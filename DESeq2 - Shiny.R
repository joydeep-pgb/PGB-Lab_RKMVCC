# Enhanced single-file Shiny app for DESeq2 differential expression analysis
# Save as app.R and run with: shiny::runApp('path/to/app.R') or open in RStudio and click "Run App"

library(shiny)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(DT)
library(vsn)
# Optional / enhance visuals
suppressPackageStartupMessages({
  haveEV <- requireNamespace('EnhancedVolcano', quietly=TRUE)
})

# Helper: read uploaded file whether csv/tsv
read_table_auto <- function(file){
  ext <- tools::file_ext(file$datapath)
  if (ext %in% c('csv','CSV')) read.csv(file$datapath, row.names=1, check.names=FALSE)
  else read.table(file$datapath, header=TRUE, sep='	', row.names=1, check.names=FALSE)
}

# Small example dataset for quick testing
make_example_data <- function(){
  set.seed(42)
  counts <- matrix(rnbinom(1000, mu=50, size=1/0.2), nrow=100, ncol=10)
  rownames(counts) <- paste0('gene', sprintf('%03d', 1:nrow(counts)))
  colnames(counts) <- paste0('S', 1:ncol(counts))
  # design: 2 conditions x 5 replicates
  meta <- data.frame(sample = colnames(counts), condition = rep(c('A','B'), each=5), batch = rep(c('X','Y'), length.out=10), stringsAsFactors=FALSE)
  list(counts = counts, meta = meta)
}

ui <- fluidPage(
  titlePanel("DESeq2 — Enhanced"),
  sidebarLayout(
    sidebarPanel(
      helpText("Upload a raw counts table (rows = genes, cols = samples) and a sample metadata table (rows = samples)."),
      fileInput('counts','Counts table (CSV or TSV)', accept = c('.csv','.tsv','.txt')),
      fileInput('meta',
                "Sample metadata (CSV or TSV). Must contain a column 'sample' or rownames matching counts column names.",
                accept = c('.csv', '.tsv', '.txt')),
      checkboxInput('header','Counts has header (colnames)', TRUE),
      textInput('sampleCol','Sample column name in metadata (leave blank to use rownames)', ''),
      textInput('design','Design formula (DESeq2 format)', value='~ condition'),
      numericInput('minCount','Minimum total counts per gene (filtering)', value=10, min=0),
      selectInput('transform','Transform for PCA/heatmap', choices=c('vst','rlog','none'), selected='vst'),
      checkboxInput('blind','Blind transformation for vst/rlog (TRUE = blind)', FALSE),
      hr(),
      h4('Contrast / Results'),
      uiOutput('factor_ui'),
      uiOutput('contrast_ui'),
      actionButton('run','Run DESeq2 / (Re)compute', class='btn-primary'),
      actionButton('runContrast','Get results for selected contrast', class='btn-success'),
      hr(),
      downloadButton('downloadRes','Download results (CSV)'),
      downloadButton('downloadExample','Download example data (zip)'),
      width=3
    ),
    mainPanel(
      tabsetPanel(
        tabPanel('Summary', verbatimTextOutput('summaryText')),
        tabPanel('Counts preview', DTOutput('countsTable')),
        tabPanel('Metadata preview', DTOutput('metaTable')),
        tabPanel('DE results', DTOutput('resultsTable')),
        tabPanel('MA plot', plotOutput('maPlot')),
        tabPanel('Enhanced Volcano', plotOutput('evPlot')),
        tabPanel('Volcano', plotOutput('volcanoPlot')),
        tabPanel('PCA', plotOutput('pcaPlot')),
        tabPanel('QC', plotOutput('qcPlots')),
        tabPanel('Heatmap', plotOutput('heatmapPlot')),
        tabPanel('Logs', verbatimTextOutput('logText'))
      )
    )
  )
)

server <- function(input, output, session){
  rv <- reactiveValues(
    counts = NULL,
    meta = NULL,
    dds = NULL,
    res_contrast = NULL,
    normMat = NULL,
    log = ''
  )
  
  # Provide example data download
  output$downloadExample <- downloadHandler(
    filename = function(){ paste0('example_counts_and_meta_', Sys.Date(), '.zip') },
    content = function(file){
      tmpdir <- tempdir(); oldwd <- setwd(tmpdir)
      ex <- make_example_data()
      write.csv(ex$counts, file='example_counts.csv', row.names=TRUE)
      write.csv(ex$meta, file='example_meta.csv', row.names=FALSE)
      zip(zipfile='example_data.zip', files=c('example_counts.csv','example_meta.csv'))
      setwd(oldwd)
      file.copy(file.path(tmpdir,'example_data.zip'), file)
    }
  )
  
  observeEvent(input$counts, {
    req(input$counts)
    tryCatch({
      df <- read_table_auto(input$counts)
      rv$counts <- as.matrix(df)
      rv$log <- paste0(rv$log, '
Loaded counts: ', nrow(df), ' genes x ', ncol(df), ' samples')
    }, error=function(e){ rv$log <- paste0(rv$log,'
Error reading counts: ', e$message) })
  })
  
  observeEvent(input$meta, {
    req(input$meta)
    tryCatch({
      md <- read_table_auto(input$meta)
      if (input$sampleCol != '' && input$sampleCol %in% colnames(md)){
        rownames(md) <- as.character(md[[input$sampleCol]])
      }
      rv$meta <- as.data.frame(md, stringsAsFactors=TRUE)
      rv$log <- paste0(rv$log, '
Loaded metadata: ', nrow(md), ' samples, ', ncol(md), ' columns')
    }, error=function(e){ rv$log <- paste0(rv$log,'
Error reading metadata: ', e$message) })
  })
  
  # If user hasn't uploaded data, populate with example on first load
  observe({
    if (is.null(rv$counts) && is.null(rv$meta)){
      ex <- make_example_data()
      rv$counts <- ex$counts
      rownames(ex$meta) <- ex$meta$sample
      rv$meta <- ex$meta
      rv$log <- paste0(rv$log, '
Loaded built-in example dataset (100 genes x 10 samples).')
    }
  })
  
  output$countsTable <- renderDT({
    req(rv$counts)
    dat <- as.data.frame(rv$counts)
    datatable(head(dat, n=50), options=list(scrollX=TRUE))
  })
  
  output$metaTable <- renderDT({
    req(rv$meta)
    datatable(rv$meta, options=list(scrollX=TRUE))
  })
  
  # Dynamic UI: factor columns for contrast
  output$factor_ui <- renderUI({
    req(rv$meta)
    meta <- rv$meta
    factorCols <- names(which(sapply(meta, function(x) length(unique(x)) > 1)))
    if (length(factorCols) == 0) return(NULL)
    selectInput('factor_col', 'Factor column (for contrasts)', choices=factorCols, selected=factorCols[1])
  })
  
  output$contrast_ui <- renderUI({
    req(rv$meta, input$factor_col)
    vals <- unique(as.character(rv$meta[[input$factor_col]]))
    tagList(
      selectInput('contrast_ref','Reference level', choices=vals, selected=vals[1]),
      selectInput('contrast_treat','Treated level', choices=vals, selected=vals[2])
    )
  })
  
  observeEvent(input$run, {
    rv$log <- paste0(rv$log, '
Running DESeq2...')
    tryCatch({
      req(rv$counts, rv$meta)
      counts <- rv$counts
      meta <- rv$meta
      sample_names <- colnames(counts)
      if (is.null(sample_names)) stop('Counts table must have column names matching sample names')
      
      if (input$sampleCol != '' && input$sampleCol %in% colnames(meta)){
        rownames(meta) <- as.character(meta[[input$sampleCol]])
      }
      
      if (!all(sample_names %in% rownames(meta))) stop('Not all count column names are present in metadata rownames. Make sure sample names match.')
      
      meta <- meta[sample_names, , drop=FALSE]
      
      keep <- rowSums(counts) >= input$minCount
      countsFilt <- counts[keep, , drop=FALSE]
      
      designFormula <- as.formula(input$design)
      dds <- DESeqDataSetFromMatrix(countData = countsFilt, colData = meta, design = designFormula)
      # ensure factors
      for (cname in colnames(colData(dds))){
        if (is.character(colData(dds)[[cname]])) colData(dds)[[cname]] <- factor(colData(dds)[[cname]])
      }
      dds <- DESeq(dds)
      rv$dds <- dds
      # precompute normalized matrix for transformations
      if (input$transform == 'vst') rv$normMat <- assay(vst(dds, blind=input$blind))
      else if (input$transform == 'rlog') rv$normMat <- assay(rlog(dds, blind=input$blind))
      else rv$normMat <- NULL
      
      rv$log <- paste0(rv$log, '
DESeq2 finished. Genes in object: ', nrow(dds))
    }, error=function(e){ rv$log <- paste0(rv$log,'
Error running DESeq2: ', e$message) })
  })
  
  # Compute results for selected contrast
  observeEvent(input$runContrast, {
    rv$log <- paste0(rv$log, '
Computing contrast results...')
    tryCatch({
      req(rv$dds, input$factor_col, input$contrast_ref, input$contrast_treat)
      # results using factor name
      res <- results(rv$dds, contrast = c(input$factor_col, input$contrast_treat, input$contrast_ref))
      resOrdered <- as.data.frame(res[order(res$padj), ])
      resOrdered$gene <- rownames(resOrdered)
      rv$res_contrast <- resOrdered
      rv$log <- paste0(rv$log, '
Contrast results ready: ', nrow(resOrdered), ' genes')
    }, error=function(e){ rv$log <- paste0(rv$log,'
Error computing contrast: ', e$message) })
  })
  
  output$summaryText <- renderPrint({
    if (is.null(rv$dds)){
      cat('No analysis run yet.
Upload data and click Run DESeq2.')
    } else {
      cat('DESeq2 object present.
')
      cat('Number of genes in object:', nrow(rv$dds), '
')
      cat('Design:', input$design, '
')
      cat('Filtering threshold (min total counts):', input$minCount, '
')
      if (!is.null(rv$res_contrast)){
        cat('
Last contrast:', paste(input$factor_col, input$contrast_treat, 'vs', input$contrast_ref), '
')
        cat('Significant (padj < 0.05):', sum(rv$res_contrast$padj < 0.05, na.rm=TRUE), '
')
      }
    }
  })
  
  output$resultsTable <- renderDT({
    req(rv$res_contrast)
    datatable(rv$res_contrast, options=list(pageLength=20, scrollX=TRUE))
  })
  
  output$maPlot <- renderPlot({
    req(rv$res_contrast)
    res <- rv$res_contrast
    plotMA(as.data.frame(res), main='MA plot', ylim=c(-5,5))
  })
  
  output$evPlot <- renderPlot({
    req(rv$res_contrast)
    if (!haveEV) {
      plot.new(); text(0.5,0.5,'EnhancedVolcano not installed. Install with BiocManager::install("EnhancedVolcano")')
      return()
    }
    df <- rv$res_contrast
    # EnhancedVolcano expects a dataframe/bioconductor result with log2FC and pvalue columns
    # Create a temporary frame
    tmp <- df
    tmp$logFC <- tmp$log2FoldChange
    tmp$pval <- tmp$pvalue
    EnhancedVolcano::EnhancedVolcano(tmp, lab = tmp$gene, x = 'logFC', y = 'pval', pCutoff = 0.05, FCcutoff = 1)
  })
  
  output$volcanoPlot <- renderPlot({
    req(rv$res_contrast)
    df <- rv$res_contrast
    df$logp <- -log10(df$padj + 1e-300)
    df$log2FoldChange[is.na(df$log2FoldChange)] <- 0
    ggplot(df, aes(x=log2FoldChange, y=logp)) +
      geom_point(alpha=0.4, size=1) +
      xlab('log2 fold change') + ylab('-log10 adjusted p-value') +
      theme_minimal() +
      geom_vline(xintercept=c(-1,1), linetype='dashed') +
      geom_hline(yintercept=-log10(0.05), linetype='dashed')
  })
  
  output$pcaPlot <- renderPlot({
    req(rv$dds)
    if (!is.null(rv$normMat)) mat <- rv$normMat else mat <- assay(rlog(rv$dds, blind=TRUE))
    pca <- prcomp(t(mat))
    percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
    df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], row.names=rownames(pca$x))
    meta <- as.data.frame(colData(rv$dds))
    factorCols <- sapply(meta, is.factor)
    colorBy <- if (any(factorCols)) names(which(factorCols))[1] else NULL
    if (!is.null(colorBy)) df$group <- as.factor(meta[[colorBy]]) else df$group <- 'all'
    ggplot(df, aes(x=PC1, y=PC2, color=group, label=rownames(df))) + geom_point(size=3) +
      xlab(paste0('PC1: ', percentVar[1], '%')) + ylab(paste0('PC2: ', percentVar[2], '%')) +
      theme_minimal()
  })
  
  output$qcPlots <- renderPlot({
    req(rv$dds)
    par(mfrow=c(2,2))
    # library sizes
    libs <- colSums(counts(rv$dds))
    barplot(libs, main='Library sizes', las=2)
    # meanSdPlot from vsn
    mat <- as.matrix(counts(rv$dds, normalized=TRUE))
    tryCatch({
      meanSdPlot(log2(mat+1), ranks=FALSE, main='Mean-SD plot (log2 normalized counts)')
    }, error=function(e){ plot.new(); text(0.5,0.5,'meanSdPlot failed') })
    # sample distances heatmap
    if (!is.null(rv$normMat)) {
      dists <- dist(t(rv$normMat))
      matDist <- as.matrix(dists)
      heatmap(matDist, main='Sample distance (heatmap)')
    } else {
      plot.new(); text(0.5,0.5,'No transformed matrix for distance heatmap')
    }
    # counts per gene distribution
    boxplot(log2(counts(rv$dds)+1), main='Log2 counts per sample', las=2)
  })
  
  output$heatmapPlot <- renderPlot({
    req(rv$normMat)
    mat <- rv$normMat
    rvgenes <- apply(mat, 1, var)
    top <- head(order(rvgenes, decreasing=TRUE), 50)
    matTop <- mat[top, , drop=FALSE]
    annotation_col <- as.data.frame(colData(rv$dds))
    pheatmap(matTop, scale='row', annotation_col=annotation_col, show_rownames=TRUE, fontsize_row=6)
  })
  
  output$logText <- renderText({ rv$log })
  
  output$downloadRes <- downloadHandler(
    filename = function(){ paste0('deseq2_contrast_results_', Sys.Date(), '.csv') },
    content = function(file){
      req(rv$res_contrast)
      write.csv(rv$res_contrast, file, row.names=FALSE)
    }
  )
}

shinyApp(ui, server)
