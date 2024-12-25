cnet <- function(x,
                 foldChange   = NULL,
                 layout = "kk",
                 termcolor="#E5C494",
                 palette=c("darkgreen", "#0AFF34", "#B3B3B3", "#FF6347", "red"),
                 colorEdge = FALSE,
                 Edgecolors="None",
                 circular = FALSE,
                 node_label = FALSE) {
  
  if (circular) {
    layout <- "linear"
    geom_edge <- geom_edge_arc
  } else {
    geom_edge <- geom_edge_link
  }
  
  geneSets <- lapply(x$Objects,function(x)strsplit(x,"/")[[1]]) #enrichplot:::extract_geneSets(x, showCategory)
  if(node_label) {
    names(geneSets)<-x$ID
  }else{
    names(geneSets)<-x$Term
  }
  g <- enrichplot:::list2graph(geneSets)
  
  #foldChange <- enrichplot:::fc_readable(x, foldChange)
  
  size <- sapply(geneSets, length)
  V(g)$size <- min(size)/2
  
  n <- length(geneSets)
  V(g)$size[1:n] <- size
  
  if (colorEdge) {
    E(g)$Edge <- rep(names(geneSets), sapply(geneSets, length))
    edge_layer <- geom_edge(aes_(color = ~Edge), alpha=.8)
  } else {
    edge_layer <- geom_edge(alpha=.8, colour='darkgrey')
  }
  
  if (!is.null(foldChange)) {
    fc <- foldChange[V(g)$name[(n+1):length(V(g))]]
    V(g)$color <- NA
    V(g)$color[(n+1):length(V(g))] <- fc
    palette <- palette#enrichplot:::fc_palette(fc)
    p <- ggraph(g, layout=layout, circular = circular) +
      edge_layer +
      geom_node_point(aes_(color=~as.numeric(as.character(color)), size=~size)) +
      scale_color_gradientn(name = "fold change", colors=palette, na.value = termcolor)
    if(Edgecolors[1]=="None"){
      p <- p
    }else{
      p <- p+scale_edge_colour_manual(name="Edges",values = Edgecolors)
    }
  } else {
    V(g)$color <- "#B3B3B3"
    V(g)$color[1:n] <- termcolor#"#E5C494"
    p <- ggraph(g, layout=layout, circular=circular) +
      edge_layer +
      geom_node_point(aes_(color=~I(color), size=~size))
  }
  
  p <- p + scale_size(range=c(5, 15), breaks=unique(round(seq(min(size), max(size), length.out=4)))) +
    theme_void()
  
  p <- p + geom_node_text(aes_(label=""), repel=T)
  
  return(p)
}