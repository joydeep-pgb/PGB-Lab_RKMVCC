import sys
from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                               QHBoxLayout, QGridLayout, QLabel, QPushButton, 
                               QLineEdit, QTextEdit, QComboBox, QFileDialog, 
                               QMessageBox, QFrame, QSplitter)
from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtGui import QFont, QPalette, QColor
from Bio import SeqIO

class ExtractionWorker(QThread):
    """Worker thread for sequence extraction to prevent GUI freezing"""
    finished = Signal(dict, set, int, int)
    error = Signal(str)
    
    def __init__(self, gene_ids, fasta_file, output_file, line_length):
        super().__init__()
        self.gene_ids = gene_ids
        self.fasta_file = fasta_file
        self.output_file = output_file
        self.line_length = line_length
    
    def run(self):
        try:
            # Read FASTA file
            fasta_sequences = {}
            for record in SeqIO.parse(self.fasta_file, "fasta"):
                fasta_sequences[record.id] = str(record.seq)
            
            # Extract sequences
            extracted_genes = {gene_id: fasta_sequences[gene_id] 
                             for gene_id in self.gene_ids if gene_id in fasta_sequences}
            
            # Write to output file
            with open(self.output_file, 'w') as out_file:
                for gene_id, sequence in extracted_genes.items():
                    out_file.write(f'>{gene_id}\n')
                    for i in range(0, len(sequence), self.line_length):
                        out_file.write(sequence[i:i+self.line_length] + '\n')
            
            missing = self.gene_ids - extracted_genes.keys()
            self.finished.emit(extracted_genes, missing, len(self.gene_ids), len(extracted_genes))
            
        except Exception as e:
            self.error.emit(str(e))

class GeneExtractorApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Extract Subsequences from FASTA/Q Files")
        self.setGeometry(100, 100, 900, 700)
        self.setMinimumSize(800, 600)
        
        # Set up the UI
        self.setup_ui()
        self.apply_dark_theme()
        
        # Worker thread
        self.worker = None
        
    def setup_ui(self):
        # Create central widget and main layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)
        main_layout.setSpacing(15)
        main_layout.setContentsMargins(20, 20, 20, 20)
        
        # Gene IDs section
        gene_ids_label = QLabel("Gene IDs to Extract")
        gene_ids_label.setFont(QFont("Arial", 12, QFont.Bold))
        main_layout.addWidget(gene_ids_label)
        
        self.id_text = QTextEdit()
        self.id_text.setMinimumHeight(200)
        self.id_text.setPlaceholderText("Enter gene IDs, one per line...")
        main_layout.addWidget(self.id_text)
        
        # File selection section
        file_frame = QFrame()
        file_layout = QGridLayout(file_frame)
        file_layout.setSpacing(10)
        
        # Input FASTA file
        file_layout.addWidget(QLabel("Input FASTA File:"), 0, 0)
        self.fasta_entry = QLineEdit()
        self.fasta_entry.setPlaceholderText("Select input FASTA file...")
        file_layout.addWidget(self.fasta_entry, 0, 1)
        
        browse_fasta_btn = QPushButton("Browse...")
        browse_fasta_btn.clicked.connect(self.browse_fasta)
        file_layout.addWidget(browse_fasta_btn, 0, 2)
        
        # Output file
        file_layout.addWidget(QLabel("Output File:"), 1, 0)
        self.output_entry = QLineEdit()
        self.output_entry.setPlaceholderText("Select output file location...")
        file_layout.addWidget(self.output_entry, 1, 1)
        
        browse_output_btn = QPushButton("Browse...")
        browse_output_btn.clicked.connect(self.browse_output)
        file_layout.addWidget(browse_output_btn, 1, 2)
        
        # Line length option
        file_layout.addWidget(QLabel("Line Length:"), 2, 0)
        self.line_length = QComboBox()
        self.line_length.addItems(["60", "70", "80", "100", "120"])
        self.line_length.setCurrentText("60")
        self.line_length.setMaximumWidth(100)
        file_layout.addWidget(self.line_length, 2, 1, Qt.AlignLeft)
        
        # Make the entry fields expand
        file_layout.setColumnStretch(1, 1)
        main_layout.addWidget(file_frame)
        
        # Results section
        results_frame = QFrame()
        results_layout = QVBoxLayout(results_frame)
        
        self.results_label = QLabel("Ready to extract sequences")
        self.results_label.setFont(QFont("Arial", 10, QFont.Bold))
        results_layout.addWidget(self.results_label)
        
        self.missing_ids_text = QTextEdit()
        self.missing_ids_text.setMaximumHeight(150)
        self.missing_ids_text.setReadOnly(True)
        self.missing_ids_text.setPlaceholderText("Missing IDs will be displayed here...")
        results_layout.addWidget(self.missing_ids_text)
        
        main_layout.addWidget(results_frame)
        
        # Action buttons
        button_layout = QHBoxLayout()
        
        self.extract_btn = QPushButton("Extract Sequences")
        self.extract_btn.setMinimumHeight(35)
        self.extract_btn.clicked.connect(self.extract_sequences)
        button_layout.addWidget(self.extract_btn)
        
        clear_btn = QPushButton("Clear All")
        clear_btn.setMinimumHeight(35)
        clear_btn.clicked.connect(self.clear_all)
        button_layout.addWidget(clear_btn)
        
        button_layout.addStretch()
        
        exit_btn = QPushButton("Exit")
        exit_btn.setMinimumHeight(35)
        exit_btn.clicked.connect(self.close)
        button_layout.addWidget(exit_btn)
        
        main_layout.addLayout(button_layout)
        
    def apply_dark_theme(self):
        """Apply a dark theme to the application"""
        self.setStyleSheet("""
            QMainWindow {
                background-color: #2e2e2e;
                color: white;
            }
            QWidget {
                background-color: #2e2e2e;
                color: white;
            }
            QLabel {
                color: white;
                padding: 2px;
            }
            QTextEdit {
                background-color: #3c3f41;
                color: white;
                border: 1px solid #555;
                border-radius: 4px;
                padding: 5px;
                selection-background-color: #5a5a5a;
            }
            QLineEdit {
                background-color: #3c3f41;
                color: white;
                border: 1px solid #555;
                border-radius: 4px;
                padding: 5px;
                selection-background-color: #5a5a5a;
            }
            QComboBox {
                background-color: #3c3f41;
                color: white;
                border: 1px solid #555;
                border-radius: 4px;
                padding: 5px;
                min-width: 80px;
            }
            QComboBox::drop-down {
                border: none;
                background-color: #555;
                border-radius: 2px;
            }
            QComboBox::down-arrow {
                border: 2px solid white;
                border-color: transparent transparent white transparent;
                width: 0px;
                height: 0px;
            }
            QPushButton {
                background-color: #3c3f41;
                color: white;
                border: 1px solid #555;
                border-radius: 4px;
                padding: 8px 16px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #4e5254;
                border-color: #777;
            }
            QPushButton:pressed {
                background-color: #2a2a2a;
            }
            QPushButton:disabled {
                background-color: #444;
                color: #888;
                border-color: #333;
            }
            QFrame {
                background-color: #2e2e2e;
            }
        """)
        
    def browse_fasta(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, 
            "Select Input FASTA File", 
            "", 
            "FASTA files (*.fasta *.fa);;All files (*.*)"
        )
        if file_path:
            self.fasta_entry.setText(file_path)
            
    def browse_output(self):
        file_path, _ = QFileDialog.getSaveFileName(
            self, 
            "Save Output File", 
            "", 
            "FASTA files (*.fasta);;All files (*.*)"
        )
        if file_path:
            self.output_entry.setText(file_path)
            
    def extract_sequences(self):
        # Get input values
        gene_ids_text = self.id_text.toPlainText().strip()
        gene_ids = {line.strip() for line in gene_ids_text.splitlines() if line.strip()}
        fasta_file = self.fasta_entry.text().strip()
        output_file = self.output_entry.text().strip()
        
        # Validate inputs
        if not gene_ids:
            QMessageBox.warning(self, "Input Error", "Please enter at least one Gene ID")
            return
        if not fasta_file:
            QMessageBox.warning(self, "Input Error", "Please select an input FASTA file")
            return
        if not output_file:
            QMessageBox.warning(self, "Input Error", "Please specify an output file")
            return
        
        try:
            line_length = int(self.line_length.currentText())
        except ValueError:
            QMessageBox.warning(self, "Input Error", "Line length must be an integer")
            return
            
        # Disable extract button during processing
        self.extract_btn.setEnabled(False)
        self.extract_btn.setText("Extracting...")
        self.results_label.setText("Processing...")
        
        # Start worker thread
        self.worker = ExtractionWorker(gene_ids, fasta_file, output_file, line_length)
        self.worker.finished.connect(self.on_extraction_finished)
        self.worker.error.connect(self.on_extraction_error)
        self.worker.start()
        
    def on_extraction_finished(self, extracted_genes, missing, total_requested, extracted_count):
        # Re-enable extract button
        self.extract_btn.setEnabled(True)
        self.extract_btn.setText("Extract Sequences")
        
        if not extracted_genes:
            QMessageBox.information(self, "No Results", "No matching sequences found")
            self.results_label.setText("Ready to extract sequences")
            return
            
        # Show results
        self.results_label.setText(
            f"Extraction complete! Requested: {total_requested}, "
            f"Extracted: {extracted_count}, Missing: {len(missing)}"
        )
        
        # Display missing IDs
        if missing:
            missing_text = "Missing IDs:\n" + "\n".join(sorted(missing))
        else:
            missing_text = "All requested IDs were found and extracted."
        self.missing_ids_text.setPlainText(missing_text)
        
        QMessageBox.information(
            self, 
            "Success", 
            f"Successfully extracted {extracted_count} sequences to:\n{self.output_entry.text()}"
        )
        
    def on_extraction_error(self, error_message):
        # Re-enable extract button
        self.extract_btn.setEnabled(True)
        self.extract_btn.setText("Extract Sequences")
        self.results_label.setText("Ready to extract sequences")
        
        QMessageBox.critical(self, "Error", f"Error during extraction:\n{error_message}")
        
    def clear_all(self):
        self.id_text.clear()
        self.fasta_entry.clear()
        self.output_entry.clear()
        self.line_length.setCurrentText("60")
        self.results_label.setText("Ready to extract sequences")
        self.missing_ids_text.clear()

def main():
    app = QApplication(sys.argv)
    
    # Set application properties
    app.setApplicationName("Gene Extractor")
    app.setApplicationVersion("2.0")
    app.setOrganizationName("Bioinformatics Tools")
    
    window = GeneExtractorApp()
    window.show()
    
    sys.exit(app.exec())

if __name__ == "__main__":
    main()