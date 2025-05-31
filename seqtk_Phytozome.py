from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QLabel, QLineEdit, QPushButton,
    QTextEdit, QVBoxLayout, QHBoxLayout, QFileDialog, QComboBox,
    QProgressBar, QMessageBox, QMenuBar, QMenu, QGroupBox, QGridLayout,
    QSplitter, QFrame, QStatusBar, QToolBar, QSpacerItem, QSizePolicy
)
from PySide6.QtCore import Qt, QThread, Signal, QTimer
from PySide6.QtGui import QAction, QIcon, QFont, QPixmap
from Bio import SeqIO
import sys
import os


class ExtractionWorker(QThread):
    """Worker thread for sequence extraction to prevent GUI freezing"""
    progress_updated = Signal(int)
    status_updated = Signal(str)
    extraction_completed = Signal(dict)
    extraction_failed = Signal(str)

    def __init__(self, gene_ids, fasta_file, output_file, line_len):
        super().__init__()
        self.gene_ids = gene_ids
        self.fasta_file = fasta_file
        self.output_file = output_file
        self.line_len = line_len

    def run(self):
        try:
            self.progress_updated.emit(10)
            self.status_updated.emit("Loading FASTA file...")
            
            fasta_dict = {}
            for record in SeqIO.parse(self.fasta_file, "fasta"):
                desc = record.description
                locus = next((f.split("=")[1] for f in desc.split() if f.startswith("locus=")), None)
                if locus:
                    fasta_dict[locus] = {"header": desc, "sequence": str(record.seq)}

            self.progress_updated.emit(40)
            self.status_updated.emit("Extracting sequences...")

            extracted = {gid: fasta_dict[gid] for gid in self.gene_ids if gid in fasta_dict}
            
            self.progress_updated.emit(70)
            self.status_updated.emit("Writing output file...")

            with open(self.output_file, 'w') as f:
                for data in extracted.values():
                    f.write(f">{data['header']}\n")
                    seq = data['sequence']
                    for i in range(0, len(seq), self.line_len):
                        f.write(seq[i:i+self.line_len] + '\n')

            self.progress_updated.emit(100)
            missing = self.gene_ids - extracted.keys()
            result = {
                'extracted_count': len(extracted),
                'missing_count': len(missing),
                'missing_ids': missing
            }
            self.extraction_completed.emit(result)

        except Exception as e:
            self.extraction_failed.emit(str(e))


class GeneExtractor(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Gene Sequence Extractor - Phytozome")
        self.setGeometry(100, 100, 1000, 700)
        self.setAcceptDrops(True)
        
        # Initialize worker thread
        self.worker = None
        
        self.init_ui()
        self.init_menu_bar()
        self.init_status_bar()
        self.apply_styles()

    def init_ui(self):
        # Central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # Main layout
        main_layout = QVBoxLayout()
        central_widget.setLayout(main_layout)
        
        # Create splitter for better layout
        splitter = QSplitter(Qt.Vertical)
        main_layout.addWidget(splitter)
        
        # Input section
        input_group = QGroupBox("Input Configuration")
        input_layout = QGridLayout()
        
        # Locus IDs section
        input_layout.addWidget(QLabel("Locus IDs to Extract:"), 0, 0)
        self.id_text = QTextEdit()
        self.id_text.setPlaceholderText("Enter locus IDs, one per line...")
        self.id_text.setMaximumHeight(120)
        input_layout.addWidget(self.id_text, 1, 0, 1, 3)
        
        # File inputs
        input_layout.addWidget(QLabel("Input FASTA File:"), 2, 0)
        self.fasta_input = QLineEdit()
        self.fasta_input.setPlaceholderText("Select or drag & drop FASTA file...")
        browse_fasta_btn = QPushButton("Browse...")
        browse_fasta_btn.clicked.connect(self.browse_fasta)
        input_layout.addWidget(self.fasta_input, 2, 1)
        input_layout.addWidget(browse_fasta_btn, 2, 2)
        
        input_layout.addWidget(QLabel("Output File:"), 3, 0)
        self.output_file = QLineEdit()
        self.output_file.setPlaceholderText("Select output file location...")
        browse_out_btn = QPushButton("Browse...")
        browse_out_btn.clicked.connect(self.browse_output)
        input_layout.addWidget(self.output_file, 3, 1)
        input_layout.addWidget(browse_out_btn, 3, 2)
        
        # Options
        input_layout.addWidget(QLabel("Line Length:"), 4, 0)
        self.line_length = QComboBox()
        self.line_length.addItems(["60", "70", "80", "100", "120"])
        self.line_length.setCurrentText("60")
        input_layout.addWidget(self.line_length, 4, 1)
        
        input_group.setLayout(input_layout)
        splitter.addWidget(input_group)
        
        # Progress and control section
        control_group = QGroupBox("Extraction Control")
        control_layout = QVBoxLayout()
        
        # Progress bar
        self.progress = QProgressBar()
        self.progress.setValue(0)
        control_layout.addWidget(self.progress)
        
        # Buttons
        btn_layout = QHBoxLayout()
        self.extract_btn = QPushButton("Extract Sequences")
        self.extract_btn.clicked.connect(self.extract_sequences)
        self.extract_btn.setStyleSheet("QPushButton { font-weight: bold; }")
        
        clear_btn = QPushButton("Clear All")
        clear_btn.clicked.connect(self.clear_all)
        
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.clicked.connect(self.cancel_extraction)
        self.cancel_btn.setEnabled(False)
        
        btn_layout.addWidget(self.extract_btn)
        btn_layout.addWidget(clear_btn)
        btn_layout.addWidget(self.cancel_btn)
        btn_layout.addStretch()
        
        control_layout.addLayout(btn_layout)
        control_group.setLayout(control_layout)
        splitter.addWidget(control_group)
        
        # Results section
        results_group = QGroupBox("Results")
        results_layout = QVBoxLayout()
        
        self.results_text = QTextEdit()
        self.results_text.setMaximumHeight(100)
        self.results_text.setReadOnly(True)
        self.results_text.setPlaceholderText("Extraction results will appear here...")
        results_layout.addWidget(self.results_text)
        
        results_group.setLayout(results_layout)
        splitter.addWidget(results_group)
        
        # Set splitter proportions
        splitter.setSizes([300, 100, 100])

    def init_menu_bar(self):
        menubar = self.menuBar()
        
        # File menu
        file_menu = menubar.addMenu('&File')
        
        open_fasta_action = QAction('&Open FASTA File', self)
        open_fasta_action.setShortcut('Ctrl+O')
        open_fasta_action.triggered.connect(self.browse_fasta)
        file_menu.addAction(open_fasta_action)
        
        save_results_action = QAction('&Save Results', self)
        save_results_action.setShortcut('Ctrl+S')
        save_results_action.triggered.connect(self.save_results)
        file_menu.addAction(save_results_action)
        
        file_menu.addSeparator()
        
        load_ids_action = QAction('Load &IDs from File', self)
        load_ids_action.triggered.connect(self.load_ids_from_file)
        file_menu.addAction(load_ids_action)
        
        export_ids_action = QAction('&Export IDs to File', self)
        export_ids_action.triggered.connect(self.export_ids_to_file)
        file_menu.addAction(export_ids_action)
        
        file_menu.addSeparator()
        
        exit_action = QAction('E&xit', self)
        exit_action.setShortcut('Ctrl+Q')
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)
        
        # Edit menu
        edit_menu = menubar.addMenu('&Edit')
        
        clear_action = QAction('&Clear All', self)
        clear_action.setShortcut('Ctrl+L')
        clear_action.triggered.connect(self.clear_all)
        edit_menu.addAction(clear_action)
        
        clear_ids_action = QAction('Clear &IDs', self)
        clear_ids_action.triggered.connect(self.id_text.clear)
        edit_menu.addAction(clear_ids_action)
        
        # Tools menu
        tools_menu = menubar.addMenu('&Tools')
        
        validate_action = QAction('&Validate Input', self)
        validate_action.triggered.connect(self.validate_input)
        tools_menu.addAction(validate_action)
        
        count_ids_action = QAction('&Count IDs', self)
        count_ids_action.triggered.connect(self.count_ids)
        tools_menu.addAction(count_ids_action)
        
        # Help menu
        help_menu = menubar.addMenu('&Help')
        
        about_action = QAction('&About', self)
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)
        
        help_action = QAction('&Help', self)
        help_action.setShortcut('F1')
        help_action.triggered.connect(self.show_help)
        help_menu.addAction(help_action)

    def init_status_bar(self):
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready to extract sequences")

    def apply_styles(self):
        self.setStyleSheet("""
            QMainWindow {
                background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1,
                    stop: 0 #f8f9fa, stop: 1 #e9ecef);
                color: #212529;
            }
            
            QMenuBar {
                background-color: black;
                color: white;
                border: none;
                padding: 1px;
            }
            
            QMenuBar::item {
                background-color: transparent;
                padding: 8px 12px;
                border-radius: 4px;
                margin: 2px;
            }
            
            QMenuBar::item:selected {
                background-color: #495057;
            }
            
            QMenu {
                background-color: black;
                border: none;
                border-radius: 8px;
                padding: 8px;
            }
            
            QMenu::item {
                padding: 8px 20px;
                border-radius: 4px;
                margin: 2px;
            }
            
            QMenu::item:selected {
                background-color: #007bff;
                color: White;
            }
            
            QStatusBar {
                background-color: #f8f9fa;
                border-top: 1px solid #dee2e6;
                color: #6c757d;
                font-size: 12px;
            }
            
            QGroupBox {
                font-weight: 600;
                font-size: 14px;
                color: #495057;
                border: 2px solid #dee2e6;
                border-radius: 12px;
                margin-top: 15px;
                padding-top: 15px;
                background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1,
                    stop: 0 #ffffff, stop: 1 #f8f9fa);
            }
            
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 15px;
                padding: 5px 10px;
                background-color: #007bff;
                color: white;
                border-radius: 6px;
                font-weight: bold;
            }
            
            QPushButton {
                background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1,
                    stop: 0 #007bff, stop: 1 #0056b3);
                border: 1px solid #0056b3;
                color: white;
                padding: 12px 20px;
                border-radius: 8px;
                font-size: 13px;
                font-weight: 600;
                min-height: 20px;
            }
            
            QPushButton:hover {
                background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1,
                    stop: 0 #0056b3, stop: 1 #004085);
                border-color: #004085;
            }
            
            QPushButton:pressed {
                background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1,
                    stop: 0 #004085, stop: 1 #002752);
                border-color: #002752;
            }
            
            QPushButton:disabled {
                background: #6c757d;
                color: #adb5bd;
                border-color: #6c757d;
            }
            
            QPushButton[text="Clear All"] {
                background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1,
                    stop: 0 #6c757d, stop: 1 #495057);
                border-color: #495057;
            }
            
            QPushButton[text="Clear All"]:hover {
                background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1,
                    stop: 0 #5a6268, stop: 1 #343a40);
                border-color: #343a40;
            }
            
            QPushButton[text="Cancel"] {
                background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1,
                    stop: 0 #dc3545, stop: 1 #c82333);
                border-color: #c82333;
            }
            
            QPushButton[text="Cancel"]:hover {
                background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1,
                    stop: 0 #c82333, stop: 1 #a71e2a);
                border-color: #a71e2a;
            }
            
            QLineEdit {
                padding: 12px 16px;
                border: 2px solid #ced4da;
                border-radius: 8px;
                font-size: 13px;
                background-color: white;
                selection-background-color: #007bff;
                color: #495057;
            }
            
            QLineEdit:focus {
                border-color: #007bff;
                background-color: #ffffff;
            }
            
            QLineEdit:hover {
                border-color: #adb5bd;
            }
            
            QTextEdit {
                border: 2px solid #ced4da;
                border-radius: 8px;
                font-family: 'Monaco', 'Menlo', 'Ubuntu Mono', 'Consolas', 'source-code-pro', monospace;
                font-size: 12px;
                background-color: #f8f9fa;
                color: #495057;
                selection-background-color: #007bff;
                padding: 8px;
            }
            
            QTextEdit:focus {
                border-color: #007bff;
                background-color: white;
            }
            
            QComboBox {
                border: 2px solid #ced4da;
                border-radius: 8px;
                padding: 1px 1px;
                font-size: 13px;
                background-color: white;
                color: #495057;
                min-width: 100px;
            }
            
            QComboBox:focus {
                border-color: #007bff;
            }
            
            QComboBox:hover {
                border-color: #adb5bd;
            }
            
            QComboBox::drop-down {
                border: none;
                width: 20px;
                subcontrol-origin: padding;
                subcontrol-position: top right;
            }
            
            QComboBox::down-arrow {
                image: none;
                border-left: 5px solid transparent;
                border-right: 5px solid transparent;
                border-top: 5px solid #6c757d;
                margin-right: 5px;
            }
            
            QComboBox QAbstractItemView {
                border: 2px solid #dee2e6;
                border-radius: 1px;
                background-color: white;
                selection-background-color: #007bff;
                selection-color: white;
                outline: none;
            }
            
            QProgressBar {
                border: 2px solid #e9ecef;
                border-radius: 12px;
                text-align: center;
                font-weight: bold;
                font-size: 12px;
                color: #495057;
                background-color: #f8f9fa;
                height: 24px;
            }
            
            QProgressBar::chunk {
                background: qlineargradient(x1: 0, y1: 0, x2: 1, y2: 0,
                    stop: 0 #28a745, stop: 0.5 #20c997, stop: 1 #17a2b8);
                border-radius: 10px;
                margin: 2px;
            }
            
            QLabel {
                color: #495057;
                font-size: 13px;
                font-weight: 500;
            }
            
            QSplitter::handle {
                background-color: #dee2e6;
                height: 2px;
                border-radius: 1px;
            }
            
            QSplitter::handle:hover {
                background-color: #007bff;
            }
            
            /* Scrollbars */
            QScrollBar:vertical {
                background-color: #f8f9fa;
                width: 12px;
                border-radius: 6px;
                border: 1px solid #e9ecef;
            }
            
            QScrollBar::handle:vertical {
                background-color: #ced4da;
                border-radius: 5px;
                min-height: 20px;
                margin: 1px;
            }
            
            QScrollBar::handle:vertical:hover {
                background-color: #adb5bd;
            }
            
            QScrollBar::handle:vertical:pressed {
                background-color: #6c757d;
            }
            
            QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
                height: 0px;
                border: none;
                background: none;
            }
            
            QScrollBar:horizontal {
                background-color: #f8f9fa;
                height: 12px;
                border-radius: 6px;
                border: 1px solid #e9ecef;
            }
            
            QScrollBar::handle:horizontal {
                background-color: #ced4da;
                border-radius: 5px;
                min-width: 20px;
                margin: 1px;
            }
            
            QScrollBar::handle:horizontal:hover {
                background-color: #adb5bd;
            }
            
            QScrollBar::handle:horizontal:pressed {
                background-color: #6c757d;
            }
            
            QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {
                width: 0px;
                border: none;
                background: none;
            }
        """)

    def browse_fasta(self):
        file, _ = QFileDialog.getOpenFileName(
            self, "Select FASTA File", "", 
            "FASTA files (*.fasta *.fa *.fas);;All files (*.*)"
        )
        if file:
            self.fasta_input.setText(file)
            self.status_bar.showMessage(f"Loaded FASTA file: {os.path.basename(file)}")

    def browse_output(self):
        file, _ = QFileDialog.getSaveFileName(
            self, "Save Output File", "", 
            "FASTA files (*.fasta);;All files (*.*)"
        )
        if file:
            self.output_file.setText(file)

    def load_ids_from_file(self):
        file, _ = QFileDialog.getOpenFileName(
            self, "Load IDs from File", "", 
            "Text files (*.txt);;All files (*.*)"
        )
        if file:
            try:
                with open(file, 'r') as f:
                    content = f.read()
                    self.id_text.setPlainText(content)
                self.status_bar.showMessage(f"Loaded IDs from: {os.path.basename(file)}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to load file: {str(e)}")

    def export_ids_to_file(self):
        if not self.id_text.toPlainText().strip():
            QMessageBox.information(self, "Info", "No IDs to export")
            return
            
        file, _ = QFileDialog.getSaveFileName(
            self, "Export IDs to File", "", 
            "Text files (*.txt);;All files (*.*)"
        )
        if file:
            try:
                with open(file, 'w') as f:
                    f.write(self.id_text.toPlainText())
                self.status_bar.showMessage(f"Exported IDs to: {os.path.basename(file)}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to save file: {str(e)}")

    def save_results(self):
        if not self.results_text.toPlainText().strip():
            QMessageBox.information(self, "Info", "No results to save")
            return
            
        file, _ = QFileDialog.getSaveFileName(
            self, "Save Results", "", 
            "Text files (*.txt);;All files (*.*)"
        )
        if file:
            try:
                with open(file, 'w') as f:
                    f.write(self.results_text.toPlainText())
                self.status_bar.showMessage(f"Results saved to: {os.path.basename(file)}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to save results: {str(e)}")

    def validate_input(self):
        issues = []
        
        if not self.id_text.toPlainText().strip():
            issues.append("• No locus IDs provided")
        
        if not self.fasta_input.text().strip():
            issues.append("• No input FASTA file selected")
        elif not os.path.exists(self.fasta_input.text()):
            issues.append("• Input FASTA file does not exist")
        
        if not self.output_file.text().strip():
            issues.append("• No output file specified")
        
        if issues:
            QMessageBox.warning(self, "Validation Issues", 
                              "Please fix the following issues:\n\n" + "\n".join(issues))
        else:
            QMessageBox.information(self, "Validation", "All inputs are valid!")

    def count_ids(self):
        ids = [line.strip() for line in self.id_text.toPlainText().splitlines() if line.strip()]
        unique_ids = set(ids)
        
        msg = f"Total lines: {len(ids)}\nUnique IDs: {len(unique_ids)}"
        if len(ids) != len(unique_ids):
            msg += f"\nDuplicate IDs found: {len(ids) - len(unique_ids)}"
        
        QMessageBox.information(self, "ID Count", msg)

    def show_about(self):
        QMessageBox.about(self, "About Gene Sequence Extractor",
                         "Gene Sequence Extractor v2.0\n\n"
                         "A tool for extracting gene sequences from Phytozome FASTA files "
                         "based on locus IDs.\n\n"
                         "Built with PySide6 and BioPython")

    def show_help(self):
        help_text = """
Gene Sequence Extractor Help

How to use:
1. Enter locus IDs in the text area (one per line)
2. Select your input FASTA file
3. Choose where to save the output
4. Click 'Extract Sequences'

Features:
• Drag & drop FASTA files
• Threaded extraction (non-blocking)
• Progress tracking
• Detailed results reporting
• File management tools

Keyboard shortcuts:
• Ctrl+O: Open FASTA file
• Ctrl+S: Save results
• Ctrl+L: Clear all fields
• Ctrl+Q: Exit
• F1: Show this help
        """
        QMessageBox.information(self, "Help", help_text)

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event):
        for url in event.mimeData().urls():
            path = url.toLocalFile()
            if path.lower().endswith(('.fasta', '.fa', '.fas')):
                self.fasta_input.setText(path)
                self.status_bar.showMessage(f"Dropped FASTA file: {os.path.basename(path)}")
            elif path.lower().endswith('.txt'):
                try:
                    with open(path, 'r') as f:
                        content = f.read()
                        self.id_text.setPlainText(content)
                    self.status_bar.showMessage(f"Loaded IDs from: {os.path.basename(path)}")
                except Exception as e:
                    QMessageBox.warning(self, "Warning", f"Could not load text file: {str(e)}")

    def extract_sequences(self):
        gene_ids = {line.strip() for line in self.id_text.toPlainText().splitlines() if line.strip()}
        fasta_file = self.fasta_input.text()
        output_file = self.output_file.text()
        
        try:
            line_len = int(self.line_length.currentText())
        except ValueError:
            QMessageBox.warning(self, "Input Error", "Invalid line length")
            return

        if not gene_ids or not fasta_file or not output_file:
            QMessageBox.warning(self, "Input Error", "Missing required input")
            return

        if not os.path.exists(fasta_file):
            QMessageBox.warning(self, "File Error", "Input FASTA file does not exist")
            return

        # Start extraction in worker thread
        self.worker = ExtractionWorker(gene_ids, fasta_file, output_file, line_len)
        self.worker.progress_updated.connect(self.progress.setValue)
        self.worker.status_updated.connect(self.status_bar.showMessage)
        self.worker.extraction_completed.connect(self.on_extraction_completed)
        self.worker.extraction_failed.connect(self.on_extraction_failed)
        
        self.extract_btn.setEnabled(False)
        self.cancel_btn.setEnabled(True)
        self.worker.start()

    def cancel_extraction(self):
        if self.worker and self.worker.isRunning():
            self.worker.terminate()
            self.worker.wait()
            self.progress.setValue(0)
            self.status_bar.showMessage("Extraction cancelled")
            self.extract_btn.setEnabled(True)
            self.cancel_btn.setEnabled(False)

    def on_extraction_completed(self, result):
        self.extract_btn.setEnabled(True)
        self.cancel_btn.setEnabled(False)
        
        extracted_count = result['extracted_count']
        missing_count = result['missing_count']
        missing_ids = result['missing_ids']
        
        msg = f"Extraction completed successfully!\n\n"
        msg += f"Extracted sequences: {extracted_count}\n"
        msg += f"Missing sequences: {missing_count}\n"
        
        if missing_ids:
            msg += f"\nMissing IDs:\n" + "\n".join(list(missing_ids)[:10])
            if len(missing_ids) > 10:
                msg += f"\n... and {len(missing_ids) - 10} more"
        
        self.results_text.setPlainText(msg)
        self.status_bar.showMessage(f"Extracted {extracted_count} sequences, {missing_count} missing")
        
        QMessageBox.information(self, "Extraction Complete", 
                              f"Successfully extracted {extracted_count} sequences.\n"
                              f"{missing_count} sequences were not found in the input file.")

    def on_extraction_failed(self, error_msg):
        self.extract_btn.setEnabled(True)
        self.cancel_btn.setEnabled(False)
        self.progress.setValue(0)
        self.status_bar.showMessage("Extraction failed")
        
        QMessageBox.critical(self, "Extraction Error", f"Extraction failed:\n\n{error_msg}")

    def clear_all(self):
        self.id_text.clear()
        self.fasta_input.clear()
        self.output_file.clear()
        self.results_text.clear()
        self.line_length.setCurrentText("60")
        self.progress.setValue(0)
        self.status_bar.showMessage("Ready to extract sequences")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    
    # Set application properties
    app.setApplicationName("Gene Sequence Extractor")
    app.setApplicationVersion("2.0")
    app.setOrganizationName("Bioinformatics Tools")
    
    window = GeneExtractor()
    window.show()
    sys.exit(app.exec())