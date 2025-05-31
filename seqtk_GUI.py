import sys
from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                               QHBoxLayout, QGridLayout, QLabel, QPushButton, 
                               QLineEdit, QTextEdit, QComboBox, QFileDialog, 
                               QMessageBox, QFrame, QSplitter, QMenuBar, QDialog)
from PySide6.QtCore import Qt, QThread, Signal, QUrl
from PySide6.QtGui import QFont, QPalette, QColor, QAction, QKeySequence, QDesktopServices
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

class AboutDialog(QDialog):
    """About dialog for the application"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("About Gene Extractor")
        self.setFixedSize(400, 300)
        self.setModal(True)
        
        layout = QVBoxLayout(self)
        layout.setSpacing(15)
        
        # Title
        title = QLabel("Gene Extractor")
        title.setAlignment(Qt.AlignCenter)
        title.setFont(QFont("Arial", 16, QFont.Bold))
        layout.addWidget(title)
        
        # Version
        version = QLabel("Version 2.0")
        version.setAlignment(Qt.AlignCenter)
        layout.addWidget(version)
        
        # Description
        description = QLabel(
            "A bioinformatics tool for extracting specific gene sequences "
            "from FASTA files based on gene IDs.\n\n"
            "Features:\n"
            "• Extract multiple sequences by ID\n"
            "• Customizable line length\n"
            "• Missing ID detection\n"
            "• Dark theme interface"
        )
        description.setWordWrap(True)
        description.setAlignment(Qt.AlignCenter)
        layout.addWidget(description)
        
        # Dependencies
        deps = QLabel(
            "Built with:\n"
            "• PySide6 (Qt for Python)\n"
            "• BioPython"
        )
        deps.setAlignment(Qt.AlignCenter)
        layout.addWidget(deps)
        
        # Close button
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.accept)
        layout.addWidget(close_btn)

class GeneExtractorApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Extract Subsequences from FASTA/Q Files")
        self.setGeometry(100, 100, 900, 700)
        self.setMinimumSize(800, 600)
        
        # Set up the UI
        self.setup_menubar()
        self.setup_ui()
        self.apply_dark_theme()
        
        # Worker thread
        self.worker = None
        
    def setup_menubar(self):
        """Create and configure the menu bar"""
        menubar = self.menuBar()
        
        # File menu
        file_menu = menubar.addMenu("&File")
        
        # Open FASTA file action
        open_fasta_action = QAction("&Open FASTA File...", self)
        open_fasta_action.setShortcut(QKeySequence.StandardKey.Open)
        open_fasta_action.setStatusTip("Open a FASTA file for processing")
        open_fasta_action.triggered.connect(self.browse_fasta)
        file_menu.addAction(open_fasta_action)
        
        # Save output as action
        save_output_action = QAction("&Save Output As...", self)
        save_output_action.setShortcut(QKeySequence.StandardKey.SaveAs)
        save_output_action.setStatusTip("Choose output file location")
        save_output_action.triggered.connect(self.browse_output)
        file_menu.addAction(save_output_action)
        
        file_menu.addSeparator()
        
        # Recent files submenu (placeholder for future implementation)
        recent_menu = file_menu.addMenu("Recent Files")
        no_recent_action = QAction("No recent files", self)
        no_recent_action.setEnabled(False)
        recent_menu.addAction(no_recent_action)
        
        file_menu.addSeparator()
        
        # Exit action
        exit_action = QAction("E&xit", self)
        exit_action.setShortcut(QKeySequence.StandardKey.Quit)
        exit_action.setStatusTip("Exit the application")
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)
        
        # Edit menu
        edit_menu = menubar.addMenu("&Edit")
        
        # Clear all action
        clear_action = QAction("&Clear All", self)
        clear_action.setShortcut(QKeySequence("Ctrl+L"))
        clear_action.setStatusTip("Clear all fields")
        clear_action.triggered.connect(self.clear_all)
        edit_menu.addAction(clear_action)
        
        edit_menu.addSeparator()
        
        # Copy gene IDs action
        copy_ids_action = QAction("Copy &Gene IDs", self)
        copy_ids_action.setShortcut(QKeySequence.StandardKey.Copy)
        copy_ids_action.setStatusTip("Copy gene IDs to clipboard")
        copy_ids_action.triggered.connect(self.copy_gene_ids)
        edit_menu.addAction(copy_ids_action)
        
        # Paste gene IDs action
        paste_ids_action = QAction("&Paste Gene IDs", self)
        paste_ids_action.setShortcut(QKeySequence.StandardKey.Paste)
        paste_ids_action.setStatusTip("Paste gene IDs from clipboard")
        paste_ids_action.triggered.connect(self.paste_gene_ids)
        edit_menu.addAction(paste_ids_action)
        
        # Tools menu
        tools_menu = menubar.addMenu("&Tools")
        
        # Extract sequences action
        extract_action = QAction("&Extract Sequences", self)
        extract_action.setShortcut(QKeySequence("Ctrl+E"))
        extract_action.setStatusTip("Start sequence extraction")
        extract_action.triggered.connect(self.extract_sequences)
        tools_menu.addAction(extract_action)
        
        tools_menu.addSeparator()
        
        # Validate FASTA action
        validate_action = QAction("&Validate FASTA File", self)
        validate_action.setShortcut(QKeySequence("Ctrl+V"))
        validate_action.setStatusTip("Validate the selected FASTA file")
        validate_action.triggered.connect(self.validate_fasta)
        tools_menu.addAction(validate_action)
        
        # Count sequences action
        count_action = QAction("&Count Sequences", self)
        count_action.setShortcut(QKeySequence("Ctrl+U"))
        count_action.setStatusTip("Count sequences in FASTA file")
        count_action.triggered.connect(self.count_sequences)
        tools_menu.addAction(count_action)
        
        # Settings menu
        settings_menu = menubar.addMenu("&Settings")
        
        # Theme submenu
        theme_menu = settings_menu.addMenu("&Theme")
        
        dark_theme_action = QAction("&Dark Theme", self)
        dark_theme_action.setCheckable(True)
        dark_theme_action.setChecked(True)
        dark_theme_action.triggered.connect(self.toggle_theme)
        theme_menu.addAction(dark_theme_action)
        
        light_theme_action = QAction("&Light Theme", self)
        light_theme_action.setCheckable(True)
        light_theme_action.triggered.connect(self.toggle_theme)
        theme_menu.addAction(light_theme_action)
        
        # Help menu
        help_menu = menubar.addMenu("&Help")
        
        # User guide action
        guide_action = QAction("&User Guide", self)
        guide_action.setShortcut(QKeySequence.StandardKey.HelpContents)
        guide_action.setStatusTip("Show user guide")
        guide_action.triggered.connect(self.show_user_guide)
        help_menu.addAction(guide_action)
        
        # Keyboard shortcuts action
        shortcuts_action = QAction("&Keyboard Shortcuts", self)
        shortcuts_action.setShortcut(QKeySequence("Ctrl+?"))
        shortcuts_action.setStatusTip("Show keyboard shortcuts")
        shortcuts_action.triggered.connect(self.show_shortcuts)
        help_menu.addAction(shortcuts_action)
        
        help_menu.addSeparator()
        
        # About action
        about_action = QAction("&About", self)
        about_action.setStatusTip("About this application")
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)
        
        # Store actions for later reference
        self.extract_action = extract_action
        self.dark_theme_action = dark_theme_action
        self.light_theme_action = light_theme_action
        
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
            QMenuBar {
                background-color: #3c3f41;
                color: white;
                border-bottom: 1px solid #555;
                padding: 2px;
            }
            QMenuBar::item {
                background-color: transparent;
                padding: 5px 10px;
                border-radius: 3px;
            }
            QMenuBar::item:selected {
                background-color: #4e5254;
            }
            QMenu {
                background-color: #3c3f41;
                color: white;
                border: 1px solid #555;
                padding: 2px;
            }
            QMenu::item {
                padding: 5px 20px;
                border-radius: 3px;
            }
            QMenu::item:selected {
                background-color: #4e5254;
            }
            QMenu::separator {
                height: 1px;
                background-color: #555;
                margin: 2px 5px;
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
            QDialog {
                background-color: #2e2e2e;
                color: white;
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
        self.extract_action.setEnabled(False)
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
        self.extract_action.setEnabled(True)
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
        self.extract_action.setEnabled(True)
        self.extract_btn.setText("Extract Sequences")
        self.results_label.setText("Ready to extract sequences")
        
        QMessageBox.critical(self, "Error", f"Error during extraction:\n{error_message}")
        
    def copy_gene_ids(self):
        """Copy gene IDs to clipboard"""
        clipboard = QApplication.clipboard()
        clipboard.setText(self.id_text.toPlainText())
        
    def paste_gene_ids(self):
        """Paste gene IDs from clipboard"""
        clipboard = QApplication.clipboard()
        self.id_text.setPlainText(clipboard.text())
        
    def validate_fasta(self):
        """Validate the selected FASTA file"""
        fasta_file = self.fasta_entry.text().strip()
        if not fasta_file:
            QMessageBox.warning(self, "Input Error", "Please select a FASTA file first")
            return
            
        try:
            count = 0
            for record in SeqIO.parse(fasta_file, "fasta"):
                count += 1
                if count > 10:  # Just check first 10 records for validation
                    break
            
            if count > 0:
                QMessageBox.information(self, "Validation Result", 
                                      f"FASTA file appears to be valid.\nFound {count}+ sequences.")
            else:
                QMessageBox.warning(self, "Validation Result", 
                                  "No valid sequences found in the file.")
        except Exception as e:
            QMessageBox.critical(self, "Validation Error", 
                               f"Error validating FASTA file:\n{str(e)}")
            
    def count_sequences(self):
        """Count total sequences in the FASTA file"""
        fasta_file = self.fasta_entry.text().strip()
        if not fasta_file:
            QMessageBox.warning(self, "Input Error", "Please select a FASTA file first")
            return
            
        try:
            count = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))
            QMessageBox.information(self, "Sequence Count", 
                                  f"Total sequences in file: {count}")
        except Exception as e:
            QMessageBox.critical(self, "Count Error", 
                               f"Error counting sequences:\n{str(e)}")
            
    def toggle_theme(self):
        """Toggle between dark and light themes"""
        sender = self.sender()
        if sender == self.dark_theme_action:
            self.light_theme_action.setChecked(False)
            self.apply_dark_theme()
        else:
            self.dark_theme_action.setChecked(False)
            self.apply_light_theme()
            
    def apply_light_theme(self):
        """Apply light theme to the application"""
        self.setStyleSheet("""
            QMainWindow {
                background-color: #f0f0f0;
                color: black;
            }
            QWidget {
                background-color: #f0f0f0;
                color: black;
            }
            QMenuBar {
                background-color: #e0e0e0;
                color: black;
                border-bottom: 1px solid #ccc;
                padding: 2px;
            }
            QMenuBar::item {
                background-color: transparent;
                padding: 5px 10px;
                border-radius: 3px;
            }
            QMenuBar::item:selected {
                background-color: #d0d0d0;
            }
            QMenu {
                background-color: white;
                color: black;
                border: 1px solid #ccc;
                padding: 2px;
            }
            QMenu::item {
                padding: 5px 20px;
                border-radius: 3px;
            }
            QMenu::item:selected {
                background-color: #e0e0e0;
            }
            QMenu::separator {
                height: 1px;
                background-color: #ccc;
                margin: 2px 5px;
            }
            QLabel {
                color: black;
                padding: 2px;
            }
            QTextEdit {
                background-color: white;
                color: black;
                border: 1px solid #ccc;
                border-radius: 4px;
                padding: 5px;
                selection-background-color: #b0b0ff;
            }
            QLineEdit {
                background-color: white;
                color: black;
                border: 1px solid #ccc;
                border-radius: 4px;
                padding: 5px;
                selection-background-color: #b0b0ff;
            }
            QComboBox {
                background-color: white;
                color: black;
                border: 1px solid #ccc;
                border-radius: 4px;
                padding: 5px;
                min-width: 80px;
            }
            QPushButton {
                background-color: #e0e0e0;
                color: black;
                border: 1px solid #ccc;
                border-radius: 4px;
                padding: 8px 16px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #d0d0d0;
                border-color: #aaa;
            }
            QPushButton:pressed {
                background-color: #c0c0c0;
            }
            QPushButton:disabled {
                background-color: #f5f5f5;
                color: #888;
                border-color: #ddd;
            }
            QFrame {
                background-color: #f0f0f0;
            }
            QDialog {
                background-color: #f0f0f0;
                color: black;
            }
        """)
        
    def show_user_guide(self):
        """Show user guide dialog"""
        guide_text = """
<h3>Gene Extractor User Guide</h3>

<h4>Getting Started:</h4>
<ol>
<li><b>Enter Gene IDs:</b> Type or paste gene IDs in the text area, one per line</li>
<li><b>Select Input File:</b> Choose your FASTA file using File → Open or the Browse button</li>
<li><b>Choose Output Location:</b> Specify where to save extracted sequences</li>
<li><b>Set Line Length:</b> Choose the number of characters per line (default: 60)</li>
<li><b>Extract:</b> Click "Extract Sequences" or press Ctrl+E</li>
</ol>

<h4>Menu Options:</h4>
<ul>
<li><b>File Menu:</b> Open files, save output, view recent files, exit</li>
<li><b>Edit Menu:</b> Clear fields, copy/paste gene IDs</li>
<li><b>Tools Menu:</b> Extract sequences, validate FASTA, count sequences</li>
<li><b>Settings Menu:</b> Switch between dark and light themes</li>
<li><b>Help Menu:</b> View this guide, keyboard shortcuts, about info</li>
</ul>

<h4>Tips:</h4>
<ul>
<li>Use Ctrl+V to validate your FASTA file before extraction</li>
<li>Use Ctrl+U to count total sequences in your file</li>
<li>Missing IDs will be displayed in the results area</li>
<li>The application supports both .fasta and .fa file extensions</li>
</ul>
        """
        
        msg = QMessageBox(self)
        msg.setWindowTitle("User Guide")
        msg.setTextFormat(Qt.RichText)
        msg.setText(guide_text)
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec()
        
    def show_shortcuts(self):
        """Show keyboard shortcuts dialog"""
        shortcuts_text = """
<h3>Keyboard Shortcuts</h3>

<table cellpadding="5">
<tr><td><b>Ctrl+O</b></td><td>Open FASTA file</td></tr>
<tr><td><b>Ctrl+S</b></td><td>Save output as</td></tr>
<tr><td><b>Ctrl+Q</b></td><td>Exit application</td></tr>
<tr><td><b>Ctrl+L</b></td><td>Clear all fields</td></tr>
<tr><td><b>Ctrl+C</b></td><td>Copy gene IDs</td></tr>
<tr><td><b>Ctrl+V</b></td><td>Paste gene IDs / Validate FASTA</td></tr>
<tr><td><b>Ctrl+E</b></td><td>Extract sequences</td></tr>
<tr><td><b>Ctrl+U</b></td><td>Count sequences</td></tr>
<tr><td><b>Ctrl+?</b></td><td>Show shortcuts</td></tr>
<tr><td><b>F1</b></td><td>User guide</td></tr>
</table>
        """
        
        msg = QMessageBox(self)
        msg.setWindowTitle("Keyboard Shortcuts")
        msg.setTextFormat(Qt.RichText)
        msg.setText(shortcuts_text)
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec()
        
    def show_about(self):
        """Show about dialog"""
        about_dialog = AboutDialog(self)
        about_dialog.exec()
        
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