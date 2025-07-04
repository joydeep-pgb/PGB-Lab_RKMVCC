#!/usr/bin/env python3

import sys
import os
from pathlib import Path
from typing import List, Tuple, Dict
from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                              QHBoxLayout, QLabel, QLineEdit, QPushButton, 
                              QTextEdit, QFileDialog, QMessageBox, QProgressBar,
                              QSplitter, QGroupBox, QGridLayout, QFrame,
                              QMenuBar, QMenu, QStyleFactory)
from PySide6.QtCore import Qt, QThread, Signal, QTimer
from PySide6.QtGui import QFont, QIcon, QPixmap, QAction, QActionGroup  # Fixed import
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class DomainExtractionWorker(QThread):
    """Worker thread for domain extraction to prevent GUI freezing."""
    
    progress = Signal(int)
    status = Signal(str)
    finished = Signal(str, int, int)
    error = Signal(str)
    
    def __init__(self, fasta_file: str, coordinates_text: str, output_file: str):
        super().__init__()
        self.fasta_file = fasta_file
        self.coordinates_text = coordinates_text
        self.output_file = output_file
    
    def run(self):
        try:
            # Load sequences
            self.status.emit("Loading FASTA sequences...")
            seq_dict = self.load_sequences()
            self.progress.emit(25)
            
            # Parse coordinates
            self.status.emit("Parsing coordinates...")
            coords = self.parse_coordinates()
            self.progress.emit(50)
            
            if not coords:
                self.error.emit("No valid coordinates found in the input text.")
                return
            
            # Extract domains
            self.status.emit("Extracting domains...")
            success_count, fail_count = self.extract_domains(seq_dict, coords)
            self.progress.emit(100)
            
            self.finished.emit(f"Extraction completed successfully!\n"
                             f"Output file: {self.output_file}",
                             success_count, fail_count)
            
        except Exception as e:
            self.error.emit(f"Error during extraction: {str(e)}")
    
    def load_sequences(self) -> Dict[str, SeqRecord]:
        """Load sequences from FASTA file."""
        try:
            return SeqIO.to_dict(SeqIO.parse(self.fasta_file, "fasta"))
        except Exception as e:
            raise Exception(f"Error loading FASTA file: {str(e)}")
    
    def parse_coordinates(self) -> List[Tuple[str, int, int]]:
        """Parse coordinates from text input."""
        coords = []
        lines = self.coordinates_text.strip().split('\n')
        
        for line_num, line in enumerate(lines, start=1):
            line = line.strip()
            if not line or line.startswith('#'):  # Skip empty lines and comments
                continue
            
            try:
                parts = line.split('\t')
                if len(parts) < 2:
                    parts = line.split()  # Try space separation
                
                if len(parts) < 2:
                    continue
                
                protein_id = parts[0].strip()
                coord_range = parts[1].strip()
                
                if '-' not in coord_range:
                    continue
                
                start_str, end_str = coord_range.split('-', 1)
                start, end = int(start_str), int(end_str)
                
                if start <= 0 or end <= 0 or start > end:
                    continue
                
                coords.append((protein_id, start, end))
                
            except (ValueError, IndexError):
                continue
        
        return coords
    
    def extract_domains(self, seq_dict: Dict[str, SeqRecord], 
                       coords: List[Tuple[str, int, int]]) -> Tuple[int, int]:
        """Extract domain sequences and write to output file."""
        success_count = 0
        fail_count = 0
        total = len(coords)
        
        with open(self.output_file, "w") as out:
            for i, (protein_id, start, end) in enumerate(coords):
                if protein_id not in seq_dict:
                    fail_count += 1
                    continue
                
                full_seq = seq_dict[protein_id].seq
                seq_length = len(full_seq)
                
                if end > seq_length:
                    fail_count += 1
                    continue
                
                domain_seq = full_seq[start - 1:end]
                out.write(f">{protein_id}_{start}-{end}\n")
                out.write(f"{domain_seq}\n")
                success_count += 1
                
                # Update progress
                progress = 50 + int((i + 1) / total * 50)
                self.progress.emit(progress)
        
        return success_count, fail_count

class DomainExtractorGUI(QMainWindow):
    """Main GUI application for domain extraction."""
    
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Protein Domain Extractor")
        self.setGeometry(100, 100, 1000, 700)
        self.setMinimumSize(800, 600)
        
        # Initialize variables
        self.fasta_file = ""
        self.worker = None
        
        self.setup_ui()
        self.setup_connections()
        self.setup_menu()
        self.apply_theme("Light")  # Apply default theme
    
    def setup_ui(self):
        """Setup the user interface."""
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # Main layout
        main_layout = QVBoxLayout(central_widget)
        main_layout.setSpacing(10)
        main_layout.setContentsMargins(10, 10, 10, 10)
        
        # Title
        title_label = QLabel("Protein Domain Extractor")
        title_font = QFont()
        title_font.setPointSize(16)
        title_font.setBold(True)
        title_label.setFont(title_font)
        title_label.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(title_label)
        
        # Separator
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        main_layout.addWidget(separator)
        
        # Input section
        input_group = QGroupBox("Input Files")
        input_layout = QGridLayout(input_group)
        input_layout.setColumnStretch(1, 1)
        
        # FASTA file selection
        input_layout.addWidget(QLabel("FASTA File:"), 0, 0)
        self.fasta_line_edit = QLineEdit()
        self.fasta_line_edit.setPlaceholderText("Select your FASTA file containing protein sequences...")
        input_layout.addWidget(self.fasta_line_edit, 0, 1)
        
        self.browse_fasta_btn = QPushButton("Browse")
        self.browse_fasta_btn.setMaximumWidth(100)
        input_layout.addWidget(self.browse_fasta_btn, 0, 2)
        
        # Output file selection
        input_layout.addWidget(QLabel("Output File:"), 1, 0)
        self.output_line_edit = QLineEdit()
        self.output_line_edit.setPlaceholderText("Select output file for extracted domains...")
        input_layout.addWidget(self.output_line_edit, 1, 1)
        
        self.browse_output_btn = QPushButton("Browse")
        self.browse_output_btn.setMaximumWidth(100)
        input_layout.addWidget(self.browse_output_btn, 1, 2)
        
        main_layout.addWidget(input_group)
        
        # Coordinates section
        coords_group = QGroupBox("Coordinates Data")
        coords_layout = QVBoxLayout(coords_group)
        
        coords_label = QLabel("Paste your coordinates data below (Format: ProteinID\\tStart-End):")
        coords_layout.addWidget(coords_label)
        
        self.coordinates_text = QTextEdit()
        self.coordinates_text.setPlaceholderText(
            "Example:\n"
            "Protein1\t10-50\n"
            "Protein2\t25-75\n"
            "Protein3\t100-200\n\n"
            "You can also use space separation:\n"
            "Protein1 10-50\n"
            "Protein2 25-75"
        )
        coords_layout.addWidget(self.coordinates_text)
        
        main_layout.addWidget(coords_group)
        
        # Control buttons
        button_layout = QHBoxLayout()
        button_layout.setContentsMargins(0, 10, 0, 10)
        
        self.extract_btn = QPushButton("Extract Domains")
        self.extract_btn.setMinimumHeight(40)
        button_layout.addWidget(self.extract_btn)
        
        self.clear_btn = QPushButton("Clear All")
        self.clear_btn.setMinimumHeight(40)
        button_layout.addWidget(self.clear_btn)
        
        main_layout.addLayout(button_layout)
        
        # Progress section
        progress_group = QGroupBox("Progress")
        progress_layout = QVBoxLayout(progress_group)
        
        self.status_label = QLabel("Ready to extract domains...")
        progress_layout.addWidget(self.status_label)
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        progress_layout.addWidget(self.progress_bar)
        
        main_layout.addWidget(progress_group)
        
        # Results section
        results_group = QGroupBox("Results")
        results_layout = QVBoxLayout(results_group)
        
        self.results_text = QTextEdit()
        self.results_text.setReadOnly(True)
        self.results_text.setMinimumHeight(150)
        results_layout.addWidget(self.results_text)
        
        main_layout.addWidget(results_group)
        
        # Add stretch to push everything up
        main_layout.addStretch(1)
    
    def setup_menu(self):
        """Set up the professional menu bar."""
        menu_bar = QMenuBar(self)
        self.setMenuBar(menu_bar)
        
        # File menu
        file_menu = menu_bar.addMenu("&File")
        
        open_fasta_action = QAction("&Open FASTA...", self)
        open_fasta_action.setShortcut("Ctrl+O")
        open_fasta_action.triggered.connect(self.browse_fasta_file)
        file_menu.addAction(open_fasta_action)
        
        save_output_action = QAction("&Save Output As...", self)
        save_output_action.setShortcut("Ctrl+S")
        save_output_action.triggered.connect(self.browse_output_file)
        file_menu.addAction(save_output_action)
        
        file_menu.addSeparator()
        
        exit_action = QAction("&Exit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)
        
        # Edit menu
        edit_menu = menu_bar.addMenu("&Edit")
        
        clear_action = QAction("&Clear All", self)
        clear_action.setShortcut("Ctrl+Shift+C")
        clear_action.triggered.connect(self.clear_all)
        edit_menu.addAction(clear_action)
        
        # View menu
        view_menu = menu_bar.addMenu("&View")
        
        # Theme submenu
        theme_menu = view_menu.addMenu("&Theme")
        theme_group = QActionGroup(self)
        
        light_theme = QAction("&Light", self, checkable=True)
        light_theme.setChecked(True)
        light_theme.triggered.connect(lambda: self.apply_theme("Light"))
        theme_menu.addAction(light_theme)
        theme_group.addAction(light_theme)
        
        dark_theme = QAction("&Dark", self, checkable=True)
        dark_theme.triggered.connect(lambda: self.apply_theme("Dark"))
        theme_menu.addAction(dark_theme)
        theme_group.addAction(dark_theme)
        
        blue_theme = QAction("&Blue", self, checkable=True)
        blue_theme.triggered.connect(lambda: self.apply_theme("Blue"))
        theme_menu.addAction(blue_theme)
        theme_group.addAction(blue_theme)
        
        # Help menu
        help_menu = menu_bar.addMenu("&Help")
        
        about_action = QAction("&About", self)
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)
    
    def setup_connections(self):
        """Setup signal-slot connections."""
        self.browse_fasta_btn.clicked.connect(self.browse_fasta_file)
        self.browse_output_btn.clicked.connect(self.browse_output_file)
        self.extract_btn.clicked.connect(self.extract_domains)
        self.clear_btn.clicked.connect(self.clear_all)
    
    def apply_theme(self, theme_name: str):
        """Apply a selected theme to the application."""
        if theme_name == "Light":
            # Light theme
            self.setStyleSheet("""
                QMainWindow, QWidget {
                    background-color: #f0f0f0;
                    color: #000000;
                }
                QGroupBox {
                    font-weight: bold;
                    border: 1px solid #aaaaaa;
                    border-radius: 5px;
                    margin-top: 1ex;
                    padding-top: 10px;
                    background-color: #ffffff;
                }
                QGroupBox::title {
                    subcontrol-origin: margin;
                    left: 10px;
                    padding: 0 5px;
                }
                QPushButton {
                    background-color: #e0e0e0;
                    border: 1px solid #aaaaaa;
                    border-radius: 4px;
                    padding: 5px 10px;
                }
                QPushButton:hover {
                    background-color: #d0d0d0;
                }
                QPushButton:pressed {
                    background-color: #c0c0c0;
                }
                QLineEdit, QTextEdit {
                    background-color: #ffffff;
                    border: 1px solid #aaaaaa;
                    border-radius: 3px;
                    padding: 3px;
                }
                QProgressBar {
                    border: 1px solid #aaaaaa;
                    border-radius: 3px;
                    text-align: center;
                }
                QProgressBar::chunk {
                    background-color: #4CAF50;
                    width: 10px;
                }
                QMenuBar {
                    background-color: #e0e0e0;
                    padding: 2px;
                    border-bottom: 1px solid #cccccc;
                }
                QMenuBar::item {
                    padding: 5px 10px;
                }
                QMenuBar::item:selected {
                    background-color: #d0d0d0;
                }
                QMenu {
                    background-color: #ffffff;
                    border: 1px solid #cccccc;
                }
                QMenu::item:selected {
                    background-color: #e0e0e0;
                }
            """)
            self.extract_btn.setStyleSheet("""
                QPushButton {
                    background-color: #4CAF50;
                    color: white;
                    font-weight: bold;
                }
                QPushButton:hover {
                    background-color: #45a049;
                }
                QPushButton:disabled {
                    background-color: #cccccc;
                    color: #666666;
                }
            """)
            self.clear_btn.setStyleSheet("""
                QPushButton {
                    background-color: #f44336;
                    color: white;
                    font-weight: bold;
                }
                QPushButton:hover {
                    background-color: #da190b;
                }
            """)
        
        elif theme_name == "Dark":
            # Dark theme
            self.setStyleSheet("""
                QMainWindow, QWidget {
                    background-color: #303030;
                    color: #ffffff;
                }
                QGroupBox {
                    font-weight: bold;
                    border: 1px solid #555555;
                    border-radius: 5px;
                    margin-top: 1ex;
                    padding-top: 10px;
                    background-color: #404040;
                }
                QGroupBox::title {
                    color: #ffffff;
                    subcontrol-origin: margin;
                    left: 10px;
                    padding: 0 5px;
                }
                QPushButton {
                    background-color: #505050;
                    color: #ffffff;
                    border: 1px solid #666666;
                    border-radius: 4px;
                    padding: 5px 10px;
                }
                QPushButton:hover {
                    background-color: #606060;
                }
                QPushButton:pressed {
                    background-color: #404040;
                }
                QLineEdit, QTextEdit {
                    background-color: #202020;
                    color: #ffffff;
                    border: 1px solid #555555;
                    border-radius: 3px;
                    padding: 3px;
                }
                QProgressBar {
                    border: 1px solid #555555;
                    border-radius: 3px;
                    text-align: center;
                    color: white;
                }
                QProgressBar::chunk {
                    background-color: #4CAF50;
                    width: 10px;
                }
                QMenuBar {
                    background-color: #353535;
                    color: #ffffff;
                    padding: 2px;
                    border-bottom: 1px solid #252525;
                }
                QMenuBar::item {
                    padding: 5px 10px;
                    color: #ffffff;
                }
                QMenuBar::item:selected {
                    background-color: #505050;
                }
                QMenu {
                    background-color: #404040;
                    color: #ffffff;
                    border: 1px solid #555555;
                }
                QMenu::item:selected {
                    background-color: #505050;
                }
            """)
            self.extract_btn.setStyleSheet("""
                QPushButton {
                    background-color: #2E7D32;
                    color: white;
                    font-weight: bold;
                }
                QPushButton:hover {
                    background-color: #1B5E20;
                }
                QPushButton:disabled {
                    background-color: #555555;
                    color: #aaaaaa;
                }
            """)
            self.clear_btn.setStyleSheet("""
                QPushButton {
                    background-color: #C62828;
                    color: white;
                    font-weight: bold;
                }
                QPushButton:hover {
                    background-color: #B71C1C;
                }
            """)
        
        elif theme_name == "Blue":
            # Blue theme
            self.setStyleSheet("""
                QMainWindow, QWidget {
                    background-color: #e0f0ff;
                    color: #000000;
                }
                QGroupBox {
                    font-weight: bold;
                    border: 1px solid #90c0ff;
                    border-radius: 5px;
                    margin-top: 1ex;
                    padding-top: 10px;
                    background-color: #f0f8ff;
                }
                QGroupBox::title {
                    subcontrol-origin: margin;
                    left: 10px;
                    padding: 0 5px;
                }
                QPushButton {
                    background-color: #90c0ff;
                    border: 1px solid #70a0e0;
                    border-radius: 4px;
                    padding: 5px 10px;
                }
                QPushButton:hover {
                    background-color: #70a0e0;
                }
                QPushButton:pressed {
                    background-color: #5080c0;
                }
                QLineEdit, QTextEdit {
                    background-color: #ffffff;
                    border: 1px solid #90c0ff;
                    border-radius: 3px;
                    padding: 3px;
                }
                QProgressBar {
                    border: 1px solid #90c0ff;
                    border-radius: 3px;
                    text-align: center;
                }
                QProgressBar::chunk {
                    background-color: #4CAF50;
                    width: 10px;
                }
                QMenuBar {
                    background-color: #c0e0ff;
                    padding: 2px;
                    border-bottom: 1px solid #a0c0e0;
                }
                QMenuBar::item {
                    padding: 5px 10px;
                }
                QMenuBar::item:selected {
                    background-color: #a0d0ff;
                }
                QMenu {
                    background-color: #f0f8ff;
                    border: 1px solid #90c0ff;
                }
                QMenu::item:selected {
                    background-color: #c0e0ff;
                }
            """)
            self.extract_btn.setStyleSheet("""
                QPushButton {
                    background-color: #2196F3;
                    color: white;
                    font-weight: bold;
                }
                QPushButton:hover {
                    background-color: #0b7dda;
                }
                QPushButton:disabled {
                    background-color: #90c0ff;
                    color: #666666;
                }
            """)
            self.clear_btn.setStyleSheet("""
                QPushButton {
                    background-color: #f44336;
                    color: white;
                    font-weight: bold;
                }
                QPushButton:hover {
                    background-color: #da190b;
                }
            """)
    
    def browse_fasta_file(self):
        """Browse for FASTA file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select FASTA File", "", 
            "FASTA Files (*.fasta *.fa *.fas);;All Files (*)"
        )
        if file_path:
            self.fasta_line_edit.setText(file_path)
    
    def browse_output_file(self):
        """Browse for output file."""
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Domain Sequences", "", 
            "FASTA Files (*.fasta *.fa *.fas);;All Files (*)"
        )
        if file_path:
            self.output_line_edit.setText(file_path)
    
    def extract_domains(self):
        """Start domain extraction process."""
        # Validate inputs
        if not self.fasta_line_edit.text().strip():
            QMessageBox.warning(self, "Input Error", "Please select a FASTA file.")
            return
        
        if not self.output_line_edit.text().strip():
            QMessageBox.warning(self, "Input Error", "Please specify an output file.")
            return
        
        if not self.coordinates_text.toPlainText().strip():
            QMessageBox.warning(self, "Input Error", "Please paste coordinates data.")
            return
        
        # Check if FASTA file exists
        if not Path(self.fasta_line_edit.text()).exists():
            QMessageBox.warning(self, "File Error", "FASTA file does not exist.")
            return
        
        # Create output directory if needed
        output_path = Path(self.output_line_edit.text())
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Disable UI during processing
        self.extract_btn.setEnabled(False)
        self.progress_bar.setVisible(True)
        self.progress_bar.setValue(0)
        self.results_text.clear()
        
        # Start worker thread
        self.worker = DomainExtractionWorker(
            self.fasta_line_edit.text(),
            self.coordinates_text.toPlainText(),
            self.output_line_edit.text()
        )
        
        self.worker.progress.connect(self.progress_bar.setValue)
        self.worker.status.connect(self.status_label.setText)
        self.worker.finished.connect(self.on_extraction_finished)
        self.worker.error.connect(self.on_extraction_error)
        
        self.worker.start()
    
    def on_extraction_finished(self, message: str, success_count: int, fail_count: int):
        """Handle successful extraction completion."""
        self.extract_btn.setEnabled(True)
        self.progress_bar.setVisible(False)
        self.status_label.setText("Extraction completed successfully!")
        
        result_text = f"{message}\n\nStatistics:\n"
        result_text += f"Successfully extracted: {success_count} domains\n"
        if fail_count > 0:
            result_text += f"Failed extractions: {fail_count} domains\n"
        
        self.results_text.setPlainText(result_text)
        
        QMessageBox.information(self, "Success", 
                              f"Domain extraction completed!\n\n"
                              f"Successfully extracted: {success_count} domains\n"
                              f"Failed extractions: {fail_count} domains")
    
    def on_extraction_error(self, error_message: str):
        """Handle extraction errors."""
        self.extract_btn.setEnabled(True)
        self.progress_bar.setVisible(False)
        self.status_label.setText("Extraction failed!")
        
        self.results_text.setPlainText(f"Error: {error_message}")
        QMessageBox.critical(self, "Error", error_message)
    
    def clear_all(self):
        """Clear all input fields."""
        self.fasta_line_edit.clear()
        self.output_line_edit.clear()
        self.coordinates_text.clear()
        self.results_text.clear()
        self.status_label.setText("Ready to extract domains...")
        self.progress_bar.setVisible(False)
    
    def show_about(self):
        """Show about dialog."""
        about_text = (
            "<b>Protein Domain Extractor</b><br><br>"
            "Version: 1.1.0<br>"
            "Developed by BioTools<br><br>"
            "This application extracts protein domains from FASTA files<br>"
            "based on user-provided coordinates.<br><br>"
            "© 2025 BioTools - All rights reserved"
        )
        QMessageBox.about(self, "About Protein Domain Extractor", about_text)

def main():
    """Main function to run the application."""
    app = QApplication(sys.argv)
    app.setApplicationName("Protein Domain Extractor")
    app.setOrganizationName("BioTools")
    
    # Set application style
    app.setStyle(QStyleFactory.create("Fusion"))
    
    window = DomainExtractorGUI()
    window.show()
    
    sys.exit(app.exec())

if __name__ == "__main__":
    main()