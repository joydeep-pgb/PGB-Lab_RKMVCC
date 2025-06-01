import sys
import os
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QLabel, QLineEdit, QPushButton, 
                             QTextEdit, QFileDialog, QMenuBar, QMenu, 
                             QMessageBox, QProgressBar, QFrame)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt6.QtGui import QAction, QFont, QPalette, QColor
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class FastaProcessor(QThread):
    """Worker thread for processing FASTA files"""
    progress = pyqtSignal(str)
    finished = pyqtSignal(bool, str, int)
    
    def __init__(self, input_file, output_file):
        super().__init__()
        self.input_file = input_file
        self.output_file = output_file
    
    def run(self):
        try:
            self.progress.emit("Starting FASTA processing...")
            
            # Check if input file exists
            if not os.path.exists(self.input_file):
                self.finished.emit(False, f"Input file '{self.input_file}' not found!", 0)
                return
            
            self.progress.emit("Reading sequences...")
            sequences = []
            
            for record in SeqIO.parse(self.input_file, "fasta"):
                # Extract ID (first part of the header)
                original_id = record.id
                simplified_id = original_id.split()[0] if ' ' in original_id else original_id
                
                # Create new record with simplified ID
                new_record = SeqRecord(
                    seq=record.seq,
                    id=simplified_id,
                    description=""  # Remove description to keep only ID
                )
                
                sequences.append(new_record)
            
            self.progress.emit(f"Processing {len(sequences)} sequences...")
            
            # Write sequences to output file with 80 characters per line
            with open(self.output_file, 'w') as output_handle:
                for i, record in enumerate(sequences):
                    # Write header
                    output_handle.write(f">{record.id}\n")
                    
                    # Write sequence with 80 characters per line
                    sequence = str(record.seq)
                    for j in range(0, len(sequence), 80):
                        output_handle.write(sequence[j:j+80] + "\n")
                    
                    # Update progress periodically
                    if i % 100 == 0:
                        self.progress.emit(f"Processed {i+1}/{len(sequences)} sequences...")
            
            self.finished.emit(True, f"Successfully processed {len(sequences)} sequences", len(sequences))
            
        except Exception as e:
            self.finished.emit(False, f"Error processing FASTA file: {str(e)}", 0)

class FastaHeaderFixerGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.is_dark_theme = False
        self.processor_thread = None
        self.init_ui()
        self.apply_light_theme()
        
    def init_ui(self):
        self.setWindowTitle("Phytozome FASTA Header Fixer")
        self.setGeometry(100, 100, 800, 600)
        
        # Create menu bar
        self.create_menu_bar()
        
        # Create central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # Main layout
        layout = QVBoxLayout(central_widget)
        layout.setSpacing(15)
        layout.setContentsMargins(20, 20, 20, 20)
        
        # Title
        title_label = QLabel("FASTA Header Fixer")
        title_font = QFont()
        title_font.setPointSize(18)
        title_font.setBold(True)
        title_label.setFont(title_font)
        title_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(title_label)
        
        # Description
        desc_label = QLabel("This tool simplifies FASTA headers to contain only sequence IDs and formats sequences to 80 characters per line.")
        desc_label.setWordWrap(True)
        desc_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(desc_label)
        
        # Separator
        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        layout.addWidget(line)
        
        # Input file section
        input_layout = QHBoxLayout()
        input_label = QLabel("Input FASTA File:")
        input_label.setMinimumWidth(120)
        self.input_path = QLineEdit()
        self.input_path.setPlaceholderText("Select input FASTA file...")
        self.browse_input_btn = QPushButton("Browse")
        self.browse_input_btn.clicked.connect(self.browse_input_file)
        
        input_layout.addWidget(input_label)
        input_layout.addWidget(self.input_path)
        input_layout.addWidget(self.browse_input_btn)
        layout.addLayout(input_layout)
        
        # Output file section
        output_layout = QHBoxLayout()
        output_label = QLabel("Output FASTA File:")
        output_label.setMinimumWidth(120)
        self.output_path = QLineEdit()
        self.output_path.setPlaceholderText("Select output FASTA file...")
        self.browse_output_btn = QPushButton("Browse")
        self.browse_output_btn.clicked.connect(self.browse_output_file)
        
        output_layout.addWidget(output_label)
        output_layout.addWidget(self.output_path)
        output_layout.addWidget(self.browse_output_btn)
        layout.addLayout(output_layout)
        
        # Process button
        self.process_btn = QPushButton("Process FASTA File")
        self.process_btn.setMinimumHeight(40)
        self.process_btn.clicked.connect(self.process_fasta)
        layout.addWidget(self.process_btn)
        
        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        layout.addWidget(self.progress_bar)
        
        # Log output
        log_label = QLabel("Processing Log:")
        layout.addWidget(log_label)
        
        self.log_output = QTextEdit()
        self.log_output.setMaximumHeight(200)
        self.log_output.setReadOnly(True)
        layout.addWidget(self.log_output)
        
        # Status bar
        self.statusBar().showMessage("Ready")
        
    def create_menu_bar(self):
        menubar = self.menuBar()
        
        # File menu
        file_menu = menubar.addMenu("File")
        
        exit_action = QAction("Exit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)
        
        # View menu
        view_menu = menubar.addMenu("View")
        
        # Theme submenu
        theme_menu = QMenu("Theme", self)
        view_menu.addMenu(theme_menu)
        
        light_theme_action = QAction("Light Theme", self)
        light_theme_action.triggered.connect(self.apply_light_theme)
        theme_menu.addAction(light_theme_action)
        
        dark_theme_action = QAction("Dark Theme", self)
        dark_theme_action.triggered.connect(self.apply_dark_theme)
        theme_menu.addAction(dark_theme_action)
        
        # Help menu
        help_menu = menubar.addMenu("Help")
        
        about_action = QAction("About", self)
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)
        
    def apply_light_theme(self):
        """Apply light theme to the application"""
        self.is_dark_theme = False
        palette = QPalette()
        
        # Light theme colors
        palette.setColor(QPalette.ColorRole.Window, QColor(255, 255, 255))
        palette.setColor(QPalette.ColorRole.WindowText, QColor(0, 0, 0))
        palette.setColor(QPalette.ColorRole.Base, QColor(255, 255, 255))
        palette.setColor(QPalette.ColorRole.AlternateBase, QColor(245, 245, 245))
        palette.setColor(QPalette.ColorRole.ToolTipBase, QColor(255, 255, 220))
        palette.setColor(QPalette.ColorRole.ToolTipText, QColor(0, 0, 0))
        palette.setColor(QPalette.ColorRole.Text, QColor(0, 0, 0))
        palette.setColor(QPalette.ColorRole.Button, QColor(240, 240, 240))
        palette.setColor(QPalette.ColorRole.ButtonText, QColor(0, 0, 0))
        palette.setColor(QPalette.ColorRole.BrightText, QColor(255, 0, 0))
        palette.setColor(QPalette.ColorRole.Link, QColor(42, 130, 218))
        palette.setColor(QPalette.ColorRole.Highlight, QColor(42, 130, 218))
        palette.setColor(QPalette.ColorRole.HighlightedText, QColor(255, 255, 255))
        
        self.setPalette(palette)
        self.statusBar().showMessage("Light theme applied")
        
    def apply_dark_theme(self):
        """Apply dark theme to the application"""
        self.is_dark_theme = True
        palette = QPalette()
        
        # Dark theme colors
        palette.setColor(QPalette.ColorRole.Window, QColor(53, 53, 53))
        palette.setColor(QPalette.ColorRole.WindowText, QColor(255, 255, 255))
        palette.setColor(QPalette.ColorRole.Base, QColor(25, 25, 25))
        palette.setColor(QPalette.ColorRole.AlternateBase, QColor(53, 53, 53))
        palette.setColor(QPalette.ColorRole.ToolTipBase, QColor(0, 0, 0))
        palette.setColor(QPalette.ColorRole.ToolTipText, QColor(255, 255, 255))
        palette.setColor(QPalette.ColorRole.Text, QColor(255, 255, 255))
        palette.setColor(QPalette.ColorRole.Button, QColor(53, 53, 53))
        palette.setColor(QPalette.ColorRole.ButtonText, QColor(255, 255, 255))
        palette.setColor(QPalette.ColorRole.BrightText, QColor(255, 0, 0))
        palette.setColor(QPalette.ColorRole.Link, QColor(42, 130, 218))
        palette.setColor(QPalette.ColorRole.Highlight, QColor(42, 130, 218))
        palette.setColor(QPalette.ColorRole.HighlightedText, QColor(0, 0, 0))
        
        self.setPalette(palette)
        self.statusBar().showMessage("Dark theme applied")
        
    def browse_input_file(self):
        """Open file dialog to select input FASTA file"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, 
            "Select Input FASTA File", 
            "", 
            "FASTA Files (*.fa *.fasta *.fas);;All Files (*)"
        )
        if file_path:
            self.input_path.setText(file_path)
            self.log_message(f"Selected input file: {file_path}")
            
    def browse_output_file(self):
        """Open file dialog to select output FASTA file"""
        file_path, _ = QFileDialog.getSaveFileName(
            self, 
            "Save Output FASTA File", 
            "", 
            "FASTA Files (*.fa *.fasta *.fas);;All Files (*)"
        )
        if file_path:
            self.output_path.setText(file_path)
            self.log_message(f"Selected output file: {file_path}")
            
    def log_message(self, message):
        """Add message to log output"""
        self.log_output.append(f"[{QTimer().remainingTime()}] {message}")
        
    def process_fasta(self):
        """Start FASTA processing in a separate thread"""
        input_file = self.input_path.text().strip()
        output_file = self.output_path.text().strip()
        
        # Validate inputs
        if not input_file:
            QMessageBox.warning(self, "Warning", "Please select an input FASTA file.")
            return
            
        if not output_file:
            QMessageBox.warning(self, "Warning", "Please select an output FASTA file.")
            return
            
        if not os.path.exists(input_file):
            QMessageBox.critical(self, "Error", f"Input file does not exist: {input_file}")
            return
            
        # Disable UI during processing
        self.process_btn.setEnabled(False)
        self.browse_input_btn.setEnabled(False)
        self.browse_output_btn.setEnabled(False)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)  # Indeterminate progress
        
        # Clear log
        self.log_output.clear()
        self.log_message("Starting FASTA processing...")
        
        # Start processing thread
        self.processor_thread = FastaProcessor(input_file, output_file)
        self.processor_thread.progress.connect(self.update_progress)
        self.processor_thread.finished.connect(self.processing_finished)
        self.processor_thread.start()
        
        self.statusBar().showMessage("Processing FASTA file...")
        
    def update_progress(self, message):
        """Update progress display"""
        self.log_message(message)
        
    def processing_finished(self, success, message, seq_count):
        """Handle processing completion"""
        # Re-enable UI
        self.process_btn.setEnabled(True)
        self.browse_input_btn.setEnabled(True)
        self.browse_output_btn.setEnabled(True)
        self.progress_bar.setVisible(False)
        
        if success:
            self.log_message(f"✓ {message}")
            self.log_message(f"Output saved to: {self.output_path.text()}")
            self.statusBar().showMessage(f"Processing completed successfully! ({seq_count} sequences)")
            QMessageBox.information(self, "Success", f"Processing completed successfully!\n\n{message}\nOutput saved to: {self.output_path.text()}")
        else:
            self.log_message(f"✗ {message}")
            self.statusBar().showMessage("Processing failed!")
            QMessageBox.critical(self, "Error", f"Processing failed!\n\n{message}")
            
    def show_about(self):
        """Show about dialog"""
        QMessageBox.about(
            self, 
            "About FASTA Header Fixer",
            "FASTA Header Fixer v1.0\n\n"
            "This application simplifies FASTA headers to contain only sequence IDs "
            "and formats sequences to 80 characters per line.\n\n"
            "Built with PyQt6 and BioPython"
        )
        
    def closeEvent(self, event):
        """Handle application close event"""
        if self.processor_thread and self.processor_thread.isRunning():
            reply = QMessageBox.question(
                self, 
                "Confirm Exit", 
                "Processing is still running. Are you sure you want to exit?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No
            )
            if reply == QMessageBox.StandardButton.Yes:
                self.processor_thread.terminate()
                self.processor_thread.wait()
                event.accept()
            else:
                event.ignore()
        else:
            event.accept()

def main():
    app = QApplication(sys.argv)
    app.setApplicationName("FASTA Header Fixer")
    app.setApplicationVersion("1.0")
    
    window = FastaHeaderFixerGUI()
    window.show()
    
    sys.exit(app.exec())

if __name__ == "__main__":
    main()