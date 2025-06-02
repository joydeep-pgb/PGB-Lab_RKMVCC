## For Ubuntu/Debian-based systems, run: sudo apt install libxcb-cursor0 libxcb-xinerama0 libxcb-icccm4 libxcb-image0 libxcb-keysyms1 libxcb-render-util0

from PySide6.QtWidgets import (
    QApplication, QWidget, QLabel, QLineEdit, QPushButton,
    QTextEdit, QVBoxLayout, QHBoxLayout, QFileDialog, QComboBox,
    QProgressBar, QMessageBox, QMenuBar, QMenu, QMainWindow, QStatusBar
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QPalette, QColor
from Bio import SeqIO
import sys, os

class GeneExtractor(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Extract subsequences by locus from Phytozome")
        self.setGeometry(100, 100, 800, 600)
        self.setAcceptDrops(True)
        
        # Create central widget and layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        self.main_layout = QVBoxLayout(central_widget)
        
        self.init_ui()
        self.create_menu()
        self.create_status_bar()
        
        # Initialize with light theme
        self.current_theme = "light"
        self.apply_theme()

    def init_ui(self):
        self.main_layout.addWidget(QLabel("Locus IDs to Extract"))
        self.id_text = QTextEdit()
        self.main_layout.addWidget(self.id_text)

        self.fasta_input = QLineEdit()
        browse_fasta_btn = QPushButton("Browse Input FASTA")
        browse_fasta_btn.clicked.connect(self.browse_fasta)

        fasta_layout = QHBoxLayout()
        fasta_layout.addWidget(QLabel("Input FASTA:"))
        fasta_layout.addWidget(self.fasta_input)
        fasta_layout.addWidget(browse_fasta_btn)
        self.main_layout.addLayout(fasta_layout)

        self.output_file = QLineEdit()
        browse_out_btn = QPushButton("Browse Output File")
        browse_out_btn.clicked.connect(self.browse_output)

        out_layout = QHBoxLayout()
        out_layout.addWidget(QLabel("Output File:"))
        out_layout.addWidget(self.output_file)
        out_layout.addWidget(browse_out_btn)
        self.main_layout.addLayout(out_layout)

        line_layout = QHBoxLayout()
        line_layout.addWidget(QLabel("Line Length:"))
        self.line_length = QComboBox()
        self.line_length.addItems(["60", "70", "80", "100", "120"])
        self.line_length.setCurrentText("60")
        line_layout.addWidget(self.line_length)
        self.main_layout.addLayout(line_layout)

        self.progress = QProgressBar()
        self.progress.setValue(0)
        self.main_layout.addWidget(self.progress)

        btn_layout = QHBoxLayout()
        extract_btn = QPushButton("Extract Sequences")
        extract_btn.clicked.connect(self.extract_sequences)
        clear_btn = QPushButton("Clear")
        clear_btn.clicked.connect(self.clear_all)
        btn_layout.addWidget(extract_btn)
        btn_layout.addWidget(clear_btn)
        self.main_layout.addLayout(btn_layout)

    def create_menu(self):
        # Create menu bar
        menu_bar = QMenuBar()
        self.setMenuBar(menu_bar)
        
        # File menu
        file_menu = menu_bar.addMenu("&File")
        
        open_fasta_action = file_menu.addAction("&Open FASTA...")
        open_fasta_action.triggered.connect(self.browse_fasta)
        
        open_ids_action = file_menu.addAction("Open &IDs...")
        open_ids_action.triggered.connect(self.open_ids_file)
        
        file_menu.addSeparator()
        
        set_output_action = file_menu.addAction("Set &Output File...")
        set_output_action.triggered.connect(self.browse_output)
        
        file_menu.addSeparator()
        
        exit_action = file_menu.addAction("&Exit")
        exit_action.triggered.connect(self.close)
        
        # Theme menu
        theme_menu = menu_bar.addMenu("&Theme")
        
        light_theme_action = theme_menu.addAction("&Light Mode")
        light_theme_action.triggered.connect(lambda: self.set_theme("light"))
        
        dark_theme_action = theme_menu.addAction("&Dark Mode")
        dark_theme_action.triggered.connect(lambda: self.set_theme("dark"))
        
        # Help menu
        help_menu = menu_bar.addMenu("&Help")
        
        about_action = help_menu.addAction("&About")
        about_action.triggered.connect(self.show_about)
        
        about_qt_action = help_menu.addAction("About &Qt")
        about_qt_action.triggered.connect(QApplication.aboutQt)

    def create_status_bar(self):
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready to extract sequences")

    def set_theme(self, theme_name):
        self.current_theme = theme_name
        self.apply_theme()
        self.status_bar.showMessage(f"Switched to {theme_name} theme")

    def apply_theme(self):
        if self.current_theme == "dark":
            self.apply_dark_theme()
        else:
            self.apply_light_theme()

    def apply_dark_theme(self):
        dark_palette = QPalette()
        
        # Base colors
        dark_palette.setColor(QPalette.Window, QColor(53, 53, 53))
        dark_palette.setColor(QPalette.WindowText, Qt.white)
        dark_palette.setColor(QPalette.Base, QColor(35, 35, 35))
        dark_palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
        dark_palette.setColor(QPalette.ToolTipBase, Qt.white)
        dark_palette.setColor(QPalette.ToolTipText, Qt.white)
        dark_palette.setColor(QPalette.Text, Qt.white)
        dark_palette.setColor(QPalette.Button, QColor(53, 53, 53))
        dark_palette.setColor(QPalette.ButtonText, Qt.white)
        dark_palette.setColor(QPalette.BrightText, Qt.red)
        dark_palette.setColor(QPalette.Link, QColor(42, 130, 218))
        dark_palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
        dark_palette.setColor(QPalette.HighlightedText, Qt.black)
        
        # Disabled colors
        dark_palette.setColor(QPalette.Disabled, QPalette.Text, QColor(160, 160, 160))
        dark_palette.setColor(QPalette.Disabled, QPalette.ButtonText, QColor(160, 160, 160))
        
        # Set the palette
        QApplication.instance().setPalette(dark_palette)
        
        # Additional stylesheet for fine-tuning
        self.setStyleSheet("""
            QMenuBar {
                background-color: #323232;
                color: white;
            }
            QMenuBar::item:selected {
                background-color: #505050;
            }
            QMenu {
                background-color: #424242;
                color: white;
                border: 1px solid #323232;
            }
            QMenu::item:selected {
                background-color: #505050;
            }
            QStatusBar {
                background-color: #323232;
                color: white;
            }
            QComboBox, QLineEdit, QTextEdit, QSpinBox {
                background-color: #2D2D2D;
                color: white;
                border: 1px solid #555555;
            }
            QPushButton {
                background-color: #505050;
                color: white;
                border: 1px solid #555555;
                padding: 5px;
            }
            QPushButton:hover {
                background-color: #606060;
            }
            QPushButton:pressed {
                background-color: #404040;
            }
            QProgressBar {
                border: 1px solid #555555;
                border-radius: 3px;
                text-align: center;
            }
            QProgressBar::chunk {
                background-color: #4CAF50;
            }
        """)

    def apply_light_theme(self):
        # Reset to default palette
        QApplication.instance().setPalette(QApplication.style().standardPalette())
        self.setStyleSheet("")
        
        # Apply minimal styling for consistency
        self.setStyleSheet("""
            QStatusBar {
                background-color: #F0F0F0;
                color: #333333;
                border-top: 1px solid #CCCCCC;
            }
            QPushButton {
                padding: 5px;
            }
            QProgressBar {
                border: 1px solid #CCCCCC;
                border-radius: 3px;
                text-align: center;
            }
            QProgressBar::chunk {
                background-color: #4CAF50;
            }
        """)

    def browse_fasta(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select FASTA File", "", "FASTA files (*.fasta *.fa)")
        if file:
            self.fasta_input.setText(file)
            self.status_bar.showMessage(f"Loaded FASTA file: {os.path.basename(file)}")

    def open_ids_file(self):
        file, _ = QFileDialog.getOpenFileName(self, "Open Gene IDs", "", "Text files (*.txt)")
        if file:
            try:
                with open(file, 'r') as f:
                    self.id_text.setPlainText(f.read())
                self.status_bar.showMessage(f"Loaded gene IDs from: {os.path.basename(file)}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Could not read file:\n{str(e)}")

    def browse_output(self):
        file, _ = QFileDialog.getSaveFileName(self, "Save Output File", "", "FASTA files (*.fasta)")
        if file:
            self.output_file.setText(file)
            self.status_bar.showMessage(f"Output set to: {os.path.basename(file)}")

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event):
        for url in event.mimeData().urls():
            path = url.toLocalFile()
            if path.lower().endswith(('.fasta', '.fa')):
                self.fasta_input.setText(path)
                self.status_bar.showMessage(f"Loaded FASTA file: {os.path.basename(path)}")
            elif path.lower().endswith('.txt'):
                try:
                    with open(path, 'r') as f:
                        self.id_text.setPlainText(f.read())
                    self.status_bar.showMessage(f"Loaded gene IDs from: {os.path.basename(path)}")
                except Exception as e:
                    QMessageBox.critical(self, "Error", f"Could not read file:\n{str(e)}")

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

        try:
            self.progress.setValue(10)
            self.status_bar.showMessage("Reading FASTA file...")
            QApplication.processEvents()  # Update UI
            
            fasta_dict = {}
            for record in SeqIO.parse(fasta_file, "fasta"):
                desc = record.description
                locus = next((f.split("=")[1] for f in desc.split() if f.startswith("locus=")), None)
                if locus:
                    fasta_dict[locus] = {"header": desc, "sequence": str(record.seq)}
            
            self.progress.setValue(40)
            self.status_bar.showMessage("Extracting sequences...")
            QApplication.processEvents()  # Update UI
            
            extracted = {gid: fasta_dict[gid] for gid in gene_ids if gid in fasta_dict}
            with open(output_file, 'w') as f:
                for data in extracted.values():
                    f.write(f">{data['header']}\n")
                    seq = data['sequence']
                    for i in range(0, len(seq), line_len):
                        f.write(seq[i:i+line_len] + '\n')

            self.progress.setValue(100)
            missing = gene_ids - extracted.keys()
            msg = f"Extracted {len(extracted)} sequences. Missing {len(missing)}."
            self.status_bar.showMessage(msg)
            QMessageBox.information(self, "Extraction Complete", msg)
        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))
            self.progress.setValue(0)
            self.status_bar.showMessage("Error during extraction")

    def clear_all(self):
        self.id_text.clear()
        self.fasta_input.clear()
        self.output_file.clear()
        self.line_length.setCurrentText("60")
        self.progress.setValue(0)
        self.status_bar.showMessage("Ready to extract sequences")

    def show_about(self):
        QMessageBox.about(self, "About Gene Extractor",
                         "Gene Extractor v1.0\n\n"
                         "Extract subsequences by locus from Phytozome files.\n"
                         "Developed with PySide6 and Biopython.\n\n"
                         f"Current theme: {self.current_theme.capitalize()} Mode")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = GeneExtractor()
    window.show()
    sys.exit(app.exec())