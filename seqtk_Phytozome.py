from PySide6.QtWidgets import (
    QApplication, QWidget, QLabel, QLineEdit, QPushButton,
    QTextEdit, QVBoxLayout, QHBoxLayout, QFileDialog, QComboBox,
    QProgressBar, QMessageBox
)
from PySide6.QtCore import Qt
from Bio import SeqIO
import sys, os

class GeneExtractor(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Extract subsequences by locus from Phytozome")
        self.setGeometry(100, 100, 800, 600)
        self.setAcceptDrops(True)
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()

        layout.addWidget(QLabel("Locus IDs to Extract"))
        self.id_text = QTextEdit()
        layout.addWidget(self.id_text)

        self.fasta_input = QLineEdit()
        browse_fasta_btn = QPushButton("Browse Input FASTA")
        browse_fasta_btn.clicked.connect(self.browse_fasta)

        fasta_layout = QHBoxLayout()
        fasta_layout.addWidget(QLabel("Input FASTA:"))
        fasta_layout.addWidget(self.fasta_input)
        fasta_layout.addWidget(browse_fasta_btn)
        layout.addLayout(fasta_layout)

        self.output_file = QLineEdit()
        browse_out_btn = QPushButton("Browse Output File")
        browse_out_btn.clicked.connect(self.browse_output)

        out_layout = QHBoxLayout()
        out_layout.addWidget(QLabel("Output File:"))
        out_layout.addWidget(self.output_file)
        out_layout.addWidget(browse_out_btn)
        layout.addLayout(out_layout)

        line_layout = QHBoxLayout()
        line_layout.addWidget(QLabel("Line Length:"))
        self.line_length = QComboBox()
        self.line_length.addItems(["60", "70", "80", "100", "120"])
        self.line_length.setCurrentText("60")
        line_layout.addWidget(self.line_length)
        layout.addLayout(line_layout)

        self.progress = QProgressBar()
        self.progress.setValue(0)
        layout.addWidget(self.progress)

        btn_layout = QHBoxLayout()
        extract_btn = QPushButton("Extract Sequences")
        extract_btn.clicked.connect(self.extract_sequences)
        clear_btn = QPushButton("Clear")
        clear_btn.clicked.connect(self.clear_all)
        btn_layout.addWidget(extract_btn)
        btn_layout.addWidget(clear_btn)
        layout.addLayout(btn_layout)

        self.status_label = QLabel("Ready to extract sequences")
        layout.addWidget(self.status_label)

        self.setLayout(layout)

    def browse_fasta(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select FASTA File", "", "FASTA files (*.fasta *.fa)")
        if file:
            self.fasta_input.setText(file)

    def browse_output(self):
        file, _ = QFileDialog.getSaveFileName(self, "Save Output File", "", "FASTA files (*.fasta)")
        if file:
            self.output_file.setText(file)

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event):
        for url in event.mimeData().urls():
            path = url.toLocalFile()
            if path.lower().endswith(('.fasta', '.fa')):
                self.fasta_input.setText(path)
            elif path.lower().endswith('.txt'):
                self.output_file.setText(path)

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
            fasta_dict = {}
            for record in SeqIO.parse(fasta_file, "fasta"):
                desc = record.description
                locus = next((f.split("=")[1] for f in desc.split() if f.startswith("locus=")), None)
                if locus:
                    fasta_dict[locus] = {"header": desc, "sequence": str(record.seq)}
            self.progress.setValue(40)

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
            self.status_label.setText(msg)
            QMessageBox.information(self, "Extraction Complete", msg)
        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))
            self.progress.setValue(0)

    def clear_all(self):
        self.id_text.clear()
        self.fasta_input.clear()
        self.output_file.clear()
        self.line_length.setCurrentText("60")
        self.status_label.setText("Ready to extract sequences")
        self.progress.setValue(0)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyleSheet("""{}""".format('''{}'''))
    window = GeneExtractor()
    window.show()
    sys.exit(app.exec())
".replace("'''{}'''
