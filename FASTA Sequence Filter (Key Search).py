import sys
from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QLabel, QTextEdit, QLineEdit,
    QPushButton, QFileDialog, QHBoxLayout, QScrollArea, QFrame
)
from PySide6.QtGui import QFont
from PySide6.QtCore import Qt


class FastaFilterApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("🧬 FASTA Sequence Filter")
        self.resize(900, 700)

        layout = QVBoxLayout()

        # Title
        title = QLabel("🧬 FASTA Sequence Filter")
        title.setAlignment(Qt.AlignCenter)
        title.setFont(QFont("Segoe UI", 20, QFont.Bold))
        layout.addWidget(title)

        # File upload button
        self.file_btn = QPushButton("📁 Upload FASTA File (or paste text below)")
        self.file_btn.clicked.connect(self.load_file)
        layout.addWidget(self.file_btn)

        # FASTA input box
        fasta_label = QLabel("FASTA Sequences:")
        layout.addWidget(fasta_label)
        self.fasta_input = QTextEdit()
        self.fasta_input.setFont(QFont("Courier New", 10))
        layout.addWidget(self.fasta_input)

        # Search term section
        search_layout = QHBoxLayout()
        search_label = QLabel("Search Term:")
        search_layout.addWidget(search_label)
        self.search_term = QLineEdit()
        self.search_term.setText("small nuclear ribonucleoprotein")
        search_layout.addWidget(self.search_term)
        self.search_btn = QPushButton("🔍 Filter Sequences")
        self.search_btn.clicked.connect(self.filter_sequences)
        search_layout.addWidget(self.search_btn)
        layout.addLayout(search_layout)

        # Button to copy all filtered sequences
        self.copy_all_button = QPushButton("📋 Copy All Filtered Sequences")
        self.copy_all_button.setVisible(False)
        self.copy_all_button.clicked.connect(self.copy_all_sequences)
        layout.addWidget(self.copy_all_button)

        # Results section (scrollable)
        self.results_area = QScrollArea()
        self.results_area.setWidgetResizable(True)
        self.results_container = QVBoxLayout()

        results_widget = QWidget()
        results_widget.setLayout(self.results_container)
        self.results_area.setWidget(results_widget)

        layout.addWidget(self.results_area)

        self.setLayout(layout)

        # Pre-fill with sample data
        self.fasta_input.setPlainText(self.sample_fasta())

    def load_file(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Open FASTA File", "", "FASTA Files (*.fasta *.fa *.fas *.fna *.txt)")
        if file_path:
            with open(file_path, "r") as f:
                self.fasta_input.setPlainText(f.read())

    def parse_fasta(self, fasta_text):
        sequences = []
        lines = fasta_text.strip().split("\n")
        current_header = ""
        current_sequence = ""

        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                if current_header and current_sequence:
                    sequences.append({"header": current_header, "sequence": current_sequence})
                current_header = line
                current_sequence = ""
            else:
                current_sequence += line
        if current_header and current_sequence:
            sequences.append({"header": current_header, "sequence": current_sequence})
        return sequences

    def filter_sequences(self):
        fasta_text = self.fasta_input.toPlainText()
        search_term = self.search_term.text().lower()

        # Clear previous results
        for i in reversed(range(self.results_container.count())):
            widget = self.results_container.itemAt(i).widget()
            if widget:
                widget.deleteLater()

        if not fasta_text.strip():
            self.copy_all_button.setVisible(False)
            self.add_result_message("Please enter FASTA sequences to search.", error=True)
            return

        if not search_term.strip():
            self.copy_all_button.setVisible(False)
            self.add_result_message("Please enter a search term.", error=True)
            return

        try:
            sequences = self.parse_fasta(fasta_text)
            self.filtered_sequences = [seq for seq in sequences if search_term in seq["header"].lower()]

            if self.filtered_sequences:
                self.copy_all_button.setVisible(True)
            else:
                self.copy_all_button.setVisible(False)

            stats_label = QLabel(f"📊 Found {len(self.filtered_sequences)} sequence(s) containing \"{search_term}\" out of {len(sequences)} total sequences")
            stats_label.setStyleSheet("background-color: #4facfe; color: white; padding: 10px; border-radius: 8px;")
            stats_label.setAlignment(Qt.AlignCenter)
            self.results_container.addWidget(stats_label)

            if not self.filtered_sequences:
                self.add_result_message("No sequences found containing the search term in their headers.", error=True)
            else:
                for seq in self.filtered_sequences:
                    self.add_sequence_item(seq["header"], seq["sequence"], search_term)

        except Exception as e:
            self.copy_all_button.setVisible(False)
            self.add_result_message(f"Error parsing FASTA file: {str(e)}", error=True)

    def add_sequence_item(self, header, sequence, search_term):
        frame = QFrame()
        frame.setFrameShape(QFrame.StyledPanel)
        frame.setStyleSheet("background-color: rgba(248, 249, 250, 0.8); padding: 10px; border-radius: 8px;")
        vbox = QVBoxLayout(frame)

        # Highlight search term
        highlighted_header = header.replace(
            search_term,
            f"<span style='background: yellow; padding: 2px 4px; border-radius: 3px;'>{search_term}</span>"
        )

        header_label = QLabel(highlighted_header)
        header_label.setTextFormat(Qt.RichText)
        header_label.setStyleSheet("font-weight: bold; color: #2c3e50;")
        vbox.addWidget(header_label)

        seq_label = QLabel(sequence)
        seq_label.setFont(QFont("Courier New", 9))
        seq_label.setStyleSheet("background-color: rgba(255, 255, 255, 0.7); padding: 8px; border-radius: 4px;")
        seq_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        vbox.addWidget(seq_label)

        self.results_container.addWidget(frame)

    def add_result_message(self, text, error=False):
        label = QLabel(text)
        label.setWordWrap(True)
        if error:
            label.setStyleSheet("color: #e74c3c; font-style: italic; background: rgba(231, 76, 60, 0.1); padding: 10px; border: 1px dashed #e74c3c; border-radius: 8px;")
        else:
            label.setStyleSheet("color: #2c3e50; padding: 10px;")
        label.setAlignment(Qt.AlignCenter)
        self.results_container.addWidget(label)

    def copy_all_sequences(self):
        if not hasattr(self, 'filtered_sequences'):
            return
        fasta_text = ""
        for seq in self.filtered_sequences:
            fasta_text += f"{seq['header']}\n"
            fasta_text += '\n'.join([seq['sequence'][i:i + 60] for i in range(0, len(seq['sequence']), 60)]) + "\n"
        QApplication.clipboard().setText(fasta_text.strip())

    def sample_fasta(self):
        return """>Manin00g000640.1 Similar to PHL5: Myb family transcription factor PHL5 (Arabidopsis thaliana)
MNSHNINYQGQFHQKQAMSSDCNIEFTNHYSSYFATQPPWNLRINHRQALDGASGQQNLG
PGRSSSAILGQFESPASAFYATERFMGFPQYDSQTQNSPGCDLEFPAYQSSGENFSGAVS
SLEPGENIELRNALLSKLKSQKYCGNLFHSNEILETKIFPPEQIKMFGSEAVEIPDQRMN
YNTSQQEKQSLRYSFGGTSASSGGSVSSGAVLSCKTRIRWTQDLHEKFVECVNRLGGPEK
ATPKAILKLMDSEGLTIFHVKSHLQKYRMAKYTPEFAEGKSEKREILTDASQLDVKTTLQ
IKEALQLQLDVQRRLHEQIEIQRNLQLRIEEQGKQLKMLFDQQQKTNKGNLKTTQNLDIT
PEDDPLFSLEDVEVSTAENSGNTRFPSKIS
>Manin00g000650.1 Similar to UBA2C: UBP1-associated protein 2C (Arabidopsis thaliana)
MEEMKKRKMDEMNEMNNGLQLATSSSSSSQEHLRSLLEPLSKPQLVELLSRVGSQYPSIA
EEIKSIASADPAHRKLFVRGLAWNTTSETLCAAFQVHGEIEEGAVIYDKATGKSRGYGFI
TYKHMESTQNALRAPSKLIDGRMAVCNLACEGLTGATTADLAQRKLYIGGLSPEITTEML
LNFFGRHGEIEEGSVAYDKDTNDSRGFGFVTYKTVEAAKKAIEDPQKLLGGRSITVKLAD
THKGKPVQTQVSPAMINLAGIPLAAGYPQPGKVHANTPPVGYNYPQNIASYPNSSYASPP
AAAAQYPAQPPISYPPVPIKDPLGIPSTTPVGMGGYPYYLGKQ
>Manin00g000660.1 Protein of unknown function
MLVLTQVAAMKSKTPCQGFLAENCSSCPDLQQVVHVNNAMTVQLQDRLTLSLLKIVLSMG
RVRGKGKKPTMIASHEDHGSGEELKIPAFRRKGRPQKPLKAEVKEEEAEKIVEDGKDPKL
SISSKDMGNQATITKGRKRKKSTQAKEDKDMDKHESGIGLKPSTDESMKSVGFRQNCSRR
KSKPRRAAEAGVECK"""


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = FastaFilterApp()
    window.show()
    sys.exit(app.exec())
