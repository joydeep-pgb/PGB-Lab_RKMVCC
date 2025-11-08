from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QLabel,
    QLineEdit, QPushButton, QTextEdit, QMessageBox, QFrame, QStatusBar
)
from PyQt6.QtGui import QFont, QTextCursor
from PyQt6.QtCore import Qt
import sys


class DNASequenceComplementApp(QMainWindow):
    def __init__(self):
        super().__init__()

        # === Window Properties ===
        self.setWindowTitle("DNA Sequence Complement Tool")
        self.setGeometry(300, 200, 650, 420)
        self.setStyleSheet("""
            QWidget {
                background-color: #f8fafc;
                font-family: 'Segoe UI';
                color: #202124;
            }
            QLabel {
                font-size: 16px;
                font-weight: 600;
            }
            QLineEdit {
                border: 1px solid #c0c4c7;
                border-radius: 6px;
                padding: 6px 8px;
                font-size: 15px;
            }
            QTextEdit {
                border: 1px solid #c0c4c7;
                border-radius: 6px;
                padding: 8px;
                font-size: 15px;
                background-color: #ffffff;
            }
            QPushButton {
                background-color: #0078d4;
                color: white;
                font-weight: 600;
                font-size: 15px;
                border-radius: 6px;
                padding: 8px 12px;
            }
            QPushButton:hover {
                background-color: #005a9e;
            }
        """)

        # === Central Widget ===
        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        main_layout = QVBoxLayout()
        main_layout.setSpacing(15)
        main_layout.setContentsMargins(30, 25, 30, 25)

        # --- Title Label ---
        title_label = QLabel("🧬 DNA Sequence Complement Calculator")
        title_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        title_label.setFont(QFont("Segoe UI", 18, QFont.Weight.Bold))
        main_layout.addWidget(title_label)

        # --- Separator line ---
        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        main_layout.addWidget(line)

        # --- Input Section ---
        input_label = QLabel("Enter DNA sequence:")
        main_layout.addWidget(input_label)

        self.entry = QLineEdit()
        self.entry.setPlaceholderText("Example: ATCGTAGC")
        main_layout.addWidget(self.entry)

        # --- Button ---
        calculate_button = QPushButton("Calculate Complement")
        calculate_button.clicked.connect(self.calculate_complements)
        main_layout.addWidget(calculate_button)

        # --- Result Section ---
        result_label = QLabel("Results:")
        main_layout.addWidget(result_label)

        self.result_text = QTextEdit()
        self.result_text.setReadOnly(True)
        main_layout.addWidget(self.result_text)

        central_widget.setLayout(main_layout)

        # === Status Bar ===
        self.status = QStatusBar()
        self.status.setStyleSheet("""
            QStatusBar {
                background-color: #e8eaed;
                color: #202124;
                font-size: 14px;
            }
        """)
        self.setStatusBar(self.status)
        self.status.showMessage("Ready")

    # === Core Logic ===
    def complement_base(self, base: str) -> str:
        base = base.upper()
        mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return mapping.get(base, base)

    def reverse_complement(self, seq: str) -> str:
        complement_seq = [self.complement_base(base) for base in seq]
        return ''.join(complement_seq[::-1])

    def gc_content(self, seq: str) -> float:
        """Calculate GC content percentage."""
        gc_count = seq.count('G') + seq.count('C')
        return (gc_count / len(seq)) * 100 if len(seq) > 0 else 0

    def calculate_complements(self):
        dna_sequence = self.entry.text().strip().upper()

        # --- Validation ---
        if not dna_sequence:
            self.show_error("Please enter a DNA sequence before calculating.")
            self.status.showMessage("Waiting for valid input...")
            return

        if any(base not in "ATCG" for base in dna_sequence):
            self.show_error("Invalid characters detected. Please use only A, T, C, G.")
            self.status.showMessage("Error: invalid DNA sequence entered.")
            return

        # --- Computation ---
        complement_sequence = ''.join([self.complement_base(base) for base in dna_sequence])
        reverse_complement_sequence = self.reverse_complement(dna_sequence)
        gc_percent = self.gc_content(dna_sequence)

        # --- Display Results ---
        result_text = (
            f"<b>Input Sequence:</b><br>{dna_sequence}<br><br>"
            f"<b>Complement:</b><br>{complement_sequence}<br><br>"
            f"<b>Reverse Complement:</b><br>{reverse_complement_sequence}<br><br>"
            f"<b>GC Content:</b> {gc_percent:.2f}%"
        )
        self.result_text.setHtml(result_text)
        self.result_text.moveCursor(QTextCursor.MoveOperation.Start)

        self.status.showMessage(f"Calculation complete — GC Content: {gc_percent:.2f}%")

    # === Error Message Popup ===
    def show_error(self, message):
        msg = QMessageBox(self)
        msg.setIcon(QMessageBox.Icon.Warning)
        msg.setWindowTitle("Input Error")
        msg.setText(message)
        msg.setStandardButtons(QMessageBox.StandardButton.Ok)
        msg.exec()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setApplicationDisplayName("DNA Sequence Complement Tool")
    window = DNASequenceComplementApp()
    window.show()
    sys.exit(app.exec())
