import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QLabel, QLineEdit, QTextEdit, QPushButton, QMenu
from PyQt5.QtGui import QContextMenuEvent
from PyQt5.QtCore import Qt

class DNASequenceComplementApp(QWidget):
    def __init__(self):
        super().__init__()
        
        self.setWindowTitle("DNA Sequence Complement")
        self.setGeometry(300, 300, 400, 300)
        
        layout = QVBoxLayout()

        # Label and input field
        label = QLabel("Enter DNA sequence:")
        layout.addWidget(label)
        
        self.entry = QLineEdit()
        layout.addWidget(self.entry)

        # Calculate button
        calculate_button = QPushButton("Calculate")
        calculate_button.clicked.connect(self.calculate_complements)
        layout.addWidget(calculate_button)

        # Result text field
        self.result_text = QTextEdit()
        self.result_text.setReadOnly(True)
        layout.addWidget(self.result_text)
        
        self.setLayout(layout)

    def contextMenuEvent(self, event):
        menu = QMenu(self)
        
        copy_action = menu.addAction("Copy")
        paste_action = menu.addAction("Paste")
        
        action = menu.exec_(self.mapToGlobal(event.pos()))
        
        if action == copy_action:
            if self.focusWidget() == self.entry:
                self.entry.copy()
            elif self.focusWidget() == self.result_text:
                self.result_text.copy()
        elif action == paste_action:
            if self.focusWidget() == self.entry:
                self.entry.paste()

    def complement_base(self, base):
        complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return complements.get(base, base)

    def reverse_complement(self, seq):
        complement_seq = [self.complement_base(base) for base in seq]
        reverse_complement_seq = ''.join(complement_seq[::-1])
        return reverse_complement_seq

    def calculate_complements(self):
        dna_sequence = self.entry.text().upper()

        complement_sequence = ''.join([self.complement_base(base) for base in dna_sequence])
        reverse_complement_sequence = self.reverse_complement(dna_sequence)

        result_text = f"Complement: {complement_sequence}\nReverse complement: {reverse_complement_sequence}"
        self.result_text.setPlainText(result_text)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = DNASequenceComplementApp()
    window.show()
    sys.exit(app.exec_())
