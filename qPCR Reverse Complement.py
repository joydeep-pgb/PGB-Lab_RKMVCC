import os
os.environ.setdefault("QT_QPA_PLATFORM", "xcb")

import sys
from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QTextEdit, QPushButton, QMenu, QFrame
)
from PyQt5.QtGui import QFont
from PyQt5.QtCore import Qt

# Nucleotide colors, loosely modeled on standard genome-browser base coloring
BASE_COLORS = {
    "A": "#4ADE80",  # green
    "T": "#F87171",  # red
    "U": "#F87171",  # red
    "C": "#60A5FA",  # blue
    "G": "#FBBF24",  # amber
}
AMBIGUOUS_COLOR = "#C084FC"   # IUPAC ambiguity codes (R, Y, N, etc.)
INVALID_COLOR = "#7D8590"     # anything else (whitespace, stray characters)

MONO_FONT_STACK = '"Cascadia Code", "Fira Code", "JetBrains Mono", "DejaVu Sans Mono", monospace'
SANS_FONT_STACK = '"Inter", "Segoe UI", "Noto Sans", sans-serif'

STYLE_SHEET = f"""
QWidget {{
    background-color: #1B1F23;
    color: #E6EDF3;
    font-family: {SANS_FONT_STACK};
    font-size: 13px;
}}

QLabel#eyebrow {{
    color: #7D8590;
    font-size: 11px;
    font-weight: 600;
    letter-spacing: 2px;
}}

QLabel#stats {{
    color: #7D8590;
    font-size: 11px;
}}

QLineEdit {{
    background-color: #22272E;
    border: 1px solid #333B44;
    border-radius: 8px;
    padding: 10px 12px;
    font-family: {MONO_FONT_STACK};
    font-size: 14px;
    color: #E6EDF3;
    selection-background-color: #2DD4BF;
    selection-color: #0B1215;
}}
QLineEdit:focus {{
    border: 1px solid #2DD4BF;
}}

QPushButton {{
    background-color: #2DD4BF;
    color: #0B1215;
    border: none;
    border-radius: 8px;
    padding: 10px 18px;
    font-weight: 600;
    font-size: 13px;
}}
QPushButton:hover {{
    background-color: #5EEAD4;
}}
QPushButton:pressed {{
    background-color: #14B8A6;
}}

QTextEdit {{
    background-color: #22272E;
    border: 1px solid #333B44;
    border-radius: 8px;
    padding: 12px;
    font-family: {MONO_FONT_STACK};
    font-size: 14px;
    selection-background-color: #2DD4BF;
    selection-color: #0B1215;
}}

QFrame#divider {{
    background-color: #333B44;
    max-height: 1px;
    min-height: 1px;
}}
"""


class DNASequenceComplementApp(QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Reverse Complement")
        self.setGeometry(300, 300, 460, 440)
        self.setStyleSheet(STYLE_SHEET)

        layout = QVBoxLayout()
        layout.setContentsMargins(24, 24, 24, 24)
        layout.setSpacing(10)

        title = QLabel("qPCR · Reverse Complement")
        title_font = QFont()
        title_font.setPointSize(15)
        title_font.setWeight(QFont.DemiBold)
        title.setFont(title_font)
        layout.addWidget(title)
        layout.addSpacing(4)

        input_label = QLabel("SEQUENCE")
        input_label.setObjectName("eyebrow")
        layout.addWidget(input_label)

        self.entry = QLineEdit()
        self.entry.setPlaceholderText("e.g. ATCGGCTAAGCT")
        self.entry.returnPressed.connect(self.calculate_complements)
        layout.addWidget(self.entry)

        self.stats_label = QLabel("Length: 0 nt   ·   GC content: —")
        self.stats_label.setObjectName("stats")
        layout.addWidget(self.stats_label)

        button_row = QHBoxLayout()
        button_row.addStretch()
        calculate_button = QPushButton("Calculate")
        calculate_button.clicked.connect(self.calculate_complements)
        button_row.addWidget(calculate_button)
        layout.addLayout(button_row)

        layout.addSpacing(6)
        divider = QFrame()
        divider.setObjectName("divider")
        divider.setFrameShape(QFrame.HLine)
        layout.addWidget(divider)
        layout.addSpacing(6)

        output_label = QLabel("RESULT")
        output_label.setObjectName("eyebrow")
        layout.addWidget(output_label)

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

    def colorize(self, seq):
        """Render a sequence as HTML with per-base nucleotide coloring."""
        spans = []
        for base in seq:
            if base in BASE_COLORS:
                color = BASE_COLORS[base]
            elif base.isalpha():
                color = AMBIGUOUS_COLOR
            else:
                color = INVALID_COLOR
            spans.append(f'<span style="color:{color};">{base}</span>')
        return ''.join(spans)

    def calculate_complements(self):
        dna_sequence = self.entry.text().upper().strip()

        if not dna_sequence:
            self.result_text.setPlainText("")
            self.stats_label.setText("Length: 0 nt   ·   GC content: —")
            return

        complement_sequence = ''.join(self.complement_base(b) for b in dna_sequence)
        reverse_complement_sequence = self.reverse_complement(dna_sequence)

        length = len(dna_sequence)
        gc_count = sum(1 for b in dna_sequence if b in ("G", "C"))
        gc_content = (gc_count / length * 100) if length else 0
        self.stats_label.setText(f"Length: {length} nt   ·   GC content: {gc_content:.1f}%")

        html = (
            '<div style="font-family:{font}; font-size:14px; line-height:1.6;">'
            '<b style="color:#7D8590; font-weight:600;">COMPLEMENT</b><br>{comp}<br><br>'
            '<b style="color:#7D8590; font-weight:600;">REVERSE COMPLEMENT</b><br>{rc}'
            '</div>'
        ).format(
            font=MONO_FONT_STACK,
            comp=self.colorize(complement_sequence),
            rc=self.colorize(reverse_complement_sequence),
        )
        self.result_text.setHtml(html)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = DNASequenceComplementApp()
    window.show()
    sys.exit(app.exec_())