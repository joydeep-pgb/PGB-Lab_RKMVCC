#!/usr/bin/env python3
"""
PlantTF-BLAST GUI
==================
A PyQt6 desktop application for running custom BLAST searches against a
family-annotated transcription factor (TF) protein database (e.g. a
PlantTFDB species download) and reporting PlantTFDB-style predictions:

    Query_ID    TF_Family    Best_Hit_ID    E-value    Description

Requirements:
    - Python 3.9+
    - PyQt6           (pip install PyQt6)
    - NCBI BLAST+     (makeblastdb, blastp, blastx on PATH)
    - openpyxl        (optional; only needed to export .xlsx)

Usage:
    python main.py
"""

import sys

from PyQt6.QtWidgets import QApplication

from app.main_window import MainWindow


def main():
    app = QApplication(sys.argv)
    app.setApplicationName("PlantTF-BLAST")
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
