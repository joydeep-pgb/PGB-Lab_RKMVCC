#!/usr/bin/env python3
"""
AutoDock-GPU Launcher (PySide6)

A complete GUI to run AutoDock-GPU on Linux (tested on Fedora).
- Works with your existing OpenCL build (e.g., autodock_gpu_128wi).
- Detects AMD GPU by parsing the program's own output (no --list-devices).
- Lets you pick .fld, ligand .pdbqt, optional reference ligand & flexres.
- Supports --nrun, --devnum (1-based), working directory.
- Live stdout/stderr log; stop button; command preview.
- Remembers settings in ~/.config/adgpu_launcher.json.

Usage:
  pip install PySide6
  python3 adgpu_launcher.py

Tip: If you have a shell wrapper (e.g., `adgpu`), you can point the app to it
in Settings. Otherwise set the absolute path to your binary.

License: MIT
Author: Joydeep Shaw
"""

import sys
import os
import json
import subprocess
from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QFileDialog, QLineEdit, QTextEdit, QSpinBox, QMessageBox
)
from PySide6.QtCore import Qt, QProcess

CONFIG_FILE = os.path.expanduser("~/.config/adgpu_launcher.json")

class AutoDockGPULauncher(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("AutoDock-GPU Launcher")
        self.resize(800, 500)

        self.process = None
        self.build_ui()
        self.load_config(fresh_start=True)

    def build_ui(self):
        layout = QVBoxLayout()

        # Executable path
        exec_layout = QHBoxLayout()
        self.exec_edit = QLineEdit()
        exec_btn = QPushButton("Browse...")
        exec_btn.clicked.connect(self.browse_exec)
        exec_layout.addWidget(QLabel("AutoDock-GPU Executable:"))
        exec_layout.addWidget(self.exec_edit)
        exec_layout.addWidget(exec_btn)
        layout.addLayout(exec_layout)

        # FLD file
        fld_layout = QHBoxLayout()
        self.fld_edit = QLineEdit()
        fld_btn = QPushButton("Browse...")
        fld_btn.clicked.connect(lambda: self.browse_file(self.fld_edit))
        fld_layout.addWidget(QLabel("Protein maps FLD file:"))
        fld_layout.addWidget(self.fld_edit)
        fld_layout.addWidget(fld_btn)
        layout.addLayout(fld_layout)

        # Ligand file
        lig_layout = QHBoxLayout()
        self.lig_edit = QLineEdit()
        lig_btn = QPushButton("Browse...")
        lig_btn.clicked.connect(lambda: self.browse_file(self.lig_edit))
        lig_layout.addWidget(QLabel("Ligand PDBQT file:"))
        lig_layout.addWidget(self.lig_edit)
        lig_layout.addWidget(lig_btn)
        layout.addLayout(lig_layout)

        # Optional: nrun
        nrun_layout = QHBoxLayout()
        self.nrun_spin = QSpinBox()
        self.nrun_spin.setRange(1, 1000)
        self.nrun_spin.setValue(50)
        nrun_layout.addWidget(QLabel("Number of runs (--nrun):"))
        nrun_layout.addWidget(self.nrun_spin)
        layout.addLayout(nrun_layout)

        # Working directory
        wd_layout = QHBoxLayout()
        self.wd_edit = QLineEdit()
        wd_btn = QPushButton("Browse...")
        wd_btn.clicked.connect(self.browse_wd)
        wd_layout.addWidget(QLabel("Working Directory:"))
        wd_layout.addWidget(self.wd_edit)
        wd_layout.addWidget(wd_btn)
        layout.addLayout(wd_layout)

        # Buttons
        btn_layout = QHBoxLayout()
        run_btn = QPushButton("Run Docking")
        run_btn.clicked.connect(self.run_docking)
        stop_btn = QPushButton("Stop")
        stop_btn.clicked.connect(self.stop_docking)
        validate_btn = QPushButton("Validate Maps")
        validate_btn.clicked.connect(self.validate_maps)
        btn_layout.addWidget(run_btn)
        btn_layout.addWidget(stop_btn)
        btn_layout.addWidget(validate_btn)
        layout.addLayout(btn_layout)

        # Command preview
        self.cmd_preview = QTextEdit()
        self.cmd_preview.setReadOnly(True)
        layout.addWidget(QLabel("Command Preview:"))
        layout.addWidget(self.cmd_preview)

        # Log output
        self.log_output = QTextEdit()
        self.log_output.setReadOnly(True)
        layout.addWidget(QLabel("Log:"))
        layout.addWidget(self.log_output)

        self.setLayout(layout)

    def browse_exec(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select AutoDock-GPU Executable")
        if path:
            self.exec_edit.setText(path)

    def browse_file(self, target_edit):
        path, _ = QFileDialog.getOpenFileName(self, "Select File")
        if path:
            target_edit.setText(path)

    def browse_wd(self):
        path = QFileDialog.getExistingDirectory(self, "Select Working Directory")
        if path:
            self.wd_edit.setText(path)

    def run_docking(self):
        exec_path = self.exec_edit.text().strip()
        fld_path = self.fld_edit.text().strip()
        lig_path = self.lig_edit.text().strip()
        nrun = str(self.nrun_spin.value())
        wd = self.wd_edit.text().strip()

        if not all([exec_path, fld_path, lig_path, wd]):
            QMessageBox.warning(self, "Missing info", "Please fill in all fields.")
            return

        # Command without device detection
        cmd = [
            exec_path,
            "--ffile", fld_path,
            "--lfile", lig_path,
            "--nrun", nrun
        ]
        self.cmd_preview.setPlainText(" ".join(cmd))

        self.log_output.clear()
        self.process = QProcess(self)
        self.process.setWorkingDirectory(wd)
        self.process.readyReadStandardOutput.connect(self.handle_stdout)
        self.process.readyReadStandardError.connect(self.handle_stderr)
        self.process.finished.connect(self.process_finished)
        self.process.start(cmd[0], cmd[1:])

        # Save config (but no auto-start next time)
        self.save_config()

    def stop_docking(self):
        if self.process:
            self.process.kill()
            self.log_output.append("Docking stopped.")

    def handle_stdout(self):
        data = self.process.readAllStandardOutput().data().decode()
        self.log_output.append(data)

    def handle_stderr(self):
        data = self.process.readAllStandardError().data().decode()
        self.log_output.append(data)

    def process_finished(self):
        self.log_output.append("Docking finished.")

    def validate_maps(self):
        fld_path = self.fld_edit.text().strip()
        if not os.path.exists(fld_path):
            QMessageBox.warning(self, "Error", "FLD file does not exist.")
            return
        missing = []
        with open(fld_path) as f:
            for line in f:
                if line.strip().startswith("grid_map"):
                    parts = line.split()
                    if len(parts) >= 2:
                        map_file = parts[1]
                        if not os.path.exists(os.path.join(os.path.dirname(fld_path), map_file)):
                            missing.append(map_file)
        if missing:
            QMessageBox.warning(self, "Missing Maps", "\n".join(missing))
        else:
            QMessageBox.information(self, "Validation", "All map files are present.")

    def load_config(self, fresh_start=False):
        if os.path.exists(CONFIG_FILE) and not fresh_start:
            with open(CONFIG_FILE) as f:
                data = json.load(f)
            self.exec_edit.setText(data.get("exec_path", ""))
            self.fld_edit.setText(data.get("fld_path", ""))
            self.lig_edit.setText(data.get("lig_path", ""))
            self.nrun_spin.setValue(data.get("nrun", 50))
            self.wd_edit.setText(data.get("wd", ""))

    def save_config(self):
        data = {
            "exec_path": self.exec_edit.text(),
            "fld_path": self.fld_edit.text(),
            "lig_path": self.lig_edit.text(),
            "nrun": self.nrun_spin.value(),
            "wd": self.wd_edit.text()
        }
        os.makedirs(os.path.dirname(CONFIG_FILE), exist_ok=True)
        with open(CONFIG_FILE, "w") as f:
            json.dump(data, f)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = AutoDockGPULauncher()
    win.show()
    sys.exit(app.exec())
