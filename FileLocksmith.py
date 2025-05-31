## For Ubuntu/Debian-based systems, run: sudo apt install libxcb-cursor0 libxcb-xinerama0 libxcb-icccm4 libxcb-image0 libxcb-keysyms1 libxcb-render-util0
## sudo dnf install polkit-gnome   # Gnome
## sudo dnf install polkit-kde     # KDE

import sys
import os
import re
import getpass
import subprocess
from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QComboBox,
    QPushButton, QTreeWidget, QTreeWidgetItem, QHeaderView, QMessageBox
)
from PySide6.QtCore import Qt


class DriveProcessViewer(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Drive Process Manager")
        self.setGeometry(100, 100, 840, 600)
        self.setup_ui()
        self.populate_drives()

    def setup_ui(self):
        main_layout = QVBoxLayout(self)

        # Top layout
        top_layout = QHBoxLayout()
        top_layout.addWidget(QLabel("Select Drive:"))

        self.drive_combo = QComboBox()
        top_layout.addWidget(self.drive_combo)

        refresh_button = QPushButton("↻ Refresh Drives")
        refresh_button.clicked.connect(self.populate_drives)
        top_layout.addWidget(refresh_button)

        scan_button = QPushButton("Scan Processes")
        scan_button.clicked.connect(self.scan_processes)
        top_layout.addWidget(scan_button)

        main_layout.addLayout(top_layout)

        # Tree widget
        self.tree = QTreeWidget()
        self.tree.setColumnCount(4)
        self.tree.setHeaderLabels(["PID", "User", "Access", "Command"])
        self.tree.header().setSectionResizeMode(QHeaderView.Stretch)
        self.tree.itemSelectionChanged.connect(self.on_tree_select)
        main_layout.addWidget(self.tree)

        # Bottom layout
        bottom_layout = QHBoxLayout()
        self.status_label = QLabel("Ready")
        bottom_layout.addWidget(self.status_label)

        self.kill_btn = QPushButton("Kill Selected Process")
        self.kill_btn.setEnabled(False)
        self.kill_btn.clicked.connect(self.kill_selected_process)
        bottom_layout.addWidget(self.kill_btn)

        self.unmount_btn = QPushButton("Unmount Drive")
        self.unmount_btn.setEnabled(False)
        self.unmount_btn.clicked.connect(self.unmount_drive)
        bottom_layout.addWidget(self.unmount_btn)

        main_layout.addLayout(bottom_layout)

    def run_as_root(self, command):
        """Try pkexec first, fallback to sudo."""
        try:
            return subprocess.run(["pkexec"] + command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        except Exception:
            return subprocess.run(["sudo"] + command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

    def populate_drives(self):
        try:
            username = getpass.getuser()
            result = subprocess.run(["df", "-h", "--output=target"], stdout=subprocess.PIPE, text=True, check=True)
            drives = [
                line.strip()
                for line in result.stdout.splitlines()[1:]
                if any([
                    line.startswith("/media/"),
                    line.startswith("/mnt/"),
                    line.startswith(f"/run/media/{username}/")
                ])
            ]
            self.drive_combo.clear()
            self.drive_combo.addItems(drives)
            if drives:
                self.update_status(f"Detected {len(drives)} mounted drives")
            else:
                self.update_status("No external drives detected")
        except Exception as e:
            self.update_status(f"Drive error: {e}")

    def scan_processes(self):
        drive = self.drive_combo.currentText()
        if not drive:
            QMessageBox.warning(self, "No Drive", "Select a drive first")
            return

        self.tree.clear()
        self.update_status(f"Scanning {drive}...")

        try:
            result = self.run_as_root(["fuser", "-vam", drive])
            lines = result.stdout.splitlines()
            if len(lines) < 2:
                self.update_status("No processes using this drive")
                return

            for line in lines[2:]:
                parts = re.split(r'\s{2,}', line.strip())
                if len(parts) >= 4:
                    pid = parts[1].split()[0]
                    user = parts[0]
                    access = parts[2]
                    command = parts[3]
                    item = QTreeWidgetItem([pid, user, access, command])
                    self.tree.addTopLevelItem(item)

            self.update_status(f"{self.tree.topLevelItemCount()} process(es) found")
            self.unmount_btn.setEnabled(True)

        except subprocess.CalledProcessError as e:
            self.update_status(f"Error: {e.stderr.strip()}")

    def on_tree_select(self):
        self.kill_btn.setEnabled(len(self.tree.selectedItems()) > 0)

    def kill_selected_process(self):
        selected = self.tree.selectedItems()
        if not selected:
            return
        pid = selected[0].text(0)
        reply = QMessageBox.question(self, "Confirm Kill", f"Kill process {pid}?", QMessageBox.Yes | QMessageBox.No)
        if reply == QMessageBox.Yes:
            try:
                self.run_as_root(["kill", "-9", pid])
                index = self.tree.indexOfTopLevelItem(selected[0])
                self.tree.takeTopLevelItem(index)
                self.update_status(f"Process {pid} terminated")
                self.kill_btn.setEnabled(False)
            except subprocess.CalledProcessError:
                self.update_status(f"Failed to kill process {pid}")

    def unmount_drive(self):
        drive = self.drive_combo.currentText()
        if not drive:
            return
        try:
            self.run_as_root(["umount", drive])
            self.update_status(f"Unmounted {drive}")
            self.populate_drives()
            self.unmount_btn.setEnabled(False)
        except subprocess.CalledProcessError as e:
            self.update_status(f"Unmount failed: {e.stderr.strip()}")

    def update_status(self, message):
        self.status_label.setText(message)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    viewer = DriveProcessViewer()
    viewer.show()
    sys.exit(app.exec())
