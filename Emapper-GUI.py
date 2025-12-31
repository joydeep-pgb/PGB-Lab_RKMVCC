import sys
import os
from collections import deque

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout,
    QHBoxLayout, QLabel, QLineEdit, QPushButton,
    QFileDialog, QSpinBox, QTextEdit, QGroupBox,
    QTabWidget, QProgressBar, QComboBox, QDoubleSpinBox,
    QStyleFactory, QMessageBox
)
from PySide6.QtCore import (
    QProcess, QSettings, QTimer, QElapsedTimer
)
from PySide6.QtGui import QFont, QTextCursor


class EggNogPro(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("eggNOG-mapper Professional")
        self.resize(900, 700)

        self.settings = QSettings("BioTools", "EggNogMapperPro")

        # Process & queue
        self.process = QProcess(self)
        self.job_queue = deque()
        self.is_running = False

        # Timers
        self.elapsed_timer = QElapsedTimer()
        self.ui_timer = QTimer()
        self.ui_timer.setInterval(1000)
        self.ui_timer.timeout.connect(self.update_runtime)

        self.setup_style()
        self.init_ui()
        self.load_settings()
        self.setup_process_signals()

    # ================= STYLE =================
    def setup_style(self):
        QApplication.setStyle(QStyleFactory.create("Fusion"))
        self.setStyleSheet("""
            QMainWindow { background-color: #2b2b2b; color: #ffffff; }
            QGroupBox { font-weight: bold; border: 1px solid #555;
                        margin-top: 10px; padding-top: 15px;
                        border-radius: 5px; color: #ddd; }
            QLabel { color: #e0e0e0; }
            QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {
                background-color: #3b3b3b; color: #fff;
                border: 1px solid #555; padding: 5px;
                border-radius: 3px;
            }
            QPushButton {
                background-color: #3c3c3c; color: white;
                border: 1px solid #555; padding: 6px 12px;
                border-radius: 4px;
            }
            QPushButton:hover { background-color: #4c4c4c; }
            QPushButton:pressed { background-color: #2c2c2c; }
            QProgressBar { border: 1px solid #555; border-radius: 3px; }
            QProgressBar::chunk { background-color: #2b78e4; }
        """)

    # ================= UI =================
    def init_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QVBoxLayout(central)

        self.tabs = QTabWidget()
        main_layout.addWidget(self.tabs)

        self.tab_run = QWidget()
        self.setup_run_tab()
        self.tabs.addTab(self.tab_run, "Run Configuration")

        self.tab_advanced = QWidget()
        self.setup_advanced_tab()
        self.tabs.addTab(self.tab_advanced, "Advanced Options")

        # Log
        log_group = QGroupBox("Execution Log")
        log_layout = QVBoxLayout(log_group)
        self.log_window = QTextEdit(readOnly=True)
        self.log_window.setFont(QFont("Consolas", 10))
        self.log_window.setStyleSheet("background:#1e1e1e; color:#d4d4d4;")
        log_layout.addWidget(self.log_window)
        main_layout.addWidget(log_group, stretch=1)

        # Control bar
        bar = QHBoxLayout()
        self.status_label = QLabel("Ready")

        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setFixedSize(150, 15)
        self.progress_bar.hide()

        self.timer_label = QLabel("00:00:00")
        self.timer_label.setFixedWidth(70)
        self.timer_label.hide()

        self.btn_queue = QPushButton("Add to Queue")
        self.btn_queue.clicked.connect(self.enqueue_job)

        self.btn_run = QPushButton("Start Queue")
        self.btn_run.setStyleSheet("background:#2b78e4; font-weight:bold;")
        self.btn_run.clicked.connect(self.start_queue)

        self.btn_stop = QPushButton("Stop")
        self.btn_stop.setStyleSheet("background:#c42b2b; font-weight:bold;")
        self.btn_stop.setEnabled(False)
        self.btn_stop.clicked.connect(self.kill_process)

        bar.addWidget(self.status_label)
        bar.addWidget(self.progress_bar)
        bar.addWidget(self.timer_label)
        bar.addStretch()
        bar.addWidget(self.btn_queue)
        bar.addWidget(self.btn_stop)
        bar.addWidget(self.btn_run)
        main_layout.addLayout(bar)

    # ================= TABS =================
    def setup_run_tab(self):
        layout = QVBoxLayout(self.tab_run)
        grp = QGroupBox("Required Inputs")
        l = QVBoxLayout(grp)

        self.input_file = self.create_browser_row(l, "Input FASTA:", True)
        self.db_dir = self.create_browser_row(l, "Database:", False)
        self.output_dir = self.create_browser_row(l, "Output Dir:", False)
        self.output_prefix = self.create_text_row(l, "Prefix:")

        layout.addWidget(grp)
        layout.addStretch()

    def setup_advanced_tab(self):
        layout = QVBoxLayout(self.tab_advanced)

        perf = QGroupBox("Performance")
        p = QHBoxLayout(perf)

        self.cpu_spin = QSpinBox()
        self.cpu_spin.setRange(1, 128)
        self.cpu_spin.setValue(8)

        self.mode_combo = QComboBox()
        self.mode_combo.addItems(["diamond", "mmseqs"])

        p.addWidget(QLabel("CPU:"))
        p.addWidget(self.cpu_spin)
        p.addWidget(QLabel("Mode:"))
        p.addWidget(self.mode_combo)
        p.addStretch()

        filt = QGroupBox("Filtering")
        f = QHBoxLayout(filt)

        # ✅ FIXED E-VALUE (precision restored)
        self.evalue_spin = QDoubleSpinBox()
        self.evalue_spin.setDecimals(5)
        self.evalue_spin.setRange(0.0, 10.0)
        self.evalue_spin.setSingleStep(0.0001)
        self.evalue_spin.setValue(0.001)
        self.evalue_spin.setToolTip("Default: 0.001")

        self.score_spin = QDoubleSpinBox()
        self.score_spin.setRange(0, 99999)
        self.score_spin.setValue(0)

        f.addWidget(QLabel("E-value:"))
        f.addWidget(self.evalue_spin)
        f.addWidget(QLabel("Score:"))
        f.addWidget(self.score_spin)
        f.addStretch()

        layout.addWidget(perf)
        layout.addWidget(filt)
        layout.addStretch()

    # ================= HELPERS =================
    def create_browser_row(self, layout, label, is_file):
        row = QHBoxLayout()
        row.addWidget(QLabel(label))
        le = QLineEdit()
        row.addWidget(le)
        btn = QPushButton("Browse")
        btn.clicked.connect(
            lambda: self.browse_file(le) if is_file else self.browse_dir(le)
        )
        row.addWidget(btn)
        layout.addLayout(row)
        return le

    def create_text_row(self, layout, label):
        row = QHBoxLayout()
        row.addWidget(QLabel(label))
        le = QLineEdit()
        row.addWidget(le)
        layout.addLayout(row)
        return le

    # ================= QUEUE =================
    def enqueue_job(self):
        if not self.input_file.text() or not self.db_dir.text():
            QMessageBox.warning(self, "Missing Input", "Input FASTA and DB required.")
            return

        job = {
            "input": self.input_file.text(),
            "db": self.db_dir.text(),
            "outdir": self.output_dir.text() or os.getcwd(),
            "prefix": self.output_prefix.text() or "results",
            "cpu": self.cpu_spin.value(),
            "mode": self.mode_combo.currentText(),
            "evalue": self.evalue_spin.value(),
            "score": self.score_spin.value()
        }

        self.job_queue.append(job)
        self.log_window.append(
            f"<b>✔ Job queued</b> ({len(self.job_queue)} in queue)<br>"
        )

    def start_queue(self):
        if self.is_running:
            return
        if not self.job_queue:
            QMessageBox.information(self, "Queue Empty", "No jobs in queue.")
            return
        self.run_next_job()

    def run_next_job(self):
        if not self.job_queue:
            self.status_label.setText("All jobs completed")
            return

        job = self.job_queue.popleft()
        self.is_running = True

        self.progress_bar.show()
        self.timer_label.show()
        self.elapsed_timer.start()
        self.ui_timer.start()

        self.btn_run.setEnabled(False)
        self.btn_stop.setEnabled(True)

        out_path = os.path.join(job["outdir"], job["prefix"])

        args = [
            "-i", job["input"],
            "--data_dir", job["db"],
            "-o", out_path,
            "-m", job["mode"],
            "--cpu", str(job["cpu"]),
            "--evalue", f"{job['evalue']:.5f}",
            "--score", str(job["score"]),
            "--override"
        ]

        self.log_window.append(
            f"<hr><b>▶ Running job</b><br>"
            f"E-value used: {job['evalue']}<br>"
            f"emapper.py {' '.join(args)}<br>"
        )

        self.process.start("emapper.py", args)

    # ================= PROCESS =================
    def setup_process_signals(self):
        self.process.readyReadStandardOutput.connect(self.handle_stdout)
        self.process.readyReadStandardError.connect(self.handle_stderr)
        self.process.finished.connect(self.process_finished)

    def handle_stdout(self):
        self.log_window.moveCursor(QTextCursor.End)
        self.log_window.insertPlainText(
            self.process.readAllStandardOutput().data().decode()
        )

    def handle_stderr(self):
        self.log_window.moveCursor(QTextCursor.End)
        self.log_window.insertHtml(
            f"<span style='color:#ff9d00;'>"
            f"{self.process.readAllStandardError().data().decode()}</span>"
        )

    def process_finished(self, code):
        self.ui_timer.stop()
        self.timer_label.hide()
        self.progress_bar.hide()

        self.is_running = False
        self.btn_run.setEnabled(True)
        self.btn_stop.setEnabled(False)

        status = "✔ Success" if code == 0 else f"✖ Failed ({code})"
        self.log_window.append(f"<b>{status}</b><br>")

        self.run_next_job()

    def kill_process(self):
        if self.process.state() == QProcess.Running:
            self.process.kill()
            self.job_queue.clear()
            self.log_window.append("<b>⛔ Job terminated by user</b><br>")

    # ================= TIMER =================
    def update_runtime(self):
        t = self.elapsed_timer.elapsed() // 1000
        h, m, s = t // 3600, (t % 3600) // 60, t % 60
        self.timer_label.setText(f"{h:02d}:{m:02d}:{s:02d}")

    # ================= IO =================
    def browse_file(self, le):
        f, _ = QFileDialog.getOpenFileName(self, "Select FASTA")
        if f:
            le.setText(f)

    def browse_dir(self, le):
        d = QFileDialog.getExistingDirectory(self, "Select Directory")
        if d:
            le.setText(d)

    # ================= SETTINGS =================
    def load_settings(self):
        self.input_file.setText(self.settings.value("input", ""))
        self.db_dir.setText(self.settings.value("db", ""))
        self.output_dir.setText(self.settings.value("outdir", ""))
        self.output_prefix.setText(self.settings.value("prefix", "results"))

    def closeEvent(self, e):
        self.settings.setValue("input", self.input_file.text())
        self.settings.setValue("db", self.db_dir.text())
        self.settings.setValue("outdir", self.output_dir.text())
        self.settings.setValue("prefix", self.output_prefix.text())
        e.accept()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    gui = EggNogPro()
    gui.show()
    sys.exit(app.exec())
