import sys
import os
import datetime
from pathlib import Path
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QTabWidget,
    QGroupBox, QLabel, QComboBox, QTextEdit, QLineEdit, QPushButton,
    QFileDialog, QMessageBox, QSpinBox, QCheckBox, QTableWidget,
    QTableWidgetItem, QFormLayout, QDoubleSpinBox, QScrollArea, QGridLayout
)
from PyQt6.QtCore import QThread, pyqtSignal
from Bio.Blast.Applications import (
    NcbiblastnCommandline, NcbiblastpCommandline, NcbiblastxCommandline,
    NcbitblastnCommandline, NcbitblastxCommandline, NcbimakeblastdbCommandline
)
import multiprocessing
import warnings
from Bio import BiopythonDeprecationWarning
warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)


# ===========================================================
#                    WORKER CLASSES
# ===========================================================
class BlastWorker(QThread):
    progress_signal = pyqtSignal(str)
    result_signal = pyqtSignal(object)
    error_signal = pyqtSignal(str)
    finished_signal = pyqtSignal()
    command_signal = pyqtSignal(str)

    def __init__(self, program, query, database, parameters, output_dir):
        super().__init__()
        self.program = program
        self.query = query
        self.database = database
        self.parameters = parameters
        self.output_dir = output_dir

    def run(self):
        try:
            self.run_local_blast()
        except Exception as e:
            self.error_signal.emit(str(e))
        finally:
            self.finished_signal.emit()

    def run_local_blast(self):
        self.progress_signal.emit("Preparing local BLAST...")

        # === Timestamped Output Folder ===
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        out_dir = os.path.join(self.output_dir, f"{self.program}_{timestamp}")
        os.makedirs(out_dir, exist_ok=True)

        query_path = os.path.join(out_dir, "query.fasta")
        output_path = os.path.join(out_dir, f"{self.program}_results.txt")

        with open(query_path, "w") as qf:
            qf.write(f">query\n{self.query}")

        # === Select proper BLAST command ===
        blast_class = {
            'blastn': NcbiblastnCommandline,
            'blastp': NcbiblastpCommandline,
            'blastx': NcbiblastxCommandline,
            'tblastn': NcbitblastnCommandline,
            'tblastx': NcbitblastxCommandline
        }.get(self.program, NcbiblastpCommandline)

        # === Minimal command arguments ===
        cmd_args = {
            'query': query_path,
            'db': self.database,
            'out': output_path,
            'outfmt': 6,
            'evalue': self.parameters['evalue'],
            'max_target_seqs': self.parameters['max_target_seqs'],
            'num_threads': self.parameters['num_threads']
        }

        # Add selected task
        if self.parameters.get('task'):
            cmd_args['task'] = self.parameters['task']

        # Handle filtering
        if not self.parameters['filter']:
            cmd_args['dust' if self.program == 'blastn' else 'seg'] = 'no'

        cmdline = blast_class(**cmd_args)
        cmd_str = str(cmdline)
        self.command_signal.emit(f"Executing:\n{cmd_str}")
        print(f"\n[BLAST COMMAND] {cmd_str}\n")

        # === Run BLAST ===
        self.progress_signal.emit("Running BLAST locally...")
        stdout, stderr = cmdline()

        self.progress_signal.emit("Parsing results...")
        with open(output_path, 'r') as f:
            results = self.parse_tabular_output(f.read())

        self.progress_signal.emit(f"Results saved to: {output_path}")
        self.command_signal.emit(f"Output file: {output_path}")
        self.result_signal.emit(results)

    def parse_tabular_output(self, data):
        results = []
        for line in data.strip().split('\n'):
            if not line or line.startswith('#'):
                continue
            f = line.split('\t')
            if len(f) >= 12:
                results.append({
                    'qseqid': f[0], 'sseqid': f[1], 'pident': float(f[2]),
                    'length': int(f[3]), 'mismatch': int(f[4]), 'gapopen': int(f[5]),
                    'qstart': int(f[6]), 'qend': int(f[7]), 'sstart': int(f[8]), 'send': int(f[9]),
                    'evalue': float(f[10]), 'bitscore': float(f[11])
                })
        return results


class DatabaseBuilder(QThread):
    progress_signal = pyqtSignal(str)
    finished_signal = pyqtSignal(bool, str)
    command_signal = pyqtSignal(str)

    def __init__(self, fasta_file, db_type, db_name, title, output_dir):
        super().__init__()
        self.fasta_file = fasta_file
        self.db_type = db_type
        self.db_name = db_name
        self.title = title
        self.output_dir = output_dir

    def run(self):
        try:
            if not os.path.exists(self.fasta_file):
                self.finished_signal.emit(False, "Input FASTA not found")
                return

            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            out_dir = os.path.join(self.output_dir, f"makeblastdb_{timestamp}")
            os.makedirs(out_dir, exist_ok=True)
            db_path = os.path.join(out_dir, self.db_name)

            self.progress_signal.emit("Building BLAST database...")
            cmd = NcbimakeblastdbCommandline(
                input_file=self.fasta_file,
                dbtype=self.db_type,
                out=db_path,
                title=self.title,
                parse_seqids=True,
                hash_index=True
            )
            cmd_str = str(cmd)
            self.command_signal.emit(f"Executing:\n{cmd_str}")
            print(f"\n[MAKEBLASTDB COMMAND] {cmd_str}\n")

            stdout, stderr = cmd()

            exts = ['.nhr', '.nin', '.nsq'] if self.db_type == 'nucl' else ['.phr', '.pin', '.psq']
            if all(os.path.exists(db_path + e) for e in exts):
                self.finished_signal.emit(True, f"Database '{self.db_name}' created at:\n{out_dir}")
            else:
                self.finished_signal.emit(False, f"Database build failed.\n{stderr}")
        except Exception as e:
            self.finished_signal.emit(False, f"Error: {e}")


# ===========================================================
#                    MAIN GUI
# ===========================================================
class BlastApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle("PyBLAST Pro")
        self.setGeometry(100, 100, 1050, 800)

        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)

        tabs = QTabWidget()
        layout.addWidget(tabs)
        tabs.addTab(self.create_blast_tab(), "BLAST Search")
        tabs.addTab(self.create_database_tab(), "Database Builder")
        tabs.addTab(self.create_results_tab(), "Results")

        self.status_bar = self.statusBar()

    # ------------------- BLAST TAB -------------------
    def create_blast_tab(self):
        widget = QWidget()
        layout = QVBoxLayout(widget)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        content = QWidget()
        vbox = QVBoxLayout(content)

        config = QGroupBox("BLAST Configuration")
        form = QFormLayout(config)

        # Program selection
        self.program_combo = QComboBox()
        self.program_combo.addItems(["blastn", "blastp", "blastx", "tblastn", "tblastx"])
        self.program_combo.currentTextChanged.connect(self.update_task_options)
        form.addRow("Program:", self.program_combo)

        # Task selection
        self.task_combo = QComboBox()
        form.addRow("Task:", self.task_combo)
        self.update_task_options("blastn")  # default load

        # Database
        db_layout = QHBoxLayout()
        self.db_combo = QComboBox()
        self.db_combo.setEditable(True)
        db_browse = QPushButton("Browse Folder...")
        db_browse.clicked.connect(self.browse_database)
        db_layout.addWidget(self.db_combo)
        db_layout.addWidget(db_browse)
        form.addRow("Database Folder:", db_layout)

        # Output directory
        out_layout = QHBoxLayout()
        self.output_dir = QLineEdit()
        self.output_dir.setPlaceholderText("Select output folder for results...")
        browse_out = QPushButton("Browse...")
        browse_out.clicked.connect(self.browse_output_dir)
        out_layout.addWidget(self.output_dir)
        out_layout.addWidget(browse_out)
        form.addRow("Output Folder:", out_layout)
        vbox.addWidget(config)

        # Parameters
        params = QGroupBox("Parameters")
        grid = QGridLayout(params)

        grid.addWidget(QLabel("E-value:"), 0, 0)
        self.evalue_spin = QDoubleSpinBox()
        self.evalue_spin.setDecimals(10)
        self.evalue_spin.setRange(1e-200, 1e+10)
        self.evalue_spin.setValue(1e-5)
        grid.addWidget(self.evalue_spin, 0, 1)

        grid.addWidget(QLabel("Max Target Seqs:"), 0, 2)
        self.max_target = QSpinBox()
        self.max_target.setValue(500)
        grid.addWidget(self.max_target, 0, 3)

        grid.addWidget(QLabel("Threads:"), 1, 0)
        self.threads = QSpinBox()
        self.threads.setRange(1, multiprocessing.cpu_count())
        self.threads.setValue(1)
        grid.addWidget(self.threads, 1, 1)

        self.filter_check = QCheckBox("Low Complexity Filter (SEG/DUST)")
        self.filter_check.setChecked(True)
        grid.addWidget(self.filter_check, 2, 0, 1, 4)
        vbox.addWidget(params)

        # Query input
        query = QGroupBox("Query Sequence")
        qv = QVBoxLayout(query)
        self.query_input = QTextEdit()
        self.query_input.setPlaceholderText("Paste FASTA sequence or load from file...")
        qv.addWidget(self.query_input)
        load_query_btn = QPushButton("Load FASTA File")
        load_query_btn.clicked.connect(self.load_fasta_query)
        qv.addWidget(load_query_btn)
        vbox.addWidget(query)

        scroll.setWidget(content)
        layout.addWidget(scroll)

        run = QPushButton("Run Local BLAST")
        run.clicked.connect(self.run_blast_search)
        layout.addWidget(run)
        return widget

    # ------------------- DATABASE TAB -------------------
    def create_database_tab(self):
        widget = QWidget()
        layout = QVBoxLayout(widget)
        group = QGroupBox("Create Local Database")
        form = QFormLayout(group)

        self.fasta_path = QLineEdit()
        browse = QPushButton("Browse FASTA")
        browse.clicked.connect(self.browse_fasta_file)
        h = QHBoxLayout()
        h.addWidget(self.fasta_path)
        h.addWidget(browse)
        form.addRow("FASTA File:", h)

        self.db_type_combo = QComboBox()
        self.db_type_combo.addItems(["nucl", "prot"])
        form.addRow("Type:", self.db_type_combo)

        self.db_name_input = QLineEdit()
        form.addRow("Database Name:", self.db_name_input)

        self.db_title_input = QLineEdit()
        form.addRow("Title:", self.db_title_input)

        out_layout = QHBoxLayout()
        self.db_output_dir = QLineEdit()
        self.db_output_dir.setPlaceholderText("Select output folder for database...")
        browse_db_out = QPushButton("Browse...")
        browse_db_out.clicked.connect(self.browse_db_output_dir)
        out_layout.addWidget(self.db_output_dir)
        out_layout.addWidget(browse_db_out)
        form.addRow("Output Folder:", out_layout)

        layout.addWidget(group)
        build = QPushButton("Build Database")
        build.clicked.connect(self.build_database)
        layout.addWidget(build)
        self.db_progress = QLabel("Ready.")
        layout.addWidget(self.db_progress)
        return widget

    # ------------------- RESULTS TAB -------------------
    def create_results_tab(self):
        w = QWidget()
        v = QVBoxLayout(w)
        self.hits_table = QTableWidget()
        v.addWidget(self.hits_table)
        self.log_box = QTextEdit()
        self.log_box.setReadOnly(True)
        v.addWidget(QLabel("Command Log"))
        v.addWidget(self.log_box)
        return w

    # ------------------- LOGIC -------------------
    def update_task_options(self, program):
        self.task_combo.clear()
        task_map = {
            'blastn': ['megablast', 'blastn', 'dc-megablast'],
            'blastp': ['blastp', 'blastp-fast', 'blastp-short'],
            'blastx': ['blastx', 'blastx-fast'],
            'tblastn': ['tblastn', 'tblastn-fast'],
            'tblastx': ['tblastx']
        }
        self.task_combo.addItems(task_map.get(program, ['default']))

    def load_fasta_query(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select Query FASTA File", "", "FASTA Files (*.fasta *.fa *.fna *.faa)")
        if file:
            try:
                with open(file, 'r') as f:
                    content = f.read().strip()
                if not content.startswith('>'):
                    QMessageBox.warning(self, "Invalid File", "The selected file does not appear to be a valid FASTA file.")
                    return
                self.query_input.setPlainText(content)
                QMessageBox.information(self, "FASTA Loaded", f"Loaded query sequence from:\n{file}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to load FASTA file:\n{str(e)}")

    def browse_database(self):
        folder = QFileDialog.getExistingDirectory(self, "Select Database Folder")
        if folder:
            db_files = list(Path(folder).glob("*.pin")) + list(Path(folder).glob("*.nin"))
            if db_files:
                prefix = db_files[0].with_suffix('').as_posix()
                self.db_combo.setCurrentText(prefix)
            else:
                self.db_combo.setCurrentText(folder)

    def browse_fasta_file(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select FASTA", "", "FASTA Files (*.fasta *.fa)")
        if file:
            self.fasta_path.setText(file)

    def browse_output_dir(self):
        folder = QFileDialog.getExistingDirectory(self, "Select Output Folder for BLAST Results")
        if folder:
            self.output_dir.setText(folder)

    def browse_db_output_dir(self):
        folder = QFileDialog.getExistingDirectory(self, "Select Output Folder for Database")
        if folder:
            self.db_output_dir.setText(folder)

    def run_blast_search(self):
        query = self.query_input.toPlainText().strip()
        if not query:
            QMessageBox.warning(self, "Error", "Please enter or load a query sequence!")
            return
        db = self.db_combo.currentText().strip()
        if not db:
            QMessageBox.warning(self, "Error", "Please select a database folder!")
            return
        outdir = self.output_dir.text().strip() or os.getcwd()

        params = {
            'evalue': self.evalue_spin.value(),
            'max_target_seqs': self.max_target.value(),
            'num_threads': self.threads.value(),
            'filter': self.filter_check.isChecked(),
            'task': self.task_combo.currentText()
        }

        self.worker = BlastWorker(self.program_combo.currentText(), query, db, params, outdir)
        self.worker.progress_signal.connect(self.status_bar.showMessage)
        self.worker.result_signal.connect(self.display_results)
        self.worker.error_signal.connect(self.show_error)
        self.worker.command_signal.connect(self.log_command)
        self.worker.start()

    def build_database(self):
        fasta = self.fasta_path.text().strip()
        dbname = self.db_name_input.text().strip()
        title = self.db_title_input.text().strip()
        outdir = self.db_output_dir.text().strip() or os.getcwd()
        if not fasta or not dbname:
            QMessageBox.warning(self, "Error", "Provide FASTA and database name!")
            return
        self.db_builder = DatabaseBuilder(fasta, self.db_type_combo.currentText(), dbname, title, outdir)
        self.db_builder.progress_signal.connect(self.db_progress.setText)
        self.db_builder.command_signal.connect(self.log_command)
        self.db_builder.finished_signal.connect(lambda s, m: QMessageBox.information(self, "Result", m))
        self.db_builder.start()

    def display_results(self, results):
        self.hits_table.setColumnCount(12)
        headers = ["Query", "Subject", "%Identity", "Length", "Mismatch", "GapOpen",
                   "Q.start", "Q.end", "S.start", "S.end", "Evalue", "BitScore"]
        self.hits_table.setHorizontalHeaderLabels(headers)
        self.hits_table.setRowCount(len(results))
        for i, r in enumerate(results):
            self.hits_table.setItem(i, 0, QTableWidgetItem(r['qseqid']))
            self.hits_table.setItem(i, 1, QTableWidgetItem(r['sseqid']))
            self.hits_table.setItem(i, 2, QTableWidgetItem(f"{r['pident']:.1f}"))
            self.hits_table.setItem(i, 3, QTableWidgetItem(str(r['length'])))
            self.hits_table.setItem(i, 10, QTableWidgetItem(f"{r['evalue']:.2e}"))
        self.hits_table.resizeColumnsToContents()
        self.status_bar.showMessage("Local BLAST completed.")

    def log_command(self, cmd):
        self.log_box.append(f"\n{cmd}\n")

    def show_error(self, msg):
        QMessageBox.critical(self, "Error", msg)
        self.status_bar.showMessage(msg)


# ===========================================================
#                    MAIN
# ===========================================================
def main():
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    win = BlastApp()
    win.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
