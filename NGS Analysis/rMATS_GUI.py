import os
import sys
import subprocess
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QLineEdit, QPushButton, QGroupBox, QSpinBox, QCheckBox,
    QFileDialog, QTextEdit, QComboBox, QProgressBar, QMessageBox,
    QMenuBar, QMenu, QToolBar, QStatusBar
)
from PySide6.QtCore import QProcess, Qt, QTimer, QByteArray
from PySide6.QtGui import QColor, QPalette, QTextCursor, QAction, QKeySequence, QIcon
from PySide6.QtCore import QStandardPaths

class rMATSGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("rMATS Alternative Splicing Analyzer")
        self.setGeometry(100, 100, 900, 700)

        # Create menu bar first
        self.create_menu_bar()
        
        # Create status bar
        self.statusBar = QStatusBar()
        self.setStatusBar(self.statusBar)
        self.statusBar.showMessage("Ready")

        # Central Widget and Layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)

        # Input Parameters Group
        input_group = QGroupBox("Input Parameters")
        input_layout = QVBoxLayout(input_group)

        # Batch 1 (always required)
        b1_layout = QHBoxLayout()
        self.b1_label = QLabel("Batch 1 File*:")
        b1_layout.addWidget(self.b1_label)
        self.b1_edit = QLineEdit()
        b1_layout.addWidget(self.b1_edit)
        self.b1_button = QPushButton("Browse")
        self.b1_button.clicked.connect(lambda: self.browse_file(self.b1_edit))
        b1_layout.addWidget(self.b1_button)
        input_layout.addLayout(b1_layout)

        # Batch 2 (optional when statoff enabled)
        b2_layout = QHBoxLayout()
        self.b2_label = QLabel("Batch 2 File*:")
        b2_layout.addWidget(self.b2_label)
        self.b2_edit = QLineEdit()
        b2_layout.addWidget(self.b2_edit)
        self.b2_button = QPushButton("Browse")
        self.b2_button.clicked.connect(lambda: self.browse_file(self.b2_edit))
        b2_layout.addWidget(self.b2_button)
        input_layout.addLayout(b2_layout)

        # GTF File (always required)
        gtf_layout = QHBoxLayout()
        self.gtf_label = QLabel("GTF File*:")
        gtf_layout.addWidget(self.gtf_label)
        self.gtf_edit = QLineEdit()
        gtf_layout.addWidget(self.gtf_edit)
        self.gtf_button = QPushButton("Browse")
        self.gtf_button.clicked.connect(lambda: self.browse_file(self.gtf_edit))
        gtf_layout.addWidget(self.gtf_button)
        input_layout.addLayout(gtf_layout)

        # Output Directory (always required)
        out_layout = QHBoxLayout()
        self.out_label = QLabel("Output Directory*:")
        out_layout.addWidget(self.out_label)
        self.out_edit = QLineEdit()
        out_layout.addWidget(self.out_edit)
        self.out_button = QPushButton("Browse")
        self.out_button.clicked.connect(self.browse_directory)
        out_layout.addWidget(self.out_button)
        input_layout.addLayout(out_layout)

        # Temp Directory (always required)
        temp_layout = QHBoxLayout()
        self.temp_label = QLabel("Temp Directory*:")
        temp_layout.addWidget(self.temp_label)
        self.temp_edit = QLineEdit()
        temp_layout.addWidget(self.temp_edit)
        self.temp_button = QPushButton("Browse")
        self.temp_button.clicked.connect(self.browse_temp_directory)
        temp_layout.addWidget(self.temp_button)
        input_layout.addLayout(temp_layout)

        # Options Group
        options_group = QGroupBox("Analysis Options")
        options_layout = QVBoxLayout(options_group)

        # Read Type
        type_layout = QHBoxLayout()
        type_layout.addWidget(QLabel("Read Type*:"))
        self.type_combo = QComboBox()
        self.type_combo.addItems(["paired", "single"])
        self.type_combo.setCurrentText("paired")
        type_layout.addWidget(self.type_combo)
        options_layout.addLayout(type_layout)

        # Read Length
        read_len_layout = QHBoxLayout()
        read_len_layout.addWidget(QLabel("Read Length*:"))
        self.read_len_spin = QSpinBox()
        self.read_len_spin.setRange(1, 500)
        self.read_len_spin.setValue(295)
        read_len_layout.addWidget(self.read_len_spin)
        options_layout.addLayout(read_len_layout)

        # Threads
        threads_layout = QHBoxLayout()
        threads_layout.addWidget(QLabel("Threads*:"))
        self.threads_spin = QSpinBox()
        self.threads_spin.setRange(1, 64)
        self.threads_spin.setValue(8)
        threads_layout.addWidget(self.threads_spin)
        options_layout.addLayout(threads_layout)

        # Checkboxes
        self.var_read_check = QCheckBox("Variable Read Length")
        self.var_read_check.setChecked(True)
        options_layout.addWidget(self.var_read_check)

        self.statoff_check = QCheckBox("Disable Statistical Analysis (--statoff)")
        self.statoff_check.setChecked(False)
        self.statoff_check.stateChanged.connect(self.toggle_b2_requirement)
        options_layout.addWidget(self.statoff_check)

        # Progress and Output
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        main_layout.addWidget(self.progress_bar)
        
        # Status label (keeping for backward compatibility, but status bar will be primary)
        self.status_label = QLabel("Ready")
        self.status_label.setVisible(False)  # Hide since we have status bar now
        main_layout.addWidget(self.status_label)

        self.output_text = QTextEdit()
        self.output_text.setReadOnly(True)
        main_layout.addWidget(self.output_text)

        # Run Button
        self.run_button = QPushButton("Run rMATS Analysis")
        self.run_button.clicked.connect(self.run_analysis)
        main_layout.addWidget(self.run_button)

        # Debug Button (moved to menu, keeping for compatibility)
        self.debug_button = QPushButton("Debug Information")
        self.debug_button.clicked.connect(self.show_debug_info)
        self.debug_button.setVisible(False)  # Hide since it's in menu now
        main_layout.addWidget(self.debug_button)

        # Add groups to main layout
        main_layout.addWidget(input_group)
        main_layout.addWidget(options_group)

        # QProcess for running commands
        self.process = QProcess()
        self.process.readyReadStandardOutput.connect(self.handle_stdout)
        self.process.readyReadStandardError.connect(self.handle_stderr)
        self.process.started.connect(self.process_started)
        self.process.finished.connect(self.analysis_finished)
        self.process.errorOccurred.connect(self.process_error)

        # Timer for process monitoring
        self.timer = QTimer()
        self.timer.timeout.connect(self.monitor_process)
        
        # Initial UI update
        self.toggle_b2_requirement()

    def create_menu_bar(self):
        """Create professional menu bar"""
        menubar = self.menuBar()
        
        # File Menu
        file_menu = menubar.addMenu('&File')
        
        # New Analysis
        new_action = QAction('&New Analysis', self)
        new_action.setShortcut(QKeySequence.New)
        new_action.setStatusTip('Start a new analysis (clear all fields)')
        new_action.triggered.connect(self.new_analysis)
        file_menu.addAction(new_action)
        
        file_menu.addSeparator()
        
        # Load Configuration
        load_config_action = QAction('&Load Configuration...', self)
        load_config_action.setShortcut(QKeySequence.Open)
        load_config_action.setStatusTip('Load analysis configuration from file')
        load_config_action.triggered.connect(self.load_configuration)
        file_menu.addAction(load_config_action)
        
        # Save Configuration
        save_config_action = QAction('&Save Configuration...', self)
        save_config_action.setShortcut(QKeySequence.Save)
        save_config_action.setStatusTip('Save current analysis configuration')
        save_config_action.triggered.connect(self.save_configuration)
        file_menu.addAction(save_config_action)
        
        file_menu.addSeparator()
        
        # Export Results
        export_action = QAction('&Export Results...', self)
        export_action.setStatusTip('Export analysis output to file')
        export_action.triggered.connect(self.export_results)
        file_menu.addAction(export_action)
        
        file_menu.addSeparator()
        
        # Exit
        exit_action = QAction('E&xit', self)
        exit_action.setShortcut(QKeySequence.Quit)
        exit_action.setStatusTip('Exit the application')
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)
        
        # Analysis Menu
        analysis_menu = menubar.addMenu('&Analysis')
        
        # Run Analysis
        run_action = QAction('&Run Analysis', self)
        run_action.setShortcut(QKeySequence('Ctrl+R'))
        run_action.setStatusTip('Start rMATS analysis')
        run_action.triggered.connect(self.run_analysis)
        analysis_menu.addAction(run_action)
        
        # Stop Analysis
        self.stop_action = QAction('&Stop Analysis', self)
        self.stop_action.setShortcut(QKeySequence('Ctrl+S'))
        self.stop_action.setStatusTip('Stop running analysis')
        self.stop_action.triggered.connect(self.stop_analysis)
        self.stop_action.setEnabled(False)
        analysis_menu.addAction(self.stop_action)
        
        analysis_menu.addSeparator()
        
        # Validate Inputs
        validate_action = QAction('&Validate Inputs', self)
        validate_action.setStatusTip('Check if all required inputs are valid')
        validate_action.triggered.connect(self.validate_inputs)
        analysis_menu.addAction(validate_action)
        
        # Clear Output
        clear_action = QAction('&Clear Output', self)
        clear_action.setShortcut(QKeySequence('Ctrl+L'))
        clear_action.setStatusTip('Clear the output window')
        clear_action.triggered.connect(self.clear_output)
        analysis_menu.addAction(clear_action)
        
        # Tools Menu
        tools_menu = menubar.addMenu('&Tools')
        
        # Check rMATS Installation
        check_rmats_action = QAction('Check &rMATS Installation', self)
        check_rmats_action.setStatusTip('Verify rMATS installation and version')
        check_rmats_action.triggered.connect(self.check_rmats_installation)
        tools_menu.addAction(check_rmats_action)
        
        # System Information
        debug_action = QAction('&System Information', self)
        debug_action.setStatusTip('Show system and environment information')
        debug_action.triggered.connect(self.show_debug_info)
        tools_menu.addAction(debug_action)
        
        tools_menu.addSeparator()
        
        # Set Default Directories
        set_defaults_action = QAction('Set &Default Directories...', self)
        set_defaults_action.setStatusTip('Configure default output and temp directories')
        set_defaults_action.triggered.connect(self.set_default_directories)
        tools_menu.addAction(set_defaults_action)
        
        # View Menu
        view_menu = menubar.addMenu('&View')
        
        # Toggle Status Bar
        status_bar_action = QAction('&Status Bar', self)
        status_bar_action.setCheckable(True)
        status_bar_action.setChecked(True)
        status_bar_action.setStatusTip('Show or hide the status bar')
        status_bar_action.triggered.connect(self.toggle_status_bar)
        view_menu.addAction(status_bar_action)
        
        # Help Menu
        help_menu = menubar.addMenu('&Help')
        
        # rMATS Documentation
        docs_action = QAction('rMATS &Documentation', self)
        docs_action.setStatusTip('Open rMATS documentation')
        docs_action.triggered.connect(self.open_documentation)
        help_menu.addAction(docs_action)
        
        help_menu.addSeparator()
        
        # About
        about_action = QAction('&About', self)
        about_action.setStatusTip('About this application')
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)

    def new_analysis(self):
        """Clear all fields for a new analysis"""
        self.b1_edit.clear()
        self.b2_edit.clear()
        self.gtf_edit.clear()
        self.out_edit.clear()
        self.temp_edit.clear()
        self.output_text.clear()
        self.statusBar.showMessage("New analysis - all fields cleared")

    def load_configuration(self):
        """Load configuration from a file"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, 
            "Load Configuration", 
            "", 
            "Configuration Files (*.conf *.cfg *.txt);;All Files (*)"
        )
        if file_path:
            try:
                with open(file_path, 'r') as f:
                    # Simple key=value format
                    for line in f:
                        line = line.strip()
                        if '=' in line and not line.startswith('#'):
                            key, value = line.split('=', 1)
                            key = key.strip()
                            value = value.strip()
                            
                            if key == 'batch1':
                                self.b1_edit.setText(value)
                            elif key == 'batch2':
                                self.b2_edit.setText(value)
                            elif key == 'gtf':
                                self.gtf_edit.setText(value)
                            elif key == 'output':
                                self.out_edit.setText(value)
                            elif key == 'temp':
                                self.temp_edit.setText(value)
                            elif key == 'read_type':
                                self.type_combo.setCurrentText(value)
                            elif key == 'read_length':
                                self.read_len_spin.setValue(int(value))
                            elif key == 'threads':
                                self.threads_spin.setValue(int(value))
                            elif key == 'variable_read':
                                self.var_read_check.setChecked(value.lower() == 'true')
                            elif key == 'statoff':
                                self.statoff_check.setChecked(value.lower() == 'true')
                
                self.statusBar.showMessage(f"Configuration loaded from {os.path.basename(file_path)}")
            except Exception as e:
                QMessageBox.warning(self, "Error", f"Failed to load configuration:\n{str(e)}")

    def save_configuration(self):
        """Save current configuration to a file"""
        file_path, _ = QFileDialog.getSaveFileName(
            self, 
            "Save Configuration", 
            "rmats_config.conf",
            "Configuration Files (*.conf *.cfg *.txt);;All Files (*)"
        )
        if file_path:
            try:
                with open(file_path, 'w') as f:
                    f.write("# rMATS Analysis Configuration\n")
                    f.write(f"batch1={self.b1_edit.text()}\n")
                    f.write(f"batch2={self.b2_edit.text()}\n")
                    f.write(f"gtf={self.gtf_edit.text()}\n")
                    f.write(f"output={self.out_edit.text()}\n")
                    f.write(f"temp={self.temp_edit.text()}\n")
                    f.write(f"read_type={self.type_combo.currentText()}\n")
                    f.write(f"read_length={self.read_len_spin.value()}\n")
                    f.write(f"threads={self.threads_spin.value()}\n")
                    f.write(f"variable_read={self.var_read_check.isChecked()}\n")
                    f.write(f"statoff={self.statoff_check.isChecked()}\n")
                
                self.statusBar.showMessage(f"Configuration saved to {os.path.basename(file_path)}")
            except Exception as e:
                QMessageBox.warning(self, "Error", f"Failed to save configuration:\n{str(e)}")

    def export_results(self):
        """Export analysis output to a file"""
        if not self.output_text.toPlainText().strip():
            QMessageBox.information(self, "No Results", "No analysis output to export.")
            return
            
        file_path, _ = QFileDialog.getSaveFileName(
            self, 
            "Export Results", 
            "rmats_output.txt",
            "Text Files (*.txt);;Log Files (*.log);;All Files (*)"
        )
        if file_path:
            try:
                with open(file_path, 'w') as f:
                    f.write(self.output_text.toPlainText())
                self.statusBar.showMessage(f"Results exported to {os.path.basename(file_path)}")
            except Exception as e:
                QMessageBox.warning(self, "Error", f"Failed to export results:\n{str(e)}")

    def stop_analysis(self):
        """Stop the running analysis"""
        if self.process.state() == QProcess.Running:
            self.process.kill()
            self.output_text.append("\n<font color='orange'>Analysis stopped by user</font>")
            self.statusBar.showMessage("Analysis stopped")

    def validate_inputs(self):
        """Validate all input fields"""
        errors = []
        warnings = []
        
        # Check required fields
        if not self.b1_edit.text().strip():
            errors.append("Batch 1 file is required")
        elif not os.path.exists(self.b1_edit.text().strip()):
            errors.append("Batch 1 file does not exist")
            
        if not self.statoff_check.isChecked():
            if not self.b2_edit.text().strip():
                errors.append("Batch 2 file is required when statistical analysis is enabled")
            elif not os.path.exists(self.b2_edit.text().strip()):
                errors.append("Batch 2 file does not exist")
        
        if not self.gtf_edit.text().strip():
            errors.append("GTF file is required")
        elif not os.path.exists(self.gtf_edit.text().strip()):
            errors.append("GTF file does not exist")
            
        if not self.out_edit.text().strip():
            errors.append("Output directory is required")
        elif not os.path.exists(self.out_edit.text().strip()):
            warnings.append("Output directory does not exist (will be created)")
            
        if not self.temp_edit.text().strip():
            errors.append("Temp directory is required")
        elif not os.path.exists(self.temp_edit.text().strip()):
            warnings.append("Temp directory does not exist (will be created)")
        
        # Show results
        if errors:
            QMessageBox.critical(self, "Validation Errors", 
                               f"Found {len(errors)} error(s):\n\n" + "\n".join(f"• {e}" for e in errors))
        elif warnings:
            QMessageBox.warning(self, "Validation Warnings", 
                              f"Found {len(warnings)} warning(s):\n\n" + "\n".join(f"• {w}" for w in warnings))
        else:
            QMessageBox.information(self, "Validation Success", "All inputs are valid!")

    def clear_output(self):
        """Clear the output window"""
        self.output_text.clear()
        self.statusBar.showMessage("Output cleared")

    def check_rmats_installation(self):
        """Check rMATS installation and show detailed information"""
        try:
            result = subprocess.run(
                ["rmats.py", "--version"],
                capture_output=True, 
                text=True,
                timeout=10
            )
            
            if result.returncode == 0:
                version_info = result.stdout.strip() or "Version information not available"
                QMessageBox.information(self, "rMATS Installation", 
                                      f"rMATS is installed and working!\n\n{version_info}")
            else:
                QMessageBox.warning(self, "rMATS Installation", 
                                  f"rMATS found but version check failed.\nReturn code: {result.returncode}")
        except FileNotFoundError:
            QMessageBox.critical(self, "rMATS Installation", 
                               "rMATS not found!\n\nPlease ensure rMATS is installed and 'rmats.py' is in your PATH.")
        except Exception as e:
            QMessageBox.critical(self, "rMATS Installation", f"Error checking rMATS:\n{str(e)}")

    def set_default_directories(self):
        """Set default output and temp directories"""
        # This is a placeholder - you could implement a dialog to set and save default directories
        documents_path = QStandardPaths.writableLocation(QStandardPaths.DocumentsLocation)
        default_output = os.path.join(documents_path, "rMATS_Output")
        default_temp = os.path.join(documents_path, "rMATS_Temp")
        
        msg = QMessageBox(self)
        msg.setWindowTitle("Default Directories")
        msg.setText(f"Suggested default directories:\n\nOutput: {default_output}\nTemp: {default_temp}\n\nWould you like to set these as defaults?")
        msg.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        
        if msg.exec() == QMessageBox.Yes:
            if not self.out_edit.text().strip():
                self.out_edit.setText(default_output)
            if not self.temp_edit.text().strip():
                self.temp_edit.setText(default_temp)
            self.statusBar.showMessage("Default directories set")

    def toggle_status_bar(self, checked):
        """Toggle status bar visibility"""
        self.statusBar.setVisible(checked)

    def open_documentation(self):
        """Open rMATS documentation"""
        import webbrowser
        webbrowser.open("https://rmats.sourceforge.net/")
        self.statusBar.showMessage("Opening rMATS documentation in browser")

    def show_about(self):
        """Show about dialog"""
        about_text = """
        <h3>rMATS Alternative Splicing Analyzer</h3>
        <p>A graphical user interface for rMATS (replicate Multivariate Analysis of Transcript Splicing)</p>
        
        <p><b>Features:</b></p>
        <ul>
        <li>Easy-to-use GUI for rMATS analysis</li>
        <li>Configuration save/load functionality</li>
        <li>Real-time analysis monitoring</li>
        <li>Comprehensive input validation</li>
        </ul>
        
        <p><b>rMATS:</b> A computational tool to detect differential alternative splicing events from replicate RNA-Seq data.</p>
        
        <p>Built with PySide6/Qt</p>
        """
        
        msg = QMessageBox(self)
        msg.setWindowTitle("About rMATS GUI")
        msg.setText(about_text)
        msg.setIcon(QMessageBox.Information)
        msg.exec()

    def toggle_b2_requirement(self):
        """Update UI based on statoff selection"""
        if self.statoff_check.isChecked():
            self.b2_label.setText("Batch 2 File (optional):")
            pal = self.b2_edit.palette()
            pal.setColor(QPalette.Base, QColor(240, 240, 240))
            self.b2_edit.setPalette(pal)
        else:
            self.b2_label.setText("Batch 2 File*:")
            self.b2_edit.setPalette(QApplication.palette())

    def browse_file(self, line_edit):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select File")
        if file_path:
            line_edit.setText(file_path)

    def browse_directory(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Directory")
        if dir_path:
            self.out_edit.setText(dir_path)

    def browse_temp_directory(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Temp Directory")
        if dir_path:
            self.temp_edit.setText(dir_path)

    def run_analysis(self):
        # Reset UI
        self.output_text.clear()
        self.statusBar.showMessage("Starting...")
        self.progress_bar.setValue(0)
        
        # Enable stop action, disable run
        self.stop_action.setEnabled(True)
        
        # Validate required inputs
        required_fields = [
            ("Batch 1", self.b1_edit.text()),
            ("GTF File", self.gtf_edit.text()),
            ("Output Directory", self.out_edit.text()),
            ("Temp Directory", self.temp_edit.text()),
        ]
        
        # Batch 2 is required only when statoff is disabled
        if not self.statoff_check.isChecked():
            required_fields.append(("Batch 2", self.b2_edit.text()))
        
        # Check for empty required fields
        errors = [name for name, value in required_fields if not value.strip()]
        if errors:
            error_msg = f"ERROR: Missing required fields: {', '.join(errors)}"
            self.output_text.append(f"<font color='red'>{error_msg}</font>")
            self.statusBar.showMessage(error_msg)
            self.stop_action.setEnabled(False)
            return

        # Build command with proper argument order
        cmd = [
            "rmats.py",
            "--b1", self.b1_edit.text().strip(),
        ]

        # Add --b2 immediately after --b1 if needed
        if self.b2_edit.text().strip() or not self.statoff_check.isChecked():
            cmd.extend(["--b2", self.b2_edit.text().strip()])

        # Add the rest of the parameters
        cmd.extend([
            "--gtf", self.gtf_edit.text().strip(),
            "-t", self.type_combo.currentText(),
            "--readLength", str(self.read_len_spin.value()),
            "--nthread", str(self.threads_spin.value()),
            "--od", self.out_edit.text().strip(),
            "--tmp", self.temp_edit.text().strip()
        ])

        # Add optional flags at the end
        if self.var_read_check.isChecked():
            cmd.append("--variable-read-length")
        
        if self.statoff_check.isChecked():
            cmd.append("--statoff")

        # Display command
        full_cmd = " ".join(cmd)
        self.output_text.append("Running command:\n" + full_cmd)
        self.output_text.append("\nChecking rMATS installation with 'rmats.py -h'...")
        self.statusBar.showMessage("Checking rMATS installation")
        
        # Check if rMATS is installed and functional
        try:
            # Run help command to test installation
            result = subprocess.run(
                ["rmats.py", "-h"],
                capture_output=True, 
                text=True,
                timeout=10  # Add timeout to prevent hanging
            )
            
            # Check for successful execution
            if result.returncode == 0 or "usage: rmats.py" in result.stdout:
                self.output_text.append("rMATS installation verified (help command works)")
                if result.stdout:
                    self.output_text.append("Help command output (first line):")
                    self.output_text.append(result.stdout.split('\n')[0] if result.stdout else "No output")
                self.statusBar.showMessage("rMATS found, starting analysis")
            else:
                error_msg = f"rMATS help command failed: exit code {result.returncode}"
                if result.stderr:
                    error_msg += f" | stderr: {result.stderr.strip()[:200]}"  # Show first 200 chars
                self.output_text.append(f"<font color='red'>{error_msg}</font>")
                self.statusBar.showMessage("rMATS not functional")
                self.stop_action.setEnabled(False)
                return
        except FileNotFoundError:
            self.output_text.append("<font color='red'>rMATS not found: 'rmats.py' command not available</font>")
            self.statusBar.showMessage("rMATS not found")
            self.stop_action.setEnabled(False)
            return
        except subprocess.TimeoutExpired:
            self.output_text.append("<font color='red'>rMATS check timed out (is it installed?)</font>")
            self.statusBar.showMessage("rMATS check timeout")
            self.stop_action.setEnabled(False)
            return
        except Exception as e:
            self.output_text.append(f"<font color='red'>Error checking rMATS: {str(e)}</font>")
            self.statusBar.showMessage("rMATS check failed")
            self.stop_action.setEnabled(False)
            return

        # Show process starting
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)  # Indeterminate progress
        self.run_button.setEnabled(False)
        self.statusBar.showMessage("Starting process...")
        
        # Clear previous output and prepare for new run
        self.output_text.append("\nStarting rMATS analysis...")
        self.output_text.append("-" * 80)
        
        # Start process
        try:
            # Start process monitoring
            self.timer.start(1000)  # Check every second
            
            # Start the actual process using the current environment
            self.process.setProcessChannelMode(QProcess.MergedChannels)  # Merge stdout and stderr
            self.process.start("rmats.py", cmd[1:])
                
        except Exception as e:
            self.output_text.append(f"<font color='red'>Failed to start process: {str(e)}</font>")
            self.statusBar.showMessage(f"Start failed: {str(e)}")
            self.run_button.setEnabled(True)
            self.progress_bar.setVisible(False)
            self.stop_action.setEnabled(False)

    def process_started(self):
        """Called when process successfully starts"""
        self.statusBar.showMessage("Analysis running...")
        self.output_text.append("\nProcess started successfully")
        self.output_text.append("Waiting for output...\n")
        self.output_text.moveCursor(QTextCursor.End)

    def process_error(self, error):
        """Handle process errors"""
        error_types = {
            QProcess.FailedToStart: "Process failed to start",
            QProcess.Crashed: "Process crashed",
            QProcess.Timedout: "Process timed out",
            QProcess.WriteError: "Write error",
            QProcess.ReadError: "Read error",
            QProcess.UnknownError: "Unknown error"
        }
        error_msg = error_types.get(error, f"Error code: {error}")
        self.output_text.append(f"<font color='red'>PROCESS ERROR: {error_msg}</font>")
        self.statusBar.showMessage(f"Error: {error_msg}")
        self.run_button.setEnabled(True)
        self.progress_bar.setVisible(False)
        self.stop_action.setEnabled(False)
        self.timer.stop()

    def monitor_process(self):
        """Monitor process state"""
        state = self.process.state()
        state_names = {
            QProcess.NotRunning: "Not running",
            QProcess.Starting: "Starting...",
            QProcess.Running: "Running"
        }
        self.statusBar.showMessage(f"Status: {state_names.get(state, 'Unknown')}")
        
        # If process has finished but our finished handler didn't trigger
        if state == QProcess.NotRunning and self.process.exitStatus() == QProcess.CrashExit:
            self.output_text.append("<font color='red'>Process crashed unexpectedly</font>")
            self.statusBar.showMessage("Process crashed")
            self.run_button.setEnabled(True)
            self.progress_bar.setVisible(False)
            self.stop_action.setEnabled(False)
            self.timer.stop()

    def handle_stdout(self):
        """Handle all process output (stdout and stderr)"""
        while self.process.canReadLine():
            data = self.process.readLine()
            if data:
                try:
                    # Try UTF-8 decoding first
                    text = data.data().decode('utf-8').strip()
                except UnicodeDecodeError:
                    try:
                        # Fallback to system encoding
                        text = data.data().decode(sys.getfilesystemencoding()).strip()
                    except:
                        # If all else fails, use raw bytes
                        text = str(data.data())
                
                # Append to output with preserved line breaks
                self.output_text.append(text)
                self.output_text.moveCursor(QTextCursor.End)

    def handle_stderr(self):
        # Not used since we're using merged channels
        pass

    def analysis_finished(self, exit_code, exit_status):
        self.timer.stop()
        self.progress_bar.setVisible(False)
        self.run_button.setEnabled(True)
        self.stop_action.setEnabled(False)
        
        if exit_status == QProcess.NormalExit and exit_code == 0:
            self.output_text.append("\nAnalysis completed successfully!")
            self.statusBar.showMessage("Analysis completed successfully")
        else:
            self.output_text.append(f"\n<font color='red'>Analysis failed with exit code {exit_code}</font>")
            self.statusBar.showMessage(f"Analysis failed (exit code {exit_code})")
            
            # Try to capture any final error messages
            error_data = self.process.readAllStandardError()
            if error_data:
                error_text = error_data.data().decode('utf-8', 'ignore').strip()
                if error_text:
                    self.output_text.append(f"<font color='red'>Final error: {error_text}</font>")

    def show_debug_info(self):
        """Show diagnostic information"""
        info = [
            f"Python version: {sys.version}",
            f"Platform: {sys.platform}",
            f"Working directory: {os.getcwd()}",
            f"PATH environment: {os.getenv('PATH', 'Not set')}"
        ]
        
        # Check rMATS availability
        try:
            result = subprocess.run(
                ["which", "rmats.py"],
                capture_output=True,
                text=True
            )
            if result.returncode == 0:
                info.append(f"rMATS path: {result.stdout.strip()}")
            else:
                info.append("rMATS not found in PATH")
        except Exception as e:
            info.append(f"rMATS check error: {str(e)}")
            
        # Show dialog with debug info
        msg = QMessageBox(self)
        msg.setWindowTitle("System Information")
        msg.setText("\n".join(info))
        msg.setIcon(QMessageBox.Information)
        msg.exec()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = rMATSGUI()
    window.show()
    sys.exit(app.exec())