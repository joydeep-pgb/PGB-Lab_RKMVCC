# requirements: PySide6, matplotlib, seaborn, pandas, openpyxl
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QComboBox, QCheckBox, QLabel, QDockWidget, QPushButton, QSpinBox,
    QTableWidget, QTableWidgetItem, QFileDialog, QMessageBox, QTabWidget,
    QColorDialog, QGroupBox, QGridLayout, QLineEdit, QTextEdit, QSplitter,
    QDoubleSpinBox, QSlider, QFrame, QScrollArea
)
from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QColor, QFont, QAction
import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

import seaborn as sns
import pandas as pd
import numpy as np
import sys
import os

class DataTableWidget(QWidget):
    """Data import and editing widget"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent_window = parent
        layout = QVBoxLayout(self)
        
        # Import buttons
        btn_layout = QHBoxLayout()
        self.import_csv_btn = QPushButton("Import CSV")
        self.import_excel_btn = QPushButton("Import Excel")
        self.sample_data_btn = QPushButton("Load Sample Data")
        self.clear_btn = QPushButton("Clear Data")
        
        btn_layout.addWidget(self.import_csv_btn)
        btn_layout.addWidget(self.import_excel_btn)
        btn_layout.addWidget(self.sample_data_btn)
        btn_layout.addWidget(self.clear_btn)
        btn_layout.addStretch()
        
        layout.addLayout(btn_layout)
        
        # Data table
        self.table = QTableWidget()
        self.table.setMinimumHeight(200)  # Reduced minimum height
        layout.addWidget(self.table)
        
        # Connect signals
        self.import_csv_btn.clicked.connect(self.import_csv)
        self.import_excel_btn.clicked.connect(self.import_excel)
        self.sample_data_btn.clicked.connect(self.load_sample_data)
        self.clear_btn.clicked.connect(self.clear_data)
        self.table.cellChanged.connect(self.on_data_changed)
        
        # Do NOT auto-load sample data initially
        # self.load_sample_data()
    
    def import_csv(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Import CSV File", "", "CSV Files (*.csv)"
        )
        if file_path:
            try:
                df = pd.read_csv(file_path)
                self.load_dataframe(df)
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to import CSV:\n{str(e)}")
    
    def import_excel(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Import Excel File", "", "Excel Files (*.xlsx *.xls)"
        )
        if file_path:
            try:
                df = pd.read_excel(file_path)
                self.load_dataframe(df)
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to import Excel:\n{str(e)}")
    
    def load_sample_data(self):
        """Load sample data similar to GraphPad Prism examples"""
        rng = np.random.default_rng(42)
        n = 100
        
        # Create more realistic sample data
        groups = ["Control", "Treatment A", "Treatment B", "Treatment C"]
        data = []
        
        for i, group in enumerate(groups):
            base_value = 10 + i * 2
            values = rng.normal(base_value, 2, n//4)
            for val in values:
                data.append({
                    "Group": group,
                    "Value": val,
                    "Time": rng.uniform(0, 24),
                    "Dose": [0, 5, 10, 20][i],
                    "Response": val * (1 + i * 0.2) + rng.normal(0, 1)
                })
        
        df = pd.DataFrame(data)
        self.load_dataframe(df)
    
    def clear_data(self):
        self.table.setRowCount(0)
        self.table.setColumnCount(0)
        if self.parent_window:
            self.parent_window.data_changed()
    
    def load_dataframe(self, df):
        """Load a pandas DataFrame into the table"""
        self.table.setRowCount(len(df))
        self.table.setColumnCount(len(df.columns))
        self.table.setHorizontalHeaderLabels(df.columns.tolist())
        
        for i, row in df.iterrows():
            for j, value in enumerate(row):
                item = QTableWidgetItem(str(value))
                self.table.setItem(i, j, item)
        
        if self.parent_window and hasattr(self.parent_window, "data_changed"):
            self.parent_window.data_changed()
    
    def get_dataframe(self):
        """Extract current table data as pandas DataFrame"""
        if self.table.rowCount() == 0 or self.table.columnCount() == 0:
            return pd.DataFrame()
        
        columns = [self.table.horizontalHeaderItem(i).text() 
                  for i in range(self.table.columnCount())]
        
        data = []
        for row in range(self.table.rowCount()):
            row_data = {}
            for col in range(self.table.columnCount()):
                item = self.table.item(row, col)
                value = item.text() if item else ""
                
                # Try to convert to numeric
                try:
                    value = float(value)
                except ValueError:
                    pass  # Keep as string
                
                row_data[columns[col]] = value
            data.append(row_data)
        
        return pd.DataFrame(data)
    
    def on_data_changed(self):
        """Called when table data is modified"""
        if self.parent_window:
            QTimer.singleShot(100, self.parent_window.data_changed)

class StyleControlWidget(QWidget):
    """Advanced styling controls"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent_window = parent
        
        # Create a scroll area to contain all controls
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        
        # Create a container widget for the scroll area
        container = QWidget()
        layout = QVBoxLayout(container)
        
        # Plot type group
        plot_group = QGroupBox("Plot Type")
        plot_layout = QVBoxLayout(plot_group)
        
        plot_layout.addWidget(QLabel("Kind:"))
        self.kind = QComboBox()
        self.kind.addItems([
            "scatter", "line", "hist", "kde", "box", "violin", 
            "bar", "strip", "swarm", "point", "count"
        ])
        plot_layout.addWidget(self.kind)
        
        layout.addWidget(plot_group)
        
        # Variables group
        var_group = QGroupBox("Variables")
        var_layout = QGridLayout(var_group)
        
        var_layout.addWidget(QLabel("X:"), 0, 0)
        self.xcol = QComboBox()
        var_layout.addWidget(self.xcol, 0, 1)
        
        var_layout.addWidget(QLabel("Y:"), 1, 0)
        self.ycol = QComboBox()
        var_layout.addWidget(self.ycol, 1, 1)
        
        var_layout.addWidget(QLabel("Hue:"), 2, 0)
        self.hue = QComboBox()
        var_layout.addWidget(self.hue, 2, 1)
        
        var_layout.addWidget(QLabel("Size:"), 3, 0)
        self.size_col = QComboBox()
        var_layout.addWidget(self.size_col, 3, 1)
        
        var_layout.addWidget(QLabel("Style:"), 4, 0)
        self.style_col = QComboBox()
        var_layout.addWidget(self.style_col, 4, 1)
        
        layout.addWidget(var_group)
        
        # Appearance group
        appear_group = QGroupBox("Appearance")
        appear_layout = QVBoxLayout(appear_group)
        
        # Seaborn style and context
        appear_layout.addWidget(QLabel("Context:"))
        self.context = QComboBox()
        self.context.addItems(["paper", "notebook", "talk", "poster"])
        self.context.setCurrentText("notebook")
        appear_layout.addWidget(self.context)
        
        appear_layout.addWidget(QLabel("Style:"))
        self.sns_style = QComboBox()
        self.sns_style.addItems(["white", "dark", "whitegrid", "darkgrid", "ticks"])
        appear_layout.addWidget(self.sns_style)
        
        appear_layout.addWidget(QLabel("Palette:"))
        self.palette = QComboBox()
        self.palette.addItems([
            "deep", "muted", "pastel", "bright", "dark", "colorblind",
            "husl", "rocket", "mako", "flare", "crest", "viridis", "plasma", "custom"
        ])
        appear_layout.addWidget(self.palette)
        
        # Custom color button
        self.custom_color_btn = QPushButton("Choose Custom Colors")
        self.custom_color_btn.clicked.connect(self.choose_custom_colors)
        appear_layout.addWidget(self.custom_color_btn)
        
        layout.addWidget(appear_group)
        
        # Font group
        font_group = QGroupBox("Font Settings")
        font_layout = QGridLayout(font_group)
        
        font_layout.addWidget(QLabel("Title Size:"), 0, 0)
        self.title_size = QSpinBox()
        self.title_size.setRange(8, 32)
        self.title_size.setValue(14)
        font_layout.addWidget(self.title_size, 0, 1)
        
        font_layout.addWidget(QLabel("Label Size:"), 1, 0)
        self.label_size = QSpinBox()
        self.label_size.setRange(6, 24)
        self.label_size.setValue(12)
        font_layout.addWidget(self.label_size, 1, 1)
        
        font_layout.addWidget(QLabel("Tick Size:"), 2, 0)
        self.tick_size = QSpinBox()
        self.tick_size.setRange(6, 20)
        self.tick_size.setValue(10)
        font_layout.addWidget(self.tick_size, 2, 1)
        
        font_layout.addWidget(QLabel("Legend Size:"), 3, 0)
        self.legend_size = QSpinBox()
        self.legend_size.setRange(6, 20)
        self.legend_size.setValue(10)
        font_layout.addWidget(self.legend_size, 3, 1)
        
        layout.addWidget(font_group)
        
        # Size and spacing group
        size_group = QGroupBox("Plot Settings")
        size_layout = QGridLayout(size_group)
        
        size_layout.addWidget(QLabel("Figure Width:"), 0, 0)
        self.fig_width = QDoubleSpinBox()
        self.fig_width.setRange(4, 20)
        self.fig_width.setValue(8)
        self.fig_width.setSingleStep(0.5)
        size_layout.addWidget(self.fig_width, 0, 1)
        
        size_layout.addWidget(QLabel("Figure Height:"), 1, 0)
        self.fig_height = QDoubleSpinBox()
        self.fig_height.setRange(3, 15)
        self.fig_height.setValue(6)
        self.fig_height.setSingleStep(0.5)
        size_layout.addWidget(self.fig_height, 1, 1)
        
        size_layout.addWidget(QLabel("Point Size:"), 2, 0)
        self.point_size = QSpinBox()
        self.point_size.setRange(10, 200)
        self.point_size.setValue(50)
        size_layout.addWidget(self.point_size, 2, 1)
        
        size_layout.addWidget(QLabel("Line Width:"), 3, 0)
        self.line_width = QDoubleSpinBox()
        self.line_width.setRange(0.5, 5.0)
        self.line_width.setValue(1.5)
        self.line_width.setSingleStep(0.1)
        size_layout.addWidget(self.line_width, 3, 1)
        
        layout.addWidget(size_group)
        
        # Toggles group
        toggle_group = QGroupBox("Display Options")
        toggle_layout = QVBoxLayout(toggle_group)
        
        self.show_grid = QCheckBox("Show Grid")
        self.show_grid.setChecked(True)
        toggle_layout.addWidget(self.show_grid)
        
        self.show_legend = QCheckBox("Show Legend")
        self.show_legend.setChecked(True)
        toggle_layout.addWidget(self.show_legend)
        
        self.show_spines = QCheckBox("Show Spines")
        self.show_spines.setChecked(True)
        toggle_layout.addWidget(self.show_spines)
        
        self.tight_layout = QCheckBox("Tight Layout")
        self.tight_layout.setChecked(True)
        toggle_layout.addWidget(self.tight_layout)
        
        layout.addWidget(toggle_group)
        
        # Custom title and labels
        labels_group = QGroupBox("Labels")
        labels_layout = QVBoxLayout(labels_group)
        
        labels_layout.addWidget(QLabel("Title:"))
        self.custom_title = QLineEdit()
        labels_layout.addWidget(self.custom_title)
        
        labels_layout.addWidget(QLabel("X Label:"))
        self.custom_xlabel = QLineEdit()
        labels_layout.addWidget(self.custom_xlabel)
        
        labels_layout.addWidget(QLabel("Y Label:"))
        self.custom_ylabel = QLineEdit()
        labels_layout.addWidget(self.custom_ylabel)
        
        layout.addWidget(labels_group)
        
        # Apply and export buttons
        btn_layout2 = QVBoxLayout()
        self.apply_btn = QPushButton("Apply Changes")
        self.apply_btn.setStyleSheet("QPushButton { background-color: #4CAF50; color: white; font-weight: bold; }")
        self.export_pdf_btn = QPushButton("Export to PDF")
        self.export_png_btn = QPushButton("Export to PNG")
        
        btn_layout2.addWidget(self.apply_btn)
        btn_layout2.addWidget(self.export_pdf_btn)
        btn_layout2.addWidget(self.export_png_btn)
        
        layout.addLayout(btn_layout2)
        layout.addStretch()
        
        # Set the container as the scroll area's widget
        scroll.setWidget(container)
        
        # Create a main layout for this widget and add the scroll area
        main_layout = QVBoxLayout(self)
        main_layout.addWidget(scroll)
        
        # Custom colors storage
        self.custom_colors = []
        
        # Connect export signals
        self.export_pdf_btn.clicked.connect(self.export_pdf)
        self.export_png_btn.clicked.connect(self.export_png)
    
    def choose_custom_colors(self):
        """Open color picker for custom palette"""
        colors = []
        for i in range(6):  # Allow up to 6 custom colors
            color = QColorDialog.getColor(QColor(255, 255, 255), self)
            if color.isValid():
                colors.append(color.name())
            else:
                break
        
        if colors:
            self.custom_colors = colors
            QMessageBox.information(self, "Custom Colors", f"Set {len(colors)} custom colors")
    
    def update_column_combos(self, columns):
        """Update column selection combos when data changes"""
        combos = [self.xcol, self.ycol, self.hue, self.size_col, self.style_col]
        
        for combo in combos:
            current = combo.currentText()
            combo.clear()
            combo.addItem("")  # Empty option
            combo.addItems(columns)
            
            # Try to restore previous selection
            idx = combo.findText(current)
            if idx >= 0:
                combo.setCurrentIndex(idx)
    
    def export_pdf(self):
        """Export current plot to editable PDF"""
        if not self.parent_window:
            return
            
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export to PDF", "plot.pdf", "PDF Files (*.pdf)"
        )
        
        if file_path:
            try:
                # Create a new figure with the same settings for export
                fig = Figure(figsize=(self.fig_width.value(), self.fig_height.value()))
                ax = fig.add_subplot(111)
                
                # Redraw the plot on the new figure
                self.parent_window.draw_plot(ax, for_export=True)
                
                # Save as editable PDF
                with PdfPages(file_path) as pdf:
                    pdf.savefig(fig, bbox_inches='tight', 
                               facecolor='white', edgecolor='none')
                
                QMessageBox.information(self, "Export Success", 
                                      f"Plot exported to:\n{file_path}")
                
            except Exception as e:
                QMessageBox.critical(self, "Export Error", 
                                   f"Failed to export PDF:\n{str(e)}")
    
    def export_png(self):
        """Export current plot to PNG"""
        if not self.parent_window:
            return
            
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export to PNG", "plot.png", "PNG Files (*.png)"
        )
        
        if file_path:
            try:
                # Create a new figure with high DPI for export
                fig = Figure(figsize=(self.fig_width.value(), self.fig_height.value()), 
                           dpi=300)
                ax = fig.add_subplot(111)
                
                # Redraw the plot on the new figure
                self.parent_window.draw_plot(ax, for_export=True)
                
                # Save as PNG
                fig.savefig(file_path, bbox_inches='tight', dpi=300,
                           facecolor='white', edgecolor='none')
                
                QMessageBox.information(self, "Export Success", 
                                      f"Plot exported to:\n{file_path}")
                
            except Exception as e:
                QMessageBox.critical(self, "Export Error", 
                                   f"Failed to export PNG:\n{str(e)}")

class PlotCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None):
        self.fig = Figure(constrained_layout=True)
        super().__init__(self.fig)
        self.ax = self.fig.add_subplot(111)
        self.setParent(parent)

class StatisticsWidget(QWidget):
    """Statistics and analysis panel"""
    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        
        # Statistics display
        self.stats_text = QTextEdit()
        self.stats_text.setMaximumHeight(200)
        self.stats_text.setReadOnly(True)
        layout.addWidget(QLabel("Data Statistics:"))
        layout.addWidget(self.stats_text)
        
        # Analysis options
        analysis_group = QGroupBox("Statistical Analysis")
        analysis_layout = QVBoxLayout(analysis_group)
        
        self.show_corr = QCheckBox("Show Correlation")
        self.show_regression = QCheckBox("Show Regression Line")
        self.show_ci = QCheckBox("Show Confidence Intervals")
        
        analysis_layout.addWidget(self.show_corr)
        analysis_layout.addWidget(self.show_regression)
        analysis_layout.addWidget(self.show_ci)
        
        layout.addWidget(analysis_group)
        layout.addStretch()
    
    def update_statistics(self, df):
        """Update statistics display"""
        if df.empty:
            self.stats_text.clear()
            return
        
        # Generate basic statistics
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) > 0:
            stats = df[numeric_cols].describe()
            stats_str = "Descriptive Statistics:\n\n"
            stats_str += stats.to_string()
            
            # Add correlation matrix if multiple numeric columns
            if len(numeric_cols) > 1:
                corr = df[numeric_cols].corr()
                stats_str += "\n\nCorrelation Matrix:\n"
                stats_str += corr.to_string()
            
            self.stats_text.setText(stats_str)
        else:
            self.stats_text.setText("No numeric data available for statistics.")

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Advanced Seaborn Plotter - GraphPad Prism Style")
        self.resize(1600, 1000)  # Increased initial size
        
        # Set application style
        self.setStyleSheet("""
            QMainWindow {
                background-color: #f5f5f5;
            }
            QDockWidget {
                font-weight: bold;
            }
            QGroupBox {
                font-weight: bold;
                border: 2px solid #cccccc;
                border-radius: 5px;
                margin-top: 1ex;
                padding-top: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px 0 5px;
            }
            QScrollArea {
                border: none;
            }
        """)
        
        # Central widget with splitter
        central = QWidget()
        main_layout = QHBoxLayout(central)
        
        # Create splitter for resizable panels
        splitter = QSplitter(Qt.Horizontal)
        
        # Left panel - Data table in tab widget
        left_tabs = QTabWidget()
        left_tabs.setMinimumWidth(300)  # Set minimum width
        
        # Data tab
        self.data_widget = DataTableWidget(self)
        left_tabs.addTab(self.data_widget, "Data")
        
        # Statistics tab
        self.stats_widget = StatisticsWidget(self)
        left_tabs.addTab(self.stats_widget, "Statistics")
        
        # Right panel - Plot
        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)
        
        self.canvas = PlotCanvas(self)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        
        right_layout.addWidget(self.toolbar)
        right_layout.addWidget(self.canvas)
        
        # Add panels to splitter
        splitter.addWidget(left_tabs)
        splitter.addWidget(right_widget)
        splitter.setSizes([500, 1100])  # Increased initial sizes
        
        main_layout.addWidget(splitter)
        self.setCentralWidget(central)
        
        # Style controls dock
        self.create_style_dock()
        
        # Create menu bar AFTER data_widget is initialized
        self.create_menu_bar()
        
        # Connect signals
        self.style_widget.apply_btn.clicked.connect(self.redraw)
        
        # Initial setup
        self.data_changed()
        self.redraw()
    
    def create_menu_bar(self):
        """Create menu bar like GraphPad Prism"""
        menubar = self.menuBar()
        
        # File menu
        file_menu = menubar.addMenu('File')
        
        new_action = QAction('New', self)
        new_action.triggered.connect(self.data_widget.clear_data)
        file_menu.addAction(new_action)
        
        file_menu.addSeparator()
        
        import_csv_action = QAction('Import CSV...', self)
        import_csv_action.triggered.connect(self.data_widget.import_csv)
        file_menu.addAction(import_csv_action)
        
        import_excel_action = QAction('Import Excel...', self)
        import_excel_action.triggered.connect(self.data_widget.import_excel)
        file_menu.addAction(import_excel_action)
        
        file_menu.addSeparator()
        
        export_pdf_action = QAction('Export PDF...', self)
        export_pdf_action.triggered.connect(self.style_widget.export_pdf)
        file_menu.addAction(export_pdf_action)
        
        export_png_action = QAction('Export PNG...', self)
        export_png_action.triggered.connect(self.style_widget.export_png)
        file_menu.addAction(export_png_action)
        
        # View menu
        view_menu = menubar.addMenu('View')
        
        toggle_style_dock = QAction('Toggle Style Panel', self)
        toggle_style_dock.triggered.connect(self.toggle_style_dock)
        view_menu.addAction(toggle_style_dock)
    
    def create_style_dock(self):
        """Create dockable style control panel"""
        self.style_dock = QDockWidget("Plot Controls", self)
        self.style_dock.setObjectName("styleDock")
        self.style_dock.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        self.style_dock.setMinimumWidth(400)  # Increased minimum width
        
        self.style_widget = StyleControlWidget(self)
        self.style_dock.setWidget(self.style_widget)
        self.addDockWidget(Qt.RightDockWidgetArea, self.style_dock)
    
    def toggle_style_dock(self):
        """Toggle style dock visibility"""
        self.style_dock.setVisible(not self.style_dock.isVisible())
    
    def data_changed(self):
        """Called when data table changes"""
        df = self.data_widget.get_dataframe()
        
        if not df.empty:
            # Update column combos
            numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
            all_cols = df.columns.tolist()
            categorical_cols = df.select_dtypes(include=['object', 'category']).columns.tolist()
            
            self.style_widget.update_column_combos(all_cols)
            
            # Set smart defaults
            if len(numeric_cols) >= 2:
                self.style_widget.xcol.setCurrentText(numeric_cols[0])
                self.style_widget.ycol.setCurrentText(numeric_cols[1])
            
            if categorical_cols:
                self.style_widget.hue.setCurrentText(categorical_cols[0])
            
            # Update statistics
            self.stats_widget.update_statistics(df)
        
        # Auto-redraw
        QTimer.singleShot(200, self.redraw)
    
    def get_current_palette(self):
        """Get current color palette"""
        palette_name = self.style_widget.palette.currentText()
        
        if palette_name == "custom" and self.style_widget.custom_colors:
            return self.style_widget.custom_colors
        else:
            return palette_name
    
    def draw_plot(self, ax, for_export=False):
        """Draw the plot on given axes"""
        df = self.data_widget.get_dataframe()
        
        if df.empty:
            ax.text(0.5, 0.5, 'No data available\nImport CSV/Excel or load sample data', 
                   ha='center', va='center', transform=ax.transAxes, fontsize=14)
            return
        
        # Clear the axes
        ax.clear()
        
        # Get plot parameters
        kind = self.style_widget.kind.currentText()
        x = self.style_widget.xcol.currentText() or None
        y = self.style_widget.ycol.currentText() or None
        hue = self.style_widget.hue.currentText() or None
        size_col = self.style_widget.size_col.currentText() or None
        style_col = self.style_widget.style_col.currentText() or None
        
        # Set figure size
        if not for_export:
            self.canvas.fig.set_size_inches(
                self.style_widget.fig_width.value(),
                self.style_widget.fig_height.value()
            )
        
        # Apply seaborn theme
        sns.set_theme(
            context=self.style_widget.context.currentText(),
            style=self.style_widget.sns_style.currentText(),
            palette=self.get_current_palette()
        )
        
        # Common plot parameters
        plot_kwargs = {
            'data': df,
            'x': x,
            'y': y,
            'hue': hue,
            'ax': ax
        }
        
        # Add size and style if applicable and available
        if size_col and kind in ['scatter', 'line', 'strip']:
            plot_kwargs['size'] = size_col
        if style_col and kind in ['scatter', 'line']:
            plot_kwargs['style'] = style_col
        
        try:
            # Plot based on kind
            if kind == "scatter":
                plot_kwargs['s'] = self.style_widget.point_size.value()
                if self.stats_widget.show_regression.isChecked():
                    sns.regplot(data=df, x=x, y=y, ax=ax, scatter=False, 
                               ci=95 if self.stats_widget.show_ci.isChecked() else None)
                sns.scatterplot(**plot_kwargs)
                
            elif kind == "line":
                plot_kwargs['linewidth'] = self.style_widget.line_width.value()
                sns.lineplot(**plot_kwargs)
                
            elif kind == "hist":
                plot_kwargs.pop('y', None)  # Remove y for histogram
                sns.histplot(**plot_kwargs, kde=True)
                
            elif kind == "kde":
                sns.kdeplot(**plot_kwargs, fill=True)
                
            elif kind == "box":
                if not hue:
                    plot_kwargs['x'] = hue or 'Group'
                sns.boxplot(**plot_kwargs)
                
            elif kind == "violin":
                if not hue:
                    plot_kwargs['x'] = hue or 'Group'
                sns.violinplot(**plot_kwargs)
                
            elif kind == "bar":
                if not hue:
                    plot_kwargs['x'] = hue or 'Group'
                sns.barplot(**plot_kwargs)
                
            elif kind == "strip":
                if not hue:
                    plot_kwargs['x'] = hue or 'Group'
                plot_kwargs['size'] = self.style_widget.point_size.value() / 10
                sns.stripplot(**plot_kwargs)
                
            elif kind == "swarm":
                if not hue:
                    plot_kwargs['x'] = hue or 'Group'
                plot_kwargs['size'] = self.style_widget.point_size.value() / 10
                sns.swarmplot(**plot_kwargs)
                
            elif kind == "point":
                sns.pointplot(**plot_kwargs)
                
            elif kind == "count":
                plot_kwargs.pop('y', None)
                sns.countplot(**plot_kwargs)
        
        except Exception as e:
            ax.text(0.5, 0.5, f'Error creating plot:\n{str(e)}', 
                   ha='center', va='center', transform=ax.transAxes, 
                   fontsize=12, color='red')
            return
        
        # Customize appearance
        self.apply_styling(ax)
    
    def apply_styling(self, ax):
        """Apply custom styling to the plot"""
        # Set custom labels if provided
        title = self.style_widget.custom_title.text()
        xlabel = self.style_widget.custom_xlabel.text()
        ylabel = self.style_widget.custom_ylabel.text()
        
        if not title:
            x_col = self.style_widget.xcol.currentText()
            y_col = self.style_widget.ycol.currentText()
            kind = self.style_widget.kind.currentText()
            title = f"{kind.title()} Plot"
            if y_col and x_col:
                title += f": {y_col} vs {x_col}"
        
        ax.set_title(title, fontsize=self.style_widget.title_size.value(), 
                    fontweight='bold', pad=20)
        
        if xlabel:
            ax.set_xlabel(xlabel, fontsize=self.style_widget.label_size.value(), 
                         fontweight='bold')
        if ylabel:
            ax.set_ylabel(ylabel, fontsize=self.style_widget.label_size.value(), 
                         fontweight='bold')
        
        # Apply font sizes
        ax.tick_params(axis='both', which='major', 
                      labelsize=self.style_widget.tick_size.value())
        
        # Grid settings
        ax.grid(self.style_widget.show_grid.isChecked(), alpha=0.3)
        
        # Legend settings
        legend = ax.get_legend()
        if legend:
            legend.set_visible(self.style_widget.show_legend.isChecked())
            if self.style_widget.show_legend.isChecked():
                legend.set_title(legend.get_title().get_text(), 
                               prop={'size': self.style_widget.legend_size.value(), 
                                     'weight': 'bold'})
                for text in legend.get_texts():
                    text.set_fontsize(self.style_widget.legend_size.value())
        
        # Spine settings
        if not self.style_widget.show_spines.isChecked():
            for spine in ax.spines.values():
                spine.set_visible(False)
        
        # Tight layout
        if self.style_widget.tight_layout.isChecked():
            self.canvas.fig.tight_layout()
        
        # Show correlation if requested
        if (self.stats_widget.show_corr.isChecked() and 
            hasattr(self, 'last_x_col') and hasattr(self, 'last_y_col')):
            df = self.data_widget.get_dataframe()
            x_col = self.style_widget.xcol.currentText()
            y_col = self.style_widget.ycol.currentText()
            
            if (x_col and y_col and x_col in df.columns and y_col in df.columns):
                try:
                    corr = df[x_col].corr(df[y_col])
                    ax.text(0.05, 0.95, f'Correlation: r = {corr:.3f}', 
                           transform=ax.transAxes, fontsize=10,
                           bbox=dict(boxstyle="round,pad=0.3", facecolor="white", 
                                    alpha=0.8))
                except:
                    pass
    
    def redraw(self):
        """Redraw the main plot"""
        self.draw_plot(self.canvas.ax)
        self.canvas.draw_idle()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    
    # Set application style
    app.setStyle('Fusion')  # Modern look
    
    # Apply dark theme colors similar to GraphPad Prism
    palette = app.palette()
    palette.setColor(palette.ColorRole.Window, QColor(240, 240, 240))
    palette.setColor(palette.ColorRole.WindowText, QColor(0, 0, 0))
    palette.setColor(palette.ColorRole.Base, QColor(255, 255, 255))
    palette.setColor(palette.ColorRole.AlternateBase, QColor(245, 245, 245))
    palette.setColor(palette.ColorRole.ToolTipBase, QColor(255, 255, 220))
    palette.setColor(palette.ColorRole.ToolTipText, QColor(0, 0, 0))
    palette.setColor(palette.ColorRole.Text, QColor(0, 0, 0))
    palette.setColor(palette.ColorRole.Button, QColor(240, 240, 240))
    palette.setColor(palette.ColorRole.ButtonText, QColor(0, 0, 0))
    palette.setColor(palette.ColorRole.BrightText, QColor(255, 0, 0))
    palette.setColor(palette.ColorRole.Link, QColor(42, 130, 218))
    palette.setColor(palette.ColorRole.Highlight, QColor(42, 130, 218))
    palette.setColor(palette.ColorRole.HighlightedText, QColor(255, 255, 255))
    app.setPalette(palette)
    
    window = MainWindow()
    window.show()
    
    sys.exit(app.exec())