import sys
from threading import Thread
from PySide6.QtCore import Signal, QObject
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget,
    QVBoxLayout, QLineEdit, QPushButton,
    QTextEdit, QProgressBar, QComboBox,
    QLabel, QFileDialog, QHBoxLayout, QCheckBox
)
from yt_dlp import YoutubeDL


class DownloadSignals(QObject):
    progress = Signal(float)
    log = Signal(str)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("YouTube Downloader")

        # URL input
        self.url_input = QLineEdit()
        self.url_input.setPlaceholderText("Paste YouTube URL here")

        # Quality selection
        self.quality_label = QLabel("Select Quality:")
        self.quality_box = QComboBox()
        self.quality_box.addItems(["Best Available", "1080p", "720p", "Audio Only"])

        # Playlist checkbox
        self.playlist_checkbox = QCheckBox("Download Playlist")

        # Output folder selection
        self.folder_label = QLabel("Output Folder: Not selected")
        self.browse_btn = QPushButton("Browse Folder")
        self.browse_btn.clicked.connect(self.choose_folder)
        self.output_folder = None

        folder_layout = QHBoxLayout()
        folder_layout.addWidget(self.folder_label)
        folder_layout.addWidget(self.browse_btn)

        # Download button & status
        self.download_btn = QPushButton("Download")
        self.status = QTextEdit()
        self.status.setReadOnly(True)
        self.progress_bar = QProgressBar()
        self.progress_bar.setValue(0)

        # Layout setup
        layout = QVBoxLayout()
        layout.addWidget(self.url_input)
        layout.addWidget(self.quality_label)
        layout.addWidget(self.quality_box)
        layout.addWidget(self.playlist_checkbox)
        layout.addLayout(folder_layout)
        layout.addWidget(self.download_btn)
        layout.addWidget(self.progress_bar)
        layout.addWidget(self.status)

        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)

        # Signals
        self.signals = DownloadSignals()
        self.signals.progress.connect(self.progress_bar.setValue)
        self.signals.log.connect(self.log)

        # Events
        self.download_btn.clicked.connect(self.start_download)

    def choose_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Select Download Folder")
        if folder:
            self.output_folder = folder
            self.folder_label.setText(f"Output Folder: {folder}")

    def start_download(self):
        url = self.url_input.text().strip()
        if not url:
            self.log("Please enter a valid URL.")
            return
        if not self.output_folder:
            self.log("Please choose an output folder.")
            return

        quality_choice = self.quality_box.currentText()
        is_playlist = self.playlist_checkbox.isChecked()
        Thread(target=self.download_video, args=(url, quality_choice, is_playlist), daemon=True).start()

    def download_video(self, url, quality_choice, is_playlist):
        self.signals.log.emit(f"Starting download: {url}")

        # Map quality choice to yt-dlp format string
        if quality_choice == "1080p":
            fmt = "bestvideo[height<=1080]+bestaudio/best[height<=1080]"
        elif quality_choice == "720p":
            fmt = "bestvideo[height<=720]+bestaudio/best[height<=720]"
        elif quality_choice == "Audio Only":
            fmt = "bestaudio/best"
        else:
            fmt = "best"

        ydl_opts = {
            'format': fmt,
            'outtmpl': f'{self.output_folder}/%(title)s.%(ext)s',
            'progress_hooks': [self.hook],
            'noplaylist': not is_playlist
        }

        try:
            with YoutubeDL(ydl_opts) as ydl:
                ydl.download([url])
            self.signals.log.emit("✅ Download completed!")
        except Exception as e:
            self.signals.log.emit(f"❌ Error: {e}")

    def hook(self, d):
        if d['status'] == 'downloading':
            percent_str = d.get('_percent_str', '0.0%').strip('%')
            try:
                percent = float(percent_str)
            except:
                percent = 0.0
            self.signals.progress.emit(int(percent))
        elif d['status'] == 'finished':
            self.signals.progress.emit(100)
            self.signals.log.emit("Finished downloading, processing...")

    def log(self, message):
        self.status.append(message)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    mw = MainWindow()
    mw.resize(500, 350)
    mw.show()
    sys.exit(app.exec())
