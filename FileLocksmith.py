import os
import re
import subprocess
import tkinter as tk
from tkinter import messagebox
import tkinter.font as tkFont
from ttkbootstrap import Style
from ttkbootstrap.widgets import Frame, Label, Combobox, Button, Treeview, Scrollbar

## pip install ttkbootstrap
## sudo apt install policykit-1-gnome
## sudo apt install polkit-kde-agent-1

class DriveProcessViewer:
    def __init__(self, root):
        self.root = root
        self.root.title("Drive Process Manager")
        self.root.geometry("1000x600")

        self.style = Style("flatly")  # Options: flatly, darkly, etc.
        self.root.configure(bg=self.style.colors.bg)

        # Set default font to Arial
        default_font = tkFont.nametofont("TkDefaultFont")
        default_font.configure(family="Arial", size=10)
        self.root.option_add("*Font", default_font)

        # Layout
        top = Frame(root, padding=10)
        top.pack(fill=tk.X)

        middle = Frame(root, padding=10)
        middle.pack(fill=tk.BOTH, expand=True)

        bottom = Frame(root, padding=10)
        bottom.pack(fill=tk.X)

        Label(top, text="Select Drive:", font=("Arial", 11, "bold")).pack(side=tk.LEFT, padx=(0, 10))

        self.drive_var = tk.StringVar()
        self.drive_combo = Combobox(top, textvariable=self.drive_var, width=60)
        self.drive_combo.pack(side=tk.LEFT, padx=(0, 10))

        Button(top, text="↻ Refresh Drives", command=self.populate_drives).pack(side=tk.LEFT, padx=5)
        Button(top, text="Scan Processes", bootstyle="primary", command=self.scan_processes).pack(side=tk.LEFT, padx=5)

        # Treeview
        tree_frame = Frame(middle)
        tree_frame.pack(fill=tk.BOTH, expand=True)

        self.tree = Treeview(tree_frame, columns=("PID", "User", "Access", "Command"), show="headings", bootstyle="info")
        for col in self.tree["columns"]:
            self.tree.heading(col, text=col)
        self.tree.column("PID", width=80)
        self.tree.column("User", width=120)
        self.tree.column("Access", width=100)
        self.tree.column("Command", width=500)

        vsb = Scrollbar(tree_frame, orient="vertical", command=self.tree.yview)
        hsb = Scrollbar(tree_frame, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vsb.pack(side=tk.RIGHT, fill=tk.Y)
        hsb.pack(side=tk.BOTTOM, fill=tk.X)

        # Status & Buttons
        self.status_var = tk.StringVar()
        Label(bottom, textvariable=self.status_var, anchor=tk.W, bootstyle="secondary").pack(fill=tk.X, side=tk.LEFT)

        self.kill_btn = Button(bottom, text="Kill Selected Process", command=self.kill_selected_process, state=tk.DISABLED)
        self.kill_btn.pack(side=tk.RIGHT, padx=5)

        self.unmount_btn = Button(bottom, text="Unmount Drive", command=self.unmount_drive, state=tk.DISABLED, bootstyle="danger")
        self.unmount_btn.pack(side=tk.RIGHT, padx=5)

        self.tree.bind("<<TreeviewSelect>>", self.on_tree_select)

        self.populate_drives()
        self.update_status("Ready")

    def run_as_root(self, command):
        """Run sensitive commands using pkexec or fallback to sudo"""
        try:
            return subprocess.run(["pkexec"] + command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        except Exception:
            return subprocess.run(["sudo"] + command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

    def populate_drives(self):
        try:
            result = subprocess.run(["df", "-h", "--output=target"], stdout=subprocess.PIPE, text=True, check=True)
            drives = [line.strip() for line in result.stdout.splitlines()[1:] if line.startswith("/media/") or line.startswith("/mnt/")]
            self.drive_combo["values"] = drives
            if drives:
                self.drive_var.set(drives[0])
                self.update_status(f"Detected {len(drives)} mounted drives")
            else:
                self.update_status("No external drives detected")
        except Exception as e:
            self.update_status(f"Drive error: {e}")

    def scan_processes(self):
        drive = self.drive_var.get()
        if not drive:
            messagebox.showwarning("No Drive", "Select a drive first")
            return
        self.tree.delete(*self.tree.get_children())
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
                    self.tree.insert("", "end", values=(pid, user, access, command))

            self.update_status(f"{len(self.tree.get_children())} process(es) found")
            self.unmount_btn["state"] = tk.NORMAL

        except subprocess.CalledProcessError as e:
            self.update_status(f"Error: {e.stderr.strip()}")

    def on_tree_select(self, event):
        if self.tree.selection():
            self.kill_btn["state"] = tk.NORMAL
        else:
            self.kill_btn["state"] = tk.DISABLED

    def kill_selected_process(self):
        selected = self.tree.selection()
        if not selected:
            return
        item = self.tree.item(selected[0])
        pid = str(item["values"][0])
        if messagebox.askyesno("Confirm Kill", f"Kill process {pid}?"):
            try:
                self.run_as_root(["kill", "-9", pid])
                self.tree.delete(selected[0])
                self.update_status(f"Process {pid} terminated")
                self.kill_btn["state"] = tk.DISABLED
            except subprocess.CalledProcessError:
                self.update_status(f"Failed to kill process {pid}")

    def unmount_drive(self):
        drive = self.drive_var.get()
        if not drive:
            return
        try:
            self.run_as_root(["umount", drive])
            self.update_status(f"Unmounted {drive}")
            self.populate_drives()
            self.unmount_btn["state"] = tk.DISABLED
        except subprocess.CalledProcessError as e:
            self.update_status(f"Unmount failed: {e.stderr.strip()}")

    def update_status(self, message):
        self.status_var.set(message)
        self.root.update_idletasks()


if __name__ == "__main__":
    root = tk.Tk()
    app = DriveProcessViewer(root)
    root.mainloop()
