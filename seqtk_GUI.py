import tkinter as tk
from tkinter import ttk, filedialog, scrolledtext, messagebox
from Bio import SeqIO

class GeneExtractorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("extract subsequences from FASTA/Q files")
        self.root.geometry("750x650")
        self.root.resizable(True, True)
        
        # Configure style
        self.style = ttk.Style()
        self.style.theme_use('clam')
        self.configure_styles()
        
        # Create main frames
        self.create_widgets()
        
    def configure_styles(self):
        self.root.configure(bg='#2e2e2e')
        self.style.configure('TFrame', background='#2e2e2e')
        self.style.configure('TLabel', background='#2e2e2e', foreground='white')
        self.style.configure('TButton', background='#3c3f41', foreground='white')
        self.style.configure('TEntry', fieldbackground='#3c3f41', foreground='white')
        self.style.configure('TCombobox', fieldbackground='#3c3f41', foreground='white')
        self.style.configure('Header.TLabel', font=('Arial', 12, 'bold'))
        self.style.map('TButton', background=[('active', '#4e5254')])
        
    def create_widgets(self):
        # Main container
        main_frame = ttk.Frame(self.root, padding=10)
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Gene IDs section
        ttk.Label(main_frame, text="Gene IDs to Extract", style='Header.TLabel').grid(row=0, column=0, sticky=tk.W, pady=(0, 5))
        
        self.id_text = scrolledtext.ScrolledText(main_frame, width=85, height=15, 
                                               bg='#3c3f41', fg='white', insertbackground='white')
        self.id_text.grid(row=1, column=0, columnspan=3, pady=(0, 10))
        
        # File selection section
        file_frame = ttk.Frame(main_frame)
        file_frame.grid(row=2, column=0, columnspan=3, sticky=tk.EW, pady=10)
        
        # Input FASTA
        ttk.Label(file_frame, text="Input FASTA File:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5))
        self.fasta_entry = ttk.Entry(file_frame, width=50)
        self.fasta_entry.grid(row=0, column=1, sticky=tk.EW, padx=5)
        ttk.Button(file_frame, text="Browse...", command=self.browse_fasta).grid(row=0, column=2)
        
        # Output File
        ttk.Label(file_frame, text="Output File:").grid(row=1, column=0, sticky=tk.W, padx=(0, 5), pady=(10, 0))
        self.output_entry = ttk.Entry(file_frame, width=50)
        self.output_entry.grid(row=1, column=1, sticky=tk.EW, padx=5, pady=(10, 0))
        ttk.Button(file_frame, text="Browse...", command=self.browse_output).grid(row=1, column=2, pady=(10, 0))
        
        # Line length option
        ttk.Label(file_frame, text="Line Length:").grid(row=2, column=0, sticky=tk.W, padx=(0, 5), pady=(10, 0))
        self.line_length = ttk.Combobox(file_frame, width=10, values=[60, 70, 80, 100, 120])
        self.line_length.set(60)
        self.line_length.grid(row=2, column=1, sticky=tk.W, padx=5, pady=(10, 0))
        
        # Results section
        results_frame = ttk.Frame(main_frame)
        results_frame.grid(row=3, column=0, columnspan=3, sticky=tk.EW, pady=10)
        
        self.results_var = tk.StringVar()
        self.results_var.set("Ready to extract sequences")
        ttk.Label(results_frame, textvariable=self.results_var).grid(row=0, column=0, sticky=tk.W)
        
        self.missing_ids_text = scrolledtext.ScrolledText(results_frame, width=85, height=8, 
                                                         bg='#3c3f41', fg='white', state=tk.DISABLED)
        self.missing_ids_text.grid(row=1, column=0, sticky=tk.EW, pady=(5, 0))
        
        # Action buttons
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=4, column=0, columnspan=3, pady=10)
        
        ttk.Button(button_frame, text="Extract Sequences", command=self.extract_sequences).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Clear All", command=self.clear_all).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Exit", command=self.root.destroy).pack(side=tk.RIGHT, padx=5)
        
        # Configure grid weights
        main_frame.columnconfigure(0, weight=1)
        file_frame.columnconfigure(1, weight=1)
        
    def browse_fasta(self):
        file_path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")])
        if file_path:
            self.fasta_entry.delete(0, tk.END)
            self.fasta_entry.insert(0, file_path)
            
    def browse_output(self):
        file_path = filedialog.asksaveasfilename(
            defaultextension=".fasta",
            filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")]
        )
        if file_path:
            self.output_entry.delete(0, tk.END)
            self.output_entry.insert(0, file_path)
            
    def read_fasta_file(self, file_path):
        try:
            return {record.id: str(record.seq) for record in SeqIO.parse(file_path, "fasta")}
        except Exception as e:
            messagebox.showerror("Error", f"Error reading FASTA file:\n{str(e)}")
            return {}
            
    def extract_genes_from_fasta(self, gene_ids, fasta_file):
        fasta_sequences = self.read_fasta_file(fasta_file)
        return {gene_id: fasta_sequences[gene_id] for gene_id in gene_ids if gene_id in fasta_sequences}
            
    def write_extracted_sequences_to_file(self, extracted_sequences, output_file, line_length=60):
        try:
            with open(output_file, 'w') as out_file:
                for gene_id, sequence in extracted_sequences.items():
                    out_file.write(f'>{gene_id}\n')
                    for i in range(0, len(sequence), line_length):
                        out_file.write(sequence[i:i+line_length] + '\n')
            return True
        except Exception as e:
            messagebox.showerror("Error", f"Error writing output file:\n{str(e)}")
            return False
            
    def extract_sequences(self):
        # Get input values
        gene_ids = {line.strip() for line in self.id_text.get("1.0", tk.END).splitlines() if line.strip()}
        fasta_file = self.fasta_entry.get()
        output_file = self.output_entry.get()
        
        # Validate inputs
        if not gene_ids:
            messagebox.showwarning("Input Error", "Please enter at least one Gene ID")
            return
        if not fasta_file:
            messagebox.showwarning("Input Error", "Please select an input FASTA file")
            return
        if not output_file:
            messagebox.showwarning("Input Error", "Please specify an output file")
            return
        
        try:
            line_length = int(self.line_length.get())
        except ValueError:
            messagebox.showwarning("Input Error", "Line length must be an integer")
            return
            
        # Perform extraction
        extracted_genes = self.extract_genes_from_fasta(gene_ids, fasta_file)
        
        if not extracted_genes:
            messagebox.showinfo("No Results", "No matching sequences found")
            return
            
        # Write results
        if self.write_extracted_sequences_to_file(extracted_genes, output_file, line_length):
            # Show results
            total_requested = len(gene_ids)
            extracted_count = len(extracted_genes)
            missing = gene_ids - extracted_genes.keys()
            
            self.results_var.set(
                f"Extraction complete! Requested: {total_requested}, "
                f"Extracted: {extracted_count}, Missing: {len(missing)}"
            )
            
            # Display missing IDs
            self.missing_ids_text.config(state=tk.NORMAL)
            self.missing_ids_text.delete(1.0, tk.END)
            if missing:
                self.missing_ids_text.insert(tk.END, "Missing IDs:\n" + "\n".join(sorted(missing)))
            else:
                self.missing_ids_text.insert(tk.END, "All requested IDs were found and extracted.")
            self.missing_ids_text.config(state=tk.DISABLED)
            
            messagebox.showinfo("Success", f"Successfully extracted {extracted_count} sequences to:\n{output_file}")
            
    def clear_all(self):
        self.id_text.delete(1.0, tk.END)
        self.fasta_entry.delete(0, tk.END)
        self.output_entry.delete(0, tk.END)
        self.line_length.set(60)
        self.results_var.set("Ready to extract sequences")
        self.missing_ids_text.config(state=tk.NORMAL)
        self.missing_ids_text.delete(1.0, tk.END)
        self.missing_ids_text.config(state=tk.DISABLED)

if __name__ == "__main__":
    root = tk.Tk()
    app = GeneExtractorApp(root)
    root.mainloop()