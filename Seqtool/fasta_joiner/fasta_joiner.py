#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
FASTA Sequence Joiner Tool
Concatenates all sequences from a FASTA file into one sequence
"""

import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext
from Bio import SeqIO


class FastaJoinerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("FASTA Sequence Joiner Tool")
        self.root.geometry("600x500")
        self.root.resizable(True, True)
        
        # Store input and output file paths
        self.input_file = ""
        self.output_file = ""
        
        # Create interface
        self.create_widgets()
        
    def create_widgets(self):
        """Create interface components"""
        # Title
        title_label = tk.Label(
            self.root,
            text="FASTA Sequence Joiner Tool",
            font=("Arial", 16, "bold"),
            pady=10
        )
        title_label.pack()
        
        # Description text
        info_frame = tk.Frame(self.root, padx=20)
        info_frame.pack(fill=tk.X, pady=5)
        
        info_text = """
This tool concatenates all sequences in a FASTA file into one sequence.
Suitable for merging multiple sequence fragments.
        """
        info_label = tk.Label(info_frame, text=info_text, justify=tk.LEFT, fg="gray")
        info_label.pack()
        
        # Step 1: Select input file
        step1_frame = tk.LabelFrame(
            self.root,
            text="Step 1: Select FASTA Input File",
            font=("Arial", 10, "bold"),
            padx=10, pady=10
        )
        step1_frame.pack(fill=tk.X, padx=20, pady=10)
        
        self.input_file_label = tk.Label(step1_frame, text="No file selected", fg="blue", cursor="hand2")
        self.input_file_label.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.input_file_label.bind("<Button-1>", self.browse_input_file)
        
        browse_btn = tk.Button(
            step1_frame,
            text="Browse File",
            command=self.browse_input_file
        )
        browse_btn.pack(side=tk.RIGHT)
        
        # Step 2: Set output options
        step2_frame = tk.LabelFrame(
            self.root,
            text="Step 2: Set Output Options",
            font=("Arial", 10, "bold"),
            padx=10, pady=10
        )
        step2_frame.pack(fill=tk.X, padx=20, pady=10)
        
        # Sequence name input
        name_frame = tk.Frame(step2_frame)
        name_frame.pack(fill=tk.X, pady=5)
        
        name_label = tk.Label(name_frame, text="New Sequence Name:", font=("Arial", 9))
        name_label.pack(side=tk.LEFT, padx=(0, 5))
        
        self.seq_name_entry = tk.Entry(
            name_frame,
            font=("Arial", 9),
            width=40
        )
        self.seq_name_entry.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.seq_name_entry.insert(0, "joined_sequence")  # Default value
        self.seq_name_entry.config(state="readonly")  # Default readonly
        
        # Add custom name button
        custom_name_btn = tk.Button(
            name_frame,
            text="Custom",
            command=self.enable_custom_name,
            font=("Arial", 8),
            padx=5
        )
        custom_name_btn.pack(side=tk.LEFT, padx=(5, 0))
        
        # File path selection
        path_frame = tk.Frame(step2_frame)
        path_frame.pack(fill=tk.X, pady=5)
        
        self.output_file_label = tk.Label(path_frame, text="No path selected", fg="blue", cursor="hand2")
        self.output_file_label.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.output_file_label.bind("<Button-1>", self.browse_output_file)
        
        output_btn = tk.Button(
            path_frame,
            text="Save Location",
            command=self.browse_output_file
        )
        output_btn.pack(side=tk.RIGHT)
        
        # Step 3: Execute
        step3_frame = tk.LabelFrame(
            self.root,
            text="Step 3: Execute Join Operation",
            font=("Arial", 10, "bold"),
            padx=10, pady=10
        )
        step3_frame.pack(fill=tk.X, padx=20, pady=10)
        
        run_btn = tk.Button(
            step3_frame,
            text="Start Joining Sequences",
            command=self.join_sequences,
            bg="#4CAF50",
            fg="white",
            font=("Arial", 11, "bold"),
            padx=20,
            pady=5
        )
        run_btn.pack()
        
        # Result display area
        result_frame = tk.LabelFrame(
            self.root,
            text="Processing Results",
            font=("Arial", 10, "bold"),
            padx=10, pady=10
        )
        result_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=10)
        
        self.result_text = scrolledtext.ScrolledText(
            result_frame,
            height=8,
            wrap=tk.WORD
        )
        self.result_text.pack(fill=tk.BOTH, expand=True)
        
        # Help button
        help_btn = tk.Button(
            self.root,
            text="Help",
            command=self.show_help
        )
        help_btn.pack(pady=5)
        
    def enable_custom_name(self):
        """Enable custom sequence name"""
        current_state = self.seq_name_entry.cget("state")
        if current_state == "readonly":
            self.seq_name_entry.config(state="normal")
            self.seq_name_entry.focus()
    
    def browse_input_file(self, event=None):
        """Browse input file"""
        file_path = filedialog.askopenfilename(
            title="Select FASTA File",
            filetypes=[("FASTA files", "*.fasta *.fa *.fna *.ffn *.faa *.frn"), ("All files", "*.*")]
        )
        
        if file_path:
            self.input_file = file_path
            from os.path import basename
            self.input_file_label.config(text=f"Selected: {basename(file_path)}")
            self.log(f"✓ Input file selected: {file_path}")
    
    def browse_output_file(self, event=None):
        """Browse output file"""
        file_path = filedialog.asksaveasfilename(
            title="Save Output File",
            defaultextension=".fasta",
            filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")],
            initialfile="joined_sequence.fasta"
        )
        
        if file_path:
            self.output_file = file_path
            from os.path import basename
            self.output_file_label.config(text=f"Save to: {basename(file_path)}")
            self.log(f"✓ Output file set: {file_path}")
    
    def join_sequences(self):
        """Execute sequence concatenation"""
        # Validate input
        if not self.input_file:
            messagebox.showwarning("Warning", "Please select a FASTA input file first!")
            return
        
        if not self.output_file:
            messagebox.showwarning("Warning", "Please select an output file path first!")
            return
        
        try:
            self.log("\n" + "="*60)
            self.log("Starting processing...")
            self.log("="*60 + "\n")
            
            # Read all sequences
            sequences = []
            seq_count = 0
            total_length = 0
            
            self.log(f"Reading file: {self.input_file}")
            
            for record in SeqIO.parse(self.input_file, "fasta"):
                sequences.append(record)
                seq_count += 1
                total_length += len(record.seq)
                self.log(f"  - Sequence {seq_count}: {record.id} (Length: {len(record.seq)})")
            
            if not sequences:
                self.log("❌ Error: No sequences found in file!")
                messagebox.showerror("Error", "No valid FASTA sequences found in the file!")
                return
            
            self.log(f"\nTotal {seq_count} sequences found")
            self.log(f"Total sequence length: {total_length} bp")
            
            # Concatenate all sequences
            self.log("\nConcatenating sequences...")
            joined_sequence = "".join(str(record.seq) for record in sequences)
            self.log(f"✓ Sequence concatenation complete")
            self.log(f"  Joined sequence length: {len(joined_sequence)} bp")
            
            # Get custom sequence name
            custom_name = self.seq_name_entry.get().strip()
            if not custom_name:
                custom_name = "joined_sequence"
            self.log(f"  New sequence name: {custom_name}")
            
            # Create output record (using new SeqRecord object)
            from Bio.SeqRecord import SeqRecord
            from Bio.Seq import Seq
            output_record = SeqRecord(
                Seq(joined_sequence),
                id=custom_name,
                description=f"joined_from_{seq_count}_sequences_total_length_{len(joined_sequence)}"
            )
            
            # Write output file
            self.log(f"\nWriting output file: {self.output_file}")
            SeqIO.write(output_record, self.output_file, "fasta")
            self.log("✓ Output file saved successfully")
            
            # Show success message
            from os.path import basename
            success_msg = f"""
Processing completed successfully!

Input file: {basename(self.input_file)}
Number of sequences: {seq_count}
Total sequence length: {total_length} bp

Joined sequence length: {len(joined_sequence)} bp
Output file: {basename(self.output_file)}
            """
            
            self.log("\n" + "="*60)
            self.log("Processing successful!")
            self.log("="*60)
            messagebox.showinfo("Success", success_msg.strip())
            
        except Exception as e:
            error_msg = f"Error during processing:\n{str(e)}"
            self.log(f"\n❌ {error_msg}")
            messagebox.showerror("Error", error_msg)
    
    def log(self, message):
        """Display log in result area"""
        self.result_text.insert(tk.END, message + "\n")
        self.result_text.see(tk.END)
        self.result_text.update()
    
    def show_help(self):
        """Show help information"""
        help_text = """
FASTA Sequence Joiner Tool - User Guide

Function Description:
  This tool concatenates all sequence fragments in a FASTA file into one continuous sequence.
  Suitable for merging multiple sequence fragments into a single sequence.

Usage Steps:
  1. Click the "Browse File" button to select the FASTA input file to process
  2. Click the "Save Location" button to select the output file save path
  3. Click the "Start Joining Sequences" button to execute the concatenation operation

Supported Formats:
  Input: .fasta, .fa, .fna, .ffn, .faa, .frn and other FASTA format files
  Output: FASTA format file

Notes:
  - Ensure the input file is a valid FASTA format
  - The joined sequence ID will use the first sequence's ID
  - Sequences are concatenated in the order they appear in the input file
  - Sequences are directly connected without separators

Example:
  Input file contains 3 sequences:
    Seq1: ATGC
    Seq2: GCTA
    Seq3: CGAT
  
  Output result:
    Seq1: ATGCGCTACGAT (9 bp)
        """
        messagebox.showinfo("Help", help_text)


def main():
    """Main function"""
    root = tk.Tk()
    app = FastaJoinerApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
