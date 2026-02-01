#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GenBank to mVISTA Annotation File Converter
Batch convert GenBank annotation files to mVISTA annotation file format

mVISTA format specification:
- Each gene line starts with > (plus strand) or < (minus strand)
- Followed by start position, end position, and gene name
- Exons are listed with start, end positions and "exon"
- UTRs are listed with start, end positions and "utr"
- Each gene block is separated by empty lines
"""

import os
import sys
from pathlib import Path
from typing import List
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QLineEdit, QListWidget, QTextEdit,
    QFileDialog, QMessageBox, QGroupBox, QProgressBar, QSplitter
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QFont
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature


class ConversionThread(QThread):
    """Background conversion thread"""
    progress_updated = pyqtSignal(int)
    file_completed = pyqtSignal(str, str)
    error_occurred = pyqtSignal(str, str)
    all_completed = pyqtSignal(int, int)
    log_message = pyqtSignal(str)

    def __init__(self, input_files: List[str], output_dir: str):
        super().__init__()
        self.input_files = input_files
        self.output_dir = output_dir
        self.success_count = 0
        self.fail_count = 0

    def run(self):
        """Execute conversion"""
        total_files = len(self.input_files)

        for i, input_file in enumerate(self.input_files):
            try:
                self.log_message.emit(f"Processing: {Path(input_file).name}")

                # Read GenBank file
                records = list(SeqIO.parse(input_file, "genbank"))

                if not records:
                    raise ValueError("No valid GenBank records found")

                # Convert each record
                for record in records:
                    output_file = self.convert_record(record, self.output_dir)
                    self.file_completed.emit(input_file, output_file)
                    self.log_message.emit(f"✓ Success: {output_file}")

                self.success_count += 1

            except Exception as e:
                self.fail_count += 1
                error_msg = f"Conversion failed: {str(e)}"
                self.error_occurred.emit(input_file, error_msg)
                self.log_message.emit(f"✗ {error_msg}")

            # Update progress
            progress = int((i + 1) / total_files * 100)
            self.progress_updated.emit(progress)

        self.all_completed.emit(self.success_count, self.fail_count)

    def convert_record(self, record: SeqRecord, output_dir: str) -> str:
        """
        Convert GenBank record to mVISTA format

        Args:
            record: GenBank record
            output_dir: Output directory

        Returns:
            Output file path
        """
        # Create output filename
        record_id = record.id.replace(" ", "_")
        output_filename = f"{record_id}_mvista.txt"
        output_path = os.path.join(output_dir, output_filename)

        # Extract features and group by gene
        features_by_gene = {}  # {gene_name: {'strand': +1/-1, 'gene_range': (start, end), 'exons': [], 'utrs': []}}

        for feature in record.features:
            feature_type = feature.type.lower()

            # Only process gene, exon, and UTR-related features
            if feature_type not in ['gene', 'exon', '5utr', '3utr', 'utr']:
                continue

            # Get position information
            if feature.location is None:
                continue

            # Convert to 1-based coordinates
            start = int(feature.location.start) + 1
            end = int(feature.location.end)

            # Get gene name as grouping identifier
            gene_name = self.get_gene_name(feature)

            # Create entry if gene doesn't exist yet
            if gene_name not in features_by_gene:
                features_by_gene[gene_name] = {
                    'gene_start': None,
                    'gene_end': None,
                    'gene_strand': 1,  # Default to plus strand
                    'exons': [],
                    'utrs': []
                }

            # Record gene start and end positions (from gene type feature)
            if feature_type == 'gene':
                features_by_gene[gene_name]['gene_start'] = start
                features_by_gene[gene_name]['gene_end'] = end
                # Get gene strand direction
                try:
                    features_by_gene[gene_name]['gene_strand'] = feature.strand if feature.strand else 1
                except (AttributeError, TypeError):
                    features_by_gene[gene_name]['gene_strand'] = 1

            # Record exon
            elif feature_type == 'exon':
                features_by_gene[gene_name]['exons'].append({
                    'start': start,
                    'end': end
                })

            # Record UTR (including 5'UTR and 3'UTR)
            elif feature_type in ['5utr', '3utr', 'utr']:
                features_by_gene[gene_name]['utrs'].append({
                    'start': start,
                    'end': end
                })

        # Write mVISTA format file
        with open(output_path, 'w', encoding='utf-8') as f:
            # Sort by gene start position
            sorted_genes = sorted(
                [(name, data) for name, data in features_by_gene.items()],
                key=lambda x: (x[1]['gene_start'] or float('inf'))
            )

            # Write each gene and its features
            for gene_name, gene_data in sorted_genes:
                # Get gene start and end positions
                gene_start = gene_data['gene_start']
                gene_end = gene_data['gene_end']
                gene_strand = gene_data['gene_strand']

                # If no gene type feature found, use min and max positions of all features
                if gene_start is None or gene_end is None:
                    all_coords = []
                    all_coords.extend([(e['start'], e['end']) for e in gene_data['exons']])
                    all_coords.extend([(u['start'], u['end']) for u in gene_data['utrs']])
                    if all_coords:
                        gene_start = min(all_coords)
                        gene_end = max([coord[1] for coord in all_coords])

                # Write gene line (use > for plus strand, < for minus strand)
                if gene_start and gene_end:
                    strand_symbol = '>' if gene_strand == 1 else '<'
                    f.write(f"{strand_symbol} {gene_start} {gene_end} {gene_name}\n")

                    # Write UTRs (sorted by position)
                    sorted_utrs = sorted(gene_data['utrs'], key=lambda x: x['start'])
                    for utr in sorted_utrs:
                        f.write(f"{utr['start']} {utr['end']} utr\n")

                    # Write exons (sorted by position)
                    sorted_exons = sorted(gene_data['exons'], key=lambda x: x['start'])
                    for exon in sorted_exons:
                        f.write(f"{exon['start']} {exon['end']} exon\n")

                # Add empty line between genes
                f.write("\n")

        return output_path

    def get_gene_name(self, feature: SeqFeature) -> str:
        """Get gene name (for grouping)"""
        # Try to get name from different fields
        if 'gene' in feature.qualifiers:
            return feature.qualifiers['gene'][0]
        elif 'locus_tag' in feature.qualifiers:
            return feature.qualifiers['locus_tag'][0]
        elif 'protein_id' in feature.qualifiers:
            return feature.qualifiers['protein_id'][0]
        elif 'label' in feature.qualifiers:
            return feature.qualifiers['label'][0]
        else:
            # If no name found, use feature ID or generate one
            if hasattr(feature, 'id') and feature.id:
                return str(feature.id)
            else:
                # Use position information as identifier
                if feature.location:
                    pos = int(feature.location.start) + 1
                    return f"feature_{pos}"
                return "unknown"


class GenBankToMVISTAConverter(QMainWindow):
    """Main window"""

    def __init__(self):
        super().__init__()
        self.input_files = []
        self.init_ui()
        self.conversion_thread = None

    def init_ui(self):
        """Initialize UI"""
        self.setWindowTitle("GenBank to mVISTA Converter")
        self.setGeometry(100, 100, 1000, 700)

        # Central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        # Main layout
        main_layout = QVBoxLayout(central_widget)

        # Title
        title_label = QLabel("GenBank → mVISTA Batch Converter")
        title_label.setFont(QFont("Arial", 16, QFont.Bold))
        title_label.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(title_label)

        # Create splitter
        splitter = QSplitter(Qt.Horizontal)

        # Left panel - File selection and settings
        left_panel = self.create_left_panel()
        splitter.addWidget(left_panel)

        # Right panel - Log display
        right_panel = self.create_right_panel()
        splitter.addWidget(right_panel)

        splitter.setSizes([500, 500])
        main_layout.addWidget(splitter)

        # Bottom progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(self.progress_bar)

    def create_left_panel(self) -> QWidget:
        """Create left control panel"""
        panel = QWidget()
        layout = QVBoxLayout(panel)

        # Input file selection group
        input_group = QGroupBox("Input Files (GenBank format)")
        input_layout = QVBoxLayout()

        # File list
        self.file_list = QListWidget()
        input_layout.addWidget(self.file_list)

        # Button layout
        btn_layout = QHBoxLayout()

        add_btn = QPushButton("Add Files")
        add_btn.clicked.connect(self.add_files)
        btn_layout.addWidget(add_btn)

        add_dir_btn = QPushButton("Add Folder")
        add_dir_btn.clicked.connect(self.add_directory)
        btn_layout.addWidget(add_dir_btn)

        clear_btn = QPushButton("Clear List")
        clear_btn.clicked.connect(self.clear_files)
        btn_layout.addWidget(clear_btn)

        input_layout.addLayout(btn_layout)

        # File count
        self.file_count_label = QLabel("0 files selected")
        input_layout.addWidget(self.file_count_label)

        input_group.setLayout(input_layout)
        layout.addWidget(input_group)

        # Output directory selection group
        output_group = QGroupBox("Output Settings")
        output_layout = QVBoxLayout()

        dir_layout = QHBoxLayout()
        dir_layout.addWidget(QLabel("Output directory:"))
        self.output_dir_edit = QLineEdit()
        self.output_dir_edit.setPlaceholderText("Select output directory...")
        dir_layout.addWidget(self.output_dir_edit)

        browse_btn = QPushButton("Browse...")
        browse_btn.clicked.connect(self.browse_output_dir)
        dir_layout.addWidget(browse_btn)

        output_layout.addLayout(dir_layout)
        output_group.setLayout(output_layout)
        layout.addWidget(output_group)

        # Feature type selection group (fixed to gene, exon, and UTR)
        feature_group = QGroupBox("Feature Type Info")
        feature_layout = QVBoxLayout()

        # Add info text
        info_label = QLabel("mVISTA format only includes these feature types:\n"
                            "• gene: Gene (with strand marker >/<)\n"
                            "• exon: Exon\n"
                            "• utr: Untranslated region (5'UTR and 3'UTR)")
        info_label.setWordWrap(True)
        feature_layout.addWidget(info_label)

        feature_group.setLayout(feature_layout)
        layout.addWidget(feature_group)

        # Convert button
        convert_btn = QPushButton("Start Conversion")
        convert_btn.setFont(QFont("Arial", 12, QFont.Bold))
        convert_btn.clicked.connect(self.start_conversion)
        layout.addWidget(convert_btn)

        return panel

    def create_right_panel(self) -> QWidget:
        """Create right log panel"""
        panel = QWidget()
        layout = QVBoxLayout(panel)

        # Log display
        log_label = QLabel("Conversion Log")
        log_label.setFont(QFont("Arial", 12, QFont.Bold))
        layout.addWidget(log_label)

        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setFont(QFont("Consolas", 9))
        layout.addWidget(self.log_text)

        # Clear log button
        clear_log_btn = QPushButton("Clear Log")
        clear_log_btn.clicked.connect(self.clear_log)
        layout.addWidget(clear_log_btn)

        return panel

    def add_files(self):
        """Add files"""
        files, _ = QFileDialog.getOpenFileNames(
            self, "Select GenBank Files", "",
            "GenBank files (*.gb *.gbk *.gbff);;All files (*.*)"
        )
        if files:
            for file in files:
                if file not in self.input_files:
                    self.input_files.append(file)
                    self.file_list.addItem(Path(file).name)
            self.update_file_count()

    def add_directory(self):
        """Add folder"""
        directory = QFileDialog.getExistingDirectory(
            self, "Select folder containing GenBank files"
        )
        if directory:
            for file in Path(directory).glob("*.gb*"):
                file_path = str(file)
                if file_path not in self.input_files:
                    self.input_files.append(file_path)
                    self.file_list.addItem(file.name)
            self.update_file_count()

    def clear_files(self):
        """Clear file list"""
        self.input_files.clear()
        self.file_list.clear()
        self.update_file_count()

    def update_file_count(self):
        """Update file count"""
        self.file_count_label.setText(f"{len(self.input_files)} files selected")

    def browse_output_dir(self):
        """Browse output directory"""
        directory = QFileDialog.getExistingDirectory(
            self, "Select output directory"
        )
        if directory:
            self.output_dir_edit.setText(directory)

    def start_conversion(self):
        """Start conversion"""
        # Validate input
        if not self.input_files:
            QMessageBox.warning(self, "Warning", "Please select input files first!")
            return

        output_dir = self.output_dir_edit.text()
        if not output_dir:
            QMessageBox.warning(self, "Warning", "Please select output directory!")
            return

        if not os.path.exists(output_dir):
            QMessageBox.warning(self, "Warning", "Output directory does not exist!")
            return

        # Clear log
        self.log_text.clear()
        self.log_message("=== Starting conversion ===")
        self.log_message(f"Input files: {len(self.input_files)}")
        self.log_message(f"Output directory: {output_dir}")
        self.log_message("Feature types: gene, exon, utr (fixed)")
        self.log_message("")

        # Reset progress bar
        self.progress_bar.setValue(0)

        # Create and start conversion thread
        self.conversion_thread = ConversionThread(
            self.input_files, output_dir
        )
        self.conversion_thread.progress_updated.connect(self.progress_bar.setValue)
        self.conversion_thread.file_completed.connect(self.on_file_completed)
        self.conversion_thread.error_occurred.connect(self.on_error_occurred)
        self.conversion_thread.all_completed.connect(self.on_all_completed)
        self.conversion_thread.log_message.connect(self.log_message)
        self.conversion_thread.start()

        # Disable buttons
        self.set_buttons_enabled(False)

    def on_file_completed(self, input_file: str, output_file: str):
        """File conversion complete callback"""
        self.log_message(f"Completed: {Path(input_file).name} → {Path(output_file).name}")

    def on_error_occurred(self, input_file: str, error_msg: str):
        """Error occurred callback"""
        self.log_message(f"Error: {Path(input_file).name}\n  {error_msg}")

    def on_all_completed(self, success_count: int, fail_count: int):
        """All files conversion complete callback"""
        self.log_message("")
        self.log_message("=== Conversion completed ===")
        self.log_message(f"Success: {success_count} files")
        self.log_message(f"Failed: {fail_count} files")
        self.log_message(f"Total: {success_count + fail_count} files")

        # Show result dialog
        QMessageBox.information(
            self, "Conversion Complete",
            f"Conversion complete!\n\nSuccess: {success_count} files\nFailed: {fail_count} files"
        )

        # Enable buttons
        self.set_buttons_enabled(True)
        self.progress_bar.setValue(100)

    def log_message(self, message: str):
        """Add log message"""
        self.log_text.append(message)
        # Auto scroll to bottom
        scrollbar = self.log_text.verticalScrollBar()
        scrollbar.setValue(scrollbar.maximum())

    def clear_log(self):
        """Clear log"""
        self.log_text.clear()

    def set_buttons_enabled(self, enabled: bool):
        """Set button enabled state"""
        self.file_list.setEnabled(enabled)
        self.output_dir_edit.setEnabled(enabled)


def main():
    """Main function"""
    app = QApplication(sys.argv)

    # Set application style
    app.setStyle('Fusion')

    # Create and show main window
    window = GenBankToMVISTAConverter()
    window.show()

    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
