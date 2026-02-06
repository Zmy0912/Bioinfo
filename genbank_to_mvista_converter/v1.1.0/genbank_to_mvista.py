#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GenBank/GFF3 to mVISTA Annotation File Converter
Batch convert GenBank or GFF3 format annotation files to mVISTA annotation file format

mVISTA Format Specification:
- Tab-separated text file
- Fields: Feature, Start, End, Name, Direction, Annotation
- Feature: feature type (gene, exon, CDS, rRNA, tRNA, mRNA, etc.)
- Start: start position (1-based)
- End: end position (1-based)
- Name: gene/feature name
- Direction: direction (+/- or forward/reverse)
- Annotation: annotation information

GFF3 Format Specification:
- Column 1: seqid (sequence ID)
- Column 2: source (origin)
- Column 3: type (feature type)
- Column 4: start (start position, 1-based)
- Column 5: end (end position, 1-based)
- Column 6: score (score)
- Column 7: strand (chain direction: +, -, .)
- Column 8: phase (phase: 0, 1, 2, .)
- Column 9: attributes (attributes, semicolon-separated key-value pairs)
"""

import os
import sys
from pathlib import Path
from typing import List, Dict
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QLineEdit, QListWidget, QTextEdit,
    QFileDialog, QMessageBox, QGroupBox, QProgressBar, QSplitter,
    QRadioButton, QButtonGroup
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

    def __init__(self, input_files: List[str], output_dir: str,
                 feature_types: List[str], file_format: str = 'genbank'):
        super().__init__()
        self.input_files = input_files
        self.output_dir = output_dir
        self.feature_types = feature_types
        self.file_format = file_format
        self.success_count = 0
        self.fail_count = 0

    def run(self):
        """Execute conversion"""
        total_files = len(self.input_files)
        
        for i, input_file in enumerate(self.input_files):
            try:
                self.log_message.emit(f"Processing: {Path(input_file).name}")
                
                if self.file_format == 'gff3':
                    # Read GFF3 file
                    self.convert_gff3_file(input_file)
                else:
                    # Read GenBank file
                    records = list(SeqIO.parse(input_file, "genbank"))
                    
                    if not records:
                        raise ValueError("No valid GenBank records found")
                    
                    # Convert each record
                    for record in records:
                        output_file = self.convert_record(record, self.output_dir, 
                                                         self.feature_types)
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

    def convert_record(self, record: SeqRecord, output_dir: str, 
                      feature_types: List[str]) -> str:
        """
        将GenBank记录转换为mVISTA格式
        
        Args:
            record: GenBank记录
            output_dir: 输出目录
            feature_types: 要提取的特征类型列表
            
        Returns:
            输出文件路径
        """
        # 创建输出文件名
        record_id = record.id.replace(" ", "_")
        output_filename = f"{record_id}_mvista.txt"
        output_path = os.path.join(output_dir, output_filename)
        
        # 提取特征并按基因分组
        features_by_gene = {}  # {gene_name: {'strand': +1/-1, 'gene_range': (start, end), 'exons': [], 'utrs': []}}
        
        for feature in record.features:
            feature_type = feature.type.lower()

            # 处理所有相关特征类型
            # RNA相关特征
            if feature_type in ['gene', 'exon', '5utr', '3utr', 'utr', 'rna', 'rrna', 'trna', 'snrna', 'scrna', 'mrna', 'ncrna', 'mirna', 'lnc_rna']:
                # 获取位置信息
                if feature.location is None:
                    continue

                # 转换为1-based坐标
                start = int(feature.location.start) + 1
                end = int(feature.location.end)

                # 获取基因名称作为分组标识
                gene_name = self.get_gene_name(feature)

                # RNA特征单独处理（作为独立基因处理）
                if feature_type in ['rna', 'rrna', 'trna', 'snrna', 'scrna', 'mrna', 'ncrna', 'mirna', 'lnc_rna']:
                    # 为RNA特征创建独立的条目
                    if gene_name not in features_by_gene:
                        features_by_gene[gene_name] = {
                            'gene_start': start,
                            'gene_end': end,
                            'gene_strand': 1,  # 默认为正链
                            'exons': [],
                            'utrs': [],
                            'is_rna': True,
                            'rna_type': feature_type
                        }
                    else:
                        # 更新RNA的范围
                        features_by_gene[gene_name]['gene_start'] = min(features_by_gene[gene_name]['gene_start'], start)
                        features_by_gene[gene_name]['gene_end'] = max(features_by_gene[gene_name]['gene_end'], end)

                    # 获取RNA的链方向
                    try:
                        features_by_gene[gene_name]['gene_strand'] = feature.strand if feature.strand else 1
                    except (AttributeError, TypeError):
                        features_by_gene[gene_name]['gene_strand'] = 1

                    continue

                # 普通基因相关特征的处理
                # 如果还没有这个基因的条目，创建一个
                if gene_name not in features_by_gene:
                    features_by_gene[gene_name] = {
                        'gene_start': None,
                        'gene_end': None,
                        'gene_strand': 1,  # 默认为正链
                        'exons': [],
                        'utrs': [],
                        'is_rna': False,
                        'rna_type': None
                    }

                # 记录基因的起始和结束位置（从gene类型特征中获取）
                if feature_type == 'gene':
                    features_by_gene[gene_name]['gene_start'] = start
                    features_by_gene[gene_name]['gene_end'] = end
                    # 获取基因的链方向
                    try:
                        features_by_gene[gene_name]['gene_strand'] = feature.strand if feature.strand else 1
                    except (AttributeError, TypeError):
                        features_by_gene[gene_name]['gene_strand'] = 1

                # 记录exon
                elif feature_type == 'exon':
                    features_by_gene[gene_name]['exons'].append({
                        'start': start,
                        'end': end
                    })

                # 记录UTR（包括5'UTR和3'UTR）
                elif feature_type in ['5utr', '3utr', 'utr']:
                    features_by_gene[gene_name]['utrs'].append({
                        'start': start,
                        'end': end
                    })
        
        # 写入mVISTA格式文件
        with open(output_path, 'w', encoding='utf-8') as f:
            # 按基因起始位置排序
            sorted_genes = sorted(
                [(name, data) for name, data in features_by_gene.items()],
                key=lambda x: (x[1]['gene_start'] or float('inf'))
            )
            
            # 写入每个基因及其特征
            for gene_name, gene_data in sorted_genes:
                # 获取基因的起始和结束位置
                gene_start = gene_data['gene_start']
                gene_end = gene_data['gene_end']
                gene_strand = gene_data['gene_strand']
                is_rna = gene_data.get('is_rna', False)
                rna_type = gene_data.get('rna_type', None)

                # 如果没有找到gene类型的特征，使用所有特征的最小和最大位置
                if gene_start is None or gene_end is None:
                    all_coords = []
                    all_coords.extend([(e['start'], e['end']) for e in gene_data['exons']])
                    all_coords.extend([(u['start'], u['end']) for u in gene_data['utrs']])
                    if all_coords:
                        gene_start = min(all_coords)
                        gene_end = max([coord[1] for coord in all_coords])

                # 写入基因行（使用 > 表示正链，< 表示负链）
                if gene_start and gene_end:
                    strand_symbol = '>' if gene_strand == 1 else '<'

                    # RNA特征直接写入
                    if is_rna and rna_type:
                        # 标准化RNA类型名称
                        rna_display = rna_type.replace('_', '').upper()
                        f.write(f"{strand_symbol} {gene_start} {gene_end} {gene_name} ({rna_display})\n")

                    # 普通基因特征
                    else:
                        f.write(f"{strand_symbol} {gene_start} {gene_end} {gene_name}\n")

                        # 写入UTR（按位置排序）
                        sorted_utrs = sorted(gene_data['utrs'], key=lambda x: x['start'])
                        for utr in sorted_utrs:
                            f.write(f"{utr['start']} {utr['end']} utr\n")

                        # 写入exon（按位置排序）
                        sorted_exons = sorted(gene_data['exons'], key=lambda x: x['start'])
                        for exon in sorted_exons:
                            f.write(f"{exon['start']} {exon['end']} exon\n")

                # 基因之间添加空行
                f.write("\n")
        
        return output_path

    def get_gene_name(self, feature: SeqFeature) -> str:
        """获取基因名称（用于分组）"""
        # 尝试从不同字段获取名称
        if 'gene' in feature.qualifiers:
            return feature.qualifiers['gene'][0]
        elif 'locus_tag' in feature.qualifiers:
            return feature.qualifiers['locus_tag'][0]
        elif 'protein_id' in feature.qualifiers:
            return feature.qualifiers['protein_id'][0]
        elif 'label' in feature.qualifiers:
            return feature.qualifiers['label'][0]
        else:
            # 如果没有找到名称，使用特征ID或生成一个
            if hasattr(feature, 'id') and feature.id:
                return str(feature.id)
            else:
                # 使用位置信息作为标识
                if feature.location:
                    pos = int(feature.location.start) + 1
                    return f"feature_{pos}"
                return "unknown"

    def get_feature_name(self, feature: SeqFeature) -> str:
        """获取特征名称"""
        # 尝试从不同字段获取名称
        if 'gene' in feature.qualifiers:
            return feature.qualifiers['gene'][0]
        elif 'locus_tag' in feature.qualifiers:
            return feature.qualifiers['locus_tag'][0]
        elif 'protein_id' in feature.qualifiers:
            return feature.qualifiers['protein_id'][0]
        elif 'label' in feature.qualifiers:
            return feature.qualifiers['label'][0]
        else:
            return "."

    def get_feature_annotation(self, feature: SeqFeature) -> str:
        """获取特征注释"""
        # 尝试从product字段获取注释
        if 'product' in feature.qualifiers:
            return feature.qualifiers['product'][0]
        elif 'note' in feature.qualifiers:
            return feature.qualifiers['note'][0]
        elif 'function' in feature.qualifiers:
            return feature.qualifiers['function'][0]
        else:
            return "."

    def convert_gff3_file(self, input_file: str):
        """
        将GFF3文件转换为mVISTA格式

        Args:
            input_file: GFF3文件路径
        """
        # 创建输出文件名
        input_path = Path(input_file)
        output_filename = f"{input_path.stem}_mvista.txt"
        output_path = os.path.join(self.output_dir, output_filename)

        # 读取GFF3文件
        features_by_gene = {}  # {gene_name: {'strand': '+/-', 'gene_range': (start, end), 'exons': [], 'utrs': []}}
        mrna_to_gene = {}  # 存储mRNA到gene的映射关系

        with open(input_file, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()

                # 跳过注释行和空行
                if line.startswith('#') or not line:
                    continue

                # 解析GFF3行
                parts = line.split('\t')
                if len(parts) < 9:
                    continue

                seqid = parts[0]
                source = parts[1]
                feature_type = parts[2].lower()
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attributes_str = parts[8]

                # 处理所有相关特征类型（包括RNA）
                if feature_type not in ['gene', 'mrna', 'rna', 'rrna', 'trna', 'snrna', 'scrna',
                                       'ncrna', 'mirna', 'lnc_rna', 'exon', 'cds',
                                       '5utr', '3utr', 'utr', 'five_prime_utr', 'three_prime_utr']:
                    continue

                # 解析属性
                attributes = self.parse_gff3_attributes(attributes_str)

                # 处理gene特征 - 建立mRNA到gene的映射
                if feature_type == 'gene':
                    gene_id = attributes.get('ID', [f'gene_{start}'])[0]
                    # 暂时用ID，稍后会用mRNA的Name覆盖
                    if gene_id not in features_by_gene:
                        features_by_gene[gene_id] = {
                            'gene_start': start,
                            'gene_end': end,
                            'gene_strand': strand if strand in ['+', '-'] else '+',
                            'exons': [],
                            'utrs': [],
                            'cds_regions': [],
                            'gene_name': attributes.get('Name', [gene_id])[0],
                            'is_rna': False,
                            'rna_type': None
                        }

                # 处理RNA相关特征（作为独立基因处理）
                elif feature_type in ['rna', 'rrna', 'trna', 'snrna', 'scrna', 'ncrna', 'mirna', 'lnc_rna']:
                    rna_id = attributes.get('ID', [f'{feature_type}_{start}'])[0]
                    rna_name = attributes.get('Name', [rna_id])[0]

                    if rna_id not in features_by_gene:
                        features_by_gene[rna_id] = {
                            'gene_start': start,
                            'gene_end': end,
                            'gene_strand': strand if strand in ['+', '-'] else '+',
                            'exons': [],
                            'utrs': [],
                            'cds_regions': [],
                            'gene_name': rna_name,
                            'is_rna': True,
                            'rna_type': feature_type
                        }
                    else:
                        # 更新RNA的范围
                        features_by_gene[rna_id]['gene_start'] = min(features_by_gene[rna_id]['gene_start'], start)
                        features_by_gene[rna_id]['gene_end'] = max(features_by_gene[rna_id]['gene_end'], end)

                    # 建立RNA到自身的映射（便于后续关联子特征）
                    mrna_to_gene[rna_id] = rna_id

                # 处理mRNA特征 - 提取基因名称并建立映射
                elif feature_type == 'mrna':
                    mrna_id = attributes.get('ID', [f'mrna_{start}'])[0]
                    parent_id = attributes.get('Parent', [''])[0]
                    gene_name = attributes.get('Name', [parent_id])[0]

                    # 如果mRNA有parent，更新该gene的名称
                    if parent_id and parent_id in features_by_gene:
                        features_by_gene[parent_id]['gene_name'] = gene_name
                        features_by_gene[parent_id]['gene_start'] = start
                        features_by_gene[parent_id]['gene_end'] = end
                        features_by_gene[parent_id]['gene_strand'] = strand if strand in ['+', '-'] else '+'
                        # 建立mRNA到gene的映射
                        mrna_to_gene[mrna_id] = parent_id
                    else:
                        # 如果mRNA没有parent，作为独立基因处理
                        if mrna_id not in features_by_gene:
                            features_by_gene[mrna_id] = {
                                'gene_start': start,
                                'gene_end': end,
                                'gene_strand': strand if strand in ['+', '-'] else '+',
                                'exons': [],
                                'utrs': [],
                                'cds_regions': [],
                                'gene_name': gene_name,
                                'is_rna': True,
                                'rna_type': 'mrna'
                            }
                        else:
                            # 更新mRNA的范围
                            features_by_gene[mrna_id]['gene_start'] = min(features_by_gene[mrna_id]['gene_start'], start)
                            features_by_gene[mrna_id]['gene_end'] = max(features_by_gene[mrna_id]['gene_end'], end)
                        # 建立mRNA到自身的映射
                        mrna_to_gene[mrna_id] = mrna_id

                # 处理exon
                elif feature_type == 'exon':
                    parent_id = attributes.get('Parent', [''])[0]
                    # 通过mRNA找到对应的gene
                    gene_id = mrna_to_gene.get(parent_id, '')
                    if gene_id and gene_id in features_by_gene:
                        features_by_gene[gene_id]['exons'].append({
                            'start': start,
                            'end': end
                        })
                    # 如果没有mRNA映射，尝试直接通过parent找gene
                    elif parent_id and parent_id in features_by_gene:
                        features_by_gene[parent_id]['exons'].append({
                            'start': start,
                            'end': end
                        })

                # 处理CDS（作为外显子的补充）
                elif feature_type == 'cds':
                    parent_id = attributes.get('Parent', [''])[0]
                    # 通过mRNA找到对应的gene
                    gene_id = mrna_to_gene.get(parent_id, '')
                    if gene_id and gene_id in features_by_gene:
                        features_by_gene[gene_id]['cds_regions'].append({
                            'start': start,
                            'end': end
                        })
                    # 如果没有mRNA映射，尝试直接通过parent找gene
                    elif parent_id and parent_id in features_by_gene:
                        features_by_gene[parent_id]['cds_regions'].append({
                            'start': start,
                            'end': end
                        })

                # 记录UTR（包括5'UTR和3'UTR）
                elif feature_type in ['5utr', '3utr', 'utr', 'five_prime_utr', 'three_prime_utr']:
                    parent_id = attributes.get('Parent', [''])[0]
                    gene_id = mrna_to_gene.get(parent_id, '')
                    if gene_id and gene_id in features_by_gene:
                        features_by_gene[gene_id]['utrs'].append({
                            'start': start,
                            'end': end
                        })
                    elif parent_id and parent_id in features_by_gene:
                        features_by_gene[parent_id]['utrs'].append({
                            'start': start,
                            'end': end
                        })

        # 如果没有exon但有cds，使用cds作为外显子
        for gene_id in features_by_gene:
            if not features_by_gene[gene_id]['exons'] and features_by_gene[gene_id]['cds_regions']:
                features_by_gene[gene_id]['exons'] = features_by_gene[gene_id]['cds_regions']

        # 写入mVISTA格式文件
        with open(output_path, 'w', encoding='utf-8') as f:
            # 按基因起始位置排序
            sorted_genes = sorted(
                [(gene_id, data) for gene_id, data in features_by_gene.items()],
                key=lambda x: (x[1]['gene_start'] or float('inf'))
            )

            # 写入每个基因及其特征
            for gene_id, gene_data in sorted_genes:
                # 获取基因的名称（优先使用mRNA的Name）
                gene_name = gene_data['gene_name']
                gene_start = gene_data['gene_start']
                gene_end = gene_data['gene_end']
                gene_strand = gene_data['gene_strand']
                is_rna = gene_data.get('is_rna', False)
                rna_type = gene_data.get('rna_type', None)

                # 如果没有找到gene类型的特征，使用所有特征的最小和最大位置
                if gene_start is None or gene_end is None:
                    all_coords = []
                    all_coords.extend([(e['start'], e['end']) for e in gene_data['exons']])
                    all_coords.extend([(u['start'], u['end']) for u in gene_data['utrs']])
                    all_coords.extend([(c['start'], c['end']) for c in gene_data['cds_regions']])
                    if all_coords:
                        gene_start = min(all_coords)
                        gene_end = max([coord[1] for coord in all_coords])

                # 写入基因行（使用 > 表示正链，< 表示负链）
                if gene_start and gene_end:
                    strand_symbol = '>' if gene_strand == '+' else '<'

                    # RNA特征直接写入
                    if is_rna and rna_type:
                        # 标准化RNA类型名称
                        rna_display = rna_type.replace('_', '').upper()
                        f.write(f"{strand_symbol} {gene_start} {gene_end} {gene_name} ({rna_display})\n")

                    # 普通基因特征（必须有exon或cds）
                    elif gene_data['exons'] or gene_data['cds_regions']:
                        f.write(f"{strand_symbol} {gene_start} {gene_end} {gene_name}\n")

                        # 写入UTR（按位置排序）
                        sorted_utrs = sorted(gene_data['utrs'], key=lambda x: x['start'])
                        for utr in sorted_utrs:
                            f.write(f"{utr['start']} {utr['end']} utr\n")

                        # 写入exon（按位置排序）
                        sorted_exons = sorted(gene_data['exons'], key=lambda x: x['start'])
                        for exon in sorted_exons:
                            f.write(f"{exon['start']} {exon['end']} exon\n")

                    # 基因之间添加空行
                    f.write("\n")

        self.file_completed.emit(input_file, output_path)

    def parse_gff3_attributes(self, attributes_str: str) -> Dict[str, List[str]]:
        """
        解析GFF3属性字段
        
        Args:
            attributes_str: 属性字符串（第9列）
            
        Returns:
            属性字典 {key: [value1, value2, ...]}
        """
        attributes = {}
        if not attributes_str or attributes_str == '.':
            return attributes
        
        # 分割属性键值对
        pairs = attributes_str.split(';')
        for pair in pairs:
            if '=' in pair:
                key, value = pair.split('=', 1)
                if key not in attributes:
                    attributes[key] = []
                attributes[key].append(value)
        
        return attributes

    def get_gene_name_from_gff3(self, attributes: Dict[str, List[str]], 
                                 feature_type: str, position: int) -> str:
        """
        从GFF3属性中获取基因名称
        
        Args:
            attributes: 属性字典
            feature_type: 特征类型
            position: 特征位置
            
        Returns:
            基因名称
        """
        # 尝试从不同字段获取名称
        if 'gene' in attributes:
            return attributes['gene'][0]
        elif 'gene_id' in attributes:
            return attributes['gene_id'][0]
        elif 'Name' in attributes:
            return attributes['Name'][0]
        elif 'name' in attributes:
            return attributes['name'][0]
        elif 'locus_tag' in attributes:
            return attributes['locus_tag'][0]
        elif 'ID' in attributes:
            # 对于基因，直接使用ID
            if feature_type == 'gene':
                return attributes['ID'][0]
            # 对于其他特征，尝试从Parent获取基因名
            elif 'Parent' in attributes:
                return attributes['Parent'][0]
            return attributes['ID'][0]
        elif 'Parent' in attributes:
            return attributes['Parent'][0]
        else:
            # 使用位置信息作为标识
            return f"feature_{position}"


class GenBankToMVISTAConverter(QMainWindow):
    """Main window"""

    def __init__(self):
        super().__init__()
        self.input_files = []
        self.file_format = 'genbank'  # Default to GenBank format
        self.init_ui()
        self.conversion_thread = None

    def init_ui(self):
        """Initialize UI"""
        self.setWindowTitle("GenBank/GFF3 to mVISTA Converter")
        self.setGeometry(100, 100, 1000, 700)

        # Central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        # Main layout
        main_layout = QVBoxLayout(central_widget)

        # Title
        title_label = QLabel("GenBank/GFF3 → mVISTA Batch Converter")
        title_label.setFont(QFont("Arial", 16, QFont.Bold))
        title_label.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(title_label)

        # Create splitter
        splitter = QSplitter(Qt.Horizontal)

        # Left panel - file selection and settings
        left_panel = self.create_left_panel()
        splitter.addWidget(left_panel)

        # Right panel - log display
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

        # File format selection group
        format_group = QGroupBox("File Format")
        format_layout = QHBoxLayout()

        self.format_button_group = QButtonGroup()

        self.genbank_radio = QRadioButton("GenBank Format")
        self.genbank_radio.setChecked(True)
        self.genbank_radio.toggled.connect(self.on_format_changed)
        self.format_button_group.addButton(self.genbank_radio)
        format_layout.addWidget(self.genbank_radio)

        self.gff3_radio = QRadioButton("GFF3 Format")
        self.gff3_radio.toggled.connect(self.on_format_changed)
        self.format_button_group.addButton(self.gff3_radio)
        format_layout.addWidget(self.gff3_radio)

        format_group.setLayout(format_layout)
        layout.addWidget(format_group)

        # Input file selection group
        input_group = QGroupBox("Input Files")
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
        dir_layout.addWidget(QLabel("Output Directory:"))
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
        feature_group = QGroupBox("Feature Type Description")
        feature_layout = QVBoxLayout()
        
        # Add description text
        info_label = QLabel("mVISTA format supports the following feature types:\n"
                            "• gene: Gene (with strand marker >/<)\n"
                            "• exon: Exon\n"
                            "• utr: Untranslated region (5'UTR and 3'UTR)\n"
                            "• RNA: rRNA, tRNA, mRNA, snRNA, ncRNA, miRNA, lnc_RNA, etc.")
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
        if self.file_format == 'genbank':
            files, _ = QFileDialog.getOpenFileNames(
                self, "Select GenBank Files", "",
                "GenBank files (*.gb *.gbk *.gbff);;All files (*.*)"
            )
        else:
            files, _ = QFileDialog.getOpenFileNames(
                self, "Select GFF3 Files", "",
                "GFF3 files (*.gff *.gff3);;All files (*.*)"
            )

        if files:
            for file in files:
                if file not in self.input_files:
                    self.input_files.append(file)
                    self.file_list.addItem(Path(file).name)
            self.update_file_count()

    def add_directory(self):
        """Add folder"""
        if self.file_format == 'genbank':
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
        else:
            directory = QFileDialog.getExistingDirectory(
                self, "Select folder containing GFF3 files"
            )
            if directory:
                for file in Path(directory).glob("*.gff*"):
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

    def on_format_changed(self):
        """File format changed"""
        if self.genbank_radio.isChecked():
            self.file_format = 'genbank'
        else:
            self.file_format = 'gff3'
        # Clear file list when format changes
        self.clear_files()

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
        self.log_message("=== Starting Conversion ===")
        self.log_message(f"Input files: {len(self.input_files)}")
        self.log_message(f"Output directory: {output_dir}")
        self.log_message(f"File format: {self.file_format.upper()}")
        self.log_message("Feature types: gene, exon, utr (fixed)")
        self.log_message("")

        # Reset progress bar
        self.progress_bar.setValue(0)

        # Create and start conversion thread
        self.conversion_thread = ConversionThread(
            self.input_files, output_dir, None, self.file_format
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
        """File conversion completion callback"""
        self.log_message(f"Completed: {Path(input_file).name} → {Path(output_file).name}")

    def on_error_occurred(self, input_file: str, error_msg: str):
        """Error occurrence callback"""
        self.log_message(f"Error: {Path(input_file).name}\n  {error_msg}")

    def on_all_completed(self, success_count: int, fail_count: int):
        """All files conversion completion callback"""
        self.log_message("")
        self.log_message("=== Conversion Complete ===")
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
        # Auto-scroll to bottom
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
    """主函数"""
    app = QApplication(sys.argv)
    
    # 设置应用样式
    app.setStyle('Fusion')
    
    # 创建并显示主窗口
    window = GenBankToMVISTAConverter()
    window.show()
    
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
