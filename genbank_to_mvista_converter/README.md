# GenBank to mVISTA Converter

A graphical tool for batch converting GenBank annotation files to mVISTA annotation file format.

## Features

- ✅ Graphical user interface, intuitive and easy to use
- ✅ Batch conversion of multiple GenBank files
- ✅ Import all GenBank files from a folder
- ✅ Automatic detection of gene strand direction (use `>` for plus strand, `<` for minus strand)
- ✅ Automatic grouping of exon and UTR features by gene
- ✅ Real-time conversion progress and detailed logs
- ✅ Statistics display after conversion completion
- ✅ Cross-platform support (Windows, Linux, macOS)

## Installation

Install dependencies using pip:

```bash
pip install -r requirements.txt
```

Or install individually:

```bash
pip install biopython PyQt5
```

## Usage

### Start the program

```bash
python genbank_to_mvista.py
```

### Operation steps

1. **Add input files**
   - Click the "Add Files" button to select one or multiple GenBank files
   - Or click "Add Folder" button to select a folder containing GenBank files
   - The program automatically recognizes files with `.gb`, `.gbk`, `.gbff` extensions

2. **Select output directory**
   - Click "Browse..." button to select output directory
   - Converted mVISTA format files will be saved in this directory

3. **Start conversion**
   - Click "Start Conversion" button
   - View the conversion progress in the log panel on the right
   - Wait for conversion to complete

## mVISTA Format Specification

The converted file uses the standard mVISTA annotation file format, displayed in gene groups.

### Format rules

1. **Gene line**: Starts with `>` (plus strand) or `<` (minus strand), followed by gene start position, end position, and gene name
   ```
   > start end gene_name
   ```
   or
   ```
   < start end gene_name
   ```

2. **UTR line**: Contains UTR start position, end position, and type "utr"
   ```
   start end utr
   ```

3. **exon line**: Contains exon start position, end position, and type "exon"
   ```
   start end exon
   ```

4. **Gene separation**: Each gene block is separated by an empty line

### File example

```
< 106481 116661 gene1
106481 106497 utr
107983 108069 exon
109884 110033 exon
111865 112023 exon

> 39424 42368 gene2
39424 39820 exon
41401 42368 exon

> 77817 81088 gene3
77817 78820 utr
79538 80107 exon
```

### Feature types

mVISTA format only includes the following feature types:

- **gene**: Gene (use `>` for plus strand, `<` for minus strand)
- **exon**: Exon
- **utr**: Untranslated region (including 5'UTR and 3'UTR)

## Output file naming

Output file naming convention:
```
{sequence_id}_mvista.txt
```

For example:
- Input: `chloroplast.gb` (sequence ID: NC_001665)
- Output: `NC_001665_mvista.txt`

## Coordinate system

- Uses **1-based** coordinate system (consistent with GenBank format)
- Start positions count from 1
- Start position ≤ End position

## Notes

1. **File format requirement**: Input files must be in standard GenBank format
2. **Coordinate system**: Uses 1-based coordinate system
3. **Output overwriting**: Existing output files will be overwritten
4. **Gene grouping**: The program automatically groups features by gene name and associates related exon and UTR features with their corresponding genes
5. **Strand direction**: The program automatically extracts strand direction information from the GenBank file, using `>` for plus strand and `<` for minus strand
6. **Gene line**: If no explicit gene type feature exists, the program infers gene start and end positions based on the position range of all related features

## Example Data

### GenBank input example

```
LOCUS       NC_001665             155500 bp    DNA     circular PLN 15-MAR-2024
DEFINITION  Nicotiana tabacum chloroplast, complete genome.
SOURCE      Nicotiana tabacum
  ORGANISM  Nicotiana tabacum
            Eukaryota; Viridiplantae; Streptophyta; Embryophyta;
            Tracheophyta; Spermatophyta; Magnoliopsida;
            eudicotyledons; Core eudicots; Gunneridae; Pentapetalae;
            asterids; lamiids; Solanales; Solanaceae; Nicotiana.
FEATURES             Location/Qualifiers
     gene            complement(106481..116661)
                     /gene="gene1"
                     /locus_tag="gene1"
     mRNA            complement(106481..116661)
                     /gene="gene1"
     exon            complement(107983..108069)
                     /gene="gene1"
     exon            complement(109884..110033)
                     /gene="gene1"
     5'UTR           complement(106481..106497)
                     /gene="gene1"
     gene            39424..42368
                     /gene="gene2"
                     /locus_tag="gene2"
     exon            39424..39820
                     /gene="gene2"
     exon            41401..42368
                     /gene="gene2"
```

### mVISTA output example

```
< 106481 116661 gene1
106481 106497 utr
107983 108069 exon
109884 110033 exon

> 39424 42368 gene2
39424 39820 exon
41401 42368 exon
```

## Troubleshooting

### Q1: Conversion fails with "No valid GenBank records found"

**Cause**: File is not in valid GenBank format

**Solution**:
- Verify file extension is `.gb`, `.gbk`, or `.gbff`
- Check if file content contains complete GenBank format (LOCUS, ORIGIN, etc.)
- Use a text editor to open and verify the file format

### Q2: Some genes do not appear in the output file

**Cause**: GenBank file may not have explicit gene, exon, or UTR features

**Solution**:
- Check if the GenBank file contains these feature types
- Ensure feature naming follows standard conventions (e.g., "gene", "exon", "5'UTR", "3'UTR")

### Q3: Gene strand direction is incorrect

**Cause**: Gene feature in GenBank file does not have explicit strand information

**Solution**:
- Check the strand attribute of gene features in the GenBank file
- If no strand information is available, the program defaults to plus strand (`>`)

### Q4: Program fails to start

**Cause**: Missing dependency libraries

**Solution**:
```bash
pip install -r requirements.txt
```

## Technical Implementation

### Dependencies

| Library | Usage |
|---------|-------|
| biopython | GenBank file parsing |
| PyQt5 | Graphical user interface |

### Core algorithm

1. **File reading**: Read GenBank files using BioPython's SeqIO module
2. **Feature extraction**: Iterate through all features, only extracting gene, exon, and UTR types
3. **Gene grouping**: Group exons and UTRs by gene name under their corresponding genes
4. **Strand direction detection**: Extract strand direction from gene feature's strand attribute
5. **Coordinate conversion**: Convert 0-based coordinates to 1-based coordinates
6. **Data sorting**: Sort output by gene start position and feature positions
7. **File writing**: Write to text file in mVISTA format

## mVISTA Format Specification

According to mVISTA official documentation:

> If a gene annotation of the sequence is available you can submit it in a simple plain text format to be displayed on the plot
>
> Each gene is defined by its start and end coordinates on the sequence, and the name listed on one line. A greater than (>) or less than (<) sign should be placed before this line to indicate plus or minus strand, although the numbering should be according to the plus strand. The exons are listed individually with the word "exon", after the start and end coordinates of each exon. UTRs are annotated the same way exons, with the word "utr" replacing "exon".

## Version History

### v2.0 (2026-02-01)
- Updated to correct mVISTA format specification
- Gene lines use `>` (plus strand) or `<` (minus strand) to indicate direction
- Only outputs three feature types: gene, exon, and UTR
- Removed feature type filter options (format is now fixed)
- Fixed strand attribute missing error

### v1.0 (2026-02-01)
- Initial version
- Support for basic GenBank to mVISTA format conversion
- Graphical user interface
- Batch processing functionality

## License

This tool is for learning and research purposes only.

## Contact

If you have questions or suggestions, please feel free to provide feedback.
