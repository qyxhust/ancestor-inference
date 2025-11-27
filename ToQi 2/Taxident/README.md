# TaxIdent - Taxonomic Hierarchy BAM Processor

TaxIdent is a C program that processes BAM files to produce carious forms of mismatch counts. It supports phylogenetic inferences of mismatches to inferred ancestors using NCBI taxonomy data. 

## Table of Contents

- [Key Features](#key-features)
- [Installation](#installation)
- [Required Data Files](#required-data-files)
- [Usage](#usage)
- [Input File Formats](#input-file-formats)
- [Output Formats](#output-formats)
- [Quick Reference](#quick-reference)
- [Troubleshooting](#troubleshooting)
- [Technical Details](#technical-details)
- [Algorithm Behavior](#algorithm-behavior)

## Multiple Distinct Processing Modes

TaxIdent operates in several different modes depending on whether taxonomic analysis and depending on how the data should be processed.

### **Genome-Based Mode (Default)**
- **Activation**: Default if not doing taxonomic analysis
- **Processing**: Groups alignments according to reference genome. If a genome is divided into multiple references, so will the catagorization
- **Output columns**: Represent individual genome files (G1, G2, G3...)

```bash
# Example: Genome-based analysis
./TaxIdent -i sample.bam -o output.txt
```

### **RG-Based Mode **
- **Activation**: Run with `--use-rg` flag
- **Processing**: Groups alignments according to the RG tag in the BAm file. 
- **Output columns**: Represent individual RG tags (G1, G2, G3...)

```bash
# Example: Genome-based analysis
./TaxIdent -i sample.bam -o output.txt --use-rg
```

### **Taxid-Based Mode **
- **Activation**: Run with --consolidate-by-taxid
- **Processing**: Consolidates genomes by taxonomic ID, keeps best alignment per taxid
- **Output columns**: Represent unique taxonomic IDs (T12345, T67890...)

```bash
# Example: Genome-based analysis
./TaxIdent -i sample.bam -o output.txt --consolidate-by-taxid
```

### **Taxid-Based Mode with Taxonomic Analysis**
- **Activation**: Use `--with-higher-taxa` flag with taxonomy files
- **Processing**: Automatically consolidates genomes by taxonomic ID, keeps best alignment per taxid
- **Output columns**: Represent unique taxonomic IDs (T12345, T67890...). Taxids after '|' represent hgiher taxonomic grous.
- **Use case**: Taxonomic classification

```bash
# Example: Taxid-based taxonomic analysis  
./TaxIdent -i sample.bam -o output.txt \
    --taxonomy-dir ncbi_20250530 \
    --acc2taxid namelookup/wgs_eukaryota.acc2taxid \
    --with-higher-taxa
```

**Important**: These modes produce fundamentally different outputs - choose based on your analysis goals.

### There are two different output format: Dense and sparse. Dense represents the mismatch data in a full matrix representing all genomes/RGs by all reads. Sparse uses a sparse representation where only the genomes/RGs/taxids with alignments reported are given for each read. The output also differs depending on whether damages are represented or not.

### Output Format Examples

#### **Genome-Based Mode output with damages **
```
read_id	total_count	G1	G2	G3	G4	G5...
read001	100	G1	7	1	4	G3	8	0	5	G5	6	2	3
```
-- **Key file**: Maps G1→genome_file.fasta, G2→another_genome.fasta

### **Taxid-Based Mode Output with damages** 
```
read_id	total_count	T12345	T67890	T54321	|	T100:genus	T200:family
read001	100	T12345	7	1	4	T67890	8	0	5	|	T100:genus	8	0	6	T200:family	8	0	6
```
- **Key file** (optional): Maps T1→T12345, T2→T67890 when `--key-file` used


## Features

- **Read Trimming**: Configurable trimming of bases from read ends (default: no trimming for non-damage mode, 5bp for damage mode)
- **Multiple Input Modes**: Supports both single BAM files and directories containing multiple BAM files
- **Multi-file Processing**: `--no-redundancy` mode available when reads don't overlap between files and can speed up calculations
- **Genome Auto-discovery**: Automatically discovers genome names from BAM files when no genome list is provided
- **Genome Name Flexibility**: Handles inconsistent genome naming with configurable ignore patterns
- **Quality Filtering**: Filters reads by mapping quality and minimum length requirements
- **Precise Mismatch Detection**: Uses MD tags and CIGAR strings for accurate mismatch counting
- **Best Alignment Selection**: When processing multiple alignments to same genome, automatically selects the alignment with the fewest mismatches
- **Multi-threaded I/O**: Supports parallel decompression for faster processing of large BAM files. This is only for htslib processing and does not have a major impact on computational speed.
- **Smart Defaults**: Defaults to sparse_damage format with damage analysis, precise indel handling, and short genome names enabled
- **Automatic Genome Key File**: Automatically creates a genome_key.txt file mapping short names (G1, G2, etc.) to full genome names when using short names
- ** Parsimony assignment of ancestral states** Uses a parsimony algorithm to reconstrcuct ancestral states to generat mismatch counts for higther taxonomic groups
- **NCBI Taxonomy Integration**: Uses official NCBI taxonomy files for accurate hierarchical structure
- **Flexible Taxonomic Levels**: Select specific taxonomic ranks to include in analysis
- **Taxonomy Tree Visualization**: Can output taxonomy tree structure for debugging and analysis



## Program Defaults 

The program uses these defaults:

1. **Output format**: `sparse_damage` 
2. **Damage analysis**: ENABLED (use `--no-damage` to disable)  
3. **Precise indel handling**: ENABLED for accurate complex indel processing. However, when analyzing higher taxonomic groups, reads with insertions are skipped.
4. **Short genome names**: ENABLED (use `--full-names` for full genome names in genome mode)
5. **Key files**: Automatically created when using short genome names (use `--no-key-file` to disable)
6. **Genome identification**: Reference names (BAM file column 3) not RG tags (use `--use-rg` for RG tags)
7. **Internal data structures**: Sparse structures for performance with many genomes (use `--dense` when all reads have alignments to all groups analyzed)
8. **Taxonomy**: Disabled (use `--taxonomy-dir` + `--acc2taxid` to enable, add `--with-higher-taxa` for taxonomic groups)
9. **Threading**: Single-threaded (use `-t N` for multi-threading)

For basic BAM processing (without taxonomy), most users can simply run:
```bash
./TaxIdent -i input.bam -o output.txt
```

For taxonomic analysis (leaf genomes only):
```bash
./TaxIdent -i input.bam -o output.txt --taxonomy-dir ncbi_20250530 --acc2taxid namelookup/wgs_eukaryota.acc2taxid
```

For taxonomic analysis with higher taxonomic groups (family, genus, order, etc.):
```bash
./TaxIdent -i input.bam -o output.txt --taxonomy-dir ncbi_20250530 --acc2taxid namelookup/wgs_eukaryota.acc2taxid --with-higher-taxa
```

## Installation

### Requirements

- HTSlib (for BAM file processing)
- GCC compiler
- Linux/Unix/macOS environment

### Compilation on Dandy System

The Dandy system requires loading HTSlib as a module before compilation and execution.

#### Loading Required Modules

```bash
# Load modules in this exact order (libdeflate is a prerequisite for htslib)
module load libdeflate/1.21 htslib/1.21
```

**IMPORTANT**: The modules must be loaded in every new session before compilation or execution.

#### Compilation

```bash
# After loading modules, compile with RPATH to embed library paths:
gcc -DWITH_HTSLIB -O3 -I/opt/software/htslib/1.21/include -L/opt/software/htslib/1.21/lib -Wl,-rpath,/opt/software/htslib/1.21/lib -o TaxIdent TaxIdent.c -lhts -lz -lm -lpthread
```

### Compilation on Local Systems (macOS/Linux)

For local development (e.g., on macOS), use the following compilation command:

```bash
# Compile with local HTSlib installation
gcc -o TaxIdent TaxIdent.c -O3 -D_GNU_SOURCE -DWITH_HTSLIB \
    -I/path/to/htslib/include \
    -L/path/to/htslib/lib \
    -lhts -lz -lm -lbz2 -llzma -lcurl -lpthread

# Example for Rasmus's local system:
gcc -o TaxIdent TaxIdent.c -O3 -D_GNU_SOURCE -DWITH_HTSLIB \
    -I/Users/rasmus_nielsen/local/include \
    -L/Users/rasmus_nielsen/local/lib \
    -lhts -lz -lm -lbz2 -llzma -lcurl -lpthread

# On macOS, fix library paths after compilation:
install_name_tool -change /usr/local/lib/libhts.3.dylib \
    /Users/rasmus_nielsen/local/lib/libhts.3.dylib \
    TaxIdent
```

## Required Data Files

### NCBI Taxonomy Files (if analyzing higher taxonomic groups)

Download from NCBI FTP site or use provided files in `ncbi_20250530/`:
- `nodes.dmp`: Taxonomic hierarchy structure (parent-child relationships)
- `names.dmp`: Taxonomic names and synonyms

### Accession to TaxID Mapping

Use files from `namelookup/` directory:
- `core_nt.acc2taxid`: Core nucleotide accessions
- `wgs_eukaryota.acc2taxid`: Whole genome shotgun eukaryotic accessions
- Other specialized mapping files as needed

**Note**: These files should be uncompressed for use. If provided as `.gz` files, uncompress them first:
```bash
gunzip namelookup/*.gz
```

## Usage

### Basic Syntax

```bash
./TaxIdent [OPTIONS]
```

### Required Arguments

- `-i, --input FILE/DIR`: Input BAM file or directory containing BAM files
- `-o, --output FILE`: Output file for the mismatch matrix

### Optional Arguments

#### Basic Options

- `-g, --genomes FILE`: Text file listing genome names [default: auto-discover from BAM]
- `-q, --min-mapq INT`: Minimum mapping quality [default: 0]
- `-n, --max-reads INT`: Maximum number of alignment records to process [default: unlimited]
- `-I, --ignore CHAR`: Character to ignore in genome names [default: none]
- `-t, --threads N`: Number of decompression threads for faster I/O [default: 0]
- `-v, --verbose`: Enable verbose output [default: disabled]
- `--silent`: Suppress all progress output (overrides --verbose) [default: disabled]
- `-h, --help`: Show help message

#### Damage Analysis Options

- `-s, --damage-sites N`: Number of damage-susceptible sites at read ends [default: 5 for damage mode, 0 for non-damage mode]
- `--damage`: Enable damage analysis (C->T, G->A in first/last s sites) [**ENABLED by default**]
- `--no-damage`: Disable damage analysis [default: disabled]
- `--asymmetric-damage`: Use asymmetric damage calculation (C->T in first s, G->A in last s) [default: disabled - uses symmetric]

#### Output Format Options

- `--format FORMAT`: Output format: dense, sparse, dense_damage, sparse_damage [**default: sparse_damage**]
- `--short-names`: Use short genome names (G1, G2, etc.) [**ENABLED by default**]
- `--full-names`: Use full genome names instead of short names [default: disabled]
- `--compress`: Use compressed output with short IDs [default: disabled]
- `--taxid-file FILE`: NCBI names.dmp file for taxid mapping (with --compress) [default: none]
- `--genome-map FILE`: Output file for genome ID to name mapping (with --compress) [default: none]

#### Processing Options

- `--text`: Process text mismatch matrix instead of BAM file [default: disabled - uses BAM mode]
- `--simple`: Use simple unoptimized algorithms (for testing/debugging) [default: disabled - uses optimized mode]
- `--use-rg`: Use RG (Read Group) tags for genome identification [default: disabled - uses reference names]
- `--dense`: Use dense internal data structures [default: disabled - uses sparse structures]
- `--skip-indels`: Skip alignments containing indels [default: disabled - processes indels]
- `--precise-indels`: Use precise calculation for alignments with indels [**ENABLED by default**]
- `--no-redundancy`: Assume reads don't appear in multiple BAM files [default: disabled - allows redundancy]
- `--parallel-files N`: Number of BAM files to process in parallel (directory mode) [default: 1]
- `--enforce-dense_strict`: Only output reads with alignments to ALL genomes/taxids/readgroups [default: disabled]
- `--penalty-mode`: For dense output formats, use max+1 penalty for missing alignments instead of -1 [default: disabled - uses -1]
  - When enabled: Genomes without alignments get a penalty value of (max_mismatches + 1), capped at read_length
  - When disabled (default): Genomes without alignments get -1 to indicate missing data
  - Only affects dense and dense_damage output formats

#### Taxonomy Options

- `--taxonomy-dir DIR`: Directory containing NCBI taxonomy files (nodes.dmp, names.dmp) [default: none - required for taxonomy]
- `--acc2taxid FILE`: Accession to taxid mapping file [default: none - required for taxonomy]
- `--with-higher-taxa`: Include higher taxonomic groups (family, genus, order, etc.) in output after "|" separator [default: disabled - outputs only leaf genomes]
- `--higher_taxa_with_mbs`: Use original ALL/SOME method for taxonomic groups (outputs 5 values) instead of parsimony method (3 values) [default: disabled - uses parsimony]
- `--consolidate-by-taxid`: Consolidate genomes by taxid without taxonomy tree (**sparse format only**, requires --acc2taxid) [default: disabled]
- `--consecutive`: Optimize for consecutive alignments (reduces memory usage) [default: disabled]
- `--test-tree`: Test taxonomy tree parsing and output tree structure (exits after test) [default: disabled]

**Note**: The `--tax-levels` option is parsed but not currently implemented. All taxonomic levels present in the data are included when using `--with-higher-taxa`.

## Examples

### Basic Examples

```bash
# On Dandy system - load modules first
module load libdeflate/1.21 htslib/1.21

# SIMPLIFIED: Process a single BAM file with taxonomy (uses sparse_damage format with short names by default)
./TaxIdent -i sample.bam -o output.txt \
    --taxonomy-dir ncbi_20250530 \
    --acc2taxid namelookup/wgs_eukaryota.acc2taxid \
    --with-higher-taxa

# With genome list file for specific ordering
./TaxIdent -i sample.bam -g genome_names.txt -o output.txt \
    --taxonomy-dir ncbi_20250530 \
    --acc2taxid namelookup/core_nt.acc2taxid \
    --with-higher-taxa

# With genome name processing (ignore everything before first underscore)
./TaxIdent -i sample.bam -o output.txt -I _ \
    --taxonomy-dir ncbi_20250530 \
    --acc2taxid namelookup/wgs_eukaryota.acc2taxid \
    --with-higher-taxa

# Optimize for consecutive alignments (lower memory)
./TaxIdent -i sample.bam -o output.txt \
    --taxonomy-dir ncbi_20250530 \
    --acc2taxid namelookup/wgs_eukaryota.acc2taxid \
    --with-higher-taxa \
    --consecutive

# Test taxonomy tree structure (debugging)
./TaxIdent --test-tree --taxonomy-dir ncbi_20250530 \
    -i sample.bam -o test.txt -n 5
# This will output tree_test.txt with parent-child relationships
```

### Directory Processing (Multiple BAM Files)

```bash
# Process all BAM files in a directory with taxonomy
./TaxIdent -i /path/to/bam/files/ -o combined.txt -I _ \
    --taxonomy-dir ncbi_20250530 \
    --acc2taxid namelookup/wgs_eukaryota.acc2taxid \
    --with-higher-taxa

# OPTIMIZED: Use --no-redundancy for 3.7x speedup when reads don't overlap between files
./TaxIdent -i /path/to/bam/files/ -o combined.txt -I _ \
    --taxonomy-dir ncbi_20250530 \
    --acc2taxid namelookup/wgs_eukaryota.acc2taxid \
    --with-higher-taxa \
    --no-redundancy

# With multi-threading for faster processing of large BAM files
./TaxIdent -i /path/to/bam/files/ -o combined.txt -I _ -t 8 -v \
    --taxonomy-dir ncbi_20250530 \
    --acc2taxid namelookup/wgs_eukaryota.acc2taxid \
    --with-higher-taxa

# Combining optimizations for maximum performance
./TaxIdent -i /path/to/bam/files/ -o combined.txt -I _ \
    --no-redundancy -t 8 -v \
    --taxonomy-dir ncbi_20250530 \
    --acc2taxid namelookup/wgs_eukaryota.acc2taxid \
    --with-higher-taxa
```

### Non-Taxonomy Mode (Default Behavior)

```bash
# Process without taxonomy (default behavior - like original BAMreader)
./TaxIdent -i sample.bam -o output.txt

# Disable damage analysis (if not needed)
./TaxIdent -i sample.bam -o output.txt --no-damage

# Use dense format with full genome names (old defaults)
./TaxIdent -i sample.bam -o output.txt --format dense --no-damage --full-names

# With quality filtering and read limit
./TaxIdent -i sample.bam -g genome_names.txt -o output.txt -q 30 -n 100000 -v
```

## Input File Formats

### Genome Names File

A text file containing one genome name per line. Example:

```
10TJ18_pseudomolecules_and_unplaced_contigs_CPclean.srt
Aizu6_pseudomolecules_and_unplaced_contigs_CPclean.srt
Akashinriki_pseudomolecules_and_unplaced_contigs_CPclean.srt
BCC_1451_ragtag.srt
```

### BAM File Requirements

- Reads must be aligned to multiple reference genomes
- MD tags preferred for precise mismatch detection (NM tags used as fallback)
- For RG-based identification: RG (Read Group) tags must be present
- For reference-based identification (default): Reference names should correspond to genome names

### Taxonomy Files Format

#### nodes.dmp
Tab-separated file with taxonomic hierarchy:
```
taxid | parent_taxid | rank | ...
```

#### names.dmp
Tab-separated file with taxonomic names:
```
taxid | name | unique_name | name_class |
```

#### acc2taxid
Tab-separated file mapping accessions to taxids:
```
accession accession.version taxid gi
```

## Genome Identification Methods

TaxIdent supports two methods for identifying which genome each read alignment belongs to:

### Method 1: Reference Names (Default)

By default, TaxIdent uses the reference sequence name (column 3 in SAM format) to identify the target genome:

- **How it works**: Extracts the reference name from each alignment record
- **Processing**: Applies ignore character (`-I`) processing if specified
- **Example**: Reference name `chr1_Genome_A_v1.0` with `-I _` becomes `Genome_A_v1.0`
- **Advantages**: 
  - Works with any BAM file aligned to multiple references
  - More genomes typically discovered (reference sequences > RG groups)
  - No dependency on RG tag consistency

```bash
# Default reference-based identification
./TaxIdent -i multi_ref.bam -o output.txt -I _
```

### Method 2: RG Tags (Optional)

Use the `--use-rg` flag to identify genomes via RG (Read Group) tags:

- **How it works**: Extracts genome name from the RG tag of each alignment
- **Processing**: Applies ignore character (`-I`) processing if specified  
- **Example**: RG tag `sample_Genome_A_v1.0.bam` with `-I _` becomes `Genome_A_v1.0.bam`
- **Advantages**:
  - Explicitly groups reads by intended target genome
  - Fewer genomes typically discovered (cleaner grouping)
  - Useful when reference names don't match genome identities

```bash
# RG-based identification
./TaxIdent -i multi_ref.bam -o output.txt --use-rg -I _
```

## Output Formats

The program supports eight different output formats, combining three dimensions:
1. Dense vs Sparse
2. With damage vs Without damage
3. With higher taxa vs Leaf nodes only

### Format Structure with Taxonomy

When using `--with-higher-taxa`, the output includes both individual genomes and taxonomic groups:

```
read_id total_count G1 G2 G3 | T100:genus T200:family T300:order
```

The "|" separator divides:
- **Left side**: Individual reference genomes (leaf nodes)
- **Right side**: Higher taxonomic groups

### Output Fields

#### Individual Genomes (Leaf Nodes)
- **With damage**: `nd md mb` (3 numbers)
  - nd: number of damage-susceptible sites (C/G in damage regions)
  - md: damage mismatches (C->T, G->A)
  - mb: background mismatches (all other)
- **Without damage**: Single number (mb + md)

#### Taxonomic Groups

**Default Parsimony Method** 
- **With damage**: `nd md mb` (3 numbers, same format as individual genomes)
  - nd: maximum damage sites across all member genomes
  - md: damage mismatches (always 0 for parsimony method)
  - mb: background mismatches calculated via phylogenetic parsimony
- **Without damage**: Single number (mb from parsimony)

**Legacy ALL/SOME Method** (with `--higher_taxa_with_mbs`):
- **With damage**: `nd md mds mb mbs` (5 numbers)
  - nd: maximum damage sites across all member genomes
  - md: damage mismatches to ALL member genomes
  - mds: damage mismatches to SOME member genomes
  - mb: background mismatches to ALL member genomes
  - mbs: background mismatches to SOME member genomes
- **Without damage**: Two numbers (mb+md, mbs+mds)



## Technical Details

### Processing Algorithm

TaxIdent uses a **linear-complexity algorithm** for taxonomic mismatch calculation:

**Phase 1: Data Collection**
1. Read BAM file(s) to discover reference genomes
2. Map reference names to taxids using acc2taxid file
3. Build global taxonomic tree from NCBI files
4. Store position-specific mismatch data for each read

**Phase 2: Per-Read Processing**
1. For each read, build minimal taxonomy tree containing only aligned genomes
2. Calculate taxonomic group statistics using one of two methods:
   - **Parsimony Method (default)**: Uses phylogenetic parsimony reconstruction with post-order and pre-order tree traversals
   - **ALL/SOME Method (legacy)**: Uses bit vector operations to track mismatches to ALL vs SOME member genomes
3. Output results for leaf genomes and taxonomic groups

### Parsimony Algorithm Details (New Default)

The new parsimony method uses phylogenetic parsimony reconstruction to calculate a single mismatch value (mb) for each taxonomic group:

**Algorithm Steps:**
1. **Leaf Initialization**: Each genome sets state = 1 for positions with non-damage mismatches, 0 otherwise
2. **Post-order Traversal**: Internal nodes apply parsimony rules:
   - State = 0 if at least one child = 0 and no child = 1
   - State = 1 if at least one child = 1 and no child = 0  
   - State = 0|1 (ambiguous) otherwise
3. **Pre-order Resolution**: Resolve ambiguous states:
   - Root: 0|1 → 1
   - Internal nodes: 0|1 → inherit parent's resolved state
4. **Final Count**: Count positions with state = 1 to get mb value


### Memory Management Options

**Consecutive Mode** (`--consecutive`):
- Processes all alignments for a read immediately
- Lower memory footprint
- Requires alignments grouped by read

**Default Mode**:
- Stores all mismatch positions until end
- Handles any alignment order
- More flexible for various input formats

### Best Alignment Selection

When multiple alignments to the same genome/RG/taxid exist for a read:
1. Compare NM (total mismatches) values from BAM tags
2. Keep only the alignment with minimum NM value
3. Clear any previously stored mismatch data for worse alignments
4. This ensures each genome gets at most one alignment per read

### Read Processing Details

#### Trimming Behavior
- **Damage Mode (default)**: No trimming by default (use `-s N` to set damage region size)
- **Non-damage Mode**: No trimming by default (use `-s N` to trim N bases from each end)
- **Minimum Length**: Reads shorter than 20 bp total (after any trimming) are filtered out
- **Length Reporting**: Output reports the actual analyzed read length

#### Genome Name Processing
When using the `-I` flag, genome names are processed as follows:
- **Processing**: Everything before and including the specified character is ignored
- **RG Example**: With `-I _`, RG tag `445.CGG3016660.20250324Asplund_BCC_1451_ragtag.srt` becomes `BCC_1451_ragtag.srt`
- **Reference Example**: With `-I _`, reference name `chr1_BCC_1451_ragtag.srt` becomes `BCC_1451_ragtag.srt`


## Error Handling

The program will report errors for:
- Missing or unreadable input files
- BAM files without required tags (MD/NM for mismatches)
- Genome names not found in the reference file
- Memory allocation failures
- Corrupted BAM files
- Missing taxonomy files
- Invalid acc2taxid mappings

Warnings are issued for:
- Reads missing RG tags when in RG mode(skipped)
- Unknown genome names in RG tags when in RG mode (skipped)
- Reads too short after trimming (skipped)
- Genomes without taxid assignments (assigned -1)
- Multiple alignments to same genome (best kept)


### Slow Processing
- Ensure acc2taxid file is uncompressed
- Use `--no-redundancy` for directory processing
- Enable multi-threading with `-t N`
- Use sparse formats for files with many genomes

### Missing Taxonomic Assignments
- Verify accession format matches acc2taxid file
- Check that taxonomy files are complete and uncompressed
- Some genomes may not have taxid assignments (normal)
- Try different acc2taxid files (core_nt.acc2taxid vs wgs_eukaryota.acc2taxid)

### Compilation Issues
- Ensure HTSlib is properly installed
- On Dandy: Load required modules first
- On macOS: Fix library paths with install_name_tool
- Check CLAUDE_MEMORY.md for detailed compilation instructions


### Most Common Use Cases

```bash
# Basic BAM processing (default behavior - no taxonomy)
./TaxIdent -i sample.bam -o output.txt

# Taxonomic processing (leaf genomes only)
./TaxIdent -i sample.bam -o output.txt \
    --taxonomy-dir ncbi_20250530 \
    --acc2taxid namelookup/wgs_eukaryota.acc2taxid

# Taxonomic processing with higher taxonomic groups (family, genus, order, etc.)
# Uses new parsimony method (3 values per taxonomic group)
./TaxIdent -i sample.bam -o output.txt \
    --taxonomy-dir ncbi_20250530 \
    --acc2taxid namelookup/wgs_eukaryota.acc2taxid \
    --with-higher-taxa

# Use legacy ALL/SOME method (5 values per taxonomic group)
./TaxIdent -i sample.bam -o output.txt \
    --taxonomy-dir ncbi_20250530 \
    --acc2taxid namelookup/wgs_eukaryota.acc2taxid \
    --with-higher-taxa \
    --higher_taxa_with_mbs

# With genome name processing and threading (recommended for large files)
./TaxIdent -i sample.bam -o output.txt -I _ -t 8 \
    --taxonomy-dir ncbi_20250530 \
    --acc2taxid namelookup/wgs_eukaryota.acc2taxid \
    --with-higher-taxa

# Large-scale directory processing with progress monitoring
./TaxIdent -i /bam/directory/ -o matrix.txt -I _ -t 8 -v \
    --taxonomy-dir ncbi_20250530 \
    --acc2taxid namelookup/wgs_eukaryota.acc2taxid \
    --with-higher-taxa

# OPTIMIZED: Maximum performance for directory processing (3.7x speedup)
./TaxIdent -i /bam/directory/ -o matrix.txt -I _ --no-redundancy -t 8 \
    --taxonomy-dir ncbi_20250530 \
    --acc2taxid namelookup/wgs_eukaryota.acc2taxid \
    --with-higher-taxa

# For non-damage analysis (e.g., modern DNA)
./TaxIdent -i modern.bam -o output.txt --no-damage --format sparse -I _ -t 8

# To use full genome names (old behavior)
./TaxIdent -i sample.bam -o output.txt --full-names -I _ -t 8

# Test taxonomy tree structure
./TaxIdent --test-tree --taxonomy-dir ncbi_20250530 -i sample.bam -o test.txt -n 5
```
## Citation

If you use TaxIdent in your research, please cite:
[Citation information to be added]

## Known Limitations and Incompatibilities

### Sparse Internal Mode with Dense Output Format
Sparse internal mode (default, without `--dense` flag) is **not compatible** with dense output formats (`--format dense` or `--format dense_damage`).

This incompatibility exists because the sparse internal mode uses an optimized processing path that writes data to a temporary file in sparse format and cannot be converted to dense output.

**Error message**: The program will exit with an error if you attempt to use sparse internal mode with dense output format.

**Solutions**:
1. Use dense internal mode: `./TaxIdent -i input.bam -o output.txt --format dense --dense`
2. Use sparse output format: `./TaxIdent -i input.bam -o output.txt --format sparse`

## License

[License information to be added]

## Contact

For questions or issues, please contact:
[Contact information to be added]