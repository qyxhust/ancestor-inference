# MLID - Machine Learning Identifier for Taxonomic Classification

A high-performance C implementation of EM algorithms for estimating mixture proportions from mismatch matrix data, with specialized support for per-genome error rate estimation and damage error modeling in ancient DNA.

## Installation

### Prerequisites
- GCC compiler with C99 support
- Make
- Math library (libm)

### Build
```bash
make
```

For debug build:
```bash
make debug
```

For clean build:
```bash
make clean && make
```

## Usage

### Basic Usage
```bash
# Standard Method BEfull with per-genome error rates
./mlid -i data.txt -M BEfull

# Damage Method CEfull with per-genome error rates (fixed damage rate)
./mlid -i damage_data.txt -M CEfull -r 0.01

# Damage Method CEDfull with per-genome error and damage rates
./mlid -i damage_data.txt -M CEDfull

# With filtering and verbose output
./mlid -i data.txt -M BEfull -f 0.001 -v

# Custom output file naming
./mlid -i data.txt -M A -o myresults              # Creates: myresults_out.txt
./mlid -i data.txt -M A -o /path/to/analysis      # Creates: /path/to/analysis_out.txt
./mlid -i data.txt -M A -o results/experiment1    # Creates: results/experiment1_out.txt

# With parameter constraints
./mlid -i data.txt -M CEDfull -C constraints.txt  # Apply custom bounds on error rates

# Force data format processing
./mlid -i large_data.txt -M BEfull -p             # Force sparse processing
./mlid -i small_data.txt -M BEfull -d             # Force dense processing

# Genome name translation
./mlid -i short_names.txt -M CEDfull              # Auto-detect genome key file
./mlid -i data.txt -M BEfull -k genome.key        # Use specific genome key file
```

### Command Line Options

```
Required Options:
  -i, --input FILE      Input mismatch matrix file (TSV format)
  -M, --method METHOD   Method to use: A, AE, B, BE, C, CE, CED, BEfull, CEfull, CEDfull

Standard Method Options:
  -e, --error-rate RATE Background error rate (default: 0.005)

Damage Method Options:
  -e, --error-rate RATE Background error rate (default: 0.005)
  -r, --damage-rate RATE Damage error rate (default: 0.01)
  -j, --joint-damage    Estimate single shared damage rate for all entities (CEDfull only)

General Options:
  -o, --output PREFIX   Output file prefix (default: input filename)
                        Creates PREFIX_out.txt. Can include full path.
  -t, --tolerance TOL   Convergence tolerance (default: 1e-6)
  -m, --max-iter N      Maximum iterations (default: 1000)
  -f, --filter THRESH   Dynamic pruning threshold (see below)
  -s, --squarem         Force SQUAREM acceleration (default: enabled)
  -S, --no-squarem      Disable SQUAREM, use standard EM
  -v, --verbose         Enable verbose output
  -h, --help            Show this help message
  -V, --version         Show version information

Data Format Options:
  -p, --sparse          Force sparse data processing
  -d, --dense           Force dense data processing

Name Translation Options:
  -c, --convert-names   Enable genome name conversion (requires --mapping-file)
  -g, --mapping-file FILE Genome mapping file for name conversion
  -k, --genome-key FILE Specify genome key file for short name translation
  -K, --use-genome-key  Force use of genome key file (auto-detect by default)
  -N, --no-genome-key   Disable genome key file detection and translation

Advanced Options:
  -C, --constraints FILE Parameter constraint file for bounds on error rates (see discussion below)
  -W, --weights FILE Parameter file specifying the prior of the Dirichlet distribution (see discussion below)
  -D, --debug-likelihood FILE Debug mode: load parameters from FILE and calculate likelihood
  -T, --test-likelihood-ratio FILE Likelihood ratio test: test null hypothesis using FILE

Dynamic Pruning:
  -f, --filter THRESH   Enable dynamic pruning during EM iterations
                       Default: min(0.00001, 10/n_reads)
                       Set to 0 to disable pruning
                       
                       Pruning iterations: 1,2,3,4,5,10,20,30,40,50,60,70,80,90,100,150,200,250,300,350,400,450,500,550...
                       Removes genomes with proportion < THRESH
                       Removes reads with no remaining alignments
                       Removed reads saved to: [input_file].removed_reads.txt

Parameter Constraints:
  -C, --constraints FILE Parameter constraint file for bounds on error rates
                       Allows setting custom bounds on error/damage rates for specific taxonomic groups
                       File format: one constraint per line
                       Examples:
                         species < 0.001          # Upper bound constraint
                         genus > 0.01 < 0.1       # Range constraint (lower and upper bounds)
                       Applied during M-step to enforce biological/technical constraints
```

### Methods

| Method | Description | Parameters | Data Format | R Equivalent |
|--------|-------------|------------|-------------|--------------|
| **A**  | Single error rate | Fixed | Standard | 0A |
| **AE** | Single error rate with estimation | Estimated | Standard | 0AE |
| **B**  | Single error rate (÷3 for substitutions) | Fixed | Standard | 0B |
| **BE** | Single error rate (÷3) with estimation | Estimated | Standard | 0BE |
| **C**  | Damage + background errors (fixed rates) | Fixed | Damage | - |
| **CE** | Damage (fixed) + background (estimated) | Mixed | Damage | - |
| **CED** | Damage + background (both estimated) | Estimated | Damage | - |
| **BEfull** | **Per-genome error rates (standard model)** | **2K-1** | Standard | - |
| **CEfull** | **Per-genome error rates (damage model)** | **2K-1** | Damage | - |
| **CEDfull** | **Per-genome error and damage rates** | **3K-1** | Damage | - |

### Full Models (New Features)

The **full models** (BEfull, CEfull, CEDfull) estimate individual error rates for each reference genome, providing much more detailed insights into sequencing patterns:

#### BEfull: Per-Genome Error Rates (Standard Model)
- **Parameters**: K proportions + K error rates = 2K-1 parameters (K-1 proportions + K error rates)
- **Use case**: When different reference genomes have systematically different error patterns
- **Model**: Each genome j has its own error rate εⱼ
- **Likelihood**: (εⱼ/3)^dᵢⱼ × (1-εⱼ)^(nᵢⱼ-dᵢⱼ)

#### CEfull: Per-Genome Error Rates (Damage Model)
- **Parameters**: K proportions + K error rates + 1 damage rate = 2K parameters total (K-1 + K + 1)
- **Use case**: Ancient DNA with per-genome background error variation, fixed damage rate
- **Model**: Each genome j has error rate εⱼ, shared damage rate eᵈ
- **Likelihood**: (1-eᵈ)^(ndᵢⱼ-mdᵢⱼ) × (eᵈ)^mdᵢⱼ × (1-εⱼ)^(nᵢⱼ-mdᵢⱼ-mbᵢⱼ) × (εⱼ/3)^mbᵢⱼ

#### CEDfull: Per-Genome Error and Damage Rates
- **Parameters**: K proportions + K error rates + K damage rates = 3K-1 parameters (default)
  - With `-j` flag: K proportions + K error rates + 1 shared damage rate = 2K parameters (~32% reduction)
- **Use case**: Ancient DNA with both per-genome error and damage variation
- **Model**: Each genome j has error rate εⱼ and damage rate eⱼᵈ (or shared eᵈ with `-j`)
- **Taxonomic support**: When taxonomic data present, estimates per-taxonomic group rates
- **Joint damage mode** (`-j`): Single shared damage rate for all entities, faster convergence
- **Likelihood**: (1-eⱼᵈ)^(ndᵢⱼ-mdᵢⱼ) × (eⱼᵈ)^mdᵢⱼ × (1-εⱼ)^(nᵢⱼ-mdᵢⱼ-mbᵢⱼ) × (εⱼ/3)^mbᵢⱼ

### Input Formats

The program supports **four input formats** and automatically detects which format is used.

#### Format 1: Dense Standard Format
Traditional format with mismatch counts for each genome:
```
read_id	total_count	Genome_A	Genome_B	Genome_C	Genome_D	Genome_E
read_001	50	2	1	0	3	0
read_002	75	0	4	2	1	0
read_003	100	5	2	0	0	1
```

- **Use with**: Methods A, AE, B, BE, BEfull
- **Missing alignments**: Use `0` for no mismatches, `-1` for no alignment
- **Use case**: Smaller datasets with most read-genome alignments

#### Format 2: Dense Damage Format  
Damage format with nd (damage sites), md (damage mismatches), mb (background mismatches):
```
read_id	total_count	nd_Genome_A	md_Genome_A	mb_Genome_A	nd_Genome_B	md_Genome_B	mb_Genome_B
read_001	50	15	2	1	15	0	3
read_002	75	20	3	2	20	1	1
read_003	100	25	5	0	25	0	2
```

- **Use with**: Methods C, CE, CED, CEfull, CEDfull
- **nd_X**: Number of damage-susceptible sites (e.g., C/G in first/last 5 positions)
- **md_X**: Number of damage mismatches (e.g., C→T transitions)
- **mb_X**: Number of background mismatches (all other errors)
- **Use case**: Ancient DNA analysis with dense genome coverage

#### Format 3: Sparse Standard Format
Header line with genome names, followed by data rows with alternating genome names and mismatch counts for ONLY genomes with alignments:
```
read_id	total_count	Genome_A	Genome_B	Genome_C	Genome_D	Genome_E
read_id_1	75	Genome_B	4	Genome_C	2	Genome_D	1
read_id_2	100	Genome_A	5	Genome_B	2	Genome_E	1
```

The integer after each genome name is the number of mismatches

- **Use with**: Methods A, AE, B, BE, BEfull
- **Format**: Header + alternating genome names and mismatch counts for ONLY aligned genomes
- **Use case**: Large datasets with sparse alignments (memory efficient)

#### Format 4: Sparse Damage Format
Header line with genome names, followed by data rows with genome names and three values (nd, md, mb) for ONLY genomes with alignments:
```
read_id	total_count	Genome_A	Genome_B	Genome_C	Genome_D	Genome_E
read_id_1	75	Genome_B	2	1	1	Genome_C	2	0	4	Genome_D	3	1	0
read_id_2	100	Genome_A	2	0	1	Genome_B	4	0	2	Genome_E	3	2	1
```
The three integers after each genome name are: [nd] [md] [mb]

- **Use with**: Methods C, CE, CED, CEfull, CEDfull
- **Format**: Header + genome names followed by nd, md, mb values for ONLY aligned genomes
- **Use case**: Large ancient DNA datasets with sparse alignments (memory efficient)

### Genome Name Translation

MLID supports multiple ways to handle genome name conversion for compatibility with different data sources:

#### Automatic Genome Key File Detection
When processing data from tools like BAMreader or TaxIdent that use short genome identifiers (0, 1, 2, ...), MLID can automatically detect and use a `.genome_key` file in the same directory as the input file to translate these to full genome names.

```bash
# Auto-detect genome key file (looks for input_file.genome_key)
./mlid -i data_with_short_names.txt -M BEfull

# Specify custom genome key file
./mlid -i data.txt -M BEfull -k /path/to/custom_genome_key.txt

# Force use of genome key file (error if not found)
./mlid -i data.txt -M BEfull -K

# Disable genome key detection (use short names as-is)
./mlid -i data.txt -M BEfull -N
```

**Genome key file format:**
```
0	Escherichia_coli_K12
1	Salmonella_enterica
2	Bacillus_subtilis
```

#### Legacy Name Conversion
For backward compatibility with older workflows, you can convert genome names using a mapping file:

```bash
# Convert names using a mapping file
./mlid -i data.txt -M BEfull -c -g mapping_file.txt
```

**Mapping file format:**
```
OldGenomeName1	NewGenomeName1
OldGenomeName2	NewGenomeName2
```

### Output Format

The program generates a single standardized output file:
- **Default naming**: Input `mydata.txt` → Output `mydata_out.txt`
- **Custom naming**: `-o /path/to/results` → Output `/path/to/results_out.txt`
- **Directory support**: `-o /results/myanalysis` → Output `/results/myanalysis_out.txt`

**Output Filtering**: Only genomes/taxonomic groups with proportion > 1e-10 are included in the output (pruned entries are omitted for large datasets).

#### Output File Structure

The output file contains two sections:

**1. Header Section (metadata):**
```
Model=AE
LogLikelihood=-12345.678901
SQUAREM=Yes
Iterations=25
Converged=Yes
ComputationTime=3.45
```

**2. Data Section (tab-separated):**
```
genome	proportion	error_rate	damage_rate
A	0.800000	0.000099	not_applicable
B	0.100000	0.000099	not_applicable
C	0.100000	0.000099	not_applicable
```

**Note**: Only genomes/taxonomic groups with significant proportions (> 1e-10) are included. For large datasets with dynamic pruning, this filters out pruned entries automatically.

### Parameter Constraint Files

Parameter constraint files allow you to set custom bounds on error and damage rates for specific taxonomic groups during the M-step of the EM algorithm. This is useful for incorporating biological knowledge or technical constraints into the optimization process.

#### Constraint File Format

Each line in the constraint file specifies a constraint for one taxonomic group:

```
# Comment lines start with #
species < 0.001                    # Upper bound: species error rate must be < 0.001
genus > 0.01 < 0.1                # Range constraint: genus error rate between 0.01 and 0.1
family < 0.05                     # Upper bound: family error rate must be < 0.05
```

#### Constraint Types

**Upper bound constraint:**
```
taxonomic_group < upper_value
```

**Range constraint:**
```  
taxonomic_group > lower_value < upper_value
```

#### Example Constraint File

```
# MLIDConstraints.txt - Example constraint file
# Constrain species-level error rates to be very low (high specificity)
species < 0.001

# Genus-level rates should be moderate  
genus > 0.010 < 0.100

# Family-level rates can be higher but still bounded
family < 0.200
```

#### Usage with Constraints

```bash
# Apply constraints during CEDfull analysis
./mlid -i damage_data.txt -M CEDfull -C MLIDConstraints.txt -v

# Constraints work with any method that estimates error rates
./mlid -i data.txt -M BEfull -C constraints.txt
```

#### Constraint Application

- Constraints are applied during the M-step of each EM iteration
- If an unconstrained parameter estimate violates a constraint, it is clamped to the constraint boundary
- Verbose output (`-v`) shows both unconstrained and constrained values for debugging
- Constraints are case-insensitive for taxonomic group names

### Dirichlet Prior Weights

The Dirichlet prior weights feature allows you to incorporate prior beliefs about the mixture proportions into the EM optimization. This is particularly useful for regularizing estimates when you have limited data or want to encode biological expectations about relative abundances at different taxonomic levels.

#### Mathematical Background

Without priors, the M-step updates proportions using the formula:
```
π_k = (1/N) Σᵢ γᵢₖ
```

With Dirichlet priors, the M-step becomes:
```
π_k = (Σᵢ γᵢₖ + α_k - 1) / (N + Σⱼ α_j - K)
```

where:
- `α_k` is the Dirichlet prior weight for entity k
- `N` is the number of reads
- `K` is the number of entities (genomes + taxonomic groups)
- `γᵢₖ` are the posterior weights from the E-step

Higher α values pull proportions toward larger values (stronger prior belief in presence), while α=1 is equivalent to no prior (uniform).

#### Weights File Format

Each line specifies the alpha parameter for one taxonomic level:

```
# Comment lines start with #
species 1                   # α=1 for species (equivalent to no prior)
subspecies 1                # α=1 for subspecies
genus 10                    # α=10 for genus (10x stronger prior than species)
all_other 100               # α=100 for all unspecified taxonomic levels
```

#### Taxonomic Level Assignment

- **Genomes** (entities before '|' separator): Always treated as "species" level
- **Taxonomic groups** (entities after '|' separator): Level extracted from rank field (genus, family, order, etc.)
- **Unspecified levels**: Use "all_other" weight as fallback
- **Missing "all_other"**: Defaults to α=1 (uniform prior)
- **Matching**: Case-insensitive

#### Example Weights File

```
# MLIDweights.txt - Example Dirichlet prior weights
# Species-level (individual genomes) get minimal prior weight
species 1

# Subspecies same as species
subspecies 1

# Genus-level groups get moderate prior weight
# Encourages genus-level groups to have non-zero proportions
genus 10

# All other taxonomic levels (family, order, class, etc.)
# get strong prior weight
all_other 100
```

#### Usage with Weights

```bash
# Apply Dirichlet prior weights during CEDfull analysis
./mlid -i data.txt -M CEDfull -W MLIDweights.txt -v

# Combine with parameter constraints
./mlid -i ancient_dna.txt -M CEDfull -C constraints.txt -W weights.txt

# Weights only apply to methods with taxonomy data
# (Currently only CEDfull with taxonomic groups in input)
./mlid -i taxonomic_data.txt -M CEDfull -W weights.txt
```

#### When to Use Dirichlet Priors

**Use Dirichlet priors when:**
- You have limited data and want to regularize proportion estimates
- You want to encode prior beliefs about relative abundances
- Higher taxonomic groups should be favored over individual genomes
- You want to prevent overfitting to noise in sparse data

**Effect of prior weights:**
- **α=1**: No effect (equivalent to no prior)
- **α>1**: Encourages non-zero proportions, prevents aggressive pruning
- **Large α**: Strong prior belief that the entity is present
- **α<1**: Sparsity-inducing prior (favors smaller proportions)

**Important note on α<1:** When using sparsity-inducing priors (α < 1), SQUAREM acceleration may occasionally extrapolate to negative proportion values, which violate the mathematical constraint π > 0 required by the Dirichlet distribution. The implementation handles this by clamping negative values to a small positive constant (10⁻¹⁰) and renormalizing, which ensures numerical stability but may not yield the exact MAP estimate in all cases.

**Behavior without -W flag:**
- All entities use α=1 (uniform prior)
- Equivalent to standard EM without priors
- Backward compatible with previous versions

#### Parameter Columns

**error_rate column:**
- Models A, B, C: `fixed=0.000100` (the value provided with -e)
- Models AE, BE, CE, CED: Single estimated value (same for all genomes)
- Models BEfull, CEfull, CEDfull: Different value per genome

**damage_rate column:**
- Models A, AE, B, BE, BEfull: `not_applicable`
- Model C: `fixed=0.010000` (the value provided with -r)
- Models CE, CEfull: `fixed=0.010000` (the value provided with -r)
- Model CED: Single estimated value (same for all genomes)
- Model CEDfull: Different value per genome

### Debug Mode (-D)

The debug mode allows you to load previously estimated parameters and calculate the log likelihood without running the EM algorithm. This is useful for:

- Verifying parameter loading and likelihood calculation
- Testing specific parameter combinations
- Debugging convergence issues

#### Usage

```bash
# Debug mode with CEDfull
./mlid -i data.txt -M CEDfull -D parameter_file.txt -v
```

#### Parameter File Format

The debug mode reads the same output format produced by mlid:

```
Model=CEDfull
LogLikelihood=-857500.629941
SQUAREM=Yes
Iterations=10
Converged=Yes
ComputationTime=0.78

genome	proportion	error_rate	damage_rate
A	0.603974	0.009987	0.010116
D	0.099891	0.009928	0.009964
C	0.099479	0.009982	0.019830
B	0.096143	0.009962	0.010101
T1:genus	0.100514	0.020226	0.008867
T5:genus	0.000000	0.013809	0.029617
```

#### Behavior

- Loads all parameters (proportions, error rates, damage rates) from the file
- Recalculates the log likelihood using these exact parameters
- Prints the likelihood and exits immediately
- Does NOT run EM optimization

### Likelihood Ratio Test (-T)

The likelihood ratio test allows you to test the null hypothesis that a specific genome or taxonomic group has zero proportion in the sample. This provides statistical evidence for the presence/absence of particular organisms.

#### Usage

```bash
# Likelihood ratio test for a genome marked with 'T' in the parameter file
./mlid -i data.txt -M CEDfull -T parameter_file_with_test_marker.txt -C constraints.txt
```

#### Requirements

- **Only supported with CEDfull method** (`-M CEDfull`)
- Parameter file must contain a 'T' marker identifying the test genome/group
- Program will exit with error if used with other methods

#### Parameter File Format

The parameter file is identical to the debug format, but with a 'T' marker after the genome/taxonomic group to be tested:

```
Model=CEDfull
LogLikelihood=-857500.629941
SQUAREM=Yes
Iterations=10
Converged=Yes
ComputationTime=0.78

genome	proportion	error_rate	damage_rate
A	0.603974	0.009987	0.010116
D	0.099891	0.009928	0.009964
C	0.099479	0.009982	0.019830
B	0.096143	0.009962	0.010101
T1:genus	0.100514	0.020226	0.008867 T
T5:genus	0.000000	0.013809	0.029617
```

Note the 'T' marker after `T1:genus` - this identifies it as the genome/taxonomic group to test.

#### Test Procedure

1. **Load L_A**: Reads the original log likelihood (L_A) from the parameter file header
2. **Remove test genome**: Identifies the genome/group marked with 'T' and removes it from analysis
3. **Renormalize**: Renormalizes remaining proportions to sum to 1.0
4. **Optimize L_0**: Runs full EM optimization to find optimal log likelihood under null hypothesis
5. **Calculate ratio**: Computes likelihood ratio = L_A - L_0 and displays results

#### Output

```
=== Likelihood Ratio Test Results ===
Test genome/taxonomic group: T1:genus
Log likelihood under alternative hypothesis (L_A): -857500.629941
Log likelihood under null hypothesis (L_0): -902884.797818
Log likelihood ratio: L_A - L_0 = -857500.629941 - -902884.797818 = 45384.167877
```

#### Interpretation

- **Positive likelihood ratio**: Evidence against null hypothesis (genome is present)
- **Larger values**: Stronger evidence against null hypothesis
- **Values >> 10**: Very strong evidence for presence of the tested genome/group

#### Example Applications

```bash
# Test if T1:genus has non-zero proportion
./mlid -i ancient_dna.txt -M CEDfull -T results_with_T_marker.txt -C constraints.txt

# Test a specific genome (mark it with 'T' in the parameter file first)
./mlid -i metagenomic_data.txt -M CEDfull -T test_parameters.txt
```

### Examples

#### Standard Per-Genome Analysis
```bash
# BEfull with standard format - estimates K error rates
./cem -i format1_dense_standard.txt -M BEfull -v

# BEfull with sparse format - more memory efficient
./cem -i format3_sparse_standard.txt -M BEfull -o results_BEfull
```

#### Damage Per-Genome Analysis
```bash
# CEfull: K error rates, fixed damage rate
./cem -i format2_dense_damage.txt -M CEfull -r 0.01

# CEDfull: K error rates + K damage rates (most complex model)
./cem -i format4_sparse_damage.txt -M CEDfull

# CEDfull with joint damage rate (single shared damage for all entities)
./cem -i format4_sparse_damage.txt -M CEDfull -j

# CEDfull with verbose output to see per-genome rate estimates
./cem -i format2_dense_damage.txt -M CEDfull -v
```

#### Comparative Analysis
```bash
# Compare standard BE vs BEfull
./cem -i data.txt -M BE -o results_BE
./cem -i data.txt -M BEfull -o results_BEfull

# Compare damage models CE vs CEfull vs CEDfull
./cem -i damage_data.txt -M CE -o results_CE
./cem -i damage_data.txt -M CEfull -o results_CEfull
./cem -i damage_data.txt -M CEDfull -o results_CEDfull
```

#### Data Format Control
```bash
# Force sparse processing for large datasets
./mlid -i huge_dataset.txt -M BEfull -p

# Force dense processing for small datasets with complete alignments
./mlid -i small_dataset.txt -M BEfull -d
```

#### Genome Name Translation
```bash
# Auto-detect genome key file for short names from BAMreader/TaxIdent
./mlid -i bamreader_output.txt -M CEDfull

# Use specific genome key file
./mlid -i taxident_output.txt -M BEfull -k /path/to/genome.key

# Force genome key usage (fail if not found)
./mlid -i data.txt -M BEfull -K

# Convert old genome names to new ones using mapping file
./mlid -i legacy_data.txt -M AE -c -g name_mapping.txt
```

#### Advanced Analysis
```bash
# Apply parameter constraints during optimization
./mlid -i ancient_dna.txt -M CEDfull -C constraints.txt -v

# Debug mode: load parameters and calculate likelihood
./mlid -i data.txt -M CEDfull -D previous_results_out.txt

# Likelihood ratio test for presence of specific genome
./mlid -i data.txt -M CEDfull -T test_parameters.txt
```

### Model Selection Guidelines

#### When to Use Full Models

**Use BEfull when:**
- You suspect different reference genomes have different error patterns
- Dataset is large enough to estimate K error rates reliably (recommended: ≥100 reads per genome)
- You need detailed error characterization for each reference genome

**Use CEfull when:**
- Working with ancient DNA with consistent damage patterns across genomes
- Different genomes show varying background error rates but similar damage
- You want to estimate per-genome error rates while keeping damage rate fixed

**Use CEDfull when:**
- Working with ancient DNA where both error and damage patterns vary per genome
- Different genomes may have been processed differently or have different preservation states
- You have sufficient data to estimate both K error rates and K damage rates (recommended: ≥200 reads per genome)

#### Parameter Count Considerations

For K reference genomes:
- **Standard models (A, AE, B, BE)**: K-1 to K parameters
- **BEfull**: 2K-1 parameters (K-1 proportions + K error rates)
- **Damage models (C, CE, CED)**: K-1 to K+1 parameters  
- **CEfull**: 2K parameters (K-1 proportions + K error rates + 1 damage rate)
- **CEDfull**: 3K-1 parameters (K-1 proportions + K error rates + K damage rates)

**Data requirements scale with parameter count** - ensure you have sufficient data for reliable estimation.

## Performance Considerations

### SQUAREM Acceleration
- **Default enabled**: Provides ~3-5x faster convergence for most models
- **Automatic fallback**: Uses standard EM steps when acceleration would increase objective
- **Can be disabled**: Use `-S` or `--no-squarem` if needed for debugging

### Two-Phase Optimization with Pruning
When both SQUAREM acceleration and dynamic pruning are enabled (default behavior), CEMfull automatically uses a two-phase optimization strategy:

**Phase 1: Standard EM with Pruning (first 10 iterations)**
- Uses standard EM algorithm (non-accelerated)
- Performs dynamic pruning at scheduled iterations (1,2,3,4,5,10...)
- Removes low-proportion genomes and orphaned reads
- Provides stable initial convergence with data structure changes

**Phase 2: SQUAREM without Pruning (remaining iterations)**
- Switches to SQUAREM acceleration for faster convergence
- Pruning is disabled to maintain fixed parameter space
- Uses accelerated optimization on the pruned dataset
- Provides rapid final convergence

This approach solves the mathematical incompatibility between SQUAREM's fixed parameter space requirement and dynamic pruning's changing data structures. The result is both fast convergence and effective genome pruning.

### Model Complexity
- **Standard < BEfull**: BEfull requires more computation per iteration
- **CE < CEfull < CEDfull**: Increasing complexity for damage models
- **Full models**: Benefit significantly from SQUAREM acceleration

### Memory Usage
- **Per-genome rates**: Additional K or 2K parameters stored in memory
- **Sparse formats**: Still provide significant memory savings for full models
- **Large K**: Memory usage scales linearly with number of genomes

## Troubleshooting

### Common Issues with Full Models

**Convergence problems:**
- SQUAREM acceleration is enabled by default; try `-S` to disable if experiencing issues
- Full models may require more iterations: use `-m 2000` or higher
- Lower tolerance may be needed: try `-t 1e-8`
- Insufficient data per genome can cause instability

**Parameter estimation issues:**
- Some genomes may have very few reads, leading to unreliable rate estimates
- Use filtering (`-f`) to remove low-proportion genomes before analysis
- Consider whether you have enough data for the full model complexity

**Memory issues:**
- Full models require more memory for per-genome rate storage
- Use sparse format when possible for large datasets
- Consider subsetting your data for initial exploratory analysis

### Getting Help
```bash
# Show detailed help including full models
./cem --help

# Run with verbose output for debugging
./cem -i data.txt -M BEfull -v

# Test with a small dataset first
./cem -i small_test.txt -M BEfull -m 100 -v
```

## File Structure
```
CEMfull/
├── src/                 # Source files
│   ├── main.c           # Main program with full model support
│   ├── em_algorithms.c  # EM algorithm implementations including full models
│   ├── io_utils.c       # Input/output with format detection
│   ├── memory.c         # Memory management
│   └── ...
├── include/             # Header files
│   ├── em_types.h       # Extended data types for per-genome rates
│   ├── em_algorithms.h  # Algorithm interfaces
│   └── ...
├── build/              # Build directory (generated)
├── Makefile            # Build system
└── README.md           # This file
```

## License

????