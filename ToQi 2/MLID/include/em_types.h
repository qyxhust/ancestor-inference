#ifndef EM_TYPES_H
#define EM_TYPES_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <getopt.h>

// Genome mapping structure
typedef struct genome_mapping {
    char *compressed_id;  // Short genome name (e.g., "G1")
    char *full_name;      // Full genome name (e.g., "AACE03000010.1")
} genome_mapping_t;

// Hash table for fast genome name lookup
typedef struct genome_hash_entry {
    char *compressed_id;
    char *full_name;
    struct genome_hash_entry *next;  // For collision chaining
} genome_hash_entry_t;

typedef struct genome_hash_table {
    genome_hash_entry_t **buckets;
    int bucket_count;
} genome_hash_table_t;

// Machine precision constants
#define EM_DOUBLE_MIN 2.2250738585072014e-308
#define EM_DOUBLE_MAX 1.7976931348623157e+308
#define EM_TOLERANCE 1e-6
#define EM_MAX_ITER 1000

// Method types
typedef enum {
    METHOD_A,       // A: Single-rate model (no binomial coefficient)
    METHOD_AE,      // AE: Single-rate model (estimated rate)
    METHOD_B,       // B: Single-rate model (error rate / 3)
    METHOD_BE,      // BE: Single-rate model (error rate / 3, estimated rate)
    METHOD_C,       // C: Damage model (fixed damage and background rates)
    METHOD_CE,      // CE: Damage model (fixed damage rate, estimated background rate)
    METHOD_CED,     // CED: Damage model (estimated damage and background rates)
    METHOD_BEFULL,  // BEfull: Per-genome error rates (standard model)
    METHOD_CEFULL,  // CEfull: Per-genome error rates (damage model, fixed damage)
    METHOD_CEDFULL, // CEDfull: Per-genome error and damage rates
    // METHOD_FULL removed - non-functional extended model
    METHOD_INVALID  // Invalid method (for error handling)
} em_method_t;

// Data format type
typedef enum {
    FORMAT_DENSE,   // Traditional dense matrix format
    FORMAT_SPARSE,  // New sparse format
    FORMAT_AUTO     // Auto-detect format
} data_format_t;

// Taxonomic group structure
typedef struct {
    int taxid;              // NCBI taxonomic ID (e.g., 2759)
    char *rank;             // Taxonomic rank (e.g., "domain", "species", etc.)
    char *name;             // Full taxonomic name (T{taxid}:{rank})
} taxonomic_group_t;

// Parameter constraint types
typedef enum {
    CONSTRAINT_INEQUALITY,  // Range constraint with min/max bounds
    CONSTRAINT_EQUALITY     // Fixed value constraint
} constraint_type_t;

// Parameter constraint structure
typedef struct {
    char *entity_type;      // "species" or taxonomic rank (e.g., "genus", "family")
    constraint_type_t type; // Type of constraint (INEQUALITY or EQUALITY)
    double min_bound;       // Minimum value (-1.0 if no lower bound, for INEQUALITY only)
    double max_bound;       // Maximum value (-1.0 if no upper bound, for INEQUALITY only)
    double fixed_value;     // Fixed value (for EQUALITY only)
} parameter_constraint_t;

// Dirichlet prior weight structure
typedef struct {
    char *taxonomic_level;  // Taxonomic level (e.g., "species", "genus", "family", "all_other")
    double weight;          // Alpha parameter for Dirichlet prior
} dirichlet_weight_t;

// Taxonomic data for reads (sparse format)
typedef struct {
    int n_taxa;                    // Number of taxonomic groups in this read
    int *taxon_indices;           // Array of indices into the global taxonomic groups array
    
    // For damage models (5 values per taxon): nd, md, mds, mb, mbs
    short *nd_values;             // Damage sites (nd)
    short *md_values;             // Damage mismatches to ALL genomes (md)
    short *mds_values;            // Damage mismatches to SOME genomes (mds)
    short *mb_values;             // Background mismatches to ALL genomes (mb) - integer format
    short *mbs_values;            // Background mismatches to SOME genomes (mbs) - integer format
    
    // New float format support (for new sparse formats)
    float *mb_values_float;       // Background mismatches to ALL genomes (mb) - float format
    float *mbs_values_float;      // Background mismatches to SOME genomes (mbs) - float format (unused in new format)
    
    // For standard models (2 values per taxon): mb+md, mbs+mds
    short *total_all_values;      // Total mismatches to ALL genomes (mb+md) - integer format
    short *total_some_values;     // Total mismatches to SOME genomes (mbs+mds) - integer format (unused in new format)
    float *total_all_values_float; // Total mismatches to ALL genomes (mb+md) - float format
} read_taxonomic_data_t;

// Main data structure for EM algorithm
typedef struct {
    short **n_matrix;      // Sites matrix [n_reads x n_genomes] - using short for memory efficiency
    short **d_matrix;      // Mismatches matrix [n_reads x n_genomes] - using short for memory efficiency
    double *proportions;   // Mixture proportions [n_genomes]
    double error_rate;     // Error rate parameter (background rate for damage models)
    int n_reads;          // Number of reads
    int n_genomes;        // Number of genomes
    
    // Damage model specific data
    short **nd_matrix;     // Damage sites matrix [n_reads x n_genomes] (for damage models) - using short
    short **md_matrix;     // Damage errors matrix [n_reads x n_genomes] (for damage models) - using short
    short **mb_matrix;     // Background errors matrix [n_reads x n_genomes] (for damage models) - using short
    double damage_rate;    // Damage error rate (for damage models)
    
    // Per-genome error rates (for *full models)
    double *error_rates;    // Per-genome error rates [n_genomes] (for *full models)
    double *damage_rates;   // Per-genome damage rates [n_genomes] (for CEDfull model)
    
    // Working matrices
    double **Q_matrix;     // Likelihood matrix [n_reads x n_genomes]
    double **W_matrix;     // Posterior weights [n_reads x n_genomes]
    
    // Convergence tracking
    double *prev_proportions;   // Previous iteration proportions
    double prev_error_rate;     // Previous iteration error rate
    double prev_damage_rate;    // Previous iteration damage rate
    double *prev_error_rates;   // Previous iteration per-genome error rates [n_genomes]
    double *prev_damage_rates;  // Previous iteration per-genome damage rates [n_genomes]
    double log_likelihood;     // Current log-likelihood
    double prev_log_likelihood; // Previous log-likelihood
    int iteration;             // Current iteration number
    
    // Dynamic pruning (for dense format)
    int *active_genomes;       // Boolean array: 1 if genome is active, 0 if pruned [n_genomes]
    int *active_reads;         // Boolean array: 1 if read is active, 0 if removed [n_reads]
    int n_active_genomes;      // Current number of active genomes
    int n_active_reads;        // Current number of active reads
    char **read_names;         // Read names for tracking removed reads
    FILE *removed_reads_file;  // File handle for writing removed read names
    double min_proportion;     // Threshold for pruning genomes
    
    // Taxonomic groups support (only for sparse format with higher taxa)
    int has_taxonomic_groups;         // Flag: 1 if input has taxonomic groups, 0 if not
    taxonomic_group_t *taxonomic_groups;  // Array of all taxonomic groups
    int n_taxonomic_groups;           // Number of taxonomic groups
    read_taxonomic_data_t *read_taxonomic_data; // Per-read taxonomic data [n_reads]
    char **genome_names;              // Array of genome names [n_genomes]
} em_data_t;

// Configuration structure
typedef struct {
    char *input_file;      // Input mismatch matrix file
    char *output_prefix;   // Output file prefix  
    em_method_t method;    // Method to use
    double error_rate;     // Initial error rate (background rate for damage models)
    double damage_rate;    // Initial damage rate (for damage models)
    int joint_damage_rate; // 0=separate damage rates (default), 1=joint/shared
    double tolerance;      // Convergence tolerance
    int max_iterations;    // Maximum EM iterations
    int verbose;           // Verbose output flag
    double min_proportion; // Minimum proportion threshold
    int use_squarem;       // Use SQUAREM acceleration (0=standard EM, 1=SQUAREM) [default: 1]
    data_format_t data_format; // Data format (dense, sparse, auto-detect)
    int force_sparse;      // Force sparse processing even for dense format
    int force_dense;       // Force dense processing even for sparse format
    int convert_names;     // Flag to enable name conversion (0=off, 1=on)
    char *mapping_file;    // Genome mapping file for name conversion
    void *genome_mappings; // Loaded genome mappings (genome_mapping_t*)
    int n_mappings;        // Number of mappings loaded
    char *genome_key_file; // Explicit genome key file path
    int use_genome_key;    // 0=disabled, 1=auto-detect (default), 2=force use
    char *constraints_file; // Parameter constraint file path
    parameter_constraint_t *constraints; // Parsed constraints array
    int n_constraints;     // Number of constraints
    char *weights_file;    // Dirichlet prior weights file path
    dirichlet_weight_t *weights; // Parsed weights array
    int n_weights;         // Number of weights
} em_config_t;

// Sparse alignment data structure (CSR format)
typedef struct {
    int *row_ptr;       // [n_reads + 1] - start index for each read
    int *col_indices;   // [nnz] - genome indices for each alignment  
    short *n_values;    // [nnz] - site counts for aligned pairs (using short for memory efficiency)
    short *d_values;    // [nnz] - mismatch counts for aligned pairs
    // Damage model specific data
    short *nd_values;   // [nnz] - damage site counts for aligned pairs (for damage models)
    short *md_values;   // [nnz] - damage error counts for aligned pairs (for damage models)
    short *mb_values;   // [nnz] - background error counts for aligned pairs (for damage models)
    int n_reads;        // Number of reads
    int n_genomes;      // Number of genomes  
    int nnz;            // Number of non-zero (aligned) entries
    
    // Taxonomic groups support
    int has_taxonomic_groups;         // Flag: 1 if input has taxonomic groups, 0 if not
    taxonomic_group_t *taxonomic_groups;  // Array of all taxonomic groups
    int n_taxonomic_groups;           // Number of taxonomic groups
    read_taxonomic_data_t *read_taxonomic_data; // Per-read taxonomic data [n_reads]
    
    // New format detection
    int use_float_mb_format;          // Flag: 1 if using new float mb format, 0 if using integer format
} sparse_alignment_t;

// Sparse EM data structure
typedef struct {
    sparse_alignment_t *alignment_data; // Sparse alignment matrix
    double *proportions;   // Mixture proportions [n_genomes]
    double error_rate;     // Error rate parameter (background rate for damage models)
    double damage_rate;    // Damage error rate (for damage models)
    double *error_rates;   // Per-genome error rates (for BEfull/CEfull/CEDfull)
    double *damage_rates;  // Per-genome damage rates (for CEDfull)
    double shared_damage_rate;     // Shared damage rate (used if joint_damage_rate=1)
    int joint_damage_rate;         // 0=separate, 1=joint (copied from config)
    int n_reads;          // Number of reads
    int n_genomes;        // Number of genomes
    
    // Working sparse arrays
    double *Q_sparse;      // Likelihood values [nnz]
    double *W_sparse;      // Posterior weights [nnz]
    
    // Convergence tracking
    double *prev_proportions;  // Previous iteration proportions
    
    // Dynamic pruning
    int *active_genomes;       // Boolean array: 1 if genome is active, 0 if pruned [n_genomes]
    int *active_reads;         // Boolean array: 1 if read is active, 0 if removed [n_reads]
    int n_active_genomes;      // Current number of active genomes
    int n_active_reads;        // Current number of active reads
    char **read_names;         // Read names for tracking removed reads
    FILE *removed_reads_file;  // File handle for writing removed read names
    double min_proportion;     // Threshold for pruning genomes
    double prev_error_rate;    // Previous iteration error rate
    double prev_damage_rate;   // Previous iteration damage rate
    double *prev_error_rates;  // Previous iteration per-genome error rates
    double *prev_damage_rates; // Previous iteration per-genome damage rates
    double log_likelihood;     // Current log-likelihood
    double prev_log_likelihood; // Previous log-likelihood
    int iteration;             // Current iteration number
    
    // Pre-calculated values for extended models with taxonomic groups
    double log_half;           // log(0.5) - constant for Class A and damage mixing
    double log_p_mix;          // log(p_mix) for Class B - recalculated when error_rate changes
    double *log_p_mix_per_genome; // log(p_mix) for per-genome models [n_genomes] - recalculated when error_rates change
    
    // Extended storage for taxonomic groups (sparse format: only store non-zero taxonomic alignments)
    int n_taxonomic_groups;    // Number of taxonomic groups (for easy access)
    int n_taxonomic_entries;   // Total number of non-zero taxonomic entries across all reads
    double *Q_taxonomic_sparse; // Q values for taxonomic groups [n_taxonomic_entries] - sparse storage
    double *W_taxonomic_sparse; // W values for taxonomic groups [n_taxonomic_entries] - sparse storage
    int *taxonomic_row_ptr;    // Start index for each read's taxonomic data [n_reads + 1] - CSR format
    int *taxonomic_col_indices; // Taxonomic group indices for each entry [n_taxonomic_entries] - CSR format
    double *taxonomic_proportions; // Proportions for taxonomic groups [n_taxonomic_groups]
    double *prev_taxonomic_proportions; // Previous proportions for taxonomic groups
    
    // Extended model rate arrays removed - unified approach used instead
    
    // New format detection and unified genome+taxonomic arrays
    int use_float_mb_format;          // Flag: 1 if using new float mb format, 0 if using integer format
    int n_total_entities;             // Total entities: n_genomes + n_taxonomic_groups
    double *unified_proportions;      // Unified proportions [n_genomes + n_taxonomic_groups]
    double *unified_error_rates;      // Unified error rates [n_genomes + n_taxonomic_groups] (for CEDfull)
    double *unified_damage_rates;     // Unified damage rates [n_genomes + n_taxonomic_groups] (for CEDfull)
    double *prev_unified_proportions; // Previous unified proportions
    double *prev_unified_error_rates; // Previous unified error rates
    double *prev_unified_damage_rates; // Previous unified damage rates
    int *unified_active_entities;     // Boolean array: 1 if entity is active, 0 if pruned [n_total_entities]
    int n_active_entities;            // Current number of active entities (genomes + taxonomic groups)
    
    // Parameter constraints
    parameter_constraint_t *constraints; // Pointer to constraints (from config)
    int n_constraints;                   // Number of constraints

    // Fixed parameter masks for equality constraints (error rates only)
    int *error_rate_fixed_mask;      // 1=fixed, 0=free [n_genomes + n_taxonomic_groups]
    double *error_rate_fixed_values; // Fixed values for equality-constrained error rates

    // Dirichlet prior weights for proportions
    dirichlet_weight_t *weights;     // Pointer to weights (from config)
    int n_weights;                   // Number of weights
} sparse_em_data_t;

// Results structure
typedef struct {
    double *final_proportions;  // Final estimated proportions
    double final_error_rate;    // Final estimated error rate (background rate for damage models)
    double final_damage_rate;   // Final estimated damage rate (for damage models)
    double *final_error_rates;  // Final per-genome error rates (for *full models)
    double *final_damage_rates; // Final per-genome damage rates (for CEDfull model)
    double final_log_likelihood; // Final log-likelihood
    int converged;              // Convergence flag
    int iterations;             // Number of iterations taken
    double computation_time;    // Total computation time in seconds
    double *genome_names;       // Optional: genome names for output
    
    // Taxonomic group results
    double *final_taxonomic_proportions; // Final taxonomic group proportions
    double *final_taxonomic_error_rates;  // Final per-taxonomic group error rates (CEDfull only)
    double *final_taxonomic_damage_rates; // Final per-taxonomic group damage rates (CEDfull only)
    taxonomic_group_t *taxonomic_groups; // Taxonomic group information (names, ranks)
    int n_taxonomic_groups;     // Number of taxonomic groups
} em_results_t;

#endif // EM_TYPES_H