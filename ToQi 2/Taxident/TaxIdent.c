/*
 * TaxIdent - Taxonomic Hierarchy BAM Processor
 *  
 * Author: Claude Code Assistant and Rasmus Nielsen
 * Date: 2025
 *
 *If you are actually reading this code - good luck! This is an expriment in vibe coding - and it shows. 
 * 
 * Dependencies: HTSlib
 * Compile: gcc -O3 -o TaxIdent TaxIdent.c -lhts -lz -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <pthread.h>
#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <htslib/sam.h>

#ifdef WITH_HTSLIB
#include <htslib/sam.h>
#include <htslib/hts.h>
#endif

#define INITIAL_MAX_GENOMES 200
#define GENOME_GROWTH_FACTOR 2
#define MAX_NAME_LEN 256
#define MAX_BAM_FILES 1000
#define MIN_READ_LENGTH 20
#define MAX_TAXID 10000000  // Maximum expected taxid value
#define MAX_RANK_LEN 32

// Parsimony state constants (2-bit encoding: 00=unset, 01=state_0, 10=state_1, 11=ambiguous)
#define PARSIMONY_UNSET 0
#define PARSIMONY_STATE_0 1
#define PARSIMONY_STATE_1 2
#define PARSIMONY_AMBIGUOUS 3

// Taxonomy data structures
typedef struct TaxNode {
    int taxid;
    int parent_taxid;
    char rank[MAX_RANK_LEN];
    int *leaf_genomes;      // Array of genome indices
    int n_leaves;
    int leaves_capacity;
    struct TaxNode **children;
    int n_children;
    int children_capacity;
    int is_active;          // Whether this node appears in our data
} TaxNode;

typedef struct {
    TaxNode **nodes;        // Array of taxonomy nodes indexed by taxid
    int max_taxid;
    int n_active_nodes;     // Number of nodes actually present in data
    int *active_taxids;     // List of active taxids for iteration
} TaxonomyTree;

// Per-position mismatch information
typedef struct {
    uint8_t *has_mismatch;     // Bit vector: 1 if genome has mismatch at this position
    uint8_t *is_damage;        // Bit vector: 1 if mismatch is damage type
    int n_genomes;
} PositionMismatch;

// Per-read mismatch storage
typedef struct {
    char *read_id;
    int read_length;
    int trimmed_length;
    PositionMismatch *positions;  // Array of position-specific mismatch info
    int *aligned_genomes;          // Which genomes have alignments
    int *genome_nd_values;         // nd value for each aligned genome
    int *genome_nm_values;         // NM value (total mismatches) for each aligned genome
    int n_aligned;
    int capacity;
    
} ReadMismatches;

// Hash table for read mismatches
typedef struct ReadMismatchEntry {
    ReadMismatches *data;
    struct ReadMismatchEntry *next;
} ReadMismatchEntry;

// Lightweight tree node for per-read taxonomy trees (linear algorithm)
typedef struct TreeNode {
    int taxid;
    char rank[MAX_RANK_LEN];
    struct TreeNode *parent;
    struct TreeNode **children;
    int n_children;
    int max_children;  // Allocated capacity for children array
    
    // Position-wise bit vectors for efficient combination (reusable across reads)
    uint8_t *md_vec;    // Damage mismatches to ALL descendants
    uint8_t *mds_vec;   // Damage mismatches to SOME descendants  
    uint8_t *mb_vec;    // Background mismatches to ALL descendants
    uint8_t *mbs_vec;   // Background mismatches to SOME descendants
    uint8_t nd_value;   // Max damage sites (scalar)
    
    // Parsimony state vectors (2 bits per position: 00=unset, 01=0, 10=1, 11=0|1)
    uint8_t *parsimony_state_vec;  // 2-bit parsimony states
    uint8_t *damage_parsimony_state_vec;  // 2-bit parsimony states for damage mismatches

    // Genome index for leaf nodes (-1 for internal nodes)
    int genome_idx;
} TreeNode;

// Program options
typedef struct {
    char *input_file;   // Input file (BAM or text)
    char *genome_list;  // File listing genome names (for BAM mode)
    char *output_file;
    int min_mapq;
    long long max_reads;  // Use long long to handle billions of reads
    int verbose;
    int bam_mode;       // 0 = text mode (default), 1 = BAM mode
    int use_simple_mode; // Use original unoptimized algorithms
    char ignore_char;   // Character to ignore prefixes until (-I flag)
    int damage_sites;   // Number of damage-susceptible sites at read ends (default 5)
    int enable_damage;  // Enable damage analysis
    char *output_format; // Output format: dense, sparse, dense_damage, sparse_damage
    int short_names;    // Use short genome names (G1, G2, etc.) instead of full names
    int create_key_file; // Create key file mapping short names to full names
    int num_threads;    // Number of decompression threads
    int silent;         // Suppress all progress output
    int asymmetric_damage; // Use asymmetric damage calculation (C->T in first s, G->A in last s)
    int skip_indels;    // Skip alignments containing indels
    int precise_indels; // Use precise calculation for alignments with indels (slower)
    int compress_output; // Use compressed read IDs and genome IDs
    char *taxid_file;   // NCBI names.dmp file for taxid mapping
    char *genome_map_file; // Output file for genome ID to name mapping
    int parallel_files; // Number of BAM files to process in parallel (directory mode)
    int use_rg_tag;     // Use RG tag for genome identification instead of reference name (default: 0)
    int use_dense;      // Use dense internal data structures instead of sparse (default: 0, sparse is default)
    int no_redundancy;  // Assume reads don't appear in multiple BAM files (allows optimization)
    
    // New taxonomy-related options
    char *taxonomy_dir;     // Directory containing NCBI taxonomy files
    char *acc2taxid_file;   // Accession to taxid mapping file
    int with_higher_taxa;   // Include higher taxonomic groups in output
    int higher_taxa_with_mbs; // Use old ALL/SOME method (5 values) instead of parsimony (3 values)
    int consecutive_mode;   // Optimize for consecutive alignments
    char *tax_levels;       // Comma-separated list of taxonomic levels to include
    int test_tree;          // Test tree parsing and exit
    int enforce_dense_strict; // Only output reads with alignments to ALL genomes/taxids/readgroups
    int use_penalty_mode;   // 0 = use -1 for missing (default), 1 = use max+1 penalty
    int consolidate_by_taxid; // Consolidate genomes by taxid without taxonomy tree (sparse mode only)
} options_t;

// Genome names - now dynamic
char (*genome_names)[MAX_NAME_LEN] = NULL;  // Will be allocated dynamically
int n_genomes = 0;
int max_genomes_allocated = 0;  // Track current allocation size

// Thread safety for dynamic genome operations
pthread_mutex_t genome_mutex = PTHREAD_MUTEX_INITIALIZER;

// Tracking for unmapped accessions
long long invalid_taxid_alignments = 0;
FILE *unmapped_report_file = NULL;
char *global_output_filename = NULL;

// Compressed ID mapping - now dynamic
char (*genome_compressed_ids)[32] = NULL;  // Will be allocated dynamically
typedef struct {
    char full_name[MAX_NAME_LEN];
    char taxid[32];  // NCBI taxid if available, otherwise custom ID like A1, A2, etc.
} genome_mapping_t;
genome_mapping_t *genome_mappings = NULL;  // Will be allocated dynamically
int custom_id_counter = 1;  // For generating A1, A2, etc.

// Linear algorithm: Node pool for per-read taxonomy trees
static TreeNode *node_pool = NULL;
static int node_pool_capacity = 0;
static int node_pool_used = 0;
static uint8_t *shared_bit_vectors = NULL;  // Reusable bit vector memory
static int max_read_length = 200;  // Will grow as needed

// Reusable node_map to eliminate 80MB allocation per read
static TreeNode **reusable_node_map = NULL;
static int node_map_allocated = 0;

// Tracking for used taxids to optimize clearing
static int *used_taxids = NULL;  // Array to store which taxids were used
static int n_used_taxids = 0;    // Number of taxids used in current read
static int used_taxids_capacity = 0;  // Capacity of the array

// Hash table for genome discovery (eliminate O(n²) search)
#define DISCOVERY_HASH_SIZE 8192  // 8K buckets for genome discovery
typedef struct discovery_hash_entry {
    char name[MAX_NAME_LEN];
    struct discovery_hash_entry *next;
} discovery_hash_entry_t;

// Hash table for fast genome lookup
#define HASH_TABLE_SIZE 1048576  // 1M buckets - power of 2 for fast modulo, good balance for 100K-1.5M genomes
typedef struct genome_hash_entry {
    char name[MAX_NAME_LEN];
    int index;
    struct genome_hash_entry *next;
} genome_hash_entry_t;

genome_hash_entry_t *genome_hash_table[HASH_TABLE_SIZE];

// Character type lookup table for optimized MD parsing
#define CHAR_DIGIT 1
#define CHAR_ALPHA 2
#define CHAR_OTHER 0
static unsigned char char_types[256] = {0};

// Damage statistics structure - UPDATED for taxonomic analysis
typedef struct {
    int nd;  // Number of damage-susceptible sites (C or G in first/last s positions)
    int md;  // Number of damage mismatches (C->T or G->A in damage sites)
    int mb;  // Number of background mismatches (all other errors)
    // New fields for taxonomic analysis
    int mbs; // Background mismatches to SOME (not all) genomes
    int mds; // Damage mismatches to SOME (not all) genomes
} damage_stats_t;

// Global taxonomy variables
TaxonomyTree *taxonomy_tree = NULL;
int using_taxid_mode = 0;  // Global flag set by --with-higher-taxa option
ReadMismatchEntry **read_mismatch_table = NULL;  // Hash table for read mismatches
int read_mismatch_table_size = 1048576;  // 1M buckets

// Accession to taxid mapping
typedef struct AccTaxidEntry {
    char accession[MAX_NAME_LEN];
    int taxid;
    struct AccTaxidEntry *next;
} AccTaxidEntry;

AccTaxidEntry **acc2taxid_table = NULL;
int acc2taxid_table_size = 4194304;  // 4M buckets for large mapping files

// Read data structure for tracking best alignments across BAM files
typedef struct {
    char read_id[MAX_NAME_LEN];
    int read_length;
    int *mismatches;  // Will be allocated based on max_genomes_allocated
    int *indel_count;  // Number of indels in each alignment
    int *has_alignment;
    int max_mismatches;
    damage_stats_t *damage_stats;  // Damage statistics per genome
} read_data_t;

// ============ SPARSE MODE DATA STRUCTURES ============
// Single genome alignment entry for sparse mode (linked list node)
// OPTIMIZED: Store genome name directly instead of index to avoid lookups
typedef struct genome_alignment {
    char genome_name[MAX_NAME_LEN]; // Store full genome name directly
    int mismatches;                 // Mismatch count for this genome
    int indel_count;                // Indel count
    damage_stats_t damage;          // Damage statistics
    struct genome_alignment *next;  // Next in linked list
} genome_alignment_t;

// Memory pool for sparse mode allocations
#define MEMORY_POOL_BLOCK_SIZE 1000
typedef struct sparse_memory_pool_block {
    genome_alignment_t nodes[MEMORY_POOL_BLOCK_SIZE];
    int used;
    struct sparse_memory_pool_block *next;
} sparse_memory_pool_block_t;

typedef struct sparse_memory_pool {
    sparse_memory_pool_block_t *current_block;
    sparse_memory_pool_block_t *first_block;
    pthread_mutex_t mutex;
} sparse_memory_pool_t;

// Global memory pool for sparse mode
sparse_memory_pool_t *sparse_global_pool = NULL;

// Global options pointer for access in various functions
options_t *global_opts = NULL;

// ============ SPARSE MODE MEMORY POOL FUNCTIONS ============
sparse_memory_pool_t* create_sparse_memory_pool() {
    sparse_memory_pool_t *pool = malloc(sizeof(sparse_memory_pool_t));
    if (!pool) return NULL;
    
    pool->first_block = malloc(sizeof(sparse_memory_pool_block_t));
    if (!pool->first_block) {
        free(pool);
        return NULL;
    }
    
    pool->first_block->used = 0;
    pool->first_block->next = NULL;
    pool->current_block = pool->first_block;
    pthread_mutex_init(&pool->mutex, NULL);
    
    return pool;
}

genome_alignment_t* sparse_pool_alloc_node(sparse_memory_pool_t *pool) {
    pthread_mutex_lock(&pool->mutex);
    
    // Check if current block has space
    if (pool->current_block->used >= MEMORY_POOL_BLOCK_SIZE) {
        // Need new block
        sparse_memory_pool_block_t *new_block = malloc(sizeof(sparse_memory_pool_block_t));
        if (!new_block) {
            pthread_mutex_unlock(&pool->mutex);
            return NULL;
        }
        new_block->used = 0;
        new_block->next = NULL;
        pool->current_block->next = new_block;
        pool->current_block = new_block;
    }
    
    genome_alignment_t *node = &pool->current_block->nodes[pool->current_block->used++];
    pthread_mutex_unlock(&pool->mutex);
    
    // Initialize the node
    memset(node, 0, sizeof(genome_alignment_t));
    return node;
}

void reset_sparse_memory_pool(sparse_memory_pool_t *pool) {
    pthread_mutex_lock(&pool->mutex);
    
    // Just reset the used counters, keep all blocks allocated
    // This avoids malloc/free overhead
    sparse_memory_pool_block_t *block = pool->first_block;
    while (block) {
        block->used = 0;
        block = block->next;
    }
    pool->current_block = pool->first_block;
    
    pthread_mutex_unlock(&pool->mutex);
}

void free_sparse_memory_pool(sparse_memory_pool_t *pool) {
    if (!pool) return;
    
    // Free all blocks except the first
    sparse_memory_pool_block_t *block = pool->first_block->next;
    while (block) {
        sparse_memory_pool_block_t *next = block->next;
        free(block);
        block = next;
    }
    
    free(pool->first_block);
    pthread_mutex_destroy(&pool->mutex);
    free(pool);
}

// Forward declarations for functions used in sparse mode
unsigned int hash_string(const char *str);
int add_genome_dynamically(const char *genome_name);
const char *get_output_genome_name_compressed(int genome_index, options_t *opts);
void write_genome_mapping_file(const char *filename);
void write_short_names_mapping_file(const char *filename);
const char *extract_genome_name(const char *name, char ignore_char);
#ifdef WITH_HTSLIB
int auto_discover_genomes(const char *input_path, char ignore_char, int silent, int use_rg_tag, long long max_reads);
int calculate_trimmed_mismatches(bam1_t *read, sam_hdr_t *header, options_t *opts, damage_stats_t *damage_stats, int *indel_count);
#endif

// Taxonomy function declarations
TaxonomyTree* load_taxonomy_tree(const char *nodes_file, const char *names_file);
void free_taxonomy_tree(TaxonomyTree *tree);
int load_acc2taxid_mapping(const char *acc2taxid_file);
void free_acc2taxid_mapping(void);
int get_taxid_for_accession(const char *accession);
void map_genomes_to_taxids(void);
void calculate_taxonomic_mismatches(ReadMismatches *read_data, TaxNode *node, damage_stats_t *stats);
void process_taxonomic_groups(options_t *opts);
ReadMismatches* get_or_create_read_mismatches(const char *read_id);
void store_read_mismatch_position(ReadMismatches *read_data, int position, int genome_idx, int has_mismatch, int is_damage);
int finalize_read_mismatches(ReadMismatches *read_data);
void free_read_mismatches(ReadMismatches *data);
void output_with_taxonomy(FILE *fp, options_t *opts);
void print_genome_taxonomy_tree(void);

// Linear algorithm function declarations
void init_node_pool(void);
TreeNode* get_node_from_pool(int taxid, const char *rank);
void reset_node_pool(void);
TreeNode* build_read_taxonomy_tree(ReadMismatches *read_data);
void calculate_node_vectors_postorder(TreeNode *node, ReadMismatches *read_data, int read_length);
void output_with_taxonomy_linear(FILE *fp, options_t *opts);

// Parsimony algorithm function declarations
int get_parsimony_state(uint8_t *state_vec, int pos);
void set_parsimony_state(uint8_t *state_vec, int pos, int state);
int count_parsimony_state_1(uint8_t *state_vec, int read_length);
void calculate_parsimony_postorder(TreeNode *node, ReadMismatches *read_data, int read_length);
void resolve_parsimony_preorder(TreeNode *node, int read_length);
void debug_print_parsimony_states(TreeNode *node, ReadMismatches *read_data, const char *read_id, int read_length);

// Global instrumentation for find_or_add_genome
static long find_genome_calls = 0;
static double find_genome_time = 0;
static long genome_comparisons = 0;

// Helper function for sparse mode - find or add a genome using hash table
int find_or_add_genome(const char *genome_name) {
    struct timeval start, end;
    gettimeofday(&start, NULL);
    find_genome_calls++;
    
    // Use hash table for O(1) lookup instead of O(n) linear search
    unsigned int hash = hash_string(genome_name);
    genome_hash_entry_t *entry = genome_hash_table[hash];
    
    while (entry) {
        genome_comparisons++;
        if (strcmp(genome_name, entry->name) == 0) {
            gettimeofday(&end, NULL);
            find_genome_time += (end.tv_sec - start.tv_sec) * 1000000.0 + 
                               (end.tv_usec - start.tv_usec);
            return entry->index;
        }
        entry = entry->next;
    }
    
    // Not found, add it dynamically
    int result = add_genome_dynamically(genome_name);
    gettimeofday(&end, NULL);
    find_genome_time += (end.tv_sec - start.tv_sec) * 1000000.0 + 
                       (end.tv_usec - start.tv_usec);
    return result;
}

// Include the optimized temp file processor
#include "process_temp_file.c"

// Wrapper to match function signature and return genome count
int rewrite_sparse_with_header_optimized_wrapper(const char *temp_file, const char *output_file, options_t *opts) {
    int is_damage_format = opts->enable_damage && 
                          (strcmp(opts->output_format, "dense_damage") == 0 || 
                           strcmp(opts->output_format, "sparse_damage") == 0);
    
    return rewrite_sparse_with_header_optimized(temp_file, output_file, is_damage_format);
}

// Old version kept for compatibility
void rewrite_sparse_with_header(const char *temp_file, const char *output_file, options_t *opts) {
    FILE *temp_fp = fopen(temp_file, "r");
    if (!temp_fp) {
        fprintf(stderr, "Error: Cannot open temp file %s\n", temp_file);
        return;
    }
    
    FILE *out_fp = fopen(output_file, "w");
    if (!out_fp) {
        fprintf(stderr, "Error: Cannot create output file %s\n", output_file);
        fclose(temp_fp);
        return;
    }
    
    // Write header with discovered genomes
    int is_damage_format = opts->enable_damage && 
                          (strcmp(opts->output_format, "dense_damage") == 0 || 
                           strcmp(opts->output_format, "sparse_damage") == 0);
    
    fprintf(out_fp, "read_id\ttotal_count");
    for (int i = 0; i < n_genomes; i++) {
        const char *name = get_output_genome_name_compressed(i, opts);
        if (is_damage_format) {
            fprintf(out_fp, "\t%s_nd\t%s_md\t%s_mb", name, name, name);
        } else {
            fprintf(out_fp, "\t%s", name);
        }
    }
    fprintf(out_fp, "\n");
    
    // Copy data from temp file
    // Use dynamic buffer for potentially very long lines in dense format
    size_t buffer_size = 1048576; // 1MB buffer for dense format lines
    char *line = malloc(buffer_size);
    if (!line) {
        fprintf(stderr, "Error: Cannot allocate line buffer\n");
        fclose(temp_fp);
        fclose(out_fp);
        return;
    }
    
    while (fgets(line, buffer_size, temp_fp)) {
        fputs(line, out_fp);
    }
    
    free(line);
    
    fclose(temp_fp);
    fclose(out_fp);
}

// Helper for dense mode - find genome index from name
static int find_genome_index_by_name(const char *genome_name) {
    // For dense mode, we still need to track genomes
    // This adds overhead but is required for dense format
    return find_or_add_genome(genome_name);
}

// OPTIMIZED: Find or add genome alignment using genome name directly
genome_alignment_t* find_or_add_genome_alignment_optimized(genome_alignment_t **head, 
                                                          const char *genome_name, 
                                                          sparse_memory_pool_t *pool) {
    genome_alignment_t *current = *head;
    
    // Search for existing alignment by name
    while (current) {
        if (strcmp(current->genome_name, genome_name) == 0) {
            return current;
        }
        current = current->next;
    }
    
    // Not found, add new node
    genome_alignment_t *new_node = sparse_pool_alloc_node(pool);
    if (!new_node) return NULL;
    
    strncpy(new_node->genome_name, genome_name, MAX_NAME_LEN - 1);
    new_node->genome_name[MAX_NAME_LEN - 1] = '\0';
    new_node->mismatches = -1;  // Initialize to -1 to indicate no alignment yet
    new_node->indel_count = 0;
    memset(&new_node->damage, 0, sizeof(damage_stats_t));
    new_node->next = *head;
    *head = new_node;
    
    return new_node;
}

// Removed old cached function that uses genome_idx - no longer needed
// The optimized approach stores genome names directly without lookups

#ifdef WITH_HTSLIB
// OPTIMIZED Sparse implementation - writes full genome names to avoid lookups
void process_bam_file_sparse_optimized(options_t *opts) {
    // Open BAM file
    samFile *bam_fp = sam_open(opts->input_file, "r");
    if (!bam_fp) {
        fprintf(stderr, "Error: Cannot open BAM file: %s\n", opts->input_file);
        exit(1);
    }
    
    sam_hdr_t *header = sam_hdr_read(bam_fp);
    if (!header) {
        fprintf(stderr, "Error: Cannot read BAM header from: %s\n", opts->input_file);
        sam_close(bam_fp);
        exit(1);
    }
    
    // Set up multi-threading for decompression if requested
    if (opts->num_threads > 0) {
        hts_set_threads(bam_fp, opts->num_threads);
    }
    
    // Initialize sparse memory pool
    sparse_global_pool = create_sparse_memory_pool();
    if (!sparse_global_pool) {
        fprintf(stderr, "Error: Cannot create memory pool\n");
        sam_hdr_destroy(header);
        sam_close(bam_fp);
        exit(1);
    }
    
    // Note: Auto-discovery is already done in main() before calling this function
    // So we don't need to do it again here
    
    // Open output file
    FILE *out_fp = NULL;
    char temp_filename[1024];
    int using_temp_file = 0;
    
    // Determine output format flags
    int is_damage_format = opts->enable_damage && 
                          (strcmp(opts->output_format, "dense_damage") == 0 || 
                           strcmp(opts->output_format, "sparse_damage") == 0);
    int is_sparse = strcmp(opts->output_format, "sparse") == 0 || 
                   strcmp(opts->output_format, "sparse_damage") == 0;
    
    // OPTIMIZED: Use temp file for both sparse and dense formats
    // This allows genome discovery before writing the header
    snprintf(temp_filename, sizeof(temp_filename), "%s.tmp", opts->output_file);
    out_fp = fopen(temp_filename, "w");
    using_temp_file = 1;
    
    if (!out_fp) {
        fprintf(stderr, "Error: Cannot create output file: %s\n", 
                using_temp_file ? temp_filename : opts->output_file);
        sam_hdr_destroy(header);
        sam_close(bam_fp);
        exit(1);
    }
    
    // OPTIMIZED: Don't write header when using temp file (sparse mode)
    if (!using_temp_file) {
        if (is_damage_format) {
            fprintf(out_fp, "read_id\ttotal_count");
            for (int i = 0; i < n_genomes; i++) {
                const char *name = get_output_genome_name_compressed(i, opts);
                if (is_sparse) {
                    fprintf(out_fp, "\t%s_nd\t%s_md\t%s_mb", name, name, name);
                } else {
                    fprintf(out_fp, "\t%s_nd\t%s_md\t%s_mb", name, name, name);
                }
            }
            fprintf(out_fp, "\n");
        } else {
            fprintf(out_fp, "read_id\ttotal_count");
            for (int i = 0; i < n_genomes; i++) {
                const char *name = get_output_genome_name_compressed(i, opts);
                fprintf(out_fp, "\t%s", name);
            }
            fprintf(out_fp, "\n");
        }
    }
    
    // Data structures for sparse tracking
    char current_read[MAX_NAME_LEN] = "";
    genome_alignment_t *current_alignments = NULL;
    int read_length = 0;
    int max_mismatches = 0;
    int reads_processed = 0;
    // OPTIMIZED: No cache needed since we're not doing genome lookups
    long long alignments_processed = 0;
    int unique_reads = 0;
    
    bam1_t *read = bam_init1();
    
    // Report file processing start (unless silent)
    if (!opts->silent) {
        if (opts->consolidate_by_taxid) {
            fprintf(stderr, "\n=== TAXID CONSOLIDATION MODE ===\n");
            fprintf(stderr, "Processing mode: TAXID CONSOLIDATION (--consolidate-by-taxid)\n");
            fprintf(stderr, "Note: This mode consolidates genomes by taxonomic ID.\n");
            fprintf(stderr, "      Output columns represent unique taxids (T12345).\n");
            fprintf(stderr, "      Multiple genomes mapping to same taxid are consolidated.\n");
            fprintf(stderr, "=========================================\n");
        } else {
            fprintf(stderr, "\n=== GENOME-BASED PROCESSING MODE ===\n");
            fprintf(stderr, "Processing mode: GENOME-BASED (no --with-higher-taxa)\n");
            fprintf(stderr, "Note: This mode processes individual genomes separately.\n");
            fprintf(stderr, "      Output columns represent individual genome files.\n");
            fprintf(stderr, "      Use --with-higher-taxa for taxid-based taxonomic analysis.\n");
            fprintf(stderr, "=========================================\n");
        }
        fprintf(stderr, "Processing BAM file (sparse mode): %s\n", opts->input_file);
    }
    
    while (sam_read1(bam_fp, header, read) >= 0) {
        // Check limits
        if (opts->max_reads > 0 && alignments_processed >= opts->max_reads) {
            break;
        }
        
        alignments_processed++;
        
        // Progress reporting every 1,000,000 alignments (unless silent)
        if (!opts->silent && alignments_processed > 0 && (alignments_processed % 1000000 == 0)) {
            fprintf(stderr, "Processed %lld alignments (%d unique reads discovered)\n", 
                    alignments_processed, reads_processed);
        }
        
        // Skip unmapped reads
        if (read->core.flag & BAM_FUNMAP) continue;
        
        // Skip low quality alignments
        if (read->core.qual < opts->min_mapq) continue;
        
        char *read_name = bam_get_qname(read);
        
        // Check if this is a new read
        if (strcmp(read_name, current_read) != 0) {
            // Output previous read if it exists
            if (current_read[0] != '\0' && current_alignments != NULL) {
                // Generate compressed read ID if needed
                char output_read_id[MAX_NAME_LEN];
                if (opts->compress_output) {
                    snprintf(output_read_id, MAX_NAME_LEN, "R%d", ++unique_reads);
                } else {
                    strcpy(output_read_id, current_read);
                }
                
                // Output the read with proper format
                if (is_sparse) {
                    // OPTIMIZED: Write full genome names - will be replaced in final processing
                    fprintf(out_fp, "%s\t%d", output_read_id, read_length);
                    genome_alignment_t *align = current_alignments;
                    while (align) {
                        // Write full genome name directly
                        if (is_damage_format) {
                            fprintf(out_fp, "\t%s\t%d\t%d\t%d",
                                    align->genome_name,
                                    align->damage.nd, align->damage.md, align->damage.mb);
                        } else {
                            fprintf(out_fp, "\t%s\t%d", align->genome_name, align->mismatches);
                        }
                        align = align->next;
                    }
                    fprintf(out_fp, "\n");
                } else {
                    // Dense format - output all genomes
                    fprintf(out_fp, "%s\t%d", output_read_id, read_length);
                    
                    // Build dense arrays from sparse list
                    int *dense_mismatches = calloc(n_genomes, sizeof(int));
                    damage_stats_t *dense_damage = calloc(n_genomes, sizeof(damage_stats_t));
                    int *has_alignment = calloc(n_genomes, sizeof(int));
                    
                    // Initialize mismatches to -1 for no alignment
                    for (int i = 0; i < n_genomes; i++) {
                        dense_mismatches[i] = -1;
                        dense_damage[i].nd = -1;
                        dense_damage[i].md = -1;
                        dense_damage[i].mb = -1;
                    }
                    
                    genome_alignment_t *align = current_alignments;
                    while (align) {
                        // For dense mode, we need to find genome index from name
                        int genome_idx = find_genome_index_by_name(align->genome_name);
                        if (genome_idx >= 0) {
                            dense_mismatches[genome_idx] = align->mismatches;
                            dense_damage[genome_idx] = align->damage;
                            has_alignment[genome_idx] = 1;
                        }
                        align = align->next;
                    }
                    
                    // Set penalty for non-aligned genomes
                    int penalty = max_mismatches + 1;
                    if (penalty > read_length) penalty = read_length;
                    
                    for (int i = 0; i < n_genomes; i++) {
                        if (is_damage_format) {
                            if (has_alignment[i]) {
                                fprintf(out_fp, "\t%d\t%d\t%d", 
                                        dense_damage[i].nd, dense_damage[i].md, dense_damage[i].mb);
                            } else {
                                fprintf(out_fp, "\t-1\t-1\t-1");
                            }
                        } else {
                            if (has_alignment[i]) {
                                fprintf(out_fp, "\t%d", dense_mismatches[i]);
                            } else {
                                fprintf(out_fp, "\t-1");
                            }
                        }
                    }
                    fprintf(out_fp, "\n");
                    
                    free(dense_mismatches);
                    free(dense_damage);
                    free(has_alignment);
                }
                
                reads_processed++;
            }
            
            // Reset for new read - OPTIMIZED: No genome tracking overhead
            strcpy(current_read, read_name);
            current_alignments = NULL;  // Just reset pointer, nodes will be reused from pool
            if (opts->enable_damage) {
                read_length = read->core.l_qseq;  // No trim when tracking damage
            } else if (opts->damage_sites > 0) {
                read_length = read->core.l_qseq - 2 * opts->damage_sites;  // Trim
                if (read_length < 0) read_length = 0;
            } else {
                read_length = read->core.l_qseq;  // No trim
            }
            max_mismatches = 0;
            
            // Reset memory pool for each read to prevent memory growth
            reset_sparse_memory_pool(sparse_global_pool);
        }
        
        // Get genome identifier
        const char* genome_identifier = NULL;
        if (opts->use_rg_tag) {
            uint8_t *rg_tag = bam_aux_get(read, "RG");
            if (rg_tag) {
                genome_identifier = bam_aux2Z(rg_tag);
            }
        } else {
            genome_identifier = sam_hdr_tid2name(header, read->core.tid);
        }
        
        if (!genome_identifier) continue;
        
        // Process genome name with ignore character if specified
        char processed_genome_name[MAX_NAME_LEN];
        if (opts->consolidate_by_taxid) {
            // Convert genome to taxid format for consolidation
            const char *ref_name = sam_hdr_tid2name(header, read->core.tid);

            // Extract base accession (reuse logic from find_or_add_taxid)
            const char *genome_name = extract_genome_name(ref_name, opts->ignore_char);
            if (!genome_name) continue;

            char base_accession[MAX_NAME_LEN];
            strncpy(base_accession, genome_name, MAX_NAME_LEN - 1);
            base_accession[MAX_NAME_LEN - 1] = '\0';
            char *dot = strchr(base_accession, '.');
            if (dot) *dot = '\0';

            // Look up taxid
            int taxid = get_taxid_for_accession(base_accession);
            if (taxid > 0) {
                snprintf(processed_genome_name, MAX_NAME_LEN, "T%d", taxid);
            } else {
                // No taxid found - use original genome name
                strncpy(processed_genome_name, genome_name, MAX_NAME_LEN - 1);
                processed_genome_name[MAX_NAME_LEN - 1] = '\0';
            }
        } else {
            // Original genome name processing
            if (opts->ignore_char != '\0') {
                const char *last_sep = strrchr(genome_identifier, opts->ignore_char);
                if (last_sep && *(last_sep + 1) != '\0') {
                    strncpy(processed_genome_name, last_sep + 1, MAX_NAME_LEN - 1);
                } else {
                    strncpy(processed_genome_name, genome_identifier, MAX_NAME_LEN - 1);
                }
            } else {
                strncpy(processed_genome_name, genome_identifier, MAX_NAME_LEN - 1);
            }
            processed_genome_name[MAX_NAME_LEN - 1] = '\0';
        }
        
        // OPTIMIZED: No genome lookup - just use the name directly
        genome_alignment_t *align = find_or_add_genome_alignment_optimized(&current_alignments, 
                                                                          processed_genome_name, 
                                                                          sparse_global_pool);
        if (!align) continue;
        
        // Get NM value from BAM tag for alignment scoring (same as taxonomy mode)
        uint8_t *nm_tag = bam_aux_get(read, "NM");
        int nm_value = 0;
        if (nm_tag) {
            nm_value = bam_aux2i(nm_tag);
        }

        // Calculate damage statistics (but don't use the returned value for scoring)
        damage_stats_t current_damage_stats;
        int current_indel_count = 0;
        int calculated_nm = calculate_trimmed_mismatches(read, header, opts,
                                                   opts->enable_damage ? &current_damage_stats : NULL,
                                                   &current_indel_count);

        // Skip this alignment if damage calculation failed
        if (calculated_nm < 0) continue;

        // Update alignment data (keep best/minimum mismatches)
        if (align->mismatches == -1 || nm_value < align->mismatches) {
            align->mismatches = nm_value;
            align->indel_count = current_indel_count;
            if (opts->enable_damage) {
                align->damage = current_damage_stats;
            }
        }
        
        if (nm_value > max_mismatches) {
            max_mismatches = nm_value;
        }
    }
    
    // Output last read if it exists
    if (current_read[0] != '\0' && current_alignments != NULL) {
        char output_read_id[MAX_NAME_LEN];
        if (opts->compress_output) {
            snprintf(output_read_id, MAX_NAME_LEN, "R%d", ++unique_reads);
        } else {
            strcpy(output_read_id, current_read);
        }

        if (is_sparse) {
            fprintf(out_fp, "%s\t%d", output_read_id, read_length);
            genome_alignment_t *align = current_alignments;
            while (align) {
                // OPTIMIZED: Write full genome name
                if (is_damage_format) {
                    fprintf(out_fp, "\t%s\t%d\t%d\t%d", 
                            align->genome_name,
                            align->damage.nd, align->damage.md, align->damage.mb);
                } else {
                    fprintf(out_fp, "\t%s\t%d", align->genome_name, align->mismatches);
                }
                align = align->next;
            }
            fprintf(out_fp, "\n");
        }

        reads_processed++;
    }
    
    // Clean up
    bam_destroy1(read);
    fclose(out_fp);
    sam_hdr_destroy(header);
    sam_close(bam_fp);
    free_sparse_memory_pool(sparse_global_pool);
    
    // Handle temp file if used - OPTIMIZED version
    if (using_temp_file) {
        // Sparse mode only - dense output is not supported in sparse internal mode
        n_genomes = rewrite_sparse_with_header_optimized_wrapper(temp_filename, opts->output_file, opts);
        unlink(temp_filename);
    }
    
    // Final report
    if (!opts->silent) {
        fprintf(stderr, "\nCompleted processing:\n");
        fprintf(stderr, "  Total reads processed: %d\n", reads_processed);
        fprintf(stderr, "  Total alignments processed: %lld\n", alignments_processed);
        fprintf(stderr, "  Total genomes: %d\n", n_genomes);
        
        // Report find_or_add_genome performance
        if (find_genome_calls > 0) {
            fprintf(stderr, "\nPerformance Analysis - find_or_add_genome:\n");
            fprintf(stderr, "  Function calls: %ld\n", find_genome_calls);
            fprintf(stderr, "  Total time: %.2f ms\n", find_genome_time / 1000.0);
            fprintf(stderr, "  Avg time per call: %.3f µs\n", find_genome_time / find_genome_calls);
            fprintf(stderr, "  Avg comparisons: %.2f\n", (double)genome_comparisons / find_genome_calls);
        }
    }
    
    // Write mapping files if requested
    if (opts->compress_output && opts->genome_map_file) {
        write_genome_mapping_file(opts->genome_map_file);
    }
    
    // REMOVED: Writing genome key file here would overwrite the correct file from process_temp_file.c
    // In sparse mode, genome_names array is empty, but process_temp_file.c has already written
    // the correct genome key file with actual genome names from the temp file
}
#endif

// BAM files list
char bam_files[MAX_BAM_FILES][MAX_NAME_LEN];
int n_bam_files = 0;

// Suffix tree structure for fast read name lookups (using reversed names)
typedef struct suffix_tree_node {
    struct suffix_tree_node *children[256];  // Direct indexing by character
    read_data_t *read_data;  // Pointer to read data if this is a terminal node
} suffix_tree_node_t;

suffix_tree_node_t *suffix_tree_root = NULL;

// Function to reverse a string in place
void reverse_string(char *str, int len) {
    int i, j;
    char temp;
    for (i = 0, j = len - 1; i < j; i++, j--) {
        temp = str[i];
        str[i] = str[j];
        str[j] = temp;
    }
}

// Create a new suffix tree node
suffix_tree_node_t *create_suffix_tree_node(void) {
    suffix_tree_node_t *node = calloc(1, sizeof(suffix_tree_node_t));
    if (!node) {
        fprintf(stderr, "Error: Failed to allocate suffix tree node\n");
        exit(1);
    }
    return node;
}

// Insert a read into the suffix tree using reversed read name
void suffix_tree_insert(read_data_t *read) {
    if (!suffix_tree_root) {
        suffix_tree_root = create_suffix_tree_node();
    }
    
    // Create reversed copy of read name
    char reversed_name[MAX_NAME_LEN];
    strcpy(reversed_name, read->read_id);
    int len = strlen(reversed_name);
    reverse_string(reversed_name, len);
    
    // Navigate/build the tree
    suffix_tree_node_t *current = suffix_tree_root;
    for (int i = 0; i < len; i++) {
        unsigned char c = (unsigned char)reversed_name[i];
        if (!current->children[c]) {
            current->children[c] = create_suffix_tree_node();
        }
        current = current->children[c];
    }
    
    // Store read data at the terminal node
    current->read_data = read;
}

// Find a read in the suffix tree using reversed read name
read_data_t *suffix_tree_find(const char *read_name) {
    if (!suffix_tree_root) {
        return NULL;
    }
    
    // Create reversed copy of read name
    char reversed_name[MAX_NAME_LEN];
    strcpy(reversed_name, read_name);
    int len = strlen(reversed_name);
    reverse_string(reversed_name, len);
    
    // Navigate the tree
    suffix_tree_node_t *current = suffix_tree_root;
    for (int i = 0; i < len; i++) {
        unsigned char c = (unsigned char)reversed_name[i];
        if (!current->children[c]) {
            return NULL;  // Read not found
        }
        current = current->children[c];
    }
    
    return current->read_data;
}

// Free the suffix tree (recursive)
void suffix_tree_free(suffix_tree_node_t *node) {
    if (!node) return;
    
    for (int i = 0; i < 256; i++) {
        if (node->children[i]) {
            suffix_tree_free(node->children[i]);
        }
    }
    free(node);
}

// String interning for genome names to avoid repeated strcmp
typedef struct string_intern_entry {
    char *string;
    int id;
    struct string_intern_entry *next;
} string_intern_entry_t;

#define STRING_INTERN_SIZE 1024
string_intern_entry_t *string_intern_table[STRING_INTERN_SIZE];
int next_intern_id = 0;

// Intern a string and return its unique ID
int intern_string(const char *str) {
    unsigned int hash = 0;
    const char *p = str;
    while (*p) {
        hash = hash * 31 + *p++;
    }
    hash %= STRING_INTERN_SIZE;
    
    // Check if already interned
    string_intern_entry_t *entry = string_intern_table[hash];
    while (entry) {
        if (strcmp(entry->string, str) == 0) {
            return entry->id;
        }
        entry = entry->next;
    }
    
    // Add new entry
    entry = malloc(sizeof(string_intern_entry_t));
    entry->string = strdup(str);
    entry->id = next_intern_id++;
    entry->next = string_intern_table[hash];
    string_intern_table[hash] = entry;
    
    return entry->id;
}

// Batch processing structure for parallel processing
#define BATCH_SIZE 1000
typedef struct {
    bam1_t *reads[BATCH_SIZE];
    int count;
    int file_idx;
} read_batch_t;

// Memory pool for frequent allocations
typedef struct memory_pool {
    void *blocks[100];
    int block_count;
    size_t block_size;
    int current_block;
    size_t current_offset;
} memory_pool_t;

memory_pool_t *create_memory_pool(size_t block_size) {
    memory_pool_t *pool = calloc(1, sizeof(memory_pool_t));
    pool->block_size = block_size;
    pool->blocks[0] = malloc(block_size);
    pool->block_count = 1;
    return pool;
}

void *pool_alloc(memory_pool_t *pool, size_t size) {
    if (pool->current_offset + size > pool->block_size) {
        // Need new block
        if (pool->current_block + 1 >= pool->block_count) {
            if (pool->block_count >= 100) return NULL; // Pool full
            pool->blocks[pool->block_count++] = malloc(pool->block_size);
        }
        pool->current_block++;
        pool->current_offset = 0;
    }
    
    void *ptr = (char*)pool->blocks[pool->current_block] + pool->current_offset;
    pool->current_offset += size;
    return ptr;
}

void free_memory_pool(memory_pool_t *pool) {
    for (int i = 0; i < pool->block_count; i++) {
        free(pool->blocks[i]);
    }
    free(pool);
}

// Function prototypes
void print_usage(const char *prog_name);
int parse_options(int argc, char **argv, options_t *opts);
const char *get_output_genome_name(int genome_index, int use_short_names);
int load_genome_list(const char *filename, char ignore_char, int use_simple_mode);
const char *extract_genome_name(const char *name, char ignore_char);
int find_or_add_taxid(const char *ref_name, char ignore_char);
int get_taxid_from_index(int taxid_idx);
int find_genome_index(const char *ref_name, char ignore_char);
int find_bam_files(const char *directory);
int is_directory(const char *path);
void init_optimizations(void);
void init_genome_arrays(void);
void grow_genome_arrays(void);
unsigned int hash_string(const char *str);
void add_genome_to_hash(const char *name, int index);
int add_genome_dynamically(const char *genome_name);
int find_genome_index_fast(const char *ref_name, char ignore_char);
int find_genome_index_simple(const char *ref_name, char ignore_char);
int parse_md_tag_simple(const char *md_string, int *mismatch_positions, int max_mismatches);
#ifdef WITH_HTSLIB
int calculate_mismatches_with_indels(bam1_t *read, sam_hdr_t *header, options_t *opts, 
                                    damage_stats_t *damage_stats, const char *md_string);
int calculate_trimmed_mismatches(bam1_t *read, sam_hdr_t *header, options_t *opts, damage_stats_t *damage_stats, int *indel_count);
void process_bam_file(options_t *opts);
void process_bam_file_sparse(options_t *opts);
void process_bam_file_sparse_optimized(options_t *opts);
void process_bam_file_with_taxonomy(options_t *opts);
void process_multiple_bam_files(options_t *opts);
void process_multiple_bam_files_optimized(options_t *opts);
void process_multiple_bam_files_with_taxonomy(options_t *opts);
void process_multiple_bam_files_unified(options_t *opts);
int auto_discover_genomes(const char *input_path, char ignore_char, int silent, int use_rg_tag, long long max_reads);
#endif
int find_or_add_genome(const char *genome_name);
const char *get_output_genome_name_compressed(int genome_index, options_t *opts);
void write_genome_mapping_file(const char *filename);
void write_short_names_mapping_file(const char *filename);
void rewrite_sparse_with_header(const char *temp_file, const char *output_file, options_t *opts);
int is_damage_mismatch(char ref_base, char read_base);
int count_damage_sites(const char *ref_seq, int start, int end);
void print_read_data(const read_data_t *read_data, options_t *opts);
void process_text_file(options_t *opts);
void generate_compressed_read_id(int read_number, char *buffer);
void load_taxid_mapping(const char *taxid_file);
void generate_genome_compressed_id(int genome_index, const char *genome_name);
void write_genome_mapping_file(const char *filename);
void write_short_names_mapping_file(const char *filename);

// ============================================================================
// TAXONOMY IMPLEMENTATION FUNCTIONS
// ============================================================================

// Load taxonomy tree from NCBI nodes.dmp and names.dmp files
TaxonomyTree* load_taxonomy_tree(const char *nodes_file, const char *names_file) {
    TaxonomyTree *tree = (TaxonomyTree *)calloc(1, sizeof(TaxonomyTree));
    if (!tree) return NULL;
    
    FILE *fp = fopen(nodes_file, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open nodes file %s\n", nodes_file);
        free(tree);
        return NULL;
    }
    
    // First pass: find max taxid
    char line[4096];
    int max_taxid = 0;
    while (fgets(line, sizeof(line), fp)) {
        int taxid, parent_taxid;
        char rank[MAX_RANK_LEN];
        if (sscanf(line, "%d | %d | %31[^|]", &taxid, &parent_taxid, rank) == 3) {
            if (taxid > max_taxid) max_taxid = taxid;
        }
    }
    
    tree->max_taxid = max_taxid;
    tree->nodes = (TaxNode **)calloc(max_taxid + 1, sizeof(TaxNode *));
    if (!tree->nodes) {
        fclose(fp);
        free(tree);
        return NULL;
    }
    
    // Second pass: create nodes
    rewind(fp);
    while (fgets(line, sizeof(line), fp)) {
        int taxid, parent_taxid;
        char rank[MAX_RANK_LEN];
        if (sscanf(line, "%d | %d | %31[^|]", &taxid, &parent_taxid, rank) == 3) {
            TaxNode *node = (TaxNode *)calloc(1, sizeof(TaxNode));
            if (!node) continue;
            
            node->taxid = taxid;
            node->parent_taxid = parent_taxid;
            
            // Trim whitespace from rank
            char *p = rank;
            while (*p && isspace(*p)) p++;
            char *end = p + strlen(p) - 1;
            while (end > p && isspace(*end)) *end-- = '\0';
            strncpy(node->rank, p, MAX_RANK_LEN - 1);
            
            tree->nodes[taxid] = node;
        }
    }
    fclose(fp);
    
    // Build parent-child relationships
    for (int i = 1; i <= max_taxid; i++) {
        TaxNode *node = tree->nodes[i];
        if (!node || node->parent_taxid == node->taxid) continue;
        
        TaxNode *parent = tree->nodes[node->parent_taxid];
        if (parent) {
            if (!parent->children) {
                parent->children_capacity = 10;
                parent->children = (TaxNode **)calloc(parent->children_capacity, sizeof(TaxNode *));
            }
            if (parent->n_children >= parent->children_capacity) {
                parent->children_capacity *= 2;
                parent->children = (TaxNode **)realloc(parent->children, 
                    parent->children_capacity * sizeof(TaxNode *));
            }
            if (parent->children) {
                parent->children[parent->n_children++] = node;
            }
        }
    }
    
    return tree;
}

// Free taxonomy tree
void free_taxonomy_tree(TaxonomyTree *tree) {
    if (!tree) return;
    
    if (tree->nodes) {
        for (int i = 0; i <= tree->max_taxid; i++) {
            if (tree->nodes[i]) {
                free(tree->nodes[i]->leaf_genomes);
                free(tree->nodes[i]->children);
                free(tree->nodes[i]);
            }
        }
        free(tree->nodes);
    }
    
    free(tree->active_taxids);
    free(tree);
}

// Load accession to taxid mapping
int load_acc2taxid_mapping(const char *acc2taxid_file) {
    FILE *fp = fopen(acc2taxid_file, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open acc2taxid file %s\n", acc2taxid_file);
        return -1;
    }
    
    // Initialize hash table
    acc2taxid_table = (AccTaxidEntry **)calloc(acc2taxid_table_size, sizeof(AccTaxidEntry *));
    if (!acc2taxid_table) {
        fclose(fp);
        return -1;
    }
    
    char line[1024];
    // Skip header if present
    if (fgets(line, sizeof(line), fp)) {
        if (strncmp(line, "accession", 9) != 0) {
            rewind(fp);
        }
    }
    
    int count = 0;
    while (fgets(line, sizeof(line), fp)) {
        char accession[MAX_NAME_LEN], accession_version[MAX_NAME_LEN];
        int taxid;
        
        if (sscanf(line, "%s %s %d", accession, accession_version, &taxid) == 3 ||
            sscanf(line, "%s\t%s\t%d", accession, accession_version, &taxid) == 3) {
            
            AccTaxidEntry *entry = (AccTaxidEntry *)malloc(sizeof(AccTaxidEntry));
            if (!entry) continue;
            
            strncpy(entry->accession, accession, MAX_NAME_LEN - 1);
            entry->taxid = taxid;
            
            unsigned int hash = hash_string(accession) % acc2taxid_table_size;
            entry->next = acc2taxid_table[hash];
            acc2taxid_table[hash] = entry;
            
            count++;
        }
    }
    
    fclose(fp);
    return count;
}

// Free acc2taxid mapping
void free_acc2taxid_mapping(void) {
    if (!acc2taxid_table) return;
    
    for (int i = 0; i < acc2taxid_table_size; i++) {
        AccTaxidEntry *entry = acc2taxid_table[i];
        while (entry) {
            AccTaxidEntry *next = entry->next;
            free(entry);
            entry = next;
        }
    }
    
    free(acc2taxid_table);
    acc2taxid_table = NULL;
}

// Get taxid for an accession
int get_taxid_for_accession(const char *accession) {
    if (!acc2taxid_table) return -1;
    
    unsigned int hash = hash_string(accession) % acc2taxid_table_size;
    AccTaxidEntry *entry = acc2taxid_table[hash];
    
    while (entry) {
        if (strcmp(entry->accession, accession) == 0) {
            return entry->taxid;
        }
        entry = entry->next;
    }
    
    return -1;
}

// Map genomes to taxids - DISABLED for taxid-based processing
void map_genomes_to_taxids(void) {
    return;  // Function disabled for taxid-based processing
/*
    if (!genome_to_taxid) {
        genome_to_taxid = (int *)calloc(max_genomes_allocated, sizeof(int));
    }
    
    for (int i = 0; i < n_genomes; i++) {
        // Extract accession from genome name
        char accession[MAX_NAME_LEN];
        strncpy(accession, genome_names[i], MAX_NAME_LEN - 1);
        
        // Try to extract accession (handle various formats)
        char *dot = strchr(accession, '.');
        if (dot) *dot = '\0';
        
        int taxid = get_taxid_for_accession(accession);
        if (taxid > 0) {
            get_taxid_from_index(i) = taxid;
            
            // Mark this taxid and all ancestors as active
            if (taxonomy_tree && taxid <= taxonomy_tree->max_taxid && taxonomy_tree->nodes[taxid]) {
                TaxNode *node = taxonomy_tree->nodes[taxid];
                
                // Add this genome to this node AND all ancestor nodes
                while (node) {
                    // Add genome to this node's leaf_genomes
                    if (!node->leaf_genomes) {
                        node->leaves_capacity = 10;
                        node->leaf_genomes = (int *)calloc(node->leaves_capacity, sizeof(int));
                    }
                    if (node->n_leaves >= node->leaves_capacity) {
                        node->leaves_capacity *= 2;
                        node->leaf_genomes = (int *)realloc(node->leaf_genomes,
                            node->leaves_capacity * sizeof(int));
                    }
                    if (node->leaf_genomes) {
                        node->leaf_genomes[node->n_leaves++] = i;
                    }
                    
                    // Mark as active
                    node->is_active = 1;
                    
                    // Move to parent (unless we're at root)
                    if (node->parent_taxid != node->taxid) {
                        node = taxonomy_tree->nodes[node->parent_taxid];
                    } else {
                        break;
                    }
                }
            }
        } else {
            // Unknown taxid - this case handled in find_or_add_taxid
        }
    }
*/
}

// Get or create read mismatches entry
ReadMismatches* get_or_create_read_mismatches(const char *read_id) {
    if (!read_mismatch_table) {
        read_mismatch_table = (ReadMismatchEntry **)calloc(read_mismatch_table_size, 
                                                           sizeof(ReadMismatchEntry *));
    }
    
    unsigned int hash = hash_string(read_id) % read_mismatch_table_size;
    ReadMismatchEntry *entry = read_mismatch_table[hash];
    
    // Search for existing entry
    while (entry) {
        if (strcmp(entry->data->read_id, read_id) == 0) {
            return entry->data;
        }
        entry = entry->next;
    }
    
    // Create new entry
    ReadMismatches *data = (ReadMismatches *)calloc(1, sizeof(ReadMismatches));
    data->read_id = strdup(read_id);
    
    entry = (ReadMismatchEntry *)malloc(sizeof(ReadMismatchEntry));
    entry->data = data;
    entry->next = read_mismatch_table[hash];
    read_mismatch_table[hash] = entry;
    
    return data;
}

// Store mismatch position for a read - RIGHT-SIZED BIT VECTORS
void store_read_mismatch_position(ReadMismatches *read_data, int position, 
                                 int genome_idx, int has_mismatch, int is_damage) {
    // Bounds check - report errors instead of silently skipping
    if (position < 0) {
        fprintf(stderr, "FATAL ERROR: Negative position %d in store_read_mismatch_position()\n", position);
        fprintf(stderr, "  This indicates a bug in position calculation.\n");
        abort();
    }
    if (position >= read_data->read_length) {
        fprintf(stderr, "FATAL ERROR: Position %d exceeds read length %d\n", position, read_data->read_length);
        fprintf(stderr, "  This indicates read length calculation is incorrect.\n");
        abort();
    }
    
    if (!read_data->positions) {
        read_data->positions = (PositionMismatch *)calloc(read_data->read_length, 
                                                          sizeof(PositionMismatch));
    }
    
    PositionMismatch *pos = &read_data->positions[position];
    
    // Initialize or resize bit vectors to handle current number of aligned genomes
    if (!pos->has_mismatch) {
        // First allocation for this position
        int bytes_needed = (read_data->n_aligned + 7) / 8;
        pos->has_mismatch = (uint8_t *)calloc(bytes_needed, 1);
        pos->is_damage = (uint8_t *)calloc(bytes_needed, 1);

        
        // Check for allocation failure
        if (!pos->has_mismatch || !pos->is_damage) {
            fprintf(stderr, "FATAL: Allocation failed for %d bytes\n", bytes_needed);
            abort();
        }
        pos->n_genomes = read_data->n_aligned;
        
    } else if (pos->n_genomes < read_data->n_aligned) {
        // Need to resize - read has grown since position was allocated
        
        int new_bytes_needed = (read_data->n_aligned + 7) / 8;
        int old_bytes = (pos->n_genomes + 7) / 8;
        
        if (new_bytes_needed > old_bytes) {
            // Reallocate with larger size
            pos->has_mismatch = (uint8_t *)realloc(pos->has_mismatch, new_bytes_needed);
            pos->is_damage = (uint8_t *)realloc(pos->is_damage, new_bytes_needed);


            // Clear new bytes
            memset(pos->has_mismatch + old_bytes, 0, new_bytes_needed - old_bytes);
            memset(pos->is_damage + old_bytes, 0, new_bytes_needed - old_bytes);
        }
        pos->n_genomes = read_data->n_aligned;
    }
    
    // Use the existing lookup array pattern to find local bit index
    // This should be called with the lookup already set up in the calling function
    int local_bit_idx = -1;
    
    for (int i = 0; i < read_data->n_aligned; i++) {
        if (read_data->aligned_genomes[i] == genome_idx) {
            local_bit_idx = i;  // Use aligned array index as bit position
            break;
        }
    }
    
    if (local_bit_idx == -1) {
        fprintf(stderr, "ERROR: genome_idx=%d not found in aligned_genomes array\n", genome_idx);
        return;
    }
    
    if (local_bit_idx >= 0) {
        int byte_idx = local_bit_idx / 8;
        int bit_idx = local_bit_idx % 8;
        
        if (has_mismatch) {
            pos->has_mismatch[byte_idx] |= (1 << bit_idx);
            if (is_damage) {
                pos->is_damage[byte_idx] |= (1 << bit_idx);
            }
        }
    }
}

// Static reusable lookup array to eliminate repeated allocations
static int *reusable_aligned_lookup = NULL;
static int lookup_capacity = 0;
// Track previous usage for selective clearing
static int *prev_used_genomes = NULL;
static int prev_n_used = 0;

// Calculate taxonomic mismatches for a node - MEMORY OPTIMIZED VERSION
void calculate_taxonomic_mismatches(ReadMismatches *read_data, TaxNode *node, 
                                   damage_stats_t *stats) {
    memset(stats, 0, sizeof(damage_stats_t));
    
    if (!node || !node->n_leaves) return;
    
    // MEMORY OPTIMIZATION: Reuse lookup array instead of calloc/free every time
    // This eliminates 3.4TB of allocations for 1M alignments
    if (max_genomes_allocated > lookup_capacity) {
        reusable_aligned_lookup = realloc(reusable_aligned_lookup, 
                                         max_genomes_allocated * sizeof(int));
        if (!reusable_aligned_lookup) return; // Fallback on memory allocation failure
        
        // Also allocate tracking array
        prev_used_genomes = realloc(prev_used_genomes, 
                                   max_genomes_allocated * sizeof(int));
        if (!prev_used_genomes) return;
        
        lookup_capacity = max_genomes_allocated;
    }
    
    // SELECTIVE CLEARING OPTIMIZATION: Clear only entries used in previous call
    // This eliminates 58.1% CPU time spent in memset(2MB)
    for (int i = 0; i < prev_n_used; i++) {
        int genome_idx = prev_used_genomes[i];
        if (genome_idx < 0 || genome_idx >= max_genomes_allocated) {
            fprintf(stderr, "FATAL ERROR: Invalid prev genome_idx %d in reusable_aligned_lookup clear\n", genome_idx);
            fprintf(stderr, "  max_genomes_allocated=%d, clearing entry %d/%d\n", max_genomes_allocated, i, prev_n_used);
            abort();
        }
        reusable_aligned_lookup[genome_idx] = 0;
    }
    
    // Initialize lookup array: genome_idx -> (aligned_array_index + 1)
    // Using +1 so that 0 means "not aligned"  
    for (int i = 0; i < read_data->n_aligned; i++) {
        int genome_idx = read_data->aligned_genomes[i];
        if (genome_idx < 0) {
            fprintf(stderr, "FATAL ERROR: Negative genome index %d in reusable_aligned_lookup\n", genome_idx);
            fprintf(stderr, "  Read alignment %d/%d indicates memory corruption or uninitialized data.\n", i, read_data->n_aligned);
            abort();
        }
        if (genome_idx >= max_genomes_allocated) {
            fprintf(stderr, "FATAL ERROR: Genome index %d exceeds allocation %d\n", genome_idx, max_genomes_allocated);
            fprintf(stderr, "  Read alignment %d/%d indicates race condition or corrupted genome index.\n", i, read_data->n_aligned);
            abort();
        }
        reusable_aligned_lookup[genome_idx] = i + 1;
    }
    
    // Track which entries we're using for next selective clear
    prev_n_used = read_data->n_aligned;
    for (int i = 0; i < read_data->n_aligned; i++) {
        prev_used_genomes[i] = read_data->aligned_genomes[i];
    }
    
    // Calculate nd as the MAXIMUM nd value across all member genomes
    // that actually have alignments for this read
    int max_nd = 0;
    int any_genome_aligned = 0;
    
    // OPTIMIZED: Check each genome in this taxonomic group
    for (int i = 0; i < node->n_leaves; i++) {
        int genome_idx = node->leaf_genomes[i];
        
        // OPTIMIZATION: O(1) lookup instead of O(A) linear search
        if (genome_idx < 0 || genome_idx >= max_genomes_allocated) {
            fprintf(stderr, "FATAL ERROR: Invalid genome_idx %d in reusable_aligned_lookup read access\n", genome_idx);
            fprintf(stderr, "  max_genomes_allocated=%d, indicates corrupted taxonomy tree data.\n", max_genomes_allocated);
            abort();
        }
        int lookup_idx = reusable_aligned_lookup[genome_idx];
        if (lookup_idx > 0) {
            any_genome_aligned = 1;
            int aligned_idx = lookup_idx - 1; // Convert back to 0-based index
            
            // OPTIMIZATION: Direct access instead of second linear search
            int genome_nd = read_data->genome_nd_values ? read_data->genome_nd_values[aligned_idx] : 0;
            
            // Update maximum
            if (genome_nd > max_nd) {
                max_nd = genome_nd;
            }
        }
    }
    
    stats->nd = max_nd;
    
    // If no genomes were aligned, return early
    if (!any_genome_aligned) {
        return;
    }
    
    // For each position in the read
    for (int pos = 0; pos < read_data->trimmed_length; pos++) {
        if (!read_data->positions) break;
        
        PositionMismatch *pos_data = &read_data->positions[pos];
        
        // Check if any genome in this group has data at this position
        int has_any_data = 0;
        int mismatches_all = 1;  // Start assuming all have mismatches
        int mismatches_some = 0;
        int damages_all = 1;     // Start assuming all have damage
        int damages_some = 0;
        int genomes_with_alignments = 0;
        
        // OPTIMIZED: Check each genome in this taxonomic group
        for (int i = 0; i < node->n_leaves; i++) {
            int genome_idx = node->leaf_genomes[i];
            
            // OPTIMIZATION: O(1) lookup instead of O(A) linear search
            if (genome_idx < 0 || genome_idx >= max_genomes_allocated) {
                fprintf(stderr, "FATAL ERROR: Invalid genome_idx %d in reusable_aligned_lookup read access #2\n", genome_idx);
                fprintf(stderr, "  max_genomes_allocated=%d, indicates corrupted taxonomy tree data.\n", max_genomes_allocated);
                abort();
            }
            int lookup_idx = reusable_aligned_lookup[genome_idx];
            if (lookup_idx == 0) continue; // Genome not aligned, skip
            
            genomes_with_alignments++;
            has_any_data = 1;
            
            if (!pos_data->has_mismatch) {
                // No mismatch data at this position - all match
                mismatches_all = 0;
                damages_all = 0;
                continue;
            }
            
            // Use local bit index (aligned_idx) instead of global genome_idx
            int aligned_idx = lookup_idx - 1;  // Convert from lookup format
            int byte_idx = aligned_idx / 8;
            int bit_idx = aligned_idx % 8;
            
            int has_mismatch = (pos_data->has_mismatch[byte_idx] >> bit_idx) & 1;
            int is_damage = (pos_data->is_damage[byte_idx] >> bit_idx) & 1;
            
            if (has_mismatch) {
                mismatches_some = 1;
                if (is_damage) {
                    damages_some = 1;
                } else {
                    damages_all = 0;  // At least one genome has non-damage mismatch
                }
            } else {
                mismatches_all = 0;  // At least one genome doesn't have mismatch
                damages_all = 0;
            }
        }
        
        // Only count if we have aligned genomes in this group
        if (!has_any_data || genomes_with_alignments == 0) continue;
        
        // Update statistics based on patterns
        if (mismatches_all && genomes_with_alignments > 0) {
            // All aligned genomes have mismatches
            if (damages_all) {
                stats->md++;  // Damage mismatch to ALL
            } else {
                stats->mb++;  // Background mismatch to ALL (at least one is not damage)
            }
        } else if (mismatches_some) {
            // Some but not all have mismatches
            if (damages_some && !damages_all) {
                stats->mds++; // Damage mismatch to SOME
            } else {
                stats->mbs++; // Background mismatch to SOME
            }
        }
    }
    
    // No cleanup needed - reusable array persists for next call
}

// Finalize read mismatches (placeholder for any final processing)
int finalize_read_mismatches(ReadMismatches *data) {
    // Filter aligned genomes to keep only best alignment per taxid
    if (!data || data->n_aligned <= 0) return 0;



    // Create array to track best genome for each unique taxid
    int *taxid_to_best_idx = malloc(data->n_aligned * sizeof(int));
    int *unique_taxids = malloc(data->n_aligned * sizeof(int));

    // Remove all initialization to test original uninitialized behavior

    int n_unique_taxids = 0;

    if (!taxid_to_best_idx || !unique_taxids) {
        fprintf(stderr, "FATAL ERROR: Memory allocation failed in finalize_read_mismatches()\n");
        fprintf(stderr, "  Read: %s, alignments: %d\n", data->read_id, data->n_aligned);
        fprintf(stderr, "  Requested: %zu bytes for taxid arrays\n",
                (size_t)data->n_aligned * sizeof(int) * 2);
        free(taxid_to_best_idx);
        free(unique_taxids);
        abort();  // Make allocation failures fatal instead of silent
    }

    // Find best alignment for each taxid
    for (int i = 0; i < data->n_aligned; i++) {
        int genome_idx = data->aligned_genomes[i];

        // Validate genome index before use
        if (genome_idx < 0 || genome_idx >= n_genomes) {
            fprintf(stderr, "FATAL ERROR: Invalid genome index %d in finalize_read_mismatches()\n", genome_idx);
            fprintf(stderr, "  Read: %s, aligned entry %d/%d, n_genomes=%d\n",
                    data->read_id, i, data->n_aligned, n_genomes);
            fprintf(stderr, "  This indicates a bug in taxid/genome index management.\n");
            abort();
        }

        int taxid = get_taxid_from_index(genome_idx);
        int nm_value = data->genome_nm_values[i];


        // NOTE: The following duplicate detection logic is DEAD CODE
        // This was designed for genome-level duplicate detection, but in taxonomic mode
        // each alignment corresponds to a unique taxid, so duplicates never occur.

        /*
        // Find if this taxid is already tracked
        int taxid_idx = -1;
        for (int j = 0; j < n_unique_taxids; j++) {
            if (unique_taxids[j] == taxid) {
                taxid_idx = j;
                break;
            }
        }
        */

        // In taxonomic mode, all taxids are unique, so taxid_idx is always -1
        int taxid_idx = -1;

        if (taxid_idx == -1) {
            // New taxid - add it
            unique_taxids[n_unique_taxids] = taxid;
            taxid_to_best_idx[n_unique_taxids] = i;
            n_unique_taxids++;
        } else {
            // Existing taxid - check if this alignment is better
            int current_best_idx = taxid_to_best_idx[taxid_idx];
            int current_best_nm = data->genome_nm_values[current_best_idx];

            if (nm_value < current_best_nm) {
                // This alignment is better - update best
                taxid_to_best_idx[taxid_idx] = i;
            } else {
            }
        }
    }

    // Create filtered arrays with only best alignments per taxid
    int *new_aligned_genomes = malloc(n_unique_taxids * sizeof(int));
    int *new_genome_nd_values = malloc(n_unique_taxids * sizeof(int));
    int *new_genome_nm_values = malloc(n_unique_taxids * sizeof(int));

    if (new_aligned_genomes && new_genome_nd_values && new_genome_nm_values) {
        for (int i = 0; i < n_unique_taxids; i++) {
            int best_idx = taxid_to_best_idx[i];
            new_aligned_genomes[i] = data->aligned_genomes[best_idx];
            new_genome_nd_values[i] = data->genome_nd_values[best_idx];
            new_genome_nm_values[i] = data->genome_nm_values[best_idx];
        }

        // Replace arrays with filtered versions
        free(data->aligned_genomes);
        free(data->genome_nd_values);
        free(data->genome_nm_values);

        data->aligned_genomes = new_aligned_genomes;
        data->genome_nd_values = new_genome_nd_values;
        data->genome_nm_values = new_genome_nm_values;
        data->n_aligned = n_unique_taxids;
        data->capacity = n_unique_taxids;

    }

    free(taxid_to_best_idx);
    free(unique_taxids);
    return 1;  // Success
}

// Free read mismatches data
void free_read_mismatches(ReadMismatches *data) {
    if (!data) return;
    
    free(data->read_id);
    
    if (data->positions) {
        for (int i = 0; i < data->read_length; i++) {
            free(data->positions[i].has_mismatch);
            free(data->positions[i].is_damage);
        }
        free(data->positions);
    }
    
    free(data->aligned_genomes);
    free(data->genome_nd_values);
    free(data->genome_nm_values);
    free(data);
}

// Output results with taxonomy information
void output_with_taxonomy(FILE *fp, options_t *opts) {
    FILE *out_fp = fp;
    if (!out_fp) {
        out_fp = fopen(opts->output_file, "w");
        if (!out_fp) {
            fprintf(stderr, "Error: Cannot create output file %s\n", opts->output_file);
            return;
        }
    }
    
    // Determine if we're using sparse format
    int is_sparse = (strcmp(opts->output_format, "sparse") == 0 || 
                    strcmp(opts->output_format, "sparse_damage") == 0);
    int is_damage_format = opts->enable_damage && 
                          (strcmp(opts->output_format, "dense_damage") == 0 || 
                           strcmp(opts->output_format, "sparse_damage") == 0);
    
    // Collect active taxonomic nodes
    int n_active_taxa = 0;
    int *active_taxids = (int *)calloc(taxonomy_tree->max_taxid + 1, sizeof(int));
    
    for (int i = 1; i <= taxonomy_tree->max_taxid; i++) {
        if (taxonomy_tree->nodes[i] && taxonomy_tree->nodes[i]->is_active) {
            active_taxids[n_active_taxa++] = i;
        }
    }
    
    // Write header
    fprintf(out_fp, "read_id\ttotal_count");
    
    // Leaf genomes
    for (int i = 0; i < n_genomes; i++) {
        fprintf(out_fp, "\t%s", get_output_genome_name(i, opts->short_names));
    }
    
    // Separator for higher taxa
    fprintf(out_fp, "\t|");
    
    // Higher taxonomic groups (exclude leaf nodes that correspond to genomes)
    for (int i = 0; i < n_active_taxa; i++) {
        TaxNode *node = taxonomy_tree->nodes[active_taxids[i]];
        if (node && node->n_leaves > 0) {
            // Check if this taxid directly corresponds to any genome
            int is_leaf_genome = 0;
            for (int g = 0; g < n_genomes; g++) {
                if (get_taxid_from_index(g) == node->taxid) {
                    is_leaf_genome = 1;
                    break;
                }
            }
            // Only include if it's not a leaf genome
            if (!is_leaf_genome) {
                fprintf(out_fp, "\tT%d:%s", node->taxid, node->rank);
            }
        }
    }
    fprintf(out_fp, "\n");
    
    // Static reusable lookup array for output function (same pattern as calculate_taxonomic_mismatches)
    static int *output_reusable_lookup = NULL;
    static int output_lookup_capacity = 0;
    // Track previous usage for selective clearing
    static int *output_prev_used_genomes = NULL;
    static int output_prev_n_used = 0;
    
    // Process each read
    for (int hash_idx = 0; hash_idx < read_mismatch_table_size; hash_idx++) {
        ReadMismatchEntry *entry = read_mismatch_table[hash_idx];
        
        while (entry) {
            ReadMismatches *read_data = entry->data;
        
        // OPTIMIZATION: Reuse lookup array instead of malloc/free every time
        if (max_genomes_allocated > output_lookup_capacity) {
                output_reusable_lookup = realloc(output_reusable_lookup, 
                                               max_genomes_allocated * sizeof(int));
                if (!output_reusable_lookup) {
                    fprintf(stderr, "FATAL: Failed to allocate output_reusable_lookup for %d genomes\n", 
                            max_genomes_allocated);
                    abort();
                }
                
                // Also allocate tracking array
                output_prev_used_genomes = realloc(output_prev_used_genomes,
                                                  max_genomes_allocated * sizeof(int));
                if (!output_prev_used_genomes) {
                    fprintf(stderr, "FATAL: Failed to allocate output_prev_used_genomes for %d genomes\n", 
                            max_genomes_allocated);
                    abort();
                }
                
                output_lookup_capacity = max_genomes_allocated;
            }
            
            // SELECTIVE CLEARING: Clear only entries used in previous call
            for (int i = 0; i < output_prev_n_used; i++) {
                output_reusable_lookup[output_prev_used_genomes[i]] = 0;
            }
            
            // Initialize lookup: genome_idx -> (aligned_array_index + 1)
            for (int i = 0; i < read_data->n_aligned; i++) {
                output_reusable_lookup[read_data->aligned_genomes[i]] = i + 1;
            }
            
            // Track which entries we're using for next selective clear
            output_prev_n_used = read_data->n_aligned;
            for (int i = 0; i < read_data->n_aligned; i++) {
                output_prev_used_genomes[i] = read_data->aligned_genomes[i];
            }

            // Check --enforce-dense_strict filtering before output (taxonomy mode)
            int should_output = 1;
            if (opts->enforce_dense_strict) {
                // In taxonomy mode, we need to count unique taxids instead of genomes
                // We can get the total number of discovered taxids from the taxonomy system
                int total_taxids = taxonomy_tree->n_active_nodes;

                // Skip this read if it doesn't align to ALL taxids
                if (read_data->n_aligned != total_taxids) {
                    should_output = 0;
                }
            }

            if (should_output) {
                fprintf(out_fp, "%s\t%d", read_data->read_id, read_data->trimmed_length);
            
            if (is_sparse) {
                // SPARSE FORMAT: Only output genomes with alignments
                for (int genome_idx = 0; genome_idx < n_genomes; genome_idx++) {
                    // OPTIMIZATION: O(1) lookup instead of O(A) linear search
                    int lookup_idx = output_reusable_lookup[genome_idx];
                    if (lookup_idx > 0) {
                        int aligned_idx = lookup_idx - 1;
                        
                        // Output genome name first for sparse format
                        fprintf(out_fp, "\t%s", get_output_genome_name(genome_idx, opts->short_names));
                        
                        damage_stats_t stats = {0};
                        
                        // OPTIMIZATION: Direct access instead of linear search
                        int genome_nd = read_data->genome_nd_values ? read_data->genome_nd_values[aligned_idx] : 0;
                        stats.nd = genome_nd;
                        
                        // Calculate mismatches for this genome
                        for (int pos = 0; pos < read_data->trimmed_length; pos++) {
                            if (read_data->positions && read_data->positions[pos].has_mismatch) {
                                // Use local bit index (aligned_idx) instead of global genome_idx
                                int byte_idx = aligned_idx / 8;
                                int bit_idx = aligned_idx % 8;
                                
                                if (read_data->positions[pos].has_mismatch[byte_idx] & (1 << bit_idx)) {
                                    if (read_data->positions[pos].is_damage[byte_idx] & (1 << bit_idx)) {
                                        stats.md++;
                                    } else {
                                        stats.mb++;
                                    }
                                }
                            }
                        }
                        
                        if (is_damage_format) {
                            // For individual genomes: nd md mb (3 numbers) - matching BAMreader.c order
                            fprintf(out_fp, "\t%d\t%d\t%d", stats.nd, stats.md, stats.mb);
                        } else {
                            // Without damage: mb+md (1 number for individual genomes)
                            fprintf(out_fp, "\t%d", stats.mb + stats.md);
                        }
                    }
                }
            } else {
                // DENSE FORMAT: Output all genomes
                for (int genome_idx = 0; genome_idx < n_genomes; genome_idx++) {
                    // OPTIMIZATION: O(1) lookup instead of O(A) linear search
                    int lookup_idx = output_reusable_lookup[genome_idx];
                    int is_aligned = (lookup_idx > 0);
                    
                    if (!is_aligned) {
                        // No alignment - output -1
                        if (is_damage_format) {
                            fprintf(out_fp, "\t-1\t-1\t-1");
                        } else {
                            fprintf(out_fp, "\t-1");
                        }
                    } else {
                        damage_stats_t stats = {0};
                        
                        // OPTIMIZATION: Direct access using lookup array
                        int aligned_idx = lookup_idx - 1;
                        int genome_nd = read_data->genome_nd_values ? read_data->genome_nd_values[aligned_idx] : 0;
                        stats.nd = genome_nd;
                        
                        // Calculate mismatches for this genome
                        for (int pos = 0; pos < read_data->trimmed_length; pos++) {
                            if (read_data->positions && read_data->positions[pos].has_mismatch) {
                                // Use local bit index (aligned_idx) instead of global genome_idx
                                int byte_idx = aligned_idx / 8;
                                int bit_idx = aligned_idx % 8;
                                
                                if (read_data->positions[pos].has_mismatch[byte_idx] & (1 << bit_idx)) {
                                    if (read_data->positions[pos].is_damage[byte_idx] & (1 << bit_idx)) {
                                        stats.md++;
                                    } else {
                                        stats.mb++;
                                    }
                                }
                            }
                        }
                        
                        if (is_damage_format) {
                            // For individual genomes: nd md mb (3 numbers) - matching BAMreader.c order
                            fprintf(out_fp, "\t%d\t%d\t%d", stats.nd, stats.md, stats.mb);
                        } else {
                            // Without damage: mb+md (1 number for individual genomes)
                            fprintf(out_fp, "\t%d", stats.mb + stats.md);
                        }
                    }
                }
            }

            // Separator
            fprintf(out_fp, "\t|");

            // Output higher taxonomic group mismatches (exclude leaf nodes that correspond to genomes)
            for (int i = 0; i < n_active_taxa; i++) {
                TaxNode *node = taxonomy_tree->nodes[active_taxids[i]];
                if (node && node->n_leaves > 0) {
                    // Check if this taxid directly corresponds to any genome
                    int is_leaf_genome = 0;
                    for (int g = 0; g < n_genomes; g++) {
                        if (get_taxid_from_index(g) == node->taxid) {
                            is_leaf_genome = 1;
                            break;
                        }
                    }
                    // Skip if it's a leaf genome
                    if (is_leaf_genome) {
                        continue;
                    }
                    damage_stats_t tax_stats = {0};
                    calculate_taxonomic_mismatches(read_data, node, &tax_stats);

                    // In sparse format, only output if this taxonomic group has aligned genomes
                    if (is_sparse && tax_stats.nd == 0) {
                        continue;  // Skip taxonomic groups with no aligned genomes
                    }

                    // For sparse format with taxonomy, output taxid first
                    if (is_sparse) {
                        fprintf(out_fp, "\tT%d:%s", node->taxid, node->rank);
                    }

                    if (is_damage_format) {
                        // For taxonomic groups: nd md mds mb mbs (5 numbers)
                        fprintf(out_fp, "\t%d\t%d\t%d\t%d\t%d",
                                tax_stats.nd, tax_stats.md, tax_stats.mds,
                                tax_stats.mb, tax_stats.mbs);
                    } else {
                        fprintf(out_fp, "\t%d\t%d",
                                tax_stats.mb + tax_stats.md,
                                tax_stats.mbs + tax_stats.mds);
                    }
                }
            }

            fprintf(out_fp, "\n");
        }  // End if (should_output)

        entry = entry->next;
        }
    }
    
    free(active_taxids);
    
    if (!fp) {
        fclose(out_fp);
    }
}

// ============ LINEAR ALGORITHM IMPLEMENTATION ============

// Initialize node pool for per-read taxonomy trees
void init_node_pool(void) {
    if (!node_pool) {
        node_pool_capacity = 20000;  // Start with 20K nodes
        node_pool = (TreeNode *)calloc(node_pool_capacity, sizeof(TreeNode));
        
        // Pre-allocate shared bit vector memory (4 vectors + 2 parsimony vectors per node)
        int bytes_per_vector = (max_read_length + 7) / 8;
        int parsimony_bytes_per_vector = (max_read_length * 2 + 7) / 8;  // 2 bits per position
        int total_bytes_per_node = 4 * bytes_per_vector + 2 * parsimony_bytes_per_vector;  // Now 2 parsimony vectors
        shared_bit_vectors = (uint8_t *)calloc(node_pool_capacity * total_bytes_per_node, 1);
        
        if (!node_pool || !shared_bit_vectors) {
            fprintf(stderr, "Error: Failed to allocate node pool\n");
            return;
        }
        
        // Initialize bit vector pointers for each node
        for (int i = 0; i < node_pool_capacity; i++) {
            int offset = i * total_bytes_per_node;
            node_pool[i].md_vec = shared_bit_vectors + offset;
            node_pool[i].mds_vec = shared_bit_vectors + offset + bytes_per_vector;
            node_pool[i].mb_vec = shared_bit_vectors + offset + 2 * bytes_per_vector;
            node_pool[i].mbs_vec = shared_bit_vectors + offset + 3 * bytes_per_vector;
            node_pool[i].parsimony_state_vec = shared_bit_vectors + offset + 4 * bytes_per_vector;
            node_pool[i].damage_parsimony_state_vec = shared_bit_vectors + offset + 4 * bytes_per_vector + parsimony_bytes_per_vector;
        }
    }
    
    // Initialize reusable node_map (one-time 80MB allocation)
    if (!reusable_node_map) {
        reusable_node_map = (TreeNode **)calloc(MAX_TAXID + 1, sizeof(TreeNode*));
        if (!reusable_node_map) {
            fprintf(stderr, "Error: Failed to allocate reusable node_map\n");
            return;
        }
        node_map_allocated = 1;

        // Initialize tracking array for used taxids
        used_taxids_capacity = 10000;  // Start with capacity for 10K taxids
        used_taxids = (int *)malloc(used_taxids_capacity * sizeof(int));
        if (!used_taxids) {
            fprintf(stderr, "Error: Failed to allocate used_taxids tracking array\n");
            return;
        }
        n_used_taxids = 0;

        fprintf(stderr, "Initializing linear algorithm for taxonomic analysis...\n");
    }
}

// Get a node from the pool (with bounds checking)
TreeNode* get_node_from_pool(int taxid, const char *rank) {
    if (node_pool_used >= node_pool_capacity) {
        // Grow pool by 5000 nodes
        int old_capacity = node_pool_capacity;
        node_pool_capacity += 5000;
        
        node_pool = (TreeNode *)realloc(node_pool, node_pool_capacity * sizeof(TreeNode));
        
        // Reallocate bit vector memory
        int bytes_per_vector = (max_read_length + 7) / 8;
        int parsimony_bytes_per_vector = (max_read_length * 2 + 7) / 8;  // 2 bits per position
        int total_bytes_per_node = 4 * bytes_per_vector + 2 * parsimony_bytes_per_vector;  // Now 2 parsimony vectors
        shared_bit_vectors = (uint8_t *)realloc(shared_bit_vectors,
                                               node_pool_capacity * total_bytes_per_node);
        
        if (!node_pool || !shared_bit_vectors) {
            fprintf(stderr, "Error: Failed to grow node pool to %d nodes\n", node_pool_capacity);
            return NULL;
        }
        
        // Initialize bit vector pointers for new nodes
        for (int i = old_capacity; i < node_pool_capacity; i++) {
            int offset = i * total_bytes_per_node;
            node_pool[i].md_vec = shared_bit_vectors + offset;
            node_pool[i].mds_vec = shared_bit_vectors + offset + bytes_per_vector;
            node_pool[i].mb_vec = shared_bit_vectors + offset + 2 * bytes_per_vector;
            node_pool[i].mbs_vec = shared_bit_vectors + offset + 3 * bytes_per_vector;
            node_pool[i].parsimony_state_vec = shared_bit_vectors + offset + 4 * bytes_per_vector;
            node_pool[i].damage_parsimony_state_vec = shared_bit_vectors + offset + 4 * bytes_per_vector + parsimony_bytes_per_vector;
        }
    }
    
    TreeNode *node = &node_pool[node_pool_used++];
    
    // Initialize node
    node->taxid = taxid;
    strncpy(node->rank, rank, MAX_RANK_LEN - 1);
    node->rank[MAX_RANK_LEN - 1] = '\0';
    node->parent = NULL;
    node->children = NULL;
    node->n_children = 0;
    node->max_children = 0;
    node->nd_value = 0;
    node->genome_idx = -1;  // Default to internal node
    
    // Clear bit vectors (only clear what we'll use)
    int bytes_per_vector = (max_read_length + 7) / 8;
    int parsimony_bytes_per_vector = (max_read_length * 2 + 7) / 8;  // 2 bits per position
    memset(node->md_vec, 0, bytes_per_vector);
    memset(node->mds_vec, 0, bytes_per_vector);
    memset(node->mb_vec, 0, bytes_per_vector);
    memset(node->mbs_vec, 0, bytes_per_vector);
    memset(node->parsimony_state_vec, 0, parsimony_bytes_per_vector);
    
    return node;
}

// Reset node pool for next read (reuse memory)
void reset_node_pool(void) {
    // Free any dynamically allocated children arrays before resetting
    for (int i = 0; i < node_pool_used; i++) {
        if (node_pool[i].children) {
            free(node_pool[i].children);
            node_pool[i].children = NULL;
            node_pool[i].max_children = 0;
        }
    }
    
    // Reset usage counter - memory stays allocated for reuse
    node_pool_used = 0;
}


// Comprehensive timing for profiling analysis
static struct timeval program_start_time;
static double genome_discovery_time = 0;
static double bam_processing_time = 0;
static double measured_bam_processing_time = 0;  // Backup timing variable
static double linear_algorithm_time = 0;

// BAM processing component timing
static double total_read_lookup_time = 0;
static double total_mismatch_calc_time = 0;
static double total_read_storage_time = 0;
static long timing_sample_count = 0;

// Build minimal taxonomy tree for this read's aligned genomes
TreeNode* build_read_taxonomy_tree(ReadMismatches *read_data) {
    if (!read_data || read_data->n_aligned == 0) return NULL;
    
    // Use reusable node_map to eliminate 80MB allocation per read
    if (!reusable_node_map) {
        fprintf(stderr, "Error: Reusable node_map not initialized\n");
        return NULL;
    }

    // Clear only the previously used entries instead of entire array
    for (int i = 0; i < n_used_taxids; i++) {
        reusable_node_map[used_taxids[i]] = NULL;
    }
    n_used_taxids = 0;  // Reset the counter for this read

    TreeNode **node_map = reusable_node_map;  // Use the reusable allocation
    
    TreeNode *root = NULL;
    
    // For each aligned genome, build path from leaf to root
    int valid_genomes = 0;
    for (int i = 0; i < read_data->n_aligned; i++) {
        int genome_idx = read_data->aligned_genomes[i];
        int taxid = get_taxid_from_index(genome_idx);
        
        if (taxid <= 0 || taxid > MAX_TAXID) continue;  // Skip invalid taxids
        valid_genomes++;
        
        // Build path from leaf to root
        int current_taxid = taxid;
        TreeNode *child_node = NULL;
        
        while (current_taxid > 0 && current_taxid <= MAX_TAXID) {
            TaxNode *global_node = taxonomy_tree->nodes[current_taxid];
            if (!global_node) break;
            
            // Check if we've already created this node
            if (!node_map[current_taxid]) {
                TreeNode *new_node = get_node_from_pool(current_taxid, global_node->rank);
                if (!new_node) return NULL;  // Pool exhausted

                node_map[current_taxid] = new_node;

                // Track this taxid as used
                if (n_used_taxids >= used_taxids_capacity) {
                    // Grow the tracking array
                    used_taxids_capacity *= 2;
                    used_taxids = (int *)realloc(used_taxids, used_taxids_capacity * sizeof(int));
                    if (!used_taxids) {
                        fprintf(stderr, "Error: Failed to realloc used_taxids array\n");
                        return NULL;
                    }
                }
                used_taxids[n_used_taxids++] = current_taxid;

                // Set as leaf if this is the genome's direct taxid
                if (current_taxid == taxid) {
                    new_node->genome_idx = genome_idx;
                }
            } else {
                // Node already exists - check if this is a leaf taxid collision
                if (current_taxid == taxid) {
                    TreeNode *existing_node = node_map[current_taxid];
                    
                    // Get NM values for both genomes to choose the better alignment
                    int current_genome_nm = -1;
                    int existing_genome_nm = -1;
                    
                    // Find NM values in read_data
                    for (int j = 0; j < read_data->n_aligned; j++) {
                        if (read_data->aligned_genomes[j] == genome_idx) {
                            current_genome_nm = read_data->genome_nm_values[j];
                        }
                        if (read_data->aligned_genomes[j] == existing_node->genome_idx) {
                            existing_genome_nm = read_data->genome_nm_values[j];
                        }
                    }
                    
                    // Keep the genome with better (lower) NM value
                    if (current_genome_nm >= 0 && existing_genome_nm >= 0) {
                        if (current_genome_nm < existing_genome_nm) {
                            existing_node->genome_idx = genome_idx;
                        }
                    } else {
                        // Fallback: just use the new one (original behavior)
                        existing_node->genome_idx = genome_idx;
                    }
                }
            }
            
            TreeNode *current_node = node_map[current_taxid];
            
            // Link child to parent
            if (child_node && current_node != child_node->parent) {
                // Add child to parent's children array
                if (current_node->n_children >= current_node->max_children) {
                    current_node->max_children = current_node->max_children ? 
                                                current_node->max_children * 2 : 4;
                    current_node->children = (TreeNode **)realloc(current_node->children,
                                                                 current_node->max_children * sizeof(TreeNode *));
                }
                current_node->children[current_node->n_children++] = child_node;
                child_node->parent = current_node;
            }
            
            child_node = current_node;
            
            // Move to parent
            int next_taxid = global_node->parent_taxid;
            
            // Track root (node with no parent in taxonomy)
            if (next_taxid == current_taxid || next_taxid <= 0) {
                root = current_node;
                break;
            }
            
            current_taxid = next_taxid;
        }
    }

    // No need to free - using reusable allocation
    return root;
}

// Post-order traversal with position-wise vector combination
void calculate_node_vectors_postorder(TreeNode *node, ReadMismatches *read_data, int read_length) {
    if (!node) return;
    
    if (node->n_children == 0) {
        // LEAF NODE: Initialize from alignment bit vectors
        if (node->genome_idx >= 0) {
            // Find this genome's index in aligned_genomes array
            int aligned_idx = -1;
            for (int i = 0; i < read_data->n_aligned; i++) {
                if (read_data->aligned_genomes[i] == node->genome_idx) {
                    aligned_idx = i;
                    break;
                }
            }
            
            if (aligned_idx >= 0) {
                // Copy from position data using local bit indexing
                for (int pos = 0; pos < read_length; pos++) {
                    if (read_data->positions && read_data->positions[pos].has_mismatch) {
                        int byte_idx = aligned_idx / 8;
                        int bit_idx = aligned_idx % 8;
                        
                        // Ensure position bit vector is large enough for current access
                        PositionMismatch *pos_data = &read_data->positions[pos];
                        int required_bytes = (aligned_idx / 8) + 1;
                        int allocated_bytes = (pos_data->n_genomes + 7) / 8;
                        
                        if (required_bytes > allocated_bytes) {
                            // Need to resize this position's bit vectors
                            pos_data->has_mismatch = (uint8_t *)realloc(pos_data->has_mismatch, required_bytes);
                            pos_data->is_damage = (uint8_t *)realloc(pos_data->is_damage, required_bytes);
                            
                            // Clear new bytes
                            memset(pos_data->has_mismatch + allocated_bytes, 0, required_bytes - allocated_bytes);
                            memset(pos_data->is_damage + allocated_bytes, 0, required_bytes - allocated_bytes);
                            
                            // Update capacity to match actual requirement
                            pos_data->n_genomes = required_bytes * 8;
                        }
                        
                        int pos_byte = pos / 8;
                        int pos_bit = pos % 8;
                        
                        // Check if this genome has mismatch at this position
                        int has_mismatch_bit = (read_data->positions[pos].has_mismatch[byte_idx] & (1 << bit_idx)) ? 1 : 0;


                        if (has_mismatch_bit) {
                            if (read_data->positions[pos].is_damage[byte_idx] & (1 << bit_idx)) {
                                // Damage mismatch - set md for this position
                                node->md_vec[pos_byte] |= (1 << pos_bit);
                            } else {
                                // Background mismatch - set mb for this position
                                node->mb_vec[pos_byte] |= (1 << pos_bit);

                            }
                        }
                    }
                }
                
                // Set nd value from genome data
                node->nd_value = read_data->genome_nd_values ? read_data->genome_nd_values[aligned_idx] : 0;
            }
        }
        return;
    }
    
    // INTERNAL NODE: Process all children first (post-order)
    for (int i = 0; i < node->n_children; i++) {
        calculate_node_vectors_postorder(node->children[i], read_data, read_length);
    }
    
    // Combine children vectors position by position
    for (int pos = 0; pos < read_length; pos++) {
        int pos_byte = pos / 8;
        int pos_bit = pos % 8;
        
        // Count children with damage/background at this position
        int all_children_md = 1, some_children_md = 0;
        int all_children_mb = 1, some_children_mb = 0;
        int max_nd = 0;
        
        for (int i = 0; i < node->n_children; i++) {
            TreeNode *child = node->children[i];
            
            // Check damage patterns
            int child_has_md = (child->md_vec[pos_byte] & (1 << pos_bit)) ? 1 : 0;
            int child_has_mds = (child->mds_vec[pos_byte] & (1 << pos_bit)) ? 1 : 0;
            
            if (!child_has_md) all_children_md = 0;
            if (child_has_md || child_has_mds) some_children_md = 1;
            
            // Check background patterns
            int child_has_mb = (child->mb_vec[pos_byte] & (1 << pos_bit)) ? 1 : 0;
            int child_has_mbs = (child->mbs_vec[pos_byte] & (1 << pos_bit)) ? 1 : 0;
            
            if (!child_has_mb) all_children_mb = 0;
            if (child_has_mb || child_has_mbs) some_children_mb = 1;
            
            // Track max nd
            if (child->nd_value > max_nd) max_nd = child->nd_value;
        }
        
        // Apply damage precedence logic
        if (all_children_md) {
            // ALL children have damage mismatch
            node->md_vec[pos_byte] |= (1 << pos_bit);
        } else if (some_children_md) {
            // SOME children have damage mismatch
            node->mds_vec[pos_byte] |= (1 << pos_bit);
        } else if (all_children_mb) {
            // No damage, but ALL children have background mismatch
            node->mb_vec[pos_byte] |= (1 << pos_bit);
        } else if (some_children_mb) {
            // SOME children have background mismatch
            node->mbs_vec[pos_byte] |= (1 << pos_bit);
        }
        
        // Set nd as max of children
        node->nd_value = max_nd;
    }
}

// Parsimony post-order traversal - implements the parsimony algorithm
void calculate_parsimony_postorder(TreeNode *node, ReadMismatches *read_data, int read_length) {
    if (!node) return;
    
    if (node->n_children == 0) {
        // LEAF NODE: Set state based on mismatch data
        if (node->genome_idx >= 0) {
            // Find this genome's index in aligned_genomes array
            int aligned_idx = -1;
            for (int i = 0; i < read_data->n_aligned; i++) {
                if (read_data->aligned_genomes[i] == node->genome_idx) {
                    aligned_idx = i;
                    break;
                }
            }
            
            if (aligned_idx >= 0) {
                // Clear both parsimony state vectors
                int parsimony_vec_size = (read_length * 2 + 7) / 8;
                memset(node->parsimony_state_vec, 0, parsimony_vec_size);
                memset(node->damage_parsimony_state_vec, 0, parsimony_vec_size);

                // Set parsimony state for each position
                for (int pos = 0; pos < read_length; pos++) {
                    int mb_state = PARSIMONY_STATE_0;  // Default: no background mismatch
                    int md_state = PARSIMONY_STATE_0;  // Default: no damage mismatch

                    // Check if this position has mismatch data
                    if (read_data->positions && read_data->positions[pos].has_mismatch) {
                        int byte_idx = aligned_idx / 8;
                        int bit_idx = aligned_idx % 8;

                        int allocated_bytes = (read_data->positions[pos].n_genomes + 7) / 8;

                        // Check if this genome has mismatch at this position
                        if (read_data->positions[pos].has_mismatch[byte_idx] & (1 << bit_idx)) {
                            // Check if it's damage or background
                            if (read_data->positions[pos].is_damage[byte_idx] & (1 << bit_idx)) {
                                md_state = PARSIMONY_STATE_1;  // Has damage mismatch
                            } else {
                                mb_state = PARSIMONY_STATE_1;  // Has background mismatch
                            }
                        }
                    }

                    set_parsimony_state(node->parsimony_state_vec, pos, mb_state);
                    set_parsimony_state(node->damage_parsimony_state_vec, pos, md_state);
                }
                
                // Set nd value from genome data (same as old method)
                node->nd_value = read_data->genome_nd_values ? read_data->genome_nd_values[aligned_idx] : 0;
            }
        }
        return;
    }
    
    // INTERNAL NODE: Process children first (post-order)
    for (int i = 0; i < node->n_children; i++) {
        calculate_parsimony_postorder(node->children[i], read_data, read_length);
    }
    
    // Apply parsimony rules for background mismatches
    for (int pos = 0; pos < read_length; pos++) {
        int has_pure_0 = 0, has_pure_1 = 0;

        // Check children states
        for (int i = 0; i < node->n_children; i++) {
            int child_state = get_parsimony_state(node->children[i]->parsimony_state_vec, pos);
            if (child_state == PARSIMONY_STATE_0) has_pure_0 = 1;
            if (child_state == PARSIMONY_STATE_1) has_pure_1 = 1;
        }

        // Apply parsimony rules
        int node_state;
        if (has_pure_0 && !has_pure_1) {
            node_state = PARSIMONY_STATE_0;  // At least one 0, no 1s
        } else if (has_pure_1 && !has_pure_0) {
            node_state = PARSIMONY_STATE_1;  // At least one 1, no 0s
        } else {
            node_state = PARSIMONY_AMBIGUOUS;  // Mixed or all ambiguous
        }

        set_parsimony_state(node->parsimony_state_vec, pos, node_state);
    }

    // Apply parsimony rules for damage mismatches (same algorithm)
    for (int pos = 0; pos < read_length; pos++) {
        int has_pure_0 = 0, has_pure_1 = 0;

        // Check children states
        for (int i = 0; i < node->n_children; i++) {
            int child_state = get_parsimony_state(node->children[i]->damage_parsimony_state_vec, pos);
            if (child_state == PARSIMONY_STATE_0) has_pure_0 = 1;
            if (child_state == PARSIMONY_STATE_1) has_pure_1 = 1;
        }

        // Apply parsimony rules
        int node_state;
        if (has_pure_0 && !has_pure_1) {
            node_state = PARSIMONY_STATE_0;  // At least one 0, no 1s
        } else if (has_pure_1 && !has_pure_0) {
            node_state = PARSIMONY_STATE_1;  // At least one 1, no 0s
        } else {
            node_state = PARSIMONY_AMBIGUOUS;  // Mixed or all ambiguous
        }

        set_parsimony_state(node->damage_parsimony_state_vec, pos, node_state);
    }
    
    // Set nd as max of children
    int max_nd = 0;
    for (int i = 0; i < node->n_children; i++) {
        if (node->children[i]->nd_value > max_nd) {
            max_nd = node->children[i]->nd_value;
        }
    }
    node->nd_value = max_nd;
}

// Parsimony pre-order resolution - resolves ambiguous states using parent information
void resolve_parsimony_preorder(TreeNode *node, int read_length) {
    if (!node) return;
    
    // Resolve ambiguous states for background mismatches
    for (int pos = 0; pos < read_length; pos++) {
        int current_state = get_parsimony_state(node->parsimony_state_vec, pos);

        if (current_state == PARSIMONY_AMBIGUOUS) {
            // Resolve based on parent state
            if (node->parent) {
                int parent_state = get_parsimony_state(node->parent->parsimony_state_vec, pos);
                if (parent_state == PARSIMONY_STATE_0 || parent_state == PARSIMONY_STATE_1) {
                    set_parsimony_state(node->parsimony_state_vec, pos, parent_state);
                } else {
                    // Parent is also ambiguous, default to state 0
                    set_parsimony_state(node->parsimony_state_vec, pos, PARSIMONY_STATE_0);
                }
            } else {
                // This is the root - set ambiguous to state 0
                set_parsimony_state(node->parsimony_state_vec, pos, PARSIMONY_STATE_0);
            }
        }
    }

    // Resolve ambiguous states for damage mismatches (same algorithm)
    for (int pos = 0; pos < read_length; pos++) {
        int current_state = get_parsimony_state(node->damage_parsimony_state_vec, pos);

        if (current_state == PARSIMONY_AMBIGUOUS) {
            // Resolve based on parent state
            if (node->parent) {
                int parent_state = get_parsimony_state(node->parent->damage_parsimony_state_vec, pos);
                if (parent_state == PARSIMONY_STATE_0 || parent_state == PARSIMONY_STATE_1) {
                    set_parsimony_state(node->damage_parsimony_state_vec, pos, parent_state);
                } else {
                    // Parent is also ambiguous, default to state 0
                    set_parsimony_state(node->damage_parsimony_state_vec, pos, PARSIMONY_STATE_0);
                }
            } else {
                // This is the root - set ambiguous to state 0
                set_parsimony_state(node->damage_parsimony_state_vec, pos, PARSIMONY_STATE_0);
            }
        }
    }
    
    // Process children in pre-order
    for (int i = 0; i < node->n_children; i++) {
        resolve_parsimony_preorder(node->children[i], read_length);
    }
}

// Debug function to print parsimony states for validation
void debug_print_parsimony_states(TreeNode *node, ReadMismatches *read_data, const char *read_id, int read_length) {
    if (!node) return;
    
    const char* state_names[] = {"unset", "0", "1", "0|1"};
    
    if (node->genome_idx >= 0) {
        // LEAF NODE - print genome and its states
        
    } else {
        // INTERNAL NODE - print taxonomic group and its states
    }
    
    // Recursively print children
    for (int i = 0; i < node->n_children; i++) {
        debug_print_parsimony_states(node->children[i], read_data, read_id, read_length);
    }
}

// Helper function to sum bits in a vector (count set positions)
int sum_bit_vector(uint8_t *vector, int read_length) {
    int count = 0;
    for (int pos = 0; pos < read_length; pos++) {
        int byte_idx = pos / 8;
        int bit_idx = pos % 8;
        if (vector[byte_idx] & (1 << bit_idx)) {
            count++;
        }
    }
    return count;
}

// Parsimony 2-bit state manipulation functions
// Get parsimony state at position (returns 0-3)
int get_parsimony_state(uint8_t *state_vec, int pos) {
    int byte_idx = pos / 4;  // 4 positions per byte (2 bits each)
    int bit_offset = (pos % 4) * 2;  // Position within byte
    return (state_vec[byte_idx] >> bit_offset) & 3;
}

// Set parsimony state at position (state should be 0-3)
void set_parsimony_state(uint8_t *state_vec, int pos, int state) {
    int byte_idx = pos / 4;  // 4 positions per byte (2 bits each)
    int bit_offset = (pos % 4) * 2;  // Position within byte
    uint8_t mask = ~(3 << bit_offset);  // Clear the 2 bits
    state_vec[byte_idx] = (state_vec[byte_idx] & mask) | ((state & 3) << bit_offset);
}

// Count positions with state 1 for mb calculation
int count_parsimony_state_1(uint8_t *state_vec, int read_length) {
    int count = 0;
    for (int pos = 0; pos < read_length; pos++) {
        if (get_parsimony_state(state_vec, pos) == PARSIMONY_STATE_1) {
            count++;
        }
    }
    return count;
}

// Helper function for tree traversal to collect all nodes with bounds checking
void collect_tree_nodes(TreeNode *node, TreeNode **nodes, int *count, int max_nodes) {
    if (!node) return;
    
    // Bounds checking to prevent buffer overflow
    if (*count >= max_nodes) {
        fprintf(stderr, "ERROR: Tree has more than %d nodes - increase tree node array size\n", max_nodes);
        return;  // Prevent overflow
    }
    
    // Add this node to collection
    nodes[(*count)++] = node;
    
    // Recursively collect children
    for (int i = 0; i < node->n_children; i++) {
        collect_tree_nodes(node->children[i], nodes, count, max_nodes);
    }
}

// Linear complexity output function - O(reads × tree_height) instead of O(reads × genomes × taxonomic_groups)
void output_with_taxonomy_linear(FILE *fp, options_t *opts) {

    FILE *out_fp = fp;
    if (!out_fp) {
        out_fp = fopen(opts->output_file, "w");
        if (!out_fp) {
            fprintf(stderr, "Error: Cannot create output file %s\n", opts->output_file);
            return;
        }
    }
    
    // Initialize node pool
    init_node_pool();
    
    // Determine output format
    int is_sparse = (strcmp(opts->output_format, "sparse") == 0 || 
                    strcmp(opts->output_format, "sparse_damage") == 0);
    int is_damage_format = opts->enable_damage && 
                          (strcmp(opts->output_format, "dense_damage") == 0 || 
                           strcmp(opts->output_format, "sparse_damage") == 0);
    
    // Write header (same as original)
    fprintf(out_fp, "read_id\ttotal_count");
    
    // Leaf genomes
    for (int i = 0; i < n_genomes; i++) {
        fprintf(out_fp, "\t%s", get_output_genome_name(i, opts->short_names));
    }
    
    // Separator for higher taxa
    fprintf(out_fp, "\t|");
    
    // Header for taxonomic groups (simplified for now - could be optimized)
    // Use existing taxonomy tree to write header
    for (int i = 1; i <= taxonomy_tree->max_taxid; i++) {
        TaxNode *node = taxonomy_tree->nodes[i];
        if (node && node->is_active && node->n_leaves > 0) {
            // Check if it's not a leaf genome
            int is_leaf_genome = 0;
            for (int g = 0; g < n_genomes; g++) {
                if (get_taxid_from_index(g) == node->taxid) {
                    is_leaf_genome = 1;
                    break;
                }
            }
            if (!is_leaf_genome) {
                fprintf(out_fp, "\tT%d:%s", node->taxid, node->rank);
            }
        }
    }
    fprintf(out_fp, "\n");
    
    // Process each read with linear algorithm
    int reads_processed = 0;
    
    for (int hash_idx = 0; hash_idx < read_mismatch_table_size; hash_idx++) {
        ReadMismatchEntry *entry = read_mismatch_table[hash_idx];
        
        while (entry) {
            ReadMismatches *read_data = entry->data;

            reads_processed++;


            if (!opts->silent && reads_processed % 10000 == 0) {
                fprintf(stderr, "Progress: Processed taxonomy for %d reads\n", reads_processed);
            }
            
            // Reset node pool for this read
            reset_node_pool();
            
            
            // Step 1: Build minimal taxonomy tree for this read's aligned genomes
            TreeNode *tree_root = build_read_taxonomy_tree(read_data);
            
            if (tree_root) {
                // Step 2: Choose algorithm based on options
                if (opts->higher_taxa_with_mbs) {
                    // Use old ALL/SOME method (5 values)
                    calculate_node_vectors_postorder(tree_root, read_data, read_data->trimmed_length);
                } else {
                    // Use new parsimony method (3 values - same as genomic groups)
                    calculate_parsimony_postorder(tree_root, read_data, read_data->trimmed_length);
                    resolve_parsimony_preorder(tree_root, read_data->trimmed_length);
                    
                }
                
                // Collect all tree nodes for output with dynamic allocation
                TreeNode **all_nodes = malloc(20000 * sizeof(TreeNode*));
                if (!all_nodes) {
                    fprintf(stderr, "Error: Failed to allocate tree nodes array\n");
                    continue;  // Skip this read
                }
                int node_count = 0;
                collect_tree_nodes(tree_root, all_nodes, &node_count, 20000);
                
                
                
                // Step 3: Output results

                // Check --enforce-dense_strict filtering before output (taxonomy mode)
                int should_output = 1;
                if (opts->enforce_dense_strict) {
                    // In taxonomy mode, we need to count unique taxids instead of genomes
                    // We can get the total number of discovered taxids from the taxonomy system
                    int total_taxids = taxonomy_tree->n_active_nodes;

                    // Skip this read if it doesn't align to ALL taxids
                    if (read_data->n_aligned != total_taxids) {
                        should_output = 0;
                    }
                }

                if (should_output) {
                    fprintf(out_fp, "%s\t%d", read_data->read_id, read_data->trimmed_length);
                
                // Output leaf genomes (traverse tree leaves for aligned genomes only)
                if (is_sparse) {
                    // Output leaf nodes (genomes with alignments)
                    for (int i = 0; i < node_count; i++) {
                        TreeNode *node = all_nodes[i];
                        if (node->genome_idx >= 0) {  // Leaf node
                            fprintf(out_fp, "\t%s", get_output_genome_name(node->genome_idx, opts->short_names));

                            // For parsimony method, use parsimony state vectors; for ALL/SOME, use bit vectors
                            int nd = node->nd_value;
                            int md, mb;

                            if (opts->higher_taxa_with_mbs) {
                                // ALL/SOME method - use bit vectors
                                md = sum_bit_vector(node->md_vec, read_data->trimmed_length);
                                mb = sum_bit_vector(node->mb_vec, read_data->trimmed_length);
                            } else {
                                // Parsimony method - use parsimony state vectors
                                md = count_parsimony_state_1(node->damage_parsimony_state_vec, read_data->trimmed_length);
                                mb = count_parsimony_state_1(node->parsimony_state_vec, read_data->trimmed_length);
                            }

                            if (is_damage_format) {
                                fprintf(out_fp, "\t%d\t%d\t%d", nd, md, mb);
                            } else {
                                fprintf(out_fp, "\t%d", mb + md);
                            }
                        }
                    }
                }

                // Output taxonomic groups (internal nodes)
                fprintf(out_fp, "\t|");

                for (int i = 0; i < node_count; i++) {
                    TreeNode *node = all_nodes[i];
                    if (node->genome_idx < 0 && node->n_children > 0) {  // Internal node with children

                        // For sparse format, output taxid first
                        if (is_sparse) {
                            fprintf(out_fp, "\tT%d:%s", node->taxid, node->rank);
                        }

                        if (opts->higher_taxa_with_mbs) {
                            // Old ALL/SOME method - 5 values for damage, 2 for non-damage
                            int nd = node->nd_value;
                            int md = sum_bit_vector(node->md_vec, read_data->trimmed_length);
                            int mds = sum_bit_vector(node->mds_vec, read_data->trimmed_length);
                            int mb = sum_bit_vector(node->mb_vec, read_data->trimmed_length);
                            int mbs = sum_bit_vector(node->mbs_vec, read_data->trimmed_length);

                            if (is_damage_format) {
                                // For taxonomic groups: nd md mds mb mbs (5 numbers)
                                fprintf(out_fp, "\t%d\t%d\t%d\t%d\t%d", nd, md, mds, mb, mbs);
                            } else {
                                // Without damage: mb+md mbs+mds (2 numbers)
                                fprintf(out_fp, "\t%d\t%d", mb + md, mbs + mds);
                            }
                        } else {
                            // New parsimony method - same format as genomic groups
                            int nd = node->nd_value;
                            int mb = count_parsimony_state_1(node->parsimony_state_vec, read_data->trimmed_length);


                            if (is_damage_format) {
                                // Same as individual genomes: nd md mb (3 numbers)
                                int md = count_parsimony_state_1(node->damage_parsimony_state_vec, read_data->trimmed_length);
                                fprintf(out_fp, "\t%d\t%d\t%d", nd, md, mb);
                            } else {
                                // Without damage: single value (same as individual genomes)
                                fprintf(out_fp, "\t%d", mb);
                            }
                        }
                    }
                }

                fprintf(out_fp, "\n");
            }  // End if (should_output)

                // Free the dynamically allocated nodes array (always needed)
                free(all_nodes);
            }
            
            entry = entry->next;
        }
    }


    if (!fp) {
        fclose(out_fp);
    }
}

void print_usage(const char *prog_name) {
    printf("TaxIdent - Taxonomic Hierarchy BAM Processor\n\n");
    printf("Processes BAM alignment files with taxonomic hierarchy support.\n");
    printf("Calculates mismatch matrices for higher taxonomic groups using NCBI taxonomy.\n\n");
    printf("Usage: %s [OPTIONS]\n\n", prog_name);
    printf("Required arguments:\n");
    printf("  -i, --input FILE/DIR  Input file (BAM) or directory containing BAM files\n");
    printf("  -o, --output FILE     Output file\n\n");
    printf("Optional arguments:\n");
    printf("  -g, --genomes FILE    Text file listing genome names (optional - auto-discovers if not provided)\n");
    printf("  --text                Process text mismatch matrix instead of BAM file\n");
    printf("  --simple              Use simple unoptimized algorithms (for testing/debugging)\n");
    printf("  -q, --min-mapq INT    Minimum mapping quality [0] (BAM mode only)\n");
    printf("  -n, --max-reads INT   Maximum alignments to process [0=all, supports billions] (BAM mode only)\n");
    printf("  -I, --ignore CHAR     Ignore characters until this character in genome names\n");
    printf("  -s, --damage-sites N  Number of damage-susceptible sites at read ends [5]\n");
    printf("  --damage              Enable damage analysis (C->T, G->A in first/last s sites) [ENABLED by default]\n");
    printf("  --no-damage           Disable damage analysis\n");
    printf("  --asymmetric-damage   Use asymmetric damage (C->T in first s, G->A in last s)\n");
    printf("  --skip-indels         Skip alignments containing indels\n");
    printf("  --precise-indels      Use precise calculation for alignments with indels [ENABLED by default]\n");
    printf("  --format FORMAT       Output format: dense, sparse, dense_damage, sparse_damage [sparse_damage]\n");
    printf("  --short-names         Use short genome names (G1, G2, etc.) [ENABLED by default]\n");
    printf("  --full-names          Use full genome names instead of short names\n");
    printf("  --key-file            Create key file mapping short names to full names (redundant - now automatic)\n");
    printf("  --no-key-file         Disable automatic key file creation when using short names\n");
    printf("  --compress            Use compressed output with short IDs\n");
    printf("  --taxid-file FILE     NCBI names.dmp file for taxid mapping (with --compress)\n");
    printf("  --genome-map FILE     Output file for genome ID to name mapping (with --compress)\n");
    printf("  --use-rg              Use RG tag for genome identification instead of reference name\n");
    printf("  --dense               Use dense internal data structures (default: sparse for better performance)\n");
    printf("  --no-redundancy       Assume reads don't appear in multiple BAM files (faster for directory mode)\n");
    printf("  --penalty-mode        Use max+1 penalty for missing alignments instead of -1 (dense formats only)\n");
    printf("\nTaxonomy options:\n");
    printf("  --taxonomy-dir DIR    Directory containing NCBI taxonomy files (nodes.dmp, names.dmp)\n");
    printf("  --acc2taxid FILE      Accession to taxid mapping file\n");
    printf("  --with-higher-taxa    Include higher taxonomic groups in output\n");
    printf("  --consolidate-by-taxid Consolidate genomes by taxid (sparse mode only, requires --acc2taxid)\n");
    printf("  --higher_taxa_with_mbs Use old ALL/SOME method (5 values) instead of parsimony (3 values)\n");
    printf("  --consecutive         Optimize for consecutive alignments (reduces memory usage)\n");
    printf("  --tax-levels LEVELS   Comma-separated taxonomic levels to include (e.g., genus,family,order)\n");
    printf("\nPerformance options:\n");
    printf("  -t, --threads N       Number of decompression threads [0]\n");
    printf("  --parallel-files N    Number of BAM files to process in parallel (directory mode) [1-8]\n");
    printf("  -v, --verbose         Verbose output\n");
    printf("  --silent              Suppress all progress output\n");
    printf("  -h, --help            Show this help message\n\n");
    printf("Examples:\n");
    printf("  BAM mode (default): %s -i multi_aligned.bam -g genome_names.txt -o mismatch_matrix.txt\n", prog_name);
    printf("  BAM directory mode: %s -i /path/to/bam/files/ -g genome_names.txt -o mismatch_matrix.txt\n", prog_name);
    printf("  BAM with ignore: %s -i multi_aligned.bam -g genome_names.txt -I _ -o mismatch_matrix.txt\n", prog_name);
    printf("  Text mode: %s --text -i Mismatch.txt -o processed.txt\n\n", prog_name);
    printf("Text file format: read_id, total_count, mismatch_count_per_genome...\n");
    printf("BAM file requirements:\n");
    printf("  - Reads aligned to multiple reference genomes\n");
    printf("  - Reference names should match genome names\n");
    printf("  - NM tags present for mismatch counts\n");
}

int parse_options(int argc, char **argv, options_t *opts) {
    // Set defaults
    opts->input_file = NULL;
    opts->genome_list = NULL;
    opts->output_file = NULL;
    opts->min_mapq = 0;
    opts->max_reads = 0;
    opts->verbose = 0;
    opts->bam_mode = 1;  // Default to BAM mode
    opts->use_simple_mode = 0;  // Default to optimized mode
    opts->ignore_char = '\0';  // No character to ignore by default
    opts->damage_sites = -1;  // Will be set based on damage mode
    opts->enable_damage = 1;  // Damage analysis ENABLED by default
    opts->output_format = "sparse_damage";  // Default to sparse_damage format
    opts->short_names = 1;  // Use short genome names by default
    opts->create_key_file = 1;  // Auto-create key file when using short names
    opts->silent = 0;  // Show progress by default
    opts->asymmetric_damage = 0;  // Use symmetric damage by default
    opts->skip_indels = 0;  // Process reads with indels by default
    opts->precise_indels = 1;  // Use precise indel handling by default
    opts->compress_output = 0;  // Don't compress output by default
    opts->taxid_file = NULL;
    opts->genome_map_file = NULL;
    opts->parallel_files = 1;  // Process files sequentially by default
    opts->no_redundancy = 0;  // Default: assume reads may appear in multiple files
    
    // Initialize new taxonomy options
    opts->taxonomy_dir = NULL;
    opts->acc2taxid_file = NULL;
    opts->with_higher_taxa = 0;
    opts->higher_taxa_with_mbs = 0;
    opts->consecutive_mode = 0;
    opts->tax_levels = NULL;
    opts->test_tree = 0;
    opts->enforce_dense_strict = 0;

    static struct option long_options[] = {
        {"input", required_argument, 0, 'i'},
        {"bam", no_argument, 0, 1000},  // Use high number for --bam flag (kept for compatibility)
        {"text", no_argument, 0, 1005},  // Switch to text mode
        {"simple", no_argument, 0, 1001},  // Use simple unoptimized algorithms
        {"genomes", required_argument, 0, 'g'},
        {"output", required_argument, 0, 'o'},
        {"min-mapq", required_argument, 0, 'q'},
        {"max-reads", required_argument, 0, 'n'},
        {"ignore", required_argument, 0, 'I'},
        {"damage-sites", required_argument, 0, 's'},
        {"threads", required_argument, 0, 't'},  // Multi-threading support
        {"damage", no_argument, 0, 1002},
        {"no-damage", no_argument, 0, 1013},
        {"format", required_argument, 0, 1003},
        {"short-names", no_argument, 0, 1004},
        {"full-names", no_argument, 0, 1014},
        {"key-file", no_argument, 0, 1030},
        {"no-key-file", no_argument, 0, 1031},
        {"verbose", no_argument, 0, 'v'},
        {"silent", no_argument, 0, 1006},
        {"asymmetric-damage", no_argument, 0, 1007},
        {"skip-indels", no_argument, 0, 1008},
        {"precise-indels", no_argument, 0, 1009},
        {"compress", no_argument, 0, 1010},
        {"taxid-file", required_argument, 0, 1011},
        {"genome-map", required_argument, 0, 1012},
        {"parallel-files", required_argument, 0, 1015},
        {"use-rg", no_argument, 0, 1016},  // Use RG tag instead of reference name
        {"dense", no_argument, 0, 1017},  // Use dense internal data structures instead of sparse
        {"no-redundancy", no_argument, 0, 1018},  // Assume reads don't appear in multiple BAM files
        {"taxonomy-dir", required_argument, 0, 1019},  // Directory with NCBI taxonomy files
        {"acc2taxid", required_argument, 0, 1020},  // Accession to taxid mapping file
        {"with-higher-taxa", no_argument, 0, 1021},  // Include higher taxonomic groups
        {"higher_taxa_with_mbs", no_argument, 0, 1025},  // Use old ALL/SOME method (5 values)
        {"consecutive", no_argument, 0, 1022},  // Optimize for consecutive alignments
        {"tax-levels", required_argument, 0, 1023},  // Taxonomic levels to include
        {"test-tree", no_argument, 0, 1024},  // Test tree parsing and exit
        {"enforce-dense_strict", no_argument, 0, 1026},  // Only output reads with alignments to ALL genomes/taxids/readgroups
        {"penalty-mode", no_argument, 0, 1027},  // Use max+1 penalty instead of -1 for missing alignments
        {"consolidate-by-taxid", no_argument, 0, 1028},  // Consolidate by taxid in sparse mode
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "i:g:o:q:n:I:s:t:vh", long_options, NULL)) != -1) {
        switch (c) {
            case 'i':
                opts->input_file = optarg;
                break;
            case 1000:  // --bam flag (kept for compatibility)
                opts->bam_mode = 1;
                break;
            case 1005:  // --text flag
                opts->bam_mode = 0;
                break;
            case 1001:  // --simple flag
                opts->use_simple_mode = 1;
                break;
            case 'g':
                opts->genome_list = optarg;
                break;
            case 'o':
                opts->output_file = optarg;
                break;
            case 'q':
                opts->min_mapq = atoi(optarg);
                break;
            case 'n':
                opts->max_reads = atoll(optarg);  // Use atoll for long long
                break;
            case 'I':
                if (optarg && strlen(optarg) > 0) {
                    opts->ignore_char = optarg[0];
                } else {
                    fprintf(stderr, "Error: -I option requires a character argument\n");
                    return -1;
                }
                break;
            case 's':
                opts->damage_sites = atoi(optarg);
                if (opts->damage_sites < 1 || opts->damage_sites > 50) {
                    fprintf(stderr, "Error: damage sites must be between 1 and 50\n");
                    return -1;
                }
                break;
            case 1002:  // --damage flag
                opts->enable_damage = 1;
                break;
            case 1013:  // --no-damage flag
                opts->enable_damage = 0;
                break;
            case 1003:  // --format flag
                if (strcmp(optarg, "dense") == 0 || strcmp(optarg, "sparse") == 0 ||
                    strcmp(optarg, "dense_damage") == 0 || strcmp(optarg, "sparse_damage") == 0) {
                    opts->output_format = optarg;
                } else {
                    fprintf(stderr, "Error: Invalid format. Use dense, sparse, dense_damage, or sparse_damage\n");
                    return -1;
                }
                break;
            case 1004:  // --short-names flag
                opts->short_names = 1;
                break;
            case 1014:  // --full-names flag
                opts->short_names = 0;
                opts->create_key_file = 0;  // No key file needed with full names
                break;
            case 1030:  // --key-file flag
                opts->create_key_file = 1;
                break;
            case 1031:  // --no-key-file flag
                opts->create_key_file = 0;
                break;
            case 't':
                opts->num_threads = atoi(optarg);
                break;
            case 'v':
                opts->verbose = 1;
                break;
            case 1006:  // --silent flag
                opts->silent = 1;
                break;
            case 1007:  // --asymmetric-damage flag
                opts->asymmetric_damage = 1;
                break;
            case 1008:  // --skip-indels flag
                opts->skip_indels = 1;
                break;
            case 1009:  // --precise-indels flag
                opts->precise_indels = 1;
                break;
            case 1010:  // --compress flag
                opts->compress_output = 1;
                break;
            case 1011:  // --taxid-file
                opts->taxid_file = optarg;
                break;
            case 1012:  // --genome-map
                opts->genome_map_file = optarg;
                break;
            case 1015:  // --parallel-files
                opts->parallel_files = atoi(optarg);
                if (opts->parallel_files < 1) opts->parallel_files = 1;
                if (opts->parallel_files > 8) opts->parallel_files = 8; // Cap at 8 for safety
                break;
            case 1016:  // --use-rg
                opts->use_rg_tag = 1;
                break;
            case 1017:  // --dense
                opts->use_dense = 1;
                break;
            case 1018:  // --no-redundancy
                opts->no_redundancy = 1;
                break;
            case 1019:  // --taxonomy-dir
                opts->taxonomy_dir = optarg;
                break;
            case 1020:  // --acc2taxid
                opts->acc2taxid_file = optarg;
                break;
            case 1021:  // --with-higher-taxa
                opts->with_higher_taxa = 1;
                break;
            case 1025:  // --higher_taxa_with_mbs
                opts->higher_taxa_with_mbs = 1;
                break;
            case 1022:  // --consecutive
                opts->consecutive_mode = 1;
                break;
            case 1023:  // --tax-levels
                opts->tax_levels = optarg;
                break;
            case 1024:  // --test-tree
                opts->test_tree = 1;
                break;
            case 1026:  // --enforce-dense_strict
                opts->enforce_dense_strict = 1;
                break;
            case 1027:  // --penalty-mode
                opts->use_penalty_mode = 1;
                break;
            case 1028:  // --consolidate-by-taxid
                opts->consolidate_by_taxid = 1;
                break;
            case 'h':
                print_usage(argv[0]);
                exit(0);
            default:
                return -1;
        }
    }

    // Set default damage_sites based on damage mode if not explicitly set
    if (opts->damage_sites == -1) {
        if (opts->enable_damage) {
            opts->damage_sites = 5;  // Default to 5 for damage mode
        } else {
            opts->damage_sites = 0;  // Default to 0 (no trimming) for non-damage mode
        }
    }
    
    // Check required arguments
    if (!opts->input_file || !opts->output_file) {
        fprintf(stderr, "Error: Missing required input or output file\n");
        print_usage(argv[0]);
        return -1;
    }

    // Check mode-specific requirements
    if (opts->bam_mode && !opts->genome_list) {
        // Auto-discovery will be performed
        if (opts->verbose) {
            fprintf(stderr, "Info: No genome list provided, will auto-discover genomes from BAM files\n");
        }
    }
    
    // For text mode, genome list should not be required
    if (!opts->bam_mode && opts->genome_list) {
        fprintf(stderr, "Warning: Genome list ignored in text mode\n");
    }
    
    // Handle conflicting options
    if (opts->verbose && opts->silent) {
        fprintf(stderr, "Warning: Both --verbose and --silent specified. Using --silent mode.\n");
        opts->verbose = 0;  // Silent takes precedence
    }

#ifndef WITH_HTSLIB
    if (opts->bam_mode) {
        fprintf(stderr, "Error: BAM mode not available - program compiled without HTSlib\n");
        return -1;
    }
#endif

    return 0;
}

int load_genome_list(const char *filename, char ignore_char, int use_simple_mode) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open genome list file: %s\n", filename);
        return -1;
    }

    char line[MAX_NAME_LEN];
    n_genomes = 0;

    while (fgets(line, sizeof(line), fp)) {
        // Grow arrays if needed
        if (n_genomes >= max_genomes_allocated) {
            grow_genome_arrays();
        }
        // Remove newline
        line[strcspn(line, "\n")] = 0;
        
        // Skip empty lines and comments
        if (line[0] == '\0' || line[0] == '#') continue;
        
        // Keep full genome names from file - don't process with ignore_char
        // Only process RG values from BAM file with ignore_char
        strncpy(genome_names[n_genomes], line, MAX_NAME_LEN - 1);
        genome_names[n_genomes][MAX_NAME_LEN - 1] = '\0';
        
        // Add to hash table for fast lookup (only in optimized mode)
        if (!use_simple_mode) {
            add_genome_to_hash(line, n_genomes);
        }
        
        // Generate compressed ID for this genome
        generate_genome_compressed_id(n_genomes, line);
        
        n_genomes++;
    }

    fclose(fp);

    if (n_genomes == 0) {
        fprintf(stderr, "Error: No valid genome names found\n");
        return -1;
    }

    return n_genomes;
}

const char *extract_genome_name(const char *name, char ignore_char) {
    if (ignore_char == '\0') {
        return name;  // No character to ignore
    }
    
    // Simple logic: find first occurrence of ignore_char and return everything after it
    // This matches the Python logic: rg_value.split(ignore_char, 1)[1]
    const char *pos = strchr(name, ignore_char);
    if (pos) {
        return pos + 1;  // Return everything after the first occurrence
    }
    
    return name;  // Character not found, return original name
}

// Simple replacement for find_or_add_genome() - stores taxid strings instead of genome names
int find_or_add_taxid(const char *ref_name, char ignore_char) {
    // Extract genome/accession name (same logic as current)
    const char *genome_name = extract_genome_name(ref_name, ignore_char);
    if (!genome_name) return -1;
    
    // Extract base accession (remove version)
    char base_accession[MAX_NAME_LEN];
    strncpy(base_accession, genome_name, MAX_NAME_LEN - 1);
    base_accession[MAX_NAME_LEN - 1] = '\0';
    char *dot = strchr(base_accession, '.');
    if (dot) *dot = '\0';
    
    // Look up taxid for this accession
    int taxid = get_taxid_for_accession(base_accession);
    if (taxid <= 0) {
        // Track alignments with unmapped accessions
        invalid_taxid_alignments++;
        
        // Write accession to report file (create file if needed)
        if (!unmapped_report_file) {
            char report_filename[512];
            if (global_output_filename) {
                // Remove .txt extension if present and add missing_taxid
                char base_name[512];
                strncpy(base_name, global_output_filename, sizeof(base_name) - 1);
                base_name[sizeof(base_name) - 1] = '\0';
                
                // Remove .txt extension if present
                int len = strlen(base_name);
                if (len > 4 && strcmp(base_name + len - 4, ".txt") == 0) {
                    base_name[len - 4] = '\0';
                }
                
                snprintf(report_filename, sizeof(report_filename), "%s.missing_taxid.txt", base_name);
            } else {
                snprintf(report_filename, sizeof(report_filename), "taxonomy_report.missing_taxid.txt");
            }
            unmapped_report_file = fopen(report_filename, "w");
            if (unmapped_report_file) {
                fprintf(unmapped_report_file, "# Accessions without taxonomic mapping\n");
                fprintf(unmapped_report_file, "# Found during TaxIdent processing\n");
                fprintf(unmapped_report_file, "accession\n");
            }
        }
        
        if (unmapped_report_file) {
            fprintf(unmapped_report_file, "%s\n", base_accession);
        }
        
        taxid = -1;  // Unknown taxid
    }
    
    // Create taxid string identifier
    char taxid_str[MAX_NAME_LEN];
    snprintf(taxid_str, MAX_NAME_LEN, "T%d", taxid);

    int result = find_or_add_genome(taxid_str);

    // Use existing find_or_add_genome logic but with taxid strings
    // This automatically provides best-alignment-per-taxid:
    // - Multiple genomes with same taxid → same taxid_str → same index returned
    // - Existing best alignment logic in BAM loop handles NM comparison automatically
    return result;
}

// Extract taxid number from taxid string "T12345" -> 12345
int get_taxid_from_index(int taxid_idx) {
    if (taxid_idx < 0 || taxid_idx >= n_genomes) {
        fprintf(stderr, "FATAL ERROR: get_taxid_from_index bounds violation\n");
        fprintf(stderr, "  Attempted index: %d, n_genomes: %d\n", taxid_idx, n_genomes);
        fprintf(stderr, "  This indicates stored genome/taxid index is invalid.\n");
        fprintf(stderr, "  Root cause: likely mismatch between stored indices and current genome array size.\n");
        abort();
    }
    if (genome_names[taxid_idx][0] != 'T') {
        fprintf(stderr, "ERROR: get_taxid_from_index non-taxid string: idx=%d, name='%s'\n", taxid_idx, genome_names[taxid_idx]);
        return -1;
    }
    return atoi(&genome_names[taxid_idx][1]);
}

// Get appropriate genome/taxid name for output (short or full)
const char *get_output_genome_name(int genome_index, int use_short_names) {
    if (using_taxid_mode || (global_opts && global_opts->consolidate_by_taxid)) {
        // In taxid mode or consolidate mode, always return actual taxids (ignore use_short_names)
        return genome_names[genome_index];  // Contains taxid strings like "T4565"
    } else {
        // In genome mode, use traditional short/full name logic
        static char short_name[16];
        if (use_short_names) {
            snprintf(short_name, sizeof(short_name), "G%d", genome_index + 1);
            return short_name;
        } else {
            return genome_names[genome_index];
        }
    }
}

const char *get_output_genome_name_compressed(int genome_index, options_t *opts) {
    if (opts->compress_output) {
        return genome_compressed_ids[genome_index];
    } else if (opts->short_names) {
        static char short_name[16];
        snprintf(short_name, sizeof(short_name), "G%d", genome_index + 1);
        return short_name;
    } else {
        return genome_names[genome_index];
    }
}

int find_genome_index(const char *ref_name, char ignore_char) {
    const char *genome_name = extract_genome_name(ref_name, ignore_char);
    for (int i = 0; i < n_genomes; i++) {
        if (strcmp(genome_name, genome_names[i]) == 0) {
            return i;
        }
    }
    return -1;  // Not found
}

#ifdef WITH_HTSLIB

// Fast number parsing using multiplication-free technique
static inline int fast_parse_number(const unsigned char **ptr) {
    int result = 0;
    const unsigned char *p = *ptr;
    
    // Unroll first few digits for common cases
    if (char_types[*p] == CHAR_DIGIT) {
        result = *p++ - '0';
        if (char_types[*p] == CHAR_DIGIT) {
            result = (result << 3) + (result << 1) + (*p++ - '0'); // result * 10
            if (char_types[*p] == CHAR_DIGIT) {
                result = (result << 3) + (result << 1) + (*p++ - '0');
                // Continue with loop for longer numbers
                while (char_types[*p] == CHAR_DIGIT) {
                    result = (result << 3) + (result << 1) + (*p++ - '0');
                }
            }
        }
    }
    
    *ptr = p;
    return result;
}

// Optimized MD tag parsing using lookup table and fast number parsing
int parse_md_tag(const char *md_string, int *mismatch_positions, int max_mismatches) {
    if (!md_string) return 0;
    
    int pos = 0;
    int num_mismatches = 0;
    const unsigned char *ptr = (const unsigned char *)md_string;
    
    while (*ptr && num_mismatches < max_mismatches) {
        unsigned char char_type = char_types[*ptr];
        
        if (char_type == CHAR_DIGIT) {
            // Parse number of matching bases - optimized with fast parser
            pos += fast_parse_number(&ptr);
        } else if (*ptr == '^') {
            // Deletion - skip the ^ and the deleted bases
            ptr++; // skip ^
            while (*ptr && char_types[*ptr] == CHAR_ALPHA) {
                ptr++; // skip deleted bases
            }
            // Deletions don't advance position in read
        } else if (char_type == CHAR_ALPHA) {
            // Mismatch base - store position for future position-specific analysis
            mismatch_positions[num_mismatches++] = pos;
            pos++;
            ptr++;
        } else {
            ptr++; // skip unknown characters
        }
    }
    
    return num_mismatches;
}

// Check if a mismatch is a damage mismatch (C->T or G->A)
int is_damage_mismatch(char ref_base, char read_base) {
    return ((ref_base == 'C' || ref_base == 'c') && (read_base == 'T' || read_base == 't')) ||
           ((ref_base == 'G' || ref_base == 'g') && (read_base == 'A' || read_base == 'a'));
}

// Parse MD tag and extract reference bases at mismatch positions
typedef struct {
    int pos;       // Position in read
    char ref_base; // Reference base at this position
    char read_base; // Read base at this position  
} mismatch_info_t;

int parse_md_with_bases(const char *md_string, uint8_t *read_seq, mismatch_info_t *mismatches, int max_mismatches) {
    if (!md_string) return 0;
    
    int pos = 0;
    int num_mismatches = 0;
    const char *ptr = md_string;
    
    while (*ptr && num_mismatches < max_mismatches) {
        if (isdigit(*ptr)) {
            // Parse number of matching bases
            int matches = 0;
            while (isdigit(*ptr)) {
                matches = matches * 10 + (*ptr - '0');
                ptr++;
            }
            pos += matches;
        } else if (*ptr == '^') {
            // Deletion - skip the ^ and the deleted bases
            ptr++; // skip ^
            while (*ptr && isalpha(*ptr)) {
                ptr++; // skip deleted bases
            }
            // Deletions don't advance position in read
        } else if (isalpha(*ptr)) {
            // Mismatch - the character is the reference base
            if (pos < 10000) { // Safety check
                mismatches[num_mismatches].pos = pos;
                mismatches[num_mismatches].ref_base = *ptr;
                mismatches[num_mismatches].read_base = seq_nt16_str[bam_seqi(read_seq, pos)];
                num_mismatches++;
            }
            pos++;
            ptr++;
        } else {
            ptr++; // skip unknown characters
        }
    }
    
    return num_mismatches;
}

// Count damage-susceptible sites (C or G) in a given range of the reference sequence
int count_damage_sites(const char *ref_seq, int start, int end) {
    int count = 0;
    for (int i = start; i < end; i++) {
        char base = ref_seq[i];
        if (base == 'C' || base == 'c' || base == 'G' || base == 'g') {
            count++;
        }
    }
    return count;
}

// Calculate mismatches for reads with indels using CIGAR and MD tag
int calculate_mismatches_with_indels(bam1_t *read, sam_hdr_t *header, options_t *opts, 
                                    damage_stats_t *damage_stats, const char *md_string) {
    int read_len = read->core.l_qseq;
    uint8_t *seq = bam_get_seq(read);
    uint32_t *cigar = bam_get_cigar(read);
    
    // Initialize damage stats if provided
    if (damage_stats) {
        damage_stats->nd = 0;
        damage_stats->md = 0;
        damage_stats->mb = 0;
    }
    
    // For damage analysis mode, don't trim - analyze entire read
    int trim_start = opts->enable_damage ? 0 : 5;
    int trim_end = opts->enable_damage ? 0 : 5;
    int analysis_start = trim_start;
    int analysis_end = read_len - trim_end;
    int analysis_length = analysis_end - analysis_start;
    
    // Build position mapping between read and reference coordinates
    int *read_to_ref = (int*)calloc(read_len, sizeof(int));
    int read_pos = 0;
    int ref_pos = 0;
    
    
    // Parse CIGAR to build coordinate mapping
    for (int i = 0; i < read->core.n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        
        switch (op) {
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:
                // Match/mismatch - positions advance in both
                for (int j = 0; j < len; j++) {
                    if (read_pos < read_len) {
                        read_to_ref[read_pos] = ref_pos;
                    }
                    read_pos++;
                    ref_pos++;
                }
                break;
            case BAM_CINS:
                // Insertion - read advances, ref doesn't
                for (int j = 0; j < len; j++) {
                    if (read_pos < read_len) {
                        read_to_ref[read_pos] = -1;  // No ref position
                    }
                    read_pos++;
                }
                break;
            case BAM_CDEL:
            case BAM_CREF_SKIP:
                // Deletion - ref advances, read doesn't
                ref_pos += len;
                break;
            case BAM_CSOFT_CLIP:
                // Soft clip - read advances, ref doesn't
                for (int j = 0; j < len; j++) {
                    if (read_pos < read_len) {
                        read_to_ref[read_pos] = -1;  // No ref position
                    }
                    read_pos++;
                }
                break;
            case BAM_CHARD_CLIP:
                // Hard clip - neither advances
                break;
        }
    }
    
    
    // Now parse MD tag to find mismatches
    const char *md_ptr = md_string;
    int md_ref_pos = 0;
    int total_mismatches = 0;
    
    
    while (*md_ptr) {
        if (isdigit(*md_ptr)) {
            // Parse number of matching bases
            int matches = 0;
            while (isdigit(*md_ptr)) {
                matches = matches * 10 + (*md_ptr - '0');
                md_ptr++;
            }
            
            // For damage analysis, count C/G sites in matching regions
            if (opts->enable_damage && damage_stats) {
                for (int i = 0; i < matches; i++) {
                    // Find read position for this ref position
                    int curr_read_pos = -1;
                    for (int j = 0; j < read_len; j++) {
                        if (read_to_ref[j] == md_ref_pos + i) {
                            curr_read_pos = j;
                            break;
                        }
                    }
                    
                    
                    if (curr_read_pos >= 0) {
                        char read_base = seq_nt16_str[bam_seqi(seq, curr_read_pos)];
                        
                        // Count damage-susceptible sites based on mode
                        if (opts->asymmetric_damage) {
                            // Asymmetric mode: C in first s positions, G in last s positions
                            if (curr_read_pos < opts->damage_sites && (read_base == 'C' || read_base == 'c')) {
                                damage_stats->nd++;
                            } else if (curr_read_pos >= read_len - opts->damage_sites && (read_base == 'G' || read_base == 'g')) {
                                damage_stats->nd++;
                            }
                        } else {
                            // Symmetric mode: C or G in both first s and last s positions
                            if (curr_read_pos < opts->damage_sites || curr_read_pos >= read_len - opts->damage_sites) {
                                if (read_base == 'C' || read_base == 'c' || read_base == 'G' || read_base == 'g') {
                                    damage_stats->nd++;
                                }
                            }
                        }
                    }
                }
            }
            md_ref_pos += matches;
            
        } else if (*md_ptr == '^') {
            // Deletion in reference - skip
            md_ptr++; // skip ^
            while (*md_ptr && isalpha(*md_ptr)) {
                md_ptr++;
                md_ref_pos++;
            }
            
        } else if (isalpha(*md_ptr)) {
            // Mismatch - *md_ptr is the reference base
            char ref_base = *md_ptr;
            
            // Find read position for this ref position
            int curr_read_pos = -1;
            for (int j = 0; j < read_len; j++) {
                if (read_to_ref[j] == md_ref_pos) {
                    curr_read_pos = j;
                    break;
                }
            }
            
            if (curr_read_pos >= 0) {
                char read_base = seq_nt16_str[bam_seqi(seq, curr_read_pos)];
                
                // Count mismatch if in analysis region
                if (!opts->enable_damage) {
                    if (curr_read_pos >= analysis_start && curr_read_pos < analysis_end) {
                        total_mismatches++;
                    }
                } else {
                    total_mismatches++;
                }
                
                // For damage analysis
                if (opts->enable_damage && damage_stats) {
                    if (opts->asymmetric_damage) {
                        // Asymmetric damage mode
                        if (curr_read_pos < opts->damage_sites) {
                            // First s sites - check if reference has C (damage site)
                            if (ref_base == 'C' || ref_base == 'c') {
                                damage_stats->nd++;  // Count as damage site
                                // Check if it's a C->T damage
                                if (read_base == 'T' || read_base == 't') {
                                    damage_stats->md++;
                                } else {
                                    damage_stats->mb++;
                                }
                            } else {
                                // Non-damage site, but still a mismatch in damage region
                                damage_stats->mb++;
                            }
                        } else if (curr_read_pos >= read_len - opts->damage_sites) {
                            // Last s sites - check if reference has G (damage site)
                            if (ref_base == 'G' || ref_base == 'g') {
                                damage_stats->nd++;  // Count as damage site
                                // Check if it's a G->A damage
                                if (read_base == 'A' || read_base == 'a') {
                                    damage_stats->md++;
                                } else {
                                    damage_stats->mb++;
                                }
                            } else {
                                // Non-damage site, but still a mismatch in damage region
                                damage_stats->mb++;
                            }
                        } else {
                            // Middle region - all mismatches are background
                            damage_stats->mb++;
                        }
                    } else {
                        // Symmetric damage mode
                        if (curr_read_pos < opts->damage_sites || curr_read_pos >= read_len - opts->damage_sites) {
                            // In damage-prone regions - check if reference has C or G (damage site)
                            if (ref_base == 'C' || ref_base == 'c' || ref_base == 'G' || ref_base == 'g') {
                                damage_stats->nd++;  // Count as damage site
                                // Check if it's a damage-type mismatch
                                if (is_damage_mismatch(ref_base, read_base)) {
                                    damage_stats->md++;
                                } else {
                                    damage_stats->mb++;
                                }
                            } else {
                                // Non-damage site, but still a mismatch in damage region
                                damage_stats->mb++;
                            }
                        } else {
                            // Middle region - all mismatches are background
                            damage_stats->mb++;
                        }
                    }
                }
            }
            
            md_ref_pos++;
            md_ptr++;
        } else {
            md_ptr++; // skip unknown characters
        }
    }
    
    free(read_to_ref);
    
    if (damage_stats && damage_stats->md > damage_stats->nd) {
        fprintf(stderr, "\nERROR: md > nd in calculate_mismatches_with_indels\n");
        fprintf(stderr, "  Read: %s\n", bam_get_qname(read));
        fprintf(stderr, "  nd=%d, md=%d, mb=%d\n", damage_stats->nd, damage_stats->md, damage_stats->mb);
        fprintf(stderr, "  Read length: %d\n", read->core.l_qseq);
        fprintf(stderr, "  Damage sites: %d\n", opts->damage_sites);
    }
    
    return total_mismatches;
}

int calculate_trimmed_mismatches(bam1_t *read, sam_hdr_t *header, options_t *opts, damage_stats_t *damage_stats, int *indel_count) {
    // Skip reads shorter than minimum length
    if (read->core.l_qseq < MIN_READ_LENGTH) {
        return -1;  // Indicates read should be skipped
    }
    
    int read_len = read->core.l_qseq;
    uint8_t *seq = bam_get_seq(read);
    
    // Initialize damage stats if provided
    if (damage_stats) {
        damage_stats->nd = 0;
        damage_stats->md = 0;
        damage_stats->mb = 0;
    }
    
    // Determine analysis region based on mode
    int analysis_start = 0;
    int analysis_end = read_len;
    int analysis_length = read_len;
    
    // In non-damage mode, use damage_sites for trimming (default 0 = no trimming)
    // In damage mode, analyze entire read
    if (!opts->enable_damage && opts->damage_sites > 0) {
        analysis_start = opts->damage_sites;
        analysis_end = read_len - opts->damage_sites;
        analysis_length = analysis_end - analysis_start;
        
        // If trimmed read would be too short, skip it
        if (analysis_length <= 0) {
            return -1;
        }
    }
    
    // Get MD tag which contains reference bases at mismatch positions
    uint8_t *md_tag = bam_aux_get(read, "MD");
    if (!md_tag) {
        if (opts->enable_damage && damage_stats) {
            // MD tag is required for damage analysis
            char *read_name = bam_get_qname(read);
            fprintf(stderr, "\nERROR: MD tag is required for damage analysis but not found.\n");
            fprintf(stderr, "Missing MD tag at read: %s\n", read_name);
            fprintf(stderr, "Please run 'samtools calmd -b your.bam ref.fasta > your_with_md.bam' to add MD tags.\n");
            exit(1);
        }
        // For non-damage analysis, try NM tag fallback
        uint8_t *nm_tag = bam_aux_get(read, "NM");
        if (nm_tag) {
            int total_nm = bam_aux2i(nm_tag);
            // Scale mismatches proportionally to trimmed length
            return (total_nm * analysis_length) / read_len;
        }
        return 0;
    }
    
    const char *md_string = bam_aux2Z(md_tag);
    
    // Check if CIGAR contains indels - this affects position mapping
    uint32_t *cigar = bam_get_cigar(read);
    int has_indels = 0;
    int total_indels = 0;
    
    // Fast path: If only one CIGAR operation and it's a match, no indels
    if (read->core.n_cigar == 1 && bam_cigar_op(cigar[0]) == BAM_CMATCH) {
        has_indels = 0;
        total_indels = 0;
    } else {
        // Need to check all CIGAR operations
        for (int i = 0; i < read->core.n_cigar; i++) {
            int op = bam_cigar_op(cigar[i]);
            int len = bam_cigar_oplen(cigar[i]);
            if (op == BAM_CINS || op == BAM_CDEL || op == BAM_CREF_SKIP) {
                has_indels = 1;
                total_indels += len;  // Count total indel bases
            }
        }
    }

    // Return indel count if requested
    if (indel_count) {
        *indel_count = total_indels;
    }
    
    // Skip alignments with indels if requested
    if (has_indels && opts->skip_indels) {
        return -1;  // Signal to skip this alignment
    }
    
    // If no indels, we can parse MD tag directly (fast path)
    if (!has_indels) {
        const char *ptr = md_string;
        int pos = 0;
        int total_mismatches = 0;
        
        
        while (*ptr) {
            if (isdigit(*ptr)) {
                // Parse number of matching bases
                int matches = 0;
                while (isdigit(*ptr)) {
                    matches = matches * 10 + (*ptr - '0');
                    ptr++;
                }
                
                // For damage analysis, count C/G sites in matching regions
                if (opts->enable_damage && damage_stats) {
                    for (int i = 0; i < matches && pos + i < read_len; i++) {
                        int curr_pos = pos + i;
                        char read_base = seq_nt16_str[bam_seqi(seq, curr_pos)];


                        // Count damage-susceptible sites based on mode
                        if (opts->asymmetric_damage) {
                            // Asymmetric mode: C in first s positions, G in last s positions
                            if (curr_pos < opts->damage_sites && (read_base == 'C' || read_base == 'c')) {
                                damage_stats->nd++;
                            } else if (curr_pos >= read_len - opts->damage_sites && (read_base == 'G' || read_base == 'g')) {
                                damage_stats->nd++;
                            }
                        } else {
                            // Symmetric mode: C or G in both first s and last s positions
                            if (curr_pos < opts->damage_sites || curr_pos >= read_len - opts->damage_sites) {
                                if (read_base == 'C' || read_base == 'c' || read_base == 'G' || read_base == 'g') {
                                    damage_stats->nd++;
                                }
                            }
                        }
                    }
                }
                pos += matches;
                
            } else if (*ptr == '^') {
                // Deletion in reference - skip
                ptr++; // skip ^
                while (*ptr && isalpha(*ptr)) {
                    ptr++;
                }
                
            } else if (isalpha(*ptr)) {
                // Mismatch - *ptr is the reference base
                char ref_base = *ptr;
                char read_base = seq_nt16_str[bam_seqi(seq, pos)];

                // Count mismatches in analysis region
                if (pos >= analysis_start && pos < analysis_end) {
                    total_mismatches++;
                }

                // For damage analysis
                if (opts->enable_damage && damage_stats) {
                    if (opts->asymmetric_damage) {
                        // Asymmetric damage mode
                        if (pos < opts->damage_sites) {
                            // First s sites - check if reference has C (damage site)
                            if (ref_base == 'C' || ref_base == 'c') {
                                damage_stats->nd++;  // Count as damage site
                                if (strcmp(bam_get_qname(read), "M_A00706:957:HCMMCDSXF:2:1101:5132:18834") == 0) {
                                    fprintf(stderr, "  -> nd++ at mismatch (first region, ref=C), nd now %d\n", damage_stats->nd);
                                }
                                // Check if it's a C->T damage
                                if (read_base == 'T' || read_base == 't') {
                                    damage_stats->md++;
                                } else {
                                    damage_stats->mb++;
                                }
                            } else {
                                // Non-damage site, but still a mismatch in damage region
                                damage_stats->mb++;
                            }
                        } else if (pos >= read_len - opts->damage_sites) {
                            // Last s sites - check if reference has G (damage site)
                            if (ref_base == 'G' || ref_base == 'g') {
                                damage_stats->nd++;  // Count as damage site
                                if (strcmp(bam_get_qname(read), "M_A00706:957:HCMMCDSXF:2:1101:5132:18834") == 0) {
                                    fprintf(stderr, "  -> nd++ at mismatch (last region, ref=G), nd now %d\n", damage_stats->nd);
                                }
                                // Check if it's a G->A damage
                                if (read_base == 'A' || read_base == 'a') {
                                    damage_stats->md++;
                                } else {
                                    damage_stats->mb++;
                                }
                            } else {
                                // Non-damage site, but still a mismatch in damage region
                                damage_stats->mb++;
                            }
                        } else {
                            // Middle region - all mismatches are background
                            damage_stats->mb++;
                        }
                    } else {
                        // Symmetric damage mode
                        if (pos < opts->damage_sites || pos >= read_len - opts->damage_sites) {
                            // In damage-prone regions - check if reference has C or G (damage site)
                            if (ref_base == 'C' || ref_base == 'c' || ref_base == 'G' || ref_base == 'g') {
                                damage_stats->nd++;  // Count as damage site
                                if (strcmp(bam_get_qname(read), "M_A00706:957:HCMMCDSXF:2:1101:5132:18834") == 0) {
                                    fprintf(stderr, "  -> nd++ at mismatch (symmetric, ref=%c), nd now %d\n", ref_base, damage_stats->nd);
                                }
                                // Check if it's a damage-type mismatch
                                if (is_damage_mismatch(ref_base, read_base)) {
                                    damage_stats->md++;
                                } else {
                                    damage_stats->mb++;
                                }
                            } else {
                                // Non-damage site, but still a mismatch in damage region
                                damage_stats->mb++;
                            }
                        } else {
                            // Middle region - all mismatches are background
                            damage_stats->mb++;
                        }
                    }
                }
                
                pos++;
                ptr++;
            } else {
                ptr++; // skip unknown characters
            }
        }
        
        if (damage_stats && damage_stats->md > damage_stats->nd) {
            fprintf(stderr, "\nERROR: md > nd in calculate_trimmed_mismatches\n");
            fprintf(stderr, "  Read: %s\n", bam_get_qname(read));
            fprintf(stderr, "  nd=%d, md=%d, mb=%d\n", damage_stats->nd, damage_stats->md, damage_stats->mb);
            fprintf(stderr, "  Read length: %d\n", read->core.l_qseq);
            fprintf(stderr, "  Damage sites: %d\n", opts->damage_sites);
        }


        return total_mismatches;
    }
    
    // Complex case: has indels, need to use CIGAR for position mapping
    if (opts->precise_indels) {
        // Precise calculation using CIGAR and MD tag together
        return calculate_mismatches_with_indels(read, header, opts, damage_stats, md_string);
    } else {
        // For damage analysis, we cannot proceed without proper MD tag parsing
        if (opts->enable_damage && damage_stats) {
            char *read_name = bam_get_qname(read);
            fprintf(stderr, "\nERROR: Complex indel alignment found but --precise-indels not enabled.\n");
            fprintf(stderr, "Read with indels: %s\n", read_name);
            fprintf(stderr, "For damage analysis with indel-containing reads, use --precise-indels option.\n");
            fprintf(stderr, "Alternatively, ensure all reads have proper MD tags.\n");
            exit(1);
        }
        // Fall back to simple NM-based estimate for non-damage analysis
        uint8_t *nm_tag = bam_aux_get(read, "NM");
        if (nm_tag) {
            int total_nm = bam_aux2i(nm_tag);
            // Scale mismatches proportionally to trimmed length
            return (total_nm * analysis_length) / read_len;
        }
        
        return 0; // No mismatch information available
    }
}

void process_bam_file(options_t *opts) {
    // Open BAM file
    samFile *bam_fp = sam_open(opts->input_file, "r");
    if (!bam_fp) {
        fprintf(stderr, "Error: Cannot open BAM file: %s\n", opts->input_file);
        exit(1);
    }

    // Apply I/O optimizations
    if (opts->num_threads > 0) {
        // Sanity check for thread count
        if (opts->num_threads > 256) {
            fprintf(stderr, "Warning: Invalid thread count %d, using 4 threads instead\n", opts->num_threads);
            opts->num_threads = 4;
        }
        if (hts_set_threads(bam_fp, opts->num_threads) < 0) {
            fprintf(stderr, "Warning: Failed to set %d threads\n", opts->num_threads);
        } else if (opts->verbose) {
            fprintf(stderr, "Info: Using %d decompression threads\n", opts->num_threads);
        }
    }
    
    // Set larger cache size (128MB)
    size_t cache_size = 128 * 1024 * 1024;
    hts_set_cache_size(bam_fp, cache_size);
    if (opts->verbose) {
        fprintf(stderr, "Info: Set cache size to %zu MB\n", cache_size / (1024*1024));
    }

    // Read BAM header
    sam_hdr_t *header = sam_hdr_read(bam_fp);
    if (!header) {
        fprintf(stderr, "Error: Cannot read BAM header\n");
        sam_close(bam_fp);
        exit(1);
    }

    // Determine output format
    int is_damage_format = opts->enable_damage && 
                          (strcmp(opts->output_format, "dense_damage") == 0 || 
                           strcmp(opts->output_format, "sparse_damage") == 0);
    int is_sparse = strcmp(opts->output_format, "sparse") == 0 || 
                   strcmp(opts->output_format, "sparse_damage") == 0;
    
    // For dynamic genome discovery with sparse formats, use temporary file
    char temp_filename[1024];
    FILE *out_fp;
    int using_temp_file = 0;
    
    if (n_genomes == 0 && is_sparse) {
        // Create temporary file for data (no header)
        snprintf(temp_filename, sizeof(temp_filename), "%s.tmp", opts->output_file);
        out_fp = fopen(temp_filename, "w");
        using_temp_file = 1;
        if (!out_fp) {
            fprintf(stderr, "Error: Cannot create temporary file: %s\n", temp_filename);
            sam_hdr_destroy(header);
            sam_close(bam_fp);
            exit(1);
        }
    } else {
        // Open output file normally
        out_fp = fopen(opts->output_file, "w");
        if (!out_fp) {
            fprintf(stderr, "Error: Cannot create output file: %s\n", opts->output_file);
            sam_hdr_destroy(header);
            sam_close(bam_fp);
            exit(1);
        }
    }
    
    if (opts->verbose) {
        printf("Format detection: enable_damage=%d, format='%s', is_damage_format=%d, is_sparse=%d\n", 
               opts->enable_damage, opts->output_format, is_damage_format, is_sparse);
    }
    
    // Write header for all formats (unless using temp file)
    // Both dense and sparse formats require headers according to CEMfull specification
    if (!using_temp_file) {
        if (is_damage_format) {
            // Damage format header (both dense and sparse use same header)
            fprintf(out_fp, "read_id\ttotal_count");
            for (int i = 0; i < n_genomes; i++) {
                const char *name = get_output_genome_name_compressed(i, opts);
                if (is_sparse) {
                    // Sparse formats just list the genome names in header
                    fprintf(out_fp, "\t%s", name);
                } else {
                    // Dense formats use nd_, md_, mb_ prefixes
                    fprintf(out_fp, "\tnd_%s\tmd_%s\tmb_%s", name, name, name);
                }
            }
            fprintf(out_fp, "\n");
        } else {
            // Standard format header (both dense and sparse use same header)
            fprintf(out_fp, "read_id\ttotal_count");
            for (int i = 0; i < n_genomes; i++) {
                fprintf(out_fp, "\t%s", get_output_genome_name_compressed(i, opts));
            }
            fprintf(out_fp, "\n");
        }
    }

    // Data structures for tracking reads
    char current_read[MAX_NAME_LEN] = "";
    // Allocate arrays to handle dynamic genome discovery
    int initial_local_array_size = max_genomes_allocated;  // Start with current genome capacity
    int *mismatches = calloc(initial_local_array_size, sizeof(int));
    int *indel_counts = calloc(initial_local_array_size, sizeof(int));
    int *has_alignment = calloc(initial_local_array_size, sizeof(int));
    damage_stats_t *damage_stats_array = calloc(initial_local_array_size, sizeof(damage_stats_t));
    
    if (!mismatches || !indel_counts || !has_alignment || !damage_stats_array) {
        fprintf(stderr, "Error: Failed to allocate arrays for %d genomes\n", initial_local_array_size);
        exit(1);
    }
    
    // Variables to track current allocation size for these arrays
    int current_array_size = initial_local_array_size;
    int read_length = 0;
    int max_mismatches = 0;
    int reads_processed = 0;
    long long alignments_processed = 0;  // Track total alignments (use long long for large datasets)
    int unique_reads = 0;  // For compressed read ID generation

    bam1_t *read = bam_init1();

    // Report file processing start (unless silent)
    if (!opts->silent) {
        fprintf(stderr, "Processing BAM file: %s\n", opts->input_file);
    }
    
    if (opts->verbose) {
        fprintf(stderr, "Expected genomes: %d\n", n_genomes);
        for (int i = 0; i < n_genomes; i++) {
            fprintf(stderr, "  %d: %s\n", i+1, genome_names[i]);
        }
        fprintf(stderr, "\n");
    }

    while (sam_read1(bam_fp, header, read) >= 0) {
        
        // Check read limit
        if (opts->max_reads > 0 && reads_processed >= (int)opts->max_reads) {
            break;
        }

        // Skip unmapped reads
        if (read->core.flag & BAM_FUNMAP) continue;

        // Skip reads with low mapping quality
        if (read->core.qual < opts->min_mapq) continue;
        
        // Count alignments
        alignments_processed++;

        char *read_name = bam_get_qname(read);
        const char *ref_name = sam_hdr_tid2name(header, read->core.tid);
        
        // Check if this is a new read
        if (strcmp(read_name, current_read) != 0) {
            // Output previous read if it exists
            if (current_read[0] != '\0') {
                // Only output if we had valid alignments (read wasn't too short)
                if (read_length > 0) {
                    // Assign penalty to genomes without alignments based on penalty mode
                    if (opts->use_penalty_mode) {
                        // Penalty mode: use max+1 (capped at read_length)
                        for (int i = 0; i < n_genomes; i++) {
                            if (!has_alignment[i]) {
                                int penalty = max_mismatches + 1;
                                if (penalty > read_length) penalty = read_length;
                                mismatches[i] = penalty;
                                // Also set damage stats for consistency
                                damage_stats_array[i].nd = penalty;
                                damage_stats_array[i].md = penalty;
                                damage_stats_array[i].mb = penalty;
                            }
                        }
                    } else {
                        // Default: use -1 for missing alignments
                        for (int i = 0; i < n_genomes; i++) {
                            if (!has_alignment[i]) {
                                mismatches[i] = -1;
                                damage_stats_array[i].nd = -1;
                                damage_stats_array[i].md = -1;
                                damage_stats_array[i].mb = -1;
                            }
                        }
                    }
                    
                    // Output the read with proper format
                    if (is_sparse) {
                        // Sparse format - only output genomes with alignments
                        if (opts->compress_output) {
                            char compressed_id[32];
                            generate_compressed_read_id(unique_reads++, compressed_id);
                            fprintf(out_fp, "%s\t%d", compressed_id, read_length);
                        } else {
                            fprintf(out_fp, "%s\t%d", current_read, read_length);
                        }
                        
                        if (is_damage_format) {
                            // Sparse damage format
                            for (int i = 0; i < n_genomes; i++) {
                                if (has_alignment[i]) {
                                    const char *name = get_output_genome_name_compressed(i, opts);
                                    fprintf(out_fp, "\t%s\t%d\t%d\t%d", 
                                           name,
                                           damage_stats_array[i].nd,
                                           damage_stats_array[i].md,
                                           damage_stats_array[i].mb);
                                }
                            }
                        } else {
                            // Sparse standard format
                            for (int i = 0; i < n_genomes; i++) {
                                if (has_alignment[i]) {
                                    fprintf(out_fp, "\t%s\t%d", get_output_genome_name_compressed(i, opts), mismatches[i]);
                                }
                            }
                        }
                        fprintf(out_fp, "\n");
                    } else {
                        // Dense format - output all genomes
                        if (opts->compress_output) {
                            char compressed_id[32];
                            generate_compressed_read_id(unique_reads++, compressed_id);
                            fprintf(out_fp, "%s\t%d", compressed_id, read_length);
                        } else {
                            fprintf(out_fp, "%s\t%d", current_read, read_length);
                        }
                        
                        if (is_damage_format) {
                            // Dense damage format - values already set based on penalty mode
                            for (int i = 0; i < n_genomes; i++) {
                                fprintf(out_fp, "\t%d\t%d\t%d",
                                       damage_stats_array[i].nd,
                                       damage_stats_array[i].md,
                                       damage_stats_array[i].mb);
                            }
                        } else {
                            // Dense standard format - values already set based on penalty mode
                            for (int i = 0; i < n_genomes; i++) {
                                fprintf(out_fp, "\t%d", mismatches[i]);
                            }
                        }
                        fprintf(out_fp, "\n");
                    }
                    reads_processed++;
                }
            }
            
            // Initialize for new read
            strcpy(current_read, read_name);
            // Check if we need more space (also check against global capacity)
            if (n_genomes >= current_array_size || max_genomes_allocated > current_array_size) {
                int new_size = (max_genomes_allocated > current_array_size * 2) ? max_genomes_allocated : current_array_size * 2;
                int *new_mismatches = realloc(mismatches, new_size * sizeof(int));
                int *new_indel_counts = realloc(indel_counts, new_size * sizeof(int));
                int *new_has_alignment = realloc(has_alignment, new_size * sizeof(int));
                damage_stats_t *new_damage_stats_array = realloc(damage_stats_array, new_size * sizeof(damage_stats_t));
                
                if (!new_mismatches || !new_indel_counts || !new_has_alignment || !new_damage_stats_array) {
                    fprintf(stderr, "Error: Failed to reallocate arrays for %d genomes\n", new_size);
                    exit(1);
                }
                
                // Clear the new portion of the arrays
                int old_size = current_array_size;
                
                mismatches = new_mismatches;
                indel_counts = new_indel_counts;
                has_alignment = new_has_alignment;
                damage_stats_array = new_damage_stats_array;
                current_array_size = new_size;
                
                // Initialize the new portion
                memset(&mismatches[old_size], 0, (current_array_size - old_size) * sizeof(int));
                memset(&indel_counts[old_size], 0, (current_array_size - old_size) * sizeof(int));
                memset(&has_alignment[old_size], 0, (current_array_size - old_size) * sizeof(int));
                memset(&damage_stats_array[old_size], 0, (current_array_size - old_size) * sizeof(damage_stats_t));
            }
            
            memset(mismatches, 0, n_genomes * sizeof(int));
            memset(indel_counts, 0, n_genomes * sizeof(int));
            memset(has_alignment, 0, n_genomes * sizeof(int));
            memset(damage_stats_array, 0, n_genomes * sizeof(damage_stats_t));
            // Calculate trimmed read length - will be set properly when processing alignments
            int original_length = read->core.l_qseq;
            if (original_length >= MIN_READ_LENGTH) {
                if (opts->enable_damage) {
                    read_length = original_length;  // No trim when tracking damage
                } else if (opts->damage_sites > 0) {
                    read_length = original_length - 2 * opts->damage_sites;  // Trim
                } else {
                    read_length = original_length;  // No trim
                }
            } else {
                read_length = 0;  // Mark as invalid
            }
            max_mismatches = 0;
        }

        // Determine genome from reference name (default) or RG tag (optional)
        const char *genome_identifier = NULL;
        int genome_idx = -1;
        
        if (opts->use_rg_tag) {
            // Use RG tag (old behavior)
            uint8_t *rg_tag = bam_aux_get(read, "RG");
            if (!rg_tag) {
                if (opts->verbose) {
                    fprintf(stderr, "Warning: No RG tag for read %s\n", read_name);
                }
                continue;
            }
            genome_identifier = bam_aux2Z(rg_tag);
        } else {
            // Use reference name from alignment (column 3) - new default behavior
            const char *ref_name = sam_hdr_tid2name(header, read->core.tid);
            if (!ref_name || read->core.tid < 0) {
                if (opts->verbose) {
                    fprintf(stderr, "Warning: No reference for read %s (unmapped)\n", read_name);
                }
                continue;
            }
            genome_identifier = ref_name;
        }
        
        // Look up genome index
        genome_idx = opts->use_simple_mode ? 
            find_genome_index_simple(genome_identifier, opts->ignore_char) :
            find_genome_index_fast(genome_identifier, opts->ignore_char);
            
        if (genome_idx == -1) {
            // Try to add the genome dynamically
            const char *genome_name = extract_genome_name(genome_identifier, opts->ignore_char);
            genome_idx = add_genome_dynamically(genome_name);
            if (genome_idx == -1) {
                fprintf(stderr, "Warning: Failed to add genome '%s'\n", genome_name);
                continue;
            }
        }

        // Skip this read entirely if it's too short
        if (read_length <= 0) {
            continue;
        }

        
        // Calculate trimmed mismatch count with damage stats
        damage_stats_t current_damage_stats;
        int current_indel_count = 0;
        int nm_value = calculate_trimmed_mismatches(read, header, opts, 
                                                   opts->enable_damage ? &current_damage_stats : NULL,
                                                   &current_indel_count);
        
        // Skip this alignment if calculation failed
        if (nm_value < 0) {
            continue;
        }

        // Store the best alignment for this genome
        // Priority 1: Fewer indels
        // Priority 2: Fewer mismatches
        int is_better = 0;
        if (!has_alignment[genome_idx]) {
            is_better = 1;  // First alignment for this genome
        } else if (current_indel_count < indel_counts[genome_idx]) {
            is_better = 1;  // Fewer indels
        } else if (current_indel_count == indel_counts[genome_idx] && nm_value < mismatches[genome_idx]) {
            is_better = 1;  // Same indels but fewer mismatches
        }
        
        if (is_better) {
            mismatches[genome_idx] = nm_value;
            indel_counts[genome_idx] = current_indel_count;
            has_alignment[genome_idx] = 1;
            
            // Store damage stats if enabled
            if (opts->enable_damage) {
                // Debug: Check if we're creating an invalid state
                if (current_damage_stats.md > current_damage_stats.nd) {
                    if (opts->verbose) {
                        fprintf(stderr, "ERROR: About to store invalid damage stats for genome %d: nd=%d, md=%d, mb=%d\n",
                                genome_idx, current_damage_stats.nd, current_damage_stats.md, current_damage_stats.mb);
                    }
                }
                
                
                damage_stats_array[genome_idx] = current_damage_stats;
                
            }
            
            // Update max mismatches seen so far
            if (nm_value > max_mismatches) {
                max_mismatches = nm_value;
            }
        }

        // Progress reporting every 1,000,000 alignments (unless silent)
        if (!opts->silent && alignments_processed > 0 && (alignments_processed % 1000000 == 0)) {
            fprintf(stderr, "Processed %lld alignments (%d unique reads, %d genomes discovered)\n", 
                    alignments_processed, reads_processed, n_genomes);
        }
    }

    // Don't forget the last read
    if (current_read[0] != '\0' && read_length > 0) {
        
        for (int i = 0; i < n_genomes; i++) {
            if (!has_alignment[i]) {
                int penalty = max_mismatches + 1;
                if (penalty > read_length) penalty = read_length;
                mismatches[i] = penalty;
            }
        }
        
        // Output the final read with proper format
        if (is_sparse) {
            // Sparse format - only output genomes with alignments
            if (opts->compress_output) {
                char compressed_id[32];
                generate_compressed_read_id(unique_reads++, compressed_id);
                fprintf(out_fp, "%s\t%d", compressed_id, read_length);
            } else {
                fprintf(out_fp, "%s\t%d", current_read, read_length);
            }
            
            if (is_damage_format) {
                // Sparse damage format
                for (int i = 0; i < n_genomes; i++) {
                    if (has_alignment[i]) {
                        const char *name = get_output_genome_name_compressed(i, opts);
                        
                        
                        fprintf(out_fp, "\t%s\t%d\t%d\t%d", 
                               name,
                               damage_stats_array[i].nd,
                               damage_stats_array[i].md,
                               damage_stats_array[i].mb);
                    }
                }
            } else {
                // Sparse standard format
                for (int i = 0; i < n_genomes; i++) {
                    if (has_alignment[i]) {
                        fprintf(out_fp, "\t%s\t%d", get_output_genome_name_compressed(i, opts), mismatches[i]);
                    }
                }
            }
            fprintf(out_fp, "\n");
        } else {
            // Dense format - output all genomes
            if (opts->compress_output) {
                char compressed_id[32];
                generate_compressed_read_id(unique_reads++, compressed_id);
                fprintf(out_fp, "%s\t%d", compressed_id, read_length);
            } else {
                fprintf(out_fp, "%s\t%d", current_read, read_length);
            }
            
            if (is_damage_format) {
                // Dense damage format
                for (int i = 0; i < n_genomes; i++) {
                    if (has_alignment[i]) {
                        fprintf(out_fp, "\t%d\t%d\t%d", 
                               damage_stats_array[i].nd,
                               damage_stats_array[i].md,
                               damage_stats_array[i].mb);
                    } else {
                        // No alignment - use penalty values
                        fprintf(out_fp, "\t-1\t-1\t-1");
                    }
                }
            } else {
                // Dense standard format
                for (int i = 0; i < n_genomes; i++) {
                    fprintf(out_fp, "\t%d", mismatches[i]);
                }
            }
            fprintf(out_fp, "\n");
        }
        reads_processed++;
    }

    
    // Always print final statistics unless silent
    if (!opts->silent) {
        fprintf(stderr, "\nTotal alignments processed: %lld\n", alignments_processed);
        fprintf(stderr, "Total unique reads processed: %d\n", reads_processed);
        fprintf(stderr, "Total genomes discovered: %d\n", n_genomes);
        fprintf(stderr, "Analysis completed successfully!\n");
    }

    // Cleanup
    free(mismatches);
    free(indel_counts);
    free(has_alignment);
    free(damage_stats_array);
    bam_destroy1(read);
    sam_hdr_destroy(header);
    sam_close(bam_fp);
    fclose(out_fp);
    
    // If we used a temp file, now create the final file with correct header
    if (using_temp_file) {
        FILE *final_fp = fopen(opts->output_file, "w");
        if (!final_fp) {
            fprintf(stderr, "Error: Cannot create final output file: %s\n", opts->output_file);
            exit(1);
        }
        
        // Write correct header with all discovered genomes
        if (is_damage_format) {
            fprintf(final_fp, "read_id\ttotal_count");
            for (int i = 0; i < n_genomes; i++) {
                const char *name = get_output_genome_name_compressed(i, opts);
                fprintf(final_fp, "\t%s", name);
            }
            fprintf(final_fp, "\n");
        } else {
            fprintf(final_fp, "read_id\ttotal_count");
            for (int i = 0; i < n_genomes; i++) {
                fprintf(final_fp, "\t%s", get_output_genome_name(i, opts->short_names));
            }
            fprintf(final_fp, "\n");
        }
        
        // Copy data from temp file
        FILE *temp_fp = fopen(temp_filename, "r");
        if (!temp_fp) {
            fprintf(stderr, "Error: Cannot open temporary file: %s\n", temp_filename);
            fclose(final_fp);
            exit(1);
        }
        
        char buffer[8192];
        size_t bytes_read;
        size_t total_bytes = 0;
        int chunk_count = 0;
        while ((bytes_read = fread(buffer, 1, sizeof(buffer), temp_fp)) > 0) {
            fwrite(buffer, 1, bytes_read, final_fp);
            total_bytes += bytes_read;
            chunk_count++;
        }
        
        fclose(temp_fp);
        fclose(final_fp);
        
        // Remove temp file
        remove(temp_filename);
        
        if (!opts->silent) {
            fprintf(stderr, "Created final output with correct header\n");
        }
    }
}
#endif

void process_text_file(options_t *opts) {
    FILE *input_fp = fopen(opts->input_file, "r");
    if (!input_fp) {
        fprintf(stderr, "Error: Cannot open input file: %s\n", opts->input_file);
        exit(1);
    }

    FILE *output_fp = fopen(opts->output_file, "w");
    if (!output_fp) {
        fprintf(stderr, "Error: Cannot create output file: %s\n", opts->output_file);
        fclose(input_fp);
        exit(1);
    }

    if (opts->verbose) {
        fprintf(stderr, "Processing text file: %s\n", opts->input_file);
    }

    char line[8192];
    int line_num = 0;
    int reads_processed = 0;

    while (fgets(line, sizeof(line), input_fp)) {
        line_num++;
        
        // Remove newline
        line[strcspn(line, "\n")] = 0;
        
        // Skip empty lines
        if (line[0] == '\0') continue;
        
        // Process and copy all lines (including header)
        fprintf(output_fp, "%s\n", line);
        
        if (line_num > 1) {  // Count data lines only
            reads_processed++;
        }
        
        if (opts->verbose && line_num == 1) {
            fprintf(stderr, "Header: %s\n", line);
        }
        
        if (opts->verbose && reads_processed > 0 && (reads_processed % 10000 == 0)) {
            fprintf(stderr, "Processed %d reads\n", reads_processed);
        }
    }

    if (opts->verbose) {
        fprintf(stderr, "Total reads processed: %d\n", reads_processed);
    }

    fclose(input_fp);
    fclose(output_fp);
}

// Simple (original) algorithm functions
int find_genome_index_simple(const char *ref_name, char ignore_char) {
    const char *genome_name = extract_genome_name(ref_name, ignore_char);
    
    pthread_mutex_lock(&genome_mutex);
    
    for (int i = 0; i < n_genomes; i++) {
        if (strcmp(genome_name, genome_names[i]) == 0) {
            pthread_mutex_unlock(&genome_mutex);
            return i;
        }
    }
    
    pthread_mutex_unlock(&genome_mutex);
    return -1;  // Not found
}

int parse_md_tag_simple(const char *md_string, int *mismatch_positions, int max_mismatches) {
    if (!md_string) return 0;
    
    int pos = 0;
    int num_mismatches = 0;
    const char *ptr = md_string;
    
    while (*ptr && num_mismatches < max_mismatches) {
        if (isdigit(*ptr)) {
            // Parse number of matching bases
            int matches = 0;
            while (isdigit(*ptr)) {
                matches = matches * 10 + (*ptr - '0');
                ptr++;
            }
            pos += matches;
        } else if (*ptr == '^') {
            // Deletion - skip the ^ and the deleted bases
            ptr++; // skip ^
            while (*ptr && isalpha(*ptr)) {
                ptr++; // skip deleted bases
            }
            // Deletions don't advance position in read
        } else if (isalpha(*ptr)) {
            // Mismatch base - store position for future position-specific analysis
            mismatch_positions[num_mismatches++] = pos;
            pos++;
            ptr++;
        } else {
            ptr++; // skip unknown characters
        }
    }
    
    return num_mismatches;
}

// Optimization functions
void init_optimizations(void) {
    // Initialize character type lookup table for MD parsing
    for (int i = 0; i < 256; i++) {
        if (i >= '0' && i <= '9') {
            char_types[i] = CHAR_DIGIT;
        } else if ((i >= 'A' && i <= 'Z') || (i >= 'a' && i <= 'z')) {
            char_types[i] = CHAR_ALPHA;
        } else {
            char_types[i] = CHAR_OTHER;
        }
    }
    
    // Initialize hash table
    for (int i = 0; i < HASH_TABLE_SIZE; i++) {
        genome_hash_table[i] = NULL;
    }
}

// djb2 hash function - fast and effective for strings
unsigned int hash_string(const char *str) {
    unsigned int hash = 5381;
    int c;
    while ((c = *str++)) {
        hash = ((hash << 5) + hash) + c; // hash * 33 + c
    }
    return hash & (HASH_TABLE_SIZE - 1); // Fast modulo for power of 2
}

void add_genome_to_hash(const char *name, int index) {
    unsigned int hash = hash_string(name);
    genome_hash_entry_t *entry = malloc(sizeof(genome_hash_entry_t));
    if (!entry) {
        fprintf(stderr, "Error: Memory allocation failed for hash table entry\n");
        exit(1);
    }
    
    strncpy(entry->name, name, MAX_NAME_LEN - 1);
    entry->name[MAX_NAME_LEN - 1] = '\0';
    entry->index = index;
    entry->next = genome_hash_table[hash];
    genome_hash_table[hash] = entry;
}

// Function to dynamically add a genome during processing
int add_genome_dynamically(const char *genome_name) {
    pthread_mutex_lock(&genome_mutex);
    
    // Check if genome already exists (another thread might have added it)
    for (int i = 0; i < n_genomes; i++) {
        if (strcmp(genome_names[i], genome_name) == 0) {
            pthread_mutex_unlock(&genome_mutex);
            return i;
        }
    }
    
    // Check if we need to grow the arrays
    if (n_genomes >= max_genomes_allocated) {
        grow_genome_arrays();
    }
    
    // Add the genome
    strcpy(genome_names[n_genomes], genome_name);
    add_genome_to_hash(genome_names[n_genomes], n_genomes);
    
    // Generate compressed ID for this genome
    generate_genome_compressed_id(n_genomes, genome_name);
    
    int genome_index = n_genomes;
    n_genomes++;
    
    pthread_mutex_unlock(&genome_mutex);
    return genome_index;  // Return the index of the newly added genome
}

int find_genome_index_fast(const char *ref_name, char ignore_char) {
    const char *genome_name = extract_genome_name(ref_name, ignore_char);
    unsigned int hash = hash_string(genome_name);
    
    pthread_mutex_lock(&genome_mutex);
    
    genome_hash_entry_t *entry = genome_hash_table[hash];
    while (entry) {
        if (strcmp(genome_name, entry->name) == 0) {
            int index = entry->index;
            pthread_mutex_unlock(&genome_mutex);
            return index;
        }
        entry = entry->next;
    }
    
    pthread_mutex_unlock(&genome_mutex);
    return -1; // Not found
}

int is_directory(const char *path) {
    struct stat statbuf;
    if (stat(path, &statbuf) != 0) {
        return 0;
    }
    return S_ISDIR(statbuf.st_mode);
}

int find_bam_files(const char *directory) {
    DIR *dir;
    struct dirent *entry;
    char filepath[MAX_NAME_LEN];
    
    dir = opendir(directory);
    if (dir == NULL) {
        fprintf(stderr, "Error: Cannot open directory: %s\n", directory);
        return -1;
    }
    
    n_bam_files = 0;
    while ((entry = readdir(dir)) != NULL && n_bam_files < MAX_BAM_FILES) {
        // Check if file ends with .bam
        int len = strlen(entry->d_name);
        if (len > 4 && strcmp(entry->d_name + len - 4, ".bam") == 0) {
            snprintf(filepath, sizeof(filepath), "%s/%s", directory, entry->d_name);
            strncpy(bam_files[n_bam_files], filepath, MAX_NAME_LEN - 1);
            bam_files[n_bam_files][MAX_NAME_LEN - 1] = '\0';
            n_bam_files++;
        }
    }
    
    closedir(dir);
    return n_bam_files;
}

#ifdef WITH_HTSLIB
void process_multiple_bam_files(options_t *opts) {
    // Check if we can use optimized mode
    if (opts->no_redundancy) {
        // Use optimized version that doesn't track reads across files
        process_multiple_bam_files_optimized(opts);
        return;
    }
    
    // Display warning for standard mode with many genomes
    if (!opts->silent) {
        fprintf(stderr, "\n");
        fprintf(stderr, "================================================================================\n");
        fprintf(stderr, "WARNING: RUNNING IN MULTIPLE FILES MODE ASSUMING THE SAME READ MAY APPEAR\n");
        fprintf(stderr, "IN MULTIPLE BAM FILES. THIS IS SLOW IF THERE ARE MANY REFERENCE GENOMES.\n");
        fprintf(stderr, "TO SPEED UP THE CALCULATION CHOOSE --no-redundancy IF IT CAN BE ASSUMED\n");
        fprintf(stderr, "THAT A READ ONLY WILL APPEAR IN ONE BAM FILE.\n");
        fprintf(stderr, "================================================================================\n");
        fprintf(stderr, "\n");
    }
    
    // Dynamic array for tracking reads across multiple BAM files
    int read_capacity = 100000;  // Initial capacity
    read_data_t *reads = calloc(read_capacity, sizeof(read_data_t));
    int n_reads = 0;
    long long total_alignments = 0;  // Use long long to handle > 2 billion alignments
    
    // Optimization: Track last read to take advantage of consecutive alignments
    char last_read_name[MAX_NAME_LEN] = "";
    read_data_t *last_read_entry = NULL;
    
    if (opts->verbose) {
        printf("Processing %d BAM files...\n", n_bam_files);
        for (int i = 0; i < n_bam_files; i++) {
            printf("  %s\n", bam_files[i]);
        }
        printf("\n");
    }
    
    // Process each BAM file
    for (int file_idx = 0; file_idx < n_bam_files; file_idx++) {
        char *bam_file = bam_files[file_idx];
        
        // Report file processing start (unless silent)
        if (!opts->silent) {
            printf("Processing BAM file: %s\n", bam_file);
        }
        
        // Open BAM file
        samFile *bam_fp = sam_open(bam_file, "r");
        if (!bam_fp) {
            fprintf(stderr, "Error: Cannot open BAM file: %s\n", bam_file);
            continue;
        }
        
        // Apply I/O optimizations for each BAM file
        if (opts->num_threads > 0) {
            // Sanity check for thread count
            if (opts->num_threads > 256) {
                fprintf(stderr, "Warning: Invalid thread count %d, using 4 threads instead\n", opts->num_threads);
                opts->num_threads = 4;
            }
            if (hts_set_threads(bam_fp, opts->num_threads) < 0) {
                fprintf(stderr, "Warning: Failed to set %d threads\n", opts->num_threads);
            } else if (opts->verbose) {
                fprintf(stderr, "Info: Using %d decompression threads for %s\n", opts->num_threads, bam_file);
            }
        }
        
        // Set larger cache size (128MB)
        size_t cache_size = 128 * 1024 * 1024;
        hts_set_cache_size(bam_fp, cache_size);
        if (opts->verbose) {
            fprintf(stderr, "Info: Set cache size to %zu MB for %s\n", cache_size / (1024*1024), bam_file);
        }
        
        // Read BAM header
        sam_hdr_t *header = sam_hdr_read(bam_fp);
        if (!header) {
            fprintf(stderr, "Error: Cannot read BAM header from %s\n", bam_file);
            sam_close(bam_fp);
            continue;
        }
        
        bam1_t *read = bam_init1();
        long long alignments_in_file = 0;  // Use long long for large files
        
        while (sam_read1(bam_fp, header, read) >= 0) {
            alignments_in_file++;
            total_alignments++;
            
            // Progress reporting every 1,000,000 alignments (for large dataset tracking)
            if (total_alignments > 0 && (total_alignments % 1000000 == 0)) {
                struct timeval current_time;
                gettimeofday(&current_time, NULL);
                double elapsed = (current_time.tv_sec - program_start_time.tv_sec) + (current_time.tv_usec - program_start_time.tv_usec) / 1000000.0;
                double rate = total_alignments / elapsed;
                fprintf(stderr, "Progress: %lldM alignments (%d reads) - %.1fs elapsed, %.0f align/sec\n", 
                        total_alignments/1000000, n_reads, elapsed, rate);
            }
            
            
            
            // Check read limit
            if (opts->max_reads > 0 && total_alignments >= opts->max_reads) {
                break;
            }
            
            // Skip unmapped reads
            if (read->core.flag & BAM_FUNMAP) continue;
            
            // Skip reads with low mapping quality
            if (read->core.qual < opts->min_mapq) continue;
            
            char *read_name = bam_get_qname(read);
            
            // Get read group (RG) tag to determine genome
            uint8_t *rg_tag = bam_aux_get(read, "RG");
            if (!rg_tag) {
                if (opts->verbose) {
                    fprintf(stderr, "Warning: No RG tag for read %s\n", read_name);
                }
                continue;
            }
            
            char *rg_value = bam_aux2Z(rg_tag);
            int genome_idx = opts->use_simple_mode ? 
            find_genome_index_simple(rg_value, opts->ignore_char) :
            find_genome_index_fast(rg_value, opts->ignore_char);
            if (genome_idx == -1) {
                // Try to add the genome dynamically
                const char *genome_name = extract_genome_name(rg_value, opts->ignore_char);
                genome_idx = add_genome_dynamically(genome_name);
                if (genome_idx == -1) {
                    fprintf(stderr, "Warning: Failed to add genome '%s'\n", genome_name);
                    continue;
                }
            }
            
            // Calculate trimmed mismatch count with damage stats
            damage_stats_t current_damage_stats;
            int current_indel_count = 0;
            int nm_value = calculate_trimmed_mismatches(read, header, opts, 
                                                       opts->enable_damage ? &current_damage_stats : NULL,
                                                       &current_indel_count);
            if (nm_value < 0) {
                continue;  // Skip reads that are too short
            }
            
            
            
            // Find or create read entry
            read_data_t *read_entry = NULL;
            
            // Optimization: First check if it's the same as the last read (consecutive alignments)
            if (last_read_entry != NULL && strcmp(last_read_name, read_name) == 0) {
                read_entry = last_read_entry;
            } else {
                // Different read - use suffix tree for fast lookup
                read_entry = suffix_tree_find(read_name);
                
                // Update last read tracking for next iteration
                strcpy(last_read_name, read_name);
                last_read_entry = read_entry;
            }
            
            if (read_entry == NULL) {
                // New read - check if we need to resize the array
                if (n_reads >= read_capacity) {
                    // Double the capacity
                    int new_capacity = read_capacity * 2;
                    if (!opts->silent) {
                        fprintf(stderr, "Info: Expanding read array from %d to %d entries\n", 
                                read_capacity, new_capacity);
                    }
                    
                    read_data_t *new_reads = realloc(reads, new_capacity * sizeof(read_data_t));
                    if (!new_reads) {
                        fprintf(stderr, "Error: Failed to allocate memory for %d reads\n", new_capacity);
                        break;
                    }
                    
                    // Clear the new memory
                    memset(&new_reads[read_capacity], 0, 
                           (new_capacity - read_capacity) * sizeof(read_data_t));
                    
                    reads = new_reads;
                    read_capacity = new_capacity;
                    
                    // Update suffix tree pointers after realloc
                    // Note: This requires rebuilding the suffix tree since pointers changed
                    suffix_tree_free(suffix_tree_root);
                    suffix_tree_root = NULL;
                    for (int i = 0; i < n_reads; i++) {
                        suffix_tree_insert(&reads[i]);
                    }
                }
                
                read_entry = &reads[n_reads++];
                strcpy(read_entry->read_id, read_name);
                if (opts->enable_damage) {
                    read_entry->read_length = read->core.l_qseq;  // No trim when tracking damage
                } else if (opts->damage_sites > 0) {
                    read_entry->read_length = read->core.l_qseq - 2 * opts->damage_sites;  // Trim
                } else {
                    read_entry->read_length = read->core.l_qseq;  // No trim
                }
                read_entry->max_mismatches = 0;
                
                // Allocate arrays for this read based on current genome count
                read_entry->mismatches = calloc(max_genomes_allocated, sizeof(int));
                read_entry->indel_count = calloc(max_genomes_allocated, sizeof(int));
                read_entry->has_alignment = calloc(max_genomes_allocated, sizeof(int));
                read_entry->damage_stats = calloc(max_genomes_allocated, sizeof(damage_stats_t));
                
                if (!read_entry->mismatches || !read_entry->indel_count || 
                    !read_entry->has_alignment || !read_entry->damage_stats) {
                    fprintf(stderr, "Error: Failed to allocate memory for read data arrays\n");
                    break;
                }
                
                // Initialize all mismatches to "no alignment"
                for (int i = 0; i < n_genomes; i++) {
                    read_entry->mismatches[i] = 0;
                    read_entry->has_alignment[i] = 0;
                    read_entry->damage_stats[i].nd = 0;
                    read_entry->damage_stats[i].md = 0;
                    read_entry->damage_stats[i].mb = 0;
                }
                
                // Insert the new read into the suffix tree for fast future lookups
                suffix_tree_insert(read_entry);
                
                // Update last read tracking for the new read
                strcpy(last_read_name, read_name);
                last_read_entry = read_entry;
                
                // Additional progress reporting for unique read milestones
                // (Commented out since we now report based on alignments)
                // if (!opts->silent && n_reads > 0 && (n_reads % 10000000 == 0)) {
                //     fprintf(stderr, "Processed %d unique reads\n", n_reads);
                // }
            }
            
            // Update best alignment for this genome
            // Priority 1: Fewer indels
            // Priority 2: Fewer mismatches
            int is_better = 0;
            if (!read_entry->has_alignment[genome_idx]) {
                is_better = 1;  // First alignment for this genome
            } else if (current_indel_count < read_entry->indel_count[genome_idx]) {
                is_better = 1;  // Fewer indels
            } else if (current_indel_count == read_entry->indel_count[genome_idx] && 
                       nm_value < read_entry->mismatches[genome_idx]) {
                is_better = 1;  // Same indels but fewer mismatches
            }
            
            if (is_better) {
                read_entry->mismatches[genome_idx] = nm_value;
                read_entry->indel_count[genome_idx] = current_indel_count;
                read_entry->has_alignment[genome_idx] = 1;
                
                // Store damage stats if enabled
                if (opts->enable_damage) {
                    read_entry->damage_stats[genome_idx] = current_damage_stats;
                }
                
                if (nm_value > read_entry->max_mismatches) {
                    read_entry->max_mismatches = nm_value;
                }
            }
        }
        
        if (opts->verbose) {
            printf("Processed %lld alignments\n", alignments_in_file);
        }
        
        bam_destroy1(read);
        sam_hdr_destroy(header);
        sam_close(bam_fp);
        
        if (opts->max_reads > 0 && total_alignments >= opts->max_reads) {
            break;
        }
    }
    
    // Write output
    FILE *out_fp = fopen(opts->output_file, "w");
    if (!out_fp) {
        fprintf(stderr, "Error: Cannot create output file: %s\n", opts->output_file);
        free(reads);
        suffix_tree_free(suffix_tree_root);
        suffix_tree_root = NULL;
        return;
    }
    
    // Determine output format
    int is_damage_format = opts->enable_damage && 
                          (strcmp(opts->output_format, "dense_damage") == 0 || 
                           strcmp(opts->output_format, "sparse_damage") == 0);
    int is_sparse = strcmp(opts->output_format, "sparse") == 0 || 
                   strcmp(opts->output_format, "sparse_damage") == 0;
    
    if (opts->verbose) {
        printf("Format detection: enable_damage=%d, format='%s', is_damage_format=%d, is_sparse=%d\n", 
               opts->enable_damage, opts->output_format, is_damage_format, is_sparse);
    }
    
    // Write header for all formats
    // Both dense and sparse formats require headers according to CEMfull specification
    {
        if (is_damage_format) {
            // Damage format header (both dense and sparse use same header)
            fprintf(out_fp, "read_id\ttotal_count");
            for (int i = 0; i < n_genomes; i++) {
                const char *name = get_output_genome_name_compressed(i, opts);
                if (is_sparse) {
                    // Sparse formats just list the genome names in header
                    fprintf(out_fp, "\t%s", name);
                } else {
                    // Dense formats use nd_, md_, mb_ prefixes
                    fprintf(out_fp, "\tnd_%s\tmd_%s\tmb_%s", name, name, name);
                }
            }
            fprintf(out_fp, "\n");
        } else {
            // Standard format header (both dense and sparse use same header)
            fprintf(out_fp, "read_id\ttotal_count");
            for (int i = 0; i < n_genomes; i++) {
                fprintf(out_fp, "\t%s", get_output_genome_name_compressed(i, opts));
            }
            fprintf(out_fp, "\n");
        }
    }
    
    // Write read data
    for (int i = 0; i < n_reads; i++) {
        read_data_t *read_entry = &reads[i];
        
        // Assign -1 to genomes without alignments
        for (int j = 0; j < n_genomes; j++) {
            if (!read_entry->has_alignment[j]) {
                read_entry->mismatches[j] = -1;
                read_entry->damage_stats[j].nd = -1;
                read_entry->damage_stats[j].md = -1;
                read_entry->damage_stats[j].mb = -1;
            }
        }
        
        if (is_sparse) {
            // Sparse format - only output genomes with alignments
            fprintf(out_fp, "%s\t%d", read_entry->read_id, read_entry->read_length);
            
            if (is_damage_format) {
                // Sparse damage format
                for (int j = 0; j < n_genomes; j++) {
                    if (read_entry->has_alignment[j]) {
                        const char *name = get_output_genome_name(j, opts->short_names);
                        fprintf(out_fp, "\t%s\t%d\t%d\t%d", 
                               name,
                               read_entry->damage_stats[j].nd,
                               read_entry->damage_stats[j].md,
                               read_entry->damage_stats[j].mb);
                    }
                }
            } else {
                // Sparse standard format
                for (int j = 0; j < n_genomes; j++) {
                    if (read_entry->has_alignment[j]) {
                        const char *name = get_output_genome_name(j, opts->short_names);
                        fprintf(out_fp, "\t%s\t%d", name, read_entry->mismatches[j]);
                    }
                }
            }
            fprintf(out_fp, "\n");
        } else {
            // Dense format - output all genomes
            fprintf(out_fp, "%s\t%d", read_entry->read_id, read_entry->read_length);
            
            if (is_damage_format) {
                // Dense damage format
                for (int j = 0; j < n_genomes; j++) {
                    if (read_entry->has_alignment[j]) {
                        fprintf(out_fp, "\t%d\t%d\t%d", 
                               read_entry->damage_stats[j].nd,
                               read_entry->damage_stats[j].md,
                               read_entry->damage_stats[j].mb);
                    } else {
                        // No alignment - use penalty values
                        fprintf(out_fp, "\t-1\t-1\t-1");
                    }
                }
            } else {
                // Dense standard format
                for (int j = 0; j < n_genomes; j++) {
                    fprintf(out_fp, "\t%d", read_entry->mismatches[j]);
                }
            }
            fprintf(out_fp, "\n");
        }
    }
    
    if (opts->verbose) {
        printf("Total alignments processed: %lld\n", total_alignments);
        printf("Total unique reads processed: %d\n", n_reads);
    }
    
    fclose(out_fp);
    
    // Free all allocated arrays for each read
    for (int i = 0; i < n_reads; i++) {
        if (reads[i].mismatches) free(reads[i].mismatches);
        if (reads[i].indel_count) free(reads[i].indel_count);
        if (reads[i].has_alignment) free(reads[i].has_alignment);
        if (reads[i].damage_stats) free(reads[i].damage_stats);
    }
    free(reads);
    
    // Free the suffix tree
    suffix_tree_free(suffix_tree_root);
    suffix_tree_root = NULL;
}

// Optimized version of process_multiple_bam_files for --no-redundancy mode
void process_multiple_bam_files_optimized(options_t *opts) {
    if (!opts->silent) {
        fprintf(stderr, "Info: Running in optimized no-redundancy mode\n");
        fprintf(stderr, "      (assuming reads don't appear in multiple BAM files)\n");
        fprintf(stderr, "      (n_bam_files = %d)\n\n", n_bam_files);
        if (opts->consolidate_by_taxid) {
            fprintf(stderr, "Info: Using taxid consolidation mode\n");
            fprintf(stderr, "      Multiple genomes mapping to same taxid will be consolidated\n\n");
        }
    }
    
    // Check if we have BAM files to process
    if (n_bam_files == 0) {
        fprintf(stderr, "Error: No BAM files to process (n_bam_files = 0)\n");
        return;
    }
    
    // Initialize sparse memory pool if not already initialized
    if (!sparse_global_pool) {
        sparse_global_pool = create_sparse_memory_pool();
        if (!sparse_global_pool) {
            fprintf(stderr, "Error: Failed to create sparse memory pool\n");
            return;
        }
    }
    
    // Create temp file for all alignments
    char temp_filename[MAX_NAME_LEN];
    snprintf(temp_filename, sizeof(temp_filename), "%s.tmp", opts->output_file);
    FILE *temp_fp = fopen(temp_filename, "w");
    if (!temp_fp) {
        fprintf(stderr, "Error: Cannot create temp file: %s\n", temp_filename);
        return;
    }
    
    long long total_alignments = 0;
    long long total_reads = 0;
    
    // Process each BAM file independently
    for (int file_idx = 0; file_idx < n_bam_files; file_idx++) {
        char *bam_file = bam_files[file_idx];
        
        if (!opts->silent) {
            printf("Processing BAM file %d/%d: %s\n", file_idx + 1, n_bam_files, bam_file);
        }
        
        // Open BAM file
        samFile *bam_fp = sam_open(bam_file, "r");
        if (!bam_fp) {
            fprintf(stderr, "Error: Cannot open BAM file: %s\n", bam_file);
            continue;
        }
        
        // Apply I/O optimizations
        if (opts->num_threads > 0) {
            if (hts_set_threads(bam_fp, opts->num_threads) < 0) {
                fprintf(stderr, "Warning: Failed to set %d threads\n", opts->num_threads);
            }
        }
        
        // Set larger cache size
        size_t cache_size = 128 * 1024 * 1024;
        hts_set_cache_size(bam_fp, cache_size);
        
        // Read BAM header
        sam_hdr_t *header = sam_hdr_read(bam_fp);
        if (!header) {
            fprintf(stderr, "Error: Cannot read BAM header from %s\n", bam_file);
            sam_close(bam_fp);
            continue;
        }
        
        bam1_t *read = bam_init1();
        genome_alignment_t *current_alignments = NULL;
        char current_read_name[MAX_NAME_LEN] = "";
        int current_read_length = 0;
        long long file_alignments = 0;
        
        while (sam_read1(bam_fp, header, read) >= 0) {
            file_alignments++;
            total_alignments++;
            
            // Progress reporting
            if (!opts->silent && total_alignments > 0 && (total_alignments % 1000000 == 0)) {
                fprintf(stderr, "Processed %lld alignments, %lld unique reads\n", 
                        total_alignments, total_reads);
            }
            
            // Check read limit
            if (opts->max_reads > 0 && total_alignments >= opts->max_reads) {
                break;
            }
            
            // Skip unmapped reads
            if (read->core.flag & BAM_FUNMAP) continue;
            
            // Skip reads with low mapping quality
            if (read->core.qual < opts->min_mapq) continue;
            
            char *read_name = bam_get_qname(read);
            
            // Check if this is a new read
            if (strcmp(current_read_name, read_name) != 0) {
                // Output previous read if exists
                if (current_alignments) {
                    // Write read data to temp file
                    fprintf(temp_fp, "%s\t%d", current_read_name, current_read_length);
                    genome_alignment_t *align = current_alignments;
                    while (align) {
                        if (opts->enable_damage) {
                            fprintf(temp_fp, "\t%s\t%d\t%d\t%d", 
                                    align->genome_name,
                                    align->damage.nd, align->damage.md, align->damage.mb);
                        } else {
                            fprintf(temp_fp, "\t%s\t%d", align->genome_name, align->mismatches);
                        }
                        align = align->next;
                    }
                    fprintf(temp_fp, "\n");
                    total_reads++;
                    
                    // Free alignments - just set to NULL, memory pool handles cleanup
                    current_alignments = NULL;
                }
                
                // Start new read
                strcpy(current_read_name, read_name);
                if (opts->enable_damage) {
                    current_read_length = read->core.l_qseq;  // No trim when tracking damage
                } else if (opts->damage_sites > 0) {
                    current_read_length = read->core.l_qseq - 2 * opts->damage_sites;  // Trim
                    if (current_read_length < 0) current_read_length = 0;
                } else {
                    current_read_length = read->core.l_qseq;  // No trim
                }
            }
            
            // Get genome from RG tag or reference
            const char *genome_name = NULL;
            char taxid_name[MAX_NAME_LEN];

            if (opts->consolidate_by_taxid) {
                // Get reference name for taxid lookup
                const char *ref_name = sam_hdr_tid2name(header, read->core.tid);
                if (!ref_name) continue;

                // Extract base accession
                const char *extracted_name = extract_genome_name(ref_name, opts->ignore_char);
                if (!extracted_name) continue;

                char base_accession[MAX_NAME_LEN];
                strncpy(base_accession, extracted_name, MAX_NAME_LEN - 1);
                base_accession[MAX_NAME_LEN - 1] = '\0';
                char *dot = strchr(base_accession, '.');
                if (dot) *dot = '\0';

                // Look up taxid
                int taxid = get_taxid_for_accession(base_accession);
                if (taxid > 0) {
                    snprintf(taxid_name, MAX_NAME_LEN, "T%d", taxid);
                    genome_name = taxid_name;
                } else {
                    genome_name = extracted_name;
                }
            } else {
                // Original logic
                if (opts->use_rg_tag) {
                    uint8_t *rg_tag = bam_aux_get(read, "RG");
                    if (rg_tag) {
                        char *rg_value = bam_aux2Z(rg_tag);
                        genome_name = extract_genome_name(rg_value, opts->ignore_char);
                    }
                } else {
                    const char *ref_name = sam_hdr_tid2name(header, read->core.tid);
                    genome_name = extract_genome_name(ref_name, opts->ignore_char);
                }
            }

            if (!genome_name) continue;
            
            // Calculate mismatches
            damage_stats_t damage_stats;
            int indel_count = 0;
            int nm_value = calculate_trimmed_mismatches(read, header, opts, 
                                                       opts->enable_damage ? &damage_stats : NULL,
                                                       &indel_count);
            if (nm_value < 0) continue;
            
            // Add to current read's alignments (no genome lookup!)
            genome_alignment_t *align = find_or_add_genome_alignment_optimized(
                &current_alignments, genome_name, sparse_global_pool);
            
            
            // Update if better alignment
            if (align->mismatches == -1 || nm_value < align->mismatches) {
                align->mismatches = nm_value;
                align->indel_count = indel_count;
                if (opts->enable_damage) {
                    align->damage = damage_stats;
                }
            }
        }
        
        // Output last read if exists
        if (current_alignments) {
            // Write read data to temp file
            fprintf(temp_fp, "%s\t%d", current_read_name, current_read_length);
            genome_alignment_t *align = current_alignments;
            while (align) {
                if (opts->enable_damage) {
                    fprintf(temp_fp, "\t%s\t%d\t%d\t%d", 
                            align->genome_name,
                            align->damage.nd, align->damage.md, align->damage.mb);
                } else {
                    fprintf(temp_fp, "\t%s\t%d", align->genome_name, align->mismatches);
                }
                align = align->next;
            }
            fprintf(temp_fp, "\n");
            total_reads++;
        }
        
        bam_destroy1(read);
        sam_hdr_destroy(header);
        sam_close(bam_fp);
        
        if (opts->verbose) {
            printf("  Processed %lld alignments from %s\n", file_alignments, bam_file);
        }
        
        if (opts->max_reads > 0 && total_alignments >= opts->max_reads) {
            break;
        }
    }
    
    fclose(temp_fp);
    
    if (!opts->silent) {
        printf("Total alignments processed: %lld\n", total_alignments);
        printf("Total unique reads processed: %lld\n", total_reads);
    }
    
    // Process temp file to deduplicate genomes and write final output
    int is_damage_format = opts->enable_damage && 
                          (strcmp(opts->output_format, "sparse_damage") == 0 || 
                           strcmp(opts->output_format, "dense_damage") == 0);
    process_temp_file_optimized(temp_filename, opts->output_file, is_damage_format);
    
    // Clean up temp file
    remove(temp_filename);
}

// Unified dispatcher for multiple BAM files - routes to appropriate function based on options
void process_multiple_bam_files_unified(options_t *opts) {
    if (opts->with_higher_taxa) {
        // NEW: Multiple BAM taxonomy mode (sparse only, taxid-based)
        process_multiple_bam_files_with_taxonomy(opts);
    } else if (opts->no_redundancy) {
        // EXISTING: Optimized sparse mode (no cross-file read tracking needed)
        process_multiple_bam_files_optimized(opts);
    } else {
        // EXISTING: General mode (handles both dense AND sparse, with cross-file tracking)
        process_multiple_bam_files(opts);
    }
}

// Process multiple BAM files with taxonomy support using taxid-based processing
void process_multiple_bam_files_with_taxonomy(options_t *opts) {
#ifdef WITH_HTSLIB
    // Clear messaging about mode
    if (!opts->silent) {
        fprintf(stderr, "\n=== MULTIPLE BAM TAXID-BASED TAXONOMY MODE ===\n");
        fprintf(stderr, "Processing %d BAM files with taxid consolidation\n", n_bam_files);
        fprintf(stderr, "Cross-file best alignment per taxid automatically handled\n");
        fprintf(stderr, "===============================================\n");
    }
    
    // REUSE: Same initialization as single BAM taxonomy mode
    if (!read_mismatch_table) {
        read_mismatch_table = (ReadMismatchEntry **)calloc(read_mismatch_table_size, 
                                                           sizeof(ReadMismatchEntry *));
    }
    
    long long total_alignments = 0;
    long long total_reads = 0;

    // ONLY NEW PART: Loop through multiple files
    for (int file_idx = 0; file_idx < n_bam_files; file_idx++) {
        char *bam_file = bam_files[file_idx];
        
        if (!opts->silent) {
            fprintf(stderr, "Processing file %d/%d: %s\n", file_idx + 1, n_bam_files, bam_file);
        }
        
        // Standard BAM file opening (same as single BAM)
        samFile *bam_fp = sam_open(bam_file, "r");
        if (!bam_fp) {
            fprintf(stderr, "Warning: Cannot open BAM file: %s\n", bam_file);
            continue;
        }
        
        sam_hdr_t *header = sam_hdr_read(bam_fp);
        if (!header) {
            fprintf(stderr, "Warning: Cannot read header from: %s\n", bam_file);
            sam_close(bam_fp);
            continue;
        }
        
        bam1_t *read = bam_init1();
        
        // COPY from existing single BAM taxonomy processing (lines ~5300-5570)
        while (sam_read1(bam_fp, header, read) >= 0) {
            total_alignments++;
            
            // Progress reporting every 1,000,000 alignments with memory stats
            if (!opts->silent && total_alignments % 1000000 == 0) {
                fprintf(stderr, "Processed %lld alignments, %lld reads\n", total_alignments, total_reads);
                
                // Calculate memory usage statistics
                int total_read_entries = 0;
                size_t total_positions_memory = 0;
                size_t total_aligned_arrays_memory = 0;
                
                for (int hash_idx = 0; hash_idx < read_mismatch_table_size; hash_idx++) {
                    ReadMismatchEntry *entry = read_mismatch_table[hash_idx];
                    while (entry) {
                        total_read_entries++;
                        ReadMismatches *data = entry->data;
                        
                        // Count positions memory
                        if (data->positions) {
                            for (int pos = 0; pos < data->read_length; pos++) {
                                if (data->positions[pos].has_mismatch) {
                                    int bytes = (data->positions[pos].n_genomes + 7) / 8;
                                    total_positions_memory += bytes * 2; // has_mismatch + is_damage
                                }
                            }
                        }
                        
                        // Count aligned arrays memory
                        total_aligned_arrays_memory += data->capacity * sizeof(int) * 3; // 3 arrays
                        
                        entry = entry->next;
                    }
                }
                
                fprintf(stderr, "  Memory stats: %d read entries, %.2f MB positions, %.2f MB arrays\n",
                        total_read_entries, total_positions_memory / 1048576.0, total_aligned_arrays_memory / 1048576.0);
            }
            
            // Check limits
            if (opts->max_reads > 0 && total_alignments > opts->max_reads) {
                break;
            }
            
            // Skip unmapped reads
            if (read->core.flag & BAM_FUNMAP) continue;

            // Skip low quality
            if (read->core.qual < opts->min_mapq) continue;

            // CHECK FOR INSERTIONS: Skip alignments with insertions in taxonomy mode
            // The MD tag parser doesn't handle insertions correctly
            uint32_t *cigar = bam_get_cigar(read);
            int has_insertion = 0;

            for (int i = 0; i < read->core.n_cigar; i++) {
                if (bam_cigar_op(cigar[i]) == BAM_CINS) {
                    has_insertion = 1;
                    break;
                }
            }

            if (has_insertion) {
                // Skip this alignment - taxonomy mode doesn't handle insertions
                continue;  // Skip to next alignment
            }

            char *read_name = bam_get_qname(read);
            
            // Get reference name for taxid lookup (taxonomy mode only uses reference names)
            const char *ref_name = sam_hdr_tid2name(header, read->core.tid);
            if (!ref_name) continue;
            
            // REUSE existing taxid discovery
            int genome_idx = find_or_add_taxid(ref_name, opts->ignore_char);
            if (genome_idx < 0) continue;
            
            // REUSE existing ReadMismatches logic - automatically handles cross-file reads!
            ReadMismatches *current_read_data = get_or_create_read_mismatches(read_name);


            // Set read length if not already set (needed for MD tag parsing)
            if (current_read_data->read_length == 0) {
                current_read_data->read_length = read->core.l_qseq;
                // Set trimmed_length based on damage mode
                if (opts->enable_damage) {
                    current_read_data->trimmed_length = read->core.l_qseq;  // No trim when tracking damage
                } else if (opts->damage_sites > 0) {
                    current_read_data->trimmed_length = read->core.l_qseq - 2 * opts->damage_sites;  // Trim
                    if (current_read_data->trimmed_length < 0) current_read_data->trimmed_length = 0;
                } else {
                    current_read_data->trimmed_length = read->core.l_qseq;  // No trim
                }
            }
            
            // Track that this genome has an alignment for this read
            int already_tracked = 0;
            int genome_array_idx = -1;
            for (int i = 0; i < current_read_data->n_aligned; i++) {
                if (current_read_data->aligned_genomes[i] == genome_idx) {
                    already_tracked = 1;
                    genome_array_idx = i;
                    break;
                }
            }

            // Get NM tag to check total mismatches for this alignment
            uint8_t *nm_tag = bam_aux_get(read, "NM");
            int nm_value = 0;
            if (nm_tag) {
                nm_value = bam_aux2i(nm_tag);
            }

            // Handle new or existing genome alignment
            if (!already_tracked) {
                // New genome - add it
                if (current_read_data->n_aligned >= current_read_data->capacity) {
                    current_read_data->capacity = current_read_data->capacity ? current_read_data->capacity * 2 : 100;
                    current_read_data->aligned_genomes = realloc(current_read_data->aligned_genomes, 
                                                                current_read_data->capacity * sizeof(int));
                    current_read_data->genome_nd_values = realloc(current_read_data->genome_nd_values,
                                                                 current_read_data->capacity * sizeof(int));
                    current_read_data->genome_nm_values = realloc(current_read_data->genome_nm_values,
                                                                 current_read_data->capacity * sizeof(int));
                }
                genome_array_idx = current_read_data->n_aligned;
                current_read_data->aligned_genomes[genome_array_idx] = genome_idx;
                current_read_data->genome_nd_values[genome_array_idx] = 0;  // Initialize to 0
                current_read_data->genome_nm_values[genome_array_idx] = nm_value;
                current_read_data->n_aligned++;

                // FIX: Reallocate ALL existing position bit vectors to match new n_aligned size
                int new_bytes_needed = (current_read_data->n_aligned + 7) / 8;
                for (int pos = 0; pos < current_read_data->read_length; pos++) {
                    if (current_read_data->positions && current_read_data->positions[pos].has_mismatch) {
                        int old_bytes = (current_read_data->positions[pos].n_genomes + 7) / 8;

                        if (new_bytes_needed > old_bytes) {
                            // Reallocate this position's bit vectors
                            current_read_data->positions[pos].has_mismatch =
                                (uint8_t *)realloc(current_read_data->positions[pos].has_mismatch, new_bytes_needed);
                            current_read_data->positions[pos].is_damage =
                                (uint8_t *)realloc(current_read_data->positions[pos].is_damage, new_bytes_needed);

                            // Check for realloc failure
                            if (!current_read_data->positions[pos].has_mismatch || !current_read_data->positions[pos].is_damage) {
                                fprintf(stderr, "FATAL: Failed to reallocate position %d bit vectors for n_aligned=%d\n",
                                        pos, current_read_data->n_aligned);
                                abort();
                            }

                            // Clear new bytes
                            memset(current_read_data->positions[pos].has_mismatch + old_bytes, 0, new_bytes_needed - old_bytes);
                            memset(current_read_data->positions[pos].is_damage + old_bytes, 0, new_bytes_needed - old_bytes);

                            // Update position's n_genomes to reflect new capacity
                            current_read_data->positions[pos].n_genomes = current_read_data->n_aligned;
                        }
                    }
                }
            } else {
                // Genome already tracked - check if this alignment is better
                if (nm_value >= current_read_data->genome_nm_values[genome_array_idx]) {
                    // This alignment is not better (same or worse), skip ALL processing
                    continue;
                }
                // This alignment is better - clear previous mismatch data and update
                if (current_read_data->positions) {
                    // Find the aligned array position for this genome_idx
                    int aligned_idx = -1;
                    
                    // DIRECT TEST: Log the lookup bounds and results
                    static int lookup_count = 0;
                    lookup_count++;
                    if (lookup_count <= 5) {
                        fprintf(stderr, "LOOKUP #%d: searching for genome_idx=%d in n_aligned=%d\n", 
                                lookup_count, genome_idx, current_read_data->n_aligned);
                    }
                    
                    for (int i = 0; i < current_read_data->n_aligned; i++) {
                        if (current_read_data->aligned_genomes[i] == genome_idx) {
                            aligned_idx = i;
                            if (lookup_count <= 5) {
                                fprintf(stderr, "  Found at aligned_idx=%d\n", aligned_idx);
                            }
                            break;
                        }
                    }
                    
                    if (aligned_idx >= 0) {
                        
                        for (int pos = 0; pos < current_read_data->read_length; pos++) {
                            if (current_read_data->positions[pos].has_mismatch) {
                                int byte_idx = aligned_idx / 8;
                                int bit_idx = aligned_idx % 8;
                                
                                
                                // Ensure bit vectors are properly sized before access
                                PositionMismatch *pos_data = &current_read_data->positions[pos];
                                if (pos_data->n_genomes < current_read_data->n_aligned) {
                                    int new_bytes_needed = (current_read_data->n_aligned + 7) / 8;
                                    int old_bytes = (pos_data->n_genomes + 7) / 8;

                                    if (new_bytes_needed > old_bytes) {
                                        pos_data->has_mismatch = (uint8_t *)realloc(pos_data->has_mismatch, new_bytes_needed);
                                        pos_data->is_damage = (uint8_t *)realloc(pos_data->is_damage, new_bytes_needed);

                                        // Clear new bytes
                                        memset(pos_data->has_mismatch + old_bytes, 0, new_bytes_needed - old_bytes);
                                        memset(pos_data->is_damage + old_bytes, 0, new_bytes_needed - old_bytes);
                                    }
                                    pos_data->n_genomes = current_read_data->n_aligned;
                                }

                                current_read_data->positions[pos].has_mismatch[byte_idx] &= ~(1 << bit_idx);
                                current_read_data->positions[pos].is_damage[byte_idx] &= ~(1 << bit_idx);
                            }
                        }
                    }
                }
                // Update the NM value to the better one
                current_read_data->genome_nm_values[genome_array_idx] = nm_value;
            }

            // Process MD tag and store mismatch positions
            uint8_t *md_tag = bam_aux_get(read, "MD");
            if (!md_tag) continue;

            char *md_string = bam_aux2Z(md_tag);

            // Parse MD tag to get mismatch positions and calculate nd
            int ref_pos = 0;
            int read_pos = 0;
            char *p = md_string;
            int nd_count = 0;  // Count C/G in damage-prone regions
            uint8_t *seq = bam_get_seq(read);

            while (*p) {
                if (isdigit(*p)) {
                    // Match run - here we can count C/G in matching positions
                    int match_len = 0;
                    while (isdigit(*p)) {
                        match_len = match_len * 10 + (*p - '0');
                        p++;
                    }

                    // Count C/G in damage-prone regions during matches
                    // In matches, read base == ref base
                    for (int i = 0; i < match_len && read_pos + i < current_read_data->read_length; i++) {
                        int curr_pos = read_pos + i;
                        char read_base = seq_nt16_str[bam_seqi(seq, curr_pos)];

                        // Check if in damage-prone region and has C or G
                        if (opts->asymmetric_damage) {
                            // Asymmetric: C in first s, G in last s
                            if (curr_pos < opts->damage_sites && (read_base == 'C' || read_base == 'c')) {
                                nd_count++;
                            } else if (curr_pos >= current_read_data->read_length - opts->damage_sites &&
                                       (read_base == 'G' || read_base == 'g')) {
                                nd_count++;
                            }
                        } else {
                            // Symmetric: C or G in first s and last s
                            if (curr_pos < opts->damage_sites ||
                                curr_pos >= current_read_data->read_length - opts->damage_sites) {
                                if (read_base == 'C' || read_base == 'c' || read_base == 'G' || read_base == 'g') {
                                    nd_count++;
                                }
                            }
                        }
                    }

                    ref_pos += match_len;
                    read_pos += match_len;
                } else if (*p == '^') {
                    // Deletion in read
                    p++;
                    while (*p && isalpha(*p)) {
                        ref_pos++;
                        p++;
                    }
                } else if (isalpha(*p)) {
                    // Mismatch
                    char ref_base = *p;
                    char read_base = seq_nt16_str[bam_seqi(seq, read_pos)];

                    // Determine if this is a damage type mismatch
                    int is_damage = 0;
                    if (opts->enable_damage) {
                        if (opts->asymmetric_damage) {
                            // Asymmetric: C->T in first s, G->A in last s
                            if (read_pos < opts->damage_sites &&
                                ((ref_base == 'C' && (read_base == 'T' || read_base == 't')) ||
                                 (ref_base == 'c' && (read_base == 'T' || read_base == 't')))) {
                                is_damage = 1;
                            } else if (read_pos >= current_read_data->read_length - opts->damage_sites &&
                                       ((ref_base == 'G' && (read_base == 'A' || read_base == 'a')) ||
                                        (ref_base == 'g' && (read_base == 'A' || read_base == 'a')))) {
                                is_damage = 1;
                            }
                        } else {
                            // Symmetric: C->T or G->A in first s and last s
                            if (read_pos < opts->damage_sites ||
                                read_pos >= current_read_data->read_length - opts->damage_sites) {
                                if ((ref_base == 'C' && (read_base == 'T' || read_base == 't')) ||
                                    (ref_base == 'c' && (read_base == 'T' || read_base == 't')) ||
                                    (ref_base == 'G' && (read_base == 'A' || read_base == 'a')) ||
                                    (ref_base == 'g' && (read_base == 'A' || read_base == 'a'))) {
                                    is_damage = 1;
                                }
                            }
                        }
                    }

                    // Count damage-susceptible sites and store mismatch position
                    if (opts->enable_damage) {
                        if (opts->asymmetric_damage) {
                            if (read_pos < opts->damage_sites && (ref_base == 'C' || ref_base == 'c')) {
                                nd_count++;
                            } else if (read_pos >= current_read_data->read_length - opts->damage_sites &&
                                       (ref_base == 'G' || ref_base == 'g')) {
                                nd_count++;
                            }
                        } else {
                            if (read_pos < opts->damage_sites ||
                                read_pos >= current_read_data->read_length - opts->damage_sites) {
                                if (ref_base == 'C' || ref_base == 'c' || ref_base == 'G' || ref_base == 'g') {
                                    nd_count++;
                                }
                            }
                        }
                    }

                    // Store mismatch position without any trimming adjustment
                    // We're keeping the full read, not trimming
                    store_read_mismatch_position(current_read_data, read_pos,
                                                genome_idx, 1, is_damage);
                    ref_pos++;
                    read_pos++;
                    p++;
                } else {
                    p++;
                }
            }


            // Update nd value for this genome
            if (current_read_data && genome_array_idx >= 0) {
                current_read_data->genome_nd_values[genome_array_idx] = nd_count;
            }
        }
        
        // Standard cleanup per file
        bam_destroy1(read);
        sam_hdr_destroy(header);
        sam_close(bam_fp);

        if (opts->max_reads > 0 && total_alignments >= opts->max_reads) {
            break;
        }
    }

    // Finalize all accumulated reads
    for (int hash_idx = 0; hash_idx < read_mismatch_table_size; hash_idx++) {
        ReadMismatchEntry *entry = read_mismatch_table[hash_idx];
        while (entry) {
            total_reads += finalize_read_mismatches(entry->data);
            entry = entry->next;
        }
    }

    if (!opts->silent) {
        fprintf(stderr, "Processed %lld total alignments from %lld unique reads\n",
                total_alignments, total_reads);
    }

    // REUSE: Exact same taxonomy tree population as single BAM
    for (int i = 0; i < n_genomes; i++) {
        int taxid = get_taxid_from_index(i);
        
        if (taxid > 0 && taxid <= taxonomy_tree->max_taxid && taxonomy_tree->nodes[taxid]) {
            TaxNode *node = taxonomy_tree->nodes[taxid];
            
            while (node) {
                if (!node->leaf_genomes) {
                    node->leaves_capacity = 10;
                    node->leaf_genomes = (int *)calloc(node->leaves_capacity, sizeof(int));
                }
                if (node->n_leaves >= node->leaves_capacity) {
                    node->leaves_capacity *= 2;
                    node->leaf_genomes = (int *)realloc(node->leaf_genomes,
                        node->leaves_capacity * sizeof(int));
                }
                if (node->leaf_genomes) {
                    node->leaf_genomes[node->n_leaves++] = i;
                }
                
                node->is_active = 1;
                
                if (node->parent_taxid == node->taxid) break;
                node = (node->parent_taxid <= taxonomy_tree->max_taxid) ? 
                       taxonomy_tree->nodes[node->parent_taxid] : NULL;
            }
        }
    }
    
    // REUSE: Same taxonomy output as single BAM
    struct timeval linear_start, linear_end;
    gettimeofday(&linear_start, NULL);
    fprintf(stderr, "Starting linear algorithm taxonomic calculation...\n");
    
    output_with_taxonomy_linear(NULL, opts);
    
    gettimeofday(&linear_end, NULL);
    linear_algorithm_time = (linear_end.tv_sec - linear_start.tv_sec) + (linear_end.tv_usec - linear_start.tv_usec) / 1000000.0;

    // Print comprehensive timing analysis for multiple BAM mode
    struct timeval program_end;
    gettimeofday(&program_end, NULL);
    double total_program_time = (program_end.tv_sec - program_start_time.tv_sec) + (program_end.tv_usec - program_start_time.tv_usec) / 1000000.0;

    fprintf(stderr, "\n=== PERFORMANCE ANALYSIS ===\n");
    fprintf(stderr, "Total runtime: %.2fs\n", total_program_time);
    fprintf(stderr, "BAM processing (total): %.2fs (%.1f%%)\n", bam_processing_time, bam_processing_time/total_program_time*100);
    fprintf(stderr, "  - BAM file I/O: %.2fs (%.1f%%)\n", bam_processing_time - linear_algorithm_time, (bam_processing_time - linear_algorithm_time)/total_program_time*100);
    fprintf(stderr, "  - Linear taxonomy algorithm: %.2fs (%.1f%%)\n", linear_algorithm_time, linear_algorithm_time/total_program_time*100);
    fprintf(stderr, "Other phases (startup/cleanup): %.2fs (%.1f%%)\n",
            total_program_time - bam_processing_time,
            (total_program_time - bam_processing_time)/total_program_time*100);
    fprintf(stderr, "============================\n");

    // REUSE: Same cleanup as single BAM
    for (int i = 0; i < read_mismatch_table_size; i++) {
        ReadMismatchEntry *entry = read_mismatch_table[i];
        while (entry) {
            ReadMismatchEntry *next = entry->next;
            free_read_mismatches(entry->data);
            free(entry);
            entry = next;
        }
    }
    free(read_mismatch_table);
    read_mismatch_table = NULL;

#else
    fprintf(stderr, "Error: BAM mode not available - program compiled without HTSlib\n");
#endif
}

// Auto-discover genome names from BAM file(s)
// Process BAM file with taxonomy support - stores position-specific mismatches
void process_bam_file_with_taxonomy(options_t *opts) {
#ifdef WITH_HTSLIB

    if (!opts->silent) {
        fprintf(stderr, "\n=== TAXID-BASED TAXONOMY MODE ===\n");
        fprintf(stderr, "Processing mode: TAXID-BASED (--with-higher-taxa)\n");
        fprintf(stderr, "Note: This mode processes unique taxonomic IDs instead of individual genomes.\n");
        fprintf(stderr, "      Output columns represent taxids (T12345) with best alignment per taxid.\n");
        fprintf(stderr, "      Use non-taxonomy mode for genome-level analysis.\n");
        fprintf(stderr, "=====================================\n");
        if (opts->max_reads > 0) {
            fprintf(stderr, "Alignment limit: %lld alignments\n", opts->max_reads);
        }
        fprintf(stderr, "Discovering unique taxids during BAM processing...\n");
    }
    
    // Initialize read mismatch storage
    if (!read_mismatch_table) {
        read_mismatch_table = (ReadMismatchEntry **)calloc(read_mismatch_table_size, 
                                                           sizeof(ReadMismatchEntry *));
    }
    
    // Open BAM file
    samFile *bam_fp = sam_open(opts->input_file, "r");
    if (!bam_fp) {
        fprintf(stderr, "Error: Cannot open BAM file: %s\n", opts->input_file);
        return;
    }
    
    // Read header
    sam_hdr_t *header = sam_hdr_read(bam_fp);
    if (!header) {
        fprintf(stderr, "Error: Cannot read BAM header\n");
        sam_close(bam_fp);
        return;
    }
    
    bam1_t *read = bam_init1();
    long long total_alignments = 0;
    long long total_reads = 0;
    char current_read_name[MAX_NAME_LEN] = "";
    ReadMismatches *current_read_data = NULL;
    
    // Process alignments
    while (sam_read1(bam_fp, header, read) >= 0) {
        total_alignments++;
        
        // Progress reporting every 1,000,000 alignments
        if (!opts->silent && total_alignments % 1000000 == 0) {
            fprintf(stderr, "Processed %lld alignments, %lld reads\n", total_alignments, total_reads);
        }
        
        
        // Check limits
        if (opts->max_reads > 0 && total_alignments > opts->max_reads) {
            break;
        }
        
        // Skip unmapped reads
        if (read->core.flag & BAM_FUNMAP) continue;

        // Skip low quality
        if (read->core.qual < opts->min_mapq) continue;

        // CHECK FOR INSERTIONS: Skip alignments with insertions in taxonomy mode
        // The MD tag parser doesn't handle insertions correctly
        uint32_t *cigar = bam_get_cigar(read);
        int has_insertion = 0;

        for (int i = 0; i < read->core.n_cigar; i++) {
            if (bam_cigar_op(cigar[i]) == BAM_CINS) {
                has_insertion = 1;
                break;
            }
        }

        if (has_insertion) {
            // Skip this alignment - taxonomy mode doesn't handle insertions
            continue;  // Skip to next alignment
        }

        char *read_name = bam_get_qname(read);
        
        // Check if new read
        if (strcmp(current_read_name, read_name) != 0) {
            // Finalize previous read
            if (current_read_data) {
                total_reads += finalize_read_mismatches(current_read_data);
            }
            
            // Start new read
            strcpy(current_read_name, read_name);
            
            // Get or create read data
            current_read_data = get_or_create_read_mismatches(read_name);

            current_read_data->read_length = read->core.l_qseq;
            // Set trimmed_length based on damage mode
            if (opts->enable_damage) {
                current_read_data->trimmed_length = read->core.l_qseq;  // No trim when tracking damage
            } else if (opts->damage_sites > 0) {
                current_read_data->trimmed_length = read->core.l_qseq - 2 * opts->damage_sites;  // Trim
                if (current_read_data->trimmed_length < 0) current_read_data->trimmed_length = 0;
            } else {
                current_read_data->trimmed_length = read->core.l_qseq;  // No trim
            }
        }
        
        // Get genome name
        const char *genome_name = NULL;
        if (opts->use_rg_tag) {
            uint8_t *rg_tag = bam_aux_get(read, "RG");
            if (rg_tag) {
                char *rg_value = bam_aux2Z(rg_tag);
                genome_name = extract_genome_name(rg_value, opts->ignore_char);
            }
        } else {
            genome_name = sam_hdr_tid2name(header, read->core.tid);
            if (genome_name && opts->ignore_char) {
                genome_name = extract_genome_name(genome_name, opts->ignore_char);
            }
        }
        
        if (!genome_name) continue;
        
        // Find taxid index (replaces genome lookup)
        int genome_idx = find_or_add_taxid(sam_hdr_tid2name(header, read->core.tid), opts->ignore_char);
        if (genome_idx < 0) continue;
        
        // Track that this genome has an alignment for this read
        int already_tracked = 0;
        for (int i = 0; i < current_read_data->n_aligned; i++) {
            if (current_read_data->aligned_genomes[i] == genome_idx) {
                already_tracked = 1;
                break;
            }
        }
        // Get NM tag to check total mismatches for this alignment
        uint8_t *nm_tag = bam_aux_get(read, "NM");
        int nm_value = 0;
        if (nm_tag) {
            nm_value = bam_aux2i(nm_tag);
        }


        
        // Find the index for this genome in our aligned_genomes array
        int genome_array_idx = -1;
        for (int i = 0; i < current_read_data->n_aligned; i++) {
            if (current_read_data->aligned_genomes[i] == genome_idx) {
                genome_array_idx = i;
                break;
            }
        }
        
        if (!already_tracked) {
            // New genome - add it
            if (current_read_data->n_aligned >= current_read_data->capacity) {
                current_read_data->capacity = current_read_data->capacity ? current_read_data->capacity * 2 : 100;
                current_read_data->aligned_genomes = realloc(current_read_data->aligned_genomes, 
                                                            current_read_data->capacity * sizeof(int));
                current_read_data->genome_nd_values = realloc(current_read_data->genome_nd_values,
                                                             current_read_data->capacity * sizeof(int));
                current_read_data->genome_nm_values = realloc(current_read_data->genome_nm_values,
                                                             current_read_data->capacity * sizeof(int));
            }
            genome_array_idx = current_read_data->n_aligned;
            current_read_data->aligned_genomes[genome_array_idx] = genome_idx;
            // Initialize nd to 0 - will be calculated when processing mismatches
            current_read_data->genome_nd_values[genome_array_idx] = 0;
            // Store the NM value for this alignment
            current_read_data->genome_nm_values[genome_array_idx] = nm_value;
            current_read_data->n_aligned++;

            // FIX: Reallocate ALL existing position bit vectors to match new n_aligned size
            int new_bytes_needed = (current_read_data->n_aligned + 7) / 8;
            for (int pos = 0; pos < current_read_data->read_length; pos++) {
                if (current_read_data->positions && current_read_data->positions[pos].has_mismatch) {
                    int old_bytes = (current_read_data->positions[pos].n_genomes + 7) / 8;

                    if (new_bytes_needed > old_bytes) {
                        // Reallocate this position's bit vectors
                        current_read_data->positions[pos].has_mismatch =
                            (uint8_t *)realloc(current_read_data->positions[pos].has_mismatch, new_bytes_needed);
                        current_read_data->positions[pos].is_damage =
                            (uint8_t *)realloc(current_read_data->positions[pos].is_damage, new_bytes_needed);

                        // Check for realloc failure
                        if (!current_read_data->positions[pos].has_mismatch || !current_read_data->positions[pos].is_damage) {
                            fprintf(stderr, "FATAL: Failed to reallocate position %d bit vectors for n_aligned=%d\n",
                                    pos, current_read_data->n_aligned);
                            abort();
                        }

                        // Clear new bytes
                        memset(current_read_data->positions[pos].has_mismatch + old_bytes, 0, new_bytes_needed - old_bytes);
                        memset(current_read_data->positions[pos].is_damage + old_bytes, 0, new_bytes_needed - old_bytes);

                        // Update position's n_genomes to reflect new capacity
                        current_read_data->positions[pos].n_genomes = current_read_data->n_aligned;
                    }
                }
            }
        } else {
            // Genome already tracked - check if this alignment is better
            if (nm_value >= current_read_data->genome_nm_values[genome_array_idx]) {
                // This alignment is not better (same or worse), skip MD processing
                continue;
            }
            // This alignment is better - clear previous mismatch data
            // Clear all mismatch bits for this genome
            if (current_read_data->positions) {
                // Find the aligned array position for this genome_idx
                int aligned_idx = -1;
                for (int i = 0; i < current_read_data->n_aligned; i++) {
                    if (current_read_data->aligned_genomes[i] == genome_idx) {
                        aligned_idx = i;
                        break;
                    }
                }
                
                if (aligned_idx >= 0) {
                    for (int pos = 0; pos < current_read_data->read_length; pos++) {
                        if (current_read_data->positions[pos].has_mismatch) {
                            int byte_idx = aligned_idx / 8;
                            int bit_idx = aligned_idx % 8;
                            
                            
                            // Ensure bit vectors are properly sized before access
                            PositionMismatch *pos_data = &current_read_data->positions[pos];
                            if (pos_data->n_genomes < current_read_data->n_aligned) {
                                int new_bytes_needed = (current_read_data->n_aligned + 7) / 8;
                                int old_bytes = (pos_data->n_genomes + 7) / 8;

                                if (new_bytes_needed > old_bytes) {
                                    pos_data->has_mismatch = (uint8_t *)realloc(pos_data->has_mismatch, new_bytes_needed);
                                    pos_data->is_damage = (uint8_t *)realloc(pos_data->is_damage, new_bytes_needed);

                                    // Clear new bytes
                                    memset(pos_data->has_mismatch + old_bytes, 0, new_bytes_needed - old_bytes);
                                    memset(pos_data->is_damage + old_bytes, 0, new_bytes_needed - old_bytes);
                                }
                                pos_data->n_genomes = current_read_data->n_aligned;
                            }

                            current_read_data->positions[pos].has_mismatch[byte_idx] &= ~(1 << bit_idx);
                            current_read_data->positions[pos].is_damage[byte_idx] &= ~(1 << bit_idx);
                        }
                    }
                }
            }
            // Update the NM value to the better one
            current_read_data->genome_nm_values[genome_array_idx] = nm_value;
        }
        
        // Get MD tag for mismatch positions
        uint8_t *md_tag = bam_aux_get(read, "MD");
        if (!md_tag) continue;
        
        char *md_string = bam_aux2Z(md_tag);
        
        // Parse MD tag to get mismatch positions and calculate nd
        int ref_pos = 0;
        int read_pos = 0;
        char *p = md_string;
        int nd_count = 0;  // Count C/G in damage-prone regions
        uint8_t *seq = bam_get_seq(read);
        
        while (*p) {
            if (isdigit(*p)) {
                // Match run - here we can count C/G in matching positions
                int match_len = 0;
                while (isdigit(*p)) {
                    match_len = match_len * 10 + (*p - '0');
                    p++;
                }
                
                // Count C/G in damage-prone regions during matches
                // In matches, read base == ref base
                for (int i = 0; i < match_len && read_pos + i < current_read_data->read_length; i++) {
                    int curr_pos = read_pos + i;
                    char read_base = seq_nt16_str[bam_seqi(seq, curr_pos)];
                    
                    // Check if in damage-prone region and has C or G
                    if (opts->asymmetric_damage) {
                        // Asymmetric: C in first s, G in last s
                        if (curr_pos < opts->damage_sites && (read_base == 'C' || read_base == 'c')) {
                            nd_count++;
                        } else if (curr_pos >= current_read_data->read_length - opts->damage_sites && 
                                   (read_base == 'G' || read_base == 'g')) {
                            nd_count++;
                        }
                    } else {
                        // Symmetric: C or G in first s and last s
                        if (curr_pos < opts->damage_sites || 
                            curr_pos >= current_read_data->read_length - opts->damage_sites) {
                            if (read_base == 'C' || read_base == 'c' || read_base == 'G' || read_base == 'g') {
                                nd_count++;
                            }
                        }
                    }
                }
                
                ref_pos += match_len;
                read_pos += match_len;
            } else if (*p == '^') {
                // Deletion in read
                p++;
                while (*p && isalpha(*p)) {
                    ref_pos++;
                    p++;
                }
            } else if (isalpha(*p)) {
                // Mismatch - *p is the reference base
                char ref_base = *p;
                
                // Count nd if ref base is C or G in damage region
                if (opts->asymmetric_damage) {
                    if (read_pos < opts->damage_sites && (ref_base == 'C' || ref_base == 'c')) {
                        nd_count++;
                    } else if (read_pos >= current_read_data->read_length - opts->damage_sites && 
                               (ref_base == 'G' || ref_base == 'g')) {
                        nd_count++;
                    }
                } else {
                    if (read_pos < opts->damage_sites || 
                        read_pos >= current_read_data->read_length - opts->damage_sites) {
                        if (ref_base == 'C' || ref_base == 'c' || ref_base == 'G' || ref_base == 'g') {
                            nd_count++;
                        }
                    }
                }
                
                // Store ALL mismatches, properly categorized
                // Determine if this is a damage mismatch
                int is_damage = 0;
                char read_base = seq_nt16_str[bam_seqi(seq, read_pos)];
                
                if (opts->enable_damage) {
                    // Check if in damage-prone region and is damage-type mismatch
                    if (opts->asymmetric_damage) {
                        if (read_pos < opts->damage_sites) {
                            // First s sites - check for C->T
                            if ((ref_base == 'C' || ref_base == 'c') && (read_base == 'T' || read_base == 't')) {
                                is_damage = 1;
                            }
                        } else if (read_pos >= current_read_data->read_length - opts->damage_sites) {
                            // Last s sites - check for G->A
                            if ((ref_base == 'G' || ref_base == 'g') && (read_base == 'A' || read_base == 'a')) {
                                is_damage = 1;
                            }
                        }
                    } else {
                        // Symmetric mode
                        if (read_pos < opts->damage_sites || 
                            read_pos >= current_read_data->read_length - opts->damage_sites) {
                            // In damage-prone regions - check for C->T or G->A
                            if (((ref_base == 'C' || ref_base == 'c') && (read_base == 'T' || read_base == 't')) ||
                                ((ref_base == 'G' || ref_base == 'g') && (read_base == 'A' || read_base == 'a'))) {
                                is_damage = 1;
                            }
                        }
                    }
                }
                
                // Store mismatch position without any trimming adjustment
                // We're keeping the full read, not trimming
                store_read_mismatch_position(current_read_data, read_pos, 
                                            genome_idx, 1, is_damage);
                ref_pos++;
                read_pos++;
                p++;
            }
        }
        
        // Store the calculated nd value for this genome
        current_read_data->genome_nd_values[genome_array_idx] = nd_count;

        // Progress reporting
        if (!opts->silent && total_alignments % 1000000 == 0) {
            fprintf(stderr, "Processed %lld alignments, %lld reads\n",
                    total_alignments, total_reads);
        }
    }
    
    // Finalize last read
    if (current_read_data) {
        total_reads += finalize_read_mismatches(current_read_data);
    }
    
    bam_destroy1(read);
    sam_hdr_destroy(header);
    sam_close(bam_fp);
    
    if (!opts->silent) {
        fprintf(stderr, "Processed %lld total alignments from %lld unique reads\n", 
                total_alignments, total_reads);
    }
    
    // Populate taxonomy tree with discovered unique taxids (replaces map_genomes_to_taxids functionality)
    for (int i = 0; i < n_genomes; i++) {  // n_genomes now = unique taxids
        int taxid = get_taxid_from_index(i);
        
        if (taxid > 0 && taxid <= taxonomy_tree->max_taxid && taxonomy_tree->nodes[taxid]) {
            TaxNode *node = taxonomy_tree->nodes[taxid];
            
            // Add this taxid index to node and all ancestors
            while (node) {
                if (!node->leaf_genomes) {
                    node->leaves_capacity = 10;
                    node->leaf_genomes = (int *)calloc(node->leaves_capacity, sizeof(int));
                }
                if (node->n_leaves >= node->leaves_capacity) {
                    node->leaves_capacity *= 2;
                    node->leaf_genomes = (int *)realloc(node->leaf_genomes,
                        node->leaves_capacity * sizeof(int));
                }
                if (node->leaf_genomes) {
                    node->leaf_genomes[node->n_leaves++] = i;  // i is unique taxid index
                }
                
                node->is_active = 1;
                
                // Move to parent
                if (node->parent_taxid == node->taxid) break;
                node = (node->parent_taxid <= taxonomy_tree->max_taxid) ? 
                       taxonomy_tree->nodes[node->parent_taxid] : NULL;
            }
        }
    }
    
    // Time the linear algorithm phase
    struct timeval linear_start, linear_end;
    gettimeofday(&linear_start, NULL);
    fprintf(stderr, "Starting linear algorithm taxonomic calculation...\n");
    
    // Now calculate taxonomic mismatches and output using LINEAR ALGORITHM
    output_with_taxonomy_linear(NULL, opts);
    
    gettimeofday(&linear_end, NULL);
    linear_algorithm_time = (linear_end.tv_sec - linear_start.tv_sec) + (linear_end.tv_usec - linear_start.tv_usec) / 1000000.0;
    
    // Print comprehensive timing analysis
    struct timeval program_end;
    gettimeofday(&program_end, NULL);
    double total_program_time = (program_end.tv_sec - program_start_time.tv_sec) + (program_end.tv_usec - program_start_time.tv_usec) / 1000000.0;
    
    // Close unmapped report file and alert user if it contains data
    if (unmapped_report_file) {
        fclose(unmapped_report_file);
        unmapped_report_file = NULL;
        
        if (invalid_taxid_alignments > 0) {
            fprintf(stderr, "\n=== UNMAPPED ACCESSIONS WARNING ===\n");
            fprintf(stderr, "WARNING: %lld alignments had no taxonomic mapping\n", invalid_taxid_alignments);
            fprintf(stderr, "These reads were excluded from taxonomic analysis.\n");
            if (global_output_filename) {
                char base_name[512];
                strncpy(base_name, global_output_filename, sizeof(base_name) - 1);
                base_name[sizeof(base_name) - 1] = '\0';
                int len = strlen(base_name);
                if (len > 4 && strcmp(base_name + len - 4, ".txt") == 0) {
                    base_name[len - 4] = '\0';
                }
                fprintf(stderr, "Report file created: %s.missing_taxid.txt\n", base_name);
            } else {
                fprintf(stderr, "Report file created: taxonomy_report.missing_taxid.txt\n");
            }
            fprintf(stderr, "=====================================\n");
        }
    }
    
    
    // Print the taxonomy tree for debugging (only for small test runs)
    if (n_genomes <= 10 && total_alignments <= 100) {
        print_genome_taxonomy_tree();
    }
    
    // Clean up read mismatch table
    for (int i = 0; i < read_mismatch_table_size; i++) {
        ReadMismatchEntry *entry = read_mismatch_table[i];
        while (entry) {
            ReadMismatchEntry *next = entry->next;
            free_read_mismatches(entry->data);
            free(entry);
            entry = next;
        }
    }
    free(read_mismatch_table);
    read_mismatch_table = NULL;
    
#endif
}

// Hash function for genome discovery
unsigned int discovery_hash(const char *str) {
    unsigned int hash = 5381;  // djb2 hash
    while (*str) {
        hash = ((hash << 5) + hash) + *str++;
    }
    return hash % DISCOVERY_HASH_SIZE;
}

// Check if genome name exists in discovery hash table
int discovery_hash_contains(discovery_hash_entry_t **table, const char *name) {
    unsigned int hash = discovery_hash(name);
    discovery_hash_entry_t *entry = table[hash];
    while (entry) {
        if (strcmp(entry->name, name) == 0) {
            return 1;
        }
        entry = entry->next;
    }
    return 0;
}

// Add genome name to discovery hash table
void discovery_hash_add(discovery_hash_entry_t **table, const char *name) {
    unsigned int hash = discovery_hash(name);
    discovery_hash_entry_t *entry = malloc(sizeof(discovery_hash_entry_t));
    strcpy(entry->name, name);
    entry->next = table[hash];
    table[hash] = entry;
}

// Free discovery hash table
void discovery_hash_free(discovery_hash_entry_t **table) {
    for (int i = 0; i < DISCOVERY_HASH_SIZE; i++) {
        discovery_hash_entry_t *entry = table[i];
        while (entry) {
            discovery_hash_entry_t *next = entry->next;
            free(entry);
            entry = next;
        }
        table[i] = NULL;
    }
}

int auto_discover_genomes(const char *input_path, char ignore_char, int silent, int use_rg_tag, long long max_reads) {
    // Use hash table for fast O(1) duplicate detection instead of O(n) linear search
    discovery_hash_entry_t *discovery_table[DISCOVERY_HASH_SIZE] = {NULL};
    
    // Use dynamic allocation for discovered genomes (keep for final sorted output)
    int discovered_capacity = INITIAL_MAX_GENOMES;
    char (*discovered_genomes)[MAX_NAME_LEN] = calloc(discovered_capacity, sizeof(*discovered_genomes));
    if (!discovered_genomes) {
        fprintf(stderr, "Error: Failed to allocate memory for genome discovery\n");
        return -1;
    }
    int n_discovered = 0;
    
    if (!silent) {
        fprintf(stderr, "Auto-discovering genome names from BAM file(s)...\n");
    }
    
    // Check if input is a directory or single file
    if (is_directory(input_path)) {
        // Find all BAM files in directory
        if (find_bam_files(input_path) < 0) {
            return -1;
        }
        
        if (n_bam_files == 0) {
            fprintf(stderr, "Error: No BAM files found in directory: %s\n", input_path);
            return -1;
        }
        
        // Process each BAM file
        for (int file_idx = 0; file_idx < n_bam_files; file_idx++) {
            char *bam_file = bam_files[file_idx];
            
            // Open BAM file
            samFile *bam_fp = sam_open(bam_file, "r");
            if (!bam_fp) {
                fprintf(stderr, "Warning: Cannot open BAM file: %s\n", bam_file);
                continue;
            }
            
            sam_hdr_t *header = sam_hdr_read(bam_fp);
            if (!header) {
                fprintf(stderr, "Warning: Cannot read header from: %s\n", bam_file);
                sam_close(bam_fp);
                continue;
            }
            
            bam1_t *read = bam_init1();
            int reads_checked = 0;
            
            // Sample reads to find genome names
            // Use max_reads if specified, otherwise read entire file for complete discovery
            long long discovery_limit = max_reads;  // Use max_reads if specified, otherwise unlimited
            while (sam_read1(bam_fp, header, read) >= 0 && (discovery_limit <= 0 || reads_checked < discovery_limit)) {
                reads_checked++;
                
                const char *genome_identifier = NULL;
                
                if (use_rg_tag) {
                    // Get RG tag
                    uint8_t *rg_tag = bam_aux_get(read, "RG");
                    if (!rg_tag) continue;
                    genome_identifier = bam_aux2Z(rg_tag);
                } else {
                    // Use reference name
                    if (read->core.tid < 0) continue;  // Skip unmapped
                    genome_identifier = sam_hdr_tid2name(header, read->core.tid);
                }
                
                const char *genome_name = extract_genome_name(genome_identifier, ignore_char);
                
                // Check if this genome is already discovered using hash table (O(1) vs O(n))
                if (!discovery_hash_contains(discovery_table, genome_name)) {
                    // Grow array if needed
                    if (n_discovered >= discovered_capacity) {
                        discovered_capacity *= 2;
                        char (*new_array)[MAX_NAME_LEN] = realloc(discovered_genomes, 
                                                                  discovered_capacity * sizeof(*discovered_genomes));
                        if (!new_array) {
                            fprintf(stderr, "Error: Failed to grow discovered genomes array\n");
                            free(discovered_genomes);
                            return -1;
                        }
                        discovered_genomes = new_array;
                    }
                    strcpy(discovered_genomes[n_discovered++], genome_name);
                    discovery_hash_add(discovery_table, genome_name);  // Add to hash table for fast lookup
                    
                    // Progress reporting for genome discovery
                    if (!silent && n_discovered % 1000 == 0) {
                        fprintf(stderr, "Progress: Discovered %d genomes\n", n_discovered);
                    }
                }
            }
            
            bam_destroy1(read);
            sam_hdr_destroy(header);
            sam_close(bam_fp);
        }
    } else {
        // Single BAM file
        samFile *bam_fp = sam_open(input_path, "r");
        if (!bam_fp) {
            fprintf(stderr, "Error: Cannot open BAM file: %s\n", input_path);
            return -1;
        }
        
        sam_hdr_t *header = sam_hdr_read(bam_fp);
        if (!header) {
            fprintf(stderr, "Error: Cannot read BAM header\n");
            sam_close(bam_fp);
            return -1;
        }
        
        bam1_t *read = bam_init1();
        int reads_checked = 0;
        
        // Sample reads to find genome names
        // Use max_reads if specified, otherwise read entire file for complete discovery
        long long discovery_limit = max_reads;  // Use max_reads if specified, otherwise unlimited
        while (sam_read1(bam_fp, header, read) >= 0 && (discovery_limit <= 0 || reads_checked < discovery_limit)) {
            reads_checked++;
            
            const char *genome_identifier = NULL;
            
            if (use_rg_tag) {
                // Get RG tag
                uint8_t *rg_tag = bam_aux_get(read, "RG");
                if (!rg_tag) continue;
                genome_identifier = bam_aux2Z(rg_tag);
            } else {
                // Use reference name
                if (read->core.tid < 0) continue;  // Skip unmapped
                genome_identifier = sam_hdr_tid2name(header, read->core.tid);
            }
            
            const char *genome_name = extract_genome_name(genome_identifier, ignore_char);
            
            // Check if this genome is already discovered using hash table (O(1) vs O(n))
            if (!discovery_hash_contains(discovery_table, genome_name)) {
                // Grow array if needed
                if (n_discovered >= discovered_capacity) {
                    discovered_capacity *= 2;
                    char (*new_array)[MAX_NAME_LEN] = realloc(discovered_genomes, 
                                                              discovered_capacity * sizeof(*discovered_genomes));
                    if (!new_array) {
                        fprintf(stderr, "Error: Failed to grow discovered genomes array\n");
                        free(discovered_genomes);
                        return -1;
                    }
                    discovered_genomes = new_array;
                }
                strcpy(discovered_genomes[n_discovered++], genome_name);
                discovery_hash_add(discovery_table, genome_name);  // Add to hash table for fast lookup
                
                // Progress reporting for genome discovery  
                if (!silent && n_discovered % 1000 == 0) {
                    fprintf(stderr, "Progress: Discovered %d genomes\n", n_discovered);
                }
            }
        }
        
        bam_destroy1(read);
        sam_hdr_destroy(header);
        sam_close(bam_fp);
    }
    
    if (n_discovered == 0) {
        fprintf(stderr, "Error: No genomes found in BAM file(s)\n");
        free(discovered_genomes);
        return -1;
    }
    
    // Sort genome names for consistent output using efficient qsort (O(n log n) vs O(n²))
    qsort(discovered_genomes, n_discovered, sizeof(*discovered_genomes), 
          (int(*)(const void*, const void*))strcmp);
    
    // Copy to global genome list
    n_genomes = n_discovered;
    
    // Ensure we have enough space in global arrays
    while (n_genomes > max_genomes_allocated) {
        grow_genome_arrays();
    }
    
    for (int i = 0; i < n_genomes; i++) {
        strcpy(genome_names[i], discovered_genomes[i]);
        add_genome_to_hash(genome_names[i], i);
    }
    
    if (!silent) {
        fprintf(stderr, "Auto-discovered %d genomes\n", n_genomes);
    }
    
    // Free the temporary array and hash table
    free(discovered_genomes);
    discovery_hash_free(discovery_table);

    return n_genomes;  // Return the number of genomes discovered
}
#endif

// Generate compressed read ID (A, B, C, ..., Z, AA, AB, ...)
void generate_compressed_read_id(int read_number, char *buffer) {
    if (read_number == 0) {
        strcpy(buffer, "A");
        return;
    }
    
    char temp[32];
    int pos = 0;
    int n = read_number;
    
    while (n >= 0) {
        temp[pos++] = 'A' + (n % 26);
        n = n / 26 - 1;
        if (n < 0) break;
    }
    
    // Reverse the string
    int i;
    for (i = 0; i < pos; i++) {
        buffer[i] = temp[pos - 1 - i];
    }
    buffer[pos] = '\0';
}

// Load NCBI taxid mapping from names.dmp
typedef struct {
    char name[MAX_NAME_LEN];
    char taxid[32];
} taxid_entry_t;

taxid_entry_t *taxid_table = NULL;
int taxid_table_size = 0;
int taxid_table_capacity = 0;

void load_taxid_mapping(const char *taxid_file) {
    if (!taxid_file) return;
    
    FILE *fp = fopen(taxid_file, "r");
    if (!fp) {
        fprintf(stderr, "Warning: Could not open taxid file %s\n", taxid_file);
        return;
    }
    
    char line[4096];
    taxid_table_capacity = 10000;
    taxid_table = malloc(taxid_table_capacity * sizeof(taxid_entry_t));
    
    while (fgets(line, sizeof(line), fp)) {
        char *taxid = strtok(line, "\t|");
        if (!taxid) continue;
        
        char *name = strtok(NULL, "\t|");
        if (!name) continue;
        
        // Skip non-scientific names
        strtok(NULL, "\t|"); // unique name
        char *name_class = strtok(NULL, "\t|");
        if (!name_class || strstr(name_class, "scientific name") == NULL) continue;
        
        // Trim whitespace
        while (*name == ' ') name++;
        char *end = name + strlen(name) - 1;
        while (end > name && (*end == ' ' || *end == '\t' || *end == '\n')) end--;
        *(end + 1) = '\0';
        
        // Add to table
        if (taxid_table_size >= taxid_table_capacity) {
            taxid_table_capacity *= 2;
            taxid_table = realloc(taxid_table, taxid_table_capacity * sizeof(taxid_entry_t));
        }
        
        strncpy(taxid_table[taxid_table_size].name, name, MAX_NAME_LEN - 1);
        strncpy(taxid_table[taxid_table_size].taxid, taxid, 31);
        taxid_table_size++;
    }
    
    fclose(fp);
    fprintf(stderr, "Loaded %d taxid mappings\n", taxid_table_size);
}

// Find taxid for a genome name
const char *find_taxid_for_name(const char *genome_name) {
    if (!taxid_table) return NULL;
    
    // Try to find exact species names with word boundaries
    // This prevents matching "Hua" within "Du_Li_Huang"
    int i;
    for (i = 0; i < taxid_table_size; i++) {
        // Skip very short names that might cause false matches
        if (strlen(taxid_table[i].name) < 4) continue;
        
        // Look for the species name as a whole word
        const char *pos = strstr(genome_name, taxid_table[i].name);
        if (pos != NULL) {
            // Check if it's a word boundary before
            if (pos > genome_name) {
                char before = *(pos - 1);
                if (isalnum(before)) continue; // Not a word boundary
            }
            
            // Check if it's a word boundary after
            char after = *(pos + strlen(taxid_table[i].name));
            if (after != '\0' && isalnum(after)) continue; // Not a word boundary
            
            return taxid_table[i].taxid;
        }
    }
    
    return NULL;
}

// Generate compressed genome ID
void generate_genome_compressed_id(int genome_index, const char *genome_name) {
    // First try to find NCBI taxid
    const char *taxid = find_taxid_for_name(genome_name);
    
    if (taxid) {
        strcpy(genome_compressed_ids[genome_index], taxid);
        strcpy(genome_mappings[genome_index].taxid, taxid);
    } else {
        // Generate custom ID: A1, A2, ..., A999, B1, B2, ...
        int letter = (custom_id_counter - 1) / 999;
        int number = ((custom_id_counter - 1) % 999) + 1;
        sprintf(genome_compressed_ids[genome_index], "%c%d", 'A' + letter, number);
        strcpy(genome_mappings[genome_index].taxid, genome_compressed_ids[genome_index]);
        custom_id_counter++;
    }
    
    strcpy(genome_mappings[genome_index].full_name, genome_name);
}

// Write genome mapping file
void write_genome_mapping_file(const char *filename) {
    if (!filename) return;
    
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Could not create genome mapping file %s\n", filename);
        return;
    }
    
    fprintf(fp, "compressed_id\tfull_name\n");
    int i;
    for (i = 0; i < n_genomes; i++) {
        fprintf(fp, "%s\t%s\n", genome_mappings[i].taxid, genome_mappings[i].full_name);
    }
    
    fclose(fp);
    fprintf(stderr, "Wrote genome mapping to %s\n", filename);
}

// Write short names mapping file
void write_short_names_mapping_file(const char *filename) {
    if (!filename) return;
    
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Could not create short names mapping file %s\n", filename);
        return;
    }
    
    fprintf(fp, "short_name\tfull_name\n");
    int i;
    for (i = 0; i < n_genomes; i++) {
        fprintf(fp, "G%d\t%s\n", i + 1, genome_names[i]);
    }
    
    fclose(fp);
    fprintf(stderr, "Wrote short names mapping to %s\n", filename);
}

// Initialize dynamic genome arrays
void init_genome_arrays() {
    // PRE-ALLOCATE large capacity to avoid reallocation during processing
    max_genomes_allocated = 500000;  // Pre-allocate for 500K genomes instead of starting at 200
    
    genome_names = calloc(max_genomes_allocated, sizeof(*genome_names));
    genome_compressed_ids = calloc(max_genomes_allocated, sizeof(*genome_compressed_ids));
    genome_mappings = calloc(max_genomes_allocated, sizeof(genome_mapping_t));
    
    if (!genome_names || !genome_compressed_ids || !genome_mappings) {
        fprintf(stderr, "Error: Failed to allocate memory for genome arrays\n");
        exit(1);
    }
    
    // Memory pre-allocation completed silently
}

// Grow genome arrays when needed
void grow_genome_arrays() {
    pthread_mutex_lock(&genome_mutex);
    
    // Double-check if we still need to grow (another thread might have done it)
    if (n_genomes < max_genomes_allocated) {
        pthread_mutex_unlock(&genome_mutex);
        return;
    }
    
    int new_size = max_genomes_allocated * GENOME_GROWTH_FACTOR;
    
    char (*new_genome_names)[MAX_NAME_LEN] = realloc(genome_names, new_size * sizeof(*genome_names));
    char (*new_genome_compressed_ids)[32] = realloc(genome_compressed_ids, new_size * sizeof(*genome_compressed_ids));
    genome_mapping_t *new_genome_mappings = realloc(genome_mappings, new_size * sizeof(genome_mapping_t));
    
    if (!new_genome_names || !new_genome_compressed_ids || !new_genome_mappings) {
        fprintf(stderr, "Error: Failed to grow genome arrays from %d to %d\n", max_genomes_allocated, new_size);
        exit(1);
    }
    
    // Clear the new portion
    memset(&new_genome_names[max_genomes_allocated], 0, 
           (new_size - max_genomes_allocated) * sizeof(*new_genome_names));
    memset(&new_genome_compressed_ids[max_genomes_allocated], 0, 
           (new_size - max_genomes_allocated) * sizeof(*new_genome_compressed_ids));
    memset(&new_genome_mappings[max_genomes_allocated], 0, 
           (new_size - max_genomes_allocated) * sizeof(genome_mapping_t));
    
    genome_names = new_genome_names;
    genome_compressed_ids = new_genome_compressed_ids;
    genome_mappings = new_genome_mappings;
    max_genomes_allocated = new_size;
    
    // Suppress expansion messages since we pre-allocated
    // fprintf(stderr, "Info: Expanded genome capacity to %d\n", max_genomes_allocated);
    
    pthread_mutex_unlock(&genome_mutex);
}

// Post-order tree traversal test function
void post_order_traversal(TaxNode *node, FILE *fp) {
    if (!node) return;
    
    // First, traverse all children
    for (int i = 0; i < node->n_children; i++) {
        post_order_traversal(node->children[i], fp);
    }
    
    // Then output this node's relationships with its children
    for (int i = 0; i < node->n_children; i++) {
        // Write parent-child relationship
        // Format: parent_taxid child_taxid child_rank
        fprintf(fp, "%d\t%d\t%s\n", 
                node->taxid, 
                node->children[i]->taxid,
                node->children[i]->rank);
    }
}

// Test function to verify NCBI tree parsing
void test_taxonomy_tree_traversal(TaxonomyTree *tree) {
    if (!tree || !tree->nodes) {
        fprintf(stderr, "Error: No taxonomy tree to test\n");
        return;
    }
    
    // The root in NCBI taxonomy is always taxid 1
    TaxNode *root = tree->nodes[1];
    if (!root) {
        fprintf(stderr, "Error: No root node (taxid 1) found in tree\n");
        return;
    }
    
    FILE *fp = fopen("/Users/rasmus_nielsen/Desktop/MLIdentifier/Testing/tree_test.txt", "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot create tree_test.txt\n");
        return;
    }
    
    fprintf(stderr, "Performing post-order traversal of taxonomy tree...\n");
    fprintf(stderr, "Writing parent-child relationships to tree_test.txt\n");
    
    // Write header
    fprintf(fp, "# Parent_TaxID\tChild_TaxID\tChild_Rank\n");
    
    // Start post-order traversal from root (taxid 1)
    post_order_traversal(root, fp);
    
    fclose(fp);
    
    // Print some statistics
    fprintf(stderr, "Tree traversal complete:\n");
    fprintf(stderr, "  Maximum taxid in tree: %d\n", tree->max_taxid);
    fprintf(stderr, "  Active nodes in tree: %d\n", tree->n_active_nodes);
    fprintf(stderr, "  Root taxid: %d\n", root->taxid);
    fprintf(stderr, "  Root has %d direct children\n", root->n_children);
    
    // Exit after test as requested
    fprintf(stderr, "Tree test complete. Exiting.\n");
    exit(0);
}

// Function to print the taxonomy tree for specific genomes
void print_genome_taxonomy_tree(void) {
    if (!taxonomy_tree) {
        fprintf(stderr, "Error: No taxonomy tree available\n");
        return;
    }
    
    FILE *fp = fopen("/Users/rasmus_nielsen/Desktop/MLIdentifier/Testing/tree_5_genomes.txt", "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot create tree output file\n");
        return;
    }
    
    fprintf(fp, "# Taxonomy tree for the 5 discovered genomes\n");
    fprintf(fp, "# ==========================================\n\n");
    
    // First, list the genomes and their taxids
    fprintf(fp, "## Genomes and their direct taxids:\n");
    for (int i = 0; i < n_genomes; i++) {
        int taxid = get_taxid_from_index(i);
        if (taxid > 0 && taxid <= taxonomy_tree->max_taxid && taxonomy_tree->nodes[taxid]) {
            TaxNode *node = taxonomy_tree->nodes[taxid];
            fprintf(fp, "G%d (%s) -> TaxID %d (rank: %s)\n", 
                    i+1, genome_names[i], taxid, node->rank);
        }
    }
    fprintf(fp, "\n");
    
    // Now print the ancestry path for each genome
    fprintf(fp, "## Ancestry paths (from genome to root):\n");
    for (int i = 0; i < n_genomes; i++) {
        int taxid = get_taxid_from_index(i);
        if (taxid > 0 && taxid <= taxonomy_tree->max_taxid && taxonomy_tree->nodes[taxid]) {
            fprintf(fp, "\nG%d (%s) ancestry:\n", i+1, genome_names[i]);
            
            TaxNode *node = taxonomy_tree->nodes[taxid];
            int level = 0;
            while (node) {
                // Indent based on level
                for (int j = 0; j < level; j++) fprintf(fp, "  ");
                fprintf(fp, "-> TaxID %d (rank: %s)\n", 
                        node->taxid, node->rank);
                
                // Move to parent
                if (node->parent_taxid != node->taxid && node->parent_taxid > 0 && 
                    node->parent_taxid <= taxonomy_tree->max_taxid && 
                    taxonomy_tree->nodes[node->parent_taxid]) {
                    node = taxonomy_tree->nodes[node->parent_taxid];
                    level++;
                } else {
                    break;
                }
            }
        }
    }
    
    // Finally, show which taxonomic groups would appear in output
    fprintf(fp, "\n## Taxonomic groups in output (after | separator):\n");
    for (int i = 1; i <= taxonomy_tree->max_taxid; i++) {
        if (taxonomy_tree->nodes[i] && taxonomy_tree->nodes[i]->is_active) {
            TaxNode *node = taxonomy_tree->nodes[i];
            
            // Check if this is a leaf genome
            int is_leaf_genome = 0;
            for (int g = 0; g < n_genomes; g++) {
                if (get_taxid_from_index(g) == node->taxid) {
                    is_leaf_genome = 1;
                    break;
                }
            }
            
            if (!is_leaf_genome && node->n_leaves > 0) {
                fprintf(fp, "T%d:%s - contains genomes: ", node->taxid, node->rank);
                for (int j = 0; j < node->n_leaves; j++) {
                    if (j > 0) fprintf(fp, ", ");
                    fprintf(fp, "G%d", node->leaf_genomes[j] + 1);
                }
                fprintf(fp, "\n");
            }
        }
    }
    
    fclose(fp);
    fprintf(stderr, "Wrote taxonomy tree for 5 genomes to /Users/rasmus_nielsen/Desktop/MLIdentifier/Testing/tree_5_genomes.txt\n");
}

int main(int argc, char **argv) {
    // Start program timing for profiling
    gettimeofday(&program_start_time, NULL);
    
    options_t opts;
    memset(&opts, 0, sizeof(options_t));  // Initialize all fields to 0

    // Initialize genome arrays
    init_genome_arrays();
    
    // Initialize optimizations (will be enabled based on flags)
    init_optimizations();

    // Parse command line options
    if (parse_options(argc, argv, &opts) != 0) {
        return 1;
    }
    
    // Set global options pointer for access in other functions
    global_opts = &opts;
    global_output_filename = opts.output_file;

    // Load taxid mapping if provided
    if (opts.compress_output && opts.taxid_file) {
        load_taxid_mapping(opts.taxid_file);
    }
    
    // Load taxonomy data if requested or if testing
    if (opts.with_higher_taxa || opts.test_tree) {
        if (!opts.taxonomy_dir) {
            fprintf(stderr, "Error: --with-higher-taxa or --test-tree requires --taxonomy-dir\n");
            return 1;
        }

        if (opts.with_higher_taxa && !opts.acc2taxid_file) {
            fprintf(stderr, "Error: --with-higher-taxa requires --acc2taxid\n");
            return 1;
        }
        
        // Check for invalid dense + taxonomy combination
        if (opts.with_higher_taxa && opts.use_dense) {
            fprintf(stderr, "Error: Dense output format (--dense) is not compatible with taxonomic hierarchy analysis (--with-higher-taxa)\n");
            fprintf(stderr, "       Use sparse format instead by removing the --dense flag\n");
            return 1;
        }

        // Check for invalid sparse internal + dense output combination (without genome list)
        if (!opts.use_dense && !opts.genome_list &&
            (strcmp(opts.output_format, "dense") == 0 || strcmp(opts.output_format, "dense_damage") == 0)) {
            fprintf(stderr, "Error: Dense output formats (--format dense/dense_damage) require either:\n");
            fprintf(stderr, "       1. A genome list file (-g/--genomes), or\n");
            fprintf(stderr, "       2. Dense internal mode (--dense)\n");
            fprintf(stderr, "       \n");
            fprintf(stderr, "       The sparse internal mode cannot output dense format with dynamic genome discovery\n");
            fprintf(stderr, "       because dense format requires a fixed number of columns for all rows.\n");
            return 1;
        }
        
        // Check for invalid RG + taxonomy combination
        if (opts.with_higher_taxa && opts.use_rg_tag) {
            fprintf(stderr, "Error: RG tag identification (--use-rg) is not compatible with taxonomic analysis (--with-higher-taxa)\n");
            fprintf(stderr, "       Taxonomy mode requires reference names for accession→taxid lookup\n");
            fprintf(stderr, "       Remove --use-rg flag for taxonomic analysis\n");
            return 1;
        }

        // Validate consolidate-by-taxid requirements
        if (opts.consolidate_by_taxid) {
            if (!opts.acc2taxid_file) {
                fprintf(stderr, "Error: --consolidate-by-taxid requires --acc2taxid to be specified\n");
                return 1;
            }
            if (opts.use_dense) {
                fprintf(stderr, "Error: --consolidate-by-taxid is only compatible with sparse mode\n");
                fprintf(stderr, "       Remove --dense flag to use taxid consolidation\n");
                return 1;
            }
            if (strcmp(opts.output_format, "dense") == 0 || strcmp(opts.output_format, "dense_damage") == 0) {
                fprintf(stderr, "Error: --consolidate-by-taxid is not compatible with dense output formats\n");
                fprintf(stderr, "       Use --format sparse or --format sparse_damage\n");
                return 1;
            }
        }

        // Load NCBI taxonomy tree
        char nodes_file[512], names_file[512];
        snprintf(nodes_file, sizeof(nodes_file), "%s/nodes.dmp", opts.taxonomy_dir);
        snprintf(names_file, sizeof(names_file), "%s/names.dmp", opts.taxonomy_dir);
        
        if (!opts.silent) {
            fprintf(stderr, "Loading taxonomy tree from %s...\n", opts.taxonomy_dir);
        }
        
        taxonomy_tree = load_taxonomy_tree(nodes_file, names_file);
        if (!taxonomy_tree) {
            fprintf(stderr, "Error: Failed to load taxonomy tree\n");
            return 1;
        }
        
        
        // If just testing the tree, run test and exit
        if (opts.test_tree) {
            test_taxonomy_tree_traversal(taxonomy_tree);
            // test_taxonomy_tree_traversal exits, so we won't reach here
        }
        
    }

    // Load accession to taxid mapping if needed (for taxonomy or consolidation)
    if ((opts.with_higher_taxa || opts.consolidate_by_taxid) && opts.acc2taxid_file) {
        if (opts.with_higher_taxa) {
            // Set global taxid mode flag only for full taxonomy
            using_taxid_mode = 1;
        }

        if (!opts.silent) {
            fprintf(stderr, "Loading accession to taxid mapping from %s...\n", opts.acc2taxid_file);
        }

        int n_mappings = load_acc2taxid_mapping(opts.acc2taxid_file);
        if (n_mappings < 0) {
            fprintf(stderr, "Error: Failed to load acc2taxid mapping\n");
            if (taxonomy_tree) free_taxonomy_tree(taxonomy_tree);
            return 1;
        }

        if (!opts.silent) {
            fprintf(stderr, "Loaded %d accession to taxid mappings\n", n_mappings);
        }
    }
    
    if (opts.bam_mode) {
#ifdef WITH_HTSLIB
        // BAM mode: load genome list if provided
        if (opts.genome_list) {
            // Use provided genome list
            if (load_genome_list(opts.genome_list, opts.ignore_char, opts.use_simple_mode) < 0) {
                return 1;
            }
        } else {
            // Initialize empty genome list - will be populated dynamically
            n_genomes = 0;
            if (!opts.silent) {
                fprintf(stderr, "No genome list provided - genomes will be collected dynamically\n");
            }

            // For dense mode, we must discover all genomes upfront before writing the header
            if (opts.use_dense) {
                if (!opts.silent) {
                    fprintf(stderr, "Dense mode requires genome discovery before processing...\n");
                }

                // Call auto_discover_genomes with the correct parameters
                int n_discovered = auto_discover_genomes(opts.input_file, opts.ignore_char,
                                                        opts.silent, opts.use_rg_tag, opts.max_reads);
                if (n_discovered < 0) {
                    fprintf(stderr, "Error: Failed to auto-discover genomes\n");
                    return 1;
                }

                if (n_discovered == 0) {
                    fprintf(stderr, "Error: No genomes discovered in BAM file(s)\n");
                    return 1;
                }

                if (!opts.silent) {
                    fprintf(stderr, "Successfully discovered %d genomes for dense format\n", n_discovered);
                }
                // n_genomes is set by auto_discover_genomes
            }
        }

        // Check if input is a directory or single BAM file
        if (is_directory(opts.input_file)) {
            // Directory mode - find all BAM files
            if (find_bam_files(opts.input_file) < 0) {
                return 1;
            }
            
            if (n_bam_files == 0) {
                fprintf(stderr, "No BAM files found in directory: %s\n", opts.input_file);
                return 1;
            }
            
            if (opts.verbose) {
                printf("=== BAM to Mismatch Matrix Converter (Directory Mode) ===\n");
                printf("Directory: %s\n", opts.input_file);
                printf("Genome list: %s\n", opts.genome_list);
                printf("Output file: %s\n", opts.output_file);
                printf("Min MAPQ: %d\n", opts.min_mapq);
                printf("Loaded %d genomes\n", n_genomes);
                printf("Found %d BAM files\n\n", n_bam_files);
            }
            
            // Process multiple BAM files with unified approach
            if (opts.with_higher_taxa) {
                // Time the BAM processing phase for multiple files
                struct timeval bam_start, bam_end;
                gettimeofday(&bam_start, NULL);
                fprintf(stderr, "Starting BAM file processing...\n");

                process_multiple_bam_files_unified(&opts);

                gettimeofday(&bam_end, NULL);
                bam_processing_time = (bam_end.tv_sec - bam_start.tv_sec) + (bam_end.tv_usec - bam_start.tv_usec) / 1000000.0;
            } else {
                process_multiple_bam_files_unified(&opts);
            }
        } else {
            // Single BAM file mode
            n_bam_files = 1;
            strncpy(bam_files[0], opts.input_file, MAX_NAME_LEN - 1);
            bam_files[0][MAX_NAME_LEN - 1] = '\0';
            
            if (opts.verbose) {
                printf("=== BAM to Mismatch Matrix Converter (Single File Mode) ===\n");
                printf("BAM file: %s\n", opts.input_file);
                printf("Genome list: %s\n", opts.genome_list);
                printf("Output file: %s\n", opts.output_file);
                printf("Min MAPQ: %d\n", opts.min_mapq);
                printf("Loaded %d genomes\n\n", n_genomes);
            }
            
            // Process single BAM file - choose implementation based on options
            if (opts.with_higher_taxa) {
                // Time the BAM processing phase
                struct timeval bam_start, bam_end;
                gettimeofday(&bam_start, NULL);
                fprintf(stderr, "Starting BAM file processing...\n");
                
                process_bam_file_with_taxonomy(&opts);  // Taxonomy-aware processing

                gettimeofday(&bam_end, NULL);
                bam_processing_time = (bam_end.tv_sec - bam_start.tv_sec) + (bam_end.tv_usec - bam_start.tv_usec) / 1000000.0;

                // Print comprehensive timing analysis (moved from inside function to after timing calculation)
                struct timeval program_end;
                gettimeofday(&program_end, NULL);
                double total_program_time = (program_end.tv_sec - program_start_time.tv_sec) + (program_end.tv_usec - program_start_time.tv_usec) / 1000000.0;

                fprintf(stderr, "\n=== PERFORMANCE ANALYSIS ===\n");
                fprintf(stderr, "Total runtime: %.2fs\n", total_program_time);
                fprintf(stderr, "BAM processing (total): %.2fs (%.1f%%)\n", bam_processing_time, bam_processing_time/total_program_time*100);
                fprintf(stderr, "  - BAM file I/O: %.2fs (%.1f%%)\n", bam_processing_time - linear_algorithm_time, (bam_processing_time - linear_algorithm_time)/total_program_time*100);
                fprintf(stderr, "  - Linear taxonomy algorithm: %.2fs (%.1f%%)\n", linear_algorithm_time, linear_algorithm_time/total_program_time*100);
                fprintf(stderr, "Other phases (startup/cleanup): %.2fs (%.1f%%)\n",
                        total_program_time - bam_processing_time,
                        (total_program_time - bam_processing_time)/total_program_time*100);
                fprintf(stderr, "============================\n");
            } else if (opts.use_dense) {
                process_bam_file(&opts);  // Dense implementation
            } else {
                // Check for incompatible sparse internal + dense output combination
                if ((strcmp(opts.output_format, "dense") == 0 ||
                     strcmp(opts.output_format, "dense_damage") == 0)) {
                    fprintf(stderr, "Error: Dense output formats (--format dense/dense_damage) are not compatible\n");
                    fprintf(stderr, "       with sparse internal mode.\n");
                    fprintf(stderr, "\n");
                    fprintf(stderr, "       Solutions:\n");
                    fprintf(stderr, "       1. Use dense internal mode: add --dense flag\n");
                    fprintf(stderr, "       2. Use sparse output format: --format sparse or --format sparse_damage\n");
                    return 1;
                }
                process_bam_file_sparse_optimized(&opts);  // OPTIMIZED Sparse implementation
            }
        }
#else
        fprintf(stderr, "Error: BAM mode not available - program compiled without HTSlib\n");
        return 1;
#endif
    } else {
        // Text mode: process text file
        if (opts.verbose) {
            printf("=== Mismatch Matrix Text Processor ===\n");
            printf("Input file: %s\n", opts.input_file);
            printf("Output file: %s\n", opts.output_file);
            printf("\n");
        }

        // Process text file
        process_text_file(&opts);
    }

    // Write genome mapping file if compression is enabled
    if (opts.compress_output) {
        // If no genome map file specified, create one based on output filename
        char genome_map_filename[1024];
        if (opts.genome_map_file) {
            write_genome_mapping_file(opts.genome_map_file);
        } else {
            // Remove extension from output file and add _genome_mapping.txt
            char *base_name = strdup(opts.output_file);
            char *dot = strrchr(base_name, '.');
            if (dot && dot != base_name) {
                *dot = '\0';
            }
            snprintf(genome_map_filename, sizeof(genome_map_filename), "%s_genome_mapping.txt", base_name);
            write_genome_mapping_file(genome_map_filename);
            if (opts.verbose) {
                printf("Created genome mapping file: %s\n", genome_map_filename);
            }
            free(base_name);
        }
    }
    
    // Write short names mapping file if key file creation is requested and using short names
    // BUT not in taxid mode or consolidate mode (which output actual taxids, not G1/G2 synonyms)
    if (opts.create_key_file && opts.short_names && !using_taxid_mode && !opts.consolidate_by_taxid) {
        // Only write if genome_names array is actually populated
        // In sparse mode without genome list, it's empty and process_temp_file.c already wrote the file
        if (n_genomes > 0 && strlen(genome_names[0]) > 0) {
            // Create mapping file based on output filename
            char *base_name = strdup(opts.output_file);
            char *dot = strrchr(base_name, '.');
            if (dot && dot != base_name) {
                *dot = '\0';
            }
            char short_names_map_filename[1024];
            snprintf(short_names_map_filename, sizeof(short_names_map_filename), "%s_genome_key.txt", base_name);
            write_short_names_mapping_file(short_names_map_filename);
            if (opts.verbose) {
                printf("Created genome key file: %s\n", short_names_map_filename);
            }
            free(base_name);
        }
    }
    
    if (opts.verbose) {
        printf("Processing completed successfully!\n");
        printf("Output written to: %s\n", opts.output_file);
    }

    // Clean up taxonomy data if loaded
    if (opts.with_higher_taxa) {
        if (taxonomy_tree) {
            free_taxonomy_tree(taxonomy_tree);
        }
        if (acc2taxid_table) {
            free_acc2taxid_mapping();
        }
    }
    
    // Clean up dynamic allocations
    free(genome_names);
    free(genome_compressed_ids);
    free(genome_mappings);
    
    // Clean up hash table
    for (int i = 0; i < HASH_TABLE_SIZE; i++) {
        genome_hash_entry_t *entry = genome_hash_table[i];
        while (entry) {
            genome_hash_entry_t *next = entry->next;
            free(entry);
            entry = next;
        }
    }

    return 0;
}