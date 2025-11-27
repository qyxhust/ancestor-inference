// Optimized temp file processing with hash table - replaces O(g²) with O(g log g)
// This implementation uses a hash table for genome discovery and qsort for sorting

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Hash table configuration
#define INITIAL_HASH_SIZE 1024
#define HASH_LOAD_FACTOR 0.75

typedef struct GenomeInfo {
    char genome_name[256];
    int count;
    int index;  // Position in sorted array
    struct GenomeInfo *next;  // For hash collision chaining
} GenomeInfo;

typedef struct {
    GenomeInfo **buckets;
    int size;           // Number of buckets
    int count;          // Number of entries
    int capacity;       // Current capacity for resizing
} HashTable;

// Simple string hash function (djb2 algorithm) for temp file processing
static unsigned int hash_string_tempfile(const char *str, int table_size) {
    unsigned int hash = 5381;
    int c;
    while ((c = *str++)) {
        hash = ((hash << 5) + hash) + c; // hash * 33 + c
    }
    return hash % table_size;
}

// Initialize hash table
static HashTable* hash_table_init(int initial_size) {
    HashTable *ht = malloc(sizeof(HashTable));
    ht->size = initial_size;
    ht->count = 0;
    ht->capacity = initial_size * HASH_LOAD_FACTOR;
    ht->buckets = calloc(initial_size, sizeof(GenomeInfo*));
    return ht;
}

// Find genome in hash table
static GenomeInfo* hash_table_find(HashTable *ht, const char *genome_name) {
    unsigned int index = hash_string_tempfile(genome_name, ht->size);
    GenomeInfo *entry = ht->buckets[index];
    
    while (entry) {
        if (strcmp(entry->genome_name, genome_name) == 0) {
            return entry;
        }
        entry = entry->next;
    }
    return NULL;
}

// Add or update genome in hash table
static GenomeInfo* hash_table_insert(HashTable *ht, const char *genome_name) {
    GenomeInfo *existing = hash_table_find(ht, genome_name);
    if (existing) {
        existing->count++;
        return existing;
    }
    
    // Create new entry
    GenomeInfo *new_entry = malloc(sizeof(GenomeInfo));
    strcpy(new_entry->genome_name, genome_name);
    new_entry->count = 1;
    new_entry->index = -1;  // Will be set during sorting phase
    new_entry->next = NULL;
    
    // Insert into hash table
    unsigned int index = hash_string_tempfile(genome_name, ht->size);
    new_entry->next = ht->buckets[index];
    ht->buckets[index] = new_entry;
    ht->count++;
    
    return new_entry;
}

// Convert hash table to sorted array
static GenomeInfo** hash_table_to_sorted_array(HashTable *ht, int *count) {
    GenomeInfo **array = malloc(ht->count * sizeof(GenomeInfo*));
    int array_index = 0;
    
    // Collect all entries from hash table
    for (int i = 0; i < ht->size; i++) {
        GenomeInfo *entry = ht->buckets[i];
        while (entry) {
            array[array_index++] = entry;
            entry = entry->next;
        }
    }
    
    *count = ht->count;
    return array;
}

// Comparison function for qsort
static int compare_genomes(const void *a, const void *b) {
    const GenomeInfo *genome_a = *(const GenomeInfo**)a;
    const GenomeInfo *genome_b = *(const GenomeInfo**)b;
    return strcmp(genome_a->genome_name, genome_b->genome_name);
}

// Free hash table
static void hash_table_free(HashTable *ht) {
    for (int i = 0; i < ht->size; i++) {
        GenomeInfo *entry = ht->buckets[i];
        while (entry) {
            GenomeInfo *temp = entry;
            entry = entry->next;
            free(temp);
        }
    }
    free(ht->buckets);
    free(ht);
}

int process_temp_file_optimized(const char *temp_file, const char *output_file, int is_damage_format) {
    // Access to global options for enforce_dense_strict flag
    extern options_t *global_opts;
    extern int using_taxid_mode;

    // Check if we should use actual genome/taxid names or abbreviations
    int use_actual_names = (global_opts && global_opts->consolidate_by_taxid) || using_taxid_mode;
    FILE *temp_fp = fopen(temp_file, "r");
    if (!temp_fp) {
        fprintf(stderr, "Error: Cannot open temp file %s\n", temp_file);
        return -1;  // Error: return -1 for failure
    }
    
    // Phase 1: Collect unique genomes using hash table - O(g) average case
    printf("Phase 1: Collecting unique genomes from temp file (hash table optimization)...\n\n");
    
    HashTable *genome_hash = hash_table_init(INITIAL_HASH_SIZE);
    
    char line[65536];
    while (fgets(line, sizeof(line), temp_fp)) {
        char *token = strtok(line, "\t");
        if (!token) continue; // Skip empty lines
        
        // Skip read_id
        token = strtok(NULL, "\t");
        if (!token) continue;
        
        // Skip total_count
        token = strtok(NULL, "\t");
        
        // Process genome entries
        while (token) {
            // O(1) average case hash table lookup and insertion
            hash_table_insert(genome_hash, token);
            
            // Skip the mismatch values based on format
            if (is_damage_format) {
                // Skip nd, md, mb values
                for (int i = 0; i < 3; i++) {
                    token = strtok(NULL, "\t");
                    if (!token) break;
                }
            } else {
                // Skip single mismatch value
                token = strtok(NULL, "\t");
            }
            
            // Get next genome name
            token = strtok(NULL, "\t");
        }
    }
    
    printf("Phase 1 complete: Found %d unique genomes\n", genome_hash->count);
    
    // Phase 2: Sort genomes alphabetically - O(g log g) using qsort
    printf("Phase 2: Sorting genomes (qsort optimization)...\n");
    
    int n_unique_genomes;
    GenomeInfo **sorted_genomes = hash_table_to_sorted_array(genome_hash, &n_unique_genomes);
    qsort(sorted_genomes, n_unique_genomes, sizeof(GenomeInfo*), compare_genomes);
    
    // Set index values for O(1) lookup during phase 3
    for (int i = 0; i < n_unique_genomes; i++) {
        sorted_genomes[i]->index = i;
    }
    
    // Phase 3: Write output with abbreviations - O(1) hash table lookup per genome
    printf("Phase 3: Writing final output with abbreviations (hash table lookup)...\n\n");
    
    FILE *out_fp = fopen(output_file, "w");
    if (!out_fp) {
        fprintf(stderr, "Error: Cannot create output file %s\n", output_file);
        fclose(temp_fp);
        hash_table_free(genome_hash);
        free(sorted_genomes);
        return -1;  // Error: return -1 for failure
    }
    
    // Write header
    fprintf(out_fp, "read_id\ttotal_count");

    if (use_actual_names) {
        // Write actual genome/taxid names for consolidate_by_taxid or taxid mode
        for (int i = 0; i < n_unique_genomes; i++) {
            fprintf(out_fp, "\t%s", sorted_genomes[i]->genome_name);
        }
    } else {
        // Write abbreviations (G1, G2, etc.) for normal mode
        for (int i = 0; i < n_unique_genomes; i++) {
            fprintf(out_fp, "\tG%d", i + 1);
        }
    }
    fprintf(out_fp, "\n");
    
    // Rewind temp file to process again
    rewind(temp_fp);
    
    // Write data rows with genome abbreviations
    int lines_written = 0;
    while (fgets(line, sizeof(line), temp_fp)) {
        char line_copy[65536];
        strcpy(line_copy, line);
        
        char *token = strtok(line_copy, "\t");
        if (!token) continue;

        char *read_id = token;  // Save read_id for potential output

        // Get total_count
        token = strtok(NULL, "\t");
        if (!token) continue;
        char *total_count = token;  // Save total_count for potential output

        // Check --enforce-dense_strict filtering before output
        if (global_opts && global_opts->enforce_dense_strict) {
            // Count how many genomes this read aligns to
            int alignment_count = 0;
            char *genome_token = strtok(NULL, "\t");
            while (genome_token) {
                // Each genome entry has either 1 value (no damage) or 3 values (damage)
                alignment_count++;

                // Skip the mismatch value(s) for this genome
                if (is_damage_format) {
                    // Skip nd, md, mb values (3 values)
                    for (int i = 0; i < 3; i++) {
                        genome_token = strtok(NULL, "\t");
                        if (!genome_token) break;
                    }
                } else {
                    // Skip single mismatch value
                    genome_token = strtok(NULL, "\t");
                }

                // Get next genome name
                genome_token = strtok(NULL, "\t");
            }

            // Skip this read if it doesn't align to ALL genomes
            if (alignment_count != n_unique_genomes) {
                continue;  // Don't output this read
            }
        }

        // Output the read (re-parse the line since strtok modified it)
        strcpy(line_copy, line);
        token = strtok(line_copy, "\t");

        // Write read_id
        fprintf(out_fp, "%s", token);

        // Write total_count
        token = strtok(NULL, "\t");
        if (!token) continue;
        fprintf(out_fp, "\t%s", token);
        
        // Process genome entries and replace with abbreviations
        token = strtok(NULL, "\t");
        while (token) {
            // O(1) average case hash table lookup instead of O(g) linear search
            GenomeInfo *genome_info = hash_table_find(genome_hash, token);
            
            if (genome_info && genome_info->index >= 0) {
                if (use_actual_names) {
                    // Write actual genome/taxid name for consolidate_by_taxid or taxid mode
                    fprintf(out_fp, "\t%s", token);
                } else {
                    // Write abbreviation for normal mode
                    fprintf(out_fp, "\tG%d", genome_info->index + 1);
                }
                
                // Write mismatch values based on format
                if (is_damage_format) {
                    // Write nd, md, mb values
                    for (int i = 0; i < 3; i++) {
                        token = strtok(NULL, "\t");
                        if (token) {
                            fprintf(out_fp, "\t%s", token);
                        }
                    }
                } else {
                    // Write single mismatch value
                    token = strtok(NULL, "\t");
                    if (token) {
                        // Remove newline if present
                        char *newline = strchr(token, '\n');
                        if (newline) *newline = '\0';
                        fprintf(out_fp, "\t%s", token);
                    }
                }
            }
            
            // Get next genome name
            token = strtok(NULL, "\t");
        }
        fprintf(out_fp, "\n");
        lines_written++;
    }
    
    printf("Phase 3 complete: Wrote %d lines\n", lines_written);

    // Only create genome key file when using abbreviations (not for taxid mode)
    if (!use_actual_names) {
        char key_filename[512];
        snprintf(key_filename, sizeof(key_filename), "%s_genome_key.txt", output_file);
        FILE *key_fp = fopen(key_filename, "w");
        if (key_fp) {
            fprintf(key_fp, "short_name\tfull_name\n");
            for (int i = 0; i < n_unique_genomes; i++) {
                fprintf(key_fp, "G%d\t%s\n", i + 1, sorted_genomes[i]->genome_name);
            }
            fclose(key_fp);
            printf("Wrote genome mapping to %s\n", key_filename);
        }
    }
    
    fclose(temp_fp);
    fclose(out_fp);
    hash_table_free(genome_hash);
    free(sorted_genomes);
    
    printf("Processing complete!\n");
    return n_unique_genomes;
}

// Wrapper for backwards compatibility
int rewrite_sparse_with_header_optimized(const char *temp_file, const char *output_file, int is_damage_format) {
    return process_temp_file_optimized(temp_file, output_file, is_damage_format);
}