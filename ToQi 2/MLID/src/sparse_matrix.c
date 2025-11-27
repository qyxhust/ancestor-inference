#define _GNU_SOURCE  // For strdup
#include "sparse_matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Define maximum line buffer size (increased to handle many genomes)
#define MAX_LINE_BUFFER 65536

// Helper function to read a line with dynamic buffer allocation (optimized for long lines)
static char* read_line_dynamic(FILE *file, size_t *buffer_size) {
    size_t local_buffer_size;
    
    if (!buffer_size || *buffer_size == 0) {
        local_buffer_size = 1048576;  // 1MB default
        buffer_size = &local_buffer_size;
    }
    
    char *line = malloc(*buffer_size);
    if (!line) return NULL;
    
    size_t total_read = 0;
    
    // Read using fgets in chunks, growing buffer as needed
    while (1) {
        size_t remaining = *buffer_size - total_read;
        
        if (!fgets(line + total_read, remaining, file)) {
            if (total_read == 0) {
                free(line);
                return NULL;  // EOF or error
            }
            // Partial line at EOF
            break;
        }
        
        size_t len = strlen(line + total_read);
        total_read += len;
        
        // Check if we got the full line
        if (total_read > 0 && line[total_read - 1] == '\n') {
            line[total_read - 1] = '\0';  // Remove newline
            return line;
        }
        
        // If we filled the buffer but didn't get a newline, grow it
        if (total_read >= *buffer_size - 1) {
            *buffer_size *= 2;
            char *new_buffer = realloc(line, *buffer_size);
            if (!new_buffer) {
                free(line);
                return NULL;
            }
            line = new_buffer;
        }
    }
    
    return line;
}

// Allocate sparse alignment structure
sparse_alignment_t* sparse_alignment_alloc(int n_reads, int n_genomes, int nnz) {
    sparse_alignment_t *sparse = malloc(sizeof(sparse_alignment_t));
    if (!sparse) return NULL;
    
    sparse->n_reads = n_reads;
    sparse->n_genomes = n_genomes;
    sparse->nnz = nnz;
    
    // Allocate arrays
    sparse->row_ptr = calloc(n_reads + 1, sizeof(int));
    sparse->col_indices = malloc(nnz * sizeof(int));
    sparse->n_values = malloc(nnz * sizeof(short));
    sparse->d_values = malloc(nnz * sizeof(short));
    
    // Initialize damage arrays to NULL - will be allocated if needed
    sparse->nd_values = NULL;
    sparse->md_values = NULL;
    sparse->mb_values = NULL;
    
    // Initialize taxonomic fields to NULL/default values
    sparse->has_taxonomic_groups = 0;
    sparse->taxonomic_groups = NULL;
    sparse->n_taxonomic_groups = 0;
    sparse->read_taxonomic_data = NULL;
    sparse->use_float_mb_format = 0;
    
    if (!sparse->row_ptr || !sparse->col_indices || 
        !sparse->n_values || !sparse->d_values) {
        sparse_alignment_free(sparse);
        return NULL;
    }
    
    return sparse;
}

// Free sparse alignment structure
void sparse_alignment_free(sparse_alignment_t *sparse) {
    if (!sparse) return;
    
    free(sparse->row_ptr);
    free(sparse->col_indices);
    free(sparse->n_values);
    free(sparse->d_values);
    
    // Free damage arrays if allocated
    free(sparse->nd_values);
    free(sparse->md_values);
    free(sparse->mb_values);
    
    // Free taxonomic groups data
    if (sparse->taxonomic_groups) {
        for (int i = 0; i < sparse->n_taxonomic_groups; i++) {
            free(sparse->taxonomic_groups[i].rank);
            free(sparse->taxonomic_groups[i].name);
        }
        free(sparse->taxonomic_groups);
    }
    
    // Free per-read taxonomic data
    if (sparse->read_taxonomic_data) {
        for (int i = 0; i < sparse->n_reads; i++) {
            free(sparse->read_taxonomic_data[i].taxon_indices);
            free(sparse->read_taxonomic_data[i].nd_values);
            free(sparse->read_taxonomic_data[i].md_values);
            free(sparse->read_taxonomic_data[i].mds_values);
            free(sparse->read_taxonomic_data[i].mb_values);
            free(sparse->read_taxonomic_data[i].mbs_values);
            free(sparse->read_taxonomic_data[i].total_all_values);
            free(sparse->read_taxonomic_data[i].total_some_values);
            // Free new float arrays
            free(sparse->read_taxonomic_data[i].mb_values_float);
            free(sparse->read_taxonomic_data[i].mbs_values_float);
            free(sparse->read_taxonomic_data[i].total_all_values_float);
        }
        free(sparse->read_taxonomic_data);
    }
    
    free(sparse);
}

// Allocate damage arrays for sparse alignment structure
int sparse_alignment_alloc_damage(sparse_alignment_t *sparse) {
    if (!sparse || sparse->nnz <= 0) return -1;
    
    // Allocate damage arrays if not already allocated
    if (!sparse->nd_values) {
        sparse->nd_values = malloc(sparse->nnz * sizeof(short));
        if (!sparse->nd_values) return -1;
    }
    
    if (!sparse->md_values) {
        sparse->md_values = malloc(sparse->nnz * sizeof(short));
        if (!sparse->md_values) {
            free(sparse->nd_values);
            sparse->nd_values = NULL;
            return -1;
        }
    }
    
    if (!sparse->mb_values) {
        sparse->mb_values = malloc(sparse->nnz * sizeof(short));
        if (!sparse->mb_values) {
            free(sparse->nd_values);
            free(sparse->md_values);
            sparse->nd_values = NULL;
            sparse->md_values = NULL;
            return -1;
        }
    }
    
    return 0;
}

// Get number of non-zero entries in a row
int sparse_get_row_nnz(const sparse_alignment_t *sparse, int row) {
    if (row < 0 || row >= sparse->n_reads) return 0;
    return sparse->row_ptr[row + 1] - sparse->row_ptr[row];
}

// Find column index for a specific row-column pair
// Returns -1 if not found
int sparse_find_col_index(const sparse_alignment_t *sparse, int row, int col) {
    if (row < 0 || row >= sparse->n_reads || col < 0 || col >= sparse->n_genomes) {
        return -1;
    }
    
    int start = sparse->row_ptr[row];
    int end = sparse->row_ptr[row + 1];
    
    for (int idx = start; idx < end; idx++) {
        if (sparse->col_indices[idx] == col) {
            return idx;
        }
    }
    
    return -1;
}

// Get n_value for a specific row-column pair (returns 0 if not aligned)
double sparse_get_n_value(const sparse_alignment_t *sparse, int row, int col) {
    int idx = sparse_find_col_index(sparse, row, col);
    return (idx >= 0) ? sparse->n_values[idx] : 0.0;
}

// Get d_value for a specific row-column pair (returns 0 if not aligned)
double sparse_get_d_value(const sparse_alignment_t *sparse, int row, int col) {
    int idx = sparse_find_col_index(sparse, row, col);
    return (idx >= 0) ? sparse->d_values[idx] : 0.0;
}

// Convert dense matrices to sparse format
sparse_alignment_t* dense_to_sparse(short **n_matrix, short **d_matrix, 
                                   int n_reads, int n_genomes) {
    // First pass: count non-zero entries (alignments where n_matrix[i][j] != -1)
    int nnz = 0;
    for (int i = 0; i < n_reads; i++) {
        for (int j = 0; j < n_genomes; j++) {
            if (n_matrix[i][j] >= 0) {  // Valid alignment
                nnz++;
            }
        }
    }
    
    // Allocate sparse structure
    sparse_alignment_t *sparse = sparse_alignment_alloc(n_reads, n_genomes, nnz);
    if (!sparse) return NULL;
    
    // Second pass: fill sparse structure
    int current_idx = 0;
    sparse->row_ptr[0] = 0;
    
    for (int i = 0; i < n_reads; i++) {
        for (int j = 0; j < n_genomes; j++) {
            if (n_matrix[i][j] >= 0) {  // Valid alignment
                sparse->col_indices[current_idx] = j;
                sparse->n_values[current_idx] = n_matrix[i][j];
                sparse->d_values[current_idx] = d_matrix[i][j];
                current_idx++;
            }
        }
        sparse->row_ptr[i + 1] = current_idx;
    }
    
    return sparse;
}

// Convert sparse format back to dense matrices
void sparse_to_dense(const sparse_alignment_t *sparse, 
                     short **n_matrix, short **d_matrix) {
    // Initialize all entries to -1 (no alignment)
    for (int i = 0; i < sparse->n_reads; i++) {
        for (int j = 0; j < sparse->n_genomes; j++) {
            n_matrix[i][j] = -1;
            d_matrix[i][j] = -1;
        }
    }
    
    // Fill in aligned entries
    for (int i = 0; i < sparse->n_reads; i++) {
        int start = sparse->row_ptr[i];
        int end = sparse->row_ptr[i + 1];
        
        for (int idx = start; idx < end; idx++) {
            int j = sparse->col_indices[idx];
            n_matrix[i][j] = sparse->n_values[idx];
            d_matrix[i][j] = sparse->d_values[idx];
        }
    }
}

// Multiply sparse matrix by vector (for proportion calculations)
void sparse_matrix_multiply_vector(const sparse_alignment_t *sparse, 
                                  const double *vector, double *result) {
    for (int i = 0; i < sparse->n_reads; i++) {
        result[i] = 0.0;
        
        int start = sparse->row_ptr[i];
        int end = sparse->row_ptr[i + 1];
        
        for (int idx = start; idx < end; idx++) {
            int j = sparse->col_indices[idx];
            result[i] += sparse->n_values[idx] * vector[j];
        }
    }
}

// Sum all values in sparse matrix
double sparse_matrix_sum(const sparse_alignment_t *sparse) {
    double sum = 0.0;
    for (int idx = 0; idx < sparse->nnz; idx++) {
        sum += sparse->n_values[idx];
    }
    return sum;
}

// Calculate row sums for sparse matrix
void sparse_matrix_row_sums(const sparse_alignment_t *sparse, double *row_sums) {
    for (int i = 0; i < sparse->n_reads; i++) {
        row_sums[i] = 0.0;
        
        int start = sparse->row_ptr[i];
        int end = sparse->row_ptr[i + 1];
        
        for (int idx = start; idx < end; idx++) {
            row_sums[i] += sparse->n_values[idx];
        }
    }
}

// Detect file format by examining first data line
data_format_t detect_file_format(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) return FORMAT_AUTO;
    
    char line[MAX_LINE_BUFFER];
    int lines_read = 0;
    
    // Skip header line and read first data line
    while (fgets(line, sizeof(line), file) && lines_read < 2) {
        // Check if line was truncated
        size_t line_len = strlen(line);
        if (line_len > 0 && line[line_len - 1] != '\n') {
            fprintf(stderr, "Warning: Line %d exceeds buffer size of %d characters.\n", lines_read + 1, MAX_LINE_BUFFER);
            fclose(file);
            return FORMAT_AUTO;
        }
        lines_read++;
        if (lines_read == 2) {  // First data line
            // Always return FORMAT_AUTO to let the proper format detection in io_utils.c handle it
            fclose(file);
            return FORMAT_AUTO;
        }
    }
    
    fclose(file);
    return FORMAT_AUTO;
}

// Parse sparse damage format file
sparse_alignment_t* parse_sparse_damage_file(const char *filename, char ***genome_names) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Cannot open file %s\n", filename);
        return NULL;
    }
    
    // Use dynamic line reader for header
    size_t header_buffer_size = 0;
    char *header_line = read_line_dynamic(file, &header_buffer_size);
    if (!header_line) {
        fprintf(stderr, "Error: Cannot read header from %s\n", filename);
        fclose(file);
        return NULL;
    }
    
    // Parse header to get genome names and detect taxonomic groups
    char *header_copy = strdup(header_line);
    char *token = strtok(header_copy, "\t\n");  // Skip read_id
    token = strtok(NULL, "\t\n");  // Skip total_count
    
    // First, check if this header has taxonomic groups (contains "|")
    int has_taxonomic_groups = 0;
    char *pipe_pos = strchr(header_line, '|');
    if (pipe_pos != NULL) {
        has_taxonomic_groups = 1;
    }
    
    // Count and collect genome names from header (before "|" if present)
    char **temp_names = malloc(10000 * sizeof(char*));  // Start with 10K genomes
    int names_capacity = 10000;
    int n_genomes = 0;
    
    while ((token = strtok(NULL, "\t\n")) != NULL) {
        // Stop parsing genome names when we hit the "|" separator
        if (strcmp(token, "|") == 0) {
            break;
        }
        
        if (n_genomes >= names_capacity) {
            names_capacity *= 2;
            char **new_names = realloc(temp_names, names_capacity * sizeof(char*));
            if (!new_names) {
                for (int i = 0; i < n_genomes; i++) {
                    free(temp_names[i]);
                }
                free(temp_names);
                free(header_copy);
                free(header_line);
                fclose(file);
                return NULL;
            }
            temp_names = new_names;
        }
        temp_names[n_genomes] = strdup(token);
        n_genomes++;
    }
    
    // Parse taxonomic groups if present (after "|")
    taxonomic_group_t *taxonomic_groups = NULL;
    int n_taxonomic_groups = 0;
    
    if (has_taxonomic_groups) {
        // Continue parsing tokens after "|" for taxonomic groups
        char **temp_taxa = malloc(10000 * sizeof(char*));
        int taxa_capacity = 10000;
        
        while ((token = strtok(NULL, "\t\n")) != NULL) {
            if (n_taxonomic_groups >= taxa_capacity) {
                taxa_capacity *= 2;
                char **new_taxa = realloc(temp_taxa, taxa_capacity * sizeof(char*));
                if (!new_taxa) {
                    for (int i = 0; i < n_taxonomic_groups; i++) {
                        free(temp_taxa[i]);
                    }
                    free(temp_taxa);
                    break;  // Continue without taxonomic groups
                }
                temp_taxa = new_taxa;
            }
            temp_taxa[n_taxonomic_groups] = strdup(token);
            n_taxonomic_groups++;
        }
        
        // Create taxonomic_group_t structures
        if (n_taxonomic_groups > 0) {
            taxonomic_groups = malloc(n_taxonomic_groups * sizeof(taxonomic_group_t));
            if (taxonomic_groups) {
                for (int i = 0; i < n_taxonomic_groups; i++) {
                    // Parse T{taxid}:{rank} format
                    char *taxon_str = temp_taxa[i];
                    if (taxon_str[0] == 'T') {
                        char *colon_pos = strchr(taxon_str, ':');
                        if (colon_pos) {
                            // Store full name BEFORE modifying the string
                            taxonomic_groups[i].name = strdup(taxon_str);  // Full name (T1:genus)
                            
                            // Extract taxid and rank
                            *colon_pos = '\0';  // Split string at colon
                            taxonomic_groups[i].taxid = atoi(taxon_str + 1);  // Skip 'T'
                            taxonomic_groups[i].rank = strdup(colon_pos + 1);
                            *colon_pos = ':';  // Restore original string
                            
                        } else {
                            // Malformed taxonomic group
                            taxonomic_groups[i].taxid = -1;
                            taxonomic_groups[i].rank = strdup("unknown");
                            taxonomic_groups[i].name = strdup(taxon_str);
                        }
                    } else {
                        // Not a valid taxonomic group
                        taxonomic_groups[i].taxid = -1;
                        taxonomic_groups[i].rank = strdup("unknown");
                        taxonomic_groups[i].name = strdup(taxon_str);
                    }
                }
            }
        }
        
        // Free temporary taxa names
        for (int i = 0; i < n_taxonomic_groups; i++) {
            free(temp_taxa[i]);
        }
        free(temp_taxa);
    }
    free(header_copy);
    free(header_line);  // Free the dynamically allocated header
    
    if (n_genomes == 0) {
        fprintf(stderr, "Error: No genomes found in header\n");
        fclose(file);
        return NULL;
    }
    
    // Allocate genome names array
    *genome_names = malloc(n_genomes * sizeof(char*));
    for (int i = 0; i < n_genomes; i++) {
        (*genome_names)[i] = temp_names[i];  // Transfer ownership
    }
    free(temp_names);  // Free the array of pointers
    
    // Print progress for header parsing
    if (has_taxonomic_groups) {
        printf("Read header: %d genomes and %d taxonomic groups\n", n_genomes, n_taxonomic_groups);
    } else {
        printf("Read header: %d genomes\n", n_genomes);
    }
    fflush(stdout);
    
    // Helper function to read dynamic lines
    size_t data_buffer_size = 1048576;  // Start with 1MB for data lines
    
    // First pass: count reads and total non-zero entries
    long file_pos = ftell(file);
    int n_reads = 0;
    int total_nnz = 0;
    
    char *line;
    
    while ((line = read_line_dynamic(file, &data_buffer_size)) != NULL) {
        n_reads++;
        
        // Print progress every 100,000 reads during counting phase
        if (n_reads % 100000 == 0) {
            printf("Counting reads: %d processed...\n", n_reads);
            fflush(stdout);
        }
        
        // Count quadruplets (genome_name, nd, md, mb) after read_id and total_count
        char *line_copy2 = strdup(line);
        char *token2 = strtok(line_copy2, "\t\n");  // Skip read_id
        token2 = strtok(NULL, "\t\n");  // Skip total_count
        
        int col = 0;
        while ((token2 = strtok(NULL, "\t\n")) != NULL) {
            if (col % 4 == 0) {  // Every 4th token is a genome name
                total_nnz++;
            }
            col++;
        }
        free(line_copy2);
        free(line);  // Free the dynamic line
    }
    
    printf("Finished counting: %d reads, %d total alignments\n", n_reads, total_nnz);
    printf("Allocating sparse matrix (%.1f MB for alignments)...\n", 
           (total_nnz * (sizeof(int) + sizeof(short) * 4)) / 1024.0 / 1024.0);
    fflush(stdout);
    
    // Reset file position
    fseek(file, file_pos, SEEK_SET);
    
    // Allocate sparse structure
    sparse_alignment_t *sparse = sparse_alignment_alloc(n_reads, n_genomes, total_nnz);
    if (!sparse) {
        printf("ERROR: Failed to allocate sparse matrix - out of memory!\n");
        fclose(file);
        return NULL;
    }
    printf("Sparse matrix allocated successfully.\n");
    fflush(stdout);
    
    // Allocate damage arrays
    if (sparse_alignment_alloc_damage(sparse) != 0) {
        sparse_alignment_free(sparse);
        fclose(file);
        return NULL;
    }
    
    // Store taxonomic groups information
    sparse->has_taxonomic_groups = has_taxonomic_groups;
    sparse->taxonomic_groups = taxonomic_groups;
    sparse->n_taxonomic_groups = n_taxonomic_groups;
    
    // Allocate per-read taxonomic data if needed
    if (has_taxonomic_groups) {
        sparse->read_taxonomic_data = calloc(n_reads, sizeof(read_taxonomic_data_t));
        if (!sparse->read_taxonomic_data) {
            sparse_alignment_free(sparse);
            fclose(file);
            return NULL;
        }
    }
    
    // Build hash table for fast genome name lookup
    // Simple hash table implementation
    typedef struct genome_hash_entry {
        char *name;
        int index;
        struct genome_hash_entry *next;
    } genome_hash_entry_t;
    
    // Use efficient hash table size - twice the number of genomes for good performance
    int hash_size = n_genomes * 2;
    printf("Allocating hash table for %d genomes (size: %d entries = %.1f MB)...\n", 
           n_genomes, hash_size, (hash_size * sizeof(genome_hash_entry_t*)) / 1024.0 / 1024.0);
    fflush(stdout);
    
    genome_hash_entry_t **hash_table = calloc(hash_size, sizeof(genome_hash_entry_t*));
    
    if (!hash_table) {
        printf("ERROR: Failed to allocate hash table - out of memory!\n");
        sparse_alignment_free(sparse);
        fclose(file);
        return NULL;
    }
    printf("Hash table allocated successfully.\n");
    printf("Populating hash table with %d genome names...\n", n_genomes);
    fflush(stdout);
    
    // Populate hash table
    for (int i = 0; i < n_genomes; i++) {
        // Progress every 100,000 genomes during hash table population
        if (i % 100000 == 0 && i > 0) {
            printf("Hash table: processed %d/%d genomes...\n", i, n_genomes);
            fflush(stdout);
        }
        unsigned int hash = 0;
        for (char *p = (*genome_names)[i]; *p; p++) {
            hash = hash * 31 + *p;
        }
        hash %= hash_size;
        
        genome_hash_entry_t *entry = malloc(sizeof(genome_hash_entry_t));
        entry->name = (*genome_names)[i];
        entry->index = i;
        entry->next = hash_table[hash];
        hash_table[hash] = entry;
    }
    
    // Build taxonomic group hash table if needed
    typedef struct taxonomic_hash_entry {
        char *name;
        int index;
        struct taxonomic_hash_entry *next;
    } taxonomic_hash_entry_t;
    
    taxonomic_hash_entry_t **taxonomic_hash_table = NULL;
    int taxonomic_hash_size = 0;
    
    if (has_taxonomic_groups && n_taxonomic_groups > 0) {
        taxonomic_hash_size = n_taxonomic_groups * 2;
        printf("Creating taxonomic hash table for %d groups (size: %d entries)...\n", 
               n_taxonomic_groups, taxonomic_hash_size);
        fflush(stdout);
        
        taxonomic_hash_table = calloc(taxonomic_hash_size, sizeof(taxonomic_hash_entry_t*));
        if (!taxonomic_hash_table) {
            printf("ERROR: Failed to allocate taxonomic hash table!\n");
            // Clean up and continue without taxonomic groups
            for (int i = 0; i < hash_size; i++) {
                genome_hash_entry_t *entry = hash_table[i];
                while (entry) {
                    genome_hash_entry_t *next = entry->next;
                    free(entry);
                    entry = next;
                }
            }
            free(hash_table);
            sparse_alignment_free(sparse);
            fclose(file);
            return NULL;
        }
        
        // Populate taxonomic hash table
        for (int i = 0; i < n_taxonomic_groups; i++) {
            unsigned int hash = 0;
            for (char *p = taxonomic_groups[i].name; *p; p++) {
                hash = hash * 31 + *p;
            }
            hash %= taxonomic_hash_size;
            
            
            taxonomic_hash_entry_t *entry = malloc(sizeof(taxonomic_hash_entry_t));
            entry->name = taxonomic_groups[i].name;
            entry->index = i;
            entry->next = taxonomic_hash_table[hash];
            taxonomic_hash_table[hash] = entry;
        }
        printf("Taxonomic hash table populated successfully.\n");
        fflush(stdout);
    }
    
    printf("Hash tables populated successfully. Starting second pass to read data...\n");
    fflush(stdout);
    
    // Second pass: fill sparse structure
    int read_idx = 0;
    int current_nnz = 0;
    sparse->row_ptr[0] = 0;
    
    while ((line = read_line_dynamic(file, &data_buffer_size)) != NULL && read_idx < n_reads) {
        // Skip empty lines (which shouldn't be there but might be due to file format issues)
        if (strlen(line) == 0 || line[0] == '\0') {
            free(line);
            continue;  // Skip this line but don't increment read_idx
        }

        // Print progress every 100,000 reads during data parsing (less frequent for cleaner output)
        if ((read_idx + 1) % 100000 == 0) {
            printf("Read data for %d reads...\n", read_idx + 1);
            fflush(stdout);
        }

        char *line_copy3 = strdup(line);
        if (!line_copy3) {
            fprintf(stderr, "ERROR: Memory allocation failed for line %d (strdup failed)\n", read_idx + 1);
            fprintf(stderr, "  Line length: %zu bytes\n", strlen(line));
            free(line);
            sparse_alignment_free(sparse);
            fclose(file);
            return NULL;
        }

        char *token3 = strtok(line_copy3, "\t\n");  // Skip read_id
        if (!token3) {
            fprintf(stderr, "ERROR: Failed to parse read_id at data line %d\n", read_idx + 1);
            fprintf(stderr, "  This might be caused by empty lines in the file.\n");
            fprintf(stderr, "  Please check that your file doesn't have blank lines between data rows.\n");
            free(line_copy3);
            free(line);
            sparse_alignment_free(sparse);
            fclose(file);
            return NULL;
        }

        token3 = strtok(NULL, "\t\n");  // Get total_count
        if (!token3) {
            fprintf(stderr, "ERROR: Failed to parse total_count at data line %d\n", read_idx + 1);
            free(line_copy3);
            free(line);
            sparse_alignment_free(sparse);
            fclose(file);
            return NULL;
        }
        short total_count = (short)atoi(token3);
        
        // Process genome quadruplets (before "|" separator)
        while ((token3 = strtok(NULL, "\t\n")) != NULL) {
            // Check if we've hit the taxonomic groups separator
            if (strcmp(token3, "|") == 0) {
                break;  // Stop processing genomes, start processing taxonomic groups
            }
            
            char *genome_name = token3;
            
            // Get nd, md, mb values
            token3 = strtok(NULL, "\t\n");
            if (!token3) break;
            short nd_val = (short)atoi(token3);
            
            token3 = strtok(NULL, "\t\n");
            if (!token3) break;
            short md_val = (short)atoi(token3);
            
            token3 = strtok(NULL, "\t\n");
            if (!token3) break;
            short mb_val = (short)atoi(token3);
            
            // Find genome index using hash table
            unsigned int hash = 0;
            for (char *p = genome_name; *p; p++) {
                hash = hash * 31 + *p;
            }
            hash %= hash_size;
            
            int genome_idx = -1;
            genome_hash_entry_t *entry = hash_table[hash];
            while (entry) {
                if (strcmp(genome_name, entry->name) == 0) {
                    genome_idx = entry->index;
                    break;
                }
                entry = entry->next;
            }
            
            if (genome_idx >= 0 && current_nnz < total_nnz) {
                sparse->col_indices[current_nnz] = genome_idx;
                sparse->n_values[current_nnz] = total_count;
                sparse->d_values[current_nnz] = md_val + mb_val;  // Total mismatches
                sparse->nd_values[current_nnz] = nd_val;
                sparse->md_values[current_nnz] = md_val;
                sparse->mb_values[current_nnz] = mb_val;
                current_nnz++;
            }
        }
        
        // Process taxonomic groups if present (after "|" separator)
        if (has_taxonomic_groups && sparse->read_taxonomic_data) {
            // Add occasional progress for taxonomic group processing
            if ((read_idx + 1) % 50000 == 0) {
                printf("  Processing taxonomic groups for read %d...\n", read_idx + 1);
                fflush(stdout);
            }
            // Initialize taxonomic data for this read
            read_taxonomic_data_t *tax_data = &sparse->read_taxonomic_data[read_idx];
            
            // Count taxonomic groups in this read first
            int tax_count = 0;
            char *remaining_tokens[10000];  // Temporary storage for counting
            
            while ((token3 = strtok(NULL, "\t\n")) != NULL) {
                remaining_tokens[tax_count] = strdup(token3);
                tax_count++;
                if (tax_count >= 10000) break;  // Safety limit
            }
            
            // Detect format by analyzing first taxonomic entry pattern
            // Possibilities:
            // Old sparse_dmg: name nd md mds mb mbs (6 values)
            // New sparse_dmg: name nd md mb (4 values, mb is float)
            // Old sparse_std: name mb mbs (3 values)
            // New sparse_std: name mb (2 values, mb is float)
            
            int values_per_entry = 0;
            int n_tax_entries = 0;
            int use_float_format = 0;
            int is_damage_format = 0;
            
            // Auto-detect format by checking first taxonomic entry
            if (tax_count >= 2) {
                // Check new formats first (higher priority), then old formats
                // Priority: new damage (4) > new standard (2) > old damage (6) > old standard (3)
                
                if (tax_count % 4 == 0) {
                    values_per_entry = 4;
                    n_tax_entries = tax_count / 4;
                    use_float_format = 1;
                    is_damage_format = 1;
                } else if (tax_count % 2 == 0) {
                    values_per_entry = 2;
                    n_tax_entries = tax_count / 2;
                    use_float_format = 1;
                    is_damage_format = 0;
                } else if (tax_count % 6 == 0) {
                    values_per_entry = 6;
                    n_tax_entries = tax_count / 6;
                    use_float_format = 0;
                    is_damage_format = 1;
                } else if (tax_count % 3 == 0) {
                    values_per_entry = 3;
                    n_tax_entries = tax_count / 3;
                    use_float_format = 0;
                    is_damage_format = 0;
                } else {
                    printf("Warning: Cannot determine taxonomic format for %d tokens\n", tax_count);
                    n_tax_entries = 0;
                }
            }
            
            // Set format flag on first read (assuming consistent format across file)
            if (read_idx == 0 && n_tax_entries > 0) {
                sparse->use_float_mb_format = use_float_format;
                if (use_float_format) {
                    printf("Detected new float mb format for taxonomic groups\n");
                } else {
                    printf("Detected traditional integer format for taxonomic groups\n");
                }
            }
            
            if (n_tax_entries > 0) {
                tax_data->n_taxa = n_tax_entries;
                tax_data->taxon_indices = malloc(n_tax_entries * sizeof(int));
                
                // Allocate arrays based on detected format
                if (is_damage_format) {
                    tax_data->nd_values = malloc(n_tax_entries * sizeof(short));
                    tax_data->md_values = malloc(n_tax_entries * sizeof(short));
                    if (!use_float_format) {
                        // Old format: nd md mds mb mbs
                        tax_data->mds_values = malloc(n_tax_entries * sizeof(short));
                        tax_data->mb_values = malloc(n_tax_entries * sizeof(short));
                        tax_data->mbs_values = malloc(n_tax_entries * sizeof(short));
                    } else {
                        // New format: nd md mb (mb is float)
                        tax_data->mds_values = NULL;
                        tax_data->mb_values = NULL;
                        tax_data->mbs_values = NULL;
                        tax_data->mb_values_float = malloc(n_tax_entries * sizeof(float));
                        tax_data->mbs_values_float = NULL;
                    }
                } else {
                    // Standard format
                    tax_data->nd_values = NULL;
                    tax_data->md_values = NULL;
                    tax_data->mds_values = NULL;
                    
                    if (!use_float_format) {
                        // Old format: mb mbs
                        tax_data->mb_values = malloc(n_tax_entries * sizeof(short));
                        tax_data->mbs_values = malloc(n_tax_entries * sizeof(short));
                        tax_data->total_all_values = malloc(n_tax_entries * sizeof(short));
                        tax_data->total_some_values = malloc(n_tax_entries * sizeof(short));
                    } else {
                        // New format: mb (mb is float)
                        tax_data->mb_values = NULL;
                        tax_data->mbs_values = NULL;
                        tax_data->total_all_values = NULL;
                        tax_data->total_some_values = NULL;
                        tax_data->mb_values_float = malloc(n_tax_entries * sizeof(float));
                        tax_data->total_all_values_float = malloc(n_tax_entries * sizeof(float));
                    }
                }
                
                // Parse entries according to detected format
                for (int t = 0; t < n_tax_entries; t++) {
                    int base_idx = t * values_per_entry;
                    
                    // Find taxonomic group index by name using hash table (O(1) lookup)
                    char *tax_name = remaining_tokens[base_idx];
                    int tax_idx = -1;
                    
                    if (taxonomic_hash_table) {
                        unsigned int hash = 0;
                        for (char *p = tax_name; *p; p++) {
                            hash = hash * 31 + *p;
                        }
                        hash %= taxonomic_hash_size;
                        
                        taxonomic_hash_entry_t *entry = taxonomic_hash_table[hash];
                        while (entry) {
                            if (strcmp(tax_name, entry->name) == 0) {
                                tax_idx = entry->index;
                                break;
                            }
                            entry = entry->next;
                        }
                    }
                    
                    tax_data->taxon_indices[t] = tax_idx;
                    
                    if (is_damage_format) {
                        // Parse damage format
                        tax_data->nd_values[t] = (short)atoi(remaining_tokens[base_idx + 1]);
                        tax_data->md_values[t] = (short)atoi(remaining_tokens[base_idx + 2]);
                        
                        if (!use_float_format) {
                            // Old damage format: nd md mds mb mbs
                            tax_data->mds_values[t] = (short)atoi(remaining_tokens[base_idx + 3]);
                            tax_data->mb_values[t] = (short)atoi(remaining_tokens[base_idx + 4]);
                            tax_data->mbs_values[t] = (short)atoi(remaining_tokens[base_idx + 5]);
                        } else {
                            // New damage format: nd md mb (mb is float)
                            tax_data->mb_values_float[t] = (float)atof(remaining_tokens[base_idx + 3]);
                        }
                    } else {
                        // Parse standard format
                        if (!use_float_format) {
                            // Old standard format: mb mbs
                            tax_data->mb_values[t] = (short)atoi(remaining_tokens[base_idx + 1]);
                            tax_data->mbs_values[t] = (short)atoi(remaining_tokens[base_idx + 2]);
                            // For compatibility
                            tax_data->total_all_values[t] = tax_data->mb_values[t];
                            tax_data->total_some_values[t] = tax_data->mbs_values[t];
                        } else {
                            // New standard format: mb (mb is float)
                            tax_data->mb_values_float[t] = (float)atof(remaining_tokens[base_idx + 1]);
                            tax_data->total_all_values_float[t] = tax_data->mb_values_float[t];
                        }
                    }
                }
            }
            
            // Free temporary tokens
            for (int i = 0; i < tax_count; i++) {
                free(remaining_tokens[i]);
            }
        }
        
        free(line_copy3);
        free(line);  // Free the dynamic line
        read_idx++;
        sparse->row_ptr[read_idx] = current_nnz;
    }
    
    sparse->nnz = current_nnz;

    // Check if we read the expected number of reads
    if (read_idx != n_reads) {
        fprintf(stderr, "\nERROR: Line count mismatch detected!\n");
        fprintf(stderr, "  Expected %d data lines but only processed %d\n", n_reads, read_idx);
        fprintf(stderr, "  This often happens when the file contains empty lines between data rows.\n");
        fprintf(stderr, "  Please remove any blank lines from your input file.\n");
        fprintf(stderr, "  You can use: grep -v '^$' input.txt > cleaned.txt\n");
        sparse_alignment_free(sparse);
        fclose(file);
        return NULL;
    }

    // Debug: Print parsed sparse matrix (commented out)
    /*
    printf("Reads: %d, Genomes: %d, Non-zeros: %d\n", n_reads, n_genomes, current_nnz);
    for (int i = 0; i < n_genomes; i++) {
        printf("Genome %d: %s\n", i, (*genome_names)[i]);
    }
    for (int i = 0; i < sparse->n_reads && i < 5; i++) {
        printf("Read %d: ", i);
        int start = sparse->row_ptr[i];
        int end = sparse->row_ptr[i + 1];
        for (int idx = start; idx < end; idx++) {
            int j = sparse->col_indices[idx];
            printf("[%s: n=%.0f,nd=%.0f,md=%.0f,mb=%.0f] ", 
                   (*genome_names)[j],
                   sparse->n_values[idx], 
                   sparse->nd_values[idx],
                   sparse->md_values[idx],
                   sparse->mb_values[idx]);
        }
        printf("\n");
    }
    */
    
    // Clean up hash table
    for (int i = 0; i < hash_size; i++) {
        genome_hash_entry_t *entry = hash_table[i];
        while (entry) {
            genome_hash_entry_t *next = entry->next;
            free(entry);
            entry = next;
        }
    }
    free(hash_table);
    
    // Clean up taxonomic hash table
    if (taxonomic_hash_table) {
        for (int i = 0; i < taxonomic_hash_size; i++) {
            taxonomic_hash_entry_t *entry = taxonomic_hash_table[i];
            while (entry) {
                taxonomic_hash_entry_t *next = entry->next;
                free(entry);
                entry = next;
            }
        }
        free(taxonomic_hash_table);
    }
    
    // Print completion message
    if (has_taxonomic_groups) {
        printf("Completed parsing: %d reads, %d genomes, %d taxonomic groups (%d alignments)\n", 
               n_reads, n_genomes, n_taxonomic_groups, current_nnz);
    } else {
        printf("Completed parsing: %d reads, %d genomes (%d alignments)\n", 
               n_reads, n_genomes, current_nnz);
    }
    fflush(stdout);
    
    fclose(file);
    return sparse;
}

// Parse sparse standard format file (alternating genome names and counts)
sparse_alignment_t* parse_sparse_standard_file(const char *filename, char ***genome_names) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Cannot open file %s\n", filename);
        return NULL;
    }
    
    char line[8192];
    
    // Parse header to get genome names
    if (!fgets(line, sizeof(line), file)) {
        fprintf(stderr, "Error: Cannot read header from %s\n", filename);
        fclose(file);
        return NULL;
    }
    
    // Count genome names and taxonomic groups in header (skip read_id and total_count)
    char *line_copy = strdup(line);
    char *token = strtok(line_copy, "\t\n");
    int n_genomes = 0;
    int n_taxonomic_groups = 0;
    int found_separator = 0;
    
    // Skip read_id and total_count
    token = strtok(NULL, "\t\n");  // Skip read_id
    token = strtok(NULL, "\t\n");  // Skip total_count
    
    // Count genomes until "|" separator, then count taxonomic groups
    while (token) {
        if (strcmp(token, "|") == 0) {
            found_separator = 1;
        } else if (!found_separator) {
            n_genomes++;
        } else {
            n_taxonomic_groups++;
        }
        token = strtok(NULL, "\t\n");
    }
    free(line_copy);
    
    if (n_genomes <= 0) {
        fprintf(stderr, "Error: No genomes found in header\n");
        fclose(file);
        return NULL;
    }
    
    // Print detected counts
    if (found_separator) {
        printf("Read header: %d genomes and %d taxonomic groups\n", n_genomes, n_taxonomic_groups);
    } else {
        printf("Read header: %d genomes\n", n_genomes);
    }
    
    // Allocate genome names array
    *genome_names = malloc(n_genomes * sizeof(char*));
    if (!*genome_names) {
        fclose(file);
        return NULL;
    }
    
    // Parse genome names from header (stop at "|" separator)
    line_copy = strdup(line);
    token = strtok(line_copy, "\t\n");
    token = strtok(NULL, "\t\n");  // Skip read_id
    token = strtok(NULL, "\t\n");  // Skip total_count
    
    int genome_idx = 0;
    while (token && genome_idx < n_genomes) {
        if (strcmp(token, "|") == 0) {
            break;  // Stop at taxonomic group separator
        }
        (*genome_names)[genome_idx] = strdup(token);
        genome_idx++;
        token = strtok(NULL, "\t\n");
    }
    free(line_copy);
    
    // Parse taxonomic groups if present
    taxonomic_group_t *taxonomic_groups = NULL;
    if (found_separator && n_taxonomic_groups > 0) {
        // Parse taxonomic group names from header
        line_copy = strdup(line);
        token = strtok(line_copy, "\t\n");
        token = strtok(NULL, "\t\n");  // Skip read_id
        token = strtok(NULL, "\t\n");  // Skip total_count
        
        // Skip genome names to get to "|" separator
        int skipped_genomes = 0;
        while (token && skipped_genomes < n_genomes) {
            token = strtok(NULL, "\t\n");
            skipped_genomes++;
        }
        
        // Now should be at "|" separator
        if (token && strcmp(token, "|") == 0) {
            // Parse taxonomic group names
            taxonomic_groups = malloc(n_taxonomic_groups * sizeof(taxonomic_group_t));
            if (taxonomic_groups) {
                int tax_idx = 0;
                while ((token = strtok(NULL, "\t\n")) != NULL && tax_idx < n_taxonomic_groups) {
                    // Parse T{taxid}:{rank} format (same logic as damage parser)
                    char *taxon_str = token;
                    if (taxon_str[0] == 'T') {
                        char *colon_pos = strchr(taxon_str, ':');
                        if (colon_pos) {
                            taxonomic_groups[tax_idx].name = strdup(taxon_str);  // Full name (T1:genus)
                            
                            // Extract taxid and rank
                            *colon_pos = '\0';  // Split string at colon
                            taxonomic_groups[tax_idx].taxid = atoi(taxon_str + 1);  // Skip 'T'
                            taxonomic_groups[tax_idx].rank = strdup(colon_pos + 1);
                            *colon_pos = ':';  // Restore original string
                        } else {
                            // Malformed taxonomic group
                            taxonomic_groups[tax_idx].taxid = -1;
                            taxonomic_groups[tax_idx].rank = strdup("unknown");
                            taxonomic_groups[tax_idx].name = strdup(taxon_str);
                        }
                    } else {
                        // Not a valid taxonomic group
                        taxonomic_groups[tax_idx].taxid = -1;
                        taxonomic_groups[tax_idx].rank = strdup("unknown");
                        taxonomic_groups[tax_idx].name = strdup(taxon_str);
                    }
                    tax_idx++;
                }
            }
        }
        free(line_copy);
    }
    
    
    // First pass: count reads and non-zeros
    long file_pos = ftell(file);
    int n_reads = 0;
    int total_nnz = 0;
    
    while (fgets(line, sizeof(line), file)) {
        n_reads++;
        
        // Count genome entries in this line (stop at "|" separator for taxonomic groups)
        char temp_line[8192];
        strcpy(temp_line, line);
        char *count_token = strtok(temp_line, "\t\n");  // read_id
        count_token = strtok(NULL, "\t\n");  // total_count
        count_token = strtok(NULL, "\t\n");  // First genome or NULL
        
        // Count alternating genome names and values (only before "|" separator)
        int field_count = 0;
        while (count_token) {
            if (strcmp(count_token, "|") == 0) {
                break;  // Stop at taxonomic group separator
            }
            if (field_count % 2 == 0) {  // Genome name
                total_nnz++;
            }
            field_count++;
            count_token = strtok(NULL, "\t\n");
        }
    }
    
    
    // Allocate sparse structure
    sparse_alignment_t *sparse = sparse_alignment_alloc(n_reads, n_genomes, total_nnz);
    if (!sparse) {
        for (int i = 0; i < n_genomes; i++) {
            free((*genome_names)[i]);
        }
        free(*genome_names);
        fclose(file);
        return NULL;
    }
    
    // Allocate taxonomic data arrays before second pass
    if (found_separator && n_taxonomic_groups > 0) {
        printf("Allocating taxonomic data for %d reads...\n", n_reads);
        sparse->read_taxonomic_data = calloc(n_reads, sizeof(read_taxonomic_data_t));
        if (!sparse->read_taxonomic_data) {
            sparse_alignment_free(sparse);
            for (int i = 0; i < n_genomes; i++) {
                free((*genome_names)[i]);
            }
            free(*genome_names);
            fclose(file);
            return NULL;
        }
    } else {
        sparse->read_taxonomic_data = NULL;
    }
    
    // Reset file position
    fseek(file, file_pos, SEEK_SET);
    
    // Second pass: fill sparse structure
    int read_idx = 0;
    int current_nnz = 0;
    sparse->row_ptr[0] = 0;
    
    while (fgets(line, sizeof(line), file) && read_idx < n_reads) {
        char temp_line[8192];
        strcpy(temp_line, line);
        
        strtok(temp_line, "\t\n");  // Skip read_id
        char *total_count_str = strtok(NULL, "\t\n");  // total_count
        double total_count = total_count_str ? atof(total_count_str) : 0.0;
        
        // if (read_idx < 5) {
        // }
        
        // Parse alternating genome names and mismatch counts (stop at "|" for taxonomic groups)
        char *genome_name = strtok(NULL, "\t\n");
        while (genome_name) {
            // Check for taxonomic group separator
            if (strcmp(genome_name, "|") == 0) {
                // Parse taxonomic groups for this read (standard format: copy from damage format logic)
                if (found_separator && sparse->read_taxonomic_data) {
                    read_taxonomic_data_t *tax_data = &sparse->read_taxonomic_data[read_idx];
                    
                    // Collect remaining tokens after "|" (same as damage format)
                    char *remaining_tokens[10000];
                    int tax_count = 0;
                    
                    char *token3 = strtok(NULL, "\t\n");
                    while (token3 && tax_count < 10000) {
                        remaining_tokens[tax_count] = strdup(token3);
                        tax_count++;
                        token3 = strtok(NULL, "\t\n");
                    }
                    
                    // Parse taxonomic group entries (each entry has: name + mb + mbs = 3 values)
                    int n_tax_entries = tax_count / 3;  // Each taxonomic entry has 3 values (name + 2 numbers)
                    
                    if (n_tax_entries > 0) {
                        tax_data->n_taxa = n_tax_entries;
                        tax_data->taxon_indices = malloc(n_tax_entries * sizeof(int));
                        
                        // For standard format: allocate same arrays as damage format
                        tax_data->nd_values = malloc(n_tax_entries * sizeof(short));   // Total sites (equivalent to damage format)
                        tax_data->mb_values = malloc(n_tax_entries * sizeof(short));   // Background mismatches (all)
                        tax_data->mbs_values = malloc(n_tax_entries * sizeof(short));  // Background mismatches (some)
                        
                        // Set damage-only arrays to NULL for standard format  
                        tax_data->md_values = NULL;
                        tax_data->mds_values = NULL;
                        
                        // For compatibility with existing code, also set total_all_values/total_some_values
                        tax_data->total_all_values = malloc(n_tax_entries * sizeof(short));
                        tax_data->total_some_values = malloc(n_tax_entries * sizeof(short));
                        
                        if (tax_data->taxon_indices && tax_data->nd_values && tax_data->mb_values && tax_data->mbs_values) {
                            
                            for (int t = 0; t < n_tax_entries; t++) {
                                int base_idx = t * 3;
                                
                                // Find taxonomic group index by name (copy exact logic from damage format)
                                char *tax_name = remaining_tokens[base_idx];
                                int tax_idx = -1;
                                
                                // Use hash table lookup if available (same as damage format)
                                // Note: Would need to implement taxonomic hash table for standard format
                                // For now, use linear search
                                for (int i = 0; i < n_taxonomic_groups; i++) {
                                    if (strcmp(taxonomic_groups[i].name, tax_name) == 0) {
                                        tax_idx = i;
                                        break;
                                    }
                                }
                                
                                tax_data->taxon_indices[t] = tax_idx;
                                tax_data->nd_values[t] = (short)total_count;  // Total sites (equivalent to damage format)
                                tax_data->mb_values[t] = (short)atoi(remaining_tokens[base_idx + 1]);   // Background mismatches (all)
                                tax_data->mbs_values[t] = (short)atoi(remaining_tokens[base_idx + 2]);  // Background mismatches (some)
                                
                                
                                // For compatibility, also populate total arrays
                                if (tax_data->total_all_values && tax_data->total_some_values) {
                                    tax_data->total_all_values[t] = (short)total_count;  // Total sites
                                    tax_data->total_some_values[t] = tax_data->mb_values[t] + tax_data->mbs_values[t];  // Total mismatches
                                }
                            }
                        }
                    }
                    
                    // Free temporary tokens
                    for (int i = 0; i < tax_count; i++) {
                        free(remaining_tokens[i]);
                    }
                }
                break;
            }
            
            char *count_str = strtok(NULL, "\t\n");
            if (count_str) {
                int mismatch_count = atoi(count_str);
                
                // Find genome index
                int found_idx = -1;
                for (int i = 0; i < n_genomes; i++) {
                    if (strcmp((*genome_names)[i], genome_name) == 0) {
                        found_idx = i;
                        break;
                    }
                }
                
                if (found_idx >= 0 && current_nnz < total_nnz) {
                    sparse->col_indices[current_nnz] = found_idx;
                    sparse->n_values[current_nnz] = total_count;
                    sparse->d_values[current_nnz] = mismatch_count;
                    current_nnz++;
                }
            }
            genome_name = strtok(NULL, "\t\n");
        }
        
        read_idx++;
        sparse->row_ptr[read_idx] = current_nnz;
    }
    
    sparse->nnz = current_nnz;
    
    // Set up taxonomic group data in sparse structure
    sparse->has_taxonomic_groups = found_separator;
    sparse->taxonomic_groups = taxonomic_groups;
    sparse->n_taxonomic_groups = found_separator ? n_taxonomic_groups : 0;
    
    // Debug output
    /*
    printf("DEBUG parse_sparse_standard_file: n_reads=%d, n_genomes=%d, nnz=%d\n", 
           sparse->n_reads, sparse->n_genomes, sparse->nnz);
    for (int i = 0; i < n_genomes; i++) {
        printf("  Genome %d: %s\n", i, (*genome_names)[i]);
    }
    printf("First 5 entries: ");
    for (int i = 0; i < 5 && i < current_nnz; i++) {
        printf("[col=%d,n=%.0f,d=%.0f] ", 
               sparse->col_indices[i], sparse->n_values[i], sparse->d_values[i]);
    }
    printf("\n");
    */
    
    fclose(file);
    return sparse;
}

// Parse sparse format file (OLD FORMAT - kept for backward compatibility)
sparse_alignment_t* parse_sparse_file(const char *filename, char ***genome_names) {
    // Check if this is actually sparse standard format
    FILE *test_file = fopen(filename, "r");
    if (!test_file) {
        fprintf(stderr, "Error: Cannot open file %s\n", filename);
        return NULL;
    }
    
    // Use dynamic line reading for header
    size_t buffer_size = 0;
    char *test_line = read_line_dynamic(test_file, &buffer_size);
    if (!test_line) {
        fclose(test_file);
        return NULL;
    }
    free(test_line);  // Don't need header for format detection
    
    // Check first data line
    test_line = read_line_dynamic(test_file, &buffer_size);
    if (test_line) {
        // Count tabs to determine format
        int tab_count = 0;
        for (char *p = test_line; *p; p++) {
            if (*p == '\t') tab_count++;
        }
        
        // If we have many fields and no colons, it's likely sparse standard format
        if (tab_count > 3 && !strchr(test_line, ':')) {
            free(test_line);
            fclose(test_file);
            return parse_sparse_standard_file(filename, genome_names);
        }
        free(test_line);
    }
    fclose(test_file);
    
    // Otherwise use old parser (for backward compatibility)
    // But with dynamic allocation for large headers
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Cannot open file %s\n", filename);
        return NULL;
    }
    
    // Skip header line using dynamic reader
    buffer_size = 0;
    char *header_line = read_line_dynamic(file, &buffer_size);
    if (!header_line) {
        fprintf(stderr, "Error: Cannot read header from %s\n", filename);
        fclose(file);
        return NULL;
    }
    free(header_line);  // Free header after reading
    
    // First pass: collect all unique genome names and count reads/nnz
    long file_pos = ftell(file);
    
    // Use dynamic array for genome names
    int genome_capacity = 10000;
    char **unique_genomes = malloc(genome_capacity * sizeof(char*));
    if (!unique_genomes) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        fclose(file);
        return NULL;
    }
    int n_genomes = 0;
    int n_reads = 0;
    int total_nnz = 0;
    
    char line[MAX_LINE_BUFFER];  // Buffer for data lines
    
    while (fgets(line, sizeof(line), file)) {
        n_reads++;
        
        // Find the aligned_genomes field (last tab-separated field)
        char line_copy[4096];
        strcpy(line_copy, line);
        
        char *token = strtok(line_copy, "\t\n");
        char *last_field = NULL;
        
        while (token) {
            last_field = token;
            token = strtok(NULL, "\t\n");
        }
        
        if (last_field) {
            // Parse aligned_genomes: "genome2:2,genome5:1,genome7:3"
            char field_copy[1024];
            strcpy(field_copy, last_field);
            
            char *entry = strtok(field_copy, ",");
            while (entry) {
                total_nnz++;
                
                char *colon = strchr(entry, ':');
                if (colon) {
                    *colon = '\0';
                    char *genome_name = entry;
                    
                    // Check if this genome is already in our list
                    int found = 0;
                    for (int i = 0; i < n_genomes; i++) {
                        if (strcmp(unique_genomes[i], genome_name) == 0) {
                            found = 1;
                            break;
                        }
                    }
                    
                    if (!found && n_genomes < 1000) {
                        strcpy(unique_genomes[n_genomes], genome_name);
                        n_genomes++;
                    }
                }
                entry = strtok(NULL, ",");
            }
        }
    }
    
    // Allocate genome names array
    *genome_names = malloc(n_genomes * sizeof(char*));
    for (int i = 0; i < n_genomes; i++) {
        (*genome_names)[i] = malloc(strlen(unique_genomes[i]) + 1);
        strcpy((*genome_names)[i], unique_genomes[i]);
    }
    
    // Reset file position and skip header
    fseek(file, file_pos, SEEK_SET);
    
    // Allocate sparse structure
    sparse_alignment_t *sparse = sparse_alignment_alloc(n_reads, n_genomes, total_nnz);
    if (!sparse) {
        fclose(file);
        return NULL;
    }
    
    // Second pass: fill sparse structure
    int read_idx = 0;
    int current_nnz = 0;
    sparse->row_ptr[0] = 0;
    
    while (fgets(line, sizeof(line), file) && read_idx < n_reads) {
        char *fields[4];
        int field_count = 0;
        
        // Split line by tabs
        char *token = strtok(line, "\t\n");
        while (token && field_count < 4) {
            fields[field_count++] = token;
            token = strtok(NULL, "\t\n");
        }
        
        if (field_count >= 3) {
            double total_count = atof(fields[1]);
            char *aligned_genomes = fields[2];
            
            // Parse aligned_genomes field: "genome2:2,genome5:1,genome7:3"
            char aligned_copy[1024];
            strcpy(aligned_copy, aligned_genomes);
            
            char *entry = strtok(aligned_copy, ",");
            while (entry) {
                char *colon = strchr(entry, ':');
                if (colon) {
                    *colon = '\0';
                    char *genome_name = entry;
                    int mismatch_count = atoi(colon + 1);
                    
                    // Find genome index
                    int genome_idx = -1;
                    for (int i = 0; i < n_genomes; i++) {
                        if (strcmp(genome_name, (*genome_names)[i]) == 0) {
                            genome_idx = i;
                            break;
                        }
                    }
                    
                    if (genome_idx >= 0 && current_nnz < total_nnz) {
                        sparse->col_indices[current_nnz] = genome_idx;
                        sparse->n_values[current_nnz] = total_count;
                        sparse->d_values[current_nnz] = mismatch_count;
                        current_nnz++;
                    }
                }
                entry = strtok(NULL, ",");
            }
        }
        
        read_idx++;
        sparse->row_ptr[read_idx] = current_nnz;
    }
    
    sparse->nnz = current_nnz;
    fclose(file);
    return sparse;
}

// Parse dense format file and convert to sparse
sparse_alignment_t* parse_dense_file_to_sparse(const char *filename, char ***genome_names) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Cannot open file %s\n", filename);
        return NULL;
    }
    
    char line[MAX_LINE_BUFFER];
    
    // Read header line to get genome names
    if (!fgets(line, sizeof(line), file)) {
        fprintf(stderr, "Error: Cannot read header from %s\n", filename);
        fclose(file);
        return NULL;
    }
    
    // Check if line was truncated
    size_t line_len = strlen(line);
    if (line_len > 0 && line[line_len - 1] != '\n') {
        fprintf(stderr, "Error: Header line exceeds maximum buffer size of %d characters.\n", MAX_LINE_BUFFER);
        fprintf(stderr, "This typically happens with many genomes (>100).\n");
        fprintf(stderr, "Please use --short-names option in BAMreader to reduce header size.\n");
        fclose(file);
        return NULL;
    }
    
    // Parse header to extract genome names (skip read_id and total_count)
    char *token = strtok(line, "\t\n");
    token = strtok(NULL, "\t\n");  // Skip total_count
    
    // Count genomes and store names
    int n_genomes = 0;
    char temp_names[1000][256];  // Temporary storage
    
    while ((token = strtok(NULL, "\t\n")) != NULL) {
        strcpy(temp_names[n_genomes], token);
        n_genomes++;
    }
    
    // Allocate genome names array
    *genome_names = malloc(n_genomes * sizeof(char*));
    for (int i = 0; i < n_genomes; i++) {
        (*genome_names)[i] = malloc(strlen(temp_names[i]) + 1);
        strcpy((*genome_names)[i], temp_names[i]);
    }
    
    // Count reads and non-zero entries
    long file_pos = ftell(file);
    int n_reads = 0;
    int total_nnz = 0;
    
    while (fgets(line, sizeof(line), file)) {
        n_reads++;
        
        // Count non-negative values (skip read_id and total_count)
        char *token = strtok(line, "\t\n");
        token = strtok(NULL, "\t\n");  // Skip total_count
        
        for (int j = 0; j < n_genomes; j++) {
            token = strtok(NULL, "\t\n");
            if (token && atof(token) >= 0) {
                total_nnz++;
            }
        }
    }
    
    // Reset file position
    fseek(file, file_pos, SEEK_SET);
    
    // Allocate sparse structure
    sparse_alignment_t *sparse = sparse_alignment_alloc(n_reads, n_genomes, total_nnz);
    if (!sparse) {
        fclose(file);
        return NULL;
    }
    
    // Parse data lines
    int read_idx = 0;
    int current_nnz = 0;
    sparse->row_ptr[0] = 0;
    
    while (fgets(line, sizeof(line), file) && read_idx < n_reads) {
        char *token = strtok(line, "\t\n");  // Skip read_id
        token = strtok(NULL, "\t\n");        // Get total_count
        double total_count = atof(token);
        
        // Process genome columns
        for (int j = 0; j < n_genomes; j++) {
            token = strtok(NULL, "\t\n");
            if (token) {
                double value = atof(token);
                if (value >= 0) {  // Valid alignment (not -1)
                    sparse->col_indices[current_nnz] = j;
                    sparse->n_values[current_nnz] = total_count;
                    sparse->d_values[current_nnz] = value;
                    current_nnz++;
                }
            }
        }
        
        read_idx++;
        sparse->row_ptr[read_idx] = current_nnz;
    }
    
    sparse->nnz = current_nnz;
    fclose(file);
    return sparse;
}

// Debug: Print sparse matrix structure
void sparse_matrix_print(const sparse_alignment_t *sparse) {
    printf("Sparse Matrix: %d reads x %d genomes, %d non-zeros\n", 
           sparse->n_reads, sparse->n_genomes, sparse->nnz);
    
    for (int i = 0; i < sparse->n_reads && i < 10; i++) {  // Limit output
        printf("Read %d: ", i);
        int start = sparse->row_ptr[i];
        int end = sparse->row_ptr[i + 1];
        
        for (int idx = start; idx < end; idx++) {
            printf("(%d:%hd/%hd) ", sparse->col_indices[idx], 
                   sparse->n_values[idx], sparse->d_values[idx]);
        }
        printf("\n");
    }
}

// Validate sparse matrix structure
int sparse_matrix_validate(const sparse_alignment_t *sparse) {
    if (!sparse) return 0;
    
    // Check basic structure
    if (sparse->n_reads <= 0 || sparse->n_genomes <= 0 || sparse->nnz < 0) {
        return 0;
    }
    
    // Check row_ptr monotonicity
    for (int i = 0; i < sparse->n_reads; i++) {
        if (sparse->row_ptr[i] > sparse->row_ptr[i + 1]) {
            return 0;
        }
    }
    
    // Check column indices bounds
    for (int idx = 0; idx < sparse->nnz; idx++) {
        if (sparse->col_indices[idx] < 0 || 
            sparse->col_indices[idx] >= sparse->n_genomes) {
            return 0;
        }
    }
    
    // Check final row_ptr value
    if (sparse->row_ptr[sparse->n_reads] != sparse->nnz) {
        return 0;
    }
    
    return 1;  // Valid
}

// Print statistics about sparse matrix
void sparse_matrix_stats(const sparse_alignment_t *sparse) {
    if (!sparse) {
        printf("Sparse matrix is NULL\n");
        return;
    }
    
    double total_dense = (double)sparse->n_reads * sparse->n_genomes;
    double sparsity = 1.0 - (double)sparse->nnz / total_dense;
    
    printf("=== Sparse Matrix Statistics ===\n");
    printf("Dimensions: %d reads x %d genomes\n", sparse->n_reads, sparse->n_genomes);
    printf("Non-zero entries: %d / %.0f (%.2f%% sparse)\n", 
           sparse->nnz, total_dense, sparsity * 100.0);
    printf("Memory usage: %.2f MB (vs %.2f MB dense)\n",
           (sparse->nnz * 2 * sizeof(double) + sparse->nnz * sizeof(int) + 
            (sparse->n_reads + 1) * sizeof(int)) / (1024.0 * 1024.0),
           (total_dense * 2 * sizeof(double)) / (1024.0 * 1024.0));
    
    // Calculate average alignments per read
    double avg_alignments = (double)sparse->nnz / sparse->n_reads;
    printf("Average alignments per read: %.2f\n", avg_alignments);
}