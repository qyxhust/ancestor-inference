#include "io_utils.h"
#include "memory.h"
#include "em_algorithms.h"
#include <ctype.h>
#include <stdarg.h>

// File I/O utilities implementation

// Helper to get method name as string
static const char* get_method_string(em_method_t method) {
    switch(method) {
        case METHOD_A: return "A";
        case METHOD_AE: return "AE";
        case METHOD_B: return "B";
        case METHOD_BE: return "BE";
        case METHOD_BEFULL: return "BEfull";
        case METHOD_C: return "C";
        case METHOD_CE: return "CE";
        case METHOD_CED: return "CED";
        case METHOD_CEFULL: return "CEfull";
        case METHOD_CEDFULL: return "CEDfull";
        // METHOD_FULL removed - non-functional extended model
        case METHOD_INVALID: return "Invalid";
        default: return "Unknown";
    }
}

// Global pointer for qsort comparison in this file
static const double *io_qsort_proportions = NULL;

// Comparison function for sorting genome indices by proportion (descending)
static int io_compare_indices_desc(const void *a, const void *b) {
    int idx_a = *(const int *)a;
    int idx_b = *(const int *)b;
    
    // Sort in descending order by proportion
    if (io_qsort_proportions[idx_b] > io_qsort_proportions[idx_a]) return 1;
    if (io_qsort_proportions[idx_b] < io_qsort_proportions[idx_a]) return -1;
    return 0;
}

// Helper function to create sorted index array by proportion (descending)
static int* create_proportion_sorted_index(const double *proportions, int n_genomes) {
    // Create array of indices
    int *indices = malloc(n_genomes * sizeof(int));
    if (!indices) return NULL;
    
    for (int i = 0; i < n_genomes; i++) {
        indices[i] = i;
    }
    
    // Use quicksort for O(n log n) sorting instead of O(n²) selection sort
    io_qsort_proportions = proportions;
    qsort(indices, n_genomes, sizeof(int), io_compare_indices_desc);
    io_qsort_proportions = NULL;  // Clear for safety
    
    return indices;
}

// Write standardized output file
int write_standardized_output(const char *filename, em_results_t *results, char **genome_names,
                             int n_genomes, em_config_t *config) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        print_error("Cannot create output file: %s\n", filename);
        return -1;
    }
    
    // Write header section with metadata
    fprintf(file, "Model=%s\n", get_method_string(config->method));
    fprintf(file, "LogLikelihood=%.16f\n", results->final_log_likelihood);
    fprintf(file, "SQUAREM=%s\n", config->use_squarem ? "Yes" : "No");
    fprintf(file, "Iterations=%d\n", results->iterations);
    fprintf(file, "Converged=%s\n", results->converged ? "Yes" : "No");
    fprintf(file, "ComputationTime=%.2f\n", results->computation_time);
    fprintf(file, "\n");
    
    // Write column headers
    // Standard header (METHOD_FULL some_error_rate column removed)
    fprintf(file, "genome\tproportion\terror_rate\tdamage_rate\n");
    
    // Create sorted index array (highest proportion first)
    int *sorted_indices = create_proportion_sorted_index(results->final_proportions, n_genomes);
    if (!sorted_indices) {
        print_error("Failed to allocate memory for sorting\n");
        fclose(file);
        return -1;
    }
    
    // Write data for each genome in sorted order (only non-pruned genomes)
    int output_count = 0;
    double min_output_proportion = 1e-10;  // Only output genomes with proportion > 1e-10
    
    for (int i = 0; i < n_genomes; i++) {
        int genome_idx = sorted_indices[i];  // Use sorted index
        
        // Skip pruned genomes (proportion near zero)
        if (results->final_proportions[genome_idx] < min_output_proportion) {
            continue;
        }
        
        output_count++;
        const char *name = genome_names[genome_idx] ? genome_names[genome_idx] : "Unknown";
        
        // Translate genome name if mapping is available
        if (config->genome_mappings && config->n_mappings > 0) {
            char *translated_name = translate_genome_name(name, config->genome_mappings, config->n_mappings);
            if (translated_name) {
                name = translated_name;
            }
        }
        
        fprintf(file, "%s\t%.16f\t", name, results->final_proportions[genome_idx]);
        
        // Write error rate column
        if (config->method == METHOD_A || config->method == METHOD_B || config->method == METHOD_C) {
            // Fixed error rate models
            fprintf(file, "fixed=%.16f", config->error_rate);
        } else if (config->method == METHOD_AE || config->method == METHOD_BE || 
                   config->method == METHOD_CE || config->method == METHOD_CED) {
            // Single estimated error rate for all genomes
            fprintf(file, "%.16f", results->final_error_rate);
        } else if (config->method == METHOD_BEFULL || config->method == METHOD_CEFULL || 
                   config->method == METHOD_CEDFULL) {
            // Per-genome error rates
            if (results->final_error_rates) {
                fprintf(file, "%.16f", results->final_error_rates[genome_idx]);
            } else {
                fprintf(file, "NA");
            }
        } else {
            fprintf(file, "NA");
        }
        
        fprintf(file, "\t");
        
        // Write damage rate column
        if (config->method == METHOD_A || config->method == METHOD_AE || 
            config->method == METHOD_B || config->method == METHOD_BE || 
            config->method == METHOD_BEFULL) {
            // No damage models
            fprintf(file, "not_applicable");
        } else if (config->method == METHOD_C) {
            // Fixed damage rate
            fprintf(file, "fixed=%.16f", config->damage_rate);
        } else if (config->method == METHOD_CE || config->method == METHOD_CEFULL) {
            // Fixed damage rate (given as parameter)
            fprintf(file, "fixed=%.16f", config->damage_rate);
        } else if (config->method == METHOD_CED) {
            // Single estimated damage rate for all genomes
            fprintf(file, "%.16f", results->final_damage_rate);
        } else if (config->method == METHOD_CEDFULL) {
            // Per-genome damage rates
            if (results->final_damage_rates) {
                fprintf(file, "%.16f", results->final_damage_rates[genome_idx]);
            } else {
                fprintf(file, "NA");
            }
        } else {
            fprintf(file, "NA");
        }
        
        fprintf(file, "\n");
    }
    
    // Write taxonomic group results if present - WITH SORTING
    int taxonomic_output_count = 0;
    if (results->final_taxonomic_proportions && results->taxonomic_groups && results->n_taxonomic_groups > 0) {
        // Create sorted index array for taxonomic groups (highest proportion first)
        int *sorted_tax_indices = create_proportion_sorted_index(results->final_taxonomic_proportions, results->n_taxonomic_groups);
        if (!sorted_tax_indices) {
            print_error("Failed to allocate memory for sorting taxonomic groups\n");
            // Continue without sorting rather than failing completely
            sorted_tax_indices = malloc(results->n_taxonomic_groups * sizeof(int));
            if (sorted_tax_indices) {
                for (int i = 0; i < results->n_taxonomic_groups; i++) {
                    sorted_tax_indices[i] = i;
                }
            }
        }

        if (sorted_tax_indices) {
            for (int i = 0; i < results->n_taxonomic_groups; i++) {
                int tax_idx = sorted_tax_indices[i];  // Use sorted index

                // Skip taxonomic groups with very small proportions
                if (results->final_taxonomic_proportions[tax_idx] < min_output_proportion) {
                    continue;
                }

                taxonomic_output_count++;
                const char *name = results->taxonomic_groups[tax_idx].name ? results->taxonomic_groups[tax_idx].name : "Unknown";

                fprintf(file, "%s\t%.16f\t", name, results->final_taxonomic_proportions[tax_idx]);

                // For taxonomic groups, show the same fixed rates used in the analysis
                // (they use the same global rates as individual genomes in the extended models)
                // Output taxonomic rates for CEDfull (restored for unified approach)
                if (config->method == METHOD_CEDFULL && results->final_taxonomic_error_rates && results->final_taxonomic_damage_rates) {
                    // CEDfull - show per-taxonomic group error and damage rates
                    fprintf(file, "%.16f\t%.16f\n", results->final_taxonomic_error_rates[tax_idx], results->final_taxonomic_damage_rates[tax_idx]);
                } else if (config->method == METHOD_A || config->method == METHOD_B || config->method == METHOD_C) {
                    // Fixed rate methods - show the fixed rates used
                    fprintf(file, "fixed=%.16f\tfixed=%.16f\n", results->final_error_rate, results->final_damage_rate);
                } else if (config->method == METHOD_AE || config->method == METHOD_BE || config->method == METHOD_CE) {
                    // Estimated error rate, fixed damage rate
                    fprintf(file, "%.16f\tfixed=%.16f\n", results->final_error_rate, results->final_damage_rate);
                } else if (config->method == METHOD_CED) {
                    // Both rates estimated
                    fprintf(file, "%.16f\t%.16f\n", results->final_error_rate, results->final_damage_rate);
                } else {
                    // Other per-genome rate methods - taxonomic groups use global rates
                    fprintf(file, "not_applicable\tnot_applicable\n");
                }
            }

            // Free the sorted taxonomic indices array
            free(sorted_tax_indices);
        }
    }
    
    if (taxonomic_output_count > 0) {
        printf("Output written: %d genomes + %d taxonomic groups (filtered from %d + %d total, threshold: %.0e)\n", 
               output_count, taxonomic_output_count, n_genomes, results->n_taxonomic_groups, min_output_proportion);
    } else {
        printf("Output written: %d genomes (filtered from %d total, threshold: %.0e)\n", 
               output_count, n_genomes, min_output_proportion);
    }
    
    // Free the sorted indices array
    free(sorted_indices);
    
    fclose(file);
    return 0;
}

// Write results to files
int write_results(const char *output_prefix, em_results_t *results, char **genome_names, 
                 int n_genomes, em_config_t *config) {
    if (!results || !config) return -1;
    
    // Create output filename using the provided output_prefix
    char output_filename[512];
    if (output_prefix) {
        // Use the user-specified output prefix (from -o option)
        snprintf(output_filename, sizeof(output_filename), "%s_out.txt", output_prefix);
    } else if (config->input_file) {
        // Fallback: Remove directory path and extension from input file
        const char *base = strrchr(config->input_file, '/');
        base = base ? base + 1 : config->input_file;
        
        char base_copy[256];
        strncpy(base_copy, base, sizeof(base_copy) - 1);
        base_copy[sizeof(base_copy) - 1] = '\0';
        
        // Remove extension
        char *dot = strrchr(base_copy, '.');
        if (dot) *dot = '\0';
        
        snprintf(output_filename, sizeof(output_filename), "%s_out.txt", base_copy);
    } else {
        strcpy(output_filename, "output_out.txt");
    }
    
    // Write single standardized output file
    if (write_standardized_output(output_filename, results, genome_names, n_genomes, config) != 0) {
        print_error("Failed to write output file: %s\n", output_filename);
        return -1;
    }
    
    printf("Results written to: %s\n", output_filename);
    return 0;
}

// Write proportions CSV file
int write_proportions_csv(const char *filename, em_results_t *results, char **genome_names, int n_genomes) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        print_error("Cannot create output file: %s\n", filename);
        return -1;
    }
    
    // Write header - include per-genome rates if available
    if (results->final_error_rates || results->final_damage_rates) {
        // Full models - include parameter columns
        fprintf(file, "Genome,Proportion");
        if (results->final_error_rates) {
            fprintf(file, ",Error_Rate");
        }
        if (results->final_damage_rates) {
            fprintf(file, ",Damage_Rate");
        }
        fprintf(file, "\n");
        
        // Write data with per-genome parameters
        for (int i = 0; i < n_genomes; i++) {
            const char *name = (genome_names && genome_names[i]) ? genome_names[i] : "Unknown";
            fprintf(file, "%s,%.16f", name, results->final_proportions[i]);
            if (results->final_error_rates) {
                fprintf(file, ",%.16f", results->final_error_rates[i]);
            }
            if (results->final_damage_rates) {
                fprintf(file, ",%.16f", results->final_damage_rates[i]);
            }
            fprintf(file, "\n");
        }
    } else {
        // Simple models - just proportions
        fprintf(file, "Genome,Proportion\n");
        for (int i = 0; i < n_genomes; i++) {
            const char *name = (genome_names && genome_names[i]) ? genome_names[i] : "Unknown";
            fprintf(file, "%s,%.16f\n", name, results->final_proportions[i]);
        }
    }
    
    fclose(file);
    return 0;
}

// Write rates CSV file  
int write_rates_csv(const char *filename, em_results_t *results, char **genome_names, int n_genomes, em_config_t *config) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        print_error("Cannot create output file: %s\n", filename);
        return -1;
    }
    
    // Write header based on method
    switch (config->method) {
        case METHOD_AE:
        case METHOD_BE:
            fprintf(file, "Parameter,Value\n");
            fprintf(file, "Error_Rate,%.16f\n", results->final_error_rate);
            break;
            
        case METHOD_CE:
            fprintf(file, "Parameter,Value\n");
            fprintf(file, "Background_Rate,%.16f\n", results->final_error_rate);
            fprintf(file, "Damage_Rate,%.16f\n", config->damage_rate);
            break;
            
        case METHOD_CED:
            fprintf(file, "Parameter,Value\n");
            fprintf(file, "Background_Rate,%.16f\n", results->final_error_rate);
            fprintf(file, "Damage_Rate,%.16f\n", results->final_damage_rate);
            break;
            
        case METHOD_BEFULL:
            fprintf(file, "Genome,Error_Rate\n");
            for (int i = 0; i < n_genomes; i++) {
                const char *name = (genome_names && genome_names[i]) ? genome_names[i] : "Unknown";
                fprintf(file, "%s,%.16f\n", name, results->final_error_rates[i]);
            }
            break;
            
        case METHOD_CEFULL:
            fprintf(file, "Genome,Error_Rate,Damage_Rate\n");
            for (int i = 0; i < n_genomes; i++) {
                const char *name = (genome_names && genome_names[i]) ? genome_names[i] : "Unknown";
                fprintf(file, "%s,%.16f,%.16f\n", name, results->final_error_rates[i], config->damage_rate);
            }
            break;
            
        case METHOD_CEDFULL:
            fprintf(file, "Genome,Error_Rate,Damage_Rate\n");
            for (int i = 0; i < n_genomes; i++) {
                const char *name = (genome_names && genome_names[i]) ? genome_names[i] : "Unknown";
                fprintf(file, "%s,%.16f,%.16f\n", name, results->final_error_rates[i], results->final_damage_rates[i]);
            }
            break;
            
        default:
            break;
    }
    
    fclose(file);
    return 0;
}

// Error handling functions
void print_error(const char *format, ...) {
    va_list args;
    va_start(args, format);
    fprintf(stderr, "Error: ");
    vfprintf(stderr, format, args);
    va_end(args);
}

void print_warning(const char *format, ...) {
    va_list args;
    va_start(args, format);
    fprintf(stderr, "Warning: ");
    vfprintf(stderr, format, args);
    va_end(args);
}

void print_verbose(int verbose, const char *format, ...) {
    if (!verbose) return;
    
    va_list args;
    va_start(args, format);
    vprintf(format, args);
    va_end(args);
    fflush(stdout);
}

// String utilities
char* string_duplicate(const char *src) {
    if (!src) return NULL;
    
    size_t len = strlen(src);
    char *dst = malloc(len + 1);
    if (!dst) return NULL;
    
    strcpy(dst, src);
    return dst;
}

char** string_array_create(int size) {
    if (size <= 0) return NULL;
    
    char **array = calloc(size, sizeof(char*));
    return array;
}

void string_array_destroy(char **array, int size) {
    if (!array) return;
    
    for (int i = 0; i < size; i++) {
        free(array[i]);
    }
    free(array);
}

// Trim whitespace from string (in-place)
char* trim_whitespace(char *str) {
    if (!str) return NULL;
    
    // Trim leading whitespace
    while (isspace(*str)) str++;
    
    // Empty string
    if (*str == '\0') return str;
    
    // Trim trailing whitespace
    char *end = str + strlen(str) - 1;
    while (end > str && isspace(*end)) end--;
    
    // Write new null terminator
    *(end + 1) = '\0';
    
    return str;
}

// Count columns in TSV line
int count_columns(const char *line, char delimiter) {
    if (!line) return 0;
    
    int count = 1;  // At least one column
    for (const char *p = line; *p; p++) {
        if (*p == delimiter) count++;
    }
    
    return count;
}

// Check if string is numeric
int is_numeric_string(const char *str) {
    if (!str || !*str) return 0;
    
    // Skip leading whitespace
    while (isspace(*str)) str++;
    
    // Check for sign
    if (*str == '+' || *str == '-') str++;
    
    // Must have at least one digit
    int has_digit = 0;
    int has_dot = 0;
    
    while (*str) {
        if (isdigit(*str)) {
            has_digit = 1;
        } else if (*str == '.' && !has_dot) {
            has_dot = 1;
        } else if (!isspace(*str)) {
            return 0;
        } else {
            // Trailing whitespace is OK
            while (*str && isspace(*str)) str++;
            return *str == '\0' && has_digit;
        }
        str++;
    }
    
    return has_digit;
}

// Read mismatch matrix from file
int read_mismatch_matrix(const char *filename, em_data_t **data, char ***genome_names) {
    
    // First detect format using the new file-based detection
    input_format_t format = detect_input_format_from_file(filename);
    
    FILE *file = fopen(filename, "r");
    if (!file) {
        print_error("Cannot open input file: %s\n", filename);
        return -1;
    }
    
    // Use dynamic buffer for potentially large headers
    size_t line_buffer_size = 1048576;  // Start with 1MB
    char *line = malloc(line_buffer_size);
    if (!line) {
        print_error("Memory allocation failed\n");
        fclose(file);
        return -1;
    }
    
    int n_genomes = 0;
    int n_reads = 0;
    
    // Read header with dynamic growth if needed
    size_t header_len = 0;
    int c;
    while ((c = fgetc(file)) != EOF && c != '\n') {
        if (header_len >= line_buffer_size - 1) {
            // Need to grow buffer
            line_buffer_size *= 2;
            char *new_buffer = realloc(line, line_buffer_size);
            if (!new_buffer) {
                print_error("Memory allocation failed while reading header\n");
                free(line);
                fclose(file);
                return -1;
            }
            line = new_buffer;
        }
        line[header_len++] = c;
    }
    line[header_len] = '\0';
    
    if (header_len == 0) {
        print_error("Empty file or cannot read header\n");
        free(line);
        fclose(file);
        return -1;
    }
    
    // Handle different formats
    if (format == FORMAT_DENSE_STANDARD) {
        // Parse header to get genome names
        if (parse_tsv_header(line, genome_names, &n_genomes) != 0) {
            print_error("Failed to parse header\n");
            fclose(file);
            return -1;
        }
        
        // Count reads
        int line_count = 0;
        while (fgets(line, line_buffer_size, file)) {
            trim_whitespace(line);
            if (strlen(line) > 0) line_count++;
        }
        n_reads = line_count;
        
        // Allocate data structure
        *data = em_data_create(n_reads, n_genomes);
        if (!*data) {
            print_error("Failed to allocate data structure\n");
            string_array_destroy(*genome_names, n_genomes);
            fclose(file);
            return -1;
        }
        
        // Rewind and skip header
        rewind(file);
        fgets(line, line_buffer_size, file);
        
        // Read data
        int read_idx = 0;
        while (fgets(line, line_buffer_size, file) && read_idx < n_reads) {
            int total_count;
            int *mismatches = malloc(n_genomes * sizeof(int));
            
            if (parse_tsv_data_line(line, &total_count, &mismatches, n_genomes) != 0) {
                print_error("Failed to parse data line %d\n", read_idx + 1);
                free(mismatches);
                em_data_destroy(*data);
                string_array_destroy(*genome_names, n_genomes);
                fclose(file);
                return -1;
            }
            
            // Store data
            for (int j = 0; j < n_genomes; j++) {
                (*data)->n_matrix[read_idx][j] = (short)total_count;
                (*data)->d_matrix[read_idx][j] = (short)mismatches[j];
            }
            
            free(mismatches);
            read_idx++;
        }
        
    } else if (format == FORMAT_DENSE_DAMAGE) {
        // Parse damage format header
        if (parse_damage_tsv_header(line, genome_names, &n_genomes) != 0) {
            print_error("Failed to parse damage format header\n");
            fclose(file);
            return -1;
        }
        
        // Count reads
        int line_count = 0;
        while (fgets(line, line_buffer_size, file)) {
            trim_whitespace(line);
            if (strlen(line) > 0) line_count++;
        }
        n_reads = line_count;
        
        // Allocate data structure with damage matrices
        *data = em_data_create(n_reads, n_genomes);
        if (!*data) {
            print_error("Failed to allocate data structure\n");
            string_array_destroy(*genome_names, n_genomes);
            fclose(file);
            return -1;
        }
        
        // Allocate damage matrices
        if (em_data_allocate_damage_matrices(*data) != 0) {
            print_error("Failed to allocate damage matrices\n");
            em_data_destroy(*data);
            string_array_destroy(*genome_names, n_genomes);
            fclose(file);
            return -1;
        }
        
        // Rewind and skip header
        rewind(file);
        fgets(line, line_buffer_size, file);
        
        // Read data
        int read_idx = 0;
        while (fgets(line, line_buffer_size, file) && read_idx < n_reads) {
            int total_count;
            int *nd_values = malloc(n_genomes * sizeof(int));
            int *md_values = malloc(n_genomes * sizeof(int));
            int *mb_values = malloc(n_genomes * sizeof(int));
            
            if (parse_damage_tsv_data_line(line, &total_count, &nd_values, &md_values, &mb_values, n_genomes) != 0) {
                print_error("Failed to parse damage data line %d\n", read_idx + 1);
                free(nd_values);
                free(md_values);
                free(mb_values);
                em_data_destroy(*data);
                string_array_destroy(*genome_names, n_genomes);
                fclose(file);
                return -1;
            }
            
            // Store data (casting to short for memory efficiency)
            for (int j = 0; j < n_genomes; j++) {
                (*data)->n_matrix[read_idx][j] = (short)total_count;
                (*data)->nd_matrix[read_idx][j] = (short)nd_values[j];
                (*data)->md_matrix[read_idx][j] = (short)md_values[j];
                (*data)->mb_matrix[read_idx][j] = (short)mb_values[j];
                (*data)->d_matrix[read_idx][j] = (short)(md_values[j] + mb_values[j]);  // Total mismatches
            }
            
            free(nd_values);
            free(md_values);
            free(mb_values);
            read_idx++;
        }
    } else if (format == FORMAT_SPARSE_DAMAGE) {
        // Parse sparse damage format header
        if (parse_tsv_header(line, genome_names, &n_genomes) != 0) {
            print_error("Failed to parse sparse damage format header\n");
            fclose(file);
            return -1;
        }
        
        // First pass: count reads
        int line_count = 0;
        while (fgets(line, line_buffer_size, file)) {
            trim_whitespace(line);
            if (strlen(line) > 0) line_count++;
        }
        n_reads = line_count;
        
        // Allocate data structure with damage matrices
        *data = em_data_create(n_reads, n_genomes);
        if (!*data) {
            print_error("Failed to allocate data structure\n");
            string_array_destroy(*genome_names, n_genomes);
            fclose(file);
            return -1;
        }
        
        // Allocate damage matrices (using short for memory efficiency)
        (*data)->nd_matrix = matrix_create_short(n_reads, n_genomes);
        (*data)->md_matrix = matrix_create_short(n_reads, n_genomes);
        (*data)->mb_matrix = matrix_create_short(n_reads, n_genomes);
        
        if (!(*data)->nd_matrix || !(*data)->md_matrix || !(*data)->mb_matrix) {
            print_error("Failed to allocate damage matrices\n");
            em_data_destroy(*data);
            string_array_destroy(*genome_names, n_genomes);
            fclose(file);
            return -1;
        }
        
        // Initialize all matrices to -1 (missing data)
        for (int i = 0; i < n_reads; i++) {
            for (int j = 0; j < n_genomes; j++) {
                (*data)->n_matrix[i][j] = -1;
                (*data)->d_matrix[i][j] = -1;
                (*data)->nd_matrix[i][j] = -1;
                (*data)->md_matrix[i][j] = -1;
                (*data)->mb_matrix[i][j] = -1;
            }
        }
        
        // Rewind and skip header
        rewind(file);
        fgets(line, line_buffer_size, file);
        
        // Read sparse damage data
        int read_idx = 0;
        while (fgets(line, line_buffer_size, file) && read_idx < n_reads) {
            
            // Use the parse_sparse_damage_tsv_data_line function
            int total_count;
            char **alignment_genome_names = NULL;
            int *nd_values = NULL;
            int *md_values = NULL;
            int *mb_values = NULL;
            int n_alignments = 0;
            
            if (parse_sparse_damage_tsv_data_line(line, &total_count, &alignment_genome_names, &nd_values, &md_values, &mb_values, &n_alignments, *genome_names, n_genomes) != 0) {
                print_error("Failed to parse sparse damage data line %d\n", read_idx + 1);
                em_data_destroy(*data);
                string_array_destroy(*genome_names, n_genomes);
                fclose(file);
                return -1;
            }
            
            // Fill in actual alignment data
            for (int a = 0; a < n_alignments; a++) {
                // Find genome index
                for (int j = 0; j < n_genomes; j++) {
                    if (strcmp(alignment_genome_names[a], (*genome_names)[j]) == 0) {
                        if (nd_values[a] == -1 || md_values[a] == -1 || mb_values[a] == -1) {
                            // Explicitly marked as missing - keep -1 values
                        } else {
                            // Valid alignment - use actual values
                            (*data)->n_matrix[read_idx][j] = (short)total_count;
                            (*data)->nd_matrix[read_idx][j] = (short)nd_values[a];
                            (*data)->md_matrix[read_idx][j] = (short)md_values[a];
                            (*data)->mb_matrix[read_idx][j] = (short)mb_values[a];
                            (*data)->d_matrix[read_idx][j] = (short)(md_values[a] + mb_values[a]);
                        }
                        break;
                    }
                }
            }
            
            // Clean up
            string_array_destroy(alignment_genome_names, n_alignments);
            free(nd_values);
            free(md_values);
            free(mb_values);
            read_idx++;
        }
    } else {
        print_error("Unsupported format detected\n");
        free(line);
        fclose(file);
        return -1;
    }
    
    free(line);  // Free the dynamically allocated buffer
    fclose(file);
    return 0;
}

// Parse sparse damage TSV data line
int parse_sparse_damage_tsv_data_line(const char *line, int *total_count, char ***genome_names,
                                      int **nd_values, int **md_values, int **mb_values, int *n_alignments,
                                      char **all_genome_names, int n_total_genomes) {
    if (!line || !total_count || !genome_names || !nd_values || !md_values || !mb_values || !n_alignments) return -1;
    
    // Mark unused parameters to avoid warnings
    (void)all_genome_names;
    (void)n_total_genomes;
    
    char *line_copy = string_duplicate(line);
    if (!line_copy) return -1;
    
    char *token = strtok(line_copy, "\t\n\r");
    int col = 0;
    int alignment_count = 0;
    
    // Temporary storage for this read's alignments
    char **temp_genome_names = NULL;
    int *temp_nd_values = NULL;
    int *temp_md_values = NULL;
    int *temp_mb_values = NULL;
    
    while (token != NULL) {
        if (col == 1) {  // total_count column
            *total_count = atoi(trim_whitespace(token));
            if (*total_count <= 0) {
                free(line_copy);
                return -1;
            }
        } else if (col >= 2 && (col - 2) % 4 == 0) {  // genome name columns (every 4th position after first 2)
            char *genome_name = trim_whitespace(token);
            
            // Get the corresponding nd, md, mb values (next 3 tokens)
            token = strtok(NULL, "\t\n\r");
            if (!token) break;
            int nd_val = atoi(trim_whitespace(token));
            
            token = strtok(NULL, "\t\n\r");
            if (!token) break;
            int md_val = atoi(trim_whitespace(token));
            
            token = strtok(NULL, "\t\n\r");
            if (!token) break;
            int mb_val = atoi(trim_whitespace(token));
            
            if (nd_val < -1 || md_val < -1 || mb_val < -1) {  // Allow -1 (missing alignment)
                string_array_destroy(temp_genome_names, alignment_count);
                free(temp_nd_values);
                free(temp_md_values);
                free(temp_mb_values);
                free(line_copy);
                return -1;
            }
            
            // Store this alignment
            temp_genome_names = realloc(temp_genome_names, (alignment_count + 1) * sizeof(char*));
            temp_nd_values = realloc(temp_nd_values, (alignment_count + 1) * sizeof(int));
            temp_md_values = realloc(temp_md_values, (alignment_count + 1) * sizeof(int));
            temp_mb_values = realloc(temp_mb_values, (alignment_count + 1) * sizeof(int));
            
            if (!temp_genome_names || !temp_nd_values || !temp_md_values || !temp_mb_values) {
                string_array_destroy(temp_genome_names, alignment_count);
                free(temp_nd_values);
                free(temp_md_values);
                free(temp_mb_values);
                free(line_copy);
                return -1;
            }
            
            temp_genome_names[alignment_count] = string_duplicate(genome_name);
            if (!temp_genome_names[alignment_count]) {
                string_array_destroy(temp_genome_names, alignment_count);
                free(temp_nd_values);
                free(temp_md_values);
                free(temp_mb_values);
                free(line_copy);
                return -1;
            }
            
            temp_nd_values[alignment_count] = nd_val;
            temp_md_values[alignment_count] = md_val;
            temp_mb_values[alignment_count] = mb_val;
            alignment_count++;
            
            col += 3; // Skip the values we just processed
        }
        token = strtok(NULL, "\t\n\r");
        col++;
    }
    
    *genome_names = temp_genome_names;
    *nd_values = temp_nd_values;
    *md_values = temp_md_values;
    *mb_values = temp_mb_values;
    *n_alignments = alignment_count;
    
    free(line_copy);
    return 0;
}

// Parse TSV header line to extract genome names
int parse_tsv_header(const char *line, char ***genome_names, int *n_genomes) {
    if (!line || !genome_names || !n_genomes) return -1;
    
    char *line_copy = string_duplicate(line);
    if (!line_copy) return -1;
    
    // Count columns (skip read_id and total_count)
    int total_cols = count_columns(line, '\t');
    *n_genomes = total_cols - 2;  // Subtract read_id and total_count columns
    
    if (*n_genomes <= 0) {
        free(line_copy);
        return -1;
    }
    
    if (*n_genomes > 5000) {
        print_error("Number of genomes (%d) exceeds maximum limit of 5000\n", *n_genomes);
        free(line_copy);
        return -1;
    }
    
    // Allocate genome names array
    *genome_names = string_array_create(*n_genomes);
    if (!*genome_names) {
        free(line_copy);
        return -1;
    }
    
    // Parse genome names
    char *token = strtok(line_copy, "\t\n\r");
    int col = 0;
    int genome_idx = 0;
    
    while (token != NULL) {
        if (col >= 2) {  // Skip read_id and total_count
            (*genome_names)[genome_idx] = string_duplicate(trim_whitespace(token));
            if (!(*genome_names)[genome_idx]) {
                string_array_destroy(*genome_names, *n_genomes);
                free(line_copy);
                return -1;
            }
            genome_idx++;
        }
        token = strtok(NULL, "\t\n\r");
        col++;
    }
    
    free(line_copy);
    return 0;
}

// Parse damage TSV header
int parse_damage_tsv_header(const char *line, char ***genome_names, int *n_genomes) {
    if (!line || !genome_names || !n_genomes) return -1;
    
    char *line_copy = string_duplicate(line);
    if (!line_copy) return -1;
    
    // Count genome names by looking for nd_ prefixes
    int genome_count = 0;
    char *temp_copy = string_duplicate(line);
    char *token = strtok(temp_copy, "\t\n\r");
    
    while (token != NULL) {
        char *trimmed = trim_whitespace(token);
        if (strncmp(trimmed, "nd_", 3) == 0) {
            genome_count++;
        }
        token = strtok(NULL, "\t\n\r");
    }
    free(temp_copy);
    
    *n_genomes = genome_count;
    
    if (*n_genomes <= 0) {
        free(line_copy);
        return -1;
    }
    
    if (*n_genomes > 5000) {
        print_error("Number of genomes (%d) exceeds maximum limit of 5000\n", *n_genomes);
        free(line_copy);
        return -1;
    }
    
    // Allocate genome names array
    *genome_names = string_array_create(*n_genomes);
    if (!*genome_names) {
        free(line_copy);
        return -1;
    }
    
    // Parse genome names from nd_ columns
    token = strtok(line_copy, "\t\n\r");
    int genome_idx = 0;
    
    while (token != NULL && genome_idx < *n_genomes) {
        char *trimmed = trim_whitespace(token);
        if (strncmp(trimmed, "nd_", 3) == 0) {
            // Extract genome name (remove nd_ prefix)
            (*genome_names)[genome_idx] = string_duplicate(trimmed + 3);
            if (!(*genome_names)[genome_idx]) {
                string_array_destroy(*genome_names, *n_genomes);
                free(line_copy);
                return -1;
            }
            genome_idx++;
        }
        token = strtok(NULL, "\t\n\r");
    }
    
    free(line_copy);
    return 0;
}

// Parse TSV data line
int parse_tsv_data_line(const char *line, int *total_count, int **mismatches, int n_genomes) {
    if (!line || !total_count || !mismatches || n_genomes <= 0) return -1;
    
    char *line_copy = string_duplicate(line);
    if (!line_copy) return -1;
    
    char *token = strtok(line_copy, "\t\n\r");
    
    // Skip read_id (first column)
    if (!token) {
        free(line_copy);
        return -1;
    }
    
    // Get total_count (second column)
    token = strtok(NULL, "\t\n\r");
    if (!token) {
        free(line_copy);
        return -1;
    }
    *total_count = atoi(token);
    
    // Get mismatches for each genome
    for (int i = 0; i < n_genomes; i++) {
        token = strtok(NULL, "\t\n\r");
        if (!token) {
            free(line_copy);
            return -1;
        }
        (*mismatches)[i] = atoi(token);
    }
    
    free(line_copy);
    return 0;
}

// Parse damage TSV data line
int parse_damage_tsv_data_line(const char *line, int *total_count, int **nd_values, 
                              int **md_values, int **mb_values, int n_genomes) {
    if (!line || !total_count || !nd_values || !md_values || !mb_values || n_genomes <= 0) {
        return -1;
    }
    
    char *line_copy = string_duplicate(line);
    if (!line_copy) return -1;
    
    char *token = strtok(line_copy, "\t\n\r");
    
    // Skip read_id (first column)
    if (!token) {
        free(line_copy);
        return -1;
    }
    
    // Get total_count (second column)
    token = strtok(NULL, "\t\n\r");
    if (!token) {
        free(line_copy);
        return -1;
    }
    *total_count = atoi(token);
    
    // Get damage values for each genome (3 values per genome)
    for (int i = 0; i < n_genomes; i++) {
        // nd value
        token = strtok(NULL, "\t\n\r");
        if (!token) {
            free(line_copy);
            return -1;
        }
        (*nd_values)[i] = atoi(token);
        
        // md value
        token = strtok(NULL, "\t\n\r");
        if (!token) {
            free(line_copy);
            return -1;
        }
        (*md_values)[i] = atoi(token);
        
        // mb value
        token = strtok(NULL, "\t\n\r");
        if (!token) {
            free(line_copy);
            return -1;
        }
        (*mb_values)[i] = atoi(token);
    }
    
    free(line_copy);
    return 0;
}

// Validate input data
int validate_input_data(em_data_t *data) {
    if (!data) return -1;
    
    // Check for valid dimensions
    if (data->n_reads <= 0 || data->n_genomes <= 0) {
        print_error("Invalid data dimensions: %d reads, %d genomes\n", 
                   data->n_reads, data->n_genomes);
        return -1;
    }
    
    // Check for valid matrices
    if (!data->n_matrix || !data->d_matrix) {
        print_error("Missing required data matrices\n");
        return -1;
    }
    
    // Validate mismatch counts
    for (int i = 0; i < data->n_reads; i++) {
        for (int j = 0; j < data->n_genomes; j++) {
            if (data->d_matrix[i][j] > data->n_matrix[i][j]) {
                print_error("Invalid data: mismatches (%g) > sites (%g) for read %d, genome %d\n",
                           data->d_matrix[i][j], data->n_matrix[i][j], i, j);
                return -1;
            }
            
            if (data->d_matrix[i][j] < 0 || data->n_matrix[i][j] < 0) {
                print_error("Invalid data: negative values for read %d, genome %d\n", i, j);
                return -1;
            }
        }
    }
    
    // Validate damage data if present
    if (data->nd_matrix && data->md_matrix && data->mb_matrix) {
        for (int i = 0; i < data->n_reads; i++) {
            for (int j = 0; j < data->n_genomes; j++) {
                if (data->md_matrix[i][j] > data->nd_matrix[i][j]) {
                    print_error("Invalid damage data: md (%g) > nd (%g) for read %d, genome %d\n",
                               data->md_matrix[i][j], data->nd_matrix[i][j], i, j);
                    return -1;
                }
                
                double total_mismatches = data->md_matrix[i][j] + data->mb_matrix[i][j];
                double expected_mismatches = data->d_matrix[i][j];
                
                if (fabs(total_mismatches - expected_mismatches) > 1e-9) {
                    print_error("Invalid damage data: md + mb (%g) != d (%g) for read %d, genome %d\n",
                               total_mismatches, expected_mismatches, i, j);
                    return -1;
                }
            }
        }
    }
    
    return 0;
}

// Detect input format based on header and first data line
input_format_t detect_input_format_from_file(const char *filename) {
    
    FILE *file = fopen(filename, "r");
    if (!file) {
        return FORMAT_DENSE_STANDARD;
    }
    
    // Use dynamic allocation with larger initial size
    size_t buffer_size = 1048576;  // Start with 1MB
    char *header_line = malloc(buffer_size);
    char *data_line = malloc(buffer_size);
    
    if (!header_line || !data_line) {
        if (header_line) free(header_line);
        if (data_line) free(data_line);
        fclose(file);
        return FORMAT_DENSE_STANDARD;
    }
    
    // Read header with dynamic buffer growth if needed
    size_t header_len = 0;
    int c;
    while ((c = fgetc(file)) != EOF && c != '\n') {
        if (header_len >= buffer_size - 1) {
            // Need to grow buffer
            size_t new_size = buffer_size * 2;
            char *new_header = realloc(header_line, new_size);
            char *new_data = realloc(data_line, new_size);
            if (!new_header || !new_data) {
                if (new_header) free(new_header);
                else free(header_line);
                if (new_data) free(new_data);
                else free(data_line);
                fclose(file);
                return FORMAT_DENSE_STANDARD;
            }
            header_line = new_header;
            data_line = new_data;
            buffer_size = new_size;
        }
        header_line[header_len++] = c;
    }
    header_line[header_len] = '\0';
    
    if (header_len == 0) {
        free(header_line);
        free(data_line);
        fclose(file);
        return FORMAT_DENSE_STANDARD;
    }
    
    // Read first data line
    if (!fgets(data_line, buffer_size, file)) {
        free(header_line);
        free(data_line);
        fclose(file);
        return FORMAT_DENSE_STANDARD;
    }

    // Ensure null termination
    data_line[buffer_size - 1] = '\0';

    fclose(file);
    
    // First check header for damage prefixes and taxonomic groups
    int has_damage_prefixes = 0;
    int has_taxonomic_groups = 0;
    char *header_copy = string_duplicate(header_line);
    if (header_copy) {
        char *token = strtok(header_copy, "\t\n\r");
        while (token != NULL) {
            char *trimmed = trim_whitespace(token);
            if (strncmp(trimmed, "nd_", 3) == 0 || 
                strncmp(trimmed, "md_", 3) == 0 || 
                strncmp(trimmed, "mb_", 3) == 0) {
                has_damage_prefixes = 1;
            }
            if (strcmp(trimmed, "|") == 0) {
                has_taxonomic_groups = 1;
            }
            token = strtok(NULL, "\t\n\r");
        }
        free(header_copy);
    }
    
    // Check for error condition: dense format with taxonomic groups
    // Dense format is determined by matching header and data column counts
    // If we have taxonomic groups, we should force sparse format
    if (has_taxonomic_groups) {
        // Count columns to determine if this would be classified as dense
        int header_cols = count_columns(header_line, '\t');
        int data_cols = count_columns(data_line, '\t');
        
        if (header_cols == data_cols || (has_damage_prefixes && data_cols == (header_cols - 2) * 3 + 2)) {
            // This would be classified as dense format, but has taxonomic groups
            print_error("Error: Dense format with higher taxonomic groups is not supported.\n");
            print_error("Please use sparse format input file for taxonomic group analysis.\n");
            print_error("Dense format files cannot contain the '|' separator.\n");
            free(header_line);
            free(data_line);
            return FORMAT_UNSUPPORTED;
        }
    }
    
    // Count columns in header and data
    int header_cols = count_columns(header_line, '\t');
    int data_cols = count_columns(data_line, '\t');

    // Dense formats: 
    // - Standard: header columns = data columns
    // - Damage: data has 3x values compared to genome columns in header
    // Sparse formats: header has all genome names, but data has variable columns
    
    // Count genome columns (header columns minus metadata columns)
    int genome_cols_in_header = header_cols - 2;  // Subtract read_id and total_count
    int cols_after_metadata = data_cols - 2;      // Subtract read_id and total_count
    
    // Files with taxonomic groups are ALWAYS sparse format
    // For files with taxonomic groups, determine damage vs standard by examining the data
    if (has_taxonomic_groups) {
        // Check if the data has damage format (3 values for genomes, 5 for taxonomic groups)
        // Count tokens in first data line after metadata
        char *data_copy = string_duplicate(data_line);
        if (data_copy) {
            char *token = strtok(data_copy, "\t\n\r");
            token = strtok(NULL, "\t\n\r"); // Skip total_count
            
            // Count values between first genome name and next genome name or "|" separator
            int value_count = 0;
            int found_first_genome = 0;
            char *first_genome_name = NULL;
            
            while ((token = strtok(NULL, "\t\n\r")) != NULL) {
                if (strcmp(token, "|") == 0) break;
                
                if (!found_first_genome) {
                    // This should be the first genome name
                    first_genome_name = token;
                    found_first_genome = 1;
                    value_count = 0;
                } else {
                    // Check if this token is a number (indicating values after genome name)
                    char *endptr;
                    strtol(token, &endptr, 10);
                    if (*endptr == '\0') {
                        // This is a number (value after genome)
                        value_count++;
                    } else {
                        // This is another genome name, stop counting
                        break;
                    }
                }
            }
            free(data_copy);
            
            // If we found 3 values per genome, it's damage format
            if (value_count == 3) {
                free(header_line);
                free(data_line);
                return FORMAT_SPARSE_DAMAGE;
            } else {
                free(header_line);
                free(data_line);
                return FORMAT_SPARSE_STANDARD;
            }
        }
        
        // Fallback: assume damage format if we can't determine
        free(header_line);
        free(data_line);
        return FORMAT_SPARSE_DAMAGE;
    }
    
    if (header_cols == data_cols) {
        // Dense standard format - columns match exactly
        if (has_damage_prefixes) {
            free(header_line);
            free(data_line);
            return FORMAT_DENSE_DAMAGE;
        } else {
            free(header_line);
            free(data_line);
            return FORMAT_DENSE_STANDARD;
        }
    } else if (cols_after_metadata == genome_cols_in_header * 3) {
        // Dense damage format - data has 3x values (n, nd, md) per genome
        free(header_line);
        free(data_line);
        return FORMAT_DENSE_DAMAGE;
    } else {
        // Sparse format - check if it's damage format by looking for quadruplets
        
        // In sparse damage format, we have quadruplets: genome_name, nd, md, mb
        // So (cols_after_metadata % 4) should be 0
        if (cols_after_metadata > 0 && cols_after_metadata % 4 == 0) {
            free(header_line);
            free(data_line);
            return FORMAT_SPARSE_DAMAGE;
        } else if (cols_after_metadata > 0 && cols_after_metadata % 2 == 0) {
            // Sparse standard format has pairs: genome_name, mismatch_count
            free(header_line);
            free(data_line);
            return FORMAT_SPARSE_STANDARD;
        }
    }

    free(header_line);
    free(data_line);
    return FORMAT_DENSE_STANDARD;
}

// Legacy function - kept for compatibility but now uses file-based detection
input_format_t detect_input_format(const char *header_line) {
    if (!header_line) return FORMAT_DENSE_STANDARD;
    
    char *line_copy = string_duplicate(header_line);
    if (!line_copy) return FORMAT_DENSE_STANDARD;
    
    // Check for damage format indicators
    int has_nd_prefix = 0;
    int has_md_prefix = 0;
    int has_mb_prefix = 0;
    int damage_columns = 0;
    int total_columns = 0;
    
    char *token = strtok(line_copy, "\t\n\r");
    int col = 0;
    
    while (token != NULL) {
        total_columns++;
        char *trimmed = trim_whitespace(token);
        
        // Skip first two columns (read_id and total_count)
        if (col >= 2) {
            if (strncmp(trimmed, "nd_", 3) == 0) {
                has_nd_prefix = 1;
                damage_columns++;
            } else if (strncmp(trimmed, "md_", 3) == 0) {
                has_md_prefix = 1;
                damage_columns++;
            } else if (strncmp(trimmed, "mb_", 3) == 0) {
                has_mb_prefix = 1;
                damage_columns++;
            }
        }
        token = strtok(NULL, "\t\n\r");
        col++;
    }
    
    free(line_copy);
    
    // If we have damage prefixes, it's a damage format
    if (has_nd_prefix && has_md_prefix && has_mb_prefix) {
        // Check if it's dense (all genomes represented) or sparse
        // Dense format: each genome has exactly 3 columns (nd_, md_, mb_)
        int genome_columns = total_columns - 2;  // Subtract read_id and total_count
        if (genome_columns % 3 == 0 && damage_columns == genome_columns) {
            return FORMAT_DENSE_DAMAGE;
        } else {
            return FORMAT_SPARSE_DAMAGE;
        }
    } else {
        // Standard format (no damage prefixes)
        // Check if it's sparse by looking for alternating pattern
        
        // Sparse standard format characteristics:
        // 1. No fixed header with genome names
        // 2. First column is read_id, second is total_count
        // 3. Remaining columns alternate: genome_name, count, genome_name, count...
        
        // All formats now have headers starting with "read_id"
        // Check if this is a dense format by looking for prefixed columns
        char *line_copy_temp = string_duplicate(header_line);
        char *first_token = strtok(line_copy_temp, "\t\n\r");
        
        if (first_token && (strcmp(first_token, "read_id") == 0)) {
            // This is a header line - check if it's dense or sparse format
            // Dense formats have genome names directly in header
            // Sparse formats also have genome names in header but data is in sparse format
            
            // The key difference: we'll parse as dense first, and let the data parsing
            // determine if it's actually sparse based on the data structure
            free(line_copy_temp);
            return FORMAT_DENSE_STANDARD;  // Default assumption, will be corrected during data parsing
        }
        
        free(line_copy_temp);
        
        // Default to dense standard format
        return FORMAT_DENSE_STANDARD;
    }
}

// Check if genome names are short format (G1, G2, etc.)
int are_genome_names_short(char **genome_names, int n_genomes) {
    if (!genome_names || n_genomes <= 0) return 0;
    
    for (int i = 0; i < n_genomes; i++) {
        if (!genome_names[i]) return 0;
        
        // Check if name matches pattern G[number]
        if (genome_names[i][0] != 'G') return 0;
        
        // Check if rest is a number
        char *endptr;
        long num = strtol(genome_names[i] + 1, &endptr, 10);
        if (*endptr != '\0' || num <= 0) return 0;
    }
    
    return 1;  // All names match G[number] pattern
}

// Auto-detect genome key file based on input filename
char* detect_genome_key_file(const char *input_file) {
    if (!input_file) return NULL;
    
    // Create potential genome key filename by replacing extension with _genome_key.txt
    char *key_file = malloc(strlen(input_file) + 20);
    if (!key_file) return NULL;
    
    strcpy(key_file, input_file);
    
    // Find last dot to remove extension
    char *last_dot = strrchr(key_file, '.');
    if (last_dot && last_dot != key_file) {
        *last_dot = '\0';
    }
    
    // Append _genome_key.txt
    strcat(key_file, "_genome_key.txt");
    
    // Check if file exists
    FILE *test = fopen(key_file, "r");
    if (test) {
        fclose(test);
        return key_file;
    }
    
    free(key_file);
    return NULL;
}

// Load genome mapping from file
genome_mapping_t* load_genome_mapping(const char *filename, int *n_mappings) {
    if (!filename || !n_mappings) return NULL;
    
    FILE *file = fopen(filename, "r");
    if (!file) {
        print_error("Cannot open mapping file: %s\n", filename);
        return NULL;
    }
    
    // Count lines to determine array size
    char line[1024];
    int line_count = 0;
    while (fgets(line, sizeof(line), file)) {
        line_count++;
    }
    
    // Subtract header line
    *n_mappings = line_count - 1;
    if (*n_mappings <= 0) {
        fclose(file);
        return NULL;
    }
    
    // Allocate mapping array
    genome_mapping_t *mappings = calloc(*n_mappings, sizeof(genome_mapping_t));
    if (!mappings) {
        fclose(file);
        return NULL;
    }
    
    // Rewind and skip header
    rewind(file);
    fgets(line, sizeof(line), file);  // Skip header
    
    // Read mappings
    int idx = 0;
    while (fgets(line, sizeof(line), file) && idx < *n_mappings) {
        char *tab = strchr(line, '\t');
        if (!tab) {
            // Invalid format, cleanup and return
            free_genome_mapping(mappings, idx);
            fclose(file);
            return NULL;
        }
        
        // Split at tab
        *tab = '\0';
        char *compressed_id = trim_whitespace(line);
        char *full_name = trim_whitespace(tab + 1);
        
        // Store mapping
        mappings[idx].compressed_id = string_duplicate(compressed_id);
        mappings[idx].full_name = string_duplicate(full_name);
        
        if (!mappings[idx].compressed_id || !mappings[idx].full_name) {
            free_genome_mapping(mappings, idx + 1);
            fclose(file);
            return NULL;
        }
        
        idx++;
    }
    
    fclose(file);
    return mappings;
}

// Free genome mapping array
void free_genome_mapping(genome_mapping_t *mappings, int n_mappings) {
    if (!mappings) return;
    
    for (int i = 0; i < n_mappings; i++) {
        free(mappings[i].compressed_id);
        free(mappings[i].full_name);
    }
    free(mappings);
}

// Translate compressed genome ID to full name
// Simple hash function for genome names
static unsigned int hash_genome_name(const char *str, int bucket_count) {
    unsigned int hash = 5381;
    while (*str) {
        hash = ((hash << 5) + hash) + *str++;
    }
    return hash % bucket_count;
}

// Create hash table from genome mappings for O(1) lookups
static genome_hash_table_t* create_genome_hash_table(genome_mapping_t *mappings, int n_mappings) {
    if (!mappings || n_mappings <= 0) return NULL;
    
    // Use prime number of buckets (roughly n_mappings / 4 for good distribution)
    int bucket_count = n_mappings / 4;
    if (bucket_count < 1024) bucket_count = 1024;
    
    genome_hash_table_t *table = malloc(sizeof(genome_hash_table_t));
    if (!table) return NULL;
    
    table->bucket_count = bucket_count;
    table->buckets = calloc(bucket_count, sizeof(genome_hash_entry_t*));
    if (!table->buckets) {
        free(table);
        return NULL;
    }
    
    // Insert all mappings into hash table
    for (int i = 0; i < n_mappings; i++) {
        unsigned int bucket = hash_genome_name(mappings[i].compressed_id, bucket_count);
        
        genome_hash_entry_t *entry = malloc(sizeof(genome_hash_entry_t));
        if (!entry) continue;  // Skip on allocation failure
        
        entry->compressed_id = mappings[i].compressed_id;
        entry->full_name = mappings[i].full_name;
        entry->next = table->buckets[bucket];  // Chain to existing entries
        table->buckets[bucket] = entry;
    }
    
    return table;
}

// Fast hash table lookup for genome names
static char* hash_table_lookup(genome_hash_table_t *table, const char *compressed_id) {
    if (!table || !compressed_id) return NULL;
    
    unsigned int bucket = hash_genome_name(compressed_id, table->bucket_count);
    genome_hash_entry_t *entry = table->buckets[bucket];
    
    while (entry) {
        if (strcmp(entry->compressed_id, compressed_id) == 0) {
            return entry->full_name;
        }
        entry = entry->next;
    }
    
    return NULL;  // Not found
}

// Free hash table memory
static void free_genome_hash_table(genome_hash_table_t *table) {
    if (!table) return;
    
    for (int i = 0; i < table->bucket_count; i++) {
        genome_hash_entry_t *entry = table->buckets[i];
        while (entry) {
            genome_hash_entry_t *next = entry->next;
            free(entry);
            entry = next;
        }
    }
    
    free(table->buckets);
    free(table);
}

// Global hash table for fast lookups (initialized once)
static genome_hash_table_t *global_genome_hash_table = NULL;

// Initialize hash table from mappings (call once at startup)
void init_genome_hash_table(genome_mapping_t *mappings, int n_mappings) {
    if (global_genome_hash_table) {
        free_genome_hash_table(global_genome_hash_table);
    }
    global_genome_hash_table = create_genome_hash_table(mappings, n_mappings);
}

// Optimized translate function using hash table
char* translate_genome_name(const char *compressed_id, genome_mapping_t *mappings, int n_mappings) {
    if (!compressed_id) return NULL;
    
    // Use hash table if available, otherwise fall back to linear search
    if (global_genome_hash_table) {
        return hash_table_lookup(global_genome_hash_table, compressed_id);
    }
    
    // Fallback to linear search (for compatibility)
    if (!mappings) return NULL;
    
    for (int i = 0; i < n_mappings; i++) {
        if (strcmp(compressed_id, mappings[i].compressed_id) == 0) {
            return mappings[i].full_name;
        }
    }
    
    return NULL;  // Not found
}

// Write proportions CSV file with genome mapping and sorting
int write_proportions_csv_with_mapping(const char *filename, em_results_t *results, char **genome_names, 
                                      int n_genomes, genome_mapping_t *mappings, int n_mappings) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        print_error("Cannot create output file: %s\n", filename);
        return -1;
    }
    
    // Create array of genome indices sorted by proportion (descending)
    int *sorted_indices = malloc(n_genomes * sizeof(int));
    if (!sorted_indices) {
        fclose(file);
        return -1;
    }
    
    for (int i = 0; i < n_genomes; i++) {
        sorted_indices[i] = i;
    }
    
    // Simple bubble sort by proportion (descending)
    for (int i = 0; i < n_genomes - 1; i++) {
        for (int j = 0; j < n_genomes - i - 1; j++) {
            if (results->final_proportions[sorted_indices[j]] < results->final_proportions[sorted_indices[j+1]]) {
                int temp = sorted_indices[j];
                sorted_indices[j] = sorted_indices[j+1];
                sorted_indices[j+1] = temp;
            }
        }
    }
    
    // Write header - include per-genome rates if available
    if (results->final_error_rates || results->final_damage_rates) {
        // Full models - include parameter columns
        fprintf(file, "Genome,Proportion");
        if (results->final_error_rates) {
            fprintf(file, ",Error_Rate");
        }
        if (results->final_damage_rates) {
            fprintf(file, ",Damage_Rate");
        }
        fprintf(file, "\n");
        
        // Write data with per-genome parameters, sorted by proportion
        for (int i = 0; i < n_genomes; i++) {
            int idx = sorted_indices[i];
            const char *compressed_name = (genome_names && genome_names[idx]) ? genome_names[idx] : "Unknown";
            
            // Try to translate the name
            char *full_name = translate_genome_name(compressed_name, mappings, n_mappings);
            const char *output_name = full_name ? full_name : compressed_name;
            
            fprintf(file, "%s,%.16f", output_name, results->final_proportions[idx]);
            if (results->final_error_rates) {
                fprintf(file, ",%.16f", results->final_error_rates[idx]);
            }
            if (results->final_damage_rates) {
                fprintf(file, ",%.16f", results->final_damage_rates[idx]);
            }
            fprintf(file, "\n");
        }
    } else {
        // Simple models - just proportions
        fprintf(file, "Genome,Proportion\n");
        for (int i = 0; i < n_genomes; i++) {
            int idx = sorted_indices[i];
            const char *compressed_name = (genome_names && genome_names[idx]) ? genome_names[idx] : "Unknown";
            
            // Try to translate the name
            char *full_name = translate_genome_name(compressed_name, mappings, n_mappings);
            const char *output_name = full_name ? full_name : compressed_name;
            
            fprintf(file, "%s,%.16f\n", output_name, results->final_proportions[idx]);
        }
    }
    
    free(sorted_indices);
    fclose(file);
    return 0;
}

// Parse parameter constraints file
int parse_constraints_file(const char *filename, parameter_constraint_t **constraints, int *n_constraints) {
    if (!filename || !constraints || !n_constraints) return -1;
    
    FILE *file = fopen(filename, "r");
    if (!file) {
        print_error("Cannot open constraints file: %s\n", filename);
        return -1;
    }
    
    // First pass: count lines
    int line_count = 0;
    char buffer[1024];
    while (fgets(buffer, sizeof(buffer), file)) {
        // Skip empty lines and comments
        char *line = trim_whitespace(buffer);
        if (strlen(line) > 0 && line[0] != '#') {
            line_count++;
        }
    }
    
    if (line_count == 0) {
        printf("No constraints found in file: %s\n", filename);
        fclose(file);
        *constraints = NULL;
        *n_constraints = 0;
        return 0;
    }
    
    // Allocate constraint array
    *constraints = malloc(line_count * sizeof(parameter_constraint_t));
    if (!*constraints) {
        fclose(file);
        return -1;
    }
    
    // Second pass: parse constraints
    rewind(file);
    int constraint_idx = 0;
    
    while (fgets(buffer, sizeof(buffer), file) && constraint_idx < line_count) {
        char *line = trim_whitespace(buffer);
        if (strlen(line) == 0 || line[0] == '#') continue;
        
        // Parse constraint line - examples:
        // "species < 0.01"
        // "0.1 < genus < 0.2" 
        // "genus > 0.05"
        
        parameter_constraint_t *constraint = &(*constraints)[constraint_idx];
        constraint->entity_type = NULL;
        constraint->min_bound = -1.0;  // No lower bound
        constraint->max_bound = -1.0;  // No upper bound
        constraint->type = CONSTRAINT_INEQUALITY;  // Default to inequality

        // Parse constraints with support for spaces in entity names
        char entity[256] = {0};
        double value1, value2;
        int parsed = 0;

        // First, check for format: "value < entity < value" (handles spaces in entity)
        if (sscanf(line, "%lf < %255[^<] < %lf", &value1, entity, &value2) == 3) {
            // Trim the entity name
            char *trimmed = trim_whitespace(entity);
            constraint->entity_type = string_duplicate(trimmed);
            constraint->min_bound = value1;
            constraint->max_bound = value2;
            parsed = 1;
        } else {
            // Look for < or > operators
            char *lt_pos = strchr(line, '<');
            char *gt_pos = strchr(line, '>');

            if (lt_pos && (!gt_pos || lt_pos < gt_pos)) {
                // Format: "entity < value"
                // Extract everything before '<' as entity name
                int entity_len = lt_pos - line;
                if (entity_len > 0 && entity_len < 256) {
                    strncpy(entity, line, entity_len);
                    entity[entity_len] = '\0';
                    char *trimmed = trim_whitespace(entity);

                    // Parse the value after '<'
                    if (sscanf(lt_pos + 1, "%lf", &value1) == 1) {
                        constraint->entity_type = string_duplicate(trimmed);
                        constraint->max_bound = value1;
                        parsed = 1;
                    }
                }
            } else if (gt_pos) {
                // Format: "entity > value"
                // Extract everything before '>' as entity name
                int entity_len = gt_pos - line;
                if (entity_len > 0 && entity_len < 256) {
                    strncpy(entity, line, entity_len);
                    entity[entity_len] = '\0';
                    char *trimmed = trim_whitespace(entity);

                    // Parse the value after '>'
                    if (sscanf(gt_pos + 1, "%lf", &value1) == 1) {
                        constraint->entity_type = string_duplicate(trimmed);
                        constraint->min_bound = value1;
                        parsed = 1;
                    }
                }
            }
        }

        // NEW: Try to parse equality constraint: "entity = value"
        if (!parsed) {
            char *eq_pos = strchr(line, '=');
            if (eq_pos) {
                // Format: "entity = value"
                int entity_len = eq_pos - line;
                if (entity_len > 0 && entity_len < 256) {
                    strncpy(entity, line, entity_len);
                    entity[entity_len] = '\0';
                    char *trimmed = trim_whitespace(entity);

                    // Parse the value after '='
                    if (sscanf(eq_pos + 1, "%lf", &value1) == 1) {
                        // VALIDATE: Fixed value must be in [0, 1]
                        if (value1 < 0.0 || value1 > 1.0) {
                            printf("Warning: Equality constraint value %.6f out of valid range [0,1] for '%s'\n",
                                   value1, trimmed);
                            printf("         Constraint will be clamped to valid range.\n");
                            if (value1 < 0.0) value1 = 0.0;
                            if (value1 > 1.0) value1 = 1.0;
                        }
                        constraint->entity_type = string_duplicate(trimmed);
                        constraint->type = CONSTRAINT_EQUALITY;
                        constraint->fixed_value = value1;
                        constraint->min_bound = -1.0;
                        constraint->max_bound = -1.0;
                        parsed = 1;
                    }
                }
            }
        }

        if (!parsed) {
            printf("Warning: Cannot parse constraint line: %s\n", line);
            continue;
        }
        
        // Convert entity type to lowercase for case insensitive matching
        if (constraint->entity_type) {
            for (char *p = constraint->entity_type; *p; p++) {
                *p = tolower(*p);
            }
        }
        
        constraint_idx++;
    }
    
    fclose(file);
    *n_constraints = constraint_idx;

    // Check for conflicts between equality and inequality constraints
    for (int i = 0; i < *n_constraints; i++) {
        if ((*constraints)[i].type == CONSTRAINT_EQUALITY) {
            // Check if there's an inequality constraint for the same entity
            for (int j = 0; j < *n_constraints; j++) {
                if (i != j && (*constraints)[j].type == CONSTRAINT_INEQUALITY &&
                    strcmp((*constraints)[i].entity_type, (*constraints)[j].entity_type) == 0) {

                    // Check if equality value violates inequality bounds
                    double eq_val = (*constraints)[i].fixed_value;
                    int violates = 0;

                    if ((*constraints)[j].min_bound >= 0 && eq_val < (*constraints)[j].min_bound) {
                        violates = 1;
                    }
                    if ((*constraints)[j].max_bound >= 0 && eq_val > (*constraints)[j].max_bound) {
                        violates = 1;
                    }

                    if (violates) {
                        printf("Warning: Equality constraint '%s = %.6f' conflicts with inequality constraint ",
                               (*constraints)[i].entity_type, eq_val);
                        if ((*constraints)[j].min_bound >= 0) printf("> %.3f ", (*constraints)[j].min_bound);
                        if ((*constraints)[j].max_bound >= 0) printf("< %.3f ", (*constraints)[j].max_bound);
                        printf("\n         Equality constraint takes precedence.\n");
                    }
                }
            }
        }
    }

    if (*n_constraints > 0) {
        printf("Loaded %d parameter constraints from %s\n", *n_constraints, filename);
        for (int i = 0; i < *n_constraints; i++) {
            parameter_constraint_t *c = &(*constraints)[i];
            printf("  Constraint %d: %s", i+1, c->entity_type);
            if (c->type == CONSTRAINT_EQUALITY) {
                printf(" = %.6f (fixed)", c->fixed_value);
            } else {
                if (c->min_bound >= 0) printf(" > %.3f", c->min_bound);
                if (c->max_bound >= 0) printf(" < %.3f", c->max_bound);
            }
            printf("\n");
        }
    }

    return 0;
}

// Parse Dirichlet prior weights file
int parse_weights_file(const char *filename, dirichlet_weight_t **weights, int *n_weights) {
    if (!filename || !weights || !n_weights) return -1;

    FILE *file = fopen(filename, "r");
    if (!file) {
        print_error("Cannot open weights file: %s\n", filename);
        return -1;
    }

    // First pass: count lines
    int line_count = 0;
    char buffer[1024];
    while (fgets(buffer, sizeof(buffer), file)) {
        // Skip empty lines and comments
        char *line = trim_whitespace(buffer);
        if (strlen(line) > 0 && line[0] != '#') {
            line_count++;
        }
    }

    if (line_count == 0) {
        printf("No weights found in file: %s\n", filename);
        fclose(file);
        *weights = NULL;
        *n_weights = 0;
        return 0;
    }

    // Allocate weights array
    *weights = malloc(line_count * sizeof(dirichlet_weight_t));
    if (!*weights) {
        fclose(file);
        return -1;
    }

    // Second pass: parse weights
    rewind(file);
    int weight_idx = 0;

    while (fgets(buffer, sizeof(buffer), file) && weight_idx < line_count) {
        char *line = trim_whitespace(buffer);
        if (strlen(line) == 0 || line[0] == '#') continue;

        // Parse weight line format: "taxonomic_level weight"
        // Example: "species 1" or "genus 10" or "all_other 100"

        char taxonomic_level[256] = {0};
        double weight_value;

        // Parse using sscanf to handle spaces in taxonomic level names
        int parsed = sscanf(line, "%255s %lf", taxonomic_level, &weight_value);

        if (parsed != 2) {
            printf("Warning: Cannot parse weight line: %s\n", line);
            continue;
        }

        // Validate weight value (must be positive)
        if (weight_value <= 0.0) {
            printf("Warning: Invalid weight %.6f for '%s' (must be positive), skipping\n",
                   weight_value, taxonomic_level);
            continue;
        }

        // Store the weight
        dirichlet_weight_t *w = &(*weights)[weight_idx];
        w->taxonomic_level = string_duplicate(taxonomic_level);
        w->weight = weight_value;

        if (!w->taxonomic_level) {
            printf("Warning: Memory allocation failed for taxonomic level, skipping\n");
            continue;
        }

        // Convert taxonomic level to lowercase for case-insensitive matching
        for (char *p = w->taxonomic_level; *p; p++) {
            *p = tolower(*p);
        }

        weight_idx++;
    }

    fclose(file);
    *n_weights = weight_idx;

    if (*n_weights > 0) {
        printf("Loaded %d Dirichlet prior weights from %s\n", *n_weights, filename);
        for (int i = 0; i < *n_weights; i++) {
            dirichlet_weight_t *w = &(*weights)[i];
            printf("  Weight %d: %s = %.3f\n", i+1, w->taxonomic_level, w->weight);
        }
    }

    return 0;
}