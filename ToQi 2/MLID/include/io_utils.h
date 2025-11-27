#ifndef IO_UTILS_H
#define IO_UTILS_H

#include "em_types.h"

// File I/O functions
int read_mismatch_matrix(const char *filename, em_data_t **data, char ***genome_names);
int write_results(const char *output_prefix, em_results_t *results, char **genome_names, int n_genomes, em_config_t *config);
int write_standardized_output(const char *filename, em_results_t *results, char **genome_names, int n_genomes, em_config_t *config);
int write_proportions_csv(const char *filename, em_results_t *results, char **genome_names, int n_genomes);
int write_proportions_csv_with_mapping(const char *filename, em_results_t *results, char **genome_names, int n_genomes, genome_mapping_t *mappings, int n_mappings);
int write_rates_csv(const char *filename, em_results_t *results, char **genome_names, int n_genomes, em_config_t *config);

// Parsing utilities
int parse_tsv_header(const char *line, char ***genome_names, int *n_genomes);
int parse_tsv_data_line(const char *line, int *total_count, int **mismatches, int n_genomes);
int parse_damage_tsv_header(const char *line, char ***genome_names, int *n_genomes);
int parse_damage_tsv_data_line(const char *line, int *total_count, int **nd_values, int **md_values, int **mb_values, int n_genomes);

// Format detection
typedef enum {
    FORMAT_DENSE_STANDARD,      // genome1  genome2  genome3  ...
    FORMAT_DENSE_DAMAGE,        // nd_genome1  md_genome1  mb_genome1  nd_genome2  ...
    FORMAT_SPARSE_STANDARD,     // variable columns with genome names and values
    FORMAT_SPARSE_DAMAGE,       // variable columns with damage data
    FORMAT_UNSUPPORTED          // unsupported format (e.g., dense with taxonomic groups)
} input_format_t;

input_format_t detect_input_format(const char *header_line);
input_format_t detect_input_format_from_file(const char *filename);
input_format_t detect_input_format_with_data(const char *header_line, const char *first_data_line);

// Sparse format parsing utilities
int parse_sparse_tsv_header(const char *line, char ***genome_names, int *n_genomes);
int parse_sparse_tsv_data_line(const char *line, int *total_count, char ***genome_names, int **mismatches, int *n_alignments, char **all_genome_names, int n_total_genomes);
int parse_sparse_damage_tsv_data_line(const char *line, int *total_count, char ***genome_names, int **nd_values, int **md_values, int **mb_values, int *n_alignments, char **all_genome_names, int n_total_genomes);
char* trim_whitespace(char *str);
int count_columns(const char *line, char delimiter);
int is_numeric_string(const char *str);

// Memory management for strings
char** string_array_create(int size);
void string_array_destroy(char **array, int size);
char* string_duplicate(const char *src);

// Validation functions
int validate_input_data(em_data_t *data);
int validate_mismatch_matrix(em_data_t *data);

// Error handling
void print_error(const char *format, ...);
void print_warning(const char *format, ...);
void print_verbose(int verbose, const char *format, ...);


// Genome mapping functions
int are_genome_names_short(char **genome_names, int n_genomes);
char* detect_genome_key_file(const char *input_file);
genome_mapping_t* load_genome_mapping(const char *filename, int *n_mappings);
void free_genome_mapping(genome_mapping_t *mappings, int n_mappings);
char* translate_genome_name(const char *compressed_id, genome_mapping_t *mappings, int n_mappings);
void init_genome_hash_table(genome_mapping_t *mappings, int n_mappings);

// Parameter constraint functions
int parse_constraints_file(const char *filename, parameter_constraint_t **constraints, int *n_constraints);

// Dirichlet prior weights functions
int parse_weights_file(const char *filename, dirichlet_weight_t **weights, int *n_weights);

#endif // IO_UTILS_H