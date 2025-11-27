#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include "em_types.h"

// Sparse matrix allocation and deallocation
sparse_alignment_t* sparse_alignment_alloc(int n_reads, int n_genomes, int nnz);
void sparse_alignment_free(sparse_alignment_t *sparse);
int sparse_alignment_alloc_damage(sparse_alignment_t *sparse);

// Sparse matrix utilities
int sparse_get_row_nnz(const sparse_alignment_t *sparse, int row);
int sparse_find_col_index(const sparse_alignment_t *sparse, int row, int col);
double sparse_get_n_value(const sparse_alignment_t *sparse, int row, int col);
double sparse_get_d_value(const sparse_alignment_t *sparse, int row, int col);

// Conversion utilities
sparse_alignment_t* dense_to_sparse(short **n_matrix, short **d_matrix, 
                                   int n_reads, int n_genomes);
void sparse_to_dense(const sparse_alignment_t *sparse, 
                     short **n_matrix, short **d_matrix);

// Sparse matrix operations
void sparse_matrix_multiply_vector(const sparse_alignment_t *sparse, 
                                  const double *vector, double *result);
double sparse_matrix_sum(const sparse_alignment_t *sparse);
void sparse_matrix_row_sums(const sparse_alignment_t *sparse, double *row_sums);

// Format detection and parsing
data_format_t detect_file_format(const char *filename);
sparse_alignment_t* parse_sparse_file(const char *filename, char ***genome_names);
sparse_alignment_t* parse_sparse_standard_file(const char *filename, char ***genome_names);
sparse_alignment_t* parse_sparse_damage_file(const char *filename, char ***genome_names);
sparse_alignment_t* parse_dense_file_to_sparse(const char *filename, char ***genome_names);

// Debug and validation utilities
void sparse_matrix_print(const sparse_alignment_t *sparse);
int sparse_matrix_validate(const sparse_alignment_t *sparse);
void sparse_matrix_stats(const sparse_alignment_t *sparse);

#endif // SPARSE_MATRIX_H