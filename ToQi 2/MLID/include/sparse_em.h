#ifndef SPARSE_EM_H
#define SPARSE_EM_H

#include "em_types.h"
#include "sparse_matrix.h"

// Sparse EM data allocation and deallocation
sparse_em_data_t* sparse_em_data_alloc(const sparse_alignment_t *alignment_data);
void sparse_em_data_free(sparse_em_data_t *sparse_em_data);

// Sparse Q matrix calculation (likelihood computation)
void sparse_calcQ_method_A(sparse_em_data_t *sparse_em_data);
void sparse_calcQ_method_AE(sparse_em_data_t *sparse_em_data);
void sparse_calcQ_method_B(sparse_em_data_t *sparse_em_data);
void sparse_calcQ_method_BE(sparse_em_data_t *sparse_em_data);
void sparse_calcQ_method_BEfull(sparse_em_data_t *sparse_em_data);
void sparse_calcQ_method_C(sparse_em_data_t *sparse_em_data);
void sparse_calcQ_method_CE(sparse_em_data_t *sparse_em_data);
void sparse_calcQ_method_CED(sparse_em_data_t *sparse_em_data);
void sparse_calcQ_method_CEfull(sparse_em_data_t *sparse_em_data);
void sparse_calcQ_method_CEDfull(sparse_em_data_t *sparse_em_data);

// Sparse EM algorithm steps
void sparse_em_e_step(sparse_em_data_t *sparse_em_data);
void sparse_em_m_step(sparse_em_data_t *sparse_em_data, em_method_t method);
double sparse_calculate_log_likelihood(sparse_em_data_t *sparse_em_data);

// Extended functions removed - unified approach used instead

// Sparse EM algorithm with SQUAREM acceleration
int sparse_em_algorithm(sparse_em_data_t *sparse_em_data, em_method_t method, 
                       double tolerance, int max_iterations, int use_squarem, int verbose,
                       double min_proportion, const char *input_filename);

// Sparse fixed-point functions for SQUAREM
void sparse_fixed_point_A(double *x, double *fx, void *params);
void sparse_fixed_point_AE(double *x, double *fx, void *params);
void sparse_fixed_point_B(double *x, double *fx, void *params);
void sparse_fixed_point_BE(double *x, double *fx, void *params);
void sparse_fixed_point_BEfull(double *x, double *fx, void *params);

// Sparse objective functions for SQUAREM
double sparse_objective_A(double *x, void *params);
double sparse_objective_AE(double *x, void *params);
double sparse_objective_B(double *x, void *params);
double sparse_objective_BE(double *x, void *params);
double sparse_objective_BEfull(double *x, void *params);

// Utility functions
void sparse_normalize_proportions(double *proportions, int n_genomes);
int sparse_check_convergence(sparse_em_data_t *sparse_em_data, double tolerance);
void sparse_copy_parameters(sparse_em_data_t *sparse_em_data, double *proportions, double *error_rate);

#endif // SPARSE_EM_H