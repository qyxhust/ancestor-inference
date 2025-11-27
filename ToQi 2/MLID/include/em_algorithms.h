#ifndef EM_ALGORITHMS_H
#define EM_ALGORITHMS_H

#include "em_types.h"
#include "squarem.h"

// Memory management functions
em_data_t* em_data_create(int n_reads, int n_genomes);
void em_data_destroy(em_data_t *data);
int em_data_allocate_damage_matrices(em_data_t *data);
em_config_t* em_config_create(void);
void em_config_destroy(em_config_t *config);
em_results_t* em_results_create(int n_genomes);
void em_results_destroy(em_results_t *results);

// Matrix operations
double** matrix_create(int rows, int cols);
void matrix_destroy(double **matrix, int rows);
short** matrix_create_short(int rows, int cols);
void matrix_destroy_short(short **matrix, int rows);
void matrix_copy(double **dest, double **src, int rows, int cols);
double* vector_create(int size);
void vector_destroy(double *vector);
void vector_copy(double *dest, double *src, int size);

// Q matrix calculation functions (corresponding to R calcQ methods)
void calcQ_method_A(em_data_t *data);    // No binomial coefficient
void calcQ_method_B(em_data_t *data);    // Error rate divided by 3
void calcQ_method_AE(em_data_t *data);   // Same as A, used in rate estimation
void calcQ_method_BE(em_data_t *data);   // Same as B, used in rate estimation
void calcQ_method_C(em_data_t *data);    // Damage model (fixed rates)
void calcQ_method_CE(em_data_t *data);   // Damage model (fixed damage, estimated background)
void calcQ_method_CED(em_data_t *data);  // Damage model (estimated rates)

// Full model Q matrix calculation functions
void calcQ_method_BEfull(em_data_t *data);   // Per-genome error rates (standard model)
void calcQ_method_CEfull(em_data_t *data);   // Per-genome error rates (damage model)
void calcQ_method_CEDfull(em_data_t *data);  // Per-genome error and damage rates

// EM algorithm implementations
em_results_t* method_A_fixed_rate(em_data_t *data, em_config_t *config);
em_results_t* method_B_fixed_rate(em_data_t *data, em_config_t *config);
em_results_t* method_AE_estimated_rate(em_data_t *data, em_config_t *config);
em_results_t* method_BE_estimated_rate(em_data_t *data, em_config_t *config);
em_results_t* method_C_fixed_rates(em_data_t *data, em_config_t *config);
em_results_t* method_CE_estimated_background(em_data_t *data, em_config_t *config);
em_results_t* method_CED_estimated_rates(em_data_t *data, em_config_t *config);

// Full model implementations (per-genome rates)
em_results_t* method_BEfull_per_genome_rates(em_data_t *data, em_config_t *config);
em_results_t* method_CEfull_per_genome_rates(em_data_t *data, em_config_t *config);
em_results_t* method_CEDfull_per_genome_rates(em_data_t *data, em_config_t *config);

// Core EM steps
void em_e_step(em_data_t *data);        // E-step: Calculate posterior weights
void em_m_step_proportions(em_data_t *data);  // M-step: Update proportions
double em_m_step_error_rate(em_data_t *data); // M-step: Update error rate
double em_m_step_damage_rate(em_data_t *data); // M-step: Update damage rate
double em_m_step_background_rate(em_data_t *data); // M-step: Update background rate
void em_m_step_per_genome_error_rates(em_data_t *data); // M-step: Update per-genome error rates
void em_m_step_per_genome_damage_rates(em_data_t *data); // M-step: Update per-genome damage rates
double calculate_log_likelihood(em_data_t *data);

// Convergence checking
int check_convergence(em_data_t *data, double tolerance);
int check_convergence_with_rate(em_data_t *data, double tolerance);
int check_convergence_proportions_only(em_data_t *data, double tolerance);
int check_convergence_proportions_and_error_rate(em_data_t *data, double tolerance);
int check_convergence_damage_model(em_data_t *data, double tolerance);

// Utility functions
void initialize_proportions_uniform(em_data_t *data);
void initialize_proportions_data_driven(em_data_t *data);
void normalize_proportions(double *proportions, int n_genomes);
void filter_proportions(double *proportions, int n_genomes, double min_proportion);
double safe_log(double x);
double safe_exp(double x);

// Dynamic pruning functions for dense format
void dense_prune_genomes_and_reads(em_data_t *data, double threshold);
int dense_should_prune_at_iteration(int iter);
void setup_pruning_if_needed(em_data_t *data, em_config_t *config);
void apply_pruning_if_needed(em_data_t *data, em_config_t *config, int iteration);

// Optimization interface (placeholder for future SQUAREM implementation)
typedef struct {
    double *parameters;     // Parameters being optimized
    int n_params;          // Number of parameters
    double tolerance;      // Convergence tolerance
    int max_iter;         // Maximum iterations
    em_data_t *data;      // Data for objective function
    em_method_t method;   // Method type
} optimization_context_t;

// SQUAREM acceleration optimization
em_results_t* optimize_with_squarem(em_data_t *data, em_config_t *config);

#endif // EM_ALGORITHMS_H