#include "em_algorithms.h"
#include <stdarg.h>

// Create EM data structure
em_data_t* em_data_create(int n_reads, int n_genomes) {
    if (n_reads <= 0 || n_genomes <= 0) {
        return NULL;
    }
    
    em_data_t *data = malloc(sizeof(em_data_t));
    if (!data) return NULL;
    
    // Initialize basic fields
    data->n_reads = n_reads;
    data->n_genomes = n_genomes;
    data->error_rate = 0.005;  // Default value
    data->damage_rate = 0.05;  // Default damage rate
    data->iteration = 0;
    data->log_likelihood = -INFINITY;
    data->prev_log_likelihood = -INFINITY;
    data->prev_error_rate = 0.005;
    data->prev_damage_rate = 0.05;
    
    // Allocate standard matrices (using short for count matrices)
    data->n_matrix = matrix_create_short(n_reads, n_genomes);
    data->d_matrix = matrix_create_short(n_reads, n_genomes);
    data->Q_matrix = matrix_create(n_reads, n_genomes);
    data->W_matrix = matrix_create(n_reads, n_genomes);
    
    // Allocate damage model matrices (initialized to NULL, allocated when needed)
    data->nd_matrix = NULL;
    data->md_matrix = NULL;
    data->mb_matrix = NULL;
    
    if (!data->n_matrix || !data->d_matrix || !data->Q_matrix || !data->W_matrix) {
        em_data_destroy(data);
        return NULL;
    }
    
    // Allocate vectors
    data->proportions = vector_create(n_genomes);
    data->prev_proportions = vector_create(n_genomes);
    
    if (!data->proportions || !data->prev_proportions) {
        em_data_destroy(data);
        return NULL;
    }
    
    // Initialize proportions to uniform
    for (int i = 0; i < n_genomes; i++) {
        data->proportions[i] = 1.0 / n_genomes;
        data->prev_proportions[i] = 1.0 / n_genomes;
    }
    
    return data;
}

// Destroy EM data structure
void em_data_destroy(em_data_t *data) {
    if (!data) return;
    
    matrix_destroy_short(data->n_matrix, data->n_reads);
    matrix_destroy_short(data->d_matrix, data->n_reads);
    matrix_destroy(data->Q_matrix, data->n_reads);
    matrix_destroy(data->W_matrix, data->n_reads);
    
    // Clean up damage model matrices if allocated
    if (data->nd_matrix) matrix_destroy_short(data->nd_matrix, data->n_reads);
    if (data->md_matrix) matrix_destroy_short(data->md_matrix, data->n_reads);
    if (data->mb_matrix) matrix_destroy_short(data->mb_matrix, data->n_reads);
    
    vector_destroy(data->proportions);
    vector_destroy(data->prev_proportions);
    
    free(data);
}

// Create configuration structure
em_config_t* em_config_create(void) {
    em_config_t *config = malloc(sizeof(em_config_t));
    if (!config) return NULL;
    
    // Set defaults
    config->input_file = NULL;
    config->output_prefix = NULL;
    config->method = METHOD_A;  // Default method
    config->error_rate = 0.005;
    config->damage_rate = 0.05;  // Default damage rate as requested
    config->joint_damage_rate = 0;  // Default: separate damage rates
    config->tolerance = EM_TOLERANCE;
    config->max_iterations = EM_MAX_ITER;
    config->verbose = 0;
    config->min_proportion = -1.0;  // -1 means use default pruning threshold
    config->use_squarem = 1;  // Default to SQUAREM acceleration
    config->data_format = FORMAT_AUTO;  // Auto-detect format by default
    config->force_sparse = 0;
    config->force_dense = 0;
    config->convert_names = 0;  // Default: no name conversion
    config->mapping_file = NULL;
    config->genome_mappings = NULL;
    config->n_mappings = 0;
    config->genome_key_file = NULL;
    config->use_genome_key = 1;  // Default: auto-detect genome key files
    config->constraints_file = NULL;
    config->constraints = NULL;
    config->n_constraints = 0;
    
    return config;
}

// Allocate damage model matrices for existing em_data_t structure
int em_data_allocate_damage_matrices(em_data_t *data) {
    if (!data) return -1;
    
    // Only allocate if not already allocated (using short for memory efficiency)
    if (!data->nd_matrix) {
        data->nd_matrix = matrix_create_short(data->n_reads, data->n_genomes);
        if (!data->nd_matrix) return -1;
    }
    
    if (!data->md_matrix) {
        data->md_matrix = matrix_create_short(data->n_reads, data->n_genomes);
        if (!data->md_matrix) return -1;
    }
    
    if (!data->mb_matrix) {
        data->mb_matrix = matrix_create_short(data->n_reads, data->n_genomes);
        if (!data->mb_matrix) return -1;
    }
    
    return 0;
}

// Destroy configuration structure
void em_config_destroy(em_config_t *config) {
    if (!config) return;
    
    free(config->input_file);
    free(config->output_prefix);
    free(config->mapping_file);
    free(config->constraints_file);
    
    // Free constraints array
    if (config->constraints) {
        for (int i = 0; i < config->n_constraints; i++) {
            free(config->constraints[i].entity_type);
        }
        free(config->constraints);
    }
    
    // Note: genome_mappings is freed separately in main.c to avoid circular dependency
    free(config);
}

// Create results structure
em_results_t* em_results_create(int n_genomes) {
    if (n_genomes <= 0) return NULL;
    
    em_results_t *results = malloc(sizeof(em_results_t));
    if (!results) return NULL;
    
    results->final_proportions = vector_create(n_genomes);
    if (!results->final_proportions) {
        free(results);
        return NULL;
    }
    
    results->final_error_rate = 0.0;
    results->final_log_likelihood = -INFINITY;
    results->converged = 0;
    results->iterations = 0;
    results->genome_names = NULL;
    results->computation_time = 0.0;
    results->final_damage_rate = 0.0;
    results->final_error_rates = NULL;
    results->final_damage_rates = NULL;
    
    // Initialize taxonomic group fields
    results->final_taxonomic_proportions = NULL;
    results->final_taxonomic_error_rates = NULL;
    results->final_taxonomic_damage_rates = NULL;
    results->taxonomic_groups = NULL;
    results->n_taxonomic_groups = 0;
    
    return results;
}

// Destroy results structure
void em_results_destroy(em_results_t *results) {
    if (!results) return;
    
    vector_destroy(results->final_proportions);
    free(results->final_error_rates);
    free(results->final_damage_rates);
    
    // Free taxonomic group fields
    free(results->final_taxonomic_proportions);
    free(results->final_taxonomic_error_rates);
    free(results->final_taxonomic_damage_rates);
    if (results->taxonomic_groups) {
        for (int i = 0; i < results->n_taxonomic_groups; i++) {
            free(results->taxonomic_groups[i].rank);
            free(results->taxonomic_groups[i].name);
        }
        free(results->taxonomic_groups);
    }
    
    free(results);
}

// Create 2D matrix for counts (using short for memory efficiency)
short** matrix_create_short(int rows, int cols) {
    if (rows <= 0 || cols <= 0) return NULL;
    
    short **matrix = malloc(rows * sizeof(short*));
    if (!matrix) return NULL;
    
    for (int i = 0; i < rows; i++) {
        matrix[i] = calloc(cols, sizeof(short));  // calloc initializes to zero
        if (!matrix[i]) {
            // Clean up on failure
            for (int j = 0; j < i; j++) {
                free(matrix[j]);
            }
            free(matrix);
            return NULL;
        }
    }
    
    return matrix;
}

// Create 2D matrix
double** matrix_create(int rows, int cols) {
    if (rows <= 0 || cols <= 0) return NULL;
    
    double **matrix = malloc(rows * sizeof(double*));
    if (!matrix) return NULL;
    
    for (int i = 0; i < rows; i++) {
        matrix[i] = malloc(cols * sizeof(double));
        if (!matrix[i]) {
            // Clean up on failure
            for (int j = 0; j < i; j++) {
                free(matrix[j]);
            }
            free(matrix);
            return NULL;
        }
        
        // Initialize to zero
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = 0.0;
        }
    }
    
    return matrix;
}

// Destroy 2D matrix (short version)
void matrix_destroy_short(short **matrix, int rows) {
    if (!matrix) return;
    
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

// Destroy 2D matrix
void matrix_destroy(double **matrix, int rows) {
    if (!matrix) return;
    
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

// Copy matrix
void matrix_copy(double **dest, double **src, int rows, int cols) {
    if (!dest || !src) return;
    
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            dest[i][j] = src[i][j];
        }
    }
}

// Create vector
double* vector_create(int size) {
    if (size <= 0) return NULL;
    
    double *vector = malloc(size * sizeof(double));
    if (!vector) return NULL;
    
    // Initialize to zero
    for (int i = 0; i < size; i++) {
        vector[i] = 0.0;
    }
    
    return vector;
}

// Destroy vector
void vector_destroy(double *vector) {
    free(vector);
}

// Copy vector
void vector_copy(double *dest, double *src, int size) {
    if (!dest || !src) return;
    
    for (int i = 0; i < size; i++) {
        dest[i] = src[i];
    }
}

// Safe math functions
double safe_log(double x) {
    if (x <= 0.0) {
        return log(EM_DOUBLE_MIN);
    }
    return log(x);
}

double safe_exp(double x) {
    if (x < -700.0) {  // Prevent underflow
        return EM_DOUBLE_MIN;
    }
    if (x > 700.0) {   // Prevent overflow
        return EM_DOUBLE_MAX;
    }
    return exp(x);
}

// Check if pruning should occur at given iteration (same schedule as sparse)
int dense_should_prune_at_iteration(int iter) {
    // Prune at iterations: 1,2,3,4,5,10,20,30,40,50,60,70,80,90,100,150,200,250,300,350,400,450,500,550,...
    if (iter <= 5) return 1;
    if (iter == 10) return 1;
    if (iter >= 20 && iter <= 100 && iter % 10 == 0) return 1;  // 20,30,40,50,60,70,80,90,100
    if (iter == 150 || iter == 200 || iter == 250) return 1;
    if (iter >= 300 && iter % 50 == 0) return 1;  // 300,350,400,450,500,550,...
    return 0;
}

// Prune genomes and reads based on proportion threshold for dense format
void dense_prune_genomes_and_reads(em_data_t *data, double threshold) {
    if (!data) return;
    
    // Simply set matrix values to -1 for pruned genomes
    // The existing likelihood calculations will skip these automatically
    int genomes_pruned = 0;
    for (int j = 0; j < data->n_genomes; j++) {
        if (data->proportions[j] < threshold) {
            genomes_pruned++;
            
            // Set all matrix values to -1 for this genome column
            for (int i = 0; i < data->n_reads; i++) {
                if (data->d_matrix) {
                    data->d_matrix[i][j] = -1;
                }
                if (data->n_matrix) {
                    data->n_matrix[i][j] = -1;
                }
                if (data->nd_matrix) {
                    data->nd_matrix[i][j] = -1;
                }
                if (data->md_matrix) {
                    data->md_matrix[i][j] = -1;
                }
                if (data->mb_matrix) {
                    data->mb_matrix[i][j] = -1;
                }
            }
            
            // Set proportion to 0 to ensure it stays pruned
            data->proportions[j] = 0.0;
        }
    }
    
    if (genomes_pruned > 0) {
        printf("  Pruned %d genomes (< %.2e)\n", genomes_pruned, threshold);
        fflush(stdout);
    }
}

// Normalize proportions to sum to 1
void normalize_proportions(double *proportions, int n_genomes) {
    if (!proportions || n_genomes <= 0) return;
    
    double sum = 0.0;
    for (int i = 0; i < n_genomes; i++) {
        if (proportions[i] < 0.0) proportions[i] = 0.0;
        sum += proportions[i];
    }
    
    if (sum > 0.0) {
        for (int i = 0; i < n_genomes; i++) {
            proportions[i] /= sum;
        }
    } else {
        // If all proportions are zero, set to uniform
        for (int i = 0; i < n_genomes; i++) {
            proportions[i] = 1.0 / n_genomes;
        }
    }
}

// Initialize proportions uniformly
void initialize_proportions_uniform(em_data_t *data) {
    if (!data || !data->proportions) return;
    
    double uniform_prop = 1.0 / data->n_genomes;
    for (int i = 0; i < data->n_genomes; i++) {
        data->proportions[i] = uniform_prop;
    }
}

// Filter proportions below threshold and renormalize
void filter_proportions(double *proportions, int n_genomes, double min_proportion) {
    if (!proportions || n_genomes <= 0 || min_proportion < 0.0) return;
    
    // Set proportions below threshold to 0
    for (int i = 0; i < n_genomes; i++) {
        if (proportions[i] < min_proportion) {
            proportions[i] = 0.0;
        }
    }
    
    // Renormalize
    normalize_proportions(proportions, n_genomes);
}

// Initialize proportions using data-driven approach (similar to R implementation)
void initialize_proportions_data_driven(em_data_t *data) {
    if (!data || !data->proportions) return;
    
    // First calculate Q matrix with current error rate
    calcQ_method_A(data);  // Use method A for initialization
    
    // Calculate average log-likelihood for each genome
    double *avg_loglik = vector_create(data->n_genomes);
    if (!avg_loglik) {
        initialize_proportions_uniform(data);
        return;
    }
    
    for (int j = 0; j < data->n_genomes; j++) {
        double sum_loglik = 0.0;
        for (int i = 0; i < data->n_reads; i++) {
            sum_loglik += safe_log(data->Q_matrix[i][j] + EM_DOUBLE_MIN);
        }
        avg_loglik[j] = sum_loglik / data->n_reads;
    }
    
    // Find maximum for numerical stability
    double max_loglik = avg_loglik[0];
    for (int j = 1; j < data->n_genomes; j++) {
        if (avg_loglik[j] > max_loglik) {
            max_loglik = avg_loglik[j];
        }
    }
    
    // Convert to weights and normalize
    double sum_weights = 0.0;
    for (int j = 0; j < data->n_genomes; j++) {
        data->proportions[j] = safe_exp(avg_loglik[j] - max_loglik + 1.0);
        sum_weights += data->proportions[j];
    }
    
    // Normalize and add smoothing
    for (int j = 0; j < data->n_genomes; j++) {
        data->proportions[j] = (data->proportions[j] / sum_weights) * 0.9 + 
                              0.1 / data->n_genomes;
    }
    
    vector_destroy(avg_loglik);
}