#include "em_algorithms.h"
#include "io_utils.h"

// Q matrix calculation for Method A (no binomial coefficient)
// Q = ε^(d_ij) * (1-ε)^(n_ij-d_ij)
void calcQ_method_A(em_data_t *data) {
    if (!data || !data->Q_matrix || !data->n_matrix || !data->d_matrix) return;
    
    double log_error = safe_log(data->error_rate);
    double log_1_minus_error = safe_log(1.0 - data->error_rate);
    
    for (int i = 0; i < data->n_reads; i++) {
        for (int j = 0; j < data->n_genomes; j++) {
            // Skip if no alignment (pruned genome)
            if (data->d_matrix[i][j] < 0 || data->n_matrix[i][j] < 0) {
                data->Q_matrix[i][j] = 0.0;
                continue;
            }
            
            double d_ij = (double)data->d_matrix[i][j];
            double n_ij = (double)data->n_matrix[i][j];
            double log_Q = d_ij * log_error + 
                          (n_ij - d_ij) * log_1_minus_error;
            
            data->Q_matrix[i][j] = safe_exp(log_Q);
            
            // Handle edge cases
            if (!isfinite(data->Q_matrix[i][j]) || data->Q_matrix[i][j] <= 0.0) {
                data->Q_matrix[i][j] = EM_DOUBLE_MIN;
            }
        }
    }
}

// Q matrix calculation for Method B (error rate divided by 3)
// Q = (ε/3)^(d_ij) * (1-ε)^(n_ij-d_ij)
void calcQ_method_B(em_data_t *data) {
    if (!data || !data->Q_matrix || !data->n_matrix || !data->d_matrix) return;
    
    double error_div3 = data->error_rate / 3.0;
    if (error_div3 < EM_DOUBLE_MIN) error_div3 = EM_DOUBLE_MIN;
    
    double log_error_div3 = safe_log(error_div3);
    double log_1_minus_error = safe_log(1.0 - data->error_rate);
    
    for (int i = 0; i < data->n_reads; i++) {
        for (int j = 0; j < data->n_genomes; j++) {
            // Skip if no alignment (pruned genome)
            if (data->d_matrix[i][j] < 0 || data->n_matrix[i][j] < 0) {
                data->Q_matrix[i][j] = 0.0;
                continue;
            }
            
            double d_ij = (double)data->d_matrix[i][j];
            double n_ij = (double)data->n_matrix[i][j];
            double log_Q = d_ij * log_error_div3 + 
                          (n_ij - d_ij) * log_1_minus_error;
            
            data->Q_matrix[i][j] = safe_exp(log_Q);
            
            // Handle edge cases
            if (!isfinite(data->Q_matrix[i][j]) || data->Q_matrix[i][j] <= 0.0) {
                data->Q_matrix[i][j] = EM_DOUBLE_MIN;
            }
        }
    }
}

// Q matrix calculation for Method 0AE (same as A, used in rate estimation)
void calcQ_method_AE(em_data_t *data) {
    calcQ_method_A(data);
}

// Q matrix calculation for Method BE (same as B, used in rate estimation)
void calcQ_method_BE(em_data_t *data) {
    calcQ_method_B(data);
}

// Q matrix calculation for Damage Model C (fixed damage and background rates)
// Q = (1-ed)^(nd-md) * ed^md * (1-eb)^(n-md-mb) * (eb/3)^mb
void calcQ_method_C(em_data_t *data) {
    if (!data || !data->Q_matrix || !data->n_matrix || !data->nd_matrix || 
        !data->md_matrix || !data->mb_matrix) return;
    
    
    double log_damage = safe_log(data->damage_rate);
    double log_1_minus_damage = safe_log(1.0 - data->damage_rate);
    double background_div3 = data->error_rate / 3.0;
    if (background_div3 < EM_DOUBLE_MIN) background_div3 = EM_DOUBLE_MIN;
    double log_background_div3 = safe_log(background_div3);
    double log_1_minus_background = safe_log(1.0 - data->error_rate);
    
    for (int i = 0; i < data->n_reads; i++) {
        for (int j = 0; j < data->n_genomes; j++) {
            double nd = (double)data->nd_matrix[i][j];
            double md = (double)data->md_matrix[i][j];
            double mb = (double)data->mb_matrix[i][j];
            double n = (double)data->n_matrix[i][j];
            
            // Skip missing alignments (indicated by -1 values)
            if (nd < 0 || md < 0 || mb < 0) {
                data->Q_matrix[i][j] = 0.0;  // No contribution to likelihood
                continue;
            }
            
            // Calculate log-likelihood: log[(1-ed)^(nd-md) * ed^md * (1-eb)^(n-md-mb) * (eb/3)^mb]
            double log_Q = (nd - md) * log_1_minus_damage +
                          (n - md - mb) * log_1_minus_background +
                          mb * log_background_div3;
            
            // Handle damage error term: when md=0, ed^0 = 1 regardless of ed value
            if (md > 0) {
                log_Q += md * log_damage;
            }
            
            data->Q_matrix[i][j] = safe_exp(log_Q);
            
            // Handle edge cases
            if (!isfinite(data->Q_matrix[i][j]) || data->Q_matrix[i][j] <= 0.0) {
                data->Q_matrix[i][j] = EM_DOUBLE_MIN;
            }
            
        }
    }
}

// Q matrix calculation for Method CE (fixed damage rate, estimated background rate)
void calcQ_method_CE(em_data_t *data) {
    calcQ_method_C(data);
}

// Q matrix calculation for Method CED (estimated damage and background rates)
void calcQ_method_CED(em_data_t *data) {
    calcQ_method_C(data);
}

// E-step: Calculate posterior weights
void em_e_step(em_data_t *data) {
    if (!data || !data->Q_matrix || !data->W_matrix || !data->proportions) return;
    
    // Calculate W = Q * proportions for each read
    for (int i = 0; i < data->n_reads; i++) {
        double row_sum = 0.0;
        
        // First pass: calculate weighted likelihoods and sum
        for (int j = 0; j < data->n_genomes; j++) {
            data->W_matrix[i][j] = data->Q_matrix[i][j] * data->proportions[j];
            row_sum += data->W_matrix[i][j];
        }
        
        // Handle degenerate case
        if (row_sum <= 0.0) {
            row_sum = EM_DOUBLE_MIN;
        }
        
        // Second pass: normalize to get posterior probabilities
        for (int j = 0; j < data->n_genomes; j++) {
            data->W_matrix[i][j] /= row_sum;
            
            // Handle numerical issues
            if (!isfinite(data->W_matrix[i][j])) {
                data->W_matrix[i][j] = 0.0;
            }
        }
    }
}

// M-step: Update proportions
void em_m_step_proportions(em_data_t *data) {
    if (!data || !data->W_matrix || !data->proportions) return;
    
    // Calculate column means of W matrix
    for (int j = 0; j < data->n_genomes; j++) {
        double sum = 0.0;
        for (int i = 0; i < data->n_reads; i++) {
            sum += data->W_matrix[i][j];
        }
        data->proportions[j] = sum / data->n_reads;
    }
    
    // Normalize proportions (safety check)
    normalize_proportions(data->proportions, data->n_genomes);
}

// M-step: Update error rate (for methods Ae and BE)
double em_m_step_error_rate(em_data_t *data) {
    if (!data || !data->W_matrix || !data->n_matrix || !data->d_matrix) {
        return data->error_rate;
    }
    
    double total_weighted_mismatches = 0.0;
    double total_weighted_sites = 0.0;
    
    for (int i = 0; i < data->n_reads; i++) {
        for (int j = 0; j < data->n_genomes; j++) {
            total_weighted_mismatches += data->W_matrix[i][j] * (double)data->d_matrix[i][j];
            total_weighted_sites += data->W_matrix[i][j] * (double)data->n_matrix[i][j];
        }
    }
    
    if (total_weighted_sites > 0.0) {
        double new_rate = total_weighted_mismatches / total_weighted_sites;
        
        // Bound error rate to reasonable range
        if (new_rate < 0.000001) new_rate = 0.000001;
        if (new_rate > 0.99) new_rate = 0.99;
        
        return new_rate;
    }
    
    return data->error_rate;  // Fallback
}

// M-step: Update damage rate (for methods CED)
double em_m_step_damage_rate(em_data_t *data) {
    if (!data || !data->W_matrix || !data->nd_matrix || !data->md_matrix) {
        return data->damage_rate;
    }
    
    double total_weighted_damage_errors = 0.0;
    double total_weighted_damage_sites = 0.0;
    
    for (int i = 0; i < data->n_reads; i++) {
        for (int j = 0; j < data->n_genomes; j++) {
            // Skip missing data (indicated by -1 values)
            if (data->nd_matrix[i][j] < 0 || data->md_matrix[i][j] < 0) {
                continue;
            }
            
            total_weighted_damage_errors += data->W_matrix[i][j] * (double)data->md_matrix[i][j];
            total_weighted_damage_sites += data->W_matrix[i][j] * (double)data->nd_matrix[i][j];
        }
    }
    
    if (total_weighted_damage_sites > 0.0) {
        double new_rate = total_weighted_damage_errors / total_weighted_damage_sites;
        
        // Bound damage rate to reasonable range
        // Use minimum of 1e-6 to avoid numerical issues with log(0)
        if (new_rate < 1e-6) new_rate = 1e-6;
        if (new_rate > 0.99) new_rate = 0.99;
        
        return new_rate;
    }
    
    return data->damage_rate;  // Fallback
}

// M-step: Update background error rate for damage models (for methods CE and CED)
double em_m_step_background_rate(em_data_t *data) {
    if (!data || !data->W_matrix || !data->n_matrix || !data->nd_matrix || !data->mb_matrix) {
        return data->error_rate;
    }
    
    double total_weighted_background_errors = 0.0;
    double total_weighted_background_sites = 0.0;
    
    for (int i = 0; i < data->n_reads; i++) {
        for (int j = 0; j < data->n_genomes; j++) {
            // Skip missing data (indicated by -1 values)
            if (data->n_matrix[i][j] < 0 || data->md_matrix[i][j] < 0 || data->mb_matrix[i][j] < 0) {
                continue;
            }
            
            double background_sites = (double)data->n_matrix[i][j] - (double)data->md_matrix[i][j];
            total_weighted_background_errors += data->W_matrix[i][j] * (double)data->mb_matrix[i][j];
            total_weighted_background_sites += data->W_matrix[i][j] * background_sites;
        }
    }
    
    if (total_weighted_background_sites > 0.0) {
        double new_rate = total_weighted_background_errors / total_weighted_background_sites;
        
        // Bound background error rate to reasonable range
        if (new_rate < 0.000001) new_rate = 0.000001;
        if (new_rate > 0.99) new_rate = 0.99;
        
        return new_rate;
    }
    
    return data->error_rate;  // Fallback
}

// Calculate log-likelihood
double calculate_log_likelihood(em_data_t *data) {
    if (!data || !data->Q_matrix || !data->proportions) return -INFINITY;
    
    double log_likelihood = 0.0;
    
    for (int i = 0; i < data->n_reads; i++) {
        double row_sum = 0.0;
        
        for (int j = 0; j < data->n_genomes; j++) {
            row_sum += data->Q_matrix[i][j] * data->proportions[j];
        }
        
        if (row_sum > 0.0) {
            log_likelihood += safe_log(row_sum);
        } else {
            log_likelihood += safe_log(EM_DOUBLE_MIN);
        }
    }
    
    return log_likelihood;
}

// Setup pruning for any method (call this before the EM loop)
void setup_pruning_if_needed(em_data_t *data, em_config_t *config) {
    if (!data || !config) return;
    
    // Determine the pruning threshold
    if (config->min_proportion < 0) {
        // Use default threshold: min(0.00001, 10/n_reads)
        data->min_proportion = fmin(0.00001, 10.0 / data->n_reads);
        printf("Using dynamic pruning threshold: %.2e\n", data->min_proportion);
    } else if (config->min_proportion == 0) {
        // Pruning explicitly disabled
        data->min_proportion = 0;
    } else {
        // Use user-specified threshold
        data->min_proportion = config->min_proportion;
        printf("Using pruning threshold: %.2e\n", data->min_proportion);
    }
    
    // Set up pruning if threshold is positive
    if (data->min_proportion > 0.0) {
        // Open removed reads file
        if (config->input_file) {
            char removed_file[1024];
            snprintf(removed_file, sizeof(removed_file), "%s.removed_reads.txt", config->input_file);
            data->removed_reads_file = fopen(removed_file, "w");
            if (data->removed_reads_file) {
                printf("Removed reads will be saved to: %s\n", removed_file);
            }
        }
    }
}

// Generic pruning function that works for any method
void apply_pruning_if_needed(em_data_t *data, em_config_t *config, int iteration) {
    if (!data || !config || data->min_proportion <= 0.0) return;
    
    // Check if we should prune at this iteration
    if (!dense_should_prune_at_iteration(iteration)) return;
    
    // Apply pruning for dense format
    dense_prune_genomes_and_reads(data, data->min_proportion);
}

// Check convergence based on proportions only
int check_convergence(em_data_t *data, double tolerance) {
    if (!data || !data->proportions || !data->prev_proportions) return 0;
    
    double max_diff = 0.0;
    for (int i = 0; i < data->n_genomes; i++) {
        double diff = fabs(data->proportions[i] - data->prev_proportions[i]);
        if (diff > max_diff) {
            max_diff = diff;
        }
    }
    
    return (max_diff < tolerance);
}

// Check convergence including error rate
int check_convergence_with_rate(em_data_t *data, double tolerance) {
    if (!check_convergence(data, tolerance)) return 0;
    
    double rate_diff = fabs(data->error_rate - data->prev_error_rate);
    return (rate_diff < tolerance);
}

// Check convergence (proportions only) - alias for existing function
int check_convergence_proportions_only(em_data_t *data, double tolerance) {
    return check_convergence(data, tolerance);
}

// Check convergence including background error rate - alias for existing function
int check_convergence_proportions_and_error_rate(em_data_t *data, double tolerance) {
    return check_convergence_with_rate(data, tolerance);
}

// Check convergence for damage models (proportions, damage rate, and background rate)
int check_convergence_damage_model(em_data_t *data, double tolerance) {
    if (!check_convergence(data, tolerance)) return 0;
    
    double error_rate_diff = fabs(data->error_rate - data->prev_error_rate);
    double damage_rate_diff = fabs(data->damage_rate - data->prev_damage_rate);
    
    return (error_rate_diff < tolerance && damage_rate_diff < tolerance);
}

// Method A: Fixed rate EM algorithm
em_results_t* method_A_fixed_rate(em_data_t *data, em_config_t *config) {
    if (!data || !config) return NULL;
    
    print_verbose(config->verbose, "Running Method A (fixed rate, no binomial coefficient)\n");
    
    // Initialize proportions using data-driven approach
    initialize_proportions_data_driven(data);
    
    // Setup pruning if enabled
    setup_pruning_if_needed(data, config);
    
    em_results_t *results = em_results_create(data->n_genomes);
    if (!results) return NULL;
    
    data->iteration = 0;
    data->log_likelihood = -INFINITY;
    data->prev_log_likelihood = -INFINITY;
    
    // EM iterations
    for (int iter = 0; iter < config->max_iterations; iter++) {
        data->iteration = iter;
        
        // Save previous state
        vector_copy(data->prev_proportions, data->proportions, data->n_genomes);
        data->prev_log_likelihood = data->log_likelihood;
        
        // Calculate Q matrix
        calcQ_method_A(data);
        
        // E-step
        em_e_step(data);
        
        // M-step
        em_m_step_proportions(data);
        
        // Apply pruning if enabled
        apply_pruning_if_needed(data, config, data->iteration + 1);
        
        // Calculate log-likelihood
        data->log_likelihood = calculate_log_likelihood(data);
        
        print_verbose(config->verbose, "Iteration %d: Log-likelihood = %.6f\n", 
                     iter + 1, data->log_likelihood);
        
        // Check convergence
        if (check_convergence(data, config->tolerance)) {
            print_verbose(config->verbose, "Converged after %d iterations\n", iter + 1);
            results->converged = 1;
            break;
        }
    }
    
    // Store results
    vector_copy(results->final_proportions, data->proportions, data->n_genomes);
    results->final_error_rate = data->error_rate;  // Fixed rate
    results->final_log_likelihood = data->log_likelihood;
    results->iterations = data->iteration + 1;
    
    return results;
}

// Method B: Fixed rate EM algorithm (error rate / 3)
em_results_t* method_B_fixed_rate(em_data_t *data, em_config_t *config) {
    if (!data || !config) return NULL;
    
    print_verbose(config->verbose, "Running Method B (fixed rate / 3)\n");
    
    // Initialize proportions using uniform approach
    initialize_proportions_uniform(data);
    
    em_results_t *results = em_results_create(data->n_genomes);
    if (!results) return NULL;
    
    data->iteration = 0;
    data->log_likelihood = -INFINITY;
    data->prev_log_likelihood = -INFINITY;
    
    // EM iterations
    for (int iter = 0; iter < config->max_iterations; iter++) {
        data->iteration = iter;
        
        // Save previous state
        vector_copy(data->prev_proportions, data->proportions, data->n_genomes);
        data->prev_log_likelihood = data->log_likelihood;
        
        // Calculate Q matrix
        calcQ_method_B(data);
        
        // E-step
        em_e_step(data);
        
        // M-step
        em_m_step_proportions(data);
        
        // Apply pruning if enabled
        apply_pruning_if_needed(data, config, data->iteration + 1);
        
        // Calculate log-likelihood
        data->log_likelihood = calculate_log_likelihood(data);
        
        print_verbose(config->verbose, "Iteration %d: Log-likelihood = %.6f\n", 
                     iter + 1, data->log_likelihood);
        
        // Check convergence
        if (check_convergence(data, config->tolerance)) {
            print_verbose(config->verbose, "Converged after %d iterations\n", iter + 1);
            results->converged = 1;
            break;
        }
    }
    
    // Store results
    vector_copy(results->final_proportions, data->proportions, data->n_genomes);
    results->final_error_rate = data->error_rate;  // Fixed rate
    results->final_log_likelihood = data->log_likelihood;
    results->iterations = data->iteration + 1;
    
    return results;
}

// Method 0AE: Estimated rate EM algorithm
em_results_t* method_AE_estimated_rate(em_data_t *data, em_config_t *config) {
    if (!data || !config) return NULL;
    
    print_verbose(config->verbose, "Running Method 0AE (estimated rate, no binomial coefficient)\n");
    
    // Initialize proportions using data-driven approach
    initialize_proportions_data_driven(data);
    
    // Setup pruning if enabled
    setup_pruning_if_needed(data, config);
    
    em_results_t *results = em_results_create(data->n_genomes);
    if (!results) return NULL;
    
    data->iteration = 0;
    data->log_likelihood = -INFINITY;
    data->prev_log_likelihood = -INFINITY;
    data->prev_error_rate = data->error_rate;
    
    // EM iterations
    for (int iter = 0; iter < config->max_iterations; iter++) {
        data->iteration = iter;
        
        // Save previous state
        vector_copy(data->prev_proportions, data->proportions, data->n_genomes);
        data->prev_log_likelihood = data->log_likelihood;
        data->prev_error_rate = data->error_rate;
        
        // Calculate Q matrix
        calcQ_method_AE(data);
        
        // E-step
        em_e_step(data);
        
        // M-step: update proportions and error rate
        em_m_step_proportions(data);
        data->error_rate = em_m_step_error_rate(data);
        
        // Apply pruning if enabled (works for both sparse and dense)
        apply_pruning_if_needed(data, config, iter + 1);
        
        // Calculate log-likelihood
        data->log_likelihood = calculate_log_likelihood(data);
        
        print_verbose(config->verbose, "Iteration %d: Log-likelihood = %.6f, Error rate = %.6f\n", 
                     iter + 1, data->log_likelihood, data->error_rate);
        
        // Check convergence
        if (check_convergence_with_rate(data, config->tolerance)) {
            print_verbose(config->verbose, "Converged after %d iterations\n", iter + 1);
            results->converged = 1;
            break;
        }
    }
    
    // Store results
    vector_copy(results->final_proportions, data->proportions, data->n_genomes);
    results->final_error_rate = data->error_rate;
    results->final_log_likelihood = data->log_likelihood;
    results->iterations = data->iteration + 1;
    
    return results;
}

// Method BE: Estimated rate EM algorithm (error rate / 3)
em_results_t* method_BE_estimated_rate(em_data_t *data, em_config_t *config) {
    if (!data || !config) return NULL;
    
    print_verbose(config->verbose, "Running Method BE (estimated rate / 3)\n");
    
    // Initialize proportions using data-driven approach
    initialize_proportions_data_driven(data);
    
    // Setup pruning if enabled
    setup_pruning_if_needed(data, config);
    
    em_results_t *results = em_results_create(data->n_genomes);
    if (!results) return NULL;
    
    data->iteration = 0;
    data->log_likelihood = -INFINITY;
    data->prev_log_likelihood = -INFINITY;
    data->prev_error_rate = data->error_rate;
    
    // EM iterations
    for (int iter = 0; iter < config->max_iterations; iter++) {
        data->iteration = iter;
        
        // Save previous state
        vector_copy(data->prev_proportions, data->proportions, data->n_genomes);
        data->prev_log_likelihood = data->log_likelihood;
        data->prev_error_rate = data->error_rate;
        
        // Calculate Q matrix
        calcQ_method_BE(data);
        
        // E-step
        em_e_step(data);
        
        // M-step: update proportions and error rate
        em_m_step_proportions(data);
        data->error_rate = em_m_step_error_rate(data);
        
        // Apply pruning if enabled
        apply_pruning_if_needed(data, config, data->iteration + 1);
        
        // Calculate log-likelihood
        data->log_likelihood = calculate_log_likelihood(data);
        
        print_verbose(config->verbose, "Iteration %d: Log-likelihood = %.6f, Error rate = %.6f\n", 
                     iter + 1, data->log_likelihood, data->error_rate);
        
        // Check convergence
        if (check_convergence_with_rate(data, config->tolerance)) {
            print_verbose(config->verbose, "Converged after %d iterations\n", iter + 1);
            results->converged = 1;
            break;
        }
    }
    
    // Store results
    vector_copy(results->final_proportions, data->proportions, data->n_genomes);
    results->final_error_rate = data->error_rate;
    results->final_log_likelihood = data->log_likelihood;
    results->iterations = data->iteration + 1;
    
    return results;
}

// Method C: Damage model with fixed damage and background rates
em_results_t* method_C_fixed_rates(em_data_t *data, em_config_t *config) {
    if (!data || !config) return NULL;
    
    // Allocate results structure
    em_results_t *results = em_results_create(data->n_genomes);
    if (!results) return NULL;
    
    // Initialize convergence tracking
    vector_copy(data->prev_proportions, data->proportions, data->n_genomes);
    data->prev_error_rate = data->error_rate;
    data->prev_damage_rate = data->damage_rate;
    data->prev_log_likelihood = -INFINITY;
    
    // EM iterations
    for (data->iteration = 0; data->iteration < config->max_iterations; data->iteration++) {
        // E-step
        calcQ_method_C(data);
        em_e_step(data);
        
        // M-step (proportions only - rates are fixed)
        em_m_step_proportions(data);
        
        // Calculate log-likelihood
        data->log_likelihood = calculate_log_likelihood(data);
        
        // Check convergence
        if (check_convergence_proportions_only(data, config->tolerance)) {
            results->converged = 1;
            break;
        }
        
        // Update previous values
        vector_copy(data->prev_proportions, data->proportions, data->n_genomes);
        data->prev_log_likelihood = data->log_likelihood;
    }
    
    // Store results
    vector_copy(results->final_proportions, data->proportions, data->n_genomes);
    results->final_error_rate = data->error_rate;
    results->final_damage_rate = data->damage_rate;
    results->final_log_likelihood = data->log_likelihood;
    results->iterations = data->iteration + 1;
    
    return results;
}

// Method CE: Damage model with fixed damage rate, estimated background rate
em_results_t* method_CE_estimated_background(em_data_t *data, em_config_t *config) {
    if (!data || !config) return NULL;
    
    print_verbose(config->verbose, "Running Method CE (damage model, estimated background rate)\n");
    
    // Allocate results structure
    em_results_t *results = em_results_create(data->n_genomes);
    if (!results) return NULL;
    
    // Initialize proportions
    initialize_proportions_data_driven(data);
    
    // Setup pruning if enabled
    setup_pruning_if_needed(data, config);
    
    // Initialize convergence tracking
    vector_copy(data->prev_proportions, data->proportions, data->n_genomes);
    data->prev_error_rate = data->error_rate;
    data->prev_damage_rate = data->damage_rate;
    data->prev_log_likelihood = -INFINITY;
    
    // EM iterations
    for (data->iteration = 0; data->iteration < config->max_iterations; data->iteration++) {
        // E-step
        calcQ_method_CE(data);
        em_e_step(data);
        
        // M-step
        em_m_step_proportions(data);
        data->error_rate = em_m_step_background_rate(data);  // Estimate background rate
        
        // Apply pruning if enabled
        apply_pruning_if_needed(data, config, data->iteration + 1);
        
        // Calculate log-likelihood
        data->log_likelihood = calculate_log_likelihood(data);
        
        // Print iteration information
        if (config->verbose && (data->iteration < 5 || (data->iteration + 1) % 10 == 0)) {
            print_verbose(config->verbose, "Iteration %d: Log-likelihood = %.6f, Background rate = %.6f\n", 
                         data->iteration + 1, data->log_likelihood, data->error_rate);
        }
        
        // Check convergence
        if (check_convergence_proportions_and_error_rate(data, config->tolerance)) {
            results->converged = 1;
            print_verbose(config->verbose, "Converged after %d iterations\n", data->iteration + 1);
            break;
        }
        
        // Update previous values
        vector_copy(data->prev_proportions, data->proportions, data->n_genomes);
        data->prev_error_rate = data->error_rate;
        data->prev_log_likelihood = data->log_likelihood;
    }
    
    // Store results
    vector_copy(results->final_proportions, data->proportions, data->n_genomes);
    results->final_error_rate = data->error_rate;
    results->final_damage_rate = data->damage_rate;
    results->final_log_likelihood = data->log_likelihood;
    results->iterations = data->iteration + 1;
    
    return results;
}

// Method CED: Damage model with estimated damage and background rates
em_results_t* method_CED_estimated_rates(em_data_t *data, em_config_t *config) {
    if (!data || !config) return NULL;
    
    print_verbose(config->verbose, "Running Method CED (damage model, estimated rates)\n");
    
    // Allocate results structure
    em_results_t *results = em_results_create(data->n_genomes);
    if (!results) return NULL;
    
    // Initialize proportions
    initialize_proportions_data_driven(data);
    
    // Setup pruning if enabled
    setup_pruning_if_needed(data, config);
    
    // Initialize convergence tracking
    vector_copy(data->prev_proportions, data->proportions, data->n_genomes);
    data->prev_error_rate = data->error_rate;
    data->prev_damage_rate = data->damage_rate;
    data->prev_log_likelihood = -INFINITY;
    
    // EM iterations
    for (data->iteration = 0; data->iteration < config->max_iterations; data->iteration++) {
        // E-step
        calcQ_method_CED(data);
        em_e_step(data);
        
        // M-step
        em_m_step_proportions(data);
        data->error_rate = em_m_step_background_rate(data);  // Estimate background rate
        
        // Apply pruning if enabled
        apply_pruning_if_needed(data, config, data->iteration + 1);
        data->damage_rate = em_m_step_damage_rate(data);     // Estimate damage rate
        
        // Calculate log-likelihood
        data->log_likelihood = calculate_log_likelihood(data);
        
        // Print iteration information
        if (config->verbose && (data->iteration < 5 || (data->iteration + 1) % 10 == 0)) {
            print_verbose(config->verbose, "Iteration %d: Log-likelihood = %.6f, Background rate = %.6f, Damage rate = %.6f\n", 
                         data->iteration + 1, data->log_likelihood, data->error_rate, data->damage_rate);
        }
        
        // Check convergence (need custom function for damage models)
        if (check_convergence_damage_model(data, config->tolerance)) {
            results->converged = 1;
            print_verbose(config->verbose, "Converged after %d iterations\n", data->iteration + 1);
            break;
        }
        
        // Update previous values
        vector_copy(data->prev_proportions, data->proportions, data->n_genomes);
        data->prev_error_rate = data->error_rate;
        data->prev_damage_rate = data->damage_rate;
        data->prev_log_likelihood = data->log_likelihood;
    }
    
    // Store results
    vector_copy(results->final_proportions, data->proportions, data->n_genomes);
    results->final_error_rate = data->error_rate;
    results->final_damage_rate = data->damage_rate;
    results->final_log_likelihood = data->log_likelihood;
    results->iterations = data->iteration + 1;
    
    return results;
}

// Helper function to count active genomes (dense format)
static int count_active_genomes(em_data_t *data) {
    if (!data || !data->d_matrix) {
        return data ? data->n_genomes : 0;
    }
    
    int count = 0;
    for (int j = 0; j < data->n_genomes; j++) {
        // Check if genome j has any valid alignments (not all -1 values)
        int has_alignment = 0;
        for (int i = 0; i < data->n_reads; i++) {
            if (data->d_matrix[i][j] >= 0) {
                has_alignment = 1;
                break;
            }
        }
        if (has_alignment) {
            count++;
        }
    }
    return count;
}

// SQUAREM acceleration optimization
em_results_t* optimize_with_squarem(em_data_t *data, em_config_t *config) {
    if (config->min_proportion > 0.0) {
        // Two-phase optimization: Standard EM with pruning first, then SQUAREM without pruning
        printf("Two-phase optimization: Standard EM (first 10 iterations with pruning) + SQUAREM (remaining iterations without pruning)\n");
        
        // Phase 1: Standard EM with pruning for first 10 iterations
        printf("Phase 1: Standard EM with pruning (10 iterations)...\n");
        fflush(stdout);
        
        // Temporarily modify config for phase 1
        double orig_min_proportion = config->min_proportion;
        int orig_max_iterations = config->max_iterations;
        int orig_use_squarem = config->use_squarem;
        
        config->max_iterations = 10;
        config->use_squarem = 0;  // Use standard EM for phase 1
        
        // Run phase 1 with standard EM
        em_results_t *phase1_results = NULL;
        switch (config->method) {
            case METHOD_A:       phase1_results = method_A_fixed_rate(data, config); break;
            case METHOD_AE:      phase1_results = method_AE_estimated_rate(data, config); break;
            case METHOD_B:       phase1_results = method_B_fixed_rate(data, config); break;
            case METHOD_BE:      phase1_results = method_BE_estimated_rate(data, config); break;
            case METHOD_C:       phase1_results = method_C_fixed_rates(data, config); break;
            case METHOD_CE:      phase1_results = method_CE_estimated_background(data, config); break;
            case METHOD_CED:     phase1_results = method_CED_estimated_rates(data, config); break;
            case METHOD_BEFULL:  phase1_results = method_BEfull_per_genome_rates(data, config); break;
            case METHOD_CEFULL:  phase1_results = method_CEfull_per_genome_rates(data, config); break;
            case METHOD_CEDFULL: phase1_results = method_CEDfull_per_genome_rates(data, config); break;
            default:             
                config->min_proportion = orig_min_proportion;
                config->max_iterations = orig_max_iterations;
                config->use_squarem = orig_use_squarem;
                return NULL;
        }
        
        if (!phase1_results) {
            config->min_proportion = orig_min_proportion;
            config->max_iterations = orig_max_iterations;
            config->use_squarem = orig_use_squarem;
            return NULL;
        }
        
        printf("Phase 1 completed after %d iterations. Active genomes: %d/%d\n", 
               phase1_results->iterations, count_active_genomes(data), data->n_genomes);
        
        // Check if we should continue to phase 2
        int remaining_iterations = orig_max_iterations - phase1_results->iterations;
        if (remaining_iterations <= 0 || phase1_results->converged) {
            printf("Optimization completed in Phase 1\n");
            config->min_proportion = orig_min_proportion;
            config->max_iterations = orig_max_iterations;
            config->use_squarem = orig_use_squarem;
            return phase1_results;
        }
        
        // Phase 2: SQUAREM without pruning
        printf("Phase 2: SQUAREM acceleration without pruning...\n");
        fflush(stdout);
        
        // Restore original config and disable pruning for phase 2
        config->min_proportion = 0.0;  // Disable pruning
        config->max_iterations = remaining_iterations;
        config->use_squarem = orig_use_squarem;
        
        // Run phase 2 with SQUAREM
        em_results_t *phase2_results = NULL;
        switch (config->method) {
            case METHOD_A:       phase2_results = em_with_squarem_A(data, config); break;
            case METHOD_AE:      phase2_results = em_with_squarem_AE(data, config); break;
            case METHOD_B:       phase2_results = em_with_squarem_B(data, config); break;
            case METHOD_BE:      phase2_results = em_with_squarem_BE(data, config); break;
            case METHOD_C:       phase2_results = em_with_squarem_C(data, config); break;
            case METHOD_CE:      phase2_results = em_with_squarem_CE(data, config); break;
            case METHOD_CED:     phase2_results = em_with_squarem_CED(data, config); break;
            case METHOD_BEFULL:  phase2_results = em_with_squarem_BEfull(data, config); break;
            case METHOD_CEFULL:  phase2_results = em_with_squarem_CEfull(data, config); break;
            case METHOD_CEDFULL: phase2_results = em_with_squarem_CEDfull(data, config); break;
            default:             
                em_results_destroy(phase1_results);
                config->min_proportion = orig_min_proportion;
                config->max_iterations = orig_max_iterations;
                return NULL;
        }
        
        // Restore original config
        config->min_proportion = orig_min_proportion;
        config->max_iterations = orig_max_iterations;
        
        if (!phase2_results) {
            return phase1_results;  // Return phase 1 results if phase 2 failed
        }
        
        // Combine results
        phase2_results->iterations += phase1_results->iterations;
        
        printf("Phase 2 completed after %d iterations\n", phase2_results->iterations - phase1_results->iterations);
        printf("Total optimization: %d iterations (%d standard EM + %d SQUAREM)\n", 
               phase2_results->iterations, phase1_results->iterations, 
               phase2_results->iterations - phase1_results->iterations);
        
        em_results_destroy(phase1_results);
        return phase2_results;
        
    } else {
        // Standard SQUAREM (no pruning case)
        switch (config->method) {
            case METHOD_A:       return em_with_squarem_A(data, config);
            case METHOD_AE:      return em_with_squarem_AE(data, config);
            case METHOD_B:       return em_with_squarem_B(data, config);
            case METHOD_BE:      return em_with_squarem_BE(data, config);
            case METHOD_C:       return em_with_squarem_C(data, config);
            case METHOD_CE:      return em_with_squarem_CE(data, config);
            case METHOD_CED:     return em_with_squarem_CED(data, config);
            case METHOD_BEFULL:  return em_with_squarem_BEfull(data, config);
            case METHOD_CEFULL:  return em_with_squarem_CEfull(data, config);
            case METHOD_CEDFULL: return em_with_squarem_CEDfull(data, config);
            default:             return NULL;
        }
    }
}

// ============================================================================
// Full Models with Per-Genome Error Rates
// ============================================================================

// Q matrix calculation for BEfull (per-genome error rates, standard model)
// Q_ij = (ε_j/3)^d_ij * (1-ε_j)^(n_ij-d_ij)
void calcQ_method_BEfull(em_data_t *data) {
    if (!data || !data->Q_matrix || !data->n_matrix || !data->d_matrix || !data->error_rates) return;
    
    for (int i = 0; i < data->n_reads; i++) {
        for (int j = 0; j < data->n_genomes; j++) {
            double error_rate_j = data->error_rates[j];
            double error_div3 = error_rate_j / 3.0;
            if (error_div3 < EM_DOUBLE_MIN) error_div3 = EM_DOUBLE_MIN;
            
            double log_error_div3 = safe_log(error_div3);
            double log_1_minus_error = safe_log(1.0 - error_rate_j);
            
            double d_ij = (double)data->d_matrix[i][j];
            double n_ij = (double)data->n_matrix[i][j];
            double log_Q = d_ij * log_error_div3 + 
                          (n_ij - d_ij) * log_1_minus_error;
            
            data->Q_matrix[i][j] = safe_exp(log_Q);
            
            // Handle edge cases
            if (!isfinite(data->Q_matrix[i][j]) || data->Q_matrix[i][j] <= 0.0) {
                data->Q_matrix[i][j] = EM_DOUBLE_MIN;
            }
        }
    }
}

// Q matrix calculation for CEfull (per-genome error rates, damage model)
// Q_ij = (1-ed)^(nd_ij-md_ij) * ed^md_ij * (1-ε_j)^(n_ij-md_ij-mb_ij) * (ε_j/3)^mb_ij
void calcQ_method_CEfull(em_data_t *data) {
    if (!data || !data->Q_matrix || !data->n_matrix || !data->d_matrix || 
        !data->nd_matrix || !data->md_matrix || !data->mb_matrix || !data->error_rates) return;
    
    double log_1_minus_damage = safe_log(1.0 - data->damage_rate);
    double log_damage = safe_log(data->damage_rate);
    
    for (int i = 0; i < data->n_reads; i++) {
        for (int j = 0; j < data->n_genomes; j++) {
            double error_rate_j = data->error_rates[j];
            double error_div3 = error_rate_j / 3.0;
            if (error_div3 < EM_DOUBLE_MIN) error_div3 = EM_DOUBLE_MIN;
            
            double log_error_div3 = safe_log(error_div3);
            double log_1_minus_error = safe_log(1.0 - error_rate_j);
            
            double nd = (double)data->nd_matrix[i][j];
            double md = (double)data->md_matrix[i][j];
            double mb = (double)data->mb_matrix[i][j];
            double n = (double)data->n_matrix[i][j];
            
            double log_Q = (nd - md) * log_1_minus_damage +
                          (n - md - mb) * log_1_minus_error +
                          mb * log_error_div3;
            if (md > 0) {
                log_Q += md * log_damage;
            }
            
            data->Q_matrix[i][j] = safe_exp(log_Q);
            
            // Handle edge cases
            if (!isfinite(data->Q_matrix[i][j]) || data->Q_matrix[i][j] <= 0.0) {
                data->Q_matrix[i][j] = EM_DOUBLE_MIN;
            }
        }
    }
}

// Q matrix calculation for CEDfull (per-genome error and damage rates)
// Q_ij = (1-ed_j)^(nd_ij-md_ij) * ed_j^md_ij * (1-ε_j)^(n_ij-md_ij-mb_ij) * (ε_j/3)^mb_ij
void calcQ_method_CEDfull(em_data_t *data) {
    if (!data || !data->Q_matrix || !data->n_matrix || !data->d_matrix || 
        !data->nd_matrix || !data->md_matrix || !data->mb_matrix || 
        !data->error_rates || !data->damage_rates) return;
    
    for (int i = 0; i < data->n_reads; i++) {
        for (int j = 0; j < data->n_genomes; j++) {
            double error_rate_j = data->error_rates[j];
            double damage_rate_j = data->damage_rates[j];
            
            double error_div3 = error_rate_j / 3.0;
            if (error_div3 < EM_DOUBLE_MIN) error_div3 = EM_DOUBLE_MIN;
            
            double log_error_div3 = safe_log(error_div3);
            double log_1_minus_error = safe_log(1.0 - error_rate_j);
            double log_1_minus_damage = safe_log(1.0 - damage_rate_j);
            double log_damage = safe_log(damage_rate_j);
            
            double nd = (double)data->nd_matrix[i][j];
            double md = (double)data->md_matrix[i][j];
            double mb = (double)data->mb_matrix[i][j];
            double n = (double)data->n_matrix[i][j];
            
            // Check for missing data (-1 values)
            if (n < 0 || nd < 0 || md < 0 || mb < 0) {
                data->Q_matrix[i][j] = 0.0;  // No alignment = 0 likelihood
                continue;
            }
            
            double log_Q = (nd - md) * log_1_minus_damage +
                          (n - md - mb) * log_1_minus_error +
                          mb * log_error_div3;
            if (md > 0) {
                log_Q += md * log_damage;
            }
            
            data->Q_matrix[i][j] = safe_exp(log_Q);
            
            // Handle edge cases
            if (!isfinite(data->Q_matrix[i][j]) || data->Q_matrix[i][j] <= 0.0) {
                data->Q_matrix[i][j] = EM_DOUBLE_MIN;
            }
        }
    }
}

// M-step for per-genome error rates (BEfull, CEfull, CEDfull)
void em_m_step_per_genome_error_rates(em_data_t *data) {
    if (!data || !data->error_rates || !data->W_matrix) return;
    
    for (int j = 0; j < data->n_genomes; j++) {
        double numerator = 0.0;
        double denominator = 0.0;
        
        for (int i = 0; i < data->n_reads; i++) {
            double weight = data->W_matrix[i][j];
            if (weight > EM_DOUBLE_MIN) {
                if (data->mb_matrix) {
                    // For damage models, use background errors only
                    double mb = (double)data->mb_matrix[i][j];
                    double n = (double)data->n_matrix[i][j];
                    double md = (double)data->md_matrix[i][j];
                    
                    // Skip missing data (-1 values)
                    if (mb >= 0 && n >= 0 && md >= 0) {
                        numerator += weight * mb;
                        denominator += weight * (n - md);
                    }
                } else {
                    // For standard models, use total errors
                    double d = (double)data->d_matrix[i][j];
                    double n = (double)data->n_matrix[i][j];
                    
                    // Skip missing data (-1 values)
                    if (d >= 0 && n >= 0) {
                        numerator += weight * d;
                        denominator += weight * n;
                    }
                }
            }
        }
        
        if (denominator > EM_DOUBLE_MIN) {
            data->error_rates[j] = numerator / denominator;
            if (data->error_rates[j] < 0.000001) data->error_rates[j] = 0.000001;
            if (data->error_rates[j] > 0.99) data->error_rates[j] = 0.99;
            
            // Debug output for first few iterations
        } else {
            data->error_rates[j] = data->error_rate;
        }
    }
}

// M-step for per-genome damage rates (CEDfull only)
void em_m_step_per_genome_damage_rates(em_data_t *data) {
    if (!data || !data->damage_rates || !data->W_matrix || !data->nd_matrix || !data->md_matrix) return;
    
    for (int j = 0; j < data->n_genomes; j++) {
        double numerator = 0.0;
        double denominator = 0.0;
        
        for (int i = 0; i < data->n_reads; i++) {
            double weight = data->W_matrix[i][j];
            double nd = (double)data->nd_matrix[i][j];
            double md = (double)data->md_matrix[i][j];
            
            // Skip missing data (-1 values)
            if (weight > EM_DOUBLE_MIN && nd >= 0 && md >= 0) {
                numerator += weight * md;
                denominator += weight * nd;
            }
        }
        
        if (denominator > EM_DOUBLE_MIN) {
            data->damage_rates[j] = numerator / denominator;
            if (data->damage_rates[j] < EM_DOUBLE_MIN) data->damage_rates[j] = EM_DOUBLE_MIN;
            if (data->damage_rates[j] > 1.0 - EM_DOUBLE_MIN) data->damage_rates[j] = 1.0 - EM_DOUBLE_MIN;
        } else {
            data->damage_rates[j] = data->damage_rate;
        }
    }
}

// Method BEfull: Per-genome error rates (standard model)  
em_results_t* method_BEfull_per_genome_rates(em_data_t *data, em_config_t *config) {
    if (!data || !config) return NULL;
    
    print_verbose(config->verbose, "Running Method BEfull (per-genome error rates, standard model)\n");
    
    em_results_t *results = em_results_create(data->n_genomes);
    if (!results) return NULL;
    
    // Initialize per-genome error rates
    if (!data->error_rates) {
        data->error_rates = vector_create(data->n_genomes);
        if (!data->error_rates) {
            em_results_destroy(results);
            return NULL;
        }
        for (int j = 0; j < data->n_genomes; j++) {
            data->error_rates[j] = data->error_rate;
        }
    }
    
    // Initialize proportions using data-driven approach
    initialize_proportions_data_driven(data);
    
    // Setup pruning if enabled
    setup_pruning_if_needed(data, config);
    
    data->iteration = 0;
    data->log_likelihood = -INFINITY;
    data->prev_log_likelihood = -INFINITY;
    
    // EM iterations with convergence checking
    for (int iter = 0; iter < config->max_iterations; iter++) {
        data->iteration = iter;
        
        // Save previous state
        vector_copy(data->prev_proportions, data->proportions, data->n_genomes);
        data->prev_log_likelihood = data->log_likelihood;
        
        // Calculate Q matrix using per-genome error rates
        calcQ_method_BEfull(data);
        
        // E-step
        em_e_step(data);
        
        // M-step: update proportions and per-genome error rates
        em_m_step_proportions(data);
        em_m_step_per_genome_error_rates(data);
        
        // Calculate log-likelihood
        data->log_likelihood = calculate_log_likelihood(data);
        
        print_verbose(config->verbose, "Iteration %d: Log-likelihood = %.6f\n", 
                     iter + 1, data->log_likelihood);
        
        // Check convergence
        if (check_convergence(data, config->tolerance)) {
            print_verbose(config->verbose, "Converged after %d iterations\n", iter + 1);
            results->converged = 1;
            break;
        }
    }
    
    // Store results
    vector_copy(results->final_proportions, data->proportions, data->n_genomes);
    results->final_log_likelihood = data->log_likelihood;
    results->iterations = data->iteration + 1;
    
    // Store per-genome error rates
    results->final_error_rates = vector_create(data->n_genomes);
    if (results->final_error_rates) {
        vector_copy(results->final_error_rates, data->error_rates, data->n_genomes);
    }
    
    return results;
}

// Method CEfull: Per-genome error rates (damage model)
em_results_t* method_CEfull_per_genome_rates(em_data_t *data, em_config_t *config) {
    if (!data || !config) return NULL;
    
    print_verbose(config->verbose, "Running Method CEfull (per-genome error rates, damage model)\n");
    
    // Check for damage data
    if (!data->nd_matrix || !data->md_matrix || !data->mb_matrix) {
        print_error("CEfull method requires damage data matrices\n");
        return NULL;
    }
    
    em_results_t *results = em_results_create(data->n_genomes);
    if (!results) return NULL;
    
    // Initialize per-genome error rates
    if (!data->error_rates) {
        data->error_rates = vector_create(data->n_genomes);
        if (!data->error_rates) {
            em_results_destroy(results);
            return NULL;
        }
        for (int j = 0; j < data->n_genomes; j++) {
            data->error_rates[j] = data->error_rate;
        }
    }
    
    // Initialize proportions using data-driven approach
    initialize_proportions_data_driven(data);
    
    // Setup pruning if enabled
    setup_pruning_if_needed(data, config);
    
    data->iteration = 0;
    data->log_likelihood = -INFINITY;
    data->prev_log_likelihood = -INFINITY;
    
    // EM iterations
    for (int iter = 0; iter < config->max_iterations; iter++) {
        data->iteration = iter;
        
        // Save previous state
        vector_copy(data->prev_proportions, data->proportions, data->n_genomes);
        data->prev_log_likelihood = data->log_likelihood;
        
        // Calculate Q matrix using damage model with per-genome error rates
        calcQ_method_CEfull(data);
        
        // E-step
        em_e_step(data);
        
        // M-step: update proportions and per-genome error rates
        em_m_step_proportions(data);
        em_m_step_per_genome_error_rates(data);
        
        // Calculate log-likelihood
        data->log_likelihood = calculate_log_likelihood(data);
        
        print_verbose(config->verbose, "Iteration %d: Log-likelihood = %.6f\n", 
                     iter + 1, data->log_likelihood);
        
        // Check convergence
        if (check_convergence(data, config->tolerance)) {
            print_verbose(config->verbose, "Converged after %d iterations\n", iter + 1);
            results->converged = 1;
            break;
        }
    }
    
    // Store results
    vector_copy(results->final_proportions, data->proportions, data->n_genomes);
    results->final_damage_rate = data->damage_rate;  // Fixed damage rate
    results->final_log_likelihood = data->log_likelihood;
    results->iterations = data->iteration + 1;
    
    // Store per-genome error rates
    results->final_error_rates = vector_create(data->n_genomes);
    if (results->final_error_rates) {
        vector_copy(results->final_error_rates, data->error_rates, data->n_genomes);
    }
    
    return results;
}

// Method CEDfull: Per-genome error and damage rates  
em_results_t* method_CEDfull_per_genome_rates(em_data_t *data, em_config_t *config) {
    if (!data || !config) return NULL;
    
    print_verbose(config->verbose, "Running Method CEDfull (per-genome error and damage rates)\n");
    
    // Check for damage data
    if (!data->nd_matrix || !data->md_matrix || !data->mb_matrix) {
        print_error("CEDfull method requires damage data matrices\n");
        return NULL;
    }
    
    em_results_t *results = em_results_create(data->n_genomes);
    if (!results) return NULL;
    
    // Initialize per-genome error rates
    if (!data->error_rates) {
        data->error_rates = vector_create(data->n_genomes);
        if (!data->error_rates) {
            em_results_destroy(results);
            return NULL;
        }
        for (int j = 0; j < data->n_genomes; j++) {
            data->error_rates[j] = data->error_rate;
        }
    }
    
    // Initialize per-genome damage rates
    if (!data->damage_rates) {
        data->damage_rates = vector_create(data->n_genomes);
        if (!data->damage_rates) {
            em_results_destroy(results);
            return NULL;
        }
        for (int j = 0; j < data->n_genomes; j++) {
            data->damage_rates[j] = data->damage_rate;
        }
    }
    
    // Initialize proportions using data-driven approach
    initialize_proportions_data_driven(data);
    
    // Setup pruning if enabled
    setup_pruning_if_needed(data, config);
    
    data->iteration = 0;
    data->log_likelihood = -INFINITY;
    data->prev_log_likelihood = -INFINITY;
    
    // EM iterations
    for (int iter = 0; iter < config->max_iterations; iter++) {
        data->iteration = iter;
        
        // Save previous state
        vector_copy(data->prev_proportions, data->proportions, data->n_genomes);
        data->prev_log_likelihood = data->log_likelihood;
        
        // Calculate Q matrix using damage model with per-genome rates
        calcQ_method_CEDfull(data);
        
        // E-step
        em_e_step(data);
        
        // M-step: update proportions, per-genome error rates, and per-genome damage rates
        em_m_step_proportions(data);
        em_m_step_per_genome_error_rates(data);
        em_m_step_per_genome_damage_rates(data);
        
        // Calculate log-likelihood
        data->log_likelihood = calculate_log_likelihood(data);
        
        print_verbose(config->verbose, "Iteration %d: Log-likelihood = %.6f\n", 
                     iter + 1, data->log_likelihood);
        
        // Check convergence
        if (check_convergence(data, config->tolerance)) {
            print_verbose(config->verbose, "Converged after %d iterations\n", iter + 1);
            results->converged = 1;
            break;
        }
    }
    
    // Store results
    vector_copy(results->final_proportions, data->proportions, data->n_genomes);
    results->final_log_likelihood = data->log_likelihood;
    results->iterations = data->iteration + 1;
    
    // Store per-genome error rates
    results->final_error_rates = vector_create(data->n_genomes);
    if (results->final_error_rates) {
        vector_copy(results->final_error_rates, data->error_rates, data->n_genomes);
    }
    
    // Store per-genome damage rates
    results->final_damage_rates = vector_create(data->n_genomes);
    if (results->final_damage_rates) {
        vector_copy(results->final_damage_rates, data->damage_rates, data->n_genomes);
    }
    
    return results;
}

