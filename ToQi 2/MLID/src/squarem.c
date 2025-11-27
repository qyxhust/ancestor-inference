#define _GNU_SOURCE  // For strdup
#include "squarem.h"
#include "em_algorithms.h"
#include "io_utils.h"
#include <math.h>
#include <string.h>

// Forward declaration of EM-specific data structure for SQUAREM
typedef struct {
    em_data_t *data;
    em_config_t *config;
} em_squarem_data_t;

// Create SQUAREM data structure
squarem_t* squarem_create(int n_params) {
    if (n_params <= 0) return NULL;
    
    squarem_t *sq = malloc(sizeof(squarem_t));
    if (!sq) return NULL;
    
    sq->n_params = n_params;
    sq->step_min = 1.0;
    sq->step_max = 1.0;
    sq->mstep = 4.0;
    sq->objfn_inc = 1.0;
    
    // Allocate parameter vectors
    sq->theta_0 = vector_create(n_params);
    sq->theta_1 = vector_create(n_params);
    sq->theta_2 = vector_create(n_params);
    sq->r_vec = vector_create(n_params);
    sq->v_vec = vector_create(n_params);
    sq->theta_new = vector_create(n_params);
    
    if (!sq->theta_0 || !sq->theta_1 || !sq->theta_2 || 
        !sq->r_vec || !sq->v_vec || !sq->theta_new) {
        squarem_destroy(sq);
        return NULL;
    }
    
    return sq;
}

// Destroy SQUAREM data structure
void squarem_destroy(squarem_t *sq) {
    if (!sq) return;
    
    vector_destroy(sq->theta_0);
    vector_destroy(sq->theta_1);
    vector_destroy(sq->theta_2);
    vector_destroy(sq->r_vec);
    vector_destroy(sq->v_vec);
    vector_destroy(sq->theta_new);
    free(sq);
}

// Create SQUAREM control structure
squarem_control_t* squarem_control_create(void) {
    squarem_control_t *control = malloc(sizeof(squarem_control_t));
    if (!control) return NULL;
    
    // Set default values based on R turboEM
    control->maxiter = 1500;
    control->tol = 1e-7;
    control->trace = 0;
    control->step_min0 = 1.0;
    control->step_max0 = 1.0;
    control->mstep = 4.0;
    control->objfn_inc = 1e-10;  // Reject ANY objective increase (monotone EM)
    control->kr = 1;  // SQUAREM-1
    
    return control;
}

void squarem_control_destroy(squarem_control_t *control) {
    free(control);
}

// Create SQUAREM result structure
squarem_result_t* squarem_result_create(int n_params) {
    squarem_result_t *result = malloc(sizeof(squarem_result_t));
    if (!result) return NULL;
    
    result->par = vector_create(n_params);
    if (!result->par) {
        free(result);
        return NULL;
    }
    
    result->value = INFINITY;
    result->iter = 0;
    result->convergence = -1;
    result->message = NULL;
    
    return result;
}

void squarem_result_destroy(squarem_result_t *result) {
    if (!result) return;
    
    vector_destroy(result->par);
    free(result->message);
    free(result);
}

// Vector utility functions
double vector_norm_squared(double *vec, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += vec[i] * vec[i];
    }
    return sum;
}

double vector_dot_product(double *vec1, double *vec2, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += vec1[i] * vec2[i];
    }
    return sum;
}

void vector_copy_squarem(double *dest, double *src, int n) {
    for (int i = 0; i < n; i++) {
        dest[i] = src[i];
    }
}

void vector_add(double *result, double *vec1, double *vec2, int n) {
    for (int i = 0; i < n; i++) {
        result[i] = vec1[i] + vec2[i];
    }
}

void vector_subtract(double *result, double *vec1, double *vec2, int n) {
    for (int i = 0; i < n; i++) {
        result[i] = vec1[i] - vec2[i];
    }
}

void vector_scale_add(double *result, double *vec, double scale, int n) {
    for (int i = 0; i < n; i++) {
        result[i] += scale * vec[i];
    }
}

// Main SQUAREM algorithm implementation
squarem_result_t* squarem(double *par, fixptfn_t fixptfn, objfn_t objfn, 
                         void *data, squarem_control_t *control) {
    if (!par || !fixptfn || !data || !control) return NULL;
    
    // Get EM data and config correctly
    em_squarem_data_t *em_sq_data = (em_squarem_data_t*)data;
    em_data_t *em_data = em_sq_data->data;
    em_config_t *config = em_sq_data->config;
    
    int n_params = 0;
    
    // Determine number of parameters based on method
    if (config->method == METHOD_AE || config->method == METHOD_BE) {
        n_params = em_data->n_genomes;  // proportions + error_rate
    } else if (config->method == METHOD_CE) {
        n_params = em_data->n_genomes;  // proportions + error_rate (damage rate is fixed)
    } else if (config->method == METHOD_CED) {
        n_params = em_data->n_genomes + 1;  // proportions + error_rate + damage_rate
    } else if (config->method == METHOD_BEFULL) {
        n_params = 2 * em_data->n_genomes - 1;  // proportions + per-genome error rates
    } else if (config->method == METHOD_CEFULL) {
        n_params = 2 * em_data->n_genomes;  // proportions + per-genome error rates + damage rate
    } else if (config->method == METHOD_CEDFULL) {
        n_params = 3 * em_data->n_genomes - 1;  // proportions + per-genome error rates + per-genome damage rates
    } else {
        n_params = em_data->n_genomes - 1;  // proportions only (last one is 1-sum)
    }
    
    squarem_result_t *result = squarem_result_create(n_params);
    if (!result) return NULL;
    
    squarem_t *sq = squarem_create(n_params);
    if (!sq) {
        squarem_result_destroy(result);
        return NULL;
    }
    
    // Initialize parameters
    vector_copy_squarem(sq->theta_0, par, n_params);
    
    // Set step size parameters
    sq->step_min = control->step_min0;
    sq->step_max = control->step_max0;
    sq->mstep = control->mstep;
    sq->objfn_inc = control->objfn_inc;
    
    double obj_old = objfn ? objfn(sq->theta_0, data) : 0.0;
    
    printf("Starting SQUAREM iterations (max %d)...\n", control->maxiter);
    fflush(stdout);
    
    for (int iter = 0; iter < control->maxiter; iter++) {
        
        // Show progress every 10 iterations or for first few iterations
        // Always show progress for first 5 iterations and every 10th iteration
        if (iter < 5 || (iter + 1) % 10 == 0) {
            printf("Iteration %d: log-likelihood = %.6f", iter + 1, -obj_old);
            if (iter > 0) {
                // Note: delta calculation would require additional objective evaluation
                printf(" (converging)");
            }
            printf("\n");
            fflush(stdout);
        } else if (control->trace) {
            print_verbose(1, "  [Verbose] Iteration %d\n", iter + 1);
        }
        
        // Standard EM step 1: θ_1 = F(θ_0)
        fixptfn(sq->theta_0, sq->theta_1, data);
        
        // Standard EM step 2: θ_2 = F(θ_1)  
        fixptfn(sq->theta_1, sq->theta_2, data);
        
        // Calculate residual vectors (like R implementation)
        vector_subtract(sq->r_vec, sq->theta_1, sq->theta_0, n_params);  // q1 = θ_1 - θ_0
        vector_subtract(sq->v_vec, sq->theta_2, sq->theta_1, n_params);  // q2 = θ_2 - θ_1
        
        // Multiple convergence criteria checks (like R turboEM)
        double q1_norm = sqrt(vector_norm_squared(sq->r_vec, n_params));
        double q2_norm = sqrt(vector_norm_squared(sq->v_vec, n_params));
        
        // Primary convergence check on first-order difference
        if (q1_norm < control->tol) {
            vector_copy_squarem(result->par, sq->theta_1, n_params);
            result->iter = iter + 1;
            result->convergence = 0;
            result->message = strdup("Converged on q1");
            if (objfn) result->value = objfn(result->par, data);
            break;
        }
        
        // Secondary convergence check on second-order difference
        if (q2_norm < control->tol) {
            vector_copy_squarem(result->par, sq->theta_2, n_params);
            result->iter = iter + 1;
            result->convergence = 0;
            result->message = strdup("Converged on q2");
            if (objfn) result->value = objfn(result->par, data);
            break;
        }
        
        // Calculate q2 - q1 
        double *q2_minus_q1 = vector_create(n_params);
        vector_subtract(q2_minus_q1, sq->v_vec, sq->r_vec, n_params);            // q2 - q1
        
        double sv2 = vector_norm_squared(q2_minus_q1, n_params);                  // sv2 = ||q2 - q1||²
        double srv = vector_dot_product(sq->r_vec, q2_minus_q1, n_params);       // srv = ⟨q1, q2-q1⟩
        
        // Calculate alpha using R SQUAREM method 1: α = -⟨q1, q2-q1⟩ / ||q2-q1||²
        double alpha;
        if (sv2 > 0.0) {
            alpha = -srv / sv2;
        } else {
            alpha = 1.0;  // Fallback
        }
        
        // Ensure alpha is positive for stability
        alpha = fabs(alpha);
        
        // Bound the step size
        alpha = fmax(sq->step_min, fmin(sq->step_max, alpha));
        
        // SQUAREM extrapolation (R formula): θ_new = θ_0 + 2*alpha*q1 + alpha²*(q2-q1)
        vector_copy_squarem(sq->theta_new, sq->theta_0, n_params);              // Start with θ_0
        vector_scale_add(sq->theta_new, sq->r_vec, 2.0 * alpha, n_params);     // + 2*alpha*q1
        vector_scale_add(sq->theta_new, q2_minus_q1, alpha * alpha, n_params); // + alpha²*(q2-q1)
        
        // Clean up
        vector_destroy(q2_minus_q1);
        
        // Safeguarding: Check for invalid extrapolated points (NaN/Inf)
        int invalid_point = 0;
        for (int i = 0; i < n_params; i++) {
            if (isnan(sq->theta_new[i]) || isinf(sq->theta_new[i])) {
                invalid_point = 1;
                break;
            }
        }
        
        if (invalid_point) {
            if (control->trace) {
                print_verbose(1, "  Invalid extrapolation, using standard EM step\n");
            }
            vector_copy_squarem(sq->theta_new, sq->theta_2, n_params);
            // Contract step bounds on failure
            sq->step_max = fmax(control->step_max0, sq->step_max / sq->mstep);
        }
        
        // Check objective function and implement safeguarding (like R turboEM)
        int extrapolation_succeeded = 1;
        if (objfn) {
            double obj_new = objfn(sq->theta_new, data);

            // If objective function increased too much, use a standard EM step instead
            if (obj_new > obj_old + control->objfn_inc) {
                if (control->trace) {
                    print_verbose(1, "  Objective increased, using standard EM step\n");
                }
                vector_copy_squarem(sq->theta_0, sq->theta_2, n_params);
                extrapolation_succeeded = 0;

                // Contract step bounds on failure (R turboEM strategy)
                sq->step_max = fmax(control->step_max0, sq->step_max / sq->mstep);
            } else {
                vector_copy_squarem(sq->theta_0, sq->theta_new, n_params);
                obj_old = obj_new;
            }
        } else {
            vector_copy_squarem(sq->theta_0, sq->theta_new, n_params);
        }
        
        // Dynamic step size adaptation (key R turboEM feature)
        // Only expand step bounds if extrapolation was successful
        if (extrapolation_succeeded) {
            sq->step_max = sq->step_max * sq->mstep;
            // Also adapt minimum step if needed
            if (sq->step_min < 0) {
                sq->step_min = sq->step_min * sq->mstep;
            }
        }
    }
    
    // If we didn't converge in the loop
    if (result->convergence != 0) {
        vector_copy_squarem(result->par, sq->theta_0, n_params);
        result->iter = control->maxiter;
        result->convergence = 1;
        result->message = strdup("Maximum iterations reached");
        if (objfn) result->value = objfn(result->par, data);
    }
    
    squarem_destroy(sq);
    return result;
}

// EM-specific wrapper functions

// Fixed-point function for Method A
void em_fixpoint_function_A(double *par, double *fpar, void *data) {
    em_squarem_data_t *em_sq_data = (em_squarem_data_t*)data;
    em_data_t *em_data = em_sq_data->data;
    
    // Set proportions from parameters with proper normalization
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        em_data->proportions[i] = fmax(0.0, fmin(1.0, par[i]));
    }
    
    // Calculate the last proportion ensuring all are positive
    double sum = 0.0;
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        sum += em_data->proportions[i];
    }
    
    if (sum >= 1.0) {
        // Normalize to ensure sum < 1
        double scale = 0.99 / sum;
        for (int i = 0; i < em_data->n_genomes - 1; i++) {
            em_data->proportions[i] *= scale;
        }
        em_data->proportions[em_data->n_genomes - 1] = 0.01;
    } else {
        em_data->proportions[em_data->n_genomes - 1] = 1.0 - sum;
    }
    
    // One EM step
    calcQ_method_A(em_data);
    em_e_step(em_data);
    em_m_step_proportions(em_data);
    
    // Return new parameters
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        fpar[i] = em_data->proportions[i];
    }
}

// Fixed-point function for Method 0AE
void em_fixpoint_function_AE(double *par, double *fpar, void *data) {
    em_squarem_data_t *em_sq_data = (em_squarem_data_t*)data;
    em_data_t *em_data = em_sq_data->data;
    
    // Set proportions and error rate from parameters
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        em_data->proportions[i] = par[i];
    }
    em_data->proportions[em_data->n_genomes - 1] = 1.0;
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        em_data->proportions[em_data->n_genomes - 1] -= par[i];
    }
    em_data->error_rate = par[em_data->n_genomes - 1];
    
    // One EM step
    calcQ_method_AE(em_data);
    em_e_step(em_data);
    em_m_step_proportions(em_data);
    em_data->error_rate = em_m_step_error_rate(em_data);
    
    // Return new parameters
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        fpar[i] = em_data->proportions[i];
    }
    fpar[em_data->n_genomes - 1] = em_data->error_rate;
}

// Fixed-point function for Method B
void em_fixpoint_function_B(double *par, double *fpar, void *data) {
    em_squarem_data_t *em_sq_data = (em_squarem_data_t*)data;
    em_data_t *em_data = em_sq_data->data;
    
    // Set proportions from parameters with proper normalization
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        em_data->proportions[i] = fmax(0.0, fmin(1.0, par[i]));
    }
    
    // Calculate the last proportion ensuring all are positive
    double sum = 0.0;
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        sum += em_data->proportions[i];
    }
    
    if (sum >= 1.0) {
        // Normalize to ensure sum < 1
        double scale = 0.99 / sum;
        for (int i = 0; i < em_data->n_genomes - 1; i++) {
            em_data->proportions[i] *= scale;
        }
        em_data->proportions[em_data->n_genomes - 1] = 0.01;
    } else {
        em_data->proportions[em_data->n_genomes - 1] = 1.0 - sum;
    }
    
    // One EM step
    calcQ_method_B(em_data);
    em_e_step(em_data);
    em_m_step_proportions(em_data);
    
    // Return new parameters
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        fpar[i] = em_data->proportions[i];
    }
}

// Fixed-point function for Method BE
void em_fixpoint_function_BE(double *par, double *fpar, void *data) {
    em_squarem_data_t *em_sq_data = (em_squarem_data_t*)data;
    em_data_t *em_data = em_sq_data->data;
    
    // Set proportions and error rate from parameters
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        em_data->proportions[i] = par[i];
    }
    em_data->proportions[em_data->n_genomes - 1] = 1.0;
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        em_data->proportions[em_data->n_genomes - 1] -= par[i];
    }
    em_data->error_rate = par[em_data->n_genomes - 1];
    
    // One EM step
    calcQ_method_BE(em_data);
    em_e_step(em_data);
    em_m_step_proportions(em_data);
    em_data->error_rate = em_m_step_error_rate(em_data);
    
    // Return new parameters
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        fpar[i] = em_data->proportions[i];
    }  
    fpar[em_data->n_genomes - 1] = em_data->error_rate;
}

// Objective functions (negative log-likelihood)
double em_objective_function_A(double *par, void *data) {
    em_squarem_data_t *em_sq_data = (em_squarem_data_t*)data;
    em_data_t *em_data = em_sq_data->data;
    
    // Set proportions from parameters with proper normalization
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        em_data->proportions[i] = fmax(0.0, fmin(1.0, par[i]));
    }
    
    // Calculate the last proportion ensuring all are positive
    double sum = 0.0;
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        sum += em_data->proportions[i];
    }
    
    if (sum >= 1.0) {
        // Normalize to ensure sum < 1
        double scale = 0.99 / sum;
        for (int i = 0; i < em_data->n_genomes - 1; i++) {
            em_data->proportions[i] *= scale;
        }
        em_data->proportions[em_data->n_genomes - 1] = 0.01;
    } else {
        em_data->proportions[em_data->n_genomes - 1] = 1.0 - sum;
    }
    
    // Check parameter bounds (like R implementation)
    for (int i = 0; i < em_data->n_genomes; i++) {
        if (em_data->proportions[i] <= 0.0 || em_data->proportions[i] >= 1.0) {
            return 1e10;  // Penalty for invalid proportions
        }
    }
    
    calcQ_method_A(em_data);
    return -calculate_log_likelihood(em_data);
}

double em_objective_function_AE(double *par, void *data) {
    em_squarem_data_t *em_sq_data = (em_squarem_data_t*)data;
    em_data_t *em_data = em_sq_data->data;
    
    // Set proportions and error rate from parameters
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        em_data->proportions[i] = par[i];
    }
    em_data->proportions[em_data->n_genomes - 1] = 1.0;
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        em_data->proportions[em_data->n_genomes - 1] -= par[i];
    }
    em_data->error_rate = par[em_data->n_genomes - 1];
    
    // Check parameter bounds (like R implementation)
    for (int i = 0; i < em_data->n_genomes; i++) {
        if (em_data->proportions[i] <= 0.0 || em_data->proportions[i] >= 1.0) {
            return 1e10;  // Penalty for invalid proportions
        }
    }
    if (em_data->error_rate <= 0.0 || em_data->error_rate >= 1.0) {
        return 1e10;  // Penalty for invalid error rate
    }
    
    calcQ_method_AE(em_data);
    return -calculate_log_likelihood(em_data);
}

double em_objective_function_B(double *par, void *data) {
    em_squarem_data_t *em_sq_data = (em_squarem_data_t*)data;
    em_data_t *em_data = em_sq_data->data;
    
    // Set proportions from parameters with proper normalization
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        em_data->proportions[i] = fmax(0.0, fmin(1.0, par[i]));
    }
    
    // Calculate the last proportion ensuring all are positive
    double sum = 0.0;
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        sum += em_data->proportions[i];
    }
    
    if (sum >= 1.0) {
        // Normalize to ensure sum < 1
        double scale = 0.99 / sum;
        for (int i = 0; i < em_data->n_genomes - 1; i++) {
            em_data->proportions[i] *= scale;
        }
        em_data->proportions[em_data->n_genomes - 1] = 0.01;
    } else {
        em_data->proportions[em_data->n_genomes - 1] = 1.0 - sum;
    }
    
    // Check parameter bounds (like R implementation)
    for (int i = 0; i < em_data->n_genomes; i++) {
        if (em_data->proportions[i] <= 0.0 || em_data->proportions[i] >= 1.0) {
            return 1e10;  // Penalty for invalid proportions
        }
    }
    
    calcQ_method_B(em_data);
    return -calculate_log_likelihood(em_data);
}

double em_objective_function_BE(double *par, void *data) {
    em_squarem_data_t *em_sq_data = (em_squarem_data_t*)data;
    em_data_t *em_data = em_sq_data->data;
    
    // Set proportions and error rate from parameters
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        em_data->proportions[i] = par[i];
    }
    em_data->proportions[em_data->n_genomes - 1] = 1.0;
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        em_data->proportions[em_data->n_genomes - 1] -= par[i];
    }
    em_data->error_rate = par[em_data->n_genomes - 1];
    
    // Check parameter bounds (like R implementation)
    for (int i = 0; i < em_data->n_genomes; i++) {
        if (em_data->proportions[i] <= 0.0 || em_data->proportions[i] >= 1.0) {
            return 1e10;  // Penalty for invalid proportions
        }
    }
    if (em_data->error_rate <= 0.0 || em_data->error_rate >= 1.0) {
        return 1e10;  // Penalty for invalid error rate
    }
    
    calcQ_method_BE(em_data);
    return -calculate_log_likelihood(em_data);
}

// ============================================================================
// Damage Model SQUAREM Functions (C, CE, CED)
// ============================================================================

// C fixpoint function: proportions only (fixed damage and error rates)
// Parameters: [prop_0, ..., prop_{K-2}]
void em_fixpoint_function_C(double *par, double *fpar, void *data) {
    em_squarem_data_t *em_sq_data = (em_squarem_data_t*)data;
    em_data_t *em_data = em_sq_data->data;
    
    // Set proportions with robust constraint handling (like method B)
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        em_data->proportions[i] = fmax(0.0, fmin(1.0, par[i]));
    }
    
    // Calculate the last proportion ensuring all are positive
    double sum = 0.0;
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        sum += em_data->proportions[i];
    }
    
    if (sum >= 1.0) {
        // Normalize to ensure sum < 1
        double scale = 0.99 / sum;
        for (int i = 0; i < em_data->n_genomes - 1; i++) {
            em_data->proportions[i] *= scale;
        }
        em_data->proportions[em_data->n_genomes - 1] = 0.01;
    } else {
        em_data->proportions[em_data->n_genomes - 1] = 1.0 - sum;
    }
    
    // Perform one EM iteration
    calcQ_method_C(em_data);
    em_e_step(em_data);
    em_m_step_proportions(em_data);
    
    // Return new parameters
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        fpar[i] = em_data->proportions[i];
    }
}

double em_objective_function_C(double *par, void *data) {
    em_squarem_data_t *em_sq_data = (em_squarem_data_t*)data;
    em_data_t *em_data = em_sq_data->data;
    
    // Set proportions from parameters (like method B)
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        em_data->proportions[i] = par[i];
    }
    em_data->proportions[em_data->n_genomes - 1] = 1.0;
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        em_data->proportions[em_data->n_genomes - 1] -= par[i];
    }
    
    // Check parameter bounds
    for (int i = 0; i < em_data->n_genomes; i++) {
        if (em_data->proportions[i] <= 0.0 || em_data->proportions[i] >= 1.0) {
            return 1e10;  // Penalty for invalid proportions
        }
    }
    
    // Calculate Q matrix and log-likelihood
    calcQ_method_C(em_data);
    return -calculate_log_likelihood(em_data);
}

// CE fixpoint function: proportions + background error rate
// Parameters: [prop_0, ..., prop_{K-2}, error_rate]
void em_fixpoint_function_CE(double *par, double *fpar, void *data) {
    em_squarem_data_t *em_sq_data = (em_squarem_data_t*)data;
    em_data_t *em_data = em_sq_data->data;
    
    // Set proportions and error rate from parameters with robust bounds (like CEfull)
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        em_data->proportions[i] = fmax(0.0, fmin(1.0, par[i]));
    }
    
    // Calculate last proportion with overflow protection
    double sum = 0.0;
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        sum += em_data->proportions[i];
    }
    
    if (sum >= 1.0) {
        // Scale down proportions if sum >= 1.0
        double scale = 0.99 / sum;
        for (int i = 0; i < em_data->n_genomes - 1; i++) {
            em_data->proportions[i] *= scale;
        }
        em_data->proportions[em_data->n_genomes - 1] = 0.01;
    } else {
        em_data->proportions[em_data->n_genomes - 1] = 1.0 - sum;
    }
    
    // Bound error rate
    em_data->error_rate = fmax(0.000001, fmin(0.99, par[em_data->n_genomes - 1]));
    
    // Perform one EM iteration
    calcQ_method_CE(em_data);
    em_e_step(em_data);
    em_m_step_proportions(em_data);
    em_data->error_rate = em_m_step_background_rate(em_data);
    
    // Apply bounds to M-step result (critical for SQUAREM stability)
    em_data->error_rate = fmax(0.000001, fmin(0.99, em_data->error_rate));
    
    // Return new parameters
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        fpar[i] = em_data->proportions[i];
    }
    fpar[em_data->n_genomes - 1] = em_data->error_rate;
}

double em_objective_function_CE(double *par, void *data) {
    em_squarem_data_t *em_sq_data = (em_squarem_data_t*)data;
    em_data_t *em_data = em_sq_data->data;
    
    // Set proportions and error rate (like method BE)
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        em_data->proportions[i] = par[i];
    }
    em_data->proportions[em_data->n_genomes - 1] = 1.0;
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        em_data->proportions[em_data->n_genomes - 1] -= par[i];
    }
    em_data->error_rate = par[em_data->n_genomes - 1];
    
    // Check parameter bounds
    for (int i = 0; i < em_data->n_genomes; i++) {
        if (em_data->proportions[i] <= 0.0 || em_data->proportions[i] >= 1.0) {
            return 1e10;  // Penalty for invalid proportions
        }
    }
    if (em_data->error_rate <= 0.0 || em_data->error_rate >= 1.0) {
        return 1e10;  // Penalty for invalid error rate
    }
    
    // Calculate Q matrix and log-likelihood
    calcQ_method_CE(em_data);
    return -calculate_log_likelihood(em_data);
}

// CED fixpoint function: proportions + error rate + damage rate
// Parameters: [prop_0, ..., prop_{K-2}, error_rate, damage_rate]
void em_fixpoint_function_CED(double *par, double *fpar, void *data) {
    em_squarem_data_t *em_sq_data = (em_squarem_data_t*)data;
    em_data_t *em_data = em_sq_data->data;
    
    // Set proportions and rates from parameters with robust bounds (like CEDfull)
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        em_data->proportions[i] = fmax(0.0, fmin(1.0, par[i]));
    }
    
    // Calculate last proportion with overflow protection
    double sum = 0.0;
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        sum += em_data->proportions[i];
    }
    
    if (sum >= 1.0) {
        // Scale down proportions if sum >= 1.0
        double scale = 0.99 / sum;
        for (int i = 0; i < em_data->n_genomes - 1; i++) {
            em_data->proportions[i] *= scale;
        }
        em_data->proportions[em_data->n_genomes - 1] = 0.01;
    } else {
        em_data->proportions[em_data->n_genomes - 1] = 1.0 - sum;
    }
    
    // Bound error and damage rates
    em_data->error_rate = fmax(0.000001, fmin(0.99, par[em_data->n_genomes - 1]));
    em_data->damage_rate = fmax(1e-6, fmin(0.99, par[em_data->n_genomes]));
    
    // Perform one EM iteration
    calcQ_method_CED(em_data);
    em_e_step(em_data);
    em_m_step_proportions(em_data);
    em_data->error_rate = em_m_step_background_rate(em_data);
    em_data->damage_rate = em_m_step_damage_rate(em_data);
    
    // Apply bounds to M-step results (critical for SQUAREM stability)
    em_data->error_rate = fmax(0.000001, fmin(0.99, em_data->error_rate));
    em_data->damage_rate = fmax(1e-6, fmin(0.99, em_data->damage_rate));
    
    // Return new parameters
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        fpar[i] = em_data->proportions[i];
    }
    fpar[em_data->n_genomes - 1] = em_data->error_rate;
    fpar[em_data->n_genomes] = em_data->damage_rate;
}

double em_objective_function_CED(double *par, void *data) {
    em_squarem_data_t *em_sq_data = (em_squarem_data_t*)data;
    em_data_t *em_data = em_sq_data->data;
    
    // Set proportions, error rate, and damage rate (like method BE)
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        em_data->proportions[i] = par[i];
    }
    em_data->proportions[em_data->n_genomes - 1] = 1.0;
    for (int i = 0; i < em_data->n_genomes - 1; i++) {
        em_data->proportions[em_data->n_genomes - 1] -= par[i];
    }
    em_data->error_rate = par[em_data->n_genomes - 1];
    em_data->damage_rate = par[em_data->n_genomes];
    // Ensure damage rate doesn't go below minimum to avoid log(0)
    if (em_data->damage_rate < 1e-6) em_data->damage_rate = 1e-6;
    
    // Check parameter bounds
    for (int i = 0; i < em_data->n_genomes; i++) {
        if (em_data->proportions[i] <= 0.0 || em_data->proportions[i] >= 1.0) {
            return 1e10;  // Penalty for invalid proportions
        }
    }
    if (em_data->error_rate <= 0.0 || em_data->error_rate >= 1.0) {
        return 1e10;  // Penalty for invalid error rate
    }
    // Use minimum of 1e-6 to avoid numerical issues with log(0)
    if (em_data->damage_rate < 1e-6 || em_data->damage_rate > 1.0) {
        return 1e10;  // Penalty for invalid damage rate
    }
    
    // Calculate Q matrix and log-likelihood
    calcQ_method_CED(em_data);
    return -calculate_log_likelihood(em_data);
}

// High-level EM with SQUAREM interface for Method A
em_results_t* em_with_squarem_A(em_data_t *data, em_config_t *config) {
    print_verbose(config->verbose, "Running Method A with SQUAREM acceleration\n");
    
    // Initialize proportions
    initialize_proportions_data_driven(data);
    
    // Setup pruning if enabled
    setup_pruning_if_needed(data, config);
    
    // Set up SQUAREM data
    em_squarem_data_t em_sq_data = {data, config};
    
    // Set up initial parameters (proportions only, excluding last one)
    int n_params = data->n_genomes - 1;
    double *initial_params = vector_create(n_params);
    for (int i = 0; i < n_params; i++) {
        initial_params[i] = data->proportions[i];
    }
    
    // Set up SQUAREM control
    squarem_control_t *sq_control = squarem_control_create();
    sq_control->maxiter = config->max_iterations;
    sq_control->tol = config->tolerance;
    sq_control->trace = config->verbose;
    
    // Run SQUAREM
    squarem_result_t *sq_result = squarem(initial_params, em_fixpoint_function_A, 
                                         em_objective_function_A, &em_sq_data, sq_control);
    
    // Create EM results
    em_results_t *results = em_results_create(data->n_genomes);
    if (results && sq_result) {
        // Set final proportions
        for (int i = 0; i < n_params; i++) {
            results->final_proportions[i] = sq_result->par[i];
        }
        results->final_proportions[data->n_genomes - 1] = 1.0;
        for (int i = 0; i < n_params; i++) {
            results->final_proportions[data->n_genomes - 1] -= sq_result->par[i];
        }
        
        results->final_error_rate = data->error_rate;  // Fixed rate
        results->final_log_likelihood = -sq_result->value;
        results->converged = (sq_result->convergence == 0);
        results->iterations = sq_result->iter;
        
        print_verbose(config->verbose, "Converged after %d iterations\n", sq_result->iter);
    }
    
    // Cleanup
    vector_destroy(initial_params);
    squarem_control_destroy(sq_control);
    squarem_result_destroy(sq_result);
    
    return results;
}

// Stub implementations for missing SQUAREM functions (required for linking)
em_results_t* em_with_squarem_AE(em_data_t *data, em_config_t *config) {
    print_verbose(config->verbose, "Running Method AE with SQUAREM acceleration\n");
    
    // Initialize proportions
    initialize_proportions_data_driven(data);
    
    // Setup pruning if enabled
    setup_pruning_if_needed(data, config);
    
    // Set up SQUAREM data
    em_squarem_data_t em_sq_data = {data, config};
    
    // Set up initial parameters (proportions + error rate)
    int n_params = data->n_genomes;  // Including error rate
    double *initial_params = vector_create(n_params);
    for (int i = 0; i < data->n_genomes - 1; i++) {
        initial_params[i] = data->proportions[i];
    }
    initial_params[data->n_genomes - 1] = data->error_rate;
    
    // Set up SQUAREM control
    squarem_control_t *sq_control = squarem_control_create();
    sq_control->maxiter = config->max_iterations;
    sq_control->tol = config->tolerance;
    sq_control->trace = config->verbose;
    
    // Run SQUAREM
    squarem_result_t *sq_result = squarem(initial_params, em_fixpoint_function_AE, 
                                         em_objective_function_AE, &em_sq_data, sq_control);
    
    // Create EM results
    em_results_t *results = em_results_create(data->n_genomes);
    if (results && sq_result) {
        // Set final proportions
        for (int i = 0; i < data->n_genomes - 1; i++) {
            results->final_proportions[i] = sq_result->par[i];
        }
        results->final_proportions[data->n_genomes - 1] = 1.0;
        for (int i = 0; i < data->n_genomes - 1; i++) {
            results->final_proportions[data->n_genomes - 1] -= sq_result->par[i];
        }
        
        results->final_error_rate = sq_result->par[data->n_genomes - 1];
        results->final_log_likelihood = -sq_result->value;
        results->converged = (sq_result->convergence == 0);
        results->iterations = sq_result->iter;
        
        print_verbose(config->verbose, "Converged after %d iterations\n", sq_result->iter);
    }
    
    // Cleanup
    vector_destroy(initial_params);
    squarem_control_destroy(sq_control);
    squarem_result_destroy(sq_result);
    
    return results;
}
em_results_t* em_with_squarem_B(em_data_t *data, em_config_t *config) {
    print_verbose(config->verbose, "Running Method B with SQUAREM acceleration\n");
    
    // Initialize proportions
    initialize_proportions_data_driven(data);
    
    // Setup pruning if enabled
    setup_pruning_if_needed(data, config);
    
    // Set up SQUAREM data
    em_squarem_data_t em_sq_data = {data, config};
    
    // Set up initial parameters (proportions only, excluding last one)
    int n_params = data->n_genomes - 1;
    double *initial_params = vector_create(n_params);
    for (int i = 0; i < n_params; i++) {
        initial_params[i] = data->proportions[i];
    }
    
    // Set up SQUAREM control
    squarem_control_t *sq_control = squarem_control_create();
    sq_control->maxiter = config->max_iterations;
    sq_control->tol = config->tolerance;
    sq_control->trace = config->verbose;
    
    // Run SQUAREM
    squarem_result_t *sq_result = squarem(initial_params, em_fixpoint_function_B, 
                                         em_objective_function_B, &em_sq_data, sq_control);
    
    // Create EM results
    em_results_t *results = em_results_create(data->n_genomes);
    if (results && sq_result) {
        // Set final proportions
        for (int i = 0; i < n_params; i++) {
            results->final_proportions[i] = sq_result->par[i];
        }
        results->final_proportions[data->n_genomes - 1] = 1.0;
        for (int i = 0; i < n_params; i++) {
            results->final_proportions[data->n_genomes - 1] -= sq_result->par[i];
        }
        
        results->final_error_rate = data->error_rate;  // Fixed rate
        results->final_log_likelihood = -sq_result->value;
        results->converged = (sq_result->convergence == 0);
        results->iterations = sq_result->iter;
        
        print_verbose(config->verbose, "Converged after %d iterations\n", sq_result->iter);
    }
    
    // Cleanup
    vector_destroy(initial_params);
    squarem_control_destroy(sq_control);
    squarem_result_destroy(sq_result);
    
    return results;
}
em_results_t* em_with_squarem_BE(em_data_t *data, em_config_t *config) {
    print_verbose(config->verbose, "Running Method BE with SQUAREM acceleration\n");
    
    // Initialize proportions
    initialize_proportions_data_driven(data);
    
    // Setup pruning if enabled
    setup_pruning_if_needed(data, config);
    
    // Set up SQUAREM data
    em_squarem_data_t em_sq_data = {data, config};
    
    // Set up initial parameters (proportions + error rate)
    int n_params = data->n_genomes;  // Including error rate
    double *initial_params = vector_create(n_params);
    for (int i = 0; i < data->n_genomes - 1; i++) {
        initial_params[i] = data->proportions[i];
    }
    initial_params[data->n_genomes - 1] = data->error_rate;
    
    // Set up SQUAREM control
    squarem_control_t *sq_control = squarem_control_create();
    sq_control->maxiter = config->max_iterations;
    sq_control->tol = config->tolerance;
    sq_control->trace = config->verbose;
    
    // Run SQUAREM
    squarem_result_t *sq_result = squarem(initial_params, em_fixpoint_function_BE, 
                                         em_objective_function_BE, &em_sq_data, sq_control);
    
    // Create EM results
    em_results_t *results = em_results_create(data->n_genomes);
    if (results && sq_result) {
        // Set final proportions
        for (int i = 0; i < data->n_genomes - 1; i++) {
            results->final_proportions[i] = sq_result->par[i];
        }
        results->final_proportions[data->n_genomes - 1] = 1.0;
        for (int i = 0; i < data->n_genomes - 1; i++) {
            results->final_proportions[data->n_genomes - 1] -= sq_result->par[i];
        }
        
        results->final_error_rate = sq_result->par[data->n_genomes - 1];
        results->final_log_likelihood = -sq_result->value;
        results->converged = (sq_result->convergence == 0);
        results->iterations = sq_result->iter;
        
        print_verbose(config->verbose, "Converged after %d iterations\n", sq_result->iter);
    }
    
    // Cleanup
    vector_destroy(initial_params);
    squarem_control_destroy(sq_control);
    squarem_result_destroy(sq_result);
    
    return results;
}
em_results_t* em_with_squarem_BEfull(em_data_t *data, em_config_t *config) {
    print_verbose(config->verbose, "Running Method BEfull with SQUAREM acceleration\n");
    
    // Initialize proportions and per-genome error rates
    initialize_proportions_data_driven(data);
    
    // Initialize per-genome error rates array
    if (!data->error_rates) {
        data->error_rates = vector_create(data->n_genomes);
        if (!data->error_rates) {
            return NULL;
        }
    }
    
    // Initialize per-genome error rates to default value
    for (int j = 0; j < data->n_genomes; j++) {
        data->error_rates[j] = config->error_rate;
    }
    
    // Set up SQUAREM data
    em_squarem_data_t em_sq_data = {data, config};
    
    // Set up initial parameters: [prop_0, ..., prop_{K-2}, error_0, ..., error_{K-1}]
    int n_params = 2 * data->n_genomes - 1;
    double *initial_params = vector_create(n_params);
    
    // Pack proportions (first K-1)
    for (int i = 0; i < data->n_genomes - 1; i++) {
        initial_params[i] = data->proportions[i];
    }
    
    // Pack per-genome error rates (next K)
    for (int j = 0; j < data->n_genomes; j++) {
        initial_params[data->n_genomes - 1 + j] = data->error_rates[j];
    }
    
    // Set up SQUAREM control with default settings
    squarem_control_t *sq_control = squarem_control_create();
    sq_control->maxiter = config->max_iterations;
    sq_control->tol = config->tolerance;
    sq_control->trace = config->verbose;
    
    // Run SQUAREM
    squarem_result_t *sq_result = squarem(initial_params, em_fixpoint_function_BEfull, 
                                         em_objective_function_BEfull, &em_sq_data, sq_control);
    
    // Create EM results
    em_results_t *results = em_results_create(data->n_genomes);
    if (results && sq_result) {
        // Unpack final proportions
        for (int i = 0; i < data->n_genomes - 1; i++) {
            results->final_proportions[i] = sq_result->par[i];
        }
        results->final_proportions[data->n_genomes - 1] = 1.0;
        for (int i = 0; i < data->n_genomes - 1; i++) {
            results->final_proportions[data->n_genomes - 1] -= sq_result->par[i];
        }
        
        // Unpack final per-genome error rates
        if (!results->final_error_rates) {
            results->final_error_rates = vector_create(data->n_genomes);
        }
        for (int j = 0; j < data->n_genomes; j++) {
            results->final_error_rates[j] = sq_result->par[data->n_genomes - 1 + j];
        }
        
        results->final_log_likelihood = -sq_result->value;
        results->converged = (sq_result->convergence == 0);
        results->iterations = sq_result->iter;
        
        print_verbose(config->verbose, "Converged after %d iterations\n", sq_result->iter);
    }
    
    // Cleanup
    vector_destroy(initial_params);
    squarem_control_destroy(sq_control);
    squarem_result_destroy(sq_result);
    
    return results;
}
em_results_t* em_with_squarem_C(em_data_t *data, em_config_t *config) {
    printf("DEBUG: em_with_squarem_C() called - dense SQUAREM implementation reached\n");
    print_verbose(config->verbose, "Running Method C with SQUAREM acceleration\n");
    
    // Initialize proportions
    initialize_proportions_data_driven(data);
    
    // Setup pruning if enabled
    setup_pruning_if_needed(data, config);
    
    // Set up SQUAREM data
    em_squarem_data_t em_sq_data = {data, config};
    
    // Set up initial parameters (proportions only, excluding last one)
    int n_params = data->n_genomes - 1;
    double *initial_params = vector_create(n_params);
    for (int i = 0; i < n_params; i++) {
        initial_params[i] = data->proportions[i];
    }
    
    // Set up SQUAREM control
    squarem_control_t *sq_control = squarem_control_create();
    sq_control->maxiter = config->max_iterations;
    sq_control->tol = config->tolerance;
    sq_control->trace = config->verbose;
    
    // Run SQUAREM
    squarem_result_t *sq_result = squarem(initial_params, em_fixpoint_function_C, 
                                         em_objective_function_C, &em_sq_data, sq_control);
    
    // Create EM results
    em_results_t *results = em_results_create(data->n_genomes);
    if (results && sq_result) {
        // Set final proportions
        for (int i = 0; i < n_params; i++) {
            results->final_proportions[i] = sq_result->par[i];
        }
        results->final_proportions[data->n_genomes - 1] = 1.0;
        for (int i = 0; i < n_params; i++) {
            results->final_proportions[data->n_genomes - 1] -= sq_result->par[i];
        }
        
        results->final_error_rate = data->error_rate;  // Fixed rate
        results->final_damage_rate = data->damage_rate;  // Fixed rate
        results->final_log_likelihood = -sq_result->value;
        results->converged = (sq_result->convergence == 0);
        results->iterations = sq_result->iter;
        
        print_verbose(config->verbose, "Converged after %d iterations\n", sq_result->iter);
    }
    
    // Cleanup
    vector_destroy(initial_params);
    squarem_control_destroy(sq_control);
    squarem_result_destroy(sq_result);
    
    return results;
}
em_results_t* em_with_squarem_CE(em_data_t *data, em_config_t *config) {
    print_verbose(config->verbose, "Running Method CE with SQUAREM acceleration\n");
    
    // Initialize proportions
    initialize_proportions_data_driven(data);
    
    // Setup pruning if enabled
    setup_pruning_if_needed(data, config);
    
    // For CE method, ensure we have a reasonable starting error rate
    // Check if error_rate wasn't set or is invalid
    if (data->error_rate <= 0.0 || data->error_rate >= 1.0) {
        // Use config error rate if valid, otherwise default
        if (config->error_rate > 0.0 && config->error_rate < 1.0) {
            data->error_rate = config->error_rate;
        } else {
            data->error_rate = 0.001;  // Default starting value for background error rate
        }
    }
    
    // Set up SQUAREM data
    em_squarem_data_t em_sq_data = {data, config};
    
    // Set up initial parameters (proportions + error rate)
    int n_params = data->n_genomes;  // Including error rate
    double *initial_params = vector_create(n_params);
    for (int i = 0; i < data->n_genomes - 1; i++) {
        initial_params[i] = data->proportions[i];
    }
    initial_params[data->n_genomes - 1] = data->error_rate;
    
    // Set up SQUAREM control
    squarem_control_t *sq_control = squarem_control_create();
    sq_control->maxiter = config->max_iterations;
    sq_control->tol = config->tolerance;
    sq_control->trace = config->verbose;
    
    // Run SQUAREM
    squarem_result_t *sq_result = squarem(initial_params, em_fixpoint_function_CE, 
                                         em_objective_function_CE, &em_sq_data, sq_control);
    
    // Create EM results
    em_results_t *results = em_results_create(data->n_genomes);
    if (results && sq_result) {
        // Set final proportions
        for (int i = 0; i < data->n_genomes - 1; i++) {
            results->final_proportions[i] = sq_result->par[i];
        }
        results->final_proportions[data->n_genomes - 1] = 1.0;
        for (int i = 0; i < data->n_genomes - 1; i++) {
            results->final_proportions[data->n_genomes - 1] -= sq_result->par[i];
        }
        
        results->final_error_rate = sq_result->par[data->n_genomes - 1];
        results->final_damage_rate = data->damage_rate;  // Fixed rate
        results->final_log_likelihood = -sq_result->value;
        results->converged = (sq_result->convergence == 0);
        results->iterations = sq_result->iter;
        
        print_verbose(config->verbose, "Converged after %d iterations\n", sq_result->iter);
    }
    
    // Cleanup
    vector_destroy(initial_params);
    squarem_control_destroy(sq_control);
    squarem_result_destroy(sq_result);
    
    return results;
}
em_results_t* em_with_squarem_CED(em_data_t *data, em_config_t *config) {
    print_verbose(config->verbose, "Running Method CED with SQUAREM acceleration\n");
    
    // Initialize proportions
    initialize_proportions_data_driven(data);
    
    // Setup pruning if enabled
    setup_pruning_if_needed(data, config);
    
    // For CED method, ensure we have reasonable starting rates
    if (data->error_rate <= 0.0 || data->error_rate >= 1.0) {
        data->error_rate = 0.001;  // Default starting value for background error rate
    }
    // Use minimum of 1e-6 to avoid numerical issues with log(0)
    if (data->damage_rate < 1e-6 || data->damage_rate >= 1.0) {
        data->damage_rate = 0.01;  // Default starting value for damage rate
    }
    
    // Set up SQUAREM data
    em_squarem_data_t em_sq_data = {data, config};
    
    // Set up initial parameters (proportions + error rate + damage rate)
    int n_params = data->n_genomes + 1;  // Including error and damage rates
    double *initial_params = vector_create(n_params);
    for (int i = 0; i < data->n_genomes - 1; i++) {
        initial_params[i] = data->proportions[i];
    }
    initial_params[data->n_genomes - 1] = data->error_rate;
    initial_params[data->n_genomes] = data->damage_rate;
    
    // Set up SQUAREM control
    squarem_control_t *sq_control = squarem_control_create();
    sq_control->maxiter = config->max_iterations;
    sq_control->tol = config->tolerance;
    sq_control->trace = config->verbose;
    
    // Run SQUAREM
    squarem_result_t *sq_result = squarem(initial_params, em_fixpoint_function_CED, 
                                         em_objective_function_CED, &em_sq_data, sq_control);
    
    // Create EM results
    em_results_t *results = em_results_create(data->n_genomes);
    if (results && sq_result) {
        // Set final proportions
        for (int i = 0; i < data->n_genomes - 1; i++) {
            results->final_proportions[i] = sq_result->par[i];
        }
        results->final_proportions[data->n_genomes - 1] = 1.0;
        for (int i = 0; i < data->n_genomes - 1; i++) {
            results->final_proportions[data->n_genomes - 1] -= sq_result->par[i];
        }
        
        results->final_error_rate = sq_result->par[data->n_genomes - 1];
        results->final_damage_rate = sq_result->par[data->n_genomes];
        results->final_log_likelihood = -sq_result->value;
        results->converged = (sq_result->convergence == 0);
        results->iterations = sq_result->iter;
        
        print_verbose(config->verbose, "Converged after %d iterations\n", sq_result->iter);
    }
    
    // Cleanup
    vector_destroy(initial_params);
    squarem_control_destroy(sq_control);
    squarem_result_destroy(sq_result);
    
    return results;
}
em_results_t* em_with_squarem_CEDfull(em_data_t *data, em_config_t *config) {
    print_verbose(config->verbose, "Running Method CEDfull with SQUAREM acceleration\n");
    
    // Initialize proportions, per-genome error rates, and per-genome damage rates
    initialize_proportions_data_driven(data);
    
    // Initialize per-genome error rates array
    if (!data->error_rates) {
        data->error_rates = vector_create(data->n_genomes);
        if (!data->error_rates) {
            return NULL;
        }
    }
    
    // Initialize per-genome damage rates array
    if (!data->damage_rates) {
        data->damage_rates = vector_create(data->n_genomes);
        if (!data->damage_rates) {
            return NULL;
        }
    }
    
    // Initialize per-genome error rates and damage rates to default values
    for (int j = 0; j < data->n_genomes; j++) {
        data->error_rates[j] = config->error_rate;
        data->damage_rates[j] = config->damage_rate;
    }
    
    // Set up SQUAREM data
    em_squarem_data_t em_sq_data = {data, config};
    
    // Set up initial parameters: [prop_0, ..., prop_{K-2}, error_0, ..., error_{K-1}, damage_0, ..., damage_{K-1}]
    int n_params = 3 * data->n_genomes - 1;
    double *initial_params = vector_create(n_params);
    
    // Pack proportions (first K-1)
    for (int i = 0; i < data->n_genomes - 1; i++) {
        initial_params[i] = data->proportions[i];
    }
    
    // Pack per-genome error rates (next K)
    for (int j = 0; j < data->n_genomes; j++) {
        initial_params[data->n_genomes - 1 + j] = data->error_rates[j];
    }
    
    // Pack per-genome damage rates (last K)
    for (int j = 0; j < data->n_genomes; j++) {
        initial_params[2 * data->n_genomes - 1 + j] = data->damage_rates[j];
    }
    
    // Set up SQUAREM control with default settings
    squarem_control_t *sq_control = squarem_control_create();
    sq_control->maxiter = config->max_iterations;
    sq_control->tol = config->tolerance;
    sq_control->trace = config->verbose;
    
    // Run SQUAREM
    squarem_result_t *sq_result = squarem(initial_params, em_fixpoint_function_CEDfull, 
                                         em_objective_function_CEDfull, &em_sq_data, sq_control);
    
    // Create EM results
    em_results_t *results = em_results_create(data->n_genomes);
    if (results && sq_result) {
        // Unpack final proportions
        for (int i = 0; i < data->n_genomes - 1; i++) {
            results->final_proportions[i] = sq_result->par[i];
        }
        results->final_proportions[data->n_genomes - 1] = 1.0;
        for (int i = 0; i < data->n_genomes - 1; i++) {
            results->final_proportions[data->n_genomes - 1] -= sq_result->par[i];
        }
        
        // Unpack final per-genome error rates
        if (!results->final_error_rates) {
            results->final_error_rates = vector_create(data->n_genomes);
        }
        for (int j = 0; j < data->n_genomes; j++) {
            results->final_error_rates[j] = sq_result->par[data->n_genomes - 1 + j];
        }
        
        // Unpack final per-genome damage rates
        if (!results->final_damage_rates) {
            results->final_damage_rates = vector_create(data->n_genomes);
        }
        for (int j = 0; j < data->n_genomes; j++) {
            results->final_damage_rates[j] = sq_result->par[2 * data->n_genomes - 1 + j];
        }
        
        results->final_log_likelihood = -sq_result->value;
        results->converged = (sq_result->convergence == 0);
        results->iterations = sq_result->iter;
        
        print_verbose(config->verbose, "Converged after %d iterations\n", sq_result->iter);
    }
    
    // Cleanup
    vector_destroy(initial_params);
    squarem_control_destroy(sq_control);
    squarem_result_destroy(sq_result);
    
    return results;
}
em_results_t* em_with_squarem_CEfull(em_data_t *data, em_config_t *config) {
    print_verbose(config->verbose, "Running Method CEfull with SQUAREM acceleration\n");
    
    // Initialize proportions and per-genome error rates
    initialize_proportions_data_driven(data);
    
    // Initialize per-genome error rates array
    if (!data->error_rates) {
        data->error_rates = vector_create(data->n_genomes);
        if (!data->error_rates) {
            return NULL;
        }
    }
    
    // Initialize per-genome error rates to default value
    for (int j = 0; j < data->n_genomes; j++) {
        data->error_rates[j] = config->error_rate;
    }
    
    // Set up SQUAREM data
    em_squarem_data_t em_sq_data = {data, config};
    
    // Set up initial parameters: [prop_0, ..., prop_{K-2}, error_0, ..., error_{K-1}, damage_rate]
    int n_params = 2 * data->n_genomes;
    double *initial_params = vector_create(n_params);
    
    // Pack proportions (first K-1)
    for (int i = 0; i < data->n_genomes - 1; i++) {
        initial_params[i] = data->proportions[i];
    }
    
    // Pack per-genome error rates (next K)
    for (int j = 0; j < data->n_genomes; j++) {
        initial_params[data->n_genomes - 1 + j] = data->error_rates[j];
    }
    
    // Pack damage rate (last parameter)
    initial_params[2 * data->n_genomes - 1] = data->damage_rate;
    
    // Set up SQUAREM control with default settings
    squarem_control_t *sq_control = squarem_control_create();
    sq_control->maxiter = config->max_iterations;
    sq_control->tol = config->tolerance;
    sq_control->trace = config->verbose;
    
    // Run SQUAREM
    squarem_result_t *sq_result = squarem(initial_params, em_fixpoint_function_CEfull, 
                                         em_objective_function_CEfull, &em_sq_data, sq_control);
    
    // Create EM results
    em_results_t *results = em_results_create(data->n_genomes);
    if (results && sq_result) {
        // Unpack final proportions
        for (int i = 0; i < data->n_genomes - 1; i++) {
            results->final_proportions[i] = sq_result->par[i];
        }
        results->final_proportions[data->n_genomes - 1] = 1.0;
        for (int i = 0; i < data->n_genomes - 1; i++) {
            results->final_proportions[data->n_genomes - 1] -= sq_result->par[i];
        }
        
        // Unpack final per-genome error rates
        if (!results->final_error_rates) {
            results->final_error_rates = vector_create(data->n_genomes);
        }
        for (int j = 0; j < data->n_genomes; j++) {
            results->final_error_rates[j] = sq_result->par[data->n_genomes - 1 + j];
        }
        
        // Unpack final damage rate
        results->final_damage_rate = sq_result->par[2 * data->n_genomes - 1];
        
        results->final_log_likelihood = -sq_result->value;
        results->converged = (sq_result->convergence == 0);
        results->iterations = sq_result->iter;
        
        print_verbose(config->verbose, "Converged after %d iterations\n", sq_result->iter);
    }
    
    // Cleanup
    vector_destroy(initial_params);
    squarem_control_destroy(sq_control);
    squarem_result_destroy(sq_result);
    
    return results;
}


// CEDfull fixpoint function
void em_fixpoint_function_CEDfull(double *par, double *fpar, void *data) {
    printf("DEBUG: em_fixpoint_function_CEDfull() called - dense CEDfull SQUAREM reached\n");
    em_squarem_data_t *em_sq_data = (em_squarem_data_t*)data;
    em_data_t *em_data = em_sq_data->data;
    int K = em_data->n_genomes;
    
    // Unpack parameters: first K-1 are proportions, next K are error rates, last K are damage rates
    for (int i = 0; i < K - 1; i++) {
        em_data->proportions[i] = fmax(0.0, fmin(1.0, par[i]));
    }
    
    // Calculate last proportion with overflow protection
    double sum = 0.0;
    for (int i = 0; i < K - 1; i++) {
        sum += em_data->proportions[i];
    }
    
    if (sum >= 1.0) {
        double scale = 0.99 / sum;
        for (int i = 0; i < K - 1; i++) {
            em_data->proportions[i] *= scale;
        }
        em_data->proportions[K - 1] = 0.01;
    } else {
        em_data->proportions[K - 1] = 1.0 - sum;
    }
    
    // Set per-genome error and damage rates
    for (int j = 0; j < K; j++) {
        em_data->error_rates[j] = fmax(0.000001, fmin(0.99, par[K - 1 + j]));
        em_data->damage_rates[j] = fmax(0.0, fmin(0.99, par[2 * K - 1 + j]));
    }
    
    // One EM step
    calcQ_method_CEDfull(em_data);
    em_e_step(em_data);
    em_m_step_proportions(em_data);
    em_m_step_per_genome_error_rates(em_data);
    em_m_step_per_genome_damage_rates(em_data);
    
    // Pack new parameters
    for (int i = 0; i < K - 1; i++) {
        fpar[i] = em_data->proportions[i];
    }
    for (int j = 0; j < K; j++) {
        fpar[K - 1 + j] = em_data->error_rates[j];
        fpar[2 * K - 1 + j] = em_data->damage_rates[j];
    }
}

// CEDfull objective function
double em_objective_function_CEDfull(double *par, void *data) {
    em_squarem_data_t *em_sq_data = (em_squarem_data_t*)data;
    em_data_t *em_data = em_sq_data->data;
    int K = em_data->n_genomes;
    
    // Unpack parameters
    for (int i = 0; i < K - 1; i++) {
        em_data->proportions[i] = fmax(0.0, fmin(1.0, par[i]));
    }
    
    // Calculate last proportion with overflow protection
    double sum = 0.0;
    for (int i = 0; i < K - 1; i++) {
        sum += em_data->proportions[i];
    }
    
    if (sum >= 1.0) {
        double scale = 0.99 / sum;
        for (int i = 0; i < K - 1; i++) {
            em_data->proportions[i] *= scale;
        }
        em_data->proportions[K - 1] = 0.01;
    } else {
        em_data->proportions[K - 1] = 1.0 - sum;
    }
    
    // Set per-genome error rates and damage rates
    for (int j = 0; j < K; j++) {
        em_data->error_rates[j] = fmax(0.000001, fmin(0.99, par[K - 1 + j]));
        em_data->damage_rates[j] = fmax(0.0, fmin(0.99, par[2 * K - 1 + j]));
    }
    
    // Calculate Q matrix and log-likelihood
    calcQ_method_CEDfull(em_data);
    return -calculate_log_likelihood(em_data);
}

// CEfull fixpoint function
void em_fixpoint_function_CEfull(double *par, double *fpar, void *data) {
    em_squarem_data_t *em_sq_data = (em_squarem_data_t*)data;
    em_data_t *em_data = em_sq_data->data;
    int K = em_data->n_genomes;
    
    // Unpack parameters: first K-1 are proportions, next K are error rates, last is damage rate
    for (int i = 0; i < K - 1; i++) {
        em_data->proportions[i] = fmax(0.0, fmin(1.0, par[i]));
    }
    
    // Calculate last proportion with overflow protection
    double sum = 0.0;
    for (int i = 0; i < K - 1; i++) {
        sum += em_data->proportions[i];
    }
    
    if (sum >= 1.0) {
        double scale = 0.99 / sum;
        for (int i = 0; i < K - 1; i++) {
            em_data->proportions[i] *= scale;
        }
        em_data->proportions[K - 1] = 0.01;
    } else {
        em_data->proportions[K - 1] = 1.0 - sum;
    }
    
    // Set per-genome error rates
    for (int j = 0; j < K; j++) {
        em_data->error_rates[j] = fmax(0.000001, fmin(0.99, par[K - 1 + j]));
    }
    
    // Set damage rate (last parameter)
    em_data->damage_rate = fmax(0.0, fmin(0.99, par[2 * K - 1]));
    
    // One EM step
    calcQ_method_CEfull(em_data);
    em_e_step(em_data);
    em_m_step_proportions(em_data);
    em_m_step_per_genome_error_rates(em_data);
    
    // Pack new parameters
    for (int i = 0; i < K - 1; i++) {
        fpar[i] = em_data->proportions[i];
    }
    for (int j = 0; j < K; j++) {
        fpar[K - 1 + j] = em_data->error_rates[j];
    }
    fpar[2 * K - 1] = em_data->damage_rate;
}

// CEfull objective function
double em_objective_function_CEfull(double *par, void *data) {
    em_squarem_data_t *em_sq_data = (em_squarem_data_t*)data;
    em_data_t *em_data = em_sq_data->data;
    int K = em_data->n_genomes;
    
    // Unpack parameters
    for (int i = 0; i < K - 1; i++) {
        em_data->proportions[i] = fmax(0.0, fmin(1.0, par[i]));
    }
    
    // Calculate last proportion with overflow protection
    double sum = 0.0;
    for (int i = 0; i < K - 1; i++) {
        sum += em_data->proportions[i];
    }
    
    if (sum >= 1.0) {
        double scale = 0.99 / sum;
        for (int i = 0; i < K - 1; i++) {
            em_data->proportions[i] *= scale;
        }
        em_data->proportions[K - 1] = 0.01;
    } else {
        em_data->proportions[K - 1] = 1.0 - sum;
    }
    
    // Set per-genome error rates and damage rate
    for (int j = 0; j < K; j++) {
        em_data->error_rates[j] = fmax(0.000001, fmin(0.99, par[K - 1 + j]));
    }
    em_data->damage_rate = fmax(0.0, fmin(0.99, par[2 * K - 1]));
    
    // Calculate Q matrix and log-likelihood
    calcQ_method_CEfull(em_data);
    return -calculate_log_likelihood(em_data);
}


// BEfull fixpoint function
void em_fixpoint_function_BEfull(double *par, double *fpar, void *data) {
    em_squarem_data_t *em_sq_data = (em_squarem_data_t*)data;
    em_data_t *em_data = em_sq_data->data;
    int K = em_data->n_genomes;
    
    // Set proportions: first K-1 from parameters, last one calculated  
    for (int i = 0; i < K - 1; i++) {
        em_data->proportions[i] = par[i];
    }
    em_data->proportions[K - 1] = 1.0;
    for (int i = 0; i < K - 1; i++) {
        em_data->proportions[K - 1] -= par[i];
    }
    
    // Set per-genome error rates
    for (int j = 0; j < K; j++) {
        em_data->error_rates[j] = par[K - 1 + j];
    }
    
    // One EM step
    calcQ_method_BEfull(em_data);
    em_e_step(em_data);
    em_m_step_proportions(em_data);
    em_m_step_per_genome_error_rates(em_data);
    
    // Pack new parameters
    for (int i = 0; i < K - 1; i++) {
        fpar[i] = em_data->proportions[i];
    }
    for (int j = 0; j < K; j++) {
        fpar[K - 1 + j] = em_data->error_rates[j];
    }
}

// BEfull objective function
double em_objective_function_BEfull(double *par, void *data) {
    em_squarem_data_t *em_sq_data = (em_squarem_data_t*)data;
    em_data_t *em_data = em_sq_data->data;
    int K = em_data->n_genomes;
    
    // Set proportions and per-genome error rates
    for (int i = 0; i < K - 1; i++) {
        em_data->proportions[i] = par[i];
    }
    em_data->proportions[K - 1] = 1.0;
    for (int i = 0; i < K - 1; i++) {
        em_data->proportions[K - 1] -= par[i];
    }
    
    for (int j = 0; j < K; j++) {
        em_data->error_rates[j] = par[K - 1 + j];
    }
    
    // Check bounds
    for (int i = 0; i < K; i++) {
        if (em_data->proportions[i] <= 0.0 || em_data->proportions[i] >= 1.0) {
            return 1e10;
        }
    }
    for (int j = 0; j < K; j++) {
        if (em_data->error_rates[j] <= 0.0 || em_data->error_rates[j] >= 1.0) {
            return 1e10;
        }
    }
    
    calcQ_method_BEfull(em_data);
    return -calculate_log_likelihood(em_data);
}
