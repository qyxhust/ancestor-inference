#ifndef SQUAREM_H
#define SQUAREM_H

#include "em_types.h"

// SQUAREM acceleration parameters
typedef struct {
    double *theta_0;        // Starting parameters
    double *theta_1;        // After first EM step
    double *theta_2;        // After second EM step
    double *r_vec;          // θ_1 - θ_0
    double *v_vec;          // θ_2 - θ_1 - r
    double *theta_new;      // Extrapolated parameters
    int n_params;           // Number of parameters
    double step_min;        // Minimum step size
    double step_max;        // Maximum step size
    double mstep;           // Step size parameter
    double objfn_inc;       // Objective function increase tolerance
} squarem_t;

// SQUAREM control parameters
typedef struct {
    int maxiter;            // Maximum iterations
    double tol;             // Convergence tolerance
    int trace;              // Verbose output flag
    double step_min0;       // Initial minimum step size
    double step_max0;       // Initial maximum step size
    double mstep;           // Step size multiplier
    double objfn_inc;       // Objective function increase tolerance
    int kr;                 // Number of steps before extrapolation
} squarem_control_t;

// SQUAREM result structure
typedef struct {
    double *par;            // Final parameters
    double value;           // Final objective function value
    int iter;               // Number of iterations
    int convergence;        // Convergence flag (0 = success)
    char *message;          // Convergence message
} squarem_result_t;

// Function pointer types for SQUAREM
typedef void (*fixptfn_t)(double *par, double *fpar, void *data);
typedef double (*objfn_t)(double *par, void *data);

// SQUAREM core functions
squarem_t* squarem_create(int n_params);
void squarem_destroy(squarem_t *sq);
squarem_control_t* squarem_control_create(void);
void squarem_control_destroy(squarem_control_t *control);
squarem_result_t* squarem_result_create(int n_params);
void squarem_result_destroy(squarem_result_t *result);

// Main SQUAREM algorithm
squarem_result_t* squarem(double *par, fixptfn_t fixptfn, objfn_t objfn, 
                         void *data, squarem_control_t *control);

// Utility functions
double vector_norm_squared(double *vec, int n);
double vector_dot_product(double *vec1, double *vec2, int n);
void vector_copy_squarem(double *dest, double *src, int n);
void vector_add(double *result, double *vec1, double *vec2, int n);
void vector_subtract(double *result, double *vec1, double *vec2, int n);
void vector_scale_add(double *result, double *vec, double scale, int n);

// EM-specific SQUAREM wrappers
void em_fixpoint_function_A(double *par, double *fpar, void *data);
void em_fixpoint_function_AE(double *par, double *fpar, void *data);
void em_fixpoint_function_B(double *par, double *fpar, void *data);
void em_fixpoint_function_BE(double *par, double *fpar, void *data);

// Damage model SQUAREM wrappers
void em_fixpoint_function_C(double *par, double *fpar, void *data);
void em_fixpoint_function_CE(double *par, double *fpar, void *data);
void em_fixpoint_function_CED(double *par, double *fpar, void *data);

// Full model SQUAREM wrappers
void em_fixpoint_function_BEfull(double *par, double *fpar, void *data);
void em_fixpoint_function_CEfull(double *par, double *fpar, void *data);
void em_fixpoint_function_CEDfull(double *par, double *fpar, void *data);

double em_objective_function_A(double *par, void *data);
double em_objective_function_AE(double *par, void *data);
double em_objective_function_B(double *par, void *data);
double em_objective_function_BE(double *par, void *data);

// Damage model objective functions
double em_objective_function_C(double *par, void *data);
double em_objective_function_CE(double *par, void *data);
double em_objective_function_CED(double *par, void *data);

// Full model objective functions
double em_objective_function_BEfull(double *par, void *data);
double em_objective_function_CEfull(double *par, void *data);
double em_objective_function_CEDfull(double *par, void *data);

// High-level EM with SQUAREM interface
em_results_t* em_with_squarem_A(em_data_t *data, em_config_t *config);
em_results_t* em_with_squarem_AE(em_data_t *data, em_config_t *config);
em_results_t* em_with_squarem_B(em_data_t *data, em_config_t *config);
em_results_t* em_with_squarem_BE(em_data_t *data, em_config_t *config);

// Damage model high-level SQUAREM interfaces
em_results_t* em_with_squarem_C(em_data_t *data, em_config_t *config);
em_results_t* em_with_squarem_CE(em_data_t *data, em_config_t *config);
em_results_t* em_with_squarem_CED(em_data_t *data, em_config_t *config);

// Full model high-level SQUAREM interfaces
em_results_t* em_with_squarem_BEfull(em_data_t *data, em_config_t *config);
em_results_t* em_with_squarem_CEfull(em_data_t *data, em_config_t *config);
em_results_t* em_with_squarem_CEDfull(em_data_t *data, em_config_t *config);

#endif // SQUAREM_H