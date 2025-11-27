#include "em_algorithms.h"
#include "io_utils.h"
#include "sparse_matrix.h"
#include "sparse_em.h"
#include <time.h>
#include <sys/time.h>

// External debug likelihood variables
extern int debug_likelihood_mode;
extern char *debug_parameter_file;
extern char **debug_genome_names;
extern int debug_n_genomes;
extern void load_debug_parameters(sparse_em_data_t *sparse_em_data);

// External likelihood ratio test variables  
extern int likelihood_ratio_test_mode;
extern char *test_parameter_file;
extern char *test_genome_name;
extern double original_log_likelihood;
extern void load_test_parameters_and_setup(sparse_em_data_t *sparse_em_data);

// Utility function to get current time in seconds
double get_time() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

// Global pointer for qsort comparison function
static double *qsort_proportions = NULL;

// Comparison function for sorting genome indices by proportion (descending)
int compare_genome_indices_desc(const void *a, const void *b) {
    int idx_a = *(const int *)a;
    int idx_b = *(const int *)b;
    
    // Sort in descending order by proportion
    if (qsort_proportions[idx_b] > qsort_proportions[idx_a]) return 1;
    if (qsort_proportions[idx_b] < qsort_proportions[idx_a]) return -1;
    return 0;
}

// Command line options structure
static struct option long_options[] = {
    {"input",       required_argument, 0, 'i'},
    {"method",      required_argument, 0, 'M'},
    {"output",      required_argument, 0, 'o'},
    {"error-rate",  required_argument, 0, 'e'},
    {"damage-rate", required_argument, 0, 'r'},
    {"filter",      required_argument, 0, 'f'},
    {"tolerance",   required_argument, 0, 't'},
    {"max-iter",    required_argument, 0, 'm'},
    {"squarem",     no_argument,       0, 's'},
    {"no-squarem",  no_argument,       0, 'S'},
    {"sparse",      no_argument,       0, 'p'},
    {"dense",       no_argument,       0, 'd'},
    {"verbose",     no_argument,       0, 'v'},
    {"convert-names", no_argument,     0, 'c'},
    {"mapping-file",  required_argument, 0, 'g'},
    {"genome-key",    required_argument, 0, 'k'},
    {"use-genome-key", no_argument,    0, 'K'},
    {"no-genome-key", no_argument,     0, 'N'},
    {"constraints", required_argument, 0, 'C'},
    {"weights", required_argument,     0, 'W'},
    {"joint-damage", no_argument,      0, 'j'},
    {"debug-likelihood", required_argument, 0, 'D'},
    {"help",        no_argument,       0, 'h'},
    {"version",     no_argument,       0, 'V'},
    {0, 0, 0, 0}
};

// Print usage information
void print_usage(const char *program_name) {
    printf("CEM - C EM Algorithm for Mixture Proportion Estimation\n");
    printf("Usage: %s [OPTIONS]\n\n", program_name);
    
    printf("Required Options:\n");
    printf("  -i, --input FILE      Input mismatch matrix file (TSV format)\n");
    printf("  -M, --method METHOD   Method to use: A, AE, B, BE, C, CE, CED, BEfull, CEfull, CEDfull\n\n");
    
    printf("Optional Options:\n");
    printf("  -o, --output PREFIX   Output file prefix (default: input filename)\n");
    printf("  -e, --error-rate RATE Initial error rate (default: 0.005)\n");
    printf("  -r, --damage-rate RATE Initial damage rate (default: 0.05, for damage models)\n");
    printf("  -f, --filter THRESH   Minimum proportion threshold (default: auto, 0 to disable)\n");
    printf("  -t, --tolerance TOL   Convergence tolerance (default: 1e-6)\n");
    printf("  -m, --max-iter N      Maximum iterations (default: 1000)\n");
    printf("  -s, --squarem         Force SQUAREM acceleration (default: enabled)\n");
    printf("  -S, --no-squarem      Disable SQUAREM, use standard EM\n");
    printf("  -p, --sparse          Force sparse data processing\n");
    printf("  -d, --dense           Force dense data processing\n");
    printf("  -v, --verbose         Enable verbose output\n");
    printf("  -c, --convert-names   Enable genome name conversion (requires --mapping-file)\n");
    printf("  -g, --mapping-file FILE Genome mapping file for name conversion\n");
    printf("  -k, --genome-key FILE Specify genome key file for short name translation\n");
    printf("  -K, --use-genome-key  Force use of genome key file (auto-detect by default)\n");
    printf("  -N, --no-genome-key   Disable genome key file detection and translation\n");
    printf("  -C, --constraints FILE Parameter constraint file for bounds on error rates\n");
    printf("  -W, --weights FILE    Dirichlet prior weights file for mixture proportions\n");
    printf("  -j, --joint-damage    Estimate single shared damage rate for all entities (CEDfull only)\n");
    printf("  -h, --help            Show this help message\n");
    printf("  -V, --version         Show version information\n\n");
    
    printf("Methods:\n");
    printf("  A   - Single-rate model (no binomial coefficient)\n");
    printf("  AE  - Single-rate model with estimated rate\n");
    printf("  B   - Single-rate model (error rate / 3)\n");
    printf("  BE  - Single-rate model (error rate / 3, estimated)\n");
    printf("  C   - Damage model (fixed damage and background rates)\n");
    printf("  CE  - Damage model (fixed damage rate, estimated background rate)\n");
    printf("  CED - Damage model (estimated damage and background rates)\n");
    printf("  BEfull  - Per-genome error rates (standard model, K error rates)\n");
    printf("  CEfull  - Per-genome error rates (damage model, K error rates, fixed damage)\n");
    printf("  CEDfull - Per-genome error rates + per-genome or shared damage rates (use -j for shared)\n\n");
    
    printf("Note: SQUAREM acceleration is enabled by default. Use --no-squarem to disable.\n");
    printf("      Damage models (C, CE, CED) require input data with damage information.\n\n");
    
    printf("Examples:\n");
    printf("  %s -i data.txt -M A                     # Method A with SQUAREM (default)\n", program_name);
    printf("  %s -i data.txt -M AE -v                 # Method AE with verbose output\n", program_name);
    printf("  %s -i data.txt -M B -e 0.01 -o results  # Method B with custom error rate\n", program_name);
    printf("  %s -i data.txt -M C -r 0.05 -e 0.005    # Damage model C with custom rates\n", program_name);
    printf("  %s -i data.txt -M CED -v                # Damage model with estimated rates\n", program_name);
}

// Print version information
void print_version(void) {
    printf("CEM version 1.0.0\n");
    printf("C EM Algorithm for Mixture Proportion Estimation\n");
    printf("Based on unified_EM_proportions.R implementation\n\n");
    
    printf("Supported methods:\n");
    printf("  Method A (0A):  Single-rate model (no binomial coefficient)\n");
    printf("  Method 0AE:      Single-rate model (estimated rate)\n");
    printf("  Method B (0B):  Single-rate model (error rate / 3)\n"); 
    printf("  Method BE (0BE): Single-rate model (error rate / 3, estimated)\n\n");
    
    printf("SQUAREM acceleration: Enabled by default\n");

#ifdef HAVE_OPENMP
    printf("OpenMP parallelization: Enabled\n");
#else
    printf("OpenMP parallelization: Disabled\n");
#endif
}

// Parse method string
em_method_t parse_method(const char *method_str) {
    if (!method_str) return METHOD_A;
    
    // Standard methods
    if (strcmp(method_str, "A") == 0) return METHOD_A;
    if (strcmp(method_str, "AE") == 0) return METHOD_AE;
    if (strcmp(method_str, "B") == 0) return METHOD_B;
    if (strcmp(method_str, "BE") == 0) return METHOD_BE;
    
    // Damage models
    if (strcmp(method_str, "C") == 0) return METHOD_C;
    if (strcmp(method_str, "CE") == 0) return METHOD_CE;
    if (strcmp(method_str, "CED") == 0) return METHOD_CED;
    
    // Full models with per-genome rates
    if (strcmp(method_str, "BEfull") == 0) return METHOD_BEFULL;
    if (strcmp(method_str, "CEfull") == 0) return METHOD_CEFULL;
    if (strcmp(method_str, "CEDfull") == 0) return METHOD_CEDFULL;
    // Method Full removed - non-functional extended model
    
    // Also accept original R names for compatibility
    if (strcmp(method_str, "0A") == 0) return METHOD_A;
    if (strcmp(method_str, "0AE") == 0) return METHOD_AE;
    if (strcmp(method_str, "0B") == 0) return METHOD_B;
    if (strcmp(method_str, "0BE") == 0) return METHOD_BE;
    
    return METHOD_INVALID;  // Invalid method
}

// Get method name string
const char* get_method_name(em_method_t method) {
    switch (method) {
        case METHOD_A:       return "A (Single-rate model, no binomial coefficient)";
        case METHOD_AE:      return "AE (Single-rate model, estimated rate)";
        case METHOD_B:       return "B (Single-rate model, error rate / 3)";
        case METHOD_BE:      return "BE (Single-rate model, error rate / 3, estimated)";
        case METHOD_C:       return "C (Damage model, fixed rates)";
        case METHOD_CE:      return "CE (Damage model, estimated background rate)";
        case METHOD_CED:     return "CED (Damage model, estimated rates)";
        case METHOD_BEFULL:  return "BEfull (Per-genome error rates, standard model)";
        case METHOD_CEFULL:  return "CEfull (Per-genome error rates, damage model)";
        case METHOD_CEDFULL: return "CEDfull (Per-genome error and damage rates)";
        // METHOD_FULL removed - non-functional extended model
        case METHOD_INVALID: return "Invalid";
        default:             return "Unknown";
    }
}

// Generate output prefix from input filename
char* generate_output_prefix(const char *input_file) {
    if (!input_file) return NULL;
    
    // Find the last occurrence of '.' to remove extension
    const char *dot = strrchr(input_file, '.');
    const char *slash = strrchr(input_file, '/');
    
    // Use basename (part after last slash)
    const char *basename = slash ? slash + 1 : input_file;
    
    size_t len;
    if (dot && dot > basename) {
        len = dot - basename;
    } else {
        len = strlen(basename);
    }
    
    char *prefix = malloc(len + 1);
    if (!prefix) return NULL;
    
    strncpy(prefix, basename, len);
    prefix[len] = '\0';
    
    return prefix;
}

// Main function
int main(int argc, char *argv[]) {
    
    // Create configuration with defaults
    em_config_t *config = em_config_create();
    if (!config) {
        print_error("Failed to allocate configuration\n");
        return EXIT_FAILURE;
    }
    
    // Parse command line arguments
    int opt;
    int option_index = 0;
    
    while ((opt = getopt_long(argc, argv, "i:M:o:e:r:f:t:m:sSpdvcg:k:KNC:W:D:T:jhV", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'i':
                config->input_file = string_duplicate(optarg);
                break;
            case 'M':
                config->method = parse_method(optarg);
                if (config->method == METHOD_INVALID) {
                    print_error("Invalid method: %s\n", optarg);
                    print_error("Valid methods are: A, AE, B, BE, C, CE, CED, BEfull, CEfull, CEDfull\n");
                    print_error("Use -h for help.\n");
                    em_config_destroy(config);
                    return EXIT_FAILURE;
                }
                break;
            case 'o':
                config->output_prefix = string_duplicate(optarg);
                break;
            case 'e':
                config->error_rate = atof(optarg);
                if (config->error_rate <= 0.0 || config->error_rate >= 1.0) {
                    print_error("Error rate must be between 0 and 1\n");
                    em_config_destroy(config);
                    return EXIT_FAILURE;
                }
                break;
            case 'r':
                config->damage_rate = atof(optarg);
                // Use minimum of 1e-6 to avoid numerical issues with log(0)
                if (config->damage_rate < 1e-6) config->damage_rate = 1e-6;
                if (config->damage_rate >= 1.0) {
                    print_error("Damage rate must be less than 1.0\n");
                    em_config_destroy(config);
                    return EXIT_FAILURE;
                }
                break;
            case 'f':
                config->min_proportion = atof(optarg);
                if (config->min_proportion < 0.0 || config->min_proportion >= 1.0) {
                    print_error("Filter threshold must be between 0 and 1\n");
                    em_config_destroy(config);
                    return EXIT_FAILURE;
                }
                break;
            case 't':
                config->tolerance = atof(optarg);
                if (config->tolerance <= 0.0) {
                    print_error("Tolerance must be positive\n");
                    em_config_destroy(config);
                    return EXIT_FAILURE;
                }
                break;
            case 'm':
                config->max_iterations = atoi(optarg);
                if (config->max_iterations <= 0) {
                    print_error("Maximum iterations must be positive\n");
                    em_config_destroy(config);
                    return EXIT_FAILURE;
                }
                break;
            case 's':
                config->use_squarem = 1;
                break;
            case 'S':
                config->use_squarem = 0;
                break;
            case 'p':
                config->force_sparse = 1;
                config->force_dense = 0;
                break;
            case 'd':
                config->force_dense = 1;
                config->force_sparse = 0;
                break;
            case 'v':
                config->verbose = 1;
                break;
            case 'c':
                config->convert_names = 1;
                break;
            case 'g':
                config->mapping_file = string_duplicate(optarg);
                break;
            case 'k':
                config->genome_key_file = string_duplicate(optarg);
                config->use_genome_key = 2;  // Force use
                break;
            case 'K':
                config->use_genome_key = 2;  // Force use
                break;
            case 'N':
                config->use_genome_key = 0;  // Disable
                break;
            case 'C':  // -C, --constraints, --pt
                config->constraints_file = string_duplicate(optarg);
                break;
            case 'W':  // -W, --weights
                config->weights_file = string_duplicate(optarg);
                break;
            case 'j':
                config->joint_damage_rate = 1;
                break;
            case 'D':  // -D, --debug-likelihood
                debug_likelihood_mode = 1;
                debug_parameter_file = strdup(optarg);
                break;
            case 'T':  // -T, --test-likelihood-ratio
                likelihood_ratio_test_mode = 1;
                test_parameter_file = strdup(optarg);
                break;
            case 'h':
                print_usage(argv[0]);
                em_config_destroy(config);
                return EXIT_SUCCESS;
            case 'V':
                print_version();
                em_config_destroy(config);
                return EXIT_SUCCESS;
            case '?':
                print_error("Unknown option. Use -h for help.\n");
                em_config_destroy(config);
                return EXIT_FAILURE;
            default:
                print_error("Unexpected option: %c\n", opt);
                em_config_destroy(config);
                return EXIT_FAILURE;
        }
    }
    
    // Validate required arguments
    if (!config->input_file) {
        print_error("Input file is required. Use -i option.\n");
        print_error("Use -h for help.\n");
        em_config_destroy(config);
        return EXIT_FAILURE;
    }
    
    // Parse parameter constraints if specified
    if (config->constraints_file) {
        if (parse_constraints_file(config->constraints_file, &config->constraints, &config->n_constraints) != 0) {
            print_error("Failed to parse constraints file: %s\n", config->constraints_file);
            em_config_destroy(config);
            return EXIT_FAILURE;
        }
    }

    // Parse Dirichlet prior weights if specified
    if (config->weights_file) {
        if (parse_weights_file(config->weights_file, &config->weights, &config->n_weights) != 0) {
            print_error("Failed to parse weights file: %s\n", config->weights_file);
            em_config_destroy(config);
            return EXIT_FAILURE;
        }
    }

    // Validate likelihood ratio test mode - only works with CEDfull
    if (likelihood_ratio_test_mode) {
        if (config->method != METHOD_CEDFULL) {
            print_error("Likelihood ratio test (-T flag) is only supported with CEDfull method.\n");
            print_error("Please use -M CEDfull with -T option.\n");
            em_config_destroy(config);
            return EXIT_FAILURE;
        }
    }
    
    // Validate name conversion options
    if (config->convert_names && !config->mapping_file) {
        print_error("Name conversion enabled but no mapping file specified.\n");
        print_error("Use --mapping-file option to specify the genome mapping file.\n");
        em_config_destroy(config);
        return EXIT_FAILURE;
    }
    
    // Generate output prefix if not provided
    if (!config->output_prefix) {
        config->output_prefix = generate_output_prefix(config->input_file);
        if (!config->output_prefix) {
            print_error("Failed to generate output prefix\n");
            em_config_destroy(config);
            return EXIT_FAILURE;
        }
    }
    
    // Print basic configuration (always shown)
    printf("=== Configuration ===\n");
    fflush(stdout);
    printf("Input file: %s\n", config->input_file);
    fflush(stdout);
    printf("Method: %s\n", get_method_name(config->method));
    fflush(stdout);
    printf("Acceleration: %s\n", config->use_squarem ? "SQUAREM" : "Standard EM");
    fflush(stdout);
    if (config->method == METHOD_CEDFULL) {
        printf("Damage estimation: %s\n", config->joint_damage_rate ? "Joint" : "Separate");
        fflush(stdout);
    }

    // Print detailed configuration if verbose
    print_verbose(config->verbose, "Output prefix: %s\n", config->output_prefix);
    print_verbose(config->verbose, "Error rate: %.6f\n", config->error_rate);
    if (config->method == METHOD_C || config->method == METHOD_CE || config->method == METHOD_CED || 
        config->method == METHOD_CEFULL || config->method == METHOD_CEDFULL) {
        print_verbose(config->verbose, "Damage rate: %.6f\n", config->damage_rate);
    }
    print_verbose(config->verbose, "Tolerance: %.2e\n", config->tolerance);
    print_verbose(config->verbose, "Max iterations: %d\n", config->max_iterations);
    print_verbose(config->verbose, "Filter threshold: %.3f\n", config->min_proportion);
    printf("\n");
    fflush(stdout);
    
    // Detect data format
    printf("Detecting data format...\n");
    fflush(stdout);
    data_format_t detected_format = FORMAT_AUTO;
    input_format_t detailed_format = detect_input_format_from_file(config->input_file);
    
    // Check for unsupported format (dense with taxonomic groups)
    if (detailed_format == FORMAT_UNSUPPORTED) {
        print_error("Unsupported format detected\n");
        em_config_destroy(config);
        return EXIT_FAILURE;
    }
    
    if (config->force_sparse) {
        detected_format = FORMAT_SPARSE;
    } else if (config->force_dense) {
        detected_format = FORMAT_DENSE;
    } else {
        // Map detailed formats to sparse/dense
        if (detailed_format == FORMAT_SPARSE_STANDARD || detailed_format == FORMAT_SPARSE_DAMAGE) {
            detected_format = FORMAT_SPARSE;
        } else {
            detected_format = FORMAT_DENSE;
        }
    }
    
    config->data_format = detected_format;
    print_verbose(config->verbose, "Data format: %s\n", 
                  detected_format == FORMAT_SPARSE ? "sparse" : "dense");
    printf("Data format: %s\n", detected_format == FORMAT_SPARSE ? "sparse" : "dense");
    fflush(stdout);
    
    // Validate format-model compatibility
    int is_damage_format = (detailed_format == FORMAT_DENSE_DAMAGE || 
                            detailed_format == FORMAT_SPARSE_DAMAGE);
    int is_damage_model = (config->method == METHOD_C || 
                          config->method == METHOD_CE || 
                          config->method == METHOD_CED ||
                          config->method == METHOD_CEFULL || 
                          config->method == METHOD_CEDFULL);
    
    if (is_damage_format && !is_damage_model) {
        print_error("ERROR: Format-model mismatch detected!\n");
        print_error("Damage format (dense_dmg or sparse_dmg) detected but using standard model.\n");
        print_error("Standard models (A, AE, B, BE, BEfull) require standard format files.\n");
        print_error("Please use a standard format file or switch to a damage model (C, CE, CED, CEfull, CEDfull).\n");
        em_config_destroy(config);
        return EXIT_FAILURE;
    }
    
    if (!is_damage_format && is_damage_model) {
        print_error("ERROR: Format-model mismatch detected!\n");
        print_error("Standard format (dense_std or sparse_std) detected but using damage model.\n");
        print_error("Damage models (C, CE, CED, CEfull, CEDfull) require damage format files.\n");
        print_error("Please use a damage format file or switch to a standard model (A, AE, B, BE, BEfull).\n");
        em_config_destroy(config);
        return EXIT_FAILURE;
    }
    
    // Read input data (sparse or dense)
    printf("Reading input data...\n");
    fflush(stdout);
    double start_time = get_time();
    em_data_t *data = NULL;
    sparse_alignment_t *sparse_data = NULL;
    sparse_em_data_t *sparse_em_data = NULL;
    char **genome_names = NULL;
    
    if (detected_format == FORMAT_SPARSE) {
        // Parse sparse format
        input_format_t detailed_format = detect_input_format_from_file(config->input_file);
        if (detailed_format == FORMAT_SPARSE_DAMAGE && !getenv("FORCE_DENSE")) {
            sparse_data = parse_sparse_damage_file(config->input_file, &genome_names);
        } else if (detailed_format == FORMAT_SPARSE_STANDARD) {
            sparse_data = parse_sparse_standard_file(config->input_file, &genome_names);
        } else {
            sparse_data = parse_sparse_file(config->input_file, &genome_names);
        }
        if (!sparse_data) {
            print_error("Failed to read sparse input file\n");
            em_config_destroy(config);
            return EXIT_FAILURE;
        }
        
        // Validate taxonomic data compatibility
        if (sparse_data->n_taxonomic_groups > 0 && config->method != METHOD_CEDFULL) {
            print_error("ERROR: Taxonomic groups detected but model %s does not support taxonomic analysis.\n", 
                       get_method_name(config->method));
            print_error("Only CEDfull supports taxonomic group estimation.\n");
            print_error("Please use CEDfull model (-M CEDfull) or remove taxonomic data from input file.\n");
            sparse_alignment_free(sparse_data);
            em_config_destroy(config);
            return EXIT_FAILURE;
        }
        
        // Validate constraint file usage
        if (config->constraints_file && config->n_constraints > 0) {
            if (config->method != METHOD_CEDFULL) {
                printf("Warning: Constraint file loaded but model %s does not support parameter constraints.\n", 
                       get_method_name(config->method));
                printf("Parameter constraints only apply to taxonomic group rates in CEDfull model.\n");
            } else if (sparse_data->n_taxonomic_groups == 0) {
                printf("Warning: Constraint file loaded but no taxonomic groups detected in input data.\n");
                printf("Parameter constraints only apply to taxonomic group error rates.\n");
            }
        }
        
        // Set debug genome names for dynamic lookup if in debug mode
        if (debug_likelihood_mode && genome_names) {
            debug_genome_names = genome_names;
            debug_n_genomes = sparse_data->n_genomes;
        }
        
        // Set test genome names for dynamic lookup if in likelihood ratio test mode
        if (likelihood_ratio_test_mode && genome_names) {
            debug_genome_names = genome_names;  // Reuse the same mechanism
            debug_n_genomes = sparse_data->n_genomes;
        }
        
        // Create sparse EM data structure
        sparse_em_data = sparse_em_data_alloc(sparse_data);
        if (!sparse_em_data) {
            print_error("Failed to create sparse EM data structure\n");
            sparse_alignment_free(sparse_data);
            em_config_destroy(config);
            return EXIT_FAILURE;
        }
        
        
        // Allocate taxonomic sparse arrays after sparse_data is available
        if (sparse_data->has_taxonomic_groups && sparse_data->n_taxonomic_groups > 0) {
            int n_taxonomic_entries = 0;
            // Count total taxonomic entries across all reads
            for (int i = 0; i < sparse_data->n_reads; i++) {
                if (sparse_data->read_taxonomic_data && sparse_data->read_taxonomic_data[i].n_taxa > 0) {
                    n_taxonomic_entries += sparse_data->read_taxonomic_data[i].n_taxa;
                }
            }
            
            sparse_em_data->n_taxonomic_entries = n_taxonomic_entries;
            printf("Allocating sparse taxonomic EM storage: %d entries across %d reads for %d groups...\n", 
                   n_taxonomic_entries, sparse_data->n_reads, sparse_data->n_taxonomic_groups);
            
            // Allocate taxonomic CSR arrays (only for taxonomic functionality)
            sparse_em_data->taxonomic_row_ptr = calloc(sparse_data->n_reads + 1, sizeof(int));
            sparse_em_data->taxonomic_col_indices = malloc(n_taxonomic_entries * sizeof(int));
            
            // Allocate Q and W arrays for taxonomic groups
            sparse_em_data->Q_taxonomic_sparse = malloc(n_taxonomic_entries * sizeof(double));
            sparse_em_data->W_taxonomic_sparse = malloc(n_taxonomic_entries * sizeof(double));
            
            if (!sparse_em_data->Q_taxonomic_sparse || !sparse_em_data->W_taxonomic_sparse ||
                !sparse_em_data->taxonomic_row_ptr || !sparse_em_data->taxonomic_col_indices) {
                print_error("Failed to allocate sparse taxonomic EM storage\n");
                sparse_em_data_free(sparse_em_data);
                sparse_alignment_free(sparse_data);
                em_config_destroy(config);
                return EXIT_FAILURE;
            }
            
            // Populate taxonomic CSR arrays from per-read taxonomic data
            int current_tax_idx = 0;
            sparse_em_data->taxonomic_row_ptr[0] = 0;
            
            for (int i = 0; i < sparse_data->n_reads; i++) {
                if (sparse_data->read_taxonomic_data && sparse_data->read_taxonomic_data[i].n_taxa > 0) {
                    read_taxonomic_data_t *tax_data = &sparse_data->read_taxonomic_data[i];
                    for (int t = 0; t < tax_data->n_taxa; t++) {
                        if (current_tax_idx < n_taxonomic_entries) {
                            sparse_em_data->taxonomic_col_indices[current_tax_idx] = tax_data->taxon_indices[t];
                            current_tax_idx++;
                        }
                    }
                }
                sparse_em_data->taxonomic_row_ptr[i + 1] = current_tax_idx;
            }
            
            printf("Sparse taxonomic EM storage allocated successfully.\n");
        }
        
        sparse_em_data->error_rate = config->error_rate;
        sparse_em_data->damage_rate = config->damage_rate;
        sparse_em_data->joint_damage_rate = config->joint_damage_rate;

        // Copy constraints for M-step enforcement
        sparse_em_data->constraints = config->constraints;
        sparse_em_data->n_constraints = config->n_constraints;

        // Copy Dirichlet prior weights for M-step
        sparse_em_data->weights = config->weights;
        sparse_em_data->n_weights = config->n_weights;

        double end_time = get_time();
        printf("Data loaded: %d reads, %d genomes (%d alignments, %.1f%% sparse, %.2f seconds)\n", 
               sparse_data->n_reads, sparse_data->n_genomes, sparse_data->nnz,
               (1.0 - (double)sparse_data->nnz / ((double)sparse_data->n_reads * sparse_data->n_genomes)) * 100.0,
               end_time - start_time);
        fflush(stdout);  // Force flush
        
    } else {
        // Parse dense format
        if (read_mismatch_matrix(config->input_file, &data, &genome_names) != 0) {
            print_error("Failed to read dense input file\n");
            em_config_destroy(config);
            return EXIT_FAILURE;
        }
        
        data->error_rate = config->error_rate;
        data->damage_rate = config->damage_rate;
        
        double end_time = get_time();
        printf("Data loaded: %d reads, %d genomes (%.2f seconds)\n", 
               data->n_reads, data->n_genomes, end_time - start_time);
    }
    
    
    // Run EM algorithm based on method and acceleration choice
    printf("Starting EM optimization...\n");
    fflush(stdout);  // Force flush to see the output immediately
    double optimization_start_time = get_time();
    em_results_t *results = NULL;
    int iterations = 0;
    
    if (detected_format == FORMAT_SPARSE) {
        
        // Run sparse EM algorithm
        iterations = sparse_em_algorithm(sparse_em_data, config->method, 
                                       config->tolerance, config->max_iterations, 
                                       config->use_squarem, config->verbose,
                                       config->min_proportion, config->input_file);
        
        // Create results structure from sparse EM data
        results = em_results_create(sparse_em_data->n_genomes);
        if (results) {
            results->converged = (iterations < config->max_iterations);
            results->iterations = iterations;
            results->final_log_likelihood = sparse_em_data->log_likelihood;
            results->final_error_rate = sparse_em_data->error_rate;
            results->final_damage_rate = sparse_em_data->damage_rate;
            
            // BUGFIX: Ensure joint normalization before copying results
            if (sparse_em_data->taxonomic_proportions && sparse_data->n_taxonomic_groups > 0) {
                // Re-normalize genome + taxonomic proportions jointly to sum to 1.0
                double total = 0.0;
                
                // Sum all proportions
                for (int j = 0; j < sparse_em_data->n_genomes; j++) {
                    total += sparse_em_data->proportions[j];
                }
                for (int t = 0; t < sparse_data->n_taxonomic_groups; t++) {
                    total += sparse_em_data->taxonomic_proportions[t];
                }
                
                
                // Normalize both arrays by the same total
                if (total > 1e-12 && fabs(total - 1.0) > 1e-9) {
                    for (int j = 0; j < sparse_em_data->n_genomes; j++) {
                        sparse_em_data->proportions[j] /= total;
                    }
                    for (int t = 0; t < sparse_data->n_taxonomic_groups; t++) {
                        sparse_em_data->taxonomic_proportions[t] /= total;
                    }
                }
            }
            
            // Copy final proportions
            for (int i = 0; i < sparse_em_data->n_genomes; i++) {
                results->final_proportions[i] = sparse_em_data->proportions[i];
            }
            
            // Copy taxonomic group results if present
            const sparse_alignment_t *sparse_data = sparse_em_data->alignment_data;
            if (sparse_em_data->taxonomic_proportions && sparse_data->n_taxonomic_groups > 0) {
                results->n_taxonomic_groups = sparse_data->n_taxonomic_groups;
                results->final_taxonomic_proportions = malloc(sparse_data->n_taxonomic_groups * sizeof(double));
                results->taxonomic_groups = malloc(sparse_data->n_taxonomic_groups * sizeof(taxonomic_group_t));
                
                if (results->final_taxonomic_proportions && results->taxonomic_groups) {
                    for (int i = 0; i < sparse_data->n_taxonomic_groups; i++) {
                        results->final_taxonomic_proportions[i] = sparse_em_data->taxonomic_proportions[i];
                        
                        // Copy taxonomic group info
                        results->taxonomic_groups[i].taxid = sparse_data->taxonomic_groups[i].taxid;
                        results->taxonomic_groups[i].rank = strdup(sparse_data->taxonomic_groups[i].rank);
                        results->taxonomic_groups[i].name = strdup(sparse_data->taxonomic_groups[i].name);
                    }
                    
                    // Copy taxonomic rates for CEDfull (from unified arrays)
                    if (config->method == METHOD_CEDFULL && sparse_em_data->error_rates && sparse_em_data->damage_rates) {
                        // Allocate taxonomic rate arrays
                        results->final_taxonomic_error_rates = malloc(sparse_data->n_taxonomic_groups * sizeof(double));
                        results->final_taxonomic_damage_rates = malloc(sparse_data->n_taxonomic_groups * sizeof(double));
                        
                        if (results->final_taxonomic_error_rates && results->final_taxonomic_damage_rates) {
                            for (int i = 0; i < sparse_data->n_taxonomic_groups; i++) {
                                // Copy from unified arrays (taxonomic groups start at index n_genomes)
                                int unified_idx = sparse_em_data->n_genomes + i;
                                results->final_taxonomic_error_rates[i] = sparse_em_data->error_rates[unified_idx];
                                results->final_taxonomic_damage_rates[i] = sparse_em_data->damage_rates[unified_idx];
                            }
                        }
                    }
                }
            }
            
            // Copy per-genome error rates for BEfull, CEfull, CEDfull, Model Full
            if ((config->method == METHOD_BEFULL || config->method == METHOD_CEFULL || 
                 config->method == METHOD_CEDFULL) && sparse_em_data->error_rates) {
                if (!results->final_error_rates) {
                    results->final_error_rates = malloc(sparse_em_data->n_genomes * sizeof(double));
                }
                if (results->final_error_rates) {
                    for (int i = 0; i < sparse_em_data->n_genomes; i++) {
                        results->final_error_rates[i] = sparse_em_data->error_rates[i];
                    }
                }
            }
            
            // Copy per-genome damage rates for CEDfull and Model Full
            if ((config->method == METHOD_CEDFULL) && sparse_em_data->damage_rates) {
                if (!results->final_damage_rates) {
                    results->final_damage_rates = malloc(sparse_em_data->n_genomes * sizeof(double));
                }
                if (results->final_damage_rates) {
                    for (int i = 0; i < sparse_em_data->n_genomes; i++) {
                        results->final_damage_rates[i] = sparse_em_data->damage_rates[i];
                    }
                }
            }
        }
        
    } else {
        // Run dense EM algorithm
        // All models now support SQUAREM acceleration
        if (config->use_squarem) {
            results = optimize_with_squarem(data, config);
        } else {
            switch (config->method) {
                case METHOD_A:
                    results = method_A_fixed_rate(data, config);
                    break;
                case METHOD_AE:
                    results = method_AE_estimated_rate(data, config);
                    break;
                case METHOD_B:
                    results = method_B_fixed_rate(data, config);
                    break;
                case METHOD_BE:
                    results = method_BE_estimated_rate(data, config);
                    break;
                case METHOD_C:
                    results = method_C_fixed_rates(data, config);
                    break;
                case METHOD_CE:
                    results = method_CE_estimated_background(data, config);
                    break;
                case METHOD_CED:
                    results = method_CED_estimated_rates(data, config);
                    break;
                case METHOD_BEFULL:
                    results = method_BEfull_per_genome_rates(data, config);
                    break;
                case METHOD_CEFULL:
                    results = method_CEfull_per_genome_rates(data, config);
                    break;
                case METHOD_CEDFULL:
                    results = method_CEDfull_per_genome_rates(data, config);
                    break;
                default:
                    print_error("Unknown method\n");
                    em_data_destroy(data);
                    string_array_destroy(genome_names, data ? data->n_genomes : 0);
                    em_config_destroy(config);
                    return EXIT_FAILURE;
            }
        }
    }
    
    double optimization_end_time = get_time();
    
    if (!results) {
        print_error("EM algorithm failed\n");
        if (detected_format == FORMAT_SPARSE) {
            sparse_em_data_free(sparse_em_data);
            sparse_alignment_free(sparse_data);
            string_array_destroy(genome_names, sparse_data ? sparse_data->n_genomes : 0);
        } else {
            em_data_destroy(data);
            string_array_destroy(genome_names, data ? data->n_genomes : 0);
        }
        em_config_destroy(config);
        return EXIT_FAILURE;
    }
    
    // Apply filtering to final results
    int total_genomes = detected_format == FORMAT_SPARSE ? sparse_data->n_genomes : data->n_genomes;
    print_verbose(config->verbose, "Applying proportion filter (threshold: %.3f)...\n", config->min_proportion);
    filter_proportions(results->final_proportions, total_genomes, config->min_proportion);
    
    // Create array of genome indices sorted by proportion (descending)
    int n_genomes = detected_format == FORMAT_SPARSE ? sparse_data->n_genomes : data->n_genomes;
    int *sorted_indices = malloc(n_genomes * sizeof(int));
    if (!sorted_indices) {
        print_error("Failed to allocate memory for sorted indices\n");
        em_results_destroy(results);
        if (detected_format == FORMAT_SPARSE) {
            sparse_em_data_free(sparse_em_data);
            sparse_alignment_free(sparse_data);
        } else {
            em_data_destroy(data);
        }
        string_array_destroy(genome_names, n_genomes);
        em_config_destroy(config);
        return EXIT_FAILURE;
    }
    
    // Initialize indices
    for (int i = 0; i < n_genomes; i++) {
        sorted_indices[i] = i;
    }
    
    // Use efficient O(n log n) sorting for large datasets
    if (n_genomes > 1000) {
        printf("Sorting %d genomes by proportion (using quicksort)...\n", n_genomes);
        fflush(stdout);
    }
    
    // Set global proportions pointer for comparison function
    qsort_proportions = results->final_proportions;
    
    // Sort the indices using qsort (O(n log n) instead of O(n²))
    qsort(sorted_indices, n_genomes, sizeof(int), compare_genome_indices_desc);
    
    // Clear global pointer for safety
    qsort_proportions = NULL;
    
    // Load genome mapping if specified or auto-detect genome key file
    genome_mapping_t *genome_mappings = NULL;
    int n_mappings = 0;
    char *mapping_file_to_use = NULL;
    
    // Check if genome key functionality should be used
    int should_use_genome_key = 0;
    
    if (config->use_genome_key == 0) {
        // User explicitly disabled genome key
        should_use_genome_key = 0;
    } else if (config->use_genome_key == 2) {
        // User explicitly enabled genome key
        should_use_genome_key = 1;
    } else {
        // Auto-detect: only use if genome names are G1, G2, etc.
        should_use_genome_key = are_genome_names_short(genome_names, n_genomes);
        if (should_use_genome_key) {
            print_verbose(config->verbose, "Detected short genome names (G1, G2, ...), will look for genome key file\n");
        }
    }
    
    // Handle old-style name conversion (--convert-names with --mapping-file)
    if (config->convert_names && config->mapping_file) {
        mapping_file_to_use = config->mapping_file;
    } else if (should_use_genome_key) {
        // Try to find genome key file
        if (config->genome_key_file) {
            // User specified a genome key file
            mapping_file_to_use = config->genome_key_file;
        } else {
            // Auto-detect genome key file
            char *auto_detected_key_file = detect_genome_key_file(config->input_file);
            if (auto_detected_key_file) {
                print_verbose(config->verbose, "Auto-detected genome key file: %s\n", auto_detected_key_file);
                mapping_file_to_use = auto_detected_key_file;
            } else if (config->use_genome_key == 2) {
                // User forced genome key but we can't find one
                print_warning("Genome key file not found and --use-genome-key was specified\n");
            }
        }
    }
    
    if (mapping_file_to_use) {
        print_verbose(config->verbose, "Loading genome mapping from %s...\n", mapping_file_to_use);
        genome_mappings = load_genome_mapping(mapping_file_to_use, &n_mappings);
        if (!genome_mappings) {
            if (config->use_genome_key == 2 || config->genome_key_file) {
                // User explicitly requested genome key, warn if it fails
                print_warning("Could not load genome mapping file: %s\n", mapping_file_to_use);
            }
            // Don't fail - just continue without mapping
        } else {
            print_verbose(config->verbose, "Loaded %d genome mappings\n", n_mappings);
            // Store in config for use by write functions
            config->genome_mappings = genome_mappings;
            config->n_mappings = n_mappings;
            
            // Initialize hash table for fast genome name lookups
            print_verbose(config->verbose, "Initializing genome name hash table...\n");
            init_genome_hash_table(genome_mappings, n_mappings);
        }
        
        // Free auto-detected filename if we allocated it
        if (mapping_file_to_use != config->mapping_file && mapping_file_to_use != config->genome_key_file) {
            free(mapping_file_to_use);
        }
    }
    
    // Calculate timing and store in results BEFORE writing output
    double total_time = optimization_end_time - start_time;
    double optimization_time = optimization_end_time - optimization_start_time;
    results->computation_time = total_time;
    
    // Write output files
    print_verbose(config->verbose, "Writing output files...\n");
    if (write_results(config->output_prefix, results, genome_names, n_genomes, config) != 0) {
        print_error("Failed to write output files\n");
        if (genome_mappings) free_genome_mapping(genome_mappings, n_mappings);
        em_results_destroy(results);
        if (detected_format == FORMAT_SPARSE) {
            sparse_em_data_free(sparse_em_data);
            sparse_alignment_free(sparse_data);
        } else {
            em_data_destroy(data);
        }
        string_array_destroy(genome_names, n_genomes);
        em_config_destroy(config);
        return EXIT_FAILURE;
    }
    
    fflush(stdout);  // Ensure previous output is flushed
    printf("\n=== Final Summary ===\n");
    fflush(stdout);  // Ensure this line is printed
    printf("Method: %s\n", get_method_name(config->method));
    printf("Acceleration: %s\n", config->use_squarem ? "SQUAREM" : "Standard EM");
    printf("Converged: %s after %d iterations\n", 
           results->converged ? "Yes" : "No", results->iterations);
    printf("Final log-likelihood: %.6f\n", results->final_log_likelihood);
    
    // Print parameter estimates
    if (config->method == METHOD_AE || config->method == METHOD_BE) {
        printf("Estimated error rate: %.6f\n", results->final_error_rate);
    } else if (config->method == METHOD_CE || config->method == METHOD_CED) {
        printf("Estimated background error rate: %.6f\n", results->final_error_rate);
        if (config->method == METHOD_CED) {
            printf("Estimated damage rate: %.6f\n", results->final_damage_rate);
        }
    } else if (config->method == METHOD_BEFULL || config->method == METHOD_CEFULL || config->method == METHOD_CEDFULL) {
        printf("Per-genome error rates estimated for %d genomes\n", n_genomes);
        if (config->method == METHOD_CEDFULL) {
            printf("Per-genome damage rates estimated for %d genomes\n", n_genomes);
        }
    }
    
    // Show top reference genomes
    printf("Top reference genomes:\n");
    int top_count = n_genomes < 3 ? n_genomes : 3;
    for (int i = 0; i < top_count; i++) {
        int idx = sorted_indices[i];
        const char *name = genome_names[idx] ? genome_names[idx] : "Unknown";
        
        // Translate genome name if mapping is available
        if (genome_mappings && n_mappings > 0) {
            char *translated_name = translate_genome_name(name, genome_mappings, n_mappings);
            if (translated_name) {
                name = translated_name;
            }
        }
        
        printf("  %d. %s: %.6f (%.2f%%)\n", i+1, name, 
               results->final_proportions[idx], 
               results->final_proportions[idx] * 100.0);
    }
    
    printf("Optimization time: %.2f seconds\n", optimization_time);
    printf("Total analysis time: %.2f seconds\n", total_time);
    
    // Machine-readable summary for easy parsing
    printf("\n=== PARSE_SUMMARY_START ===\n");
    printf("METHOD=%s\n", get_method_name(config->method));
    printf("ACCELERATION=%s\n", config->use_squarem ? "SQUAREM" : "REGULAR_EM");
    printf("CONVERGED=%s\n", results->converged ? "YES" : "NO");
    printf("ITERATIONS=%d\n", results->iterations);
    printf("LOG_LIKELIHOOD=%.6f\n", results->final_log_likelihood);
    
    // Only print global error rate if not using per-genome error rates
    if (config->method != METHOD_BEFULL && config->method != METHOD_CEFULL && config->method != METHOD_CEDFULL) {
        printf("ERROR_RATE=%.6f\n", results->final_error_rate);
    }
    
    // Only print global damage rate if not using per-genome damage rates and if it's a damage model
    if (config->method == METHOD_C || config->method == METHOD_CE || config->method == METHOD_CED || config->method == METHOD_CEFULL) {
        if (config->method != METHOD_CEDFULL) {  // CEDfull has per-genome damage rates
            printf("DAMAGE_RATE=%.6f\n", results->final_damage_rate);
        }
    }
    
    printf("OPTIMIZATION_TIME=%.2f\n", optimization_time);
    printf("TOTAL_TIME=%.2f\n", total_time);
    printf("=== PARSE_SUMMARY_END ===\n");
    
    printf("\nAnalysis completed successfully\n");
    
    // Cleanup
    free(sorted_indices);
    if (genome_mappings) free_genome_mapping(genome_mappings, n_mappings);
    em_results_destroy(results);
    if (detected_format == FORMAT_SPARSE) {
        sparse_em_data_free(sparse_em_data);
        sparse_alignment_free(sparse_data);
        string_array_destroy(genome_names, sparse_data->n_genomes);
    } else {
        em_data_destroy(data);
        string_array_destroy(genome_names, data->n_genomes);
    }
    em_config_destroy(config);
    
    return EXIT_SUCCESS;
}