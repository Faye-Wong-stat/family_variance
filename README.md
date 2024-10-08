
# Why usefulness is rarely useful


This document lists the scripts used in this project. The order to run scripts follows the order of scripts listed in this document, as later scripts use data generated from previous scripts. 

### Prepare data for phasing
"remove_end_progenies.R": 
    input data: 
        "pedigree_genotyped.txt"
        "marker.txt"
        "marker_info.txt"

"generate_vcffiles.R"

### Phase
"run_phase_parent.sh"; "run_phase_valid.sh"; "run_phase_all.sh"

### Calculate impuation and phasing error
"genotype_error.R"; "imputation_error.R"; "view_errors.R"; "get_errors.R"; "get_genotype_errors.R"; "get_imputed_errors.R": 
    archived

"get_errors_imputation.R" 

"get_markers_cM.R", "get_markers_dist.R":
    input data:
        "20A823_Genetic_Map.csv"

"phase_correct_funcs.R": 
    this contains the important functions used to calculate phasing errors

"get_phase_error.R" 

"view_phase_error.R": 
    archived

### Simulate offspring based on real genotypes
"create_marker_list.R"

"simulate_crosses.R"

### Simulate BV and phenotypes on real parents and simulated offspring
"simulate_phenotypes.R"

"simulate_phenotypes_crosses.R"

"simulate_phenotypes_crosses2.R": 
    archived

### Calculate the predictive ability of family mean, variance and usefulness
"view_usefulness.R"
    Figures 2-4: 6.5*4.33; 
    Figures S1-S2: 6.5*4.33

"examine_usefulness.R": 
    archived

### Test if family variance is greater in best parents
"select_best_parents.R"

"simulate_crosses_best_parents.R"

"view_crosses_best_parents.R": 
    archived

"view_usefulness_best_parents.R": 
    Figure 6: 6.5*4.33

"plot_BV_family_mean_sd.R": 
    Figure 5: 6.5*3.3
    Figure S5: 6.5*4.33

### Introduce phasing error
"create_marker_list2.R": 
    archived

"introduce_error.R"

"examine_error.R"

"view_examine_error.R":
    Figure 1: 6.5*4

### housekeeping: archived
"view_phenotypes.R": Figure S3-S4
"simulate_phenotypes_crosses3.R"
"view_correlation.R"

### Calculate the mean and variance by gametes
"simulate_gametes.R"

"simulate_phenotypes_gametes.R"

"view_correlation_gametes.R":
    Figure S6: 6.5*3.3

"extract_gamete_info.R"

"select_best_parents_gametes.R"

"simulate_crosses_best_parents_gametes.R"; "simulate_crosses_best_parents_gametes_use.R"; "view_usefulness_best_parents_gametes.R": 
    archived

"view_usefulness_best_parents_gametes.R":
    Figure 7: 6.5*4.33