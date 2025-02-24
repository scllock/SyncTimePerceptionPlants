# Summary of scripts

Scripts in this folder are numbered by steps in the analysis.

### List of scripts

-   `0_load_packages_and_functions.R` loads the packages and custom functions used throughout analysis
-   `1_flowering_time_experiment_and_data_processing.R` loads in the two flowering time experiments, compares experiments, produces plots exploring data, saves necessary outputs
-   `2_TPM_tables_DGE_prep_for_Elasticnet.R` filters TPM tables, performs differential gene expression, saves output of relevant results
-   `3_Post_Elasticnet_gene_analysis.R` looks at elastic net results and identifies gene lists of important genes from different Elastic Net models
-   `4_filtering_and_smoothing_topcount_data.R` loads raw topcount ciradian data and filters and smooths by fitting B-splines to curves based on chosen parameters
-   `5_FPCA_fullcurves_functional_traits.R` performs FPCA on smooth curves and identifies different traits associated 
-   `6_FPCA_pre_post_functional_traits.R` performs FPCA on smooth curves from pre and post perceived missed light cue and identifies different traits associated
-   `7_FPCA_sd_curves_functional_traits.R` performs FPCA on smooth standard deviation curves for each genotype and identifies different traits associated
-   `8_time_warping.R` performs time warping on smooth curves, taking each individual curve and warping to the mean of each genotype  
-   `9_FPCA_time_warping_functional_traits.R` performs FPCA on smooth time warping curves for each genotype and identifies different traits associated
-   `10_QTL_analysis.R` estimates genetic map by sequencing. Then uses this to perform single QTL analysis on all functional and developmental traits. Filters results for QTLs associated with both developmental and functional. Produces a network showing these connections. Saves appropriate outputs. 
-   `11_Transcriptome_plots.R` produces plots associated with transcriptome data 
-   `12_FDA_plots.R` produces plots associated with functional data analysis section 
-   `13_QTL_plots.R` produces plots associated with QTL analysis  
-   `14_ecotypes.R` produces plots associated with mutational profiles of Arabidopsis ecotypes
-   `15_splicing.R` produces plots related to alternative splicing