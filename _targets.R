##########################################################################
#                                                                        #
#                          Run allometry analyses                        #
#                                                                        #
##########################################################################

######## Options and packages

# Loading targets
library(targets)

# Loading functions
lapply(list.files("R", pattern = ".R$", full.names = TRUE), source)

# Installing if needed and loading packages
packages.in <- c("dplyr", "plyr", "stringr", "tidyverse",
                 "ggplot2", "cowplot", "ggrepel", "ggfortify",
                 "factoextra", "MuMIn", "corrplot")


for (pkg in packages.in) if(!(pkg %in% rownames(installed.packages()))) install.packages(pkg)

# Specifying target options
options(tidyverse.quiet = TRUE, clustermq.scheduler = "multiprocess")
tar_option_set(packages = packages.in,
               memory = "transient", garbage_collection = TRUE)


list(

  ### 1. NO COMPETITION ANALYSES
  # compute mean crown characteristics using selected allometry models and extracted parameters
  # perform PCA analysis: co-variations between crown characteristics
  # investigate effects of key functional traits and mean climate niche on crown tree allometry independently for each functional group
    # first for angiosperms,
    # then for gymnosperms

  tar_target(mean_characteristics_nocomp, compute_mean_crown_characteristics_nocomp()),
  tar_target(res_pca_example, pca_ecological_strategies(mean_characteristics_nocomp)), # do not return anything, PCA results are exported directly as figures in the manuscript folder

  tar_target(angiosperms_scaled, traits_climate_crown_model_selection_angiosperms(mean_characteristics_nocomp)),
  tar_target(gymnosperms_scaled, traits_climate_crown_model_selection_gymnosperms(mean_characteristics_nocomp)),
  tar_target(angio_nocomp_results, traits_climate_crown_model_running_angiosperms(angiosperms_scaled)),
  tar_target(gymno_nocomp_results, traits_climate_crown_model_running_gymnosperms(gymnosperms_scaled)),



  ### 2. COMPETITION ANALYSES
  # define species list based on data availability and competition models' robustness
  # compute crown characteristics response to competition and assess effects of species' shade tolerance on the latter
  # do the same on unscaled data for visual representation of the results
  # assess the relative importance of crown characteristics on light interception
    # for an isolated tree
    # for a codominant tree or a supressed tree

  tar_target(species_list_comp, species_selection_competition()),
  tar_target(response_to_competition, response_to_competition(species_list_comp)),
  tar_target(response_to_competition_unscaled, response_to_competition_unscaled(species_list_comp)),
  tar_target(isolated_tree_results, samsara_isolated_tree(species_list_comp)),
  tar_target(codominant_tree_results, samsara_codominated_tree(species_list_comp)),



  ### 3. SUPPORTING INFORMATION
  # determine whether weighted or non-weighted protocol effects are the best fit
  # determine whether BAT or BAL is the best proxy of the competitive abilities of tree species
  # compute mean parameters +/- SD for each crown characteristic with and without competition
  # perform pearson correlations between all considered crown characteristics, climate variables and functional traits
  # implement univariate model to investigate effects of key functional traits and mean climate niche on crown tree allometry
  # compute bivariate scatterplots

  tar_target(protocol_selection, protocol_parameter(mean_characteristics_nocomp)),
  tar_target(CI_selection, competition_index(mean_characteristics_nocomp)),
  tar_target(mean_parameters_nocomp, compute_mean_parameters_nocomp(mean_characteristics_nocomp)),
  tar_target(mean_parameters_comp, compute_mean_parameters_comp(species_list_comp)),
  tar_target(cor_matrix_values, pearson_correlation(mean_characteristics_nocomp)),
  tar_target(nothing_useful, univariate_models()), # only compute figures, does not return anything useful
  tar_target(supp_file, bivariate_scatterplots()))












