##########################################################################
#                                                                        #
#                       Run allometry relationships                      #
#                                                                        #
##########################################################################

######## Options and packages

# Loading targets
library(targets)

# Loading functions
lapply(list.files("R", pattern = ".R$", full.names = TRUE), source)

# Installing if needed and loading packages
packages.in <- c("stringr", "dplyr", "plyr", "sp", "rworldmap", "rgdal", "measurements",
                 "sf", "readxl", "stringi", "lubridate", "tidyr", "parzer", "TNRS",
                 "ggplot2", "clustermq", "nlme", "lme4", "betareg", "pals", "gdata", "Metrics",
                 "purrr", "readr", "magrittr", "rgbif", "data.table",
                 "CoordinateCleaner", "RCurl", "httr", "archive", "terra", "cowplot", "R.utils",
                 "ggspatial", "rnaturalearth", "rnaturalearthdata",
                 "ggstatsplot", "truncnorm")


for (pkg in packages.in) if(!(pkg %in% rownames(installed.packages()))) install.packages(pkg)


# Specifying target options
options(tidyverse.quiet = TRUE, clustermq.scheduler = "multiprocess")
tar_option_set(packages = packages.in,
               memory = "transient", garbage_collection = TRUE)


list(

  ### 1. Extracting data for all explored allometric relationships

  tar_target(global_species_list, get_species_list()),
  tar_target(data_allometry, get_data_allometry(global_species_list)),

  tar_target(height_data, get_data_height(data_allometry)),
  tar_target(height_species, get_species_height(height_data)),
  tar_target(height_species_comp, get_species_height_comp(height_data)),

  tar_target(diameter_data, get_data_diameter(data_allometry)),
  tar_target(diameter_species, get_species_diameter(diameter_data)),
  tar_target(diameter_species_comp, get_species_diameter_comp(diameter_data)),

  tar_target(ratio_data, get_data_ratio(data_allometry)),
  tar_target(ratio_species, get_species_ratio(ratio_data)),
  tar_target(ratio_species_comp, get_species_ratio_comp(ratio_data)),


  ### 2. Fitting all allometric relationships (use of set.seed)
  tar_target(height_asympt_nocomp, height_models_asympt(height_data, height_species), pattern = map(height_species)),
  tar_target(height_power_nocomp, height_models_power_nocomp(height_data, height_species), pattern = map(height_species)),
  tar_target(height_power_comp, height_models_power_comp(height_data, height_species_comp), pattern = map(height_species_comp)),

  tar_target(diameter_nocomp, diameter_models_nocomp(diameter_data, diameter_species), pattern = map(diameter_species)),
  tar_target(diameter_comp, diameter_models_comp(diameter_data, diameter_species_comp), pattern = map(diameter_species_comp)),

  tar_target(beta_resampling_nocomp, ratio_models_nocomp(ratio_data, ratio_species), pattern = map(ratio_species)),
  tar_target(beta_resampling_comp, ratio_models_comp(ratio_data, ratio_species_comp), pattern = map(ratio_species_comp)),


  ### 3. Extracting parameters for all considered models
  tar_target(nocomp_parameters, extracting_nocomp_parameters()),
  tar_target(comp_parameters, extracting_comp_parameters()))







