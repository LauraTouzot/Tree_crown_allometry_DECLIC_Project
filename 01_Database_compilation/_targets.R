##########################################################################
#                                                                        #
#                Compile and clean the allometry database                #
#                                                                        #
##########################################################################

######## Options and packages

# Loading targets
library(targets)

# Loading functions
lapply(list.files("R", pattern = ".R$", full.names = TRUE), source)

# Installing if needed and loading packages
packages.in <- c("baad.data", "stringr", "dplyr", "plyr", "sp", "rworldmap", "rgdal", "measurements",
                 "sf", "readxl", "stringi", "lubridate", "tidyr", "parzer", "TNRS",
                 "ggplot2", "clustermq", "nlme", "lme4", "betareg", "pals", "gdata", "Metrics",
                 "purrr", "readr", "magrittr", "rgbif", "data.table",
                 "CoordinateCleaner", "RCurl", "httr", "archive", "terra", "cowplot", "R.utils",
                 "ggspatial", "rnaturalearth", "rnaturalearthdata", "taxize",
                 "ggstatsplot", "ade4", "factoextra", "FactoMineR",
                 "truncnorm")


for (pkg in packages.in) if(!(pkg %in% rownames(installed.packages()))) install.packages(pkg)


# Specifying target options
options(tidyverse.quiet = TRUE, clustermq.scheduler = "multiprocess")
tar_option_set(packages = packages.in,
               memory = "transient", garbage_collection = TRUE)


list(

  ### 1. Compiling the database

  # baad database
  tar_target(BAAD_crown_sp_file, "data/BAAD/res_crown_BAAD.csv", format = "file"),
  tar_target(BAAD_crown_sp, read.csv(BAAD_crown_sp_file)),
  tar_target(baad, baad_data()),
  tar_target(baad_f, prepare_dataset_baad(baad)),
  tar_target(baad_crown, data_formatted_baad(baad_f)),
  tar_target(baad_crown_complete, extract_continents_baad(baad_crown))

  # FHM database
  tar_target(FHM_data_tree, read_FHM_data()),
  tar_target(FHM_data, get_FHM_weight(FHM_data_tree)),

  # FIA database
  tar_target(FIA_sp_file, "data/FIA/Crown_FIA_Tree_Plot_Coord.csv", format = "file"),
  tar_target(FIA_data, read.csv(FIA_sp_file)),

  # Spanish NFI data
  tar_target(Spain_NFI, read_spanish_data()),

  # FUNDIV data
  tar_target(FUNDIV_data, read_FUNDIV_data()),

  # FUNDIV crown data
  tar_target(FUNDIV_crown_data, prepare_FUNDIV_explore_crown()),
  tar_target(FUNDIV_crown, get_FUNDIV_coordinates(FUNDIV_crown_data)),

  # MONTANE data (unpublished data from Belledonne)
  tar_target(MONTANE_data, read_MONTANE()),

  # data from Fuhr et al. (2017)
  tar_target(data_paper_crown, crown_data_paper()),

  # data from Evans et al. (2015)
  tar_target(evan_crown_data, data_evans_crown()),

  # data from Dettmann and MacFarlane (2018)
  tar_target(dettmann_data, data_dettmann()),

  # data from Heym et al. (2017)
  tar_target(heym_2017, read_heym_crown_data()),

  # GenTree database
  tar_target(gen_tree_data, GenTree_data()),

  # Legacy Tree database
  tar_target(legacytree, Legacy_Tree_crown_data()),

  # ICP database
  tar_target(ICP_data, read_ICP()),

  # data from Anderson-Teixeira et al. (2015)
  tar_target(anderson2015, Anderson_2015_crown_data()),

  # data from Dalponte and Coomes (2016)
  tar_target(dalponte2016, Dalponte_2016_crown_data()),

  # French NFI data
  tar_target(FrenchNFI_data, french_NFI_crown_data()),

  # Quebec NFI data
  tar_target(Quebec, quebec_NFI_crown_data()),

  # Canada NFI data
  tar_target(Canada_data, canada_NFI_crown_data()),

  # data from Usoltsev database
  tar_target(Usoltsev, Usoltsev_data()),

  # data from Sullivan et al. (2018)
  tar_target(Sullivan_data, read_Sullivan()),

  # Tallo database (2022)
  tar_target(tallo_database, read_Tallo()),

  # merge all data collected on tree crowns in a single file
  tar_target(all_crown_a, merge_crown_data(baad_crown_complete, FHM_data, FIA_data, Spain_NFI, FUNDIV_data, FUNDIV_crown,
                                 MONTANE_data, data_paper_crown, evan_crown_data, dettmann_data, heym_2017, gen_tree_data,
                                 legacytree, ICP_data, anderson2015, dalponte2016, FrenchNFI_data, Quebec, Canada_data,
                                 Usoltsev, Sullivan_data, tallo_database)),

  # cleaning the database
  tar_target(all_crown, remove_duplicated_ref(all_crown_a)),
  tar_target(all_crown_clean, last_cleaning(all_crown)),
  tar_target(all_crown_checked, checking_tree_data(all_crown_clean)),
  tar_target(all_crown_checked_bis, cleaning_after_checking(all_crown_checked)),



  ### 2. Checking for taxonomy and computing supplementary variables

  tar_target(allometry_database, check_for_taxonomy_allometry(all_crown_checked_bis)),
  tar_target(allometry_supp_variables, compute_supplementary_variables(allometry_database)),
  tar_target(location_variables, extract_supplementary_variables(allometry_supp_variables)),
  tar_target(location_variables_checked, checking_plot_data(location_variables)),
  tar_target(allometry_complete_database, complete_allometry(allometry_supp_variables, location_variables_checked)),


  ### 3. Summary of the allometry database for exploration and threshold definitions

  tar_target(summary_dataset, summarizing_allometry_dataset(allometry_supp_variables)),
  tar_target(summary_species, summarizing_allometry_species(allometry_supp_variables)),


  ### 4. Comparing data availability: allometry database vs. NFI

  tar_target(species_all_NFI, read.csv("data/species_NFI_all_climate_exo.csv")),
  tar_target(NFI_data, check_for_taxonomy_NFI(species_all_NFI)),
  tar_target(data_availability_comparison, binding_databases(NFI_data, summary_species))

  )
