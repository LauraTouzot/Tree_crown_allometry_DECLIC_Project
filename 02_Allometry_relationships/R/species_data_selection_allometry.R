get_species_list <- function() {
  
  ### Loading data
  allometry_complete_database <- readRDS(file = "data/allometry_complete_database.RDS")
  NFI_data = readRDS(file = "data/NFI_TNRS_check.rds")
   
  ## extracting species list from NFI data (191 species)
  sampling <- NFI_data %>% dplyr::filter(continent == "E_U" & nplot >= 100 & ntree >= 1000 | continent == "N_A" & nplot >= 150 & ntree >= 3000)
  
  ## extracting species list in the allometry database (180 species)
  data <- allometry_complete_database
  data <- data %>% dplyr::ungroup() # just to be sure :)
  species <- unique(data$checked_name)

  data_summary <- data %>% dplyr::group_by(checked_name) %>% 
                           dplyr::summarise(nplot_crown = length(unique(location_ID)),
                                            ntree_crown = length(location_ID)) %>% 
                           dplyr::ungroup()
  
  sampling <- left_join(sampling, data_summary, by = "checked_name")
  
  selected_sp <- sampling %>% filter(continent == "E_U" & nplot_crown >= 100 & ntree_crown >= 1000 | continent == "N_A" & nplot_crown >= 150 & ntree_crown >= 3000)
  
  species_list <- unique(selected_sp$checked_name)
  species_list <- sort(species_list) # do not forget to order species list so that the rest of the code makes sense
  species_list <- species_list[-1]
  
  rm(allometry_complete_database)
  rm(data)
  rm(NFI_data)
  rm(selected_sp)
  
  return(species_list)
  
}





get_data_allometry <- function(global_species_list) {
  
  ### Loading data
  allometry_complete_database <- readRDS(file = "data/allometry_complete_database.RDS")
  data_ok <- allometry_complete_database %>% filter(checked_name %in% global_species_list) 
  
  rm(allometry_complete_database)
  
  return (data_ok)
  
}





get_data_height <- function(data_allometry) {
  
  # selecting data within the global species list and filtering individual observations
  data_height_a <- data_allometry %>% dplyr::filter(!is.na(DBH_cm) & !is.na(HT_m) & HT_m > 1.3) %>%
                                      dplyr::select(checked_name, DBH_cm, HT_m, location_ID, data, ba_plot, ba_larger_trees) %>%
                                      dplyr::rename(sp_name = checked_name, x = DBH_cm, y = HT_m, location = location_ID, protocol = data,
                                                    ba_plot = ba_plot, ba_larger = ba_larger_trees) %>% 
                                      dplyr::mutate(sp_name = as.character(sp_name), x = as.numeric(x), y = as.numeric(y), 
                                                    location = as.factor(location), protocol = as.factor(protocol),
                                                    ba_plot = as.numeric(ba_plot), ba_larger = as.numeric(ba_larger), id = as.numeric(1:n())) 
  
  
  # removing all plots with less than 2 observations and protocols with less than 9 observations from the data from the data
  sel_location_height <- names(table(data_height_a$location))[table(data_height_a$location) > 2]
  sel_protocol_height <- names(table(data_height_a$protocol))[table(data_height_a$protocol) > 9]
  data_height_b <- data_height_a[data_height_a$location %in% sel_location_height & data_height_a$protocol %in% sel_protocol_height, ]
  
  
  # selecting species with more than 300 observations for both the asymptotic and power-law relationships
  species_height <- data_height_b %>% dplyr::group_by(sp_name) %>% 
                                      dplyr::summarise(nobs_HT = sum(!is.na(y))) %>% 
                                      dplyr::filter(nobs_HT >= 200) %>%
                                      dplyr::ungroup()
  
  
  # extracting height data 
  species_height_list <- unique(species_height$sp_name)
  data_height <- data_height_b[data_height_b$sp_name %in% species_height_list,]
  
  # removing unused files
  rm(data_height_a, data_height_b, sel_location_height, species_height, species_height_list)
  gc()
  
  # returning data and species list
  return (data_height)
  
}


get_species_height <- function(height_data) {
  
  species_height_list <- unique(height_data$sp_name)
  return(species_height_list)
  
}


get_species_height_comp <- function(height_data) {
  
  data <- height_data
  
  summary <- data %>% dplyr::filter(!is.na(ba_plot) & !is.na(ba_larger)) %>%
                      dplyr::group_by(sp_name) %>% 
                      dplyr::summarise(comp_count = n()) %>%
                      dplyr::filter(comp_count > 200)
  
  new_sp_list_height <- summary$sp_name
  
  return(new_sp_list_height)
  
  
}






get_data_diameter <- function(data_allometry) {
  
  # selecting data within the global species list and filtering individual observations
  data_diameter_a <- data_allometry %>% dplyr::filter(!is.na(DBH_cm) & !is.na(C_diam_m) & C_diam_m > 0) %>%
                                        dplyr::select(checked_name, DBH_cm, C_diam_m, location_ID, data, ba_plot, ba_larger_trees) %>%
                                        dplyr::rename(sp_name = checked_name, x = DBH_cm, y = C_diam_m, location = location_ID, protocol = data, 
                                                      ba_plot = ba_plot, ba_larger = ba_larger_trees) %>% 
                                        dplyr::mutate(sp_name = as.character(sp_name), x = as.numeric(x), y = as.numeric(y), 
                                                      location = as.factor(location), protocol = as.factor(protocol),
                                                      ba_plot = as.numeric(ba_plot), ba_larger = as.numeric(ba_larger), id = as.numeric(1:n())) 
  
  
  # removing all plots with less than 2 observations and protocols with less than 9 observations from the data
  sel_location_diameter <- names(table(data_diameter_a$location))[table(data_diameter_a$location) > 2]
  sel_protocol_diameter <- names(table(data_diameter_a$protocol))[table(data_diameter_a$protocol) > 9]
  data_diameter_b <- data_diameter_a[data_diameter_a$location %in% sel_location_diameter & data_diameter_a$protocol %in% sel_protocol_diameter, ]
  
  # selecting species with more than 300 observations
  species_diameter <- data_diameter_b %>% dplyr::group_by(sp_name) %>% 
                                          dplyr::summarise(nobs_diam = sum(!is.na(y))) %>% 
                                          dplyr::filter(nobs_diam >= 200) %>%
                                          dplyr::ungroup()
  
  # extracting diameter data and diameter species list
  species_diameter_list <- unique(species_diameter$sp_name)
  data_diameter <- data_diameter_b[data_diameter_b$sp_name %in% species_diameter_list,]
  
  # removing unused files
  rm(data_diameter_a, data_diameter_b, sel_location_diameter, species_diameter, species_diameter_list)
  gc()
  
  # returning data and species list
  return (data_diameter)
  
}


get_species_diameter <- function(diameter_data) {
  
  species_diameter_list <- unique(diameter_data$sp_name)
  return(species_diameter_list)
  
}


get_species_diameter_comp <- function(diameter_data) {
  
  data <- diameter_data
  
  summary <- data %>% dplyr::filter(!is.na(ba_plot) & !is.na(ba_larger)) %>%
                      dplyr::group_by(sp_name) %>% 
                      dplyr::summarise(comp_count = n()) %>%
                      dplyr::filter(comp_count > 200)
  
  new_sp_list_diameter <- summary$sp_name
  
  return(new_sp_list_diameter)
  
  
}






get_data_ratio <- function(data_allometry) {
  
  # selecting data within the global species list and filtering individual observations
  data_ratio_a <- data_allometry %>% dplyr::filter(!is.na(DBH_cm) & !is.na(CR) & CR > 0 & CR < 1) %>%
                                     dplyr::select(checked_name, DBH_cm, CR, location_ID, data, ba_plot, ba_larger_trees) %>%
                                     dplyr::rename(sp_name = checked_name, x = DBH_cm, y = CR, location = location_ID, protocol = data, 
                                                   ba_plot = ba_plot, ba_larger = ba_larger_trees) %>% 
                                     dplyr::mutate(sp_name = as.character(sp_name), x = as.numeric(x), y = as.numeric(y), 
                                                   location = as.factor(location), protocol = as.factor(protocol),
                                                   ba_plot = as.numeric(ba_plot), ba_larger = as.numeric(ba_larger), id = as.numeric(1:n())) 
  
  
  # removing all plots with less than 2 observations and protocols with less than 9 observations from the data
  sel_location_ratio <- names(table(data_ratio_a$location))[table(data_ratio_a$location) > 2]
  sel_protocol_ratio <- names(table(data_ratio_a$protocol))[table(data_ratio_a$protocol) > 9]
  data_ratio_b <- data_ratio_a[data_ratio_a$location %in% sel_location_ratio & data_ratio_a$protocol %in% sel_protocol_ratio, ]  
  
  # selecting species with more than 300 observations
  species_ratio <- data_ratio_b %>% dplyr::group_by(sp_name) %>% 
                                    dplyr::summarise(nobs_ratio = sum(!is.na(y))) %>% 
                                    dplyr::filter(nobs_ratio >= 200) %>%
                                    dplyr::ungroup()
                                  
  # extracting ratio data and ratio species list
  species_ratio_list <- unique(species_ratio$sp_name)
  data_ratio <- data_ratio_b[data_ratio_b$sp_name %in% species_ratio_list,]
  
  # removing unused files
  rm(data_ratio_a, data_ratio_b, sel_location_ratio, species_ratio, species_ratio_list)
  gc()
  
  # returning data and species list
  return (data_ratio)
  
}


get_species_ratio <- function(ratio_data) {
  
  species_ratio_list <- unique(ratio_data$sp_name)
  return(species_ratio_list)
  
}


get_species_ratio_comp <- function(ratio_data) {
  
  data <- ratio_data
  
  summary <- data %>% dplyr::filter(!is.na(ba_plot) & !is.na(ba_larger)) %>%
                      dplyr::group_by(sp_name) %>% 
                      dplyr::summarise(comp_count = n()) %>%
                      dplyr::filter(comp_count > 200)
  
  new_sp_list_ratio <- summary$sp_name
  
  return(new_sp_list_ratio)
  
  
}













  
