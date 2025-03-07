###### Running height models #######

height_models_asympt <- function(height_data, height_species) {

  ## loading data and species list
  data_ok <- height_data
  species_list <- height_species
  
  for (i in 1:length(species_list)) {

    height_species <- species_list[i]
    print(i)
  
    ## defining number of repetitions
    n_repetition = 200
  
    ## creating file to store model parameters (1 file per species)
    asymptot_resampling_nocomp <- as.data.frame(matrix(nrow = n_repetition, ncol = length(unique(data_ok$protocol)) + 7)) 
    asymptot_resampling_nocomp[,1] <- rep(height_species, n_repetition)
    names(asymptot_resampling_nocomp) <- c("species", paste0("protocol", unique(data_ok$protocol)), "b1", "b2", "b3", "AIC", "RMSE", "weighted")
  
    asymptot_resampling_nocomp_w <- as.data.frame(matrix(nrow = n_repetition, ncol = length(unique(data_ok$protocol)) + 7)) 
    asymptot_resampling_nocomp_w[,1] <- rep(height_species, n_repetition)
    names(asymptot_resampling_nocomp_w) <- c("species", paste0("protocol", unique(data_ok$protocol)),"b1", "b2", "b3", "AIC", "RMSE", "weighted")
  
  
    ## selecting data
    data <- data_ok %>% dplyr::filter(sp_name == species_list[i]) %>%
                        dplyr::filter(!is.na(x) & !is.na(y) & x >= 10 & y > 0) %>%
                        dplyr::select(x, y, location, protocol, id) %>%
                        dplyr::mutate(x = as.numeric(x), y = as.numeric(y), 
                                      location = as.factor(droplevels.factor(location)), 
                                      id = as.numeric(id), protocol = as.factor(droplevels.factor(protocol)))
  
  
  ## classifying data based on dbh classes 
  ranged_data <- data_in_class(data)
  
  ## defining sample size for resampling
  sample_size <- what_sample_size(ranged_data)
  
  ## computing nb of datasets in which the species was surveyed
  nb_datasets_all <- length(unique(ranged_data$protocol))
  
  ## sampling data and running the models for each repetition (no competition - resampling)
  output_asympt_height_nocomp <- mod_height_asympt(ranged_data, nb_datasets_all, sample_size, 
                                                      asymptot_resampling_nocomp, asymptot_resampling_nocomp_w, 
                                                      n_repetition)
  
  ## exporting results in .csv files
  write.csv(output_asympt_height_nocomp, file =  paste0("output/height_asympt_nocomp_", height_species, ".csv"))
  
  }
  
}
  
  
height_models_power_nocomp <- function(height_data, height_species) {
  
  data_ok <- height_data
  species_list <- height_species
  
  for (i in 1:length(species_list)) {
    
  height_species <- species_list[i]
  print(i)
  
  ## defining number of repetitions
  n_repetition = 200
  
  ## creating storage files
  power_resampling_nocomp <- as.data.frame(matrix(nrow = n_repetition, ncol = length(unique(data_ok$protocol)) + 6)) 
  power_resampling_nocomp[,1] <- rep(height_species, n_repetition)
  names(power_resampling_nocomp) <- c("species", paste0("a1.protocol", unique(data_ok$protocol)), "a1", "a2", "AIC", "RMSE", "weighted")
  
  power_resampling_nocomp_w <- as.data.frame(matrix(nrow = n_repetition, ncol = length(unique(data_ok$protocol)) + 6)) 
  power_resampling_nocomp_w[,1] <- rep(height_species, n_repetition)
  names(power_resampling_nocomp_w) <- c("species", paste0("a1.protocol", unique(data_ok$protocol)), "a1", "a2", "AIC", "RMSE", "weighted")
  
  ## selecting data
  data <- data_ok %>% dplyr::filter(sp_name == species_list[i]) %>%
                      dplyr::filter(!is.na(x) & !is.na(y) & x >= 10 & y > 0) %>%
                      dplyr::select(x, y, location, protocol, id) %>%
                      dplyr::mutate(x = as.numeric(x), y = as.numeric(y), 
                                    location = as.factor(droplevels.factor(location)), 
                                    id = as.numeric(id), protocol = as.factor(droplevels.factor(protocol)))
  
  
  ## classifying data based on dbh classes 
  ranged_data <- data_in_class(data)
  
  ## defining sample size for resampling
  sample_size <- what_sample_size(ranged_data)
  
  ## computing nb of datasets in which the species was surveyed
  nb_datasets_all <- length(unique(ranged_data$protocol))
  
  ## sampling data and running the models for each repetition (no competition - resampling)
  output_power_height_nocomp <- mod_height_power_nocomp(ranged_data, nb_datasets_all, sample_size, 
                                                        power_resampling_nocomp, power_resampling_nocomp_w, 
                                                        n_repetition)
  
  ## exporting results in .csv files
  write.csv(output_power_height_nocomp, file =  paste0("output/height_power_nocomp_", height_species, ".csv"))
  
  }
  
}
  
  
  
  
  
height_models_power_comp <- function(height_data, height_species_comp) {
    
    ## loading data and species list
    data_ok <- height_data

    ## defining i
    for (i in 1:length(height_species_comp)) {
      
      new_sp_list <- height_species_comp[i]
      print(i)
      
      
      ## defining number of repetitions
      n_repetition = 200
      
      ## creating storage files
      power_resampling_c1 <- as.data.frame(matrix(nrow = n_repetition, ncol = length(unique(data_ok$protocol)) + 8)) 
      power_resampling_c1[,1] <- rep(new_sp_list, n_repetition)
      names(power_resampling_c1) <- c("species", paste0("a1.protocol", unique(data_ok$protocol)), "a1", "a2", "comp", "AIC", "RMSE", "weighted", "condition")
      
      power_resampling_c1_w <- as.data.frame(matrix(nrow = n_repetition, ncol = length(unique(data_ok$protocol)) + 8)) 
      power_resampling_c1_w[,1] <- rep(new_sp_list, n_repetition)
      names(power_resampling_c1_w) <- c("species", paste0("a1.protocol", unique(data_ok$protocol)), "a1", "a2", "comp", "AIC", "RMSE", "weighted", "condition")
      
      power_resampling_c2 <- as.data.frame(matrix(nrow = n_repetition, ncol = length(unique(data_ok$protocol)) + 8)) 
      power_resampling_c2[,1] <- rep(new_sp_list, n_repetition)
      names(power_resampling_c2) <- c("species", paste0("a1.protocol", unique(data_ok$protocol)), "a1", "a2", "comp", "AIC", "RMSE", "weighted", "condition")
      
      power_resampling_c2_w <- as.data.frame(matrix(nrow = n_repetition, ncol = length(unique(data_ok$protocol)) + 8)) 
      power_resampling_c2_w[,1] <- rep(new_sp_list, n_repetition)
      names(power_resampling_c2_w) <- c("species", paste0("a1.protocol", unique(data_ok$protocol)), "a1", "a2", "comp", "AIC", "RMSE", "weighted", "condition")
      
      
      ## selecting data
      data <- data_ok %>% dplyr::filter(sp_name == height_species_comp[i]) %>%
                          dplyr::filter(!is.na(x) & !is.na(y) & x >= 10 & y > 0 & !is.na(ba_plot) & !is.na(ba_larger) & ba_plot > 0 & ba_larger > 0) %>%
                          dplyr::select(x, y, location, protocol, id, ba_plot, ba_larger) %>%
                          dplyr::mutate(x = as.numeric(x), y = as.numeric(y), 
                                        location = as.factor(droplevels.factor(location)), 
                                        protocol = as.factor(droplevels.factor(protocol)), 
                                        id = as.numeric(id), ba_plot = as.numeric(ba_plot), ba_larger = as.numeric(ba_larger))
      
      
      ## classifying data based on dbh classes 
      ranged_data <- data_in_class(data)
      
      ## defining sample size for resampling
      sample_size <- what_sample_size(ranged_data)
      
      ## computing nb of datasets in which the species was surveyed
      nb_datasets_all <- length(unique(ranged_data$protocol))
      
      ## sampling data and running the models for each repetition (no competition - resampling)
      output_power_height_comp <- mod_height_power_comp(ranged_data, nb_datasets_all, sample_size, 
                                                        power_resampling_c1, power_resampling_c1_w, 
                                                        power_resampling_c2, power_resampling_c2_w, 
                                                        n_repetition)
      
      ## exporting results in .csv files
      write.csv(output_power_height_comp, file =  paste0("output/height_power_comp_", height_species_comp[i], ".csv"))
      
    }
    
}
    
    


  
  
  

