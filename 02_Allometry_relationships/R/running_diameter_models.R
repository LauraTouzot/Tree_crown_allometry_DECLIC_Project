###### Running diameter models #######

diameter_models_nocomp <- function(diameter_data, diameter_species) {

  ## loading data and species list
  data_ok <- diameter_data
  species_list <- diameter_species

  ## defining i
  for (i in 1:length(species_list)) {
    
    diameter_species <- species_list[i]
    print(i)
    
    
    ## defining number of repetitions
    n_repetition = 200
    
    ## creating storage files
    power_resampling_nocomp <- as.data.frame(matrix(nrow = n_repetition, ncol = length(unique(data_ok$protocol)) + 6)) 
    power_resampling_nocomp[,1] <- rep(diameter_species, n_repetition)
    names(power_resampling_nocomp) <- c("species", paste0("a1.protocol", unique(data_ok$protocol)), "a1", "a2", "AIC", "RMSE", "weighted")
    
    power_resampling_nocomp_w <- as.data.frame(matrix(nrow = n_repetition, ncol = length(unique(data_ok$protocol)) + 6)) 
    power_resampling_nocomp_w[,1] <- rep(diameter_species, n_repetition)
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
    output_power_diameter_nocomp <- mod_diameter_nocomp(ranged_data, nb_datasets_all, sample_size, 
                                                        power_resampling_nocomp, power_resampling_nocomp_w, 
                                                        n_repetition)
    
    ## exporting results in .csv files
    write.csv(output_power_diameter_nocomp, file =  paste0("output/diameter_nocomp_", diameter_species, ".csv"))
    
  }
  
}
  
  
  
  
  
  
diameter_models_comp <- function(diameter_data, diameter_species_comp) {
    
    ## loading data and species list
    data_ok <- diameter_data
    
    ## defining i
    for (i in 1:length(diameter_species_comp)) {
      
      new_sp_list <- diameter_species_comp[i]
      print(i)
      
      
      ## defining number of repetitions
      n_repetition = 200
      
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
      data <- data_ok %>% dplyr::filter(sp_name == diameter_species_comp[i]) %>%
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
      output_power_diameter_comp <- mod_diameter_comp(ranged_data, nb_datasets_all, sample_size, 
                                                      power_resampling_c1, power_resampling_c1_w, 
                                                      power_resampling_c2, power_resampling_c2_w, 
                                                      n_repetition)
      
      ## exporting results in .csv files
      write.csv(output_power_diameter_comp, file =  paste0("output/diameter_comp_", diameter_species_comp[i], ".csv"))
      
    }
    
}
    






