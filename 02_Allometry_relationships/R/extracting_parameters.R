extracting_nocomp_parameters <- function() {
  
  ### define names of files to combine to obtain one file per crown caracteristic
  file_list <- c("height_asympt_nocomp_",
                 "height_power_nocomp_",
                 "diameter_nocomp_",
                 "ratio_nocomp_")

  for (i in 1:length(file_list)) {
    
    nrep = 200
    
    ### combine outputs into one file
    dd <- list.files(path = "output/", pattern = file_list[i])
    data_list <- lapply(paste0("output/", dd), utils::read.table,
                        header = TRUE, sep = ",", dec = ".",
                        encoding = "UTF-8",
                        stringsAsFactors = FALSE)
    
    df <- data.table::rbindlist(data_list, use.names = TRUE, fill = TRUE)
    
    ### filter data based on the number of models that successfully converged
    ##  and on the way the parameters were computed (i.e. weighted vs. non weighted) 
    
    df_summary <- df %>% dplyr::filter(weighted == "yes") %>%
                         dplyr::group_by(species) %>%
                         dplyr::summarise(nobs = n()) %>%
                         dplyr::filter(nobs >= 100) %>%
                         dplyr::ungroup()
    
    if (file_list[i] == "height_asympt_nocomp_") {
    
    weighted_parameters <- df %>% dplyr::filter(species %in% df_summary$species) %>%
                                  dplyr::filter(weighted == "yes") %>%
                                  dplyr::group_by(species) %>%
                                  dplyr::slice_sample(n = 100) %>%
                                  dplyr::ungroup() %>%
                                  dplyr::select(X, species, b1, b2, b3, AIC, RMSE)
    
    ### select based on rows_id to keep parameters obtained from the same model 
    ## (i.e. same subset of data used to compute weighted and unweighted parameters)
    ## 150 repetitions have been performed for each species and each model
    
    rows_id <- weighted_parameters %>% dplyr::group_by(species) %>%
                                       dplyr::mutate(to_select = X - nrep) %>%
                                       dplyr::mutate(to_select = paste0(species, "_", to_select)) %>%
                                       dplyr::ungroup()
      
    
    no_weighted_parameters <- df %>% dplyr::group_by(species) %>%
                                     dplyr::filter(weighted == "no") %>%
                                     dplyr::mutate(to_select = paste0(species, "_", X)) %>%
                                     dplyr::select(species, b1, b2, b3, AIC, RMSE, to_select) %>%
                                     dplyr::ungroup()
    
    no_weighted_parameters <- no_weighted_parameters[no_weighted_parameters$to_select %in% rows_id$to_select,]
    
    } else {
      
      weighted_parameters <- df %>% dplyr::filter(species %in% df_summary$species) %>%
                                    dplyr::filter(weighted == "yes") %>%
                                    dplyr::group_by(species) %>%
                                    dplyr::slice_sample(n = 100) %>%
                                    dplyr::ungroup() %>%
                                    dplyr::select(X, species, a1, a2, AIC, RMSE)
      
      rows_id <- weighted_parameters %>% dplyr::group_by(species) %>%
                                         dplyr::mutate(to_select = X - nrep) %>%
                                         dplyr::mutate(to_select = paste0(species, "_", to_select)) %>%
                                         dplyr::ungroup()
      
      no_weighted_parameters <- df %>% dplyr::group_by(species) %>%
                                       dplyr::filter(weighted == "no") %>%
                                       dplyr::mutate(to_select = paste0(species, "_", X)) %>%
                                       dplyr::select(species, a1, a2, AIC, RMSE, to_select) %>%
                                       dplyr::ungroup()
      
      no_weighted_parameters <- no_weighted_parameters[no_weighted_parameters$to_select %in% rows_id$to_select,]
      
    }
    
    ### create files to store results regarding AIC and RMSE at the species and model levels
    RMSE_AIC <- as.data.frame(matrix (nrow = length(unique(weighted_parameters$species)), ncol = 7))
    colnames(RMSE_AIC) <- c("species", 
                            "total_RMSE_weighted", "mean_RMSE_weighted", 
                            "total_RMSE_no_weighted", "mean_RMSE_no_weighted",
                            "total_AIC", "mean_AIC")
    
    for (j in 1:length(unique(weighted_parameters$species))) {
      
      subset_weighted <- weighted_parameters %>% dplyr::filter(species == unique(weighted_parameters$species)[j])
      subset_no_weighted <- no_weighted_parameters %>% dplyr::filter(species == unique(weighted_parameters$species)[j])
      
      RMSE_AIC[j, "species"] <- unique(subset_weighted$species)
      RMSE_AIC[j, "total_RMSE_weighted"] <- sum(subset_weighted$RMSE)
      RMSE_AIC[j, "mean_RMSE_weighted"] <- mean(subset_weighted$RMSE)
      RMSE_AIC[j, "total_RMSE_no_weighted"] <- sum(subset_no_weighted$RMSE)
      RMSE_AIC[j, "mean_RMSE_no_weighted"] <- mean(subset_no_weighted$RMSE)
      RMSE_AIC[j, "total_AIC"] <- sum(subset_weighted$AIC)
      RMSE_AIC[j, "mean_AIC"] <- mean(subset_weighted$AIC)
      
    }
    
    ### removing unnecessary columns before saving files
    weighted_parameters <- weighted_parameters %>% dplyr::select(-X, -AIC, -RMSE)
    no_weighted_parameters <- no_weighted_parameters %>% dplyr::select(-AIC, -RMSE, -to_select)
    
    ### saving files
    write.csv(file = paste0("output/", file_list[i], "weighted_parameters.csv"), weighted_parameters)
    write.csv(file = paste0("output/", file_list[i], "no_weighted_parameters.csv"), no_weighted_parameters)
    write.csv(file = paste0("output/", file_list[i], "RMSE_AIC_results.csv"), RMSE_AIC)
    
  }
  
}




extracting_comp_parameters <- function() {
  
  ### define names of files to combine to obtain one file per crown caracteristic
  file_list <- c("height_power_comp_",
                 "diameter_comp_",
                 "ratio_comp")
  
  for (i in 1:length(file_list)) {
    
    nrep = 200
    
    ### combine outputs into one file
    dd <- list.files(path = "output/", pattern = file_list[i])
    data_list <- lapply(paste0("output/", dd), utils::read.table,
                        header = TRUE, sep = ",", dec = ".",
                        encoding = "UTF-8",
                        stringsAsFactors = FALSE)
    
    df <- data.table::rbindlist(data_list, use.names = TRUE, fill = TRUE)
    
    ### filter data based on the number of models that successfully converged
    ## the indice of competition that was included into the model (i.e. BAT or BAL)
    ##  and on the way the parameters were computed (i.e. weighted vs. non weighted) 
    
    df_summary <- df %>% dplyr::filter(weighted == "yes" & condition == "c1") %>%
                         dplyr::group_by(species) %>%
                         dplyr::summarise(nobs = n()) %>%
                         dplyr::filter(nobs >= 100) %>%
                         dplyr::ungroup()
    
      weighted_parameters_c1 <- df %>% dplyr::filter(species %in% df_summary$species) %>%
                                       dplyr::filter(weighted == "yes" & condition == "c1") %>%
                                       dplyr::group_by(species) %>%
                                       dplyr::slice_sample(n = 100) %>%
                                       dplyr::ungroup() %>%
                                       dplyr::select(X, species, a1, a2, comp, AIC, RMSE)

      rows_id_c1_no_weighted <- weighted_parameters_c1 %>% dplyr::group_by(species) %>%
                                                           dplyr::mutate(to_select = X - nrep) %>%
                                                           dplyr::mutate(to_select = paste0(species, "_", to_select)) %>%
                                                           dplyr::ungroup()
      
      rows_id_c2_no_weighted <- weighted_parameters_c1 %>% dplyr::group_by(species) %>%
                                                           dplyr::mutate(to_select = X + nrep) %>%
                                                           dplyr::mutate(to_select = paste0(species, "_", to_select)) %>%
                                                           dplyr::ungroup()
      
      rows_id_c2_weighted <- weighted_parameters_c1 %>% dplyr::group_by(species) %>%
                                                        dplyr::mutate(to_select = X + (2 * nrep)) %>%
                                                        dplyr::mutate(to_select = paste0(species, "_", to_select)) %>%
                                                        dplyr::ungroup()
      
      
      no_weighted_parameters_c1 <- df %>% dplyr::group_by(species) %>%
                                          dplyr::filter(weighted == "no" & condition == "c1") %>%
                                          dplyr::mutate(to_select = paste0(species, "_", X)) %>%
                                          dplyr::select(species, a1, a2, comp, AIC, RMSE, to_select) %>%
                                          dplyr::ungroup()
      
      no_weighted_parameters_c2 <- df %>% dplyr::group_by(species) %>%
                                          dplyr::filter(weighted == "no" & condition == "c2") %>%
                                          dplyr::mutate(to_select = paste0(species, "_", X)) %>%
                                          dplyr::select(species, a1, a2, comp, AIC, RMSE, to_select) %>%
                                          dplyr::ungroup()
      
      weighted_parameters_c2 <- df %>% dplyr::group_by(species) %>%
                                       dplyr::filter(weighted == "yes" & condition == "c2") %>%
                                       dplyr::mutate(to_select = paste0(species, "_", X)) %>%
                                       dplyr::select(species, a1, a2, comp, AIC, RMSE, to_select) %>%
                                       dplyr::ungroup()
                                      
      no_weighted_parameters_c1 <- no_weighted_parameters_c1[no_weighted_parameters_c1$to_select %in% rows_id_c1_no_weighted$to_select,]
      no_weighted_parameters_c2 <- no_weighted_parameters_c2[no_weighted_parameters_c2$to_select %in% rows_id_c2_no_weighted$to_select,]
      weighted_parameters_c2 <- weighted_parameters_c2[weighted_parameters_c2$to_select %in% rows_id_c2_weighted$to_select,]
      
    
      ### create files to store results regarding AIC and RMSE at the species and model levels
      RMSE_AIC <- as.data.frame(matrix (nrow = length(unique(weighted_parameters_c1$species)), ncol = 13))
      colnames(RMSE_AIC) <- c("species", 
                              "c1_total_RMSE_weighted", "c1_mean_RMSE_weighted", 
                              "c2_total_RMSE_weighted", "c2_mean_RMSE_weighted", 
                              "c1_total_RMSE_no_weighted", "c1_mean_RMSE_no_weighted",
                              "c2_total_RMSE_no_weighted", "c2_mean_RMSE_no_weighted",
                              "c1_total_AIC", "c1_mean_AIC",
                              "c2_total_AIC", "c2_mean_AIC")
    
    for (j in 1:length(unique(weighted_parameters_c1$species))) {
      
      subset_weighted_c1 <- weighted_parameters_c1 %>% dplyr::filter(species == unique(weighted_parameters_c1$species)[j])
      subset_no_weighted_c1 <- no_weighted_parameters_c1 %>% dplyr::filter(species == unique(weighted_parameters_c1$species)[j])
      subset_weighted_c2 <- weighted_parameters_c2 %>% dplyr::filter(species == unique(weighted_parameters_c1$species)[j])
      subset_no_weighted_c2 <- no_weighted_parameters_c2 %>% dplyr::filter(species == unique(weighted_parameters_c1$species)[j])
      
      RMSE_AIC[j, "species"] <- unique(subset_weighted_c1$species)
      
      RMSE_AIC[j, "c1_total_RMSE_weighted"] <- sum(subset_weighted_c1$RMSE)
      RMSE_AIC[j, "c1_mean_RMSE_weighted"] <- mean(subset_weighted_c1$RMSE)
      RMSE_AIC[j, "c1_total_RMSE_no_weighted"] <- sum(subset_no_weighted_c1$RMSE)
      RMSE_AIC[j, "c1_mean_RMSE_no_weighted"] <- mean(subset_no_weighted_c1$RMSE)
      RMSE_AIC[j, "c1_total_AIC"] <- sum(subset_weighted_c1$AIC)
      RMSE_AIC[j, "c1_mean_AIC"] <- mean(subset_weighted_c1$AIC)
      
      RMSE_AIC[j, "c2_total_RMSE_weighted"] <- sum(subset_weighted_c2$RMSE)
      RMSE_AIC[j, "c2_mean_RMSE_weighted"] <- mean(subset_weighted_c2$RMSE)
      RMSE_AIC[j, "c2_total_RMSE_no_weighted"] <- sum(subset_no_weighted_c2$RMSE)
      RMSE_AIC[j, "c2_mean_RMSE_no_weighted"] <- mean(subset_no_weighted_c2$RMSE)
      RMSE_AIC[j, "c2_total_AIC"] <- sum(subset_weighted_c2$AIC)
      RMSE_AIC[j, "c2_mean_AIC"] <- mean(subset_weighted_c2$AIC)
      
    }
    
    ### removing unnecessary columns before saving files
    weighted_parameters_c1 <- weighted_parameters_c1 %>% dplyr::select(-X, -AIC, -RMSE)
    no_weighted_parameters_c1 <- no_weighted_parameters_c1 %>% dplyr::select(-AIC, -RMSE, -to_select)
    weighted_parameters_c2 <- weighted_parameters_c2 %>% dplyr::select(-AIC, -RMSE, -to_select)
    no_weighted_parameters_c2 <- no_weighted_parameters_c2 %>% dplyr::select(-AIC, -RMSE, -to_select)
    
    ### saving files
    write.csv(file = paste0("output/", file_list[i], "weighted_parameters_c1.csv"), weighted_parameters_c1)
    write.csv(file = paste0("output/", file_list[i], "no_weighted_parameters_c1.csv"), no_weighted_parameters_c1)
    write.csv(file = paste0("output/", file_list[i], "weighted_parameters_c2.csv"), weighted_parameters_c2)
    write.csv(file = paste0("output/", file_list[i], "no_weighted_parameters_c2.csv"), no_weighted_parameters_c2)
    write.csv(file = paste0("output/", file_list[i], "RMSE_AIC_results.csv"), RMSE_AIC)
    
  }
  
}

