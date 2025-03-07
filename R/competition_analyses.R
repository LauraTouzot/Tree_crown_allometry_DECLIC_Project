### ANALYSES INCLUDING LOCAL COMPETITION

species_selection_competition <- function() {
  
  ### Defining species to remove (i.e., for which models were not robust - abnormal response parameters)
  to_remove <- c("Alnus rubra", "Pinus nigra", "Pinus resinosa", "Quercus stellata", "Quercus chrysolepis", "Pinus pinaster")
  
  ### Selecting species for which the three studied crown characteristics were available in presence of competition
  height <- read.csv(file = "data/height_power_comp_weighted_parameters_c2.csv")
  height <- height %>% dplyr::group_by(species) %>% 
                       dplyr::filter(!species %in% to_remove) %>%
                       dplyr::slice_sample(n = 1) %>%
                       dplyr::ungroup() %>%
                       dplyr::select(species, a1) 
  
  ratio <- read.csv(file = "data/ratio_compweighted_parameters_c2.csv")
  ratio <- ratio %>% dplyr::group_by(species) %>%
                     dplyr::slice_sample(n = 1) %>%
                     dplyr::ungroup() %>%
                     dplyr::select(species, a1)
  
  diameter <- read.csv(file = "data/diameter_comp_weighted_parameters_c2.csv")
  diameter <- diameter %>% dplyr::group_by(species) %>%
                           dplyr::slice_sample(n = 1) %>%
                           dplyr::ungroup() %>%
                           dplyr::select(species, a1)
  
  species_competition <- join_all(list(height, ratio, diameter),
                                  by = "species",
                                  type = "left")
  
  species_competition <- na.omit(species_competition)
  species_list_competition <- unique(species_competition$species)
  
  return(species_list_competition)
  
}
  
 

response_to_competition <- function(species_list_comp) { 
  
  ## BAL has been selected as the best index to account for local competition based on AIC values
  ## Note that weighted parameters have been selected based on RMSE values
  
  ## Height 
  height <- read.csv(file = "data/height_power_comp_weighted_parameters_c2.csv")
  height <- height %>% dplyr::select(-X) %>%
                       dplyr::rename(a1_height = a1,
                                     a2_height = a2,
                                     comp_height = comp) %>%
                       dplyr::filter(species %in% species_list_comp) %>%
                       dplyr::arrange(species)
    
  
  ## Relative crown depth
  ratio <- read.csv(file = "data/ratio_compweighted_parameters_c2.csv")
  ratio <- ratio %>% dplyr::select(-X) %>%
                     dplyr::rename(a1_ratio = a1,
                                   a2_ratio = a2,
                                   comp_ratio = comp) %>%
                     dplyr::filter(species %in% species_list_comp) %>%
                     dplyr::arrange(species) %>%
                     dplyr::select(-species)
  
  
  ## Crown diameter
  diameter <- read.csv(file = "data/diameter_comp_weighted_parameters_c2.csv")
  diameter <- diameter %>% dplyr::select(-X) %>%
                           dplyr::rename(a1_diameter = a1,
                                         a2_diameter = a2,
                                         comp_diameter = comp) %>%
                           dplyr::filter(species %in% species_list_comp) %>%
                           dplyr::arrange(species) %>%
                           dplyr::select(-species)
  
  
  ## All crown characteristics
  comp_parameters <- cbind(height, ratio, diameter)
  
  
  ### Scaling and centering crown characteristics at the functional group scale
  ## Summarizing crown characteristics at the species scale
  
  functional_groups <- read.csv(file = "data/functional_groups.csv")
  functional_groups <- functional_groups %>% dplyr::select(sp, new_gp)

  comp_parameters <- left_join(comp_parameters, functional_groups, by = c("species" = "sp"))  

  angio_response <- comp_parameters %>% dplyr::filter(new_gp == "angio") %>%
    
                                        dplyr::mutate(height_resp = scale(comp_height, center = TRUE, scale = TRUE)[, 1],
                                                      ratio_resp = scale(comp_ratio, center = TRUE, scale = TRUE)[, 1],
                                                      diameter_resp = scale(comp_diameter, center = TRUE, scale = TRUE)[, 1]) %>%
    
                                        dplyr::group_by(species) %>%
                                        dplyr::summarise(mean_height_par = mean(height_resp),
                                                         cv_height_par = abs(sd(height_resp)/mean(height_resp)),
                                                         se_height_par = sd(height_resp) / sqrt(n()),
                                                         
                                                         mean_ratio_par = mean(ratio_resp),
                                                         cv_ratio_par = abs(sd(ratio_resp)/mean(ratio_resp)),
                                                         se_ratio_par = sd(ratio_resp) / sqrt(n()),
                                                         
                                                         mean_diam_par = mean(diameter_resp),
                                                         cv_diam_par = abs(sd(diameter_resp)/mean(diameter_resp)),
                                                         se_diam_par = sd(diameter_resp) / sqrt(n())) %>%
                                        
                                        dplyr::mutate(lower_ci_height = mean_height_par - qt(1 - (0.05 / 2), n() - 1) * se_height_par,
                                                      upper_ci_height = mean_height_par + qt(1 - (0.05 / 2), n() - 1) * se_height_par,
                                                      
                                                      lower_ci_ratio = mean_ratio_par - qt(1 - (0.05 / 2), n() - 1) * se_ratio_par,
                                                      upper_ci_ratio = mean_ratio_par + qt(1 - (0.05 / 2), n() - 1) * se_ratio_par,
                                                      
                                                      lower_ci_diam = mean_diam_par - qt(1 - (0.05 / 2), n() - 1) * se_diam_par,
                                                      upper_ci_diam = mean_diam_par + qt(1 - (0.05 / 2), n() - 1) * se_diam_par) %>%
                                        dplyr::ungroup()
  
  
  gymno_response <- comp_parameters %>% dplyr::filter(new_gp == "gymno") %>%
    
                                        dplyr::mutate(height_resp = scale(comp_height, center = TRUE, scale = TRUE)[, 1],
                                                      ratio_resp = scale(comp_ratio, center = TRUE, scale = TRUE)[, 1],
                                                      diameter_resp = scale(comp_diameter, center = TRUE, scale = TRUE)[, 1]) %>%
                                        
                                        dplyr::group_by(species) %>%
                                        dplyr::summarise(mean_height_par = mean(height_resp),
                                                         cv_height_par = abs(sd(height_resp)/mean(height_resp)),
                                                         se_height_par = sd(height_resp) / sqrt(n()),
                                                         
                                                         mean_ratio_par = mean(ratio_resp),
                                                         cv_ratio_par = abs(sd(ratio_resp)/mean(ratio_resp)),
                                                         se_ratio_par = sd(ratio_resp) / sqrt(n()),
                                                         
                                                         mean_diam_par = mean(diameter_resp),
                                                         cv_diam_par = abs(sd(diameter_resp)/mean(diameter_resp)),
                                                         se_diam_par = sd(diameter_resp) / sqrt(n())) %>%
                                        
                                        dplyr::mutate(lower_ci_height = mean_height_par - qt(1 - (0.05 / 2), n() - 1) * se_height_par,
                                                      upper_ci_height = mean_height_par + qt(1 - (0.05 / 2), n() - 1) * se_height_par,
                                                      
                                                      lower_ci_ratio = mean_ratio_par - qt(1 - (0.05 / 2), n() - 1) * se_ratio_par,
                                                      upper_ci_ratio = mean_ratio_par + qt(1 - (0.05 / 2), n() - 1) * se_ratio_par,
                                                      
                                                      lower_ci_diam = mean_diam_par - qt(1 - (0.05 / 2), n() - 1) * se_diam_par,
                                                      upper_ci_diam = mean_diam_par + qt(1 - (0.05 / 2), n() - 1) * se_diam_par) %>%
                                        dplyr::ungroup()

  
  ## Loading and merging trait data
  traits_data <- read.csv(file = "data/sp_traits_complete.csv")
  shade_tolerance <- traits_data %>% dplyr::select(species, shade_tol_mean)
  
  angio_response <- left_join(angio_response, shade_tolerance, by = "species") %>%
                    dplyr::mutate(shade_tol_scaled = scale(shade_tol_mean, center = TRUE, scale = TRUE)[, 1]) %>%
                    dplyr::select(-shade_tol_mean) %>%
                    drop_na()
  
  gymno_response <- left_join(gymno_response, shade_tolerance, by = "species") %>%
                              dplyr::mutate(shade_tol_scaled = scale(shade_tol_mean, center = TRUE, scale = TRUE)[, 1]) %>%
                              dplyr::select(-shade_tol_mean) %>%
                              drop_na()

  
  ## Running models on scaled and centered data
  mod_h_sh_a <- lm(mean_height_par ~ shade_tol_scaled, weights = (1/cv_height_par), data = angio_response)
  mod_r_sh_a <- lm(mean_ratio_par ~ shade_tol_scaled, weights = (1/cv_ratio_par), data = angio_response)
  mod_d_sh_a <- lm(mean_diam_par ~ shade_tol_scaled, weights = (1/cv_diam_par), data = angio_response)
  
  mod_h_sh_g <- lm(mean_height_par ~ shade_tol_scaled, weights = (1/cv_height_par), data = gymno_response)
  mod_r_sh_g <- lm(mean_ratio_par ~ shade_tol_scaled, weights = (1/cv_ratio_par), data = gymno_response)
  mod_d_sh_g <- lm(mean_diam_par ~ shade_tol_scaled, weights = (1/cv_diam_par), data = gymno_response)

  
  ## Saving results from models
  angio_models <- list(mod_h_sh_a,
                       mod_r_sh_a, 
                       mod_d_sh_a)
  
  gymno_models <- list(mod_h_sh_g,
                       mod_r_sh_g, 
                       mod_d_sh_g)
  
  ### Creating file to store models' results while computing the figures
  angio_response_to_competition <- as.data.frame(matrix(nrow = 3, ncol = 3))
  colnames(angio_response_to_competition) <- c("intercept", "b", "p.value")
  rownames(angio_response_to_competition) <- c("height", "ratio", "diameter")
  
  gymno_response_to_competition <- as.data.frame(matrix(nrow = 3, ncol = 3))
  colnames(gymno_response_to_competition) <- c("intercept", "b", "p.value")
  rownames(gymno_response_to_competition) <- c("height", "ratio", "diameter")

  
  ### Computing associated figures
  data_y_angio <- list(angio_response[,"mean_height_par"],
                       angio_response[,"mean_ratio_par"],
                       angio_response[,"mean_diam_par"])
  
  data_y_gymno <- list(gymno_response[,"mean_height_par"],
                       gymno_response[,"mean_ratio_par"],
                       gymno_response[,"mean_diam_par"])
  

  ### Defining variables necessary to produce the plots
  data_sh_angio <- as.data.frame(angio_response[,"shade_tol_scaled"])
  data_sh_gymno <- as.data.frame(gymno_response[,"shade_tol_scaled"])
  
  models <- c("height", "crown_ratio", "crown_diameter")
  
  
  # for (i in 1:length(angio_models)) {
  #   
  #   model_angio <- angio_models[[i]]
  #   model_gymno <- gymno_models[[i]]
  #   
  #   pdf(file = paste0("manuscript/response_to_comp_", models[[i]], ".pdf"), width = 9, height = 9, pointsize = 10)
  #   
  #   x <- seq(-2,2,0.1)
  #   
  #   ## Saving models' results in the dedicated file
  #   angio_response_to_competition[i, "intercept"] <-  summary(model_angio)$coefficients["(Intercept)", "Estimate"]
  #   angio_response_to_competition[i, "b"] <-  summary(model_angio)$coefficients["shade_tol_scaled", "Estimate"]
  #   angio_response_to_competition[i, "p.value"] <-  summary(model_angio)$coefficients["shade_tol_scaled", "Pr(>|t|)"]
  #   
  #   gymno_response_to_competition[i, "intercept"] <-  summary(model_gymno)$coefficients["(Intercept)", "Estimate"]
  #   gymno_response_to_competition[i, "b"] <-  summary(model_gymno)$coefficients["shade_tol_scaled", "Estimate"]
  #   gymno_response_to_competition[i, "p.value"] <-  summary(model_gymno)$coefficients["shade_tol_scaled", "Pr(>|t|)"]
    
    
  #   ## Creating required information 
  #   fitted_angio <- as.data.frame(predict(model_angio, newdata = data.frame(shade_tol_scaled = seq(-2,2,0.1)), interval = "confidence"))
  #   y_angio <- fitted_angio$fit
  #   ylow_angio <- fitted_angio$lwr
  #   yhigh_angio <- fitted_angio$upr
  #   
  #   fitted_gymno <- as.data.frame(predict(model_gymno, newdata = data.frame(shade_tol_scaled = seq(-2,2,0.1)), interval = "confidence"))
  #   y_gymno <- fitted_gymno$fit
  #   ylow_gymno <- fitted_gymno$lwr
  #   yhigh_gymno <- fitted_gymno$upr
  #   
  #   y_points_angio <- as.data.frame(data_y_angio[[i]])
  #   y_points_gymno <- as.data.frame(data_y_gymno[[i]])
  #   
  #   
  #   
  #   if (summary(model_angio)$coefficients["shade_tol_scaled", "Pr(>|t|)"] <= 0.05 & summary(model_gymno)$coefficients["shade_tol_scaled", "Pr(>|t|)"] <= 0.05) {
  #     
  #     plot(y_angio ~ x, lwd = 2.5, type = "l", col = "white", las = 1, xlab = "shade tolerance", ylab = models[[i]], xlim = c(-2,2), ylim = c(-3,3))
  #     
  #     polygon(c(x,rev(x)),c(yhigh_angio, rev(ylow_angio)), border = NA, col = rgb(r = 1, g = 0.4, b = 0, max = 1, alpha = 0.5))
  #     polygon(c(x,rev(x)),c(yhigh_gymno, rev(ylow_gymno)), border = NA, col = rgb(r = 0, g = 0.6, b = 0.2, max = 1, alpha = 0.5))
  #     
  #     lines(y_angio ~ x, lwd = 4, col = "#c05a00")
  #     lines(y_gymno ~ x, lwd = 4, col = "forestgreen")
  #     
  #     points(data_sh_angio[,1], y_points_angio[,1], pch = 16, col = "#c05a00", cex = 2.5)
  #     points(data_sh_gymno[,1], y_points_gymno[,1], pch = 15, col = "forestgreen", cex = 2.5)
  #     
  #     
  #   } else if (summary(model_angio)$coefficients["shade_tol_scaled", "Pr(>|t|)"] > 0.05 & summary(model_gymno)$coefficients["shade_tol_scaled", "Pr(>|t|)"] > 0.05) {
  #     
  #     plot(y_angio ~ x, lwd = 2.5, type = "l", col = "white", las = 1, xlab = "shade tolerance", ylab = models[[i]], xlim = c(-2,2), ylim = c(-3,3))
  # 
  #     y_points_angio <- as.data.frame(data_y_angio[[i]])
  #     y_points_gymno <- as.data.frame(data_y_gymno[[i]])
  #     
  #     points(data_sh_angio[,1], y_points_angio[,1], pch = 16, col = "#c05a00", cex = 1.5)
  #     points(data_sh_gymno[,1], y_points_gymno[,1], pch = 15, col = "forestgreen", cex = 1.5)
  #     
  #     
  #   } else if (summary(model_angio)$coefficients["shade_tol_scaled", "Pr(>|t|)"] <= 0.05 & summary(model_gymno)$coefficients["shade_tol_scaled", "Pr(>|t|)"] > 0.05) {
  #   
  #     plot(y_angio ~ x, lwd = 2.5, type = "l", col = "white", las = 1, xlab = "shade tolerance", ylab = models[[i]], xlim = c(-2,2), ylim = c(-3,3))
  #     
  #     polygon(c(x,rev(x)),c(yhigh_angio, rev(ylow_angio)), border = NA, col = rgb(r = 1, g = 0.4, b = 0, max = 1, alpha = 0.5))
  # 
  #     lines(y_angio ~ x, lwd = 4, col = "#c05a00")
  # 
  #     points(data_sh_angio[,1], y_points_angio[,1], pch = 16, col = "#c05a00", cex = 2.5)
  #     points(data_sh_gymno[,1], y_points_gymno[,1], pch = 15, col = "forestgreen", cex = 1.5)
  #     
  #   
  #   }  else {
  #   
  #     plot(y_angio ~ x, lwd = 2.5, type = "l", col = "white", las = 1, xlab = "shade tolerance", ylab = models[[i]], xlim = c(-2,2), ylim = c(-3,3))
  #     
  #     polygon(c(x,rev(x)),c(yhigh_gymno, rev(ylow_gymno)), border = NA, col = rgb(r = 0, g = 0.6, b = 0.2, max = 1, alpha = 0.5))
  #     
  #     lines(y_gymno ~ x, lwd = 4.5, col = "forestgreen")
  #     
  #     points(data_sh_angio[,1], y_points_angio[,1], pch = 16, col = "#c05a00", cex = 1.5)
  #     points(data_sh_gymno[,1], y_points_gymno[,1], pch = 15, col = "forestgreen", cex = 2.5)
  #     
  #   }
  #   
  #   dev.off
  #   if (!is.null(dev.list())) dev.off()
  #   
  # }
  
  angio_response_to_competition[, "group"] <- "angio"
  gymno_response_to_competition[, "group"] <- "gymno"

  response_to_competition <- rbind(angio_response_to_competition, gymno_response_to_competition)  
  write.csv(file = "output/response_to_competition.csv", response_to_competition)

  return(response_to_competition)
  
}



response_to_competition_unscaled <- function(species_list_comp) { 
  
  ## BAL has been selected as the best index to account for local competition based on AIC values
  ## Note that weighted parameters have been selected based on RMSE values
  ## See Table X and X provided in Supplementary Materials X and X, respectively
  
  ## Height 
  height <- read.csv(file = "data/height_power_comp_weighted_parameters_c2.csv")
  height <- height %>% dplyr::select(-X) %>%
                       dplyr::rename(a1_height = a1,
                                     a2_height = a2,
                                     comp_height = comp) %>%
                       dplyr::filter(species %in% species_list_comp) %>%
                       dplyr::arrange(species)
  
  
  ## Relative crown depth
  ratio <- read.csv(file = "data/ratio_compweighted_parameters_c2.csv")
  ratio <- ratio %>% dplyr::select(-X) %>%
                     dplyr::rename(a1_ratio = a1,
                                   a2_ratio = a2,
                                   comp_ratio = comp) %>%
                     dplyr::filter(species %in% species_list_comp) %>%
                     dplyr::arrange(species) %>%
                     dplyr::select(-species)
  
  
  ## Crown diameter
  diameter <- read.csv(file = "data/diameter_comp_weighted_parameters_c2.csv")
  diameter <- diameter %>% dplyr::select(-X) %>%
                           dplyr::rename(a1_diameter = a1,
                                         a2_diameter = a2,
                                         comp_diameter = comp) %>%
                           dplyr::filter(species %in% species_list_comp) %>%
                           dplyr::arrange(species) %>%
                           dplyr::select(-species)
  
  
  ## All crown characteristics
  comp_parameters <- cbind(height, ratio, diameter)
  
  
  ### Summarizing crown characteristics at the species scale
  functional_groups <- read.csv(file = "data/functional_groups.csv")
  functional_groups <- functional_groups %>% dplyr::select(sp, new_gp)
  
  comp_parameters <- left_join(comp_parameters, functional_groups, by = c("species" = "sp"))  
  
  angio_response <- comp_parameters %>% dplyr::filter(new_gp == "angio") %>%
    
                                        dplyr::rename(height_resp = comp_height,
                                                      ratio_resp = comp_ratio,
                                                      diameter_resp = comp_diameter) %>%
                                        
                                        dplyr::group_by(species) %>%
                                        dplyr::summarise(mean_height_par = mean(height_resp),
                                                         cv_height_par = abs(sd(height_resp)/mean(height_resp)),
                                                         se_height_par = sd(height_resp) / sqrt(n()),
                                                         
                                                         mean_ratio_par = mean(ratio_resp),
                                                         cv_ratio_par = abs(sd(ratio_resp)/mean(ratio_resp)),
                                                         se_ratio_par = sd(ratio_resp) / sqrt(n()),
                                                         
                                                         mean_diam_par = mean(diameter_resp),
                                                         cv_diam_par = abs(sd(diameter_resp)/mean(diameter_resp)),
                                                         se_diam_par = sd(diameter_resp) / sqrt(n())) %>%
                                        
                                        dplyr::mutate(lower_ci_height = mean_height_par - qt(1 - (0.05 / 2), n() - 1) * se_height_par,
                                                      upper_ci_height = mean_height_par + qt(1 - (0.05 / 2), n() - 1) * se_height_par,
                                                      
                                                      lower_ci_ratio = mean_ratio_par - qt(1 - (0.05 / 2), n() - 1) * se_ratio_par,
                                                      upper_ci_ratio = mean_ratio_par + qt(1 - (0.05 / 2), n() - 1) * se_ratio_par,
                                                      
                                                      lower_ci_diam = mean_diam_par - qt(1 - (0.05 / 2), n() - 1) * se_diam_par,
                                                      upper_ci_diam = mean_diam_par + qt(1 - (0.05 / 2), n() - 1) * se_diam_par) %>%
                                        dplyr::ungroup() 
                                      
  
  gymno_response <- comp_parameters %>% dplyr::filter(new_gp == "gymno") %>%
    
                                        dplyr::rename(height_resp = comp_height,
                                                      ratio_resp = comp_ratio,
                                                      diameter_resp = comp_diameter) %>%
                                        
                                        dplyr::group_by(species) %>%
                                        dplyr::summarise(mean_height_par = mean(height_resp),
                                                         cv_height_par = abs(sd(height_resp)/mean(height_resp)),
                                                         se_height_par = sd(height_resp) / sqrt(n()),
                                                         
                                                         mean_ratio_par = mean(ratio_resp),
                                                         cv_ratio_par = abs(sd(ratio_resp)/mean(ratio_resp)),
                                                         se_ratio_par = sd(ratio_resp) / sqrt(n()),
                                                         
                                                         mean_diam_par = mean(diameter_resp),
                                                         cv_diam_par = abs(sd(diameter_resp)/mean(diameter_resp)),
                                                         se_diam_par = sd(diameter_resp) / sqrt(n())) %>%
                                        
                                        dplyr::mutate(lower_ci_height = mean_height_par - qt(1 - (0.05 / 2), n() - 1) * se_height_par,
                                                      upper_ci_height = mean_height_par + qt(1 - (0.05 / 2), n() - 1) * se_height_par,
                                                      
                                                      lower_ci_ratio = mean_ratio_par - qt(1 - (0.05 / 2), n() - 1) * se_ratio_par,
                                                      upper_ci_ratio = mean_ratio_par + qt(1 - (0.05 / 2), n() - 1) * se_ratio_par,
                                                      
                                                      lower_ci_diam = mean_diam_par - qt(1 - (0.05 / 2), n() - 1) * se_diam_par,
                                                      upper_ci_diam = mean_diam_par + qt(1 - (0.05 / 2), n() - 1) * se_diam_par) %>%
                                        dplyr::ungroup() 
  
  
  ## Loading and merging trait data
  traits_data <- read.csv(file = "data/sp_traits_complete.csv")
  shade_tolerance <- traits_data %>% dplyr::select(species, shade_tol_mean)
  
  angio_response <- left_join(angio_response, shade_tolerance, by = "species") %>% drop_na()
  gymno_response <- left_join(gymno_response, shade_tolerance, by = "species") %>% drop_na()
  
  
  ## Running models 
  mod_h_sh_a <- lm(mean_height_par ~ shade_tol_mean, weights = (1/cv_height_par), data = angio_response)
  mod_r_sh_a <- lm(mean_ratio_par ~ shade_tol_mean, weights = (1/cv_ratio_par), data = angio_response)
  mod_d_sh_a <- lm(mean_diam_par ~ shade_tol_mean, weights = (1/cv_diam_par), data = angio_response)
  
  mod_h_sh_g <- lm(mean_height_par ~ shade_tol_mean, weights = (1/cv_height_par), data = gymno_response)
  mod_r_sh_g <- lm(mean_ratio_par ~ shade_tol_mean, weights = (1/cv_ratio_par), data = gymno_response)
  mod_d_sh_g <- lm(mean_diam_par ~ shade_tol_mean, weights = (1/cv_diam_par), data = gymno_response)
  
  
  ## Saving results from models
  angio_models <- list(mod_h_sh_a,
                       mod_r_sh_a, 
                       mod_d_sh_a)
  
  gymno_models <- list(mod_h_sh_g,
                       mod_r_sh_g, 
                       mod_d_sh_g)
  
  ### Creating file to store models' results while computing the figures
  angio_response_to_competition <- as.data.frame(matrix(nrow = 3, ncol = 3))
  colnames(angio_response_to_competition) <- c("intercept", "b", "p.value")
  rownames(angio_response_to_competition) <- c("height", "ratio", "diameter")
  
  gymno_response_to_competition <- as.data.frame(matrix(nrow = 3, ncol = 3))
  colnames(gymno_response_to_competition) <- c("intercept", "b", "p.value")
  rownames(gymno_response_to_competition) <- c("height", "ratio", "diameter")
  
  
  ### Computing associated figures
  data_y_angio <- list(angio_response[,"mean_height_par"],
                       angio_response[,"mean_ratio_par"],
                       angio_response[,"mean_diam_par"])
  
  data_y_gymno <- list(gymno_response[,"mean_height_par"],
                       gymno_response[,"mean_ratio_par"],
                       gymno_response[,"mean_diam_par"])
  
  
  ### Defining variables necessary to produce the plots
  data_sh_angio <- as.data.frame(angio_response[,"shade_tol_mean"])
  data_sh_gymno <- as.data.frame(gymno_response[,"shade_tol_mean"])
  
  models <- c("height", "crown_ratio", "crown_diameter")
  
  
  for (i in 1:length(angio_models)) {
    
    model_angio <- angio_models[[i]]
    model_gymno <- gymno_models[[i]]
    
    pdf(file = paste0("manuscript/response_to_comp_unscaled_", models[[i]], ".pdf"), width = 9, height = 9, pointsize = 10)
    
    x <- seq(1,5.5,0.1)
    
    ## Saving models' results in the dedicated file
    angio_response_to_competition[i, "intercept"] <-  summary(model_angio)$coefficients["(Intercept)", "Estimate"]
    angio_response_to_competition[i, "b"] <-  summary(model_angio)$coefficients["shade_tol_mean", "Estimate"]
    angio_response_to_competition[i, "p.value"] <-  summary(model_angio)$coefficients["shade_tol_mean", "Pr(>|t|)"]
    
    gymno_response_to_competition[i, "intercept"] <-  summary(model_gymno)$coefficients["(Intercept)", "Estimate"]
    gymno_response_to_competition[i, "b"] <-  summary(model_gymno)$coefficients["shade_tol_mean", "Estimate"]
    gymno_response_to_competition[i, "p.value"] <-  summary(model_gymno)$coefficients["shade_tol_mean", "Pr(>|t|)"]
    
    
    ## Creating required information 
    fitted_angio <- as.data.frame(predict(model_angio, newdata = data.frame(shade_tol_mean = seq(1,5.5,0.1)), interval = "confidence"))
    y_angio <- fitted_angio$fit
    ylow_angio <- fitted_angio$lwr
    yhigh_angio <- fitted_angio$upr
    
    fitted_gymno <- as.data.frame(predict(model_gymno, newdata = data.frame(shade_tol_mean = seq(1,5.5,0.1)), interval = "confidence"))
    y_gymno <- fitted_gymno$fit
    ylow_gymno <- fitted_gymno$lwr
    yhigh_gymno <- fitted_gymno$upr
    
    y_points_angio <- as.data.frame(data_y_angio[[i]])
    y_points_gymno <- as.data.frame(data_y_gymno[[i]])
    
    sub_x <- seq(1, 5.5, 0.1)
    sub_y <- rep(0, length(sub_x))
    
    
    
    if (summary(model_angio)$coefficients["shade_tol_mean", "Pr(>|t|)"] <= 0.05 & summary(model_gymno)$coefficients["shade_tol_mean", "Pr(>|t|)"] <= 0.05) {
      
      plot(y_angio ~ x, lwd = 2.5, type = "l", col = "white", las = 1, xlab = "shade tolerance", ylab = models[[i]], xlim = c(1,5.5), ylim = c(-0.03,0.03))
      lines(sub_x, sub_y, col = "gray88", lty = 8, lwd = 3.5) 
      
      polygon(c(x,rev(x)),c(yhigh_angio, rev(ylow_angio)), border = NA, col = rgb(r = 1, g = 0.4, b = 0, max = 1, alpha = 0.5))
      polygon(c(x,rev(x)),c(yhigh_gymno, rev(ylow_gymno)), border = NA, col = rgb(r = 0, g = 0.6, b = 0.2, max = 1, alpha = 0.5))
      
      lines(y_angio ~ x, lwd = 4, col = "#c05a00")
      lines(y_gymno ~ x, lwd = 4, col = "forestgreen")
      
      points(data_sh_angio[,1], y_points_angio[,1], pch = 16, col = "#c05a00", cex = 2.5)
      points(data_sh_gymno[,1], y_points_gymno[,1], pch = 15, col = "forestgreen", cex = 2.5)
      
      
    } else if (summary(model_angio)$coefficients["shade_tol_mean", "Pr(>|t|)"] > 0.05 & summary(model_gymno)$coefficients["shade_tol_mean", "Pr(>|t|)"] > 0.05) {
      
      plot(y_angio ~ x, lwd = 2.5, type = "l", col = "white", las = 1, xlab = "shade tolerance", ylab = models[[i]], xlim = c(1,5.5), ylim = c(-0.03,0.03))
      lines(sub_x, sub_y, col = "gray88", lty = 8, lwd = 2) 
      
      y_points_angio <- as.data.frame(data_y_angio[[i]])
      y_points_gymno <- as.data.frame(data_y_gymno[[i]])
      
      points(data_sh_angio[,1], y_points_angio[,1], pch = 16, col = "#c05a00", cex = 1.5)
      points(data_sh_gymno[,1], y_points_gymno[,1], pch = 15, col = "forestgreen", cex = 1.5)
      
      
    } else if (summary(model_angio)$coefficients["shade_tol_mean", "Pr(>|t|)"] <= 0.05 & summary(model_gymno)$coefficients["shade_tol_mean", "Pr(>|t|)"] > 0.05) {
      
      plot(y_angio ~ x, lwd = 2.5, type = "l", col = "white", las = 1, xlab = "shade tolerance", ylab = models[[i]], xlim = c(1,5.5), ylim = c(-0.03,0.03))
      lines(sub_x, sub_y, col = "gray88", lty = 8, lwd = 2) 
      
      polygon(c(x,rev(x)),c(yhigh_angio, rev(ylow_angio)), border = NA, col = rgb(r = 1, g = 0.4, b = 0, max = 1, alpha = 0.5))
      
      lines(y_angio ~ x, lwd = 4, col = "#c05a00")
      
      points(data_sh_angio[,1], y_points_angio[,1], pch = 16, col = "#c05a00", cex = 2.5)
      points(data_sh_gymno[,1], y_points_gymno[,1], pch = 15, col = "forestgreen", cex = 1.5)
      
      
    }  else {
      
      plot(y_angio ~ x, lwd = 2.5, type = "l", col = "white", las = 1, xlab = "shade tolerance", ylab = models[[i]], xlim = c(1,5.5), ylim = c(-0.03,0.03))
      lines(sub_x, sub_y, col = "gray88", lty = 8, lwd = 2) 
      
      polygon(c(x,rev(x)),c(yhigh_gymno, rev(ylow_gymno)), border = NA, col = rgb(r = 0, g = 0.6, b = 0.2, max = 1, alpha = 0.5))
      
      lines(y_gymno ~ x, lwd = 4, col = "forestgreen")
      
      points(data_sh_angio[,1], y_points_angio[,1], pch = 16, col = "#c05a00", cex = 1.5)
      points(data_sh_gymno[,1], y_points_gymno[,1], pch = 15, col = "forestgreen", cex = 2.5)
      
    }
    
    dev.off
    if (!is.null(dev.list())) dev.off()
    
  }
  
  angio_response_to_competition[, "group"] <- "angio"
  gymno_response_to_competition[, "group"] <- "gymno"
  
  response_to_competition <- rbind(angio_response_to_competition, gymno_response_to_competition)  
  write.csv(file = "output/response_to_competition_unscaled.csv", response_to_competition)
  
  return(response_to_competition)
  
}



samsara_isolated_tree <- function(species_list_comp) {
  
  samsara <- read.csv(file = "data/out_light_alone.csv", sep = ";")

  ## Selecting only data useful for further analyses
  # DBH fixed to 15 cm
  # Local competition fixed to null
  samsara_clean <- samsara %>% dplyr::filter(dbh_cm == 15, 
                                             bal_m2ha == 0)

  
  ## Treating species by taxonomic group
  # Scaling and centering crown characteristics before computing mean values
  angiosperms <- samsara_clean %>% dplyr::filter(order == "A" & species %in% species_list_comp) %>%
    
                                   dplyr::mutate(height_15_s = scale(h_m, center = TRUE, scale = TRUE)[, 1],
                                                 diameter_15_s = scale(cradius_m, center = TRUE, scale = TRUE)[, 1],
                                                 ratio_15_s = scale(cratio, center = TRUE, scale = TRUE)[, 1],
                                                 epot_s = scale(epot, center = TRUE, scale = TRUE)[, 1]) %>%
    
                                   dplyr::select(-h_m, -cradius_m, 
                                                 -cratio, -crown_lad, 
                                                 -epot) %>%
    
                                   dplyr::group_by(species) %>%
                                   dplyr::summarise(height = mean(height_15_s), cv_h = abs(sd(height_15_s)/mean(height_15_s)),
                                                    ratio = mean(ratio_15_s), cv_r = abs(sd(ratio_15_s)/mean(ratio_15_s)),
                                                    diam = mean(diameter_15_s), cv_d = abs(sd(diameter_15_s)/mean(diameter_15_s)),
                                                    epot = mean(epot_s), cv_e = abs(sd(epot_s)/mean(epot_s))) %>%
                                                      
                                   dplyr::ungroup()
  
  
  gymnosperms <- samsara_clean %>% dplyr::filter(order == "G" & species %in% species_list_comp) %>%
    
                                   dplyr::mutate(height_15_s = scale(h_m, center = TRUE, scale = TRUE)[, 1],
                                                 diameter_15_s = scale(cradius_m, center = TRUE, scale = TRUE)[, 1],
                                                 ratio_15_s = scale(cratio, center = TRUE, scale = TRUE)[, 1],
                                                 epot_s = scale(epot, center = TRUE, scale = TRUE)[, 1]) %>%
                                  
                                   dplyr::select(-h_m, -cradius_m, 
                                                 -cratio, -crown_lad, 
                                                 -epot) %>%
                                  
                                   dplyr::group_by(species) %>%
                                   dplyr::summarise(height = mean(height_15_s), cv_h = abs(sd(height_15_s)/mean(height_15_s)),
                                                    ratio = mean(ratio_15_s), cv_r = abs(sd(ratio_15_s)/mean(ratio_15_s)),
                                                    diam = mean(diameter_15_s), cv_d = abs(sd(diameter_15_s)/mean(diameter_15_s)),
                                                    epot = mean(epot_s), cv_e = abs(sd(epot_s)/mean(epot_s))) %>%
                                  
                                   dplyr::ungroup()
  
  
  
  ### Running linear models to explore the relationships between the light potentially intercepted by an isolated tree and crown characteristics 
  ## Performing stepAIC selection using the dredge function
  ## Running models independently for angiosperms and gymnosperms
  
  options(na.action = "na.fail")
  
  # angiosperms
  model_angio <- lm(epot ~ height + ratio + diam, weights = (1/cv_e), data = angiosperms) 
  angio_all_models <- dredge(model_angio) 
  best_angio_model <- get.models(angio_all_models, 1)[[1]]
  summary_angio_best_model <- as.data.frame(angio_all_models)
  write.csv(file = "output/modelselection_isolatedtree_angiosperms.csv", summary_angio_best_model)
  
  
  # gymnosperms
  model_gymno <- lm(epot ~ height + ratio + diam, weights = (1/cv_e), data = gymnosperms) 
  gymno_all_models <- dredge(model_gymno) 
  best_gymno_model <- get.models(gymno_all_models, 1)[[1]]
  summary_gymno_best_model <- as.data.frame(gymno_all_models)
  write.csv(file = "output/modelselection_isolatedtree_gymnosperms.csv", summary_gymno_best_model)
  
  
  ## Running selected models
  angio <- lm(epot ~ height + ratio + diam, weights = (1/cv_e), data = angiosperms)
  gymno <- lm(epot ~ height + ratio + diam, weights = (1/cv_e), data = gymnosperms)

  ## Storing results
  mod_angio <- as.data.frame(confint(angio))
  mod_angio$mean <- summary(angio)$coefficients[, "Estimate"]
  mod_angio$p.value <- summary(angio)$coefficients[, "Pr(>|t|)"]
  mod_angio <- mod_angio[-1,] # remove the Intercept line of results
  colnames(mod_angio) <- c("est.inf", "est.sup", "est", "p.value")
  mod_angio$mod <- c("height", "ratio", "diam")
  
  
  mod_gymno <- as.data.frame(confint(gymno))
  mod_gymno$mean <- summary(gymno)$coefficients[, "Estimate"]
  mod_gymno$p.value <- summary(gymno)$coefficients[, "Pr(>|t|)"]
  mod_gymno <- mod_gymno[-1,] # remove the Intercept line of results
  colnames(mod_gymno) <- c("est.inf", "est.sup", "est", "p.value")
  mod_gymno$mod <- c("height", "ratio", "diam")

  
  ## Transform p.values into levels of significance
  final_angio_results <- mod_angio %>% dplyr::mutate(sign = case_when(p.value <= 0.001 ~ "***",
                                                                      p.value <= 0.01 ~ "**",
                                                                      p.value <= 0.05 ~ "*",
                                                                      TRUE ~ ""))
  
  final_gymno_results <- mod_gymno %>% dplyr::mutate(sign = case_when(p.value <= 0.001 ~ "***",
                                                                      p.value <= 0.01 ~ "**",
                                                                      p.value <= 0.05 ~ "*",
                                                                      TRUE ~ ""))

  
  ## Ploting results
  plot_angio <- final_angio_results %>% dplyr::mutate(est.inf = as.numeric(est.inf),
                                                      est.sup = as.numeric(est.sup),
                                                      est = as.numeric(est)) %>%
    
                                        dplyr::mutate(label.pos = max(est.sup, na.rm = TRUE)) %>%
                                        dplyr::mutate(label.pos = label.pos + 0.5*(max(est.sup, na.rm = TRUE) - min(est.inf, na.rm = TRUE))) %>%
                                        dplyr::mutate(mod = factor(mod, levels = c("diam", "ratio", "height"))) %>%
    
                                        ggplot(aes(x = mod, y = est)) + ylim(0, 1.5) +
                                        geom_errorbar(aes(ymin = est.inf, ymax = est.sup), colour = "black", width = 0, size = 0.5) + 
                                        geom_point(size = 2.5, fill = "chocolate2", shape = 21) +
                                        geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed", size = 0.6) +  
                                        geom_text(aes(y = 0, label = sign), nudge_y = 1.2, size = 8, color = "chocolate2") +
    
                                        xlab("") + ylab("") +
                                        theme(rect = element_rect(fill = "transparent"), 
                                              panel.background = element_rect(fill = "transparent"), 
                                              panel.grid = element_blank(), 
                                              plot.background = element_rect(fill = "transparent", color = NA),
                                              strip.background = element_blank(), 
                                              strip.text.x = element_text(size = 1),
                                              strip.text.y = element_text(size = 1, angle = 360),
                                              axis.title = element_text(size = 1),
                                              legend.position = "none",
                                              axis.text.x = element_text(size = 1, angle = 360),
                                              axis.text.y = element_text(size = 1, angle = 360)) + 
                                        coord_flip() 
  
  
  plot_gymno <- final_gymno_results %>% dplyr::mutate(est.inf = as.numeric(est.inf),
                                                      est.sup = as.numeric(est.sup),
                                                      est = as.numeric(est)) %>%
    
                                        dplyr::mutate(label.pos = max(est.sup, na.rm = TRUE)) %>%
                                        dplyr::mutate(label.pos = label.pos + 0.5*(max(est.sup, na.rm = TRUE) - min(est.inf, na.rm = TRUE))) %>%
                                        dplyr::mutate(mod = factor(mod, levels = c("diam", "ratio", "height"))) %>%
                                        
                                        ggplot(aes(x = mod, y = est)) + ylim(0, 1.5) +
                                        geom_errorbar(aes(ymin = est.inf, ymax = est.sup), colour = "black", width = 0, size = 0.5) + 
                                        geom_point(size = 2.5, fill = "forestgreen", shape = 22) +
                                        geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed", size = 0.6) +  
                                        geom_text(aes(y = 0, label = sign), nudge_y = 1.2, size = 8, color = "forestgreen") +
                                        
                                        xlab("") + ylab("") +
                                        theme(rect = element_rect(fill = "transparent"), 
                                              panel.background = element_rect(fill = "transparent"), 
                                              panel.grid = element_blank(), 
                                              plot.background = element_rect(fill = "transparent", color = NA),
                                              strip.background = element_blank(), 
                                              strip.text.x = element_text(size = 1),
                                              strip.text.y = element_text(size = 1, angle = 360),
                                              axis.title = element_text(size = 1),
                                              legend.position = "none",
                                              axis.text.x = element_text(size = 1, angle = 360),
                                              axis.text.y = element_text(size = 1, angle = 360)) + 
                                        coord_flip() 
  
  ggsave("manuscript/angiosperms_isolatedtree_results.png", plot_angio, width = 5, height = 6, units = "cm", bg = "transparent") 
  ggsave("manuscript/gymnosperms_isolatedtree_results.png", plot_gymno, width = 5, height = 6, units = "cm", bg = "transparent") 
  
  final_angio_results$group <- "angio"
  final_gymno_results$group <- "gymno"                                                    
  final_isolatedtree_results <- rbind(final_angio_results, final_gymno_results)  
  
  write.csv(file = "output/isolated_tree_results.csv", final_isolatedtree_results)

  return(final_isolatedtree_results)  
  
}


samsara_codominance_suppressed_trees <- function(species_list_comp) {
  
  ## Loading data files
  samsara <- read.csv(file = "data/out_light_compet.csv", sep = ";")
  no_comp <- read.csv(file = "data/out_light_alone.csv", sep = ";")
  
  ## Creating columns to select data according to competition contexts more easily
  samsara <- samsara %>%  tidyr::separate_wider_delim(id_pattern, "_",
                                                      names = c("dbh", "comp", "rep"))
  
  ## Selecting data in the No competition context
  no_comp <- no_comp %>% dplyr::filter(dbh_cm == 15, 
                                       bal_m2ha == 0)
  
  
  ## Creating file with responses of crown characteristics to different levels of competition
  res <- left_join(samsara, no_comp, 
                   by = c("id_allom", "species", "order"),
                   suffix = c("_comp", "_nocomp")) %>% 
         dplyr::mutate(lci = e_comp/e_nocomp,
                       rep_h = h_m_comp - h_m_nocomp,
                       rep_cd = cradius_m_comp - cradius_m_nocomp,
                       rep_cr = cratio_comp - cratio_nocomp,
                       rep_cd_r = cradius_m_comp/cradius_m_nocomp)%>%
         dplyr::filter((dbh == "15"& comp == "low") | (dbh == "30"& comp == "high")) %>%
         dplyr::mutate(comp2 = paste(dbh, comp, sep = "_")) %>%
         dplyr::filter(species %in% species_list_comp)
  
  
  ## Removing observations for which lci values were abnormally high 
  # arguing that an increase of more than a factor 3 of crown diameter is unrealistic 
  # this selection removes less than 1% of the data
  res <- res %>% dplyr::filter(rep_cd_r < 3) 
  
  
  ## Grouping species by taxonomic group and computing mean values of the responses obtained in the different competitive contexts
  res_sp <- res %>% dplyr::group_by(order) %>% 
                    dplyr::mutate(across(c(lci, h_m_nocomp,rep_h, rep_cd, rep_cr), scale)) %>% 
                    dplyr::group_by(species, order, comp2) %>% dplyr::summarise(cv_lci = abs(sd(lci)/mean(lci)),
                                                         lci = mean(lci),
                                                         h_nocomp = mean(h_m_nocomp),
                                                         rep_height = mean(rep_h),
                                                         rep_diam = mean(rep_cd),
                                                         rep_ratio = mean(rep_cr)) %>%
                    dplyr::ungroup()
  
  
  ## Running linear models to explore the relationships between the light actually intercepted by a tree under codominance and crown characteristics' responses to competition
  # Performing stepAIC selection using the dredge function
  # Running models independently for angiosperms and gymnosperms
  
  options(na.action = "na.fail")
  
  # angiosperms
  df_angio_codo <- na.omit(res_sp[res_sp$order =="A" & res_sp$comp2 == "15_low",])
  model_angio_codo <- lm(lci ~ h_nocomp + rep_height + rep_ratio + rep_diam, weights = (1/cv_lci), data = df_angio_codo) 
  angio_all_models_codo <- dredge(model_angio_codo) 
  best_angio_model_codo <- get.models(angio_all_models_codo, 1)[[1]]
  summary_angio_models_codo <- as.data.frame(angio_all_models_codo)
  summary_angio_models_codo$cond <- "codominance"
  summary_angio_models_codo$order <- "angio"
  
  
  df_angio_supp <- na.omit(res_sp[res_sp$order =="A" & res_sp$comp2 == "30_high",])
  model_angio_supp <- lm(lci ~ h_nocomp + rep_height + rep_ratio + rep_diam, weights = (1/cv_lci), data = df_angio_supp) 
  angio_all_models_supp <- dredge(model_angio_supp) 
  best_angio_model_supp <- get.models(angio_all_models_supp, 2)[[1]]
  summary_angio_models_supp <- as.data.frame(angio_all_models_supp)
  summary_angio_models_supp$cond <- "suppressed"
  summary_angio_models_supp$order <- "angio"
  
  
  # gymnosperms
  df_gymno_codo <- na.omit(res_sp[res_sp$order =="G" & res_sp$comp2 == "15_low",])
  model_gymno_codo <- lm(lci ~ h_nocomp + rep_height + rep_ratio + rep_diam, weights = (1/cv_lci), data = df_gymno_codo) 
  gymno_all_models_codo <- dredge(model_gymno_codo) 
  best_gymno_model_codo <- get.models(gymno_all_models_codo, 1)[[1]]
  summary_gymno_models_codo <- as.data.frame(gymno_all_models_codo)
  summary_gymno_models_codo$cond <- "codominance"
  summary_gymno_models_codo$order <- "gymno"
  
  
  df_gymno_supp <- na.omit(res_sp[res_sp$order =="G" & res_sp$comp2 == "30_high",])
  model_gymno_supp <- lm(lci ~ h_nocomp + rep_height + rep_ratio + rep_diam, weights = (1/cv_lci), data = df_gymno_supp) 
  gymno_all_models_supp <- dredge(model_gymno_supp) 
  best_gymno_model_supp <- get.models(gymno_all_models_supp, 1)[[1]]
  summary_gymno_models_supp <- as.data.frame(gymno_all_models_supp)
  summary_gymno_models_supp$cond <- "suppressed"
  summary_gymno_models_supp$order <- "gymno"
  
  results <- rbind(summary_angio_models_codo,
                   summary_angio_models_supp,
                   summary_gymno_models_codo,
                   summary_gymno_models_supp)
  
  write.csv(file = "output/modelselection_competition_results.csv", results)
  
  
  ## Storing models' results
  # angiosperms
  mod_angio_codo <- as.data.frame(confint(best_angio_model_codo))
  mod_angio_codo$mean <- summary(best_angio_model_codo)$coefficients[, "Estimate"]
  mod_angio_codo$p.value <- summary(best_angio_model_codo)$coefficients[, "Pr(>|t|)"]
  mod_angio_codo <- mod_angio_codo[-1,] # remove the Intercept line of results
  colnames(mod_angio_codo) <- c("est.inf", "est.sup", "est", "p.value")
  mod_angio_codo$mod <- rownames(mod_angio_codo)
  
  
  mod_angio_supp <- as.data.frame(confint(best_angio_model_supp))
  mod_angio_supp$mean <- summary(best_angio_model_supp)$coefficients[, "Estimate"]
  mod_angio_supp$p.value <- summary(best_angio_model_supp)$coefficients[, "Pr(>|t|)"]
  mod_angio_supp <- mod_angio_supp[-1,] # remove the Intercept line of results
  colnames(mod_angio_supp) <- c("est.inf", "est.sup", "est", "p.value")
  mod_angio_supp$mod <- rownames(mod_angio_supp)
  
  
  
  # gymnosperms
  mod_gymno_codo <- as.data.frame(confint(best_gymno_model_codo))
  mod_gymno_codo$mean <- summary(best_gymno_model_codo)$coefficients[, "Estimate"]
  mod_gymno_codo$p.value <- summary(best_gymno_model_codo)$coefficients[, "Pr(>|t|)"]
  mod_gymno_codo <- mod_gymno_codo[-1,] # remove the Intercept line of results
  colnames(mod_gymno_codo) <- c("est.inf", "est.sup", "est", "p.value")
  mod_gymno_codo$mod <- rownames(mod_gymno_codo)
  
  
  mod_gymno_supp <- as.data.frame(confint(best_gymno_model_supp))
  mod_gymno_supp$mean <- summary(best_gymno_model_supp)$coefficients[, "Estimate"]
  mod_gymno_supp$p.value <- summary(best_gymno_model_supp)$coefficients[, "Pr(>|t|)"]
  mod_gymno_supp <- mod_gymno_supp[-1,] # remove the Intercept line of results
  colnames(mod_gymno_supp) <- c("est.inf", "est.sup", "est", "p.value")
  mod_gymno_supp$mod <- rownames(mod_gymno_supp)
  
  
  ## Making sure to have similar files before plotting the results
  # creating list of results and of explanatory variables
  all_results <- list(best_angio_model_codo, best_angio_model_supp,
                      best_gymno_model_codo, best_gymno_model_supp)
  
  variables <- c("h_nocomp", "rep_height", "rep_ratio", "rep_diam")
  
  results <- list()
  for (i in 1:4) {
    results[[i]] <- matrix(NA, nrow = length(variables), ncol = 6) 
    colnames(results[[i]]) <- c("mod", "est.inf", "est.sup", "est", "p.value", "sign")
  }
   
  names(results) <- c("angio_codo", "angio_supp", 
                      "gymno_codo", "gymno_supp")
  
  
  for (i in 1:length(results)) {
    
    for (j in 1:length(variables)) {
      
      results[[i]][j,"mod"] <- variables[j]
      
      tryCatch({
        
        model <- all_results[[i]]
        
        results[[i]][j, "est.inf"] = confint.lm(model)[variables[j], "2.5 %"]
        results[[i]][j, "est.sup"] = confint.lm(model)[variables[j], "97.5 %"]
        results[[i]][j, "est"] = summary(model)$coefficients[variables[j], "Estimate"]
        results[[i]][j, "p.value"] = summary(model)$coefficients[variables[j], "Pr(>|t|)"]
        
      },
      error = function(e) {
        print("not included in the selected model")
      })
      
    }
    
  }
  
  
  # Creating a function to convert a matrix to a data frame, mutate, and convert back to a matrix
  mutate_matrix <- function(mat) {
    mat_df <- as.data.frame(mat)  # convert to data frame
    mat_df <- mat_df %>%
      mutate(
        sign = case_when(
          p.value <= 0.001 ~ "***",
          p.value <= 0.01 ~ "**",
          p.value <= 0.05 ~ "*",
          TRUE ~ ""
        )
      )
    as.matrix(mat_df)  # convert back to a matrix
  }
  
  # Applying the mutation to all matrices in the list of results
  final_results <- lapply(results, mutate_matrix)

 
  ## Ploting results
  # creating list of colors, shapes and names

  c <- c("chocolate2", "chocolate2",
         "forestgreen", "forestgreen")
  
  s <- c(21, 21, 22, 22)
  
  names <- c("angio_codo", "angio_supp", "gymno_codo", "gymno_supp")
  
  for (i in c(1,3)) {
    
    df <- as.data.frame(final_results[[i]])
  
    plot <- df %>% dplyr::mutate(est.inf = as.numeric(est.inf),
                                 est.sup = as.numeric(est.sup),
                                 est = as.numeric(est)) %>%
    
                   dplyr::mutate(label.pos = max(est.sup, na.rm = TRUE)) %>%
                   dplyr::mutate(label.pos = label.pos + 0.5*(max(est.sup, na.rm = TRUE) - min(est.inf, na.rm = TRUE))) %>%
                   dplyr::mutate(mod = factor(mod, levels = c("rep_diam", "rep_ratio", "rep_height", "h_nocomp"))) %>%
    
                   ggplot(aes(x = mod, y = est)) + ylim(0, 4.5) +
                   geom_errorbar(aes(ymin = est.inf, ymax = est.sup), colour = "black", width = 0, size = 0.5) + 
                   geom_point(size = 2.5, fill = c[i], shape = s[i]) +
                   geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed", size = 0.6) +  
                   geom_text(aes(y = 0, label = sign), nudge_y = 4, size = 8, color = c[i]) +
    
                   xlab("") + ylab("") +
                   theme(rect = element_rect(fill = "transparent"), 
                         panel.background = element_rect(fill = "transparent"), 
                         panel.grid = element_blank(), 
                         plot.background = element_rect(fill = "transparent", color = NA),
                         strip.background = element_blank(), 
                         strip.text.x = element_text(size = 1),
                         strip.text.y = element_text(size = 1, angle = 360),
                         axis.title = element_text(size = 1),
                         legend.position = "none",
                         axis.text.x = element_text(size = 1, angle = 360),
                         axis.text.y = element_text(size = 1, angle = 360)) + 
                   coord_flip() 
    
    ggsave(paste0("manuscript/", names[i], "_competition.png"), plot, width = 4, height = 6, units = "cm", bg = "transparent") 
    
  }
  
  
  
  for (i in c(2,4)) {
    
    df <- as.data.frame(final_results[[i]])
    
    plot <- df %>% dplyr::mutate(est.inf = as.numeric(est.inf),
                                 est.sup = as.numeric(est.sup),
                                 est = as.numeric(est)) %>%
      
                   dplyr::mutate(label.pos = max(est.sup, na.rm = TRUE)) %>%
                   dplyr::mutate(label.pos = label.pos + 0.5*(max(est.sup, na.rm = TRUE) - min(est.inf, na.rm = TRUE))) %>%
                   dplyr::mutate(mod = factor(mod, levels = c("rep_diam", "rep_ratio", "rep_height", "h_nocomp"))) %>%
                  
                   ggplot(aes(x = mod, y = est)) + ylim(0, 0.4) +
                   geom_errorbar(aes(ymin = est.inf, ymax = est.sup), colour = "black", width = 0, size = 0.5) + 
                   geom_point(size = 2.5, fill = c[i], shape = s[i]) +
                   geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed", size = 0.6) +  
                   geom_text(aes(y = 0, label = sign), nudge_y = 0.32, size = 8, color = c[i]) +
                  
                   xlab("") + ylab("") +
                   theme(rect = element_rect(fill = "transparent"), 
                         panel.background = element_rect(fill = "transparent"), 
                         panel.grid = element_blank(), 
                         plot.background = element_rect(fill = "transparent", color = NA),
                         strip.background = element_blank(), 
                         strip.text.x = element_text(size = 1),
                         strip.text.y = element_text(size = 1, angle = 360),
                         axis.title = element_text(size = 1),
                         legend.position = "none",
                         axis.text.x = element_text(size = 1, angle = 360),
                         axis.text.y = element_text(size = 1, angle = 360)) + 
                   coord_flip() 
    
    ggsave(paste0("manuscript/", names[i], "_competition.png"), plot, width = 4, height = 6, units = "cm", bg = "transparent") 
    
  }
  
  return(final_results)
  
}




samsara_codominated_tree <- function(species_list_comp) {
  
  ## Loading data files
  samsara <- read.csv(file = "data/out_light_compet.csv", sep = ";")
  no_comp <- read.csv(file = "data/out_light_alone.csv", sep = ";")

  ## Selecting useful data and renaming columns
  samsara_nocomp <- no_comp %>% dplyr::filter(dbh_cm == 15, 
                                              bal_m2ha == 0) %>%
                                dplyr::rename(h_nocomp = h_m,
                                              ratio_nocomp = cratio,
                                              diam_nocomp = cradius_m) %>%
                                dplyr::select(h_nocomp, ratio_nocomp, diam_nocomp, epot, order)

  
  samsara_codominance <- samsara %>% dplyr::mutate(parts = str_split_fixed(id_pattern, "_", 3)) %>%
                                     dplyr::mutate(dbh = parts[, 1], comp = parts[, 2], rep = parts[, 3]) %>%
                                     dplyr::select(-parts) %>%
                                     dplyr::filter(dbh == "15", comp == "low") %>%
                                     dplyr::select(-id_pattern, -bal_m2ha, -crown_lad, -epot, -dbh, -comp) %>%
                                     dplyr::group_by(species, id_allom) %>%
                                     dplyr::summarise(h_comp = mean(h_m),
                                                      ratio_comp = mean(cratio),
                                                      diam_comp = mean(cradius_m),
                                                      eint_comp = mean(e)) %>%
                                     dplyr::ungroup() %>%
                                     dplyr::select(-id_allom)  
  
  
  ## Having no competition and competition data into a single file
  samsara_complete <- cbind(samsara_codominance, samsara_nocomp) 
  samsara_complete <- samsara_complete[samsara_complete$epot > 0,]
  
  ## Calculating responses of crown characteristics to competition
  samsara_complete <- samsara_complete %>% dplyr::mutate(e = eint_comp / epot,
                                                         rep_h = h_comp - h_nocomp,
                                                         rep_r = ratio_comp - ratio_nocomp,
                                                         rep_d = diam_comp - diam_nocomp) %>%
                                           dplyr::filter(species %in% species_list_comp)
  
  ## Treating species by taxonomic group
  # Scaling and centering crown characteristics before computing mean values
  angiosperms <- samsara_complete %>% dplyr::filter(order == "A") %>%
    
                                      dplyr::mutate(height_15_s = scale(h_nocomp, center = TRUE, scale = TRUE)[, 1],
                                                    rep_height_s = scale(rep_h, center = TRUE, scale = TRUE)[, 1],
                                                    rep_ratio_15_s = scale(rep_r, center = TRUE, scale = TRUE)[, 1],
                                                    rep_diam_15_s = scale(rep_d, center = TRUE, scale = TRUE)[, 1],
                                                    e_s = scale(e, center = TRUE, scale = TRUE)[, 1]) %>%
  
                                      dplyr::group_by(species) %>%
                                      dplyr::summarise(h_no_comp = mean(height_15_s), cv_h_no_comp = abs(sd(height_15_s)/mean(height_15_s)),
                                                       height = mean(rep_height_s), cv_h = abs(sd(rep_height_s)/mean(rep_height_s)),
                                                       ratio = mean(rep_ratio_15_s), cv_r = abs(sd(rep_ratio_15_s)/mean(rep_ratio_15_s)),
                                                       diam = mean(rep_diam_15_s), cv_d = abs(sd(rep_diam_15_s)/mean(rep_diam_15_s)),
                                                       lci = mean(e_s), cv_e = abs(sd(e_s)/mean(e_s))) %>%
  
                                      dplyr::ungroup()
  
  gymnosperms <- samsara_complete %>% dplyr::filter(order == "G") %>%
    
                                      dplyr::mutate(height_15_s = scale(h_nocomp, center = TRUE, scale = TRUE)[, 1],
                                                    rep_height_s = scale(rep_h, center = TRUE, scale = TRUE)[, 1],
                                                    rep_ratio_15_s = scale(rep_r, center = TRUE, scale = TRUE)[, 1],
                                                    rep_diam_15_s = scale(rep_d, center = TRUE, scale = TRUE)[, 1],
                                                    e_s = scale(e, center = TRUE, scale = TRUE)[, 1]) %>%
                                      
                                      dplyr::group_by(species) %>%
                                      dplyr::summarise(h_no_comp = mean(height_15_s), cv_h_no_comp = abs(sd(height_15_s)/mean(height_15_s)),
                                                       height = mean(rep_height_s), cv_h = abs(sd(rep_height_s)/mean(rep_height_s)),
                                                       ratio = mean(rep_ratio_15_s), cv_r = abs(sd(rep_ratio_15_s)/mean(rep_ratio_15_s)),
                                                       diam = mean(rep_diam_15_s), cv_d = abs(sd(rep_diam_15_s)/mean(rep_diam_15_s)),
                                                       lci = mean(e_s), cv_e = abs(sd(e_s)/mean(e_s))) %>%
                                      
                                      dplyr::ungroup()
  
  
  ### Running linear models to explore the relationships between the light actually intercepted by a tree under codominance and crown characteristics' responses to competition
  ## Performing stepAIC selection using the dredge function
  ## Running models independently for angiosperms and gymnosperms
  
  options(na.action = "na.fail")
  
  # angiosperms
  model_angio <- lm(lci ~ h_no_comp + height + ratio + diam, weights = (1/cv_e), data = angiosperms) 
  angio_all_models <- dredge(model_angio) 
  best_angio_model <- get.models(angio_all_models, 1)[[1]]
  summary_angio_best_model <- as.data.frame(angio_all_models)
  write.csv(file = "output/modelselection_codominance_angiosperms.csv", summary_angio_best_model)
  
  
  # gymnosperms
  model_gymno <- lm(lci ~ h_no_comp + height + ratio + diam, weights = (1/cv_e), data = gymnosperms) 
  gymno_all_models <- dredge(model_gymno) 
  best_gymno_model <- get.models(gymno_all_models, 1)[[1]]
  summary_gymno_best_model <- as.data.frame(gymno_all_models)
  write.csv(file = "output/modelselection_codominance_gymnosperms.csv", summary_gymno_best_model)
  
  
  ## Running selected models
  angio <- lm(lci ~ h_no_comp, weights = (1/cv_e), data = angiosperms)
  gymno <- lm(lci ~ h_no_comp + height + diam, weights = (1/cv_e), data = gymnosperms)
  
  ## Storing results
  mod_angio <- as.data.frame(confint(angio))
  mod_angio$mean <- summary(angio)$coefficients[, "Estimate"]
  mod_angio$p.value <- summary(angio)$coefficients[, "Pr(>|t|)"]
  mod_angio <- mod_angio[-1,] # remove the Intercept line of results
  colnames(mod_angio) <- c("est.inf", "est.sup", "est", "p.value")
  mod_angio$mod <- c("height_nocomp")
  mod_angio <- rbind(mod_angio, c(NA, NA, NA, NA, NA), c(NA, NA, NA, NA, NA))
  
  
  mod_gymno <- as.data.frame(confint(gymno))
  mod_gymno$mean <- summary(gymno)$coefficients[, "Estimate"]
  mod_gymno$p.value <- summary(gymno)$coefficients[, "Pr(>|t|)"]
  mod_gymno <- mod_gymno[-1,] # remove the Intercept line of results
  colnames(mod_gymno) <- c("est.inf", "est.sup", "est", "p.value")
  mod_gymno$mod <- c("height_nocomp", "height", "diam")
  
  
  ## Transform p.values into levels of significance
  final_angio_results <- mod_angio %>% dplyr::mutate(sign = case_when(p.value <= 0.001 ~ "***",
                                                                      p.value <= 0.01 ~ "**",
                                                                      p.value <= 0.05 ~ "*",
                                                                      TRUE ~ ""))
  
  final_gymno_results <- mod_gymno %>% dplyr::mutate(sign = case_when(p.value <= 0.001 ~ "***",
                                                                      p.value <= 0.01 ~ "**",
                                                                      p.value <= 0.05 ~ "*",
                                                                      TRUE ~ ""))
  
  
  ## Ploting results
  angio <- final_angio_results %>% dplyr::mutate(est.inf = as.numeric(est.inf),
                                                 est.sup = as.numeric(est.sup),
                                                 est = as.numeric(est)) %>%
    
                                   dplyr::mutate(label.pos = max(est.sup, na.rm = TRUE)) %>%
                                   dplyr::mutate(label.pos = label.pos + 0.5*(max(est.sup, na.rm = TRUE) - min(est.inf, na.rm = TRUE))) %>%
                                   dplyr::mutate(mod = factor(mod, levels = c("diam", "ratio", "height", "height_nocomp"))) %>%
                      
                                   ggplot(aes(x = mod, y = est)) + ylim(-1.5, 1.5) +
                                   geom_errorbar(aes(ymin = est.inf, ymax = est.sup), colour = "black", width = 0, size = 0.5) + 
                                   geom_point(size = 2.5, fill = "chocolate2", shape = 22) +
                                   geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed", size = 0.6) +  
                                   geom_text(aes(y = 0, label = sign), nudge_y = 1.2, size = 8, color = "chocolate2") +
                      
                                   xlab("") + ylab("") +
                                   theme(rect = element_rect(fill = "transparent"), 
                                         panel.background = element_rect(fill = "transparent"), 
                                         panel.grid = element_blank(), 
                                         plot.background = element_rect(fill = "transparent", color = NA),
                                         strip.background = element_blank(), 
                                         strip.text.x = element_text(size = 1),
                                         strip.text.y = element_text(size = 1, angle = 360),
                                         axis.title = element_text(size = 1),
                                         legend.position = "none",
                                         axis.text.x = element_text(size = 1, angle = 360),
                                         axis.text.y = element_text(size = 1, angle = 360)) + 
                                   coord_flip() 
  
  ggsave("manuscript/angiosperms", plot_angio, width = 5, height = 6, units = "cm", bg = "transparent") 
  

  
  
  plot_gymno <- final_gymno_results %>% dplyr::mutate(est.inf = as.numeric(est.inf),
                                                      est.sup = as.numeric(est.sup),
                                                      est = as.numeric(est)) %>%
    
                                        dplyr::mutate(label.pos = max(est.sup, na.rm = TRUE)) %>%
                                        dplyr::mutate(label.pos = label.pos + 0.5*(max(est.sup, na.rm = TRUE) - min(est.inf, na.rm = TRUE))) %>%
                                        dplyr::mutate(mod = factor(mod, levels = c("diam", "ratio", "height", "height_nocomp"))) %>%
                                        
                                        ggplot(aes(x = mod, y = est)) + ylim(-1.5, 1.5) +
                                        geom_errorbar(aes(ymin = est.inf, ymax = est.sup), colour = "black", width = 0, size = 0.5) + 
                                        geom_point(size = 2.5, fill = "forestgreen", shape = 22) +
                                        geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed", size = 0.6) +  
                                        geom_text(aes(y = 0, label = sign), nudge_y = 1.5, size = 8, color = "forestgreen") +
                                        
                                        xlab("") + ylab("") +
                                        theme(rect = element_rect(fill = "transparent"), 
                                              panel.background = element_rect(fill = "transparent"), 
                                              panel.grid = element_blank(), 
                                              plot.background = element_rect(fill = "transparent", color = NA),
                                              strip.background = element_blank(), 
                                              strip.text.x = element_text(size = 1),
                                              strip.text.y = element_text(size = 1, angle = 360),
                                              axis.title = element_text(size = 1),
                                              legend.position = "none",
                                              axis.text.x = element_text(size = 1, angle = 360),
                                              axis.text.y = element_text(size = 1, angle = 360)) + 
                                        coord_flip() 
  
  ggsave("manuscript/gymnosperms_codominance_results.png", plot_gymno, width = 5, height = 6, units = "cm", bg = "transparent") 
  
  final_angio_results <- final_angio_results %>% dplyr::mutate(group = "angio")
  final_gymno_results$group <- "gymno"                                                    
  final_codominance_results <- rbind(final_angio_results, final_gymno_results)  
  
  write.csv(file = "output/codominant_tree_results.csv", final_codominance_results)
  
  return(final_codominance_results)  
  
}






samsara_suppressed_tree <- function() {
  
  ## Loading data files
  samsara <- read.csv(file = "data/out_light_compet.csv", sep = ";")
  no_comp <- read.csv(file = "data/out_light_alone.csv", sep = ";")
  
  ## Selecting useful data and renaming columns
  samsara_nocomp <- no_comp %>% dplyr::filter(dbh_cm == 15, 
                                              bal_m2ha == 0) %>%
                                dplyr::rename(h_nocomp = h_m,
                                              ratio_nocomp = cratio,
                                              diam_nocomp = cradius_m) %>%
                                dplyr::select(h_nocomp, ratio_nocomp, diam_nocomp, epot, order)
  
  
  samsara_suppressed <- samsara %>% dplyr::mutate(parts = str_split_fixed(id_pattern, "_", 3)) %>%
                                    dplyr::mutate(dbh = parts[, 1], comp = parts[, 2], rep = parts[, 3]) %>%
                                    dplyr::select(-parts) %>%
                                    dplyr::filter(dbh == "15", comp == "high") %>%
                                    dplyr::select(-id_pattern, -bal_m2ha, -crown_lad, -epot, -dbh, -comp) %>%
                                    dplyr::group_by(species, id_allom) %>%
                                    dplyr::summarise(h_comp = mean(h_m),
                                                     ratio_comp = mean(cratio),
                                                     diam_comp = mean(cradius_m),
                                                     eint_comp = mean(e)) %>%
                                    dplyr::ungroup() %>%
                                    dplyr::select(-id_allom)  
  
  
  ## Having no competition and competition data into a single file
  samsara_complete <- cbind(samsara_suppressed, samsara_nocomp) 
  samsara_complete <- samsara_complete[samsara_complete$epot > 0,]
  
  ## Calculating responses of crown characteristics to competition
  samsara_complete <- samsara_complete %>% dplyr::mutate(e = eint_comp / epot,
                                                         rep_h = h_comp - h_nocomp,
                                                         rep_r = ratio_comp - ratio_nocomp,
                                                         rep_d = diam_comp - diam_nocomp)
  
  ## Treating species by taxonomic group
  # Scaling and centering crown characteristics before computing mean values
  angiosperms <- samsara_complete %>% dplyr::filter(order == "A") %>%
    
                                      dplyr::mutate(height_15_s = scale(h_nocomp, center = TRUE, scale = TRUE)[, 1],
                                                    rep_height_s = scale(rep_h, center = TRUE, scale = TRUE)[, 1],
                                                    rep_ratio_15_s = scale(rep_r, center = TRUE, scale = TRUE)[, 1],
                                                    rep_diam_15_s = scale(rep_d, center = TRUE, scale = TRUE)[, 1],
                                                    e_s = scale(e, center = TRUE, scale = TRUE)[, 1]) %>%
                                      
                                      dplyr::group_by(species) %>%
                                      dplyr::summarise(h_no_comp = mean(height_15_s), cv_h_no_comp = abs(sd(height_15_s)/mean(height_15_s)),
                                                       height = mean(rep_height_s), cv_h = abs(sd(rep_height_s)/mean(rep_height_s)),
                                                       ratio = mean(rep_ratio_15_s), cv_r = abs(sd(rep_ratio_15_s)/mean(rep_ratio_15_s)),
                                                       diam = mean(rep_diam_15_s), cv_d = abs(sd(rep_diam_15_s)/mean(rep_diam_15_s)),
                                                       lci = mean(e_s), cv_e = abs(sd(e_s)/mean(e_s))) %>%
                                      
                                      dplyr::ungroup()
  
  gymnosperms <- samsara_complete %>% dplyr::filter(order == "G") %>%
    
                                      dplyr::mutate(height_15_s = scale(h_nocomp, center = TRUE, scale = TRUE)[, 1],
                                                    rep_height_s = scale(rep_h, center = TRUE, scale = TRUE)[, 1],
                                                    rep_ratio_15_s = scale(rep_r, center = TRUE, scale = TRUE)[, 1],
                                                    rep_diam_15_s = scale(rep_d, center = TRUE, scale = TRUE)[, 1],
                                                    e_s = scale(e, center = TRUE, scale = TRUE)[, 1]) %>%
                                      
                                      dplyr::group_by(species) %>%
                                      dplyr::summarise(h_no_comp = mean(height_15_s), cv_h_no_comp = abs(sd(height_15_s)/mean(height_15_s)),
                                                       height = mean(rep_height_s), cv_h = abs(sd(rep_height_s)/mean(rep_height_s)),
                                                       ratio = mean(rep_ratio_15_s), cv_r = abs(sd(rep_ratio_15_s)/mean(rep_ratio_15_s)),
                                                       diam = mean(rep_diam_15_s), cv_d = abs(sd(rep_diam_15_s)/mean(rep_diam_15_s)),
                                                       lci = mean(e_s), cv_e = abs(sd(e_s)/mean(e_s))) %>%
                                      
                                      dplyr::ungroup()
  
  
  ### Running linear models to explore the relationships between the light actually intercepted by a tree under codominance and crown characteristics' responses to competition
  ## Performing stepAIC selection using the dredge function
  ## Running models independently for angiosperms and gymnosperms
  
  options(na.action = "na.fail")
  
  # angiosperms
  model_angio <- lm(lci ~ h_no_comp + height + ratio + diam, weights = (1/cv_e), data = angiosperms) 
  angio_all_models <- dredge(model_angio) 
  best_angio_model <- get.models(angio_all_models, 1)[[1]]
  summary_angio_best_model <- as.data.frame(angio_all_models)
  write.csv(file = "output/modelselection_suppressed_angiosperms.csv", summary_angio_best_model)
  
  
  # gymnosperms
  model_gymno <- lm(lci ~ h_no_comp + height + ratio + diam, weights = (1/cv_e), data = gymnosperms) 
  gymno_all_models <- dredge(model_gymno) 
  best_gymno_model <- get.models(gymno_all_models, 1)[[1]]
  summary_gymno_best_model <- as.data.frame(gymno_all_models)
  write.csv(file = "output/modelselection_suppressed_gymnosperms.csv", summary_gymno_best_model)
  
  
  ## Running selected models
  angio <- lm(lci ~ 1, weights = (1/cv_e), data = angiosperms)
  gymno <- lm(lci ~ h_no_comp + height + ratio + diam, weights = (1/cv_e), data = gymnosperms)
  
  ## Storing results
  mod_angio <- as.data.frame(confint(angio))
  mod_angio$mean <- summary(angio)$coefficients[, "Estimate"]
  mod_angio$p.value <- summary(angio)$coefficients[, "Pr(>|t|)"]
  mod_angio <- mod_angio[-1,] # remove the Intercept line of results
  colnames(mod_angio) <- c("est.inf", "est.sup", "est", "p.value")
  mod_angio <- mod_angio %>% dplyr::mutate(mod = NA)
  
  mod_gymno <- as.data.frame(confint(gymno))
  mod_gymno$mean <- summary(gymno)$coefficients[, "Estimate"]
  mod_gymno$p.value <- summary(gymno)$coefficients[, "Pr(>|t|)"]
  mod_gymno <- mod_gymno[-1,] # remove the Intercept line of results
  colnames(mod_gymno) <- c("est.inf", "est.sup", "est", "p.value")
  mod_gymno$mod <- c("height_nocomp", "height", "ratio", "diam")
  
  
  ## Transform p.values into levels of significance
  final_angio_results <- mod_angio %>% dplyr::mutate(sign = case_when(p.value <= 0.001 ~ "***",
                                                                      p.value <= 0.01 ~ "**",
                                                                      p.value <= 0.05 ~ "*",
                                                                      TRUE ~ ""))
  
  final_gymno_results <- mod_gymno %>% dplyr::mutate(sign = case_when(p.value <= 0.001 ~ "***",
                                                                      p.value <= 0.01 ~ "**",
                                                                      p.value <= 0.05 ~ "*",
                                                                      TRUE ~ ""))
  
  
  ## Ploting results
  # no plot for angiosperms as the null model was retained
  plot_gymno <- final_gymno_results %>% dplyr::mutate(est.inf = as.numeric(est.inf),
                                                      est.sup = as.numeric(est.sup),
                                                      est = as.numeric(est)) %>%
    
                                        dplyr::mutate(label.pos = max(est.sup, na.rm = TRUE)) %>%
                                        dplyr::mutate(label.pos = label.pos + 0.5*(max(est.sup, na.rm = TRUE) - min(est.inf, na.rm = TRUE))) %>%
                                        dplyr::mutate(mod = factor(mod, levels = c("diam", "ratio", "height", "height_nocomp"))) %>%
                                        
                                        ggplot(aes(x = mod, y = est)) + ylim(-1.5, 1.5) +
                                        geom_errorbar(aes(ymin = est.inf, ymax = est.sup), colour = "black", width = 0, size = 0.2) + 
                                        geom_point(size = 1.5, fill = "forestgreen", shape = 22) +
                                        geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed", size = 0.3) +  
                                        geom_text(aes(y = 0, label = sign), nudge_y = 1.5, size = 5, color = "forestgreen") +
                                        
                                        xlab("") + ylab("") +
                                        theme(rect = element_rect(fill = "transparent"), 
                                              panel.background = element_rect(fill = "transparent"), 
                                              panel.grid = element_blank(), 
                                              plot.background = element_rect(fill = "transparent", color = NA),
                                              strip.background = element_blank(), 
                                              strip.text.x = element_text(size = 1),
                                              strip.text.y = element_text(size = 1, angle = 360),
                                              axis.title = element_text(size = 1),
                                              legend.position = "none",
                                              axis.text.x = element_text(size = 1, angle = 360),
                                              axis.text.y = element_text(size = 1, angle = 360)) + 
                                        coord_flip() 
  
  ggsave("manuscript/gymnosperms_suppressed_results.png", plot_gymno, width = 5, height = 6, units = "cm", bg = "transparent") 
  
  final_angio_results <- final_angio_results %>% dplyr::mutate(group = "angio")
  final_gymno_results$group <- "gymno"                                                    
  final_suppressed_results <- rbind(final_angio_results, final_gymno_results)  
  
  write.csv(file = "output/suppressed_tree_results.csv", final_suppressed_results)
  
  return(final_suppressed_results)  
  
}

  
  

  



    
    
    
    
  
    
  