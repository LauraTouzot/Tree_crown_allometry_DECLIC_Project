protocol_parameter <- function(mean_characteristics_nocomp) {
  
  height <- read.csv(file = "data/height_asympt_nocomp_RMSE_AIC_results.csv")
  height <- height %>% dplyr::filter(species %in% unique(mean_characteristics_nocomp$species)) %>%
                       dplyr::summarise(total_weighted = sum(mean_RMSE_weighted),
                                        total_not_weighted = sum(mean_RMSE_no_weighted))
  
  ratio <- read.csv(file = "data/ratio_nocomp_RMSE_AIC_results.csv")
  ratio <- ratio %>% dplyr::filter(species %in% unique(mean_characteristics_nocomp$species)) %>%
                     dplyr::summarise(total_weighted = sum(mean_RMSE_weighted),
                                      total_not_weighted = sum(mean_RMSE_no_weighted))
  
  diameter <- read.csv(file = "data/diameter_nocomp_RMSE_AIC_results.csv")
  diameter <- diameter %>% dplyr::filter(species %in% unique(mean_characteristics_nocomp$species)) %>%
                           dplyr::summarise(total_weighted = sum(mean_RMSE_weighted),
                                            total_not_weighted = sum(mean_RMSE_no_weighted))
  
  RMSE_protocol <- rbind(height, ratio, diameter)
  RMSE_protocol$characteristic <- c("height", "ratio", "diameter")
  
  write.csv(file = "output/protocol_parameter.csv", RMSE_protocol)
  
  return(RMSE_protocol)
  
}


competition_index <- function(mean_characteristics_nocomp) {
  
  height <- read.csv(file = "data/height_power_comp_RMSE_AIC_results.csv")
  height <- height %>% dplyr::filter(species %in% unique(mean_characteristics_nocomp$species)) %>%
                       dplyr::summarise(total_AIC_c1 = sum(c1_mean_RMSE_weighted),
                                        total_AIC_c2 = sum(c2_mean_RMSE_weighted))
  
  ratio <- read.csv(file = "data/ratio_compRMSE_AIC_results.csv")
  ratio <- ratio %>% dplyr::filter(species %in% unique(mean_characteristics_nocomp$species)) %>%
                     dplyr::summarise(total_AIC_c1 = sum(c1_mean_RMSE_weighted),
                                      total_AIC_c2 = sum(c2_mean_RMSE_weighted))
  
  diam <- read.csv(file = "data/diameter_comp_RMSE_AIC_results.csv")
  diam <- diam %>% dplyr::filter(species %in% unique(mean_characteristics_nocomp$species)) %>%
                   dplyr::summarise(total_AIC_c1 = sum(c1_mean_RMSE_weighted),
                                    total_AIC_c2 = sum(c2_mean_RMSE_weighted))
  
  CI_selection <- rbind(height, ratio, diam)
  CI_selection$characteristic <- c("height", "ratio", "diameter")
  
  write.csv(file = "output/ci_selection.csv", CI_selection)
  
  return(CI_selection)
  
}


compute_mean_parameters_nocomp <- function(mean_characteristics_nocomp) {
  
  
  ## Height 
  height_asymptot <- read.csv(file = "data/height_asympt_nocomp_weighted_parameters.csv")
  height_nocomp_parameters <- height_asymptot %>% dplyr::select(-X) %>%
                                                  dplyr::filter(species %in% mean_characteristics_nocomp$species) %>%
                                                  dplyr::group_by(species) %>%
                                                  dplyr::summarise(a1_h_mean = mean(b1),
                                                                   a2_h_mean = mean(b2),
                                                                   a3_h_mean = mean(b3),
                                                                   sd_h_a1 = sd(b1),
                                                                   sd_h_a2 = sd(b2),
                                                                   sd_h_a3 = sd(b3)) %>%
                                                 dplyr::ungroup()
  
  ### Diameter
  diameter_power <- read.csv(file = "data/diameter_nocomp_weighted_parameters.csv")
  diameter_nocomp_parameters <- diameter_power %>% dplyr::select(-X) %>%
                                                   dplyr::filter(species %in% mean_characteristics_nocomp$species) %>%
                                                   dplyr::group_by(species) %>%
                                                   dplyr::summarise(a1_cd_mean = mean(a1),
                                                                    a2_cd_mean = mean(a2),
                                                                    sd_cd_a1 = sd(a1),
                                                                    sd_cd_a2 = sd(a2)) %>%
                                                   dplyr::ungroup()
  
  ## Ratio
  ratio_beta <- read.csv(file = "data/ratio_nocomp_weighted_parameters.csv")
  ratio_nocomp_parameters <- ratio_beta %>% dplyr::select(-X) %>%
                                            dplyr::filter(species %in% mean_characteristics_nocomp$species) %>%
                                            dplyr::group_by(species) %>%
                                            dplyr::summarise(a1_rcd_mean = mean(a1),
                                                             a2_rcd_mean = mean(a2),
                                                             sd_rcd_a1 = sd(a1),
                                                             sd_rcd_a2 = sd(a2)) %>%
                                            dplyr::ungroup()
  
  
  nocomp_parameters <- join_all(list(height_nocomp_parameters, ratio_nocomp_parameters, diameter_nocomp_parameters), 
                                by = "species", type = "left")
  
  write.csv(file = "mean_parameters_nocompetition.csv", nocomp_parameters)
  
  return(nocomp_parameters)
  
}


compute_mean_parameters_comp <- function(species_list_comp) {
  
  
}


pearson_correlation <- function(mean_characteristics_nocomp) {
  
  ### Loading species' functional traits and climate niche
  ## Computing petp as map - etp
  traits_climate <- read.csv(file = "data/sp_traits_complete.csv")
  traits_climate <- traits_climate %>% dplyr::select(species, 
                                                     mat, map, etp,
                                                     sla, ssd, shade_tol_mean) %>%
                                       dplyr::mutate(petp = map - etp) %>%
                                       dplyr::select(-map, -etp)
  
  ### Creating a single file 
  df <- left_join(mean_characteristics_nocomp, traits_climate, by = "species")
  
  ### Separating angiosperms from gymnosperms
  angio <- df %>% dplyr::filter(new_gp == "angio")
  gymno <- df %>% dplyr::filter(new_gp == "gymno")
  
  ### Computing crown volumes
  angio <- angio %>% dplyr::mutate(depth_15 = height_15 * ratio_15,
                                   volume_15 = 4/3 * pi * (diameter_15/2) * (diameter_15/2) * (depth_15/2))
  
  gymno <-  gymno %>% dplyr::mutate(depth_15 = height_15 * ratio_15,
                                    volume_15 = 1/2 * pi * (diameter_15/2)*(diameter_15/2) * (depth_15))
  
  
  
  ### Performing pearson correlations
  cor_matrix_angio <- angio %>% dplyr::select(hmax,
                                              height_15,
                                              ratio_15,
                                              diameter_15,
                                              volume_15,
                                              mat,
                                              petp,
                                              sla, 
                                              ssd,
                                              shade_tol_mean)
  
  cor_matrix_values_angio <- cor(cor_matrix_angio, use = "complete.obs", method = "pearson")
  
  
  cor_matrix_gymno <- gymno %>% dplyr::select(hmax,
                                              height_15,
                                              ratio_15,
                                              diameter_15,
                                              volume_15,
                                              mat,
                                              petp,
                                              sla, 
                                              ssd,
                                              shade_tol_mean)
  
  cor_matrix_values_gymno <- cor(cor_matrix_gymno, use = "complete.obs", method = "pearson")
  
  
  ### Plotting results
  par(mfrow = c(1,2))
  corrplot(cor_matrix_values_angio, 
           type = "lower", 
           method = "number")
  
  corrplot(cor_matrix_values_gymno, 
           type = "lower", 
           method = "number")
  
  cor_matrix_values_angio <- as.data.frame(cor_matrix_values_angio)
  cor_matrix_values_angio$group <- "angio"
  cor_matrix_values_gymno <- as.data.frame(cor_matrix_values_gymno)
  cor_matrix_values_gymno$group <- "gymno"
  
  cor_matrix_values <- rbind(cor_matrix_values_angio, cor_matrix_values_gymno)
  
  return(cor_matrix_values)
  
}




univariate_models <- function() {
  
  ### Computing crown characteristics for each species based on the 100 values obtained from the resampling process
  # use species list previously defined
  
  sp_list <- read.csv(file = "output/mean_crown_characteristics_nocomp.csv")
  sp_list <- sp_list %>% dplyr::select(species)
  
  ## Height 
  height_asymptot <- read.csv(file = "data/height_asympt_nocomp_weighted_parameters.csv")
  height_nocomp_estimates <- height_asymptot %>% dplyr::select(-X) %>%
                                                 dplyr::filter(species %in% sp_list$species) %>%
                                                 dplyr::mutate(hmax = b1,
                                                               height_15 = 1.3 + b1 * (1-exp(-b2 * 15)) ^ b3,
                                                               height_30 = 1.3 + b1 * (1-exp(-b2 * 30)) ^ b3) %>%
                                                 dplyr::select(-b1, -b2, -b3) %>%
                                                 dplyr::ungroup() %>%
                                                 dplyr::arrange(species)
  
  ### Diameter
  diameter_power <- read.csv(file = "data/diameter_nocomp_weighted_parameters.csv")
  diameter_nocomp_estimates <- diameter_power %>% dplyr::select(-X) %>%
                                                  dplyr::filter(species %in% sp_list$species) %>%
                                                  dplyr::group_by(species) %>%
                                                  dplyr::mutate(diameter_15 = a1 * (15 ^ a2),
                                                                diameter_30 = a2 * (30 ^ a2)) %>%
                                                  dplyr::select(-a1, -a2) %>%
                                                  dplyr::ungroup() %>%
                                                  dplyr::arrange(species) %>%
                                                  dplyr::select(-species)
  
  ## Ratio
  ratio_beta <- read.csv(file = "data/ratio_nocomp_weighted_parameters.csv")
  ratio_nocomp_estimates <- ratio_beta %>% dplyr::select(-X) %>%
                                           dplyr::filter(species %in% sp_list$species) %>%
                                           dplyr::group_by(species) %>%
                                           dplyr::mutate(ratio_15 = (exp(a1) + a2 * 15) / (1 + exp(a1) + a2 * 15),
                                                         ratio_30 = (exp(a1) + a2 * 30) / (1 + exp(a1) + a2 * 30)) %>%
                                           dplyr::select(-a1, -a2) %>%
                                           dplyr::ungroup()  %>%
                                           dplyr::arrange(species) %>%
                                           dplyr::select(-species)
  
  
  crown_characteristics <- cbind(height_nocomp_estimates, ratio_nocomp_estimates, diameter_nocomp_estimates)
  all_crown_characteristics <- na.omit(crown_characteristics)
  
  ### Scaling and centering crown characteristics at the functional group scale
  ## Summarizing crown characteristics at the species scale
  
  functional_groups <- read.csv(file = "data/functional_groups.csv")
  functional_groups <- functional_groups %>% dplyr::select(sp, new_gp)
  
  all_crown_characteristics <- left_join(all_crown_characteristics, functional_groups, by = c("species" = "sp"))
  
  angiosperms <- all_crown_characteristics %>% dplyr::filter(new_gp == "angio") %>%
    
                                               dplyr::mutate(hmax_s = scale(hmax, center = TRUE, scale = TRUE)[, 1],
                                                             height_15_s = scale(height_15, center = TRUE, scale = TRUE)[, 1],
                                                             height_30_s = scale(height_30, center = TRUE, scale = TRUE)[, 1],
                                                             diameter_15_s = scale(diameter_15, center = TRUE, scale = TRUE)[, 1],
                                                             diameter_30_s = scale(diameter_30, center = TRUE, scale = TRUE)[, 1],
                                                             ratio_15_s = scale(ratio_15, center = TRUE, scale = TRUE)[, 1],
                                                             ratio_30_s = scale(ratio_30, center = TRUE, scale = TRUE)[, 1]) %>%
                                              
                                               dplyr::select(-hmax, -height_15, -height_30, 
                                                             -ratio_15, -ratio_30,
                                                             -diameter_15, -diameter_30) %>%
                                              
                                               dplyr::group_by(species) %>%
                                               dplyr::summarise(hmax = mean(hmax_s), cv_hmax = abs(sd(hmax_s)/mean(hmax_s)),
                                                                height_15 = mean(height_15_s), cv_height_15 = abs(sd(height_15_s)/mean(height_15_s)),
                                                                height_30 = mean(height_30_s), cv_height_30 = abs(sd(height_30_s)/mean(height_30_s)),
                                                                ratio_15 = mean(ratio_15_s), cv_ratio_15 = abs(sd(ratio_15_s)/mean(ratio_15_s)),
                                                                ratio_30 = mean(ratio_30_s), cv_ratio_30 = abs(sd(ratio_30_s)/mean(ratio_30_s)),
                                                                diameter_15 = mean(diameter_15_s), cv_diameter_15 = abs(sd(diameter_15_s)/mean(diameter_15_s)),
                                                                diameter_30 = mean(diameter_30_s), cv_diameter_30 = abs(sd(diameter_30_s)/mean(diameter_30_s))) %>%
                                               dplyr::ungroup()
  
  gymnosperms <- all_crown_characteristics %>% dplyr::filter(new_gp == "gymno") %>%
    
                                               dplyr::mutate(hmax_s = scale(hmax, center = TRUE, scale = TRUE)[, 1],
                                                             height_15_s = scale(height_15, center = TRUE, scale = TRUE)[, 1],
                                                             height_30_s = scale(height_30, center = TRUE, scale = TRUE)[, 1],
                                                             diameter_15_s = scale(diameter_15, center = TRUE, scale = TRUE)[, 1],
                                                             diameter_30_s = scale(diameter_30, center = TRUE, scale = TRUE)[, 1],
                                                             ratio_15_s = scale(ratio_15, center = TRUE, scale = TRUE)[, 1],
                                                             ratio_30_s = scale(ratio_30, center = TRUE, scale = TRUE)[, 1]) %>%
                                              
                                               dplyr::select(-hmax, -height_15, -height_30, 
                                                             -ratio_15, -ratio_30,
                                                             -diameter_15, -diameter_30) %>%
                                              
                                               dplyr::group_by(species) %>%
                                               dplyr::summarise(hmax = mean(hmax_s), cv_hmax = abs(sd(hmax_s)/mean(hmax_s)),
                                                                height_15 = mean(height_15_s), cv_height_15 = abs(sd(height_15_s)/mean(height_15_s)),
                                                                height_30 = mean(height_30_s), cv_height_30 = abs(sd(height_30_s)/mean(height_30_s)),
                                                                ratio_15 = mean(ratio_15_s), cv_ratio_15 = abs(sd(ratio_15_s)/mean(ratio_15_s)),
                                                                ratio_30 = mean(ratio_30_s), cv_ratio_30 = abs(sd(ratio_30_s)/mean(ratio_30_s)),
                                                                diameter_15 = mean(diameter_15_s), cv_diameter_15 = abs(sd(diameter_15_s)/mean(diameter_15_s)),
                                                                diameter_30 = mean(diameter_30_s), cv_diameter_30 = abs(sd(diameter_30_s)/mean(diameter_30_s))) %>%
                                               dplyr::ungroup()
  
  
  
  
  
  
  ### Loading species' functional traits and climate niche
  ## Computing petp as map - etp
  traits_climate <- read.csv(file = "data/sp_traits_complete.csv")
  traits_data <- traits_climate %>% dplyr::select(species, 
                                                     mat, map, etp,
                                                     sla, ssd, shade_tol_mean) %>%
                                    dplyr::mutate(petp = map - etp) %>%
                                    dplyr::select(-map, -etp)
  
  traits_angio <- traits_data[traits_data$species %in% unique(angiosperms$species),]
  
  mat_a <- as.data.frame(scale(traits_angio$mat, center = TRUE, scale = TRUE))
  names(mat_a) <- "mat"
  petp_a <- as.data.frame(scale(traits_angio$petp, center = TRUE, scale = TRUE))
  names(petp_a) <- "petp"
  sla_a <- as.data.frame(scale(traits_angio$sla, center = TRUE, scale = TRUE))
  names(sla_a) <- "sla"
  ssd_a <- as.data.frame(scale(traits_angio$ssd, center = TRUE, scale = TRUE))
  names(ssd_a) <- "ssd"
  sh_a <- as.data.frame(scale(traits_angio$shade_tol_mean, center = TRUE, scale = TRUE))
  names(sh_a) <- "sh"
  
  angio_scaled_summary <- cbind(angiosperms, mat_a, petp_a, sla_a, ssd_a, sh_a)
  angio_scaled_summary <- na.omit(angio_scaled_summary)
  
  traits_gymno <- traits_data[traits_data$species %in% unique(gymnosperms$species),]
  
  mat_g <- as.data.frame(scale(traits_gymno$mat, center = TRUE, scale = TRUE))
  names(mat_g) <- "mat"
  petp_g <- as.data.frame(scale(traits_gymno$petp, center = TRUE, scale = TRUE))
  names(petp_g) <- "petp"
  sla_g <- as.data.frame(scale(traits_gymno$sla, center = TRUE, scale = TRUE))
  names(sla_g) <- "sla"
  ssd_g <- as.data.frame(scale(traits_gymno$ssd, center = TRUE, scale = TRUE))
  names(ssd_g) <- "ssd"
  sh_g <- as.data.frame(scale(traits_gymno$shade_tol_mean, center = TRUE, scale = TRUE))
  names(sh_g) <- "sh"
  
  gymno_scaled_summary <- cbind(gymnosperms, mat_g, petp_g, sla_g, ssd_g, sh_g)
  gymno_scaled_summary <- na.omit(gymno_scaled_summary)
  
  
  # climatic models
  mod_hmax_mat_a <- lm(hmax ~ mat, weights = (1/cv_hmax), data = angio_scaled_summary)
  mod_h15_mat_a <- lm(height_15 ~ mat, weights = (1/cv_height_15), data = angio_scaled_summary)
  mod_r15_mat_a <- lm(ratio_15 ~ mat, weights = (1/cv_ratio_15), data = angio_scaled_summary)
  mod_d15_mat_a <- lm(diameter_15 ~ mat, weights = (1/cv_diameter_15), data = angio_scaled_summary)
  
  
  mod_hmax_mat_g <- lm(hmax ~ mat, weights = (1/cv_hmax), data = gymno_scaled_summary)
  mod_h15_mat_g <- lm(height_15 ~ mat, weights = (1/cv_height_15), data = gymno_scaled_summary)
  mod_r15_mat_g <- lm(ratio_15 ~ mat, weights = (1/cv_ratio_15), data = gymno_scaled_summary)
  mod_d15_mat_g <- lm(diameter_15 ~ mat, weights = (1/cv_diameter_15), data = gymno_scaled_summary)
  
  
  mod_hmax_petp_a <- lm(hmax ~ petp, weights = (1/cv_hmax), data = angio_scaled_summary)
  mod_h15_petp_a <- lm(height_15 ~ petp, weights = (1/cv_height_15), data = angio_scaled_summary)
  mod_r15_petp_a <- lm(ratio_15 ~ petp, weights = (1/cv_ratio_15), data = angio_scaled_summary)
  mod_d15_petp_a <- lm(diameter_15 ~ petp, weights = (1/cv_diameter_15), data = angio_scaled_summary)
  
  
  mod_hmax_petp_g <- lm(hmax ~ petp, weights = (1/cv_hmax), data = gymno_scaled_summary)
  mod_h15_petp_g <- lm(height_15 ~ petp, weights = (1/cv_height_15), data = gymno_scaled_summary)
  mod_r15_petp_g <- lm(ratio_15 ~ petp, weights = (1/cv_ratio_15), data = gymno_scaled_summary)
  mod_d15_petp_g <- lm(diameter_15 ~ petp, weights = (1/cv_diameter_15), data = gymno_scaled_summary)
  
  
  # functional traits models
  mod_hmax_sla_a <- lm(hmax ~ sla, weights = (1/cv_hmax), data = angio_scaled_summary)
  mod_h15_sla_a <- lm(height_15 ~ sla, weights = (1/cv_height_15), data = angio_scaled_summary)
  mod_r15_sla_a <- lm(ratio_15 ~ sla, weights = (1/cv_ratio_15), data = angio_scaled_summary)
  mod_d15_sla_a <- lm(diameter_15 ~ sla, weights = (1/cv_diameter_15), data = angio_scaled_summary)
  
  
  mod_hmax_sla_g <- lm(hmax ~ sla, weights = (1/cv_hmax), data = gymno_scaled_summary)
  mod_h15_sla_g <- lm(height_15 ~ sla, weights = (1/cv_height_15), data = gymno_scaled_summary)
  mod_r15_sla_g <- lm(ratio_15 ~ sla, weights = (1/cv_ratio_15), data = gymno_scaled_summary)
  mod_d15_sla_g <- lm(diameter_15 ~ sla, weights = (1/cv_diameter_15), data = gymno_scaled_summary)
  
  mod_hmax_ssd_a <- lm(hmax ~ ssd, weights = (1/cv_hmax), data = angio_scaled_summary)
  mod_h15_ssd_a <- lm(height_15 ~ ssd, weights = (1/cv_height_15), data = angio_scaled_summary)
  mod_r15_ssd_a <- lm(ratio_15 ~ ssd, weights = (1/cv_ratio_15), data = angio_scaled_summary)
  mod_d15_ssd_a <- lm(diameter_15 ~ ssd, weights = (1/cv_diameter_15), data = angio_scaled_summary)
  
  
  mod_hmax_ssd_g <- lm(hmax ~ ssd, weights = (1/cv_hmax), data = gymno_scaled_summary)
  mod_h15_ssd_g <- lm(height_15 ~ ssd, weights = (1/cv_height_15), data = gymno_scaled_summary)
  mod_r15_ssd_g <- lm(ratio_15 ~ ssd, weights = (1/cv_ratio_15), data = gymno_scaled_summary)
  mod_d15_ssd_g <- lm(diameter_15 ~ ssd, weights = (1/cv_diameter_15), data = gymno_scaled_summary)
  
  
  
  # shade tolerance models
  mod_hmax_sh_a <- lm(hmax ~ sh, weights = (1/cv_hmax), data = angio_scaled_summary)
  mod_h15_sh_a <- lm(height_15 ~ sh, weights = (1/cv_height_15), data = angio_scaled_summary)
  mod_r15_sh_a <- lm(ratio_15 ~ sh, weights = (1/cv_ratio_15), data = angio_scaled_summary)
  mod_d15_sh_a <- lm(diameter_15 ~ sh, weights = (1/cv_diameter_15), data = angio_scaled_summary)
  
  
  mod_hmax_sh_g <- lm(hmax ~ sh, weights = (1/cv_hmax), data = gymno_scaled_summary)
  mod_h15_sh_g <- lm(height_15 ~ sh, weights = (1/cv_height_15), data = gymno_scaled_summary)
  mod_r15_sh_g <- lm(ratio_15 ~ sh, weights = (1/cv_ratio_15), data = gymno_scaled_summary)
  mod_d15_sh_g <- lm(diameter_15 ~ sh, weights = (1/cv_diameter_15), data = gymno_scaled_summary)
  
  
  # storing models' results
  models <- c("max. height", "height 15", "crown ratio", "crown diameter")
  
  angio_models_mat_15 <- list(mod_hmax_mat_a, mod_h15_mat_a, mod_r15_mat_a, mod_d15_mat_a)
  angio_models_petp_15 <- list(mod_hmax_petp_a, mod_h15_petp_a, mod_r15_petp_a, mod_d15_petp_a)
  angio_models_sla_15 <- list(mod_hmax_sla_a, mod_h15_sla_a, mod_r15_sla_a, mod_d15_sla_a)
  angio_models_ssd_15 <- list(mod_hmax_ssd_a, mod_h15_ssd_a, mod_r15_ssd_a, mod_d15_ssd_a)
  angio_models_sh_15 <- list(mod_hmax_sh_a, mod_h15_sh_a, mod_r15_sh_a, mod_d15_sh_a)
  
  gymno_models_mat_15 <- list(mod_hmax_mat_g, mod_h15_mat_g, mod_r15_mat_g, mod_d15_mat_g)
  gymno_models_petp_15 <- list(mod_hmax_petp_g, mod_h15_petp_g, mod_r15_petp_g, mod_d15_petp_g)
  gymno_models_sla_15 <- list(mod_hmax_sla_g, mod_h15_sla_g, mod_r15_sla_g, mod_d15_sla_g)
  gymno_models_ssd_15 <- list(mod_hmax_ssd_g, mod_h15_ssd_g, mod_r15_ssd_g, mod_d15_ssd_g)
  gymno_models_sh_15 <- list(mod_hmax_sh_g, mod_h15_sh_g, mod_r15_sh_g, mod_d15_sh_g)
  
  
  # climatic models - results extraction
  angio_results_mat_15 <- as.data.frame(matrix(nrow = length(models), ncol = 5))
  names(angio_results_mat_15) <- c("mod", "est.inf", "est.sup", "est", "sign")
  
  gymno_results_mat_15 <- as.data.frame(matrix(nrow = length(models), ncol = 5))
  names(gymno_results_mat_15) <- c("mod", "est.inf", "est.sup", "est", "sign")
  
  
  angio_results_petp_15 <- as.data.frame(matrix(nrow = length(models), ncol = 5))
  names(angio_results_petp_15) <- c("mod", "est.inf", "est.sup", "est", "sign")
  
  gymno_results_petp_15 <- as.data.frame(matrix(nrow = length(models), ncol = 5))
  names(gymno_results_petp_15) <- c("mod", "est.inf", "est.sup", "est", "sign")
  
  
  
  for (i in 1:length(models)) {
    
    # angiosperms
    angio_results_mat_15[i, "mod"] <- models[i]
    angio_results_petp_15[i, "mod"] <- models[i]
    mod_mat_15_a <- angio_models_mat_15[[i]]
    mod_petp_15_a <- angio_models_petp_15[[i]]
    
    angio_results_mat_15[i, "est.inf"] = confint.lm(mod_mat_15_a)["mat", "2.5 %"]
    angio_results_mat_15[i, "est.sup"] = confint.lm(mod_mat_15_a)["mat", "97.5 %"]
    angio_results_mat_15[i, "est"] = summary(mod_mat_15_a)$coefficients["mat", "Estimate"]
    angio_results_mat_15[i, "sign"] = summary(mod_mat_15_a)$coefficients["mat", "Pr(>|t|)"]
    
    angio_results_petp_15[i, "est.inf"] = confint.lm(mod_petp_15_a)["petp", "2.5 %"]
    angio_results_petp_15[i, "est.sup"] = confint.lm(mod_petp_15_a)["petp", "97.5 %"]
    angio_results_petp_15[i, "est"] = summary(mod_petp_15_a)$coefficients["petp", "Estimate"]
    angio_results_petp_15[i, "sign"] = summary(mod_petp_15_a)$coefficients["petp", "Pr(>|t|)"]
    
    
    # gymnosperms
    gymno_results_mat_15[i, "mod"] <- models[i]
    gymno_results_petp_15[i, "mod"] <- models[i]
    mod_mat_15_g <- gymno_models_mat_15[[i]]
    mod_petp_15_g <- gymno_models_petp_15[[i]]
    
    gymno_results_mat_15[i, "est.inf"] = confint.lm(mod_mat_15_g)["mat", "2.5 %"]
    gymno_results_mat_15[i, "est.sup"] = confint.lm(mod_mat_15_g)["mat", "97.5 %"]
    gymno_results_mat_15[i, "est"] = summary(mod_mat_15_g)$coefficients["mat", "Estimate"]
    gymno_results_mat_15[i, "sign"] = summary(mod_mat_15_g)$coefficients["mat", "Pr(>|t|)"]
    
    gymno_results_petp_15[i, "est.inf"] = confint.lm(mod_petp_15_g)["petp", "2.5 %"]
    gymno_results_petp_15[i, "est.sup"] = confint.lm(mod_petp_15_g)["petp", "97.5 %"]
    gymno_results_petp_15[i, "est"] = summary(mod_petp_15_g)$coefficients["petp", "Estimate"]
    gymno_results_petp_15[i, "sign"] = summary(mod_petp_15_g)$coefficients["petp", "Pr(>|t|)"]
    
    
  }
  
  angio_results_mat_15 <- angio_results_mat_15 %>% mutate(sign = ifelse(sign > 0.05, "NA", 
                                                                        ifelse(sign >= 0.01, "*", 
                                                                               ifelse(sign >= 0.001, "**", "***"))))
  
  
  angio_results_petp_15 <- angio_results_petp_15 %>% mutate(sign = ifelse(sign > 0.05, "NA", 
                                                                          ifelse(sign >= 0.01, "*", 
                                                                                 ifelse(sign >= 0.001, "**", "***"))))
  
  
  gymno_results_mat_15 <- gymno_results_mat_15 %>% mutate(sign = ifelse(sign > 0.05, "NA", 
                                                                        ifelse(sign >= 0.01, "*", 
                                                                               ifelse(sign >= 0.001, "**", "***"))))
  
  
  gymno_results_petp_15 <- gymno_results_petp_15 %>% mutate(sign = ifelse(sign > 0.05, "NA", 
                                                                          ifelse(sign >= 0.01, "*", 
                                                                                 ifelse(sign >= 0.001, "**", "***"))))
  
  
  
  
  # functional traits models - results extraction
  angio_results_sla_15 <- as.data.frame(matrix(nrow = length(models), ncol = 5))
  names(angio_results_sla_15) <- c("mod", "est.inf", "est.sup", "est", "sign")
  
  gymno_results_sla_15 <- as.data.frame(matrix(nrow = length(models), ncol = 5))
  names(gymno_results_sla_15) <- c("mod", "est.inf", "est.sup", "est", "sign")
  
  
  angio_results_ssd_15 <- as.data.frame(matrix(nrow = length(models), ncol = 5))
  names(angio_results_ssd_15) <- c("mod", "est.inf", "est.sup", "est", "sign")
  
  gymno_results_ssd_15 <- as.data.frame(matrix(nrow = length(models), ncol = 5))
  names(gymno_results_ssd_15) <- c("mod", "est.inf", "est.sup", "est", "sign")
  
  
  
  for (i in 1:length(models)) {
    
    # angiosperms
    angio_results_sla_15[i, "mod"] <- models[i]
    angio_results_ssd_15[i, "mod"] <- models[i]
    mod_sla_15_a <- angio_models_sla_15[[i]]
    mod_ssd_15_a <- angio_models_ssd_15[[i]]
    
    angio_results_sla_15[i, "est.inf"] = confint.lm(mod_sla_15_a)["sla", "2.5 %"]
    angio_results_sla_15[i, "est.sup"] = confint.lm(mod_sla_15_a)["sla", "97.5 %"]
    angio_results_sla_15[i, "est"] = summary(mod_sla_15_a)$coefficients["sla", "Estimate"]
    angio_results_sla_15[i, "sign"] = summary(mod_sla_15_a)$coefficients["sla", "Pr(>|t|)"]
    
    angio_results_ssd_15[i, "est.inf"] = confint.lm(mod_ssd_15_a)["ssd", "2.5 %"]
    angio_results_ssd_15[i, "est.sup"] = confint.lm(mod_ssd_15_a)["ssd", "97.5 %"]
    angio_results_ssd_15[i, "est"] = summary(mod_ssd_15_a)$coefficients["ssd", "Estimate"]
    angio_results_ssd_15[i, "sign"] = summary(mod_ssd_15_a)$coefficients["ssd", "Pr(>|t|)"]
    
    
    # gymnosperms
    gymno_results_sla_15[i, "mod"] <- models[i]
    gymno_results_ssd_15[i, "mod"] <- models[i]
    mod_sla_15_g <- gymno_models_sla_15[[i]]
    mod_ssd_15_g <- gymno_models_ssd_15[[i]]
    
    gymno_results_sla_15[i, "est.inf"] = confint.lm(mod_sla_15_g)["sla", "2.5 %"]
    gymno_results_sla_15[i, "est.sup"] = confint.lm(mod_sla_15_g)["sla", "97.5 %"]
    gymno_results_sla_15[i, "est"] = summary(mod_sla_15_g)$coefficients["sla", "Estimate"]
    gymno_results_sla_15[i, "sign"] = summary(mod_sla_15_g)$coefficients["sla", "Pr(>|t|)"]
    
    gymno_results_ssd_15[i, "est.inf"] = confint.lm(mod_ssd_15_g)["ssd", "2.5 %"]
    gymno_results_ssd_15[i, "est.sup"] = confint.lm(mod_ssd_15_g)["ssd", "97.5 %"]
    gymno_results_ssd_15[i, "est"] = summary(mod_ssd_15_g)$coefficients["ssd", "Estimate"]
    gymno_results_ssd_15[i, "sign"] = summary(mod_ssd_15_g)$coefficients["ssd", "Pr(>|t|)"]   
    
  }
  
  angio_results_sla_15 <- angio_results_sla_15 %>% mutate(sign = ifelse(sign > 0.05, "NA", 
                                                                        ifelse(sign >= 0.01, "*", 
                                                                               ifelse(sign >= 0.001, "**", "***"))))
  
  
  angio_results_ssd_15 <- angio_results_ssd_15 %>% mutate(sign = ifelse(sign > 0.05, "NA", 
                                                                        ifelse(sign >= 0.01, "*", 
                                                                               ifelse(sign >= 0.001, "**", "***"))))
  
  
  gymno_results_sla_15 <- gymno_results_sla_15 %>% mutate(sign = ifelse(sign > 0.05, "NA", 
                                                                        ifelse(sign >= 0.01, "*", 
                                                                               ifelse(sign >= 0.001, "**", "***"))))
  
  
  gymno_results_ssd_15 <- gymno_results_ssd_15 %>% mutate(sign = ifelse(sign > 0.05, "NA", 
                                                                        ifelse(sign >= 0.01, "*", 
                                                                               ifelse(sign >= 0.001, "**", "***"))))
  
  
  
  
  
  # shade tolerance models - results extraction
  angio_results_sh_15 <- as.data.frame(matrix(nrow = length(models), ncol = 5))
  names(angio_results_sh_15) <- c("mod", "est.inf", "est.sup", "est", "sign")
  
  gymno_results_sh_15 <- as.data.frame(matrix(nrow = length(models), ncol = 5))
  names(gymno_results_sh_15) <- c("mod", "est.inf", "est.sup", "est", "sign")
  
  
  for (i in 1:length(models)) {
    
    # angiosperms
    angio_results_sh_15[i, "mod"] <- models[i]
    mod_15_a <- angio_models_sh_15[[i]]
    
    angio_results_sh_15[i, "est.inf"] = confint.lm(mod_15_a)["sh", "2.5 %"]
    angio_results_sh_15[i, "est.sup"] = confint.lm(mod_15_a)["sh", "97.5 %"]
    angio_results_sh_15[i, "est"] = summary(mod_15_a)$coefficients["sh", "Estimate"]
    angio_results_sh_15[i, "sign"] = summary(mod_15_a)$coefficients["sh", "Pr(>|t|)"]
    
    
    # gymnosperms
    gymno_results_sh_15[i, "mod"] <- models[i]
    mod_15_g <- gymno_models_sh_15[[i]]
    
    gymno_results_sh_15[i, "est.inf"] = confint.lm(mod_15_g)["sh", "2.5 %"]
    gymno_results_sh_15[i, "est.sup"] = confint.lm(mod_15_g)["sh", "97.5 %"]
    gymno_results_sh_15[i, "est"] = summary(mod_15_g)$coefficients["sh", "Estimate"]
    gymno_results_sh_15[i, "sign"] = summary(mod_15_g)$coefficients["sh", "Pr(>|t|)"]
    
    
  }
  
  angio_results_sh_15 <- angio_results_sh_15 %>% mutate(sign = ifelse(sign > 0.05, "NA", 
                                                                      ifelse(sign >= 0.01, "*", 
                                                                             ifelse(sign >= 0.001, "**", "***"))))
  
  
  gymno_results_sh_15 <- gymno_results_sh_15 %>% mutate(sign = ifelse(sign > 0.05, "NA", 
                                                                      ifelse(sign >= 0.01, "*", 
                                                                             ifelse(sign >= 0.001, "**", "***"))))
  
  plot_mat_15_a <- angio_results_mat_15 %>% dplyr::mutate(label.pos = max(est.sup, na.rm = TRUE)) %>%
                                            dplyr::mutate(label.pos = label.pos + 0.5*(max(est.sup, na.rm = T) - min(est.inf, na.rm = T))) %>%
                                            dplyr::mutate(mod = factor(models, levels = c("crown diameter", "crown ratio", "height 15", "max. height"))) %>%
                                            ggplot(aes(x = mod, y = est)) + ylim(-4, 4) +
                                            geom_errorbar(aes(ymin = est.inf, ymax = est.sup), colour = "black", width = 0, size = 0.2) + 
                                            geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", size = 0.3) + 
                                            geom_point(size = 1.5, fill = "chocolate2", shape = 21) +  
                                            geom_text(aes(y = 0, label = sign), nudge_y = 3.5, size = 5, color = "chocolate2") +
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
  
  plot_mat_15_a
  
  plot_petp_15_a <- angio_results_petp_15 %>% dplyr::mutate(label.pos = max(est.sup, na.rm = TRUE)) %>%
                                              dplyr::mutate(label.pos = label.pos + 0.5*(max(est.sup, na.rm = T) - min(est.inf, na.rm = T))) %>%
                                              dplyr::mutate(mod = factor(models, levels = c("crown diameter", "crown ratio", "height 15", "max. height"))) %>%
                                              ggplot(aes(x = mod, y = est)) + ylim(-4, 4) +
                                              geom_errorbar(aes(ymin = est.inf, ymax = est.sup), colour = "black", width = 0, size = 0.2) + 
                                              geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", size = 0.3) + 
                                              geom_point(size = 1.5, fill = "chocolate2", shape = 21)  +  
                                              geom_text(aes(y = 0, label = sign), nudge_y = 3.5, size = 5, color = "chocolate2") +
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
  plot_petp_15_a
  
  
  plot_sla_15_a <- angio_results_sla_15 %>% dplyr::mutate(label.pos = max(est.sup, na.rm = TRUE)) %>%
                                            dplyr::mutate(label.pos = label.pos + 0.5*(max(est.sup, na.rm = T) - min(est.inf, na.rm = T))) %>%
                                            dplyr::mutate(mod = factor(models, levels = c("crown diameter", "crown ratio", "height 15", "max. height"))) %>%
                                            ggplot(aes(x = mod, y = est)) + ylim(-4, 4) +
                                            geom_errorbar(aes(ymin = est.inf, ymax = est.sup), colour = "black", width = 0, size = 0.2) + 
                                            geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", size = 0.3) + 
                                            geom_point(size = 1.5, fill = "chocolate2", shape = 21)  +  
                                            geom_text(aes(y = 0, label = sign), nudge_y = 3.5, size = 5, color = "chocolate2") +
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
  plot_sla_15_a
  
  
  plot_ssd_15_a <- angio_results_ssd_15 %>% dplyr::mutate(label.pos = max(est.sup, na.rm = TRUE)) %>%
                                            dplyr::mutate(label.pos = label.pos + 0.5*(max(est.sup, na.rm = T) - min(est.inf, na.rm = T))) %>%
                                            dplyr::mutate(mod = factor(models, levels = c("crown diameter", "crown ratio", "height 15", "max. height"))) %>%
                                            ggplot(aes(x = mod, y = est)) + ylim(-4, 4) +
                                            geom_errorbar(aes(ymin = est.inf, ymax = est.sup), colour = "black", width = 0, size = 0.2) + 
                                            geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", size = 0.3) + 
                                            geom_point(size = 1.5, fill = "chocolate2", shape = 21)  +  
                                            geom_text(aes(y = 0, label = sign), nudge_y = 3.5, size = 5, color = "chocolate2") +
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
  plot_ssd_15_a
  
  
  plot_sh_15_a <- angio_results_sh_15 %>% dplyr::mutate(label.pos = max(est.sup, na.rm = TRUE)) %>%
                                          dplyr::mutate(label.pos = label.pos + 0.5*(max(est.sup, na.rm = T) - min(est.inf, na.rm = T))) %>%
                                          dplyr::mutate(mod = factor(models, levels = c("crown diameter", "crown ratio", "height 15", "max. height"))) %>%
                                          ggplot(aes(x = mod, y = est)) + ylim(-4, 4) +
                                          geom_errorbar(aes(ymin = est.inf, ymax = est.sup), colour = "black", width = 0, size = 0.2) + 
                                          geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", size = 0.3) + 
                                          geom_point(size = 1.5, fill = "chocolate2", shape = 21)  +  
                                          geom_text(aes(y = 0, label = sign), nudge_y = 3.5, size = 5, color = "chocolate2") +
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
  plot_sh_15_a
  
  
  plot_mat_15_g <- gymno_results_mat_15 %>% dplyr::mutate(label.pos = max(est.sup, na.rm = TRUE)) %>%
                                            dplyr::mutate(label.pos = label.pos + 0.5*(max(est.sup, na.rm = T) - min(est.inf, na.rm = T))) %>%
                                            dplyr::mutate(mod = factor(models, levels = c("crown diameter", "crown ratio", "height 15", "max. height"))) %>%
                                            ggplot(aes(x = mod, y = est)) + ylim(-4, 4) +
                                            geom_errorbar(aes(ymin = est.inf, ymax = est.sup), colour = "black", width = 0, size = 0.2) + 
                                            geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", size = 0.3) + 
                                            geom_point(size = 1.5, fill = "forestgreen", shape = 22)  +  
                                            geom_text(aes(y = 0, label = sign), nudge_y = 3.5, size = 5, color = "forestgreen") +
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
  plot_mat_15_g
  
  
  plot_petp_15_g <- gymno_results_petp_15 %>% dplyr::mutate(label.pos = max(est.sup, na.rm = TRUE)) %>%
                                              dplyr::mutate(label.pos = label.pos + 0.5*(max(est.sup, na.rm = T) - min(est.inf, na.rm = T))) %>%
                                              dplyr::mutate(mod = factor(models, levels = c("crown diameter", "crown ratio", "height 15", "max. height"))) %>%
                                              ggplot(aes(x = mod, y = est)) + ylim(-4, 4) +
                                              geom_errorbar(aes(ymin = est.inf, ymax = est.sup), colour = "black", width = 0, size = 0.2) + 
                                              geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", size = 0.3) + 
                                              geom_point(size = 1.5, fill = "forestgreen", shape = 22) +  
                                              geom_text(aes(y = 0, label = sign), nudge_y = 3.5, size = 5, color = "forestgreen") +
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
  plot_petp_15_g
  
  
  plot_sla_15_g <- gymno_results_sla_15 %>% dplyr::mutate(label.pos = max(est.sup, na.rm = TRUE)) %>%
                                            dplyr::mutate(label.pos = label.pos + 0.5*(max(est.sup, na.rm = T) - min(est.inf, na.rm = T))) %>%
                                            dplyr::mutate(mod = factor(models, levels = c("crown diameter", "crown ratio", "height 15", "max. height"))) %>%
                                            ggplot(aes(x = mod, y = est)) + ylim(-4, 4) +
                                            geom_errorbar(aes(ymin = est.inf, ymax = est.sup), colour = "black", width = 0, size = 0.2) + 
                                            geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", size = 0.3) + 
                                            geom_point(size = 1.5, fill = "forestgreen", shape = 22) +  
                                            geom_text(aes(y = 0, label = sign), nudge_y = 3.5, size = 5, color = "forestgreen") +
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
  plot_sla_15_g
  
  
  plot_ssd_15_g <- gymno_results_ssd_15 %>% dplyr::mutate(label.pos = max(est.sup, na.rm = TRUE)) %>%
                                            dplyr::mutate(label.pos = label.pos + 0.5*(max(est.sup, na.rm = T) - min(est.inf, na.rm = T))) %>%
                                            dplyr::mutate(mod = factor(models, levels = c("crown diameter", "crown ratio", "height 15", "max. height"))) %>%
                                            ggplot(aes(x = mod, y = est)) + ylim(-4, 4) +
                                            geom_errorbar(aes(ymin = est.inf, ymax = est.sup), colour = "black", width = 0, size = 0.2) + 
                                            geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", size = 0.3) + 
                                            geom_point(size = 1.5, fill = "forestgreen", shape = 22) +
                                            geom_text(aes(y = 0, label = sign), nudge_y = 3.5, size = 5, color = "forestgreen") +
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
  plot_ssd_15_g
  
  
  plot_sh_15_g <- gymno_results_sh_15 %>% dplyr::mutate(label.pos = max(est.sup, na.rm = TRUE)) %>%
                                          dplyr::mutate(label.pos = label.pos + 0.5*(max(est.sup, na.rm = T) - min(est.inf, na.rm = T))) %>%
                                          dplyr::mutate(mod = factor(models, levels = c("crown diameter", "crown ratio", "height 15", "max. height"))) %>%
                                          ggplot(aes(x = mod, y = est)) + ylim(-4, 4) +
                                          geom_errorbar(aes(ymin = est.inf, ymax = est.sup), colour = "black", width = 0, size = 0.2) + 
                                          geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", size = 0.3) + 
                                          geom_point(size = 1.5, fill = "forestgreen", shape = 22) +
                                          geom_text(aes(y = 0, label = sign), nudge_y = 3.5, size = 5, color = "forestgreen") +
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
  plot_sh_15_g
  
  plot_angio_15_climate <- plot_grid(plot_mat_15_a, plot_petp_15_a, nrow = 1, ncol = 2)
  ggsave("manuscript/angiosperms_climate_nocomp_15_cv.png", plot_angio_15_climate, width = 9, height = 6, units = "cm", bg = "transparent") 
  plot_angio_15_traits <- plot_grid(plot_sla_15_a, plot_ssd_15_a, nrow = 1, ncol = 2)
  ggsave("manuscript/angiosperms_traits_nocomp_15_cv.png", plot_angio_15_traits, width = 9, height = 6, units = "cm", bg = "transparent") 
  plot_angio_15_sh <- plot_grid(plot_sh_15_a, nrow = 1, ncol = 1)
  ggsave("manuscript/angiosperms_sh_nocomp_15_cv.png", plot_angio_15_sh, width = 4.5, height = 6, units = "cm", bg = "transparent") 
  
  plot_gymno_15_climate <- plot_grid(plot_mat_15_g, plot_petp_15_g, nrow = 1, ncol = 2)
  ggsave("manuscript/gymnosperms_climate_nocomp_15_cv.png", plot_gymno_15_climate, width = 9, height = 6, units = "cm", bg = "transparent") 
  plot_gymno_15_traits <- plot_grid(plot_sla_15_g, plot_ssd_15_g, nrow = 1, ncol = 2)
  ggsave("manuscript/gymnosperms_traits_nocomp_15_cv.png", plot_gymno_15_traits, width = 9, height = 6, units = "cm", bg = "transparent") 
  plot_gymno_15_sh <- plot_grid(plot_sh_15_g, nrow = 1, ncol = 1)
  ggsave("manuscript/gymnosperms_sh_nocomp_15_cv.png", plot_gymno_15_sh, width = 4.5, height = 6, units = "cm", bg = "transparent") 
  
  return(angio_scaled_summary)
  
}


bivariate_scatterplots <- function() {
  
  ### Computing crown characteristics for each species based on the 100 values obtained from the resampling process
  # use species list previously defined
  
  sp_list <- read.csv(file = "output/mean_crown_characteristics_nocomp.csv")
  sp_list <- sp_list %>% dplyr::select(species)
  
  ## Height 
  height_asymptot <- read.csv(file = "data/height_asympt_nocomp_weighted_parameters.csv")
  height_nocomp_estimates <- height_asymptot %>% dplyr::select(-X) %>%
                                                 dplyr::filter(species %in% sp_list$species) %>%
                                                 dplyr::mutate(hmax = b1,
                                                               height_15 = 1.3 + b1 * (1-exp(-b2 * 15)) ^ b3,
                                                               height_30 = 1.3 + b1 * (1-exp(-b2 * 30)) ^ b3) %>%
                                                 dplyr::select(-b1, -b2, -b3) %>%
                                                 dplyr::group_by(species) %>%
                                                 dplyr::summarise(mean_hmax = mean(hmax),
                                                                  mean_h15 = mean(height_15)) %>%
                                                 dplyr::ungroup() %>%
                                                 dplyr::arrange(species)
  
  ### Diameter
  diameter_power <- read.csv(file = "data/diameter_nocomp_weighted_parameters.csv")
  diameter_nocomp_estimates <- diameter_power %>% dplyr::select(-X) %>%
                                                  dplyr::filter(species %in% sp_list$species) %>%
                                                  dplyr::mutate(diameter_15 = a1 * (15 ^ a2),
                                                                diameter_30 = a2 * (30 ^ a2)) %>%
                                                  dplyr::select(-a1, -a2) %>%
                                                  dplyr::group_by(species) %>%
                                                  dplyr::summarise(mean_d15 = mean(diameter_15)) %>%
                                                  dplyr::ungroup() %>%
                                                  dplyr::arrange(species) %>%
                                                  dplyr::select(-species)
  
  ## Ratio
  ratio_beta <- read.csv(file = "data/ratio_nocomp_weighted_parameters.csv")
  ratio_nocomp_estimates <- ratio_beta %>% dplyr::select(-X) %>%
                                           dplyr::filter(species %in% sp_list$species) %>%
                                           dplyr::mutate(ratio_15 = (exp(a1) + a2 * 15) / (1 + exp(a1) + a2 * 15),
                                                         ratio_30 = (exp(a1) + a2 * 30) / (1 + exp(a1) + a2 * 30)) %>%
                                           dplyr::select(-a1, -a2) %>%
                                           dplyr::group_by(species) %>%
                                           dplyr::summarise(mean_r15 = mean(ratio_15)) %>%
                                           dplyr::ungroup() %>%
                                           dplyr::arrange(species) %>%
                                           dplyr::select(-species)
  
  
  crown_characteristics <- cbind(height_nocomp_estimates, ratio_nocomp_estimates, diameter_nocomp_estimates)
  all_crown_characteristics <- na.omit(crown_characteristics)
  
  ### Scaling and centering crown characteristics at the functional group scale
  ## Summarizing crown characteristics at the species scale
  
  functional_groups <- read.csv(file = "data/functional_groups.csv")
  functional_groups <- functional_groups %>% dplyr::select(sp, new_gp)
  
  all_crown_characteristics <- left_join(all_crown_characteristics, functional_groups, by = c("species" = "sp"))
  
  ### Loading species' functional traits and climate niche
  ## Computing petp as map - etp
  traits_climate <- read.csv(file = "data/sp_traits_complete.csv")
  traits_climate <- traits_climate %>% dplyr::select(species, 
                                                     mat, map, etp,
                                                     sla, ssd, shade_tol_mean) %>%
                                       dplyr::mutate(petp = map - etp) %>%
                                       dplyr::select(-map, -etp)
  
  ## Creating a single file for angiosperms and another one for gymnosperms to run the models independently
  df <- left_join(all_crown_characteristics, traits_climate, by = "species")
  df <- na.omit(df)
  
  ## Creating a new column with reduced species names
  df <- df %>% dplyr::mutate(short_name = sapply(strsplit(as.character(species), " "), 
                              function(x) paste(substr(x, 1, 2), collapse = " ")))
  
  ## Separating angiosperms from gymnosperms
  angio <-df %>% dplyr::filter(new_gp == "angio")
  gymno <-df %>% dplyr::filter(new_gp == "gymno")
  
  ## Plotting pairs of variables
  # hmax 
  plot_1 <- ggplot(gymno, aes(x = mat, y = mean_hmax)) +
            geom_point(size = 1.5, fill = "forestgreen", shape = 22) +
            geom_text_repel(aes(label = short_name), 
                            col = "forestgreen", 
                            size = 5,
                            max.overlaps = 20) + 
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
                  axis.text.x = element_text(size = 8, angle = 360),
                  axis.text.y = element_text(size = 8, angle = 360))
  
  
  plot_2 <- ggplot(gymno, aes(x = petp, y = mean_hmax)) +
            geom_point(size = 1.5, fill = "forestgreen", shape = 22) +
            geom_text_repel(aes(label = short_name), 
                            col = "forestgreen", 
                            size = 5,
                            max.overlaps = 20) + 
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
                  axis.text.x = element_text(size = 8, angle = 360),
                  axis.text.y = element_text(size = 8, angle = 360))
  
  
  plot_3 <- ggplot(gymno, aes(x = sla, y = mean_hmax)) +
            geom_point(size = 1.5, fill = "forestgreen", shape = 22) +
            geom_text_repel(aes(label = short_name), 
                            col = "forestgreen", 
                            size = 5,
                            max.overlaps = 20) + 
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
                  axis.text.x = element_text(size = 8, angle = 360),
                  axis.text.y = element_text(size = 8, angle = 360))
  
  plot_4 <- ggplot(gymno, aes(x = ssd, y = mean_hmax)) +
            geom_point(size = 1.5, fill = "forestgreen", shape = 22) +
            geom_text_repel(aes(label = short_name), 
                            col = "forestgreen", 
                            size = 5,
                            max.overlaps = 20) + 
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
                  axis.text.x = element_text(size = 8, angle = 360),
                  axis.text.y = element_text(size = 8, angle = 360))
  
  plot_5 <- ggplot(gymno, aes(x = shade_tol_mean, y = mean_hmax)) +
            geom_point(size = 1.5, fill = "forestgreen", shape = 22) +
            geom_text_repel(aes(label = short_name), 
                            col = "forestgreen", 
                            size = 5,
                            max.overlaps = 20) + 
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
                  axis.text.x = element_text(size = 8, angle = 360),
                  axis.text.y = element_text(size = 8, angle = 360))
  
  plot_global <- plot_grid(plot_1, plot_2, plot_3, plot_4, plot_5, nrow = 2, ncol = 3)
  ggsave("manuscript/gymno_hmax_scatterplot.png", plot_global, width = 30, height = 20, units = "cm", bg = "transparent") 
  
  return(df)
  
}
