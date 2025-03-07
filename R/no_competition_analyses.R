### NO COMPETITION ANALYSES

compute_mean_crown_characteristics_nocomp <- function() {
  
  ## The asymptotic model was retained as fitting the best the height ~ dbh relationship based on AIC selection
  ## Note that weighted parameters have been selected based on RMSE values
  
  ### Computing mean crown characteristics for each species based on the 100 values obtained from the resampling process
  ## Height 
  height_asymptot <- read.csv(file = "data/height_asympt_nocomp_weighted_parameters.csv")
  height_nocomp_estimates <- height_asymptot %>% dplyr::select(-X) %>%
                                                 dplyr::group_by(species) %>%
                                                 dplyr::summarise(b1_mean = mean(b1),
                                                                  b2_mean = mean(b2),
                                                                  b3_mean = mean(b3)) %>%
                                                 dplyr::mutate(hmax = b1_mean,
                                                               height_15 = 1.3 + b1_mean * (1-exp(-b2_mean * 15)) ^ b3_mean,
                                                               height_30 = 1.3 + b1_mean * (1-exp(-b2_mean * 30)) ^ b3_mean) %>%
                                                 dplyr::select(-b1_mean, -b2_mean, -b3_mean) %>%
                                                 dplyr::ungroup()
  
  ### Diameter
  diameter_power <- read.csv(file = "data/diameter_nocomp_weighted_parameters.csv")
  diameter_nocomp_estimates <- diameter_power %>% dplyr::select(-X) %>%
                                                  dplyr::group_by(species) %>%
                                                  dplyr::summarise(a1_mean = mean(a1),
                                                                   a2_mean = mean(a2)) %>%
                                                  dplyr::mutate(diameter_15 = a1_mean * (15 ^ a2_mean),
                                                                diameter_30 = a2_mean * (30 ^ a2_mean)) %>%
                                                  dplyr::select(-a1_mean, -a2_mean) %>%
                                                  dplyr::ungroup()
  
  ## Ratio
  ratio_beta <- read.csv(file = "data/ratio_nocomp_weighted_parameters.csv")
  ratio_nocomp_estimates <- ratio_beta %>% dplyr::select(-X) %>%
                                           dplyr::group_by(species) %>%
                                           dplyr::summarise(a1_mean = mean(a1),
                                                            a2_mean = mean(a2)) %>%
                                           dplyr::mutate(ratio_15 = (exp(a1_mean) + a2_mean * 15) / (1 + exp(a1_mean) + a2_mean * 15),
                                                         ratio_30 = (exp(a1_mean) + a2_mean * 30) / (1 + exp(a1_mean) + a2_mean * 30)) %>%
                                           dplyr::select(-a1_mean, -a2_mean) %>%
                                           dplyr::ungroup()
  
  ### Combine all crown characteristics so that only species for which all three metrics are available remain
  # use the characteristic with the highest number of species as the primary one
  
  crown_characteristics <- join_all(list(height_nocomp_estimates, ratio_nocomp_estimates, diameter_nocomp_estimates), 
                                    by = "species", type = "left")
  crown_characteristics <- na.omit(crown_characteristics)
  
  
  ### Computing depth as height x ratio in order to later compute volume
  mean_crown_characteristics <- crown_characteristics %>% dplyr::mutate(depth_15 = height_15 * ratio_15,
                                                                        depth_30 = height_30 * ratio_30)
  
  ### Computing volume using distinct equations according to the species' functional group
  functional_groups <- read.csv(file = "data/functional_groups.csv")
  functional_groups <- functional_groups %>% dplyr::select(sp, new_gp)
  
  mean_crown_characteristics <- left_join(mean_crown_characteristics, functional_groups, by = c("species" = "sp"))
  
  angiosperms <- mean_crown_characteristics %>% dplyr::filter(new_gp == "angio") %>%
                                                dplyr::mutate(volume_15 = 4/3 * pi * (diameter_15/2) * (diameter_15/2) * (depth_15/2),
                                                              volume_30 = 4/3 * pi * (diameter_30/2) * (diameter_30/2) * (depth_30/2))
  
  gymnosperms <-  mean_crown_characteristics %>% dplyr::filter(new_gp == "gymno") %>%
                                                 dplyr::mutate(volume_15 = 1/2 * pi * (diameter_15/2)*(diameter_15/2) * (depth_15),
                                                               volume_30 = 1/2 * pi * (diameter_30/2)*(diameter_30/2) * (depth_30))
  
  mean_characteristics <- rbind(angiosperms, gymnosperms) 
  
  ### Shortening species names to improve readability of the figures that will later be computed
  divided_names <- str_split(mean_characteristics$species, " ", simplify = TRUE)  
  short_names <- str_c(str_to_title(str_sub(divided_names[,1], 1, 2)), " ", str_to_lower(str_sub(divided_names[,2], 1, 3)))  

  mean_characteristics$short_names <- short_names
  
  write.csv(file = "output/mean_crown_characteristics_nocomp.csv", mean_characteristics)
  
  return(mean_characteristics)
  
}


pca_ecological_strategies <- function(mean_characteristics_nocomp) {
  
  ## running PCA for both functional groups combined
  rownames(mean_characteristics_nocomp) <- mean_characteristics_nocomp$short_names
  
  # keep only studied crown characteristics
  combined_groups <- mean_characteristics_nocomp %>% dplyr::select(hmax, height_15, ratio_15, diameter_15, volume_15)
  
  # run pca and centered and scaled data
  res_pca_combined <- prcomp(combined_groups, scale = TRUE, center = TRUE)
  
  # save pca figure without arrow names for aesthetic purposes
  q <- autoplot(res_pca_combined, 
                scale = 0,
                data = combined_groups, 
                label = FALSE,            
                loadings = TRUE,         
                loadings.label = FALSE,  
                loadings.label.size = 3,
                loadings.colour = "darkgrey") +
       theme(rect = element_rect(fill = "transparent"), 
             panel.background = element_rect(fill = "transparent"), 
             panel.grid = element_blank(), 
             plot.background = element_rect(fill = "transparent", color = NA),
             strip.background = element_blank()) +
    geom_segment(data = as.data.frame(res_pca_combined$rotation), 
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 col = "darkgrey") +
    geom_point(size = 3, 
               fill = c("chocolate2", "darkgreen")[unclass(factor(mean_characteristics_nocomp$new_gp))], 
               shape = 21) +
    geom_text_repel(aes(label = rownames(combined_groups)), 
                    col = c("chocolate2", "darkgreen")[unclass(factor(mean_characteristics_nocomp$new_gp))], 
                    size = 6,
                    max.overlaps = 20)
  
  ggsave("manuscript/pca_allspecies_nonames.png", q, width = 30, height = 30, units = "cm", bg = "transparent") 
  
  
  ## running PCA for angiosperms
  angiosperms <- mean_characteristics_nocomp %>% dplyr::filter(new_gp == "angio")
  
  # keep only studied crown characteristics
  angiosperms <- angiosperms %>% dplyr::select(hmax, height_15, ratio_15, diameter_15, volume_15)
  
  # run pca and centered and scaled data
  res_pca_angio <- prcomp(angiosperms, scale = TRUE, center = TRUE)
  
  # save pca figure without arrow names for aesthetic purposes
  s <- autoplot(res_pca_angio, 
                data = angiosperms, 
                scale = 0,
                label = FALSE,            
                loadings = TRUE,         
                loadings.label = FALSE,  
                loadings.label.size = 3,
                loadings.colour = "chocolate2") +
    theme(rect = element_rect(fill = "transparent"), 
          panel.background = element_rect(fill = "transparent"), 
          panel.grid = element_blank(), 
          plot.background = element_rect(fill = "transparent", color = NA),
          strip.background = element_blank()) +
    geom_point(size = 3, 
               fill = "chocolate2", 
               shape = 21) +
    geom_text_repel(aes(label = rownames(angiosperms)), 
                    size = 6, 
                    colour = "chocolate2",
                    max.overlaps = 20)
  
  ggsave("manuscript/pca_angiosperms_nonames.png", s, width = 30, height = 30, units = "cm", bg = "transparent") 
  
 
  
  ## running PCA for gymnosperms
  gymnosperms <- mean_characteristics_nocomp %>% dplyr::filter(new_gp == "gymno")
  
  # keep only studied crown characteristics
  gymnosperms <- gymnosperms %>% dplyr::select(hmax, height_15, ratio_15, diameter_15, volume_15)
  
  # run pca and centered and scaled data
  res_pca_gymno <- prcomp(gymnosperms, scale = TRUE, center = TRUE)
  
  # save pca figure without arrow names for esthetic purposes
  u <- autoplot(res_pca_gymno, 
                data = gymnosperms, 
                scale = 0,
                label = FALSE,            
                loadings = TRUE,         
                loadings.label = FALSE,  
                loadings.label.size = 3,
                loadings.colour = "forestgreen") +
    theme(rect = element_rect(fill = "transparent"), 
          panel.background = element_rect(fill = "transparent"), 
          panel.grid = element_blank(), 
          plot.background = element_rect(fill = "transparent", color = NA),
          strip.background = element_blank()) +
    geom_point(size = 3, 
               fill = "forestgreen", 
               shape = 21) +
    geom_text_repel(aes(label = rownames(gymnosperms)), 
                    size = 5, 
                    colour = "forestgreen",
                    max.overlaps = 20)
  
  ggsave("manuscript/pca_gymnosperms_nonames.png", u, width = 30, height = 30, units = "cm", bg = "transparent") 

  return(res_pca_gymno)

}



traits_climate_crown_model_selection_angiosperms <- function(mean_characteristics_nocomp) {
  
  ### Computing crown characteristics for each species based on the 100 values obtained from the resampling process
  # use species list previously defined
  
  sp_list <- mean_characteristics_nocomp
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
                                            
  
  
  ### Loading species' functional traits and climate niche
  ## Computing petp as map - etp
  traits_climate <- read.csv(file = "data/sp_traits_complete.csv")
  traits_climate <- traits_climate %>% dplyr::select(species, 
                                                     mat, map, etp,
                                                     sla, ssd, shade_tol_mean) %>%
                                       dplyr::mutate(petp = map - etp) %>%
                                       dplyr::select(-map, -etp)
  
  
  ## Creating a single file for angiosperms and another one for gymnosperms to run the models independently
  angiosperms <- left_join(angiosperms, traits_climate, by = "species")
  
  ## Scaling and centering species functional traits and climate niche at the functional group scale
  angiosperms <- angiosperms %>% dplyr::mutate(mat = scale(mat, center = TRUE, scale = TRUE)[, 1],
                                               petp = scale(petp, center = TRUE, scale = TRUE)[, 1],
                                               sla = scale(sla, center = TRUE, scale = TRUE)[, 1],
                                               ssd = scale(ssd, center = TRUE, scale = TRUE)[, 1],
                                               sh = scale(shade_tol_mean, center = TRUE, scale = TRUE)[, 1]) %>%
                                 dplyr::select(-shade_tol_mean)
  
  angiosperms <- na.omit(angiosperms)
  
  
  ## 1 Running linear models to explore the relationships between crown characteristics and functional traits / climate niche
  ## Performing stepAIC selection on (1) climate niche first, then (2) functional traits, to obtain the full and best model
  ## Running models independently for angiosperms and gymnosperms
  
  options(na.action = "na.fail")
  
  ## hmax
  model_hmax_climate <- lm(hmax ~ mat + petp, weights = (1/cv_hmax), data = angiosperms) 
  hmax_climate <- dredge(model_hmax_climate) 
  best_hmax_climate <- get.models(hmax_climate, 1)[[1]]
  summary_hmax_climate <- as.data.frame(hmax_climate)
  write.csv(file = "output/modelselection_hmax_climate_angio.csv", summary_hmax_climate)
  
  
  model_hmax_traits <- lm(hmax ~ sla + ssd + sh, weights = (1/cv_hmax), data = angiosperms)
  hmax_traits <- dredge(model_hmax_traits) 
  best_hmax_traits <- get.models(hmax_traits, 1)[[1]]
  summary_hmax_traits <- as.data.frame(hmax_traits)
  write.csv(file = "output/modelselection_hmax_traits_angio.csv", summary_hmax_traits)
  
  
  
  ## height at 15 cm DBH
  model_h15_climate <- lm(height_15 ~ mat + petp, weights = (1/cv_height_15), data = angiosperms) 
  h15_climate <- dredge(model_h15_climate) 
  best_h15_climate <- get.models(h15_climate, 1)[[1]]
  summary_h15_climate <- as.data.frame(h15_climate)
  write.csv(file = "output/modelselection_h15_climate_angio.csv", summary_h15_climate)
  
  
  model_h15_traits <- lm(height_15 ~ sla + ssd + sh, weights = (1/cv_height_15), data = angiosperms) 
  h15_traits <- dredge(model_h15_traits) 
  best_h15_traits <- get.models(h15_traits, 1)[[1]]
  summary_h15_traits <- as.data.frame(h15_traits)
  write.csv(file = "output/modelselection_h15_traits_angio.csv", summary_h15_traits)
  
  
  ## ratio_15
  model_r15_climate <- lm(ratio_15 ~ mat + petp, weights = (1/cv_ratio_15), data = angiosperms) 
  r15_climate <- dredge(model_r15_climate) 
  best_r15_climate <- get.models(r15_climate, 1)[[1]]
  summary_r15_climate <- as.data.frame(r15_climate)
  write.csv(file = "output/modelselection_r15_climate_angio.csv", summary_r15_climate)
  
  
  model_r15_traits <- lm(ratio_15 ~ sla + ssd + sh, weights = (1/cv_ratio_15), data = angiosperms) 
  r15_traits <- dredge(model_r15_traits) 
  best_r15_traits <- get.models(r15_traits, 1)[[1]]
  summary_r15_traits <- as.data.frame(r15_traits)
  write.csv(file = "output/modelselection_r15_traits_angio.csv", summary_r15_traits)
  
  
  ## diameter_15
  model_d15_climate <- lm(diameter_15 ~ mat + petp, weights = (1/cv_diameter_15), data = angiosperms) 
  d15_climate <- dredge(model_d15_climate) 
  best_d15_climate <- get.models(d15_climate, 1)[[1]]
  summary_d15_climate <- as.data.frame(d15_climate)
  write.csv(file = "output/modelselection_d15_climate_angio.csv", summary_d15_climate)
  
  
  model_d15_traits <- lm(diameter_15 ~ sla + ssd + sh, weights = (1/cv_diameter_15), data = angiosperms) 
  d15_traits <- dredge(model_d15_traits) 
  best_d15_traits <- get.models(d15_traits, 1)[[1]]
  summary_d15_traits <- as.data.frame(d15_traits)
  write.csv(file = "output/modelselection_d15_traits_angio.csv", summary_d15_traits)
  
  return(angiosperms)
  
}



traits_climate_crown_model_selection_gymnosperms <- function(mean_characteristics_nocomp) {
  
  ### Computing crown characteristics for each species based on the 100 values obtained from the resampling process
  # use species list previously defined
  
  sp_list <- mean_characteristics_nocomp
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
  traits_climate <- traits_climate %>% dplyr::select(species, 
                                                     mat, map, etp,
                                                     sla, ssd, shade_tol_mean) %>%
                                       dplyr::mutate(petp = map - etp) %>%
                                       dplyr::select(-map, -etp)
  
  ## Creating a single file for angiosperms and another one for gymnosperms to run the models independently
  gymnosperms <- left_join(gymnosperms, traits_climate, by = "species")
  
  gymnosperms <- gymnosperms %>% dplyr::mutate(mat = scale(mat, center = TRUE, scale = TRUE)[, 1],
                                               petp = scale(petp, center = TRUE, scale = TRUE)[, 1],
                                               sla = scale(sla, center = TRUE, scale = TRUE)[, 1],
                                               ssd = scale(ssd, center = TRUE, scale = TRUE)[, 1],
                                               sh = scale(shade_tol_mean, center = TRUE, scale = TRUE)[, 1]) %>%
                                 dplyr::select(-shade_tol_mean)
  
  gymnosperms <- na.omit(gymnosperms)
  
  
  ### Running linear models to explore the relationships between crown characteristics and functional traits / climate niche
  ## Performing stepAIC selection on (1) climate niche first, then (2) functional traits, to obtain the full and best model
  ## Running models independently for angiosperms and gymnosperms
  
  options(na.action = "na.fail")
  
  ## hmax
  model_hmax_climate <- lm(hmax ~ mat + petp, weights = (1/cv_hmax), data = gymnosperms) 
  hmax_climate <- dredge(model_hmax_climate) 
  best_hmax_climate <- get.models(hmax_climate, 1)[[1]]
  summary_hmax_climate <- as.data.frame(hmax_climate)
  write.csv(file = "output/modelselection_hmax_climate_gymno.csv", summary_hmax_climate)
  
  
  model_hmax_traits <- lm(hmax ~ sla + ssd + sh, weights = (1/cv_hmax), data = gymnosperms)
  hmax_traits <- dredge(model_hmax_traits) 
  best_hmax_traits <- get.models(hmax_traits, 1)[[1]]
  summary_hmax_traits <- as.data.frame(hmax_traits)
  write.csv(file = "output/modelselection_hmax_traits_gymno.csv", summary_hmax_traits)
  
  
  
  ## height at 15 cm DBH
  model_h15_climate <- lm(height_15 ~ mat + petp, weights = (1/cv_height_15), data = gymnosperms) 
  h15_climate <- dredge(model_h15_climate) 
  best_h15_climate <- get.models(h15_climate, 1)[[1]]
  summary_h15_climate <- as.data.frame(h15_climate)
  write.csv(file = "output/modelselection_h15_climate_gymno.csv", summary_h15_climate)
  
  
  model_h15_traits <- lm(height_15 ~ sla + ssd + sh, weights = (1/cv_height_15), data = gymnosperms) 
  h15_traits <- dredge(model_h15_traits) 
  best_h15_traits <- get.models(h15_traits, 1)[[1]]
  summary_h15_traits <- as.data.frame(h15_traits)
  write.csv(file = "output/modelselection_h15_traits_gymno.csv", summary_h15_traits)
  
  
  ## ratio_15
  model_r15_climate <- lm(ratio_15 ~ mat + petp, weights = (1/cv_ratio_15), data = gymnosperms) 
  r15_climate <- dredge(model_r15_climate) 
  best_r15_climate <- get.models(r15_climate, 1)[[1]]
  summary_r15_climate <- as.data.frame(r15_climate)
  write.csv(file = "output/modelselection_r15_climate_gymno.csv", summary_r15_climate)
  
  
  model_r15_traits <- lm(ratio_15 ~ sla + ssd + sh, weights = (1/cv_ratio_15), data = gymnosperms) 
  r15_traits <- dredge(model_r15_traits) 
  best_r15_traits <- get.models(r15_traits, 1)[[1]]
  summary_r15_traits <- as.data.frame(r15_traits)
  write.csv(file = "output/modelselection_r15_traits_gymno.csv", summary_r15_traits)
  
  
  ## diameter_15
  model_d15_climate <- lm(diameter_15 ~ mat + petp, weights = (1/cv_diameter_15), data = gymnosperms) 
  d15_climate <- dredge(model_d15_climate) 
  best_d15_climate <- get.models(d15_climate, 1)[[1]]
  summary_d15_climate <- as.data.frame(d15_climate)
  write.csv(file = "output/modelselection_d15_climate_gymno.csv", summary_d15_climate)
  
  
  model_d15_traits <- lm(diameter_15 ~ sla + ssd + sh, weights = (1/cv_diameter_15), data = gymnosperms) 
  d15_traits <- dredge(model_d15_traits) 
  best_d15_traits <- get.models(d15_traits, 1)[[1]]
  summary_d15_traits <- as.data.frame(d15_traits)
  write.csv(file = "output/modelselection_d15_traits_gymno.csv", summary_d15_traits)
  
  return(gymnosperms)
  
}


traits_climate_crown_model_running_angiosperms <- function(angiosperms_scaled) {
  
  angiosperms <- angiosperms_scaled
  
  ### Model have been selected based on the results obtained with the dredge function
  ## See exported files in the output folder

  # hmax ~ mat + cwd + sla
  hmax <- lm(hmax ~ mat + petp + sla, weights = (1/cv_hmax), data = angiosperms)
  
  # h15 ~ mat + cwd + sla
  h15 <- lm(height_15 ~ mat + petp + sla, weights = (1/cv_height_15), data = angiosperms)
  
  # r15 ~ sh
  r15 <- lm(ratio_15 ~ sla + sh, weights = (1/cv_ratio_15), data = angiosperms)
  
  # d15 ~ petp + ssd + sh
  d15 <- lm(diameter_15 ~ petp + ssd + sh, weights = (1/cv_diameter_15), data = angiosperms)
  
  
  ## Storing models' results
  models <- c("max. height", "height 15", "crown ratio", "crown diameter")
  angio_models <- list(hmax, h15, r15, d15)
  
  # Creating storage files 
  angio_mat <- matrix(nrow = length(models), ncol = 6)
  colnames(angio_mat) <- c("mod", "est.inf", "est.sup", "est", "p.value", "sign")

  angio_petp <- matrix(nrow = length(models), ncol = 6)
  colnames(angio_petp) <- c("mod", "est.inf", "est.sup", "est", "p.value", "sign")
  
  angio_sla <- matrix(nrow = length(models), ncol = 6)
  colnames(angio_sla) <- c("mod", "est.inf", "est.sup", "est", "p.value", "sign")
  
  angio_ssd <- matrix(nrow = length(models), ncol = 6)
  colnames(angio_ssd) <- c("mod", "est.inf", "est.sup", "est", "p.value", "sign")
  
  angio_sh <- matrix(nrow = length(models), ncol = 6)
  colnames(angio_sh) <- c("mod", "est.inf", "est.sup", "est", "p.value", "sign")
  
  angio_results <- list(angio_mat, angio_petp, angio_sla, angio_ssd, angio_sh)
  variables <- c("mat", "petp", "sla", "ssd", "sh")
  
  
  for (i in 1:length(angio_results)) {
    
    for (j in 1:length(angio_models)) {

    angio_results[[i]][j,"mod"] <- models[j]
    
    tryCatch({
      
    model <- angio_models[[j]]
    
    angio_results[[i]][j, "est.inf"] = confint.lm(model)[variables[i], "2.5 %"]
    angio_results[[i]][j, "est.sup"] = confint.lm(model)[variables[i], "97.5 %"]
    angio_results[[i]][j, "est"] = summary(model)$coefficients[variables[i], "Estimate"]
    angio_results[[i]][j, "p.value"] = summary(model)$coefficients[variables[i], "Pr(>|t|)"]
    
    },
    error = function(e) {
      print("not included in the selected model")
    })
    
    }
    
  }
  
  
  # Creating a function to convert a matrix to a data frame, mutate, and convert back to a matrix
  # as well as transform p.values into levels of significance ***
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
  final_angio_results <- lapply(angio_results, mutate_matrix)
  
  
  ### Representing the results using an effect size plot
  ## Starting with angiosperms
  variables <- c("mat", "petp", "sla", "ssd", "sh")
  
  for (i in 1:length(final_angio_results)) {
    
    df <- as.data.frame(final_angio_results[[i]])
    
    plot <- df %>% dplyr::mutate(est.inf = as.numeric(est.inf),
                                 est.sup = as.numeric(est.sup),
                                 est = as.numeric(est)) %>%
      
                   dplyr::mutate(label.pos = max(est.sup, na.rm = TRUE)) %>%
                   dplyr::mutate(label.pos = label.pos + 0.5*(max(est.sup, na.rm = TRUE) - min(est.inf, na.rm = TRUE))) %>%
                   dplyr::mutate(mod = factor(models, levels = c("crown diameter", "crown ratio", "height 15", "max. height"))) %>%
                  
                   ggplot(aes(x = mod, y = est)) + ylim(-2.5, 2.5) +
                   geom_errorbar(aes(ymin = est.inf, ymax = est.sup), colour = "black", width = 0, size = 0.5) + 
                   geom_point(size = 2.5, fill = "chocolate2", shape = 21) +
                   geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed", size = 0.6) +  
                   geom_text(aes(y = 0, label = sign), nudge_y = 2, size = 8, color = "chocolate2") +
                  
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
    
    ggsave(paste0("manuscript/angiosperms_", variables[i], ".png"), plot, width = 4, height = 6, units = "cm", bg = "transparent") 
    
  }
  
  return(final_angio_results)
  
}
  
  
  

traits_climate_crown_model_running_gymnosperms <- function(gymnosperms_scaled) {
  
  gymnosperms <- gymnosperms_scaled
  
  ### Model have been selected based on the results obtained with the dredge function
  ## See exported files in the output folder
  
  # hmax ~ 1
  hmax <- lm(hmax ~ 1, weights = (1/cv_hmax), data = gymnosperms)
  
  # h15 ~ cwd + sh
  h15 <- lm(height_15 ~ petp + ssd + sh, weights = (1/cv_height_15), data = gymnosperms)
  
  # r15 ~ mat + sh
  r15 <- lm(ratio_15 ~ mat + ssd + sh, weights = (1/cv_ratio_15), data = gymnosperms)
  
  # d15 ~ mat
  d15 <- lm(diameter_15 ~ 1, weights = (1/cv_diameter_15), data = gymnosperms)
  
  
  ## Storing models' results
  models <- c("max. height", "height 15", "crown ratio", "crown diameter")
  gymno_models <- list(hmax, h15, r15, d15)
  
  # Creating storage files 
  gymno_mat <- matrix(nrow = length(models), ncol = 6)
  colnames(gymno_mat) <- c("mod", "est.inf", "est.sup", "est", "p.value", "sign")
  
  gymno_petp <- matrix(nrow = length(models), ncol = 6)
  colnames(gymno_petp) <- c("mod", "est.inf", "est.sup", "est", "p.value", "sign")
  
  gymno_sla <- matrix(nrow = length(models), ncol = 6)
  colnames(gymno_sla) <- c("mod", "est.inf", "est.sup", "est", "p.value", "sign")
  
  gymno_ssd <- matrix(nrow = length(models), ncol = 6)
  colnames(gymno_ssd) <- c("mod", "est.inf", "est.sup", "est", "p.value", "sign")
  
  gymno_sh <- matrix(nrow = length(models), ncol = 6)
  colnames(gymno_sh) <- c("mod", "est.inf", "est.sup", "est", "p.value", "sign")
  
  gymno_results <- list(gymno_mat, gymno_petp, gymno_sla, gymno_ssd, gymno_sh)
  variables <- c("mat", "petp", "sla", "ssd", "sh")
  
  
  for (i in 1:length(gymno_results)) {
    
    for (j in 1:length(gymno_models)) {
      
      gymno_results[[i]][j,"mod"] <- models[j]
      
      tryCatch({
        
        model <- gymno_models[[j]]
        
        gymno_results[[i]][j, "est.inf"] = confint.lm(model)[variables[i], "2.5 %"]
        gymno_results[[i]][j, "est.sup"] = confint.lm(model)[variables[i], "97.5 %"]
        gymno_results[[i]][j, "est"] = summary(model)$coefficients[variables[i], "Estimate"]
        gymno_results[[i]][j, "p.value"] = summary(model)$coefficients[variables[i], "Pr(>|t|)"]
        
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
  final_gymno_results <- lapply(gymno_results, mutate_matrix)
  
  
  ### Representing the results using an effect size plot
  ## Starting with gymnosperms
  variables <- c("mat", "petp", "sla", "ssd", "sh")
  
  for (i in 1:length(final_gymno_results)) {
    
    df <- as.data.frame(final_gymno_results[[i]])
    
    plot <- df %>% dplyr::mutate(est.inf = as.numeric(est.inf),
                                 est.sup = as.numeric(est.sup),
                                 est = as.numeric(est)) %>%
      
                   dplyr::mutate(label.pos = max(est.sup, na.rm = TRUE)) %>%
                   dplyr::mutate(label.pos = label.pos + 0.5*(max(est.sup, na.rm = TRUE) - min(est.inf, na.rm = TRUE))) %>%
                   dplyr::mutate(mod = factor(models, levels = c("crown diameter", "crown ratio", "height 15", "max. height"))) %>%
                    
                   ggplot(aes(x = mod, y = est)) + ylim(-2.5, 2.5) +
                   geom_errorbar(aes(ymin = est.inf, ymax = est.sup), colour = "black", width = 0, size = 0.5) + 
                   geom_point(size = 2.5, fill = "forestgreen", shape = 22) +
                   geom_hline(yintercept = 0, color = "lightgrey", linetype = "dashed", size = 0.6) +  
                   geom_text(aes(y = 0, label = sign), nudge_y = 2, size = 8, color = "forestgreen") +
                    
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
    
    ggsave(paste0("manuscript/gymnosperms_", variables[i], ".png"), plot, width = 4, height = 6, units = "cm", bg = "transparent") 
    
  }
  
  return(final_gymno_results)
  

}


