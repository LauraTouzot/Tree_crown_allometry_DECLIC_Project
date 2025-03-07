# Merge crown data

merge_crown_data <- function(baad_crown_complete, FHM_data, FIA_data, Spain_NFI, FUNDIV_data, FUNDIV_crown,
                             MONTANE_data, data_paper_crown, evan_crown_data,
                             dettmann_data, heym_2017, gen_tree_data, legacytree,
                             ICP_data, anderson2015, dalponte2016, FrenchNFI_data, 
                             Quebec, Canada_data, Usoltsev, Sullivan_data, tallo_database) {
  
  baad_crown_complete$data <- "baad"
  baad_crown_complete$ref <- baad_crown_complete$studyName
  baad <- baad_crown_complete[!is.na(baad_crown_complete$DBH_cm), c("sp", "DBH_cm", "C_diam_m", "CR", "C_depth_m", "HT_m","W", "data", "ref", 
                                                                    "latitude_plot", "longitude_plot", "latitude_tree", "longitude_tree", "location_ID", "continent")]
  
  # crown depth data are estimated from tree total height and crown ratio
  FHM_data$data <- "FHM"
  FHM_data$ref <- "FHM"
  FHM <- FHM_data[!duplicated(FHM_data$TreeID) & !is.na(FHM_data$DBH_cm), c("sp", "DBH_cm", "C_diam_m", "CR", "C_depth_m", "HT_m", "W", "data", "ref", 
                                                                            "latitude_plot", "longitude_plot", "latitude_tree", "longitude_tree", "location_ID", "continent")]

  # crown depth data are estimated from tree total height and crown ratio  
  FIA_data$data <- "FIA"
  FIA_data$ref <- "FIA"
  FIA_data$C_diam_m <- NA
  FIA <- FIA_data[!duplicated(FIA_data$TreeID) & !is.na(FIA_data$DBH_cm), c("sp", "DBH_cm", "C_diam_m", "CR", "C_depth_m", "HT_m", "W", "data", "ref", 
                                                                            "latitude_plot", "longitude_plot", "latitude_tree", "longitude_tree", "location_ID", "continent")]
  
  Spain_NFI$data <- "Spain"
  Spain_NFI$ref <- "Spain"
  Spain <- Spain_NFI[!is.na(Spain_NFI$DBH_cm), c("sp", "DBH_cm", "C_diam_m", "CR", "C_depth_m", "HT_m", "W", "data", "ref", 
                                                 "latitude_plot", "longitude_plot", "latitude_tree", "longitude_tree", "location_ID", "continent")]
  
  FUNDIV_data$data <- "FUNDIV"
  FUNDIV_data$ref <- FUNDIV_data$country
  FUNDIV_data$C_diam_m <- NA
  FUNDIV_data$C_depth_m <- NA
  FUNDIV_data$CR <- NA
  FUNDIV <- FUNDIV_data[!is.na(FUNDIV_data$DBH_cm), c("sp", "DBH_cm", "C_diam_m", "CR", "C_depth_m", "HT_m", "W", "data", "ref", 
                                                      "latitude_plot", "longitude_plot", "latitude_tree", "longitude_tree", "location_ID", "continent")]
  
  FUNDIV_crown$data <- "FUNDIV_crown"
  FUNDIV_crown$ref <- "Jucker2015"
  FUNDIV_CROWN <- FUNDIV_crown[!is.na(FUNDIV_crown$DBH_cm), c("sp", "DBH_cm", "C_diam_m", "CR", "C_depth_m", "HT_m", "W", "data", "ref", 
                                                              "latitude_plot", "longitude_plot", "latitude_tree", "longitude_tree", "location_ID", "continent")]
  
  MONTANE_data$data <- "MONTANE"
  MONTANE_data$ref <- "MONTANE"
  MONTANE <- MONTANE_data[, c("sp", "DBH_cm", "C_diam_m", "CR", "C_depth_m", "HT_m", "W", "data", "ref", 
                              "latitude_plot", "longitude_plot", "latitude_tree", "longitude_tree", "location_ID", "continent")]
  
  data_paper_crown$data <- "Fuhr"
  data_paper_crown$ref <- "Fuhr2017"
  data_paper <- data_paper_crown[!is.na(data_paper_crown$DBH_cm), c("sp", "DBH_cm", "C_diam_m", "CR", "C_depth_m", "HT_m", "W", "data", "ref", 
                                                                    "lat", "long", "y", "x", "location_ID", "continent")]
  names(data_paper) <- c("sp", "DBH_cm", "C_diam_m", "CR", "C_depth_m", "HT_m", "W", "data", "ref", 
                         "latitude_plot", "longitude_plot", "latitude_tree", "longitude_tree", "location_ID", "continent")
  
  evan_crown_data$data <- "Evan"
  evan_crown_data$ref <- "Evan2015"
  evan_crown <- evan_crown_data[!is.na(evan_crown_data$DBH_cm), c("sp", "DBH_cm", "C_diam_m", "CR", "C_depth_m", "HT_m", "W", "data", "ref", 
                                                                  "latitude_plot", "longitude_plot", "latitude_tree", "longitude_tree", "location_ID", "continent")]
  
  dettmann_data$data <- "Dettmann"
  dettmann_data$ref <- "Dettmann2019"
  dettmann_data$HT_m <- NA
  dettmann_data$C_depth_m <- NA
  dettmann <- dettmann_data[!is.na(dettmann_data$DBH_cm), c("sp", "DBH_cm", "C_diam_m", "CR", "C_depth_m", "HT_m", "W", "data", "ref", 
                                                            "latitude_plot", "longitude_plot", "latitude_tree", "longitude_tree", "location_ID", "continent")]
  
  heym_2017$data <- "Heym"
  heym_2017$ref <- "Heym2017"
  heym <- heym_2017[!is.na(heym_2017$DBH_cm), c("sp", "DBH_cm", "C_diam_m", "CR", "C_depth_m", "HT_m", "W", "data", "ref", 
                                                "latitude_plot", "longitude_plot", "latitude_tree", "longitude_tree", "location_ID", "continent")]
  
  gen_tree_data$data <- "GenTree"
  gen_tree_data$ref <-"Opgenoorth2021"
  gen_tree_data$CR <- NA
  gen_tree_data$C_depth_m <- NA
  gentree <- gen_tree_data[!is.na(gen_tree_data$DBH_cm), c("sp", "DBH_cm", "C_diam_m", "CR", "C_depth_m", "HT_m", "W", "data", "ref", 
                                                           "latitude_plot", "longitude_plot", "latitude_tree", "longitude_tree", "location_ID", "continent")]
  
  legacytree$data <- "LegacyTree"
  legacytree$ref <- legacytree$AUTHOR
  legacytree_data <- legacytree[!is.na(legacytree$DBH_cm), c("sp", "DBH_cm", "C_diam_m", "CR", "C_depth_m", "HT_m", "data", "ref", 
                                                             "latitude_plot", "longitude_plot", "latitude_tree", "longitude_tree", "location_ID", "continent")]
  
  ICP_data$data <- "ICP"
  ICP_data$ref <- "ICP"
  ICP_data_s <- ICP_data[!is.na(ICP_data$DBH_cm), c("sp", "DBH_cm", "C_diam_m", "CR", "C_depth_m", "HT_m", "W", "data", "ref", 
                                                    "latitude_plot", "longitude_plot", "latitude_tree", "longitude_tree", "location_ID", "continent")]
  
  anderson2015$data <- "Anderson"
  anderson2015$ref <- "Anderson2015"
  anderson <- anderson2015[!is.na(anderson2015$DBH_cm), c("sp", "DBH_cm", "C_diam_m", "CR", "C_depth_m", "HT_m", "W", "data", "ref", 
                                                          "latitude_plot", "longitude_plot", "latitude_tree", "longitude_tree", "location_ID", "continent")]
  
  dalponte2016$data <- "Dalponte"
  dalponte2016$ref <- "Dalponte2016"
  dalponte2016$C_depth_m <- NA
  dalponte2016$CR <- NA
  dalponte <- dalponte2016[!is.na(dalponte2016$DBH_cm), c("sp", "DBH_cm", "C_diam_m", "CR", "C_depth_m", "HT_m", "W", "data", "ref", 
                                                          "latitude_plot", "longitude_plot", "latitude_tree", "longitude_tree", "location_ID", "continent")]
  
  FrenchNFI_data$data <- "FrenchNFI"
  FrenchNFI_data$ref <- "FrenchNFI"
  FrenchNFI_data$C_diam_m <- NA
  FrenchNFI_data$C_depth_m <- NA
  FrenchNFI_data$CR <- NA
  FrenchNFI <- FrenchNFI_data[!is.na(FrenchNFI_data$DBH_cm), c("sp", "DBH_cm", "C_diam_m", "CR", "C_depth_m", "HT_m", "W", "data", "ref", 
                                                               "latitude_plot", "longitude_plot", "latitude_tree", "longitude_tree", "location_ID", "continent")]
  
  Quebec$data <- "Quebec"
  Quebec$ref <- "QuebecNFI"
  Quebec$C_diam_m <- NA
  Quebec$C_depth_m <- NA
  Quebec$CR <- NA
  Quebec_data <- Quebec[!duplicated(Quebec$TreeID) &  !is.na(Quebec$DBH_cm), c("sp", "DBH_cm", "C_diam_m", "CR", "C_depth_m", "HT_m", "W", "data", "ref", 
                                                                               "latitude_plot", "longitude_plot", "latitude_tree", "longitude_tree", "location_ID", "continent")]
  
  Canada_data$data <- "Canada"
  Canada_data$ref <- "CanadaNFI"
  Canada_data$C_diam_m <- NA
  Canada <- Canada_data[!duplicated(Canada_data$TreeID) & !is.na(Canada_data$DBH_cm), c("sp", "DBH_cm", "C_diam_m", "CR", "C_depth_m", "HT_m", "W", "data", "ref", 
                                                                              "latitude_plot", "longitude_plot", "latitude_tree", "longitude_tree", "location_ID", "continent")]
  
  
  Usoltsev$data <- "Usoltsev"
  Usoltsev$ref <- "Usoltsev2015"
  Usoltsev$C_diam_m <- NA
  Usoltsev$C_depth_m <- NA
  Usoltsev$CR <- NA
  Usoltsev_data <- Usoltsev[!is.na(Usoltsev$DBH_cm), c("sp", "DBH_cm", "C_diam_m", "CR", "C_depth_m", "HT_m", "W", "data", "ref", 
                                                       "latitude_plot", "longitude_plot", "latitude_tree", "longitude_tree", "location_ID", "continent")]
  
  Sullivan_data$data <- "Sullivan"
  Sullivan_data$ref <- "Sullivan2017"
  Sullivan <- Sullivan_data[!is.na(Sullivan_data$DBH_cm), c("sp", "DBH_cm", "C_diam_m", "CR", "C_depth_m", "HT_m", "W", "data", "ref", 
                                                            "latitude_plot", "longitude_plot", "latitude_tree", "longitude_tree", "location_ID", "continent")]
  
  
  tallo <- tallo_database[!is.na(tallo_database$DBH_cm), c("sp", "DBH_cm", "C_diam_m", "CR", "C_depth_m", "HT_m", "W", "data", "ref", 
                                                            "latitude_plot", "longitude_plot", "latitude_tree", "longitude_tree", "location_ID", "continent")]
  
  
  df <- bind_rows(baad, FHM, FIA, Spain, FUNDIV, FUNDIV_CROWN,
                  MONTANE, data_paper, evan_crown, 
                  dettmann, heym, gentree, legacytree_data,
                  ICP_data_s, anderson2015, dalponte, 
                  FrenchNFI, Quebec_data, Canada, Usoltsev_data, Sullivan, tallo)

  return(df)
}

# REMOVE DUPLICATE DATA

remove_duplicated_ref <- function(all_crown){
  
  REMOVE_ALL <- c("ES", "FIA_SRS") 
  REMOVE_BAAD <- c("BondLamberty2002",  "Domec2012",  "Fatemi2011",  "Garber2005",  "Maguire1998",   
                   "Martin1998",  "Naidu1998",  "OHara1995",  "OHara2014",  "Reid2004",  "Roberts2003",  "Roberts2004",
                   "Roeh1997",  "Roth2007",  "Stancioiu2005",  "Wang1995",  "Wang2000",  "Whittaker1974")
  df <- all_crown[!all_crown$ref %in% REMOVE_ALL,]
  df <- df[!(df$data == "baad" & df$ref %in% REMOVE_BAAD) ,]
  
  return(df)
  
}





############################### NOT USED IN THIS VERSION OF THE CODE ###############################


# plot to check unit

# plots_check_units <- function(all_crown){

# p1 <-  ggplot(all_crown[!is.na(all_crown$HT_m), ], aes(x = DBH_cm, y = HT_m, color = data))+geom_point(alpha = 0.1)+ theme(legend.position="none")
# p2 <-  ggplot(all_crown[!is.na(all_crown$CR), ], aes(x = DBH_cm, y = CR, color = data))+geom_point(alpha = 0.1)+ theme(legend.position="none")
# p3 <-  ggplot(all_crown[!is.na(all_crown$C_depth_m), ], aes(x = DBH_cm, y = C_depth_m, color = data))+geom_point(alpha = 0.1)+ theme(legend.position="none")
# p4 <-  ggplot(all_crown[!is.na(all_crown$C_diam_m), ], aes(x = DBH_cm, y = C_diam_m, color = data))+geom_point(alpha = 0.1)
# png("figures/unit_check.png", width = 680, height = 680)
# multiplot(p1, p2, p3, p4, cols=  2)
# dev.off()

# }

# plots_check_units_perdata <- function(all_crown){
#  p1 <-  ggplot(all_crown[!is.na(all_crown$HT_m), ], aes(x = DBH_cm, y = HT_m))+geom_point()+facet_wrap(vars(data), nrow = 3)
#  p2 <-  ggplot(all_crown[!is.na(all_crown$CR), ], aes(x = DBH_cm, y = CR))+geom_point()+facet_wrap(vars(data), nrow = 3)
#  p3 <-  ggplot(all_crown[!is.na(all_crown$C_depth_m), ], aes(x = DBH_cm, y = C_depth_m))+geom_point()+facet_wrap(vars(data), nrow = 3)
#  p4 <-  ggplot(all_crown[!is.na(all_crown$C_diam_m), ], aes(x = DBH_cm, y = C_diam_m))+geom_point()+facet_wrap(vars(data), nrow = 3)
#  ggsave("figures/unit_check_perdata1.png", plot = p1, width = 20, height = 20, units = "cm")
#  ggsave("figures/unit_check_perdata2.png", plot = p2, width = 20, height = 20, units = "cm")
#  ggsave("figures/unit_check_perdata3.png", plot = p3, width = 20, height = 20, units = "cm")
#  ggsave("figures/unit_check_perdata4.png", plot = p4, width = 20, height = 20, units = "cm")
# }


#boxplots_check_units <- function(all_crown){

#  p1 <-  ggplot(all_crown, aes(x = data, y = HT_m))+geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#  p2 <-  ggplot(all_crown, aes(x = data, y = C_depth_m))+geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#  p3 <-  ggplot(all_crown, aes(x = data, y = CR))+geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# p4 <-  ggplot(all_crown, aes(x = data, y = C_diam_m))+geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#  png("figures/unit_boxplot_check.png", width = 680, height = 680)
#  multiplot(p1, p2, p3, p4, cols=  2)
#  dev.off()

# }

# compute_nobs_per_sp_per_var <- function(df){
#  df <- df[!is.na(df$sp), ]
#  df$sp <- trimws(df$sp)
#  df$sp[df$sp == "Quercus douglassi"] <- "Quercus douglasii"
  
#  # todo m
#  res <- df %>% group_by(sp) %>% summarise(Nobs_HT = sum(!is.na(HT_m)),
#                                           Nobs_C_diam = sum(!is.na(C_diam_m)),
#                                           Nobs_C_depth = sum(!is.na(C_depth_m)))  
#  write.csv(res, file.path("output", "crown_per_sp.csv"), row.names = FALSE)
#  return(res)
# }

# merge_crown_data_sp <- function(FHM_CD_per_sp, FIA_HT_per_sp, FIA_CR_per_sp,
#                                Spain_Crown_sp, FUNDIV_HT_sp, BAAD_crown_sp,
#                                MONTANE_crown_sp, data_paper_crown_sp, evan_crown_data){
#  
#  head(FHM_CD_per_sp)
#  FHM_CD_per_sp <- FHM_CD_per_sp %>% mutate(sp = paste(GENUS, SPECIES_NAME)) %>% select(sp, CD_obs) 
  
  
#  head(FIA_HT_per_sp)
#  FIA_HT_per_sp <- FIA_HT_per_sp %>% mutate(sp = paste(GENUS, SPECIES)) %>% select(sp, HT_obs) 
#  head(FIA_CR_per_sp)
#  FIA_CR_per_sp <- FIA_CR_per_sp %>% mutate(sp = paste(GENUS, SPECIES)) %>% select(sp, CR_obs) 
#  head(Spain_Crown_sp)
#  Spain_CD_sp <- Spain_Crown_sp %>% mutate(sp = Nombre.Cientifico, 
#                                           CD_obs = N_CD) %>% 
#    select(sp, CD_obs)
#  Spain_CR_sp <- Spain_Crown_sp %>% mutate(sp = Nombre.Cientifico, 
#                                           CR_obs = N_CR) %>% 
#    select(sp, CR_obs)
  
#  head(FUNDIV_HT_sp)
#  FUNDIV_HT_sp <- FUNDIV_HT_sp %>% mutate(HT_obs = N_HT) %>% 
#    select(sp, HT_obs)
#  head(BAAD_crown_sp)
#  BAAD_CR_sp <- BAAD_crown_sp %>% mutate(sp = species, CR_obs = CR_nobs) %>% 
#    select(sp, CR_obs)
#  BAAD_CD_sp <- BAAD_crown_sp %>% mutate(sp = species, CD_obs = CD_nobs) %>% 
#    select(sp, CD_obs)
#  head(MONTANE_crown_sp)
#  MONTANE_CR_sp <- MONTANE_crown_sp %>% select(sp, CR_obs)
#  MONTANE_CD_sp <- MONTANE_crown_sp %>% select(sp, CD_obs)
#  MONTANE_HT_sp <- MONTANE_crown_sp %>% select(sp, HT_obs)
#  head(data_paper_crown_sp)
#  data_paper_CR_sp <- data_paper_crown_sp %>% mutate(sp = Latin.name) %>% select(sp, CR_obs)
#  data_paper_CD_sp <- data_paper_crown_sp %>% mutate(sp = Latin.name) %>% select(sp, CD_obs)
#  data_paper_HT_sp <- data_paper_crown_sp %>% mutate(sp = Latin.name) %>% select(sp, HT_obs)
#  head(evan_crown_data)
#  evan_CR_sp <- evan_crown_data %>% mutate(sp = Species) %>% select(sp, CR_obs)
#  evan_CD_sp <- evan_crown_data %>% mutate(sp = Species) %>% select(sp, CD_obs)
#  evan_HT_sp <- evan_crown_data %>% mutate(sp = Species) %>% select(sp, HT_obs)
  
#  data_HT_sp <- bind_rows(FIA_HT_per_sp, MONTANE_HT_sp, data_paper_HT_sp, evan_HT_sp, FUNDIV_HT_sp)   
#  data_HT_sp <- data_HT_sp %>% group_by(sp) %>% summarise(HT_obs = sum(HT_obs, na.rm = TRUE))
#  data_CR_sp <- bind_rows(FIA_CR_per_sp, BAAD_CR_sp, MONTANE_CR_sp, data_paper_CR_sp, evan_CR_sp, Spain_CR_sp)   
#  data_CR_sp <- data_CR_sp %>% group_by(sp) %>% summarise(CR_obs = sum(CR_obs, na.rm = TRUE))
#  data_CD_sp <- bind_rows(FHM_CD_per_sp, BAAD_CD_sp, MONTANE_CD_sp, data_paper_CD_sp, evan_CD_sp, Spain_CD_sp)   
#  data_CD_sp <- data_CD_sp %>% group_by(sp) %>% summarise(CD_obs = sum(CD_obs, na.rm = TRUE))
#  write.csv(data_HT_sp, "output/data_HT_sp.csv", row.names = FALSE)   
#  write.csv(data_CR_sp, "output/data_CR_sp.csv", row.names = FALSE)   
#  write.csv(data_CD_sp, "output/data_CD_sp.csv", row.names = FALSE)   
# }

#####################################################################################################

