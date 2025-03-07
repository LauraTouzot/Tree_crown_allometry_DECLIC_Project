# Evans et al. 2016 Data
data_evans_crown <- function(){

  df <- read.delim(file.path("data", "doi_10.5061_dryad.2c1s7__v1", "Data record 2.txt"))
  df <- df[, 1:10]
  species <- read.csv(file.path("data", "doi_10.5061_dryad.2c1s7__v1", "Species.csv"))
  df <- left_join(df, species, by = c("Taxonomic_code" = "Species_code"))
  df <- df %>% rename(sp = Species,
                      HT_m = Height_m)
  df <- df %>% mutate(C_diam_m = Crown_radius_m*2,
                      C_depth_m = HT_m - Crown_ht_m,
                      CR = C_depth_m/ HT_m)
  df[df$CR >1 & !is.na(df$CR), c("CR", "C_depth_m")] <- NA
  
  df_W <- df %>% filter(Dataset == "ECN_W") %>% mutate(latitude_plot = 51.776789, longitude_plot = -1.335494)
  df_AH <- df %>% filter(Dataset == "ECN_AH") %>% mutate(latitude_plot = 51.161363, longitude_plot = -0.857419)
  df <- rbind(df_AH, df_W)
  
  df$location_ID <- paste(df$Dataset, "Evans", sep = "_")
  
  # checking minimum DBH
  df <- df[df$DBH_cm > 5,]
  
  df$continent <- "E_U"
  df$latitude_tree <- NA
  df$longitude_tree <- NA
  df$W <- NA
  
  df <- df %>% ungroup()
  
  # plots to check unit within the dataset
  p1 <-  ggplot(df[!is.na(df$HT_m), ], aes(x = DBH_cm, y = HT_m)) + geom_point(alpha = 0.1) 
  p2 <-  ggplot(df[!is.na(df$CR), ], aes(x = DBH_cm, y = CR)) + geom_point(alpha = 0.1) 
  p3 <-  ggplot(df[!is.na(df$C_depth_m), ], aes(x = DBH_cm, y = C_depth_m)) + geom_point(alpha = 0.1) 
  p4 <-  ggplot(df[!is.na(df$C_diam_m), ], aes(x = DBH_cm, y = C_diam_m)) + geom_point(alpha = 0.1) 
  png("figures/unit_check/unit_check_Evans.png", width = 680, height = 680)
  multiplot(p1, p2, p3, p4, cols =  2)
  dev.off()
  
  return(df)
}


############################### NOT USED IN THIS VERSION OF THE CODE ###############################

# data_evan_crown_sp <- function(df){
#   res_crown <-   df %>% group_by(sp) %>% summarise(HT_obs = sum(!is.na(HT_m)),
#                                                    CD_obs = sum(!is.na(C_diam_m)),
#                                                    CR_obs = sum(!is.na(CR)))
#   res_crown
# }

#####################################################################################################


