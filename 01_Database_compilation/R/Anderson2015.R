# Anderson-Teixeira 2015

Anderson_2015_crown_data <- function(){

  Crowns <- read.csv(file.path("data", "doi_10.5061_dryad.6nc8c__v1",
                                 "SCBI_crown_dimensions.csv"))
  Heights <- read.csv(file.path("data", "doi_10.5061_dryad.6nc8c__v1",
                               "SCBI_tree_heights.csv"))
  Species <- read.csv(file.path("data", "doi_10.5061_dryad.6nc8c__v1",
                                "SCBI_species_acronym.csv"))
  
  df <- left_join(Heights, Crowns[, c("tag", "stemtag", "height.to.base.of.crown.m", "crown.area.m2", "crown.volume.m3")], by = c("tag", "stemtag"))
  df <- left_join(df, Species, by = c("species" = "sp.acronym"))
  
  df$DBH_cm <- conv_unit(df$DBH.mm, "mm", "cm")
  df$HT_m <- df$height.m
  df$C_diam_m <- 2*sqrt(df$crown.area.m2/pi)
  df$C_depth_m <- df$HT_m - df$height.to.base.of.crown.m
  df$CR <- df$C_depth_m/df$HT_m
  
  df <- df[df$DBH_cm >= 1,]
  
  df <- df %>% rename(sp = Species) %>% dplyr::select(sp, DBH_cm, HT_m, C_diam_m, C_depth_m, CR)
  
  df$latitude_plot <- 38.8935
  df$longitude_plot <- -78.1453888888889
  df$latitude_tree <- NA
  df$longitude_tree <- NA
  df$W <- NA
  df$location_ID <- "anderson2015"
  df$continent <- "N_A"
  
  # plots to check unit within the dataset
  p1 <-  ggplot(df[!is.na(df$HT_m), ], aes(x = DBH_cm, y = HT_m)) + geom_point(alpha = 0.1) 
  p2 <-  ggplot(df[!is.na(df$CR), ], aes(x = DBH_cm, y = CR)) + geom_point(alpha = 0.1) 
  p3 <-  ggplot(df[!is.na(df$C_depth_m), ], aes(x = DBH_cm, y = C_depth_m)) + geom_point(alpha = 0.1) 
  p4 <-  ggplot(df[!is.na(df$C_diam_m), ], aes(x = DBH_cm, y = C_diam_m)) + geom_point(alpha = 0.1) 
  png("figures/unit_check/unit_check_Anderson.png", width = 680, height = 680)
  multiplot(p1, p2, p3, p4, cols =  2)
  dev.off()
  
  return(df)
}


