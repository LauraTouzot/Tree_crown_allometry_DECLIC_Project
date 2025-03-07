#Read Legacy Tree data

Legacy_Tree_crown_data <- function(){

  df <- read.csv(file.path("data", "LegacyTree","tree.csv"), sep = ",")

  # SPCD FIA SPECIES CODE
  # DBH_cm ST_OB_D_B inches
  # HT_m TR_HT in feets
  # crown_height ST_HT_BLC in feets
  # Foliage dry weight FOL_DW in pounds
  # C_diam_m CR_WID in feets
  # FOL_DW  in pounds
  # FOL_AR in square foot
  
  df <- df[df$ST_OB_D_BH >= 4,] # checking for minimum DBH
  df$DBH_cm <- conv_unit(df$ST_OB_D_BH, "inch", "cm")
  df$HT_m <- conv_unit(df$TR_HT, "ft", "m")
  df$C_diam_m <- conv_unit(df$CR_WID, "ft", "m")
  df$C_depth_m <- NA # df$HT_m - conv_unit(df$ST_HT_BLC, "ft", "m")
  df$CR <- NA #df$CROWN_RATIO/100
  df$Leaf_Biomass_Dry <-  conv_unit(df$FOL_DW, "lbs", "kg")
  df$Leaf_Area <-  conv_unit(df$FOL_AR, "ft2", "m2")
  species_code <- read.csv("data/FIA/REF_SPECIES.csv") %>% mutate(sp = paste(GENUS, SPECIES, sep = " ")) %>%
    dplyr::select("SPCD", "sp")
  # df[df$AUTHOR == "NC75_Chapman", ]$Leaf_Biomass_Dry <- NA # unit pb with this data obs
  df[df$HT_m > 150 & !is.na(df$HT_m), ]$HT_m <- NA # unit pb with this data obs
  df <- left_join(df, species_code, by = "SPCD") %>% dplyr::select(sp, DBH_cm, HT_m, C_diam_m, C_depth_m, CR, Leaf_Biomass_Dry, Leaf_Area, AUTHOR, LOC, SPCD)

  coordinates <- read.csv(file.path("data", "LegacyTree","location.csv"), sep = ",")
  coordinates$location_ID <- paste(coordinates$AUTHOR, coordinates$LOC, coordinates$SPCD, sep = "_")
  coordinates <- coordinates[,c("LAT", "LON", "location_ID")] 
  coordinates <- coordinates %>% group_by(location_ID) %>% summarize(latitude_plot = unique(LAT), longitude_plot = unique(LON))
  
  df$location_ID <- paste(df$AUTHOR, df$LOC, df$SPCD, sep = "_")
  df <- left_join(coordinates, df, by = "location_ID")
  df$location_ID <- paste(df$location_ID, "LegacyTree", sep = "_")
  
  df$continent <- "N_A"
  df$latitude_tree <- NA
  df$longitude_tree <- NA
  df$W <- NA
  
  # plots to check unit within the dataset
  p1 <-  ggplot(df[!is.na(df$HT_m), ], aes(x = DBH_cm, y = HT_m)) + geom_point(alpha = 0.1) 
  p2 <-  ggplot(df[!is.na(df$C_diam_m), ], aes(x = DBH_cm, y = C_diam_m)) + geom_point(alpha = 0.1) 
  png("figures/unit_check/unit_check_LegacyTree.png", width = 680, height = 340)
  multiplot(p1, p2, cols =  2)
  dev.off()
  
  return(df)
}

