prepare_FUNDIV_explore_crown <- function() {
  
  FUNDIV_explore_crown = read.csv("data/FUNDIV_Explo_Crown/FunDivEUROPE-treedata.csv")
  FUNDIV_explore_crown$sp = FUNDIV_explore_crown$Species
  FUNDIV_explore_crown$DBH_cm = FUNDIV_explore_crown$Diameter_cm
  FUNDIV_explore_crown$C_diam_m = 2 * sqrt((pi * (FUNDIV_explore_crown$Crown_diameter_x_m/2) * (FUNDIV_explore_crown$Crown_diameter_y_m/2)) / pi)
  FUNDIV_explore_crown$HT_m = FUNDIV_explore_crown$Height_m
  FUNDIV_explore_crown$C_depth_m = FUNDIV_explore_crown$Height_m - FUNDIV_explore_crown$Height_crown_base_m
  FUNDIV_explore_crown$CR = FUNDIV_explore_crown$C_depth_m / FUNDIV_explore_crown$HT_m
  FUNDIV_explore_crown$location_ID = paste(FUNDIV_explore_crown$Plot, "FUNDIV_Explore_Crown", sep = "_")
  FUNDIV_explore_crown = FUNDIV_explore_crown[which(FUNDIV_explore_crown$Height_m > FUNDIV_explore_crown$Height_crown_base_m),]
  # no minimum DBH
  
  return(FUNDIV_explore_crown)
}


get_FUNDIV_coordinates <- function(FUNDIV_crown_data) {
  
  FI <- read_excel("data/FunDiv_ExploratoryPlots/Plot_descriptors_-_Finland.xls", sheet = "Raw data")
  GE <- read_excel("data/FunDiv_ExploratoryPlots/Plot_descriptors_-_Germany.xls", sheet = "Raw data")
  IT <- read_excel("data/FunDiv_ExploratoryPlots/Plot_descriptors_-_Italy.xls", sheet = "Raw data")
  PO <- read_excel("data/FunDiv_ExploratoryPlots/Plot_descriptors_-_Poland.xls", sheet = "Raw data")
  RO <- read_excel("data/FunDiv_ExploratoryPlots/Plot_descriptors_-_Romania.xls", sheet = "Raw data")
  SP <- read_excel("data/FunDiv_ExploratoryPlots/Plot_descriptors_-_Spain.xls", sheet = "Raw data")
  
  FI <- FI %>% dplyr::select(PlotID, Latitude, Longitude)
  GE <- GE %>% dplyr::select(PlotID, Latitude, Longitude)
  IT <- IT %>% dplyr::select(PlotID, Latitude, Longitude)
  PO <- PO %>% dplyr::select(PlotID, Latitude, Longitude)
  RO <- RO %>% dplyr::select(PlotID, Latitude, Longitude)
  SP <- SP %>% dplyr::select(PlotID, Latitude, Longitude)
  
  df <- bind_rows(FI, GE, IT, PO, RO, SP)
  
  df_coords <- df %>% dplyr::select(PlotID, Latitude, Longitude)
  df_complete <- left_join(FUNDIV_crown_data, df_coords, by = c("Plot" = "PlotID"))
  df_complete$latitude_plot <- as.numeric(df_complete$Latitude)
  df_complete$longitude_plot <- as.numeric(df_complete$Longitude)
  
  df_complete$latitude_tree <- NA
  df_complete$longitude_tree <- NA
  df_complete$continent <- "E_U"
  df_complete$W <- NA
  
  df_complete <- df_complete %>% ungroup()
  
  # plots to check unit within the dataset
  p1 <-  ggplot(df_complete[!is.na(df_complete$HT_m), ], aes(x = DBH_cm, y = HT_m)) + geom_point(alpha = 0.1) 
  p2 <-  ggplot(df_complete[!is.na(df_complete$CR), ], aes(x = DBH_cm, y = CR)) + geom_point(alpha = 0.1) 
  p3 <-  ggplot(df_complete[!is.na(df_complete$C_depth_m), ], aes(x = DBH_cm, y = C_depth_m)) + geom_point(alpha = 0.1) 
  p4 <-  ggplot(df_complete[!is.na(df_complete$C_diam_m), ], aes(x = DBH_cm, y = C_diam_m)) + geom_point(alpha = 0.1) 
  png("figures/unit_check/unit_check_FUNDIV_crown.png", width = 680, height = 680)
  multiplot(p1, p2, p3, p4, cols =  2)
  dev.off()
  
  return(df_complete)
  
  }
