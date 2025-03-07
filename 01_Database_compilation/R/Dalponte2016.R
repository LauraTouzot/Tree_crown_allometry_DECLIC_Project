# Read Dalponte and Coomes 

Dalponte_2016_crown_data <- function(){

  Crowns <- read.csv(file.path("data", "doi_10.5061_dryad.hf5rh__v1",
                               "47_validation_plots_field_measured_trees.csv"))

  df <- Crowns %>% mutate(sp = gsub("_", " ", Species),
                          C_diam_m = 2*sqrt(Crown.area.m2/pi)) %>% 
    rename(DBH_cm = DBH.cm, HT_m = Height.m) %>% 
    dplyr::select(sp, DBH_cm, HT_m, C_diam_m, PlotID)
  
  df <- df[df$DBH_cm > 4,] # cleaning data based on min dbh of surveyed trees
  
  df$latitude_plot <- NA
  df$longitude_plot <- NA
  df$latitude_tree <- NA
  df$longitude_tree <- NA
  df$location_ID <- paste(df$PlotID, "Dalponte", sep = "_")
  df$continent <- "E_U"
  df$W <- 1/(pi * 15^2 * 10^-4)
  
  # plots to check unit within the dataset
  p1 <-  ggplot(df[!is.na(df$HT_m), ], aes(x = DBH_cm, y = HT_m)) + geom_point(alpha = 0.1) 
  p2 <-  ggplot(df[!is.na(df$C_diam_m), ], aes(x = DBH_cm, y = C_diam_m)) + geom_point(alpha = 0.1) 
  png("figures/unit_check/unit_check_Dalponte.png", width = 680, height = 340)
  multiplot(p1, p2, cols =  2)
  dev.off()
  
  return(df)
}



