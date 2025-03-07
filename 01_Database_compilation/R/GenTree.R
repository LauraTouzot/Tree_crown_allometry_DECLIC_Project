# Read GenTree data

GenTree_data <- function(){
  
  df <- read.csv(file.path("data", "GenTree", "GenTreeData.csv"))

  species_code <-  c("AA", "BP", "FS", "PA", "PC", "PH", "PN", "PO", "PP", "PS", "QP", "TB")
  species <- c("Abies Alba" , "Betula pendula", "Fagus sylvatica", "Picea abies", "Pinus cembra",
               "Pinus halepensis", "Pinus nigra", "Populus nigra", "Pinus pinaster", "Pinus sylvestris",
               "Quercus petraea", "Taxus baccata")
  names(species) <- species_code
  df$sp <- species[df$m01.spec]
  df$C_diam_m <- 2 * sqrt((pi * df$p07.canopy.1/2 * df$p08.canopy.2/2) / pi)
  df <-   df %>% rename(HT_m = p01.height,
           DBH_cm = p02.dbh,
           ) %>%
    dplyr::select(m06.tree.id, m05.tree.num,  m04.pop, m08.lat, m09.lon, m01.spec, sp, m04.pop,HT_m, DBH_cm, C_diam_m)
  
  # no minimum DBH
  
  df <- df %>% group_by(m04.pop) %>% mutate(latitude_tree = m08.lat, longitude_tree = m09.lon, 
                                            latitude_plot = mean(m08.lat, na.rm = TRUE), longitude_plot = mean(m09.lon, na.rm = TRUE))
  
  df$location_ID <- paste(df$m04.pop, "GenTree", sep = "_")
  df$continent <- "E_U"
  df$W <- NA
  
  df <- df %>% ungroup()
  
  # plots to check unit within the dataset
  p1 <-  ggplot(df[!is.na(df$HT_m), ], aes(x = DBH_cm, y = HT_m)) + geom_point(alpha = 0.1) 
  p2 <-  ggplot(df[!is.na(df$C_diam_m), ], aes(x = DBH_cm, y = C_diam_m)) + geom_point(alpha = 0.1) 
  png("figures/unit_check/unit_check_GenTree.png", width = 680, height = 340)
  multiplot(p1, p2, cols =  2)
  dev.off()

  return(df)
}


