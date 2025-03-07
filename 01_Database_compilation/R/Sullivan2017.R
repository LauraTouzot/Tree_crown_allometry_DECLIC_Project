# READ Sullivan et al. 2017 data

# plot: randomly assigned name
# species: genus/species 6 letter codes (4 genus, 2 species)
# dbh: diameter at breast height (unit: centimeter / missing value: NA)
# height: top tree height in meters, measured w/ hypsometer (unit: meter / missing value: NA)
# cbh: height in meters of the base of the crown (unit: meter / missing value: NA)

# species.code: species code
# species.name: genus and species name

read_Sullivan <- function(){
  
  df_Sullivan <- read.csv("data/Sullivan2017/hf339-01-hf-tree-2013.csv", sep = ",")
  df_species <-  read.csv("data/Sullivan2017/hf339-02-species.csv", sep = ",")
  
  colnames(df_Sullivan)[3] = "species.code"

  df <- left_join(df_Sullivan, df_species, by = "species.code")
  
  # crown data
  df <- df %>% rename(sp = species.name,
                      HT_m = height,
                      DBH_cm = dbh)
  df$C_depth_m <- df$HT_m - df$cbh
  df$CR <- df$C_depth_m / df$HT_m
  df$C_diam_m <- 2 * sqrt((pi * df$rad.in * df$rad.out) / pi)
  df <- df[df$C_diam_m > 0,]
  df <- df[df$DBH_cm >= 5,] # checking for minimum DBH
    
  # locations of trees and quadrats (from UTM zone 19 coordinates to WGS84)
  df <- df %>% group_by(plot) %>% mutate(latitude_plot = mean(sensor.est.n), longitude_plot = mean(sensor.est.e))
  
  extract_coordinates_plots <- df[,c("plot", "longitude_plot", "latitude_plot")]
  extract_coordinates_plots <- extract_coordinates_plots[!is.na(extract_coordinates_plots$longitude_plot) & !is.na(extract_coordinates_plots$latitude_plot),]
  
  unique_coordinates_plots <- extract_coordinates_plots[!duplicated(extract_coordinates_plots$plot),]
  plot <- as.data.frame(unique_coordinates_plots$plot)
  colnames(plot) <- "plot"
  
  coordinates(unique_coordinates_plots) <- ~longitude_plot + latitude_plot  
  proj4string(unique_coordinates_plots) <- CRS("+init=epsg:3726")
  unique_coordinates_plots <- spTransform(unique_coordinates_plots, CRS("+init=epsg:4326"))
  coordinatesWGS_plots <- unique_coordinates_plots@coords
  coordinatesWGS_plots <- as.data.frame(coordinatesWGS_plots)
  coordinatesWGS_plots <- cbind(plot, coordinatesWGS_plots)
  
  df <- df %>% select(-c(latitude_plot, longitude_plot))
  df <- left_join(df, coordinatesWGS_plots, by = "plot")
  
  df <- df[!is.na(df$plot),]
  
  extract_coordinates_trees <- df[,c("tag", "tree.easting", "tree.northing")]
  tag <- as.data.frame(extract_coordinates_trees$tag)
  colnames(tag) <- "tag"
  
  coordinates(extract_coordinates_trees) <- ~tree.easting + tree.northing  
  proj4string(extract_coordinates_trees) <- CRS("+init=epsg:3726")
  extract_coordinates_trees <- spTransform(extract_coordinates_trees, CRS("+init=epsg:4326"))
  coordinatesWGS_trees <- extract_coordinates_trees@coords
  coordinatesWGS_trees <- as.data.frame(coordinatesWGS_trees)
  coordinatesWGS_trees <- cbind(tag, coordinatesWGS_trees)
  colnames(coordinatesWGS_trees) <- c("tag", "longitude_tree", "latitude_tree")
  
  df <- left_join(df, coordinatesWGS_trees, by = "tag")
  
  df$location_ID <- paste(df$plot, "Sullivan", sep = "_")
  df$continent <- "N_A"
  df$W <- NA
  
  df <- df %>% ungroup()
  
  # plots to check unit within the dataset
  p1 <-  ggplot(df[!is.na(df$HT_m), ], aes(x = DBH_cm, y = HT_m)) + geom_point(alpha = 0.1) 
  p2 <-  ggplot(df[!is.na(df$CR), ], aes(x = DBH_cm, y = CR)) + geom_point(alpha = 0.1) 
  p3 <-  ggplot(df[!is.na(df$C_depth_m), ], aes(x = DBH_cm, y = C_depth_m)) + geom_point(alpha = 0.1) 
  p4 <-  ggplot(df[!is.na(df$C_diam_m), ], aes(x = DBH_cm, y = C_diam_m)) + geom_point(alpha = 0.1) 
  png("figures/unit_check/unit_check_Sullivan.png", width = 680, height = 680)
  multiplot(p1, p2, p3, p4, cols =  2)
  dev.off()
  
  return(df)
  
}
