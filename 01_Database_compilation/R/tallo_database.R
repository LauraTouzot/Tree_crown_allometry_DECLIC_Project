read_Tallo <- function() {
  
  df <- read.csv(file = "data/Tallo/Tallo.csv", sep = ",")
  
  # keep selected references only 
  ref_to_keep <- c(1, 11, 21, 47, 54, 57)
  data <- df[df$reference_id %in% ref_to_keep,]

  # selection variables and changing names
  data_b <- data %>% mutate(C_diam_m = crown_radius_m * 2) %>%
    select(species, stem_diameter_cm, height_m, C_diam_m, latitude, longitude, reference_id) %>%
    rename(sp = species, DBH_cm = stem_diameter_cm, HT_m = height_m, latitude_plot = latitude,
           longitude_plot = longitude, ref = reference_id)
  
  data_b$latitude_tree <- NA
  data_b$longitude_tree <- NA
  data_b$W <- NA
  data_b$CR <- NA
  data_b$C_depth_m <- NA
  data_b$ref <- paste(data_b$ref, "Tallo", sep = "_")
  data_b$data <- "Tallo"
  
  # adding continents based on plot coordinates
  states_map <- readOGR(dsn="data/ne_10m_admin_1_states_provinces/", 
                        layer="ne_10m_admin_1_states_provinces")
  
  countriesSP <- getMap(resolution = 'low')

  ID_unique <- as.data.frame(c(1:dim(data_b)[1]))
  colnames(ID_unique) <- "ID_unique"
  data_b <- cbind(ID_unique, data_b)

  coordinates <- data_b[,c("longitude_plot", "latitude_plot")]
  pointsSP = SpatialPoints(coordinates, proj4string = CRS(" +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
  indices = over(pointsSP, countriesSP)
  
  locations <- as.data.frame(indices$REGION)
  colnames(locations) <- "cont"
  locations <- cbind(data_b$ID_unique, locations)
  colnames(locations) <- c("ID_unique", "cont")
  
  continent <- as.data.frame(c("E_U", "N_A", "A_S", "S_A", "A_F", "A_U"))
  cont <- as.data.frame(c("Europe", "North America", "Asia", "South America", "Africa", "Australia"))
  continents <- cbind(continent, cont)
  colnames(continents) <- c("continent", "cont")
  
  locations <- left_join(locations, continents, by = "cont")
  
  tallo_data <- left_join(locations, data_b, by = "ID_unique")
  
  # adding location_ID
  plots <- as.data.frame(paste(tallo_data$latitude_plot, tallo_data$longitude_plot, sep = "_"))
  colnames(plots) <- "plots"
  tallo_data <- cbind(tallo_data, plots)
  
  locs <- as.data.frame(unique(tallo_data$plots))
  locs$location_ID <- paste(c(1:dim(locs)[1]), "Tallo", sep = "_")
  colnames(locs) <- c("plots", "location_ID")

  tallo_data <- left_join(tallo_data, locs, by = "plots")  

  # selecting required variables
  tallo_data <- tallo_data %>% select(sp, DBH_cm, C_diam_m, CR, C_depth_m, HT_m, W, data,
                                      ref, latitude_plot, longitude_plot, latitude_tree, longitude_tree,
                                      location_ID, continent)
  
  
  # adding weight for reference #1 (Aakala, Finish plots)
  tallo_data[tallo_data$ref == "1_Tallo",]$W <- 10

  return(tallo_data)
  
}
