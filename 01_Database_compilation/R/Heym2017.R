# Read Heym 2017 data

read_heym_crown_data <- function() {
  
  Crowns <- read.delim(file.path("data", "doi_10.5061_dryad.8v04m__v1",
                                 "EuMIXFOR_Scots_pine-European_beech_data", "Crowns.txt"), sep = ",")
  df <- Crowns %>% group_by(Triplet, Plot, Nr) %>% 
    mutate(distance_x1 = lead(distance, n = 1), distance_x2 = lag(distance, n = 7)) %>%
    mutate(distance_x3 = coalesce(distance_x1, distance_x2))
  
  Crowns_R <- df %>% group_by(Triplet, Plot, Nr, azimuth) %>%
    mutate(distance_a = max(distance, distance_x3), distance_c = min(distance, distance_x3)) %>% 
    dplyr::select(Triplet, Plot, Nr, azimuth, distance_a, distance_c) 
  
  Semi_axis <- Crowns_R %>% 
               rename(a = distance_a, c = distance_c) %>%
               mutate(b = sin(pi/4)*c/(sqrt(1-cos(pi/4)^2*c^2/a^2)))
  
  Area_ellipse <- Semi_axis %>%
                  mutate(area = a*b/2*(atan(a*tan(pi/4)/b)))
  
  Clean_Areas <- Area_ellipse[!is.na(Area_ellipse$area),]
  Areas <- Clean_Areas %>% group_by(Triplet, Plot, Nr) %>% summarise(crown_surface = sum(area))
  
  Trees <- read.delim(file.path("data", "doi_10.5061_dryad.8v04m__v1",
                                "EuMIXFOR_Scots_pine-European_beech_data", "Trees.txt"), sep = ",")
  
  df <- left_join(Areas, Trees, by = c("Triplet", "Plot", "Nr"))
  
  df <- df %>% rename(HT_m = h, DBH_cm = dbh, sp = species)
  df <- df %>% mutate(C_diam_m = 2 * sqrt(crown_surface/pi),
                      C_depth_m = HT_m - cbh,
                      CR = C_depth_m/HT_m)
  
  # checking for minimum DBH
  df <- df[df$DBH_cm > 7,]

  Plots <- read.delim(file.path("data", "doi_10.5061_dryad.8v04m__v1",
                                "EuMIXFOR_Scots_pine-European_beech_data", "TripletInformation.txt"), sep = ",", header = TRUE, fileEncoding ="latin1")
  
  Plots <- Plots[Plots$species != "Total",]
  
  coords.dd = Plots %>% mutate(longitude_plot = parzer::parse_lon(longitude),
           latitude_plot = parzer::parse_lat(latitude)) %>% dplyr::select(Triplet, Plot, area, longitude_plot, latitude_plot)
  
  plot_data <- coords.dd %>% group_by(Triplet, Plot) %>% mutate(W = 1/area)
  
  data_file <- left_join(df, plot_data, by = c("Triplet", "Plot"))
  data_file$location_ID <- paste(data_file$Triplet, data_file$Plot, "Heym2017", sep = "_")
  
  data_file$latitude_tree <- NA
  data_file$longitude_tree <- NA
  data_file$continent <- "E_U"

  data_file <- data_file %>% ungroup()
  
  # plots to check unit within the dataset
  p1 <-  ggplot(data_file[!is.na(data_file$HT_m), ], aes(x = DBH_cm, y = HT_m)) + geom_point(alpha = 0.1) 
  p2 <-  ggplot(data_file[!is.na(data_file$CR), ], aes(x = DBH_cm, y = CR)) + geom_point(alpha = 0.1) 
  p3 <-  ggplot(data_file[!is.na(data_file$C_depth_m), ], aes(x = DBH_cm, y = C_depth_m)) + geom_point(alpha = 0.1) 
  p4 <-  ggplot(data_file[!is.na(data_file$C_diam_m), ], aes(x = DBH_cm, y = C_diam_m)) + geom_point(alpha = 0.1) 
  png("figures/unit_check/unit_check_Heym.png", width = 680, height = 680)
  multiplot(p1, p2, p3, p4, cols =  2)
  dev.off()

  return(data_file)
  
}
