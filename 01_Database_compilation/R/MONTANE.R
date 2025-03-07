### MONTANE DATA
read_MONTANE <- function(){
  
  files <- list.files("data/MONTANE") 
  plot_names <- c("1", "1b", "2", "3", "4", "Premol")
  
  f <- function(x){
    df <-   read.csv(file.path("data/MONTANE",x))
    library(dplyr)
    df <- df %>% mutate(CR = (N+E+S+O)/4) %>% dplyr::select(NUM, ES, CR, H_Tot, H_1bv, DIA, x, y) 
    df$NUM <- as.character(df$NUM)
    return(df)
  } 
  
  list_df <- lapply(files[1:5],f)
  df <-   read.csv(file.path("data/MONTANE",files[6]))
  df <- df %>% mutate(CR = rowMeans(dplyr::select(.,N, E, S,O), na.rm = TRUE)) %>% 
    dplyr::select(NUM, ES, CR, H_Tot, H_1bv, DIA, x, y) 
  df$NUM <- as.character(df$NUM)
  df <- df[!is.na(df$ES),]
  list_df[[6]] <- df
  list_df[[1]]$plot <- plot_names[1]
  list_df[[2]]$plot <- plot_names[2]
  list_df[[3]]$plot <- plot_names[3]
  list_df[[4]]$plot <- plot_names[4]
  list_df[[5]]$plot <- plot_names[5]
  list_df[[6]]$plot <- plot_names[6]
  
  df <- bind_rows(list_df)
  
  sp_code <- c("PIAB", "FASY", "BEPE", "ABAL", "COAV", "PRAV", "POTR", "IN", 
               "SOAR", "FREX", "BEsp", "ACPS", "PIUN", "SOAU", "inconnu", NA)
  sp_names <- c("Picea abies", "Fagus sylvatica", "Betula pendula", "Abies alba", "Corylus avellana", "Prunus avium", "Populus tremula", "inconnu", 
               "Sorbus aria", "Fraxinus excelsior", "Betula", "Acer pseudoplatanus", "Pinus mugo", "Sorbus aucuparia", "inconnu", "NA")
  names(sp_names) <- sp_code
  
  df$sp <- sp_names[df$ES]
  df$C_diam_m <- df$CR * 2 
  df$DBH_cm <- df$DIA
  df$HT_m <- df$H_Tot 
  df$C_depth_m <- df$H_Tot - df$H_1bv
  df$CR <- df$C_depth_m/df$HT_m

  # remove measurement problems
  df$CR[df$C_depth_m < - 2] <- NA
  df$C_depth_m[df$C_depth_m < - 2] <- NA
  df$location_ID <- paste(df$plot, "MONTANE", sep = "_")
  
  # add weight for each location
  plot_areas <- as.data.frame(matrix(nrow = 6, ncol = 2))
  colnames(plot_areas) <- c("plot", "W")
  plot_areas$plot <- c("1","2","3","4","1b","Premol")
  plot_areas$W <- c(1/(50*50*10^-4), 1/(50*50*10^-4), 1/(50*50*10^-4), 1/(50*50*10^-4), 1/(pi * (15^2) * 10^-4), 1/0.8)
  df <- left_join(df, plot_areas, by = "plot")
  
  # add plot location
  dd <- df %>% group_by(plot) %>% summarise(latitude = unique(y, na.rm = TRUE), longitude = unique(x, na.rm = TRUE))    

  extract_coordinates <- dd[,c("latitude", "longitude")]
  
  coordinates(extract_coordinates) <- ~longitude + latitude  
  proj4string(extract_coordinates) <- CRS("+init=epsg:2154")
  extract_coordinates <- spTransform(extract_coordinates, CRS("+init=epsg:4326"))
  coordinatesWGS <- extract_coordinates@coords
  coordinatesWGS <- as.data.frame(coordinatesWGS)
  
  df$latitude_plot <- NULL
  df$longitude_plot <- NULL 
  
  df$latitude_plot <- coordinatesWGS$latitude
  df$longitude_plot <- coordinatesWGS$longitude
  
  df_TreeID <- paste(df$NUM, df$plot, sep = "_")
  coordinates_tree <- df[,c("x", "y")]
  
  coordinates(coordinates_tree) <- ~x + y
  proj4string(coordinates_tree) <- CRS("+init=epsg:2154")
  coordinates_tree <- spTransform(coordinates_tree, CRS("+init=epsg:4326"))
  coordinatesWGS_tree <- coordinates_tree@coords
  coordinatesWGS_tree <- as.data.frame(coordinatesWGS_tree)
  
  df$latitude_tree <- NULL
  df$longitude_tree <- NULL 
  
  df$latitude_tree <- coordinatesWGS_tree$y
  df$longitude_tree <- coordinatesWGS_tree$x
  
  # checking minimum DBH
  df <- df[df$DBH_cm >= 7.5,]
  
  df$continent <- "E_U"
  
  df <- df %>% ungroup()
  
  # plots to check unit within the dataset
  p1 <-  ggplot(df[!is.na(df$HT_m), ], aes(x = DBH_cm, y = HT_m)) + geom_point(alpha = 0.1) 
  p2 <-  ggplot(df[!is.na(df$CR), ], aes(x = DBH_cm, y = CR)) + geom_point(alpha = 0.1) 
  p3 <-  ggplot(df[!is.na(df$C_depth_m), ], aes(x = DBH_cm, y = C_depth_m)) + geom_point(alpha = 0.1) 
  p4 <-  ggplot(df[!is.na(df$C_diam_m), ], aes(x = DBH_cm, y = C_diam_m)) + geom_point(alpha = 0.1) 
  png("figures/unit_check/unit_check_MONTANE.png", width = 680, height = 680)
  multiplot(p1, p2, p3, p4, cols =  2)
  dev.off()
  
  return(df)
}



############################### NOT USED IN THIS VERSION OF THE CODE ###############################

# format_crown_sp_MONTANE <- function(res){
#   res_crown <-   res %>% group_by(sp) %>% summarise(CD_obs = sum(!is.na(CR)),
#                                                  CR_obs = sum(!is.na(H_1bv)),
#                                                  HT_obs = sum(!is.na(H_Tot)))
#   return(res_crown)
# }

#####################################################################################################

