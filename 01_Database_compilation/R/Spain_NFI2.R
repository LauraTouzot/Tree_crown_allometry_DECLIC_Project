
read_spanish_data <- function(){
  
  df_t <- read.table("data/SpainNFI/tree2_clean.csv",
                     header = TRUE, sep = ",", encoding="UTF-8")
  df_c <- read.table("data/SpainNFI/tree2_tipo.csv",
                     header = TRUE, sep = ",")
  df_c <- df_c[!duplicated(df_c$ID_Pma2), ] # there is some duplicated line in df_c (with same value for more 800 individuals) not sure why but there are identical; remove them. Need to check with Paloma
  df <- df_t %>%  left_join(df_c, by = c("ID_Pma2"))
  species <- read.table("data/SpainNFI/Species.txt",
                        header = TRUE,  sep = "\t")
  species$CodeIFN2 <- as.character(species$CodeIFN2)
  species <- species[as.character(species$CodeIFN2) != "xxx" &
                       !is.na(as.character(species$CodeIFN2)) &
                       as.character(species$CodeIFN2) != "", ]
  species$CodeIFN2 <- as.integer(species$CodeIFN2)
  species <- species[!duplicated(species$CodeIFN2), ]
  df <- merge(df, species, by.x = "Especie2", by.y = "CodeIFN2" )
  df$DIAMCOPA <-   2 * (sqrt((pi * df$DIAMCOPA1/2 * df$DIAMCOPA2/2) / pi)) # mean crown diameter
  df[df$h2 >60 & !is.na(df$h2), ]$h2 <- df[df$h2 > 60 & !is.na(df$h2), ]$h2/10 # Correct an error a 15cm dbh Pinus pinaster witha  height of 80m assign 8 = 80/10
  df[df$DIAMCOPA >27 & !is.na(df$DIAMCOPA), ]$DIAMCOPA <- NA  # Remove error in data DIAMCOP too large
  df[df$Forma2 == 6,]$DIAMCOPA <- NA # Forma2 == 6 means that trees have been recently managed (i.e. pruned)
  df[df$Forma2 == 6,]$h2 <- NA
  
  
  ############################### NOT USED IN THIS VERSION OF THE CODE ###############################
  
  ## NEED TO USE ALTUPV but available only for small trees (main trunk branching before 4m of height)  
#  png("figures/Spanish_NFI_Check_Provincia.png",    
#      width     = 3.25,
#      height    = 3.25,
#      units     = "in",
#      res       = 1200,
#      pointsize = 4)    
#  par(mfrow = c(2,2))
#  boxplot(df$DIAMCOPA~df$Provincia2, ylab = "crown diameter for tipo trees", xlab = "Provincia")
#  boxplot(df$ALTUPV~df$Provincia2, ylab = "Height first green branch for typo trees type 4", xlab = "Provincia")
#  boxplot(df$h2~df$Provincia2, ylab = "Tree height for all trees", xlab = "Provincia")
#  dev.off()
  
#  png("figures/Spanish_NFI_Allo.png",
#      width     = 3.25,
#      height    = 3.25,
#      units     = "in",
#      res       = 1200,
#      pointsize = 4)    
#  par(mfrow = c(2,2))
#  plot(df$dbh2, df$h2, cex = 0.1, xlab = "dbh", ylab = "height")
#  plot(df$dbh2, df$DIAMCOPA, cex = 0.1, xlab = "dbh", ylab = "Crown diameter for tipo trees")
#  plot(df$dbh2, df$ALTUPV, cex = 0.1,  xlab = "dbh", ylab = "Height first green branch for typo trees type 4")
#  plot(df$h2, df$ALTUPV, cex = 0.1, xlab = "height", ylab = "Height first green branch for typo trees type 4")
#  abline(a = 0, b = 1, col = "red")
#  dev.off()
  
  #####################################################################################################
  
  df$C_diam_m <- df$DIAMCOPA 
  df$DBH_cm <- df$dbh2/10 # There is a problem this is probably not the good dbh as the match is much better with DIAMETR4 
  df$HT_m <- df$h2 
  df$C_depth_m <- df$h2 - df$ALTUPV
  df$CR <- df$C_depth_m/df$h2
  df$C_depth_m[df$C_depth_m <0 & !is.na(df$CR) ] <- NA
  df$CR <- df$C_depth_m/df$h2
  
  df <- df %>% filter(!is.na(DBH_cm))
  
  # checking minimum DBH
  df <- df[df$DBH_cm >= 7.5,]

  species <-read.table("data/SpainNFI/Species.txt",sep = "\t", header = TRUE, encoding="UTF-8")
  species  <- species[species$CodeIFN2 != "xxx",]
  species  <- species[species$CodeIFN2 != "",]
  species$CodeIFN2  <- as.numeric(species$CodeIFN2)
  species$Nombre.Cientifico <- as.character(species$Nombre.Cientifico)
  res2 <- left_join(df, species, by = c("Especie2" = "CodeIFN2"))
  df$sp <- df$Nombre.Cientifico
  
  plots <- read.table("data/SpainNFI/all_coords2.csv", sep = ",", header = TRUE, encoding = "UTF-8")
  plots <- na.omit(plots)
  location_ID <- as.data.frame(plots$Plotcode)
  extract_coordinates <- plots[,c("CX", "CY")]
  colnames(extract_coordinates) <- c("longitude", "latitude")
  coordinates(extract_coordinates) <- ~longitude + latitude  
  proj4string(extract_coordinates) <- CRS("+init=epsg:23030")
  extract_coordinates <- spTransform(extract_coordinates, CRS("+init=epsg:4326"))
  coordinatesWGS <- extract_coordinates@coords
  coordinatesWGS <- as.data.frame(coordinatesWGS)
  coordinates_Spain <- cbind(location_ID, coordinatesWGS)
  colnames(coordinates_Spain) <- c("Plotcode2", "longitude", "latitude")
  
  df <- left_join(df, coordinates_Spain, by = "Plotcode2")
  
  Spain_data <- df %>% mutate(W = ifelse(DBH_cm >= 42.5, 1/(pi * 25^2 *10^-4), 
                                  ifelse(DBH_cm >= 22.5, 1/(pi * 15^2 *10^-4), 
                                  ifelse(DBH_cm >= 12.5, 1/(pi * 10^2 *10^-4), 1/(pi * 5^2 *10^-4)))))
  
  Spain_data$location_ID <- paste(Spain_data$Plotcode2, "SpainNFI", sep = "_")
  
  Spain_data <- Spain_data %>% ungroup()
  
  Spain_data$continent <- "E_U"
  Spain_data$latitude_plot <- Spain_data$latitude
  Spain_data$longitude_plot <- Spain_data$longitude
  Spain_data$latitude_tree <- NA
  Spain_data$longitude_tree <- NA
  
  # plots to check unit within the dataset
  p1 <-  ggplot(Spain_data[!is.na(Spain_data$HT_m), ], aes(x = DBH_cm, y = HT_m)) + geom_point(alpha = 0.1) 
  p2 <-  ggplot(Spain_data[!is.na(Spain_data$CR), ], aes(x = DBH_cm, y = CR)) + geom_point(alpha = 0.1) 
  p3 <-  ggplot(Spain_data[!is.na(Spain_data$C_depth_m), ], aes(x = DBH_cm, y = C_depth_m)) + geom_point(alpha = 0.1) 
  p4 <-  ggplot(Spain_data[!is.na(Spain_data$C_diam_m), ], aes(x = DBH_cm, y = C_diam_m)) + geom_point(alpha = 0.1) 
  png("figures/unit_check/unit_check_Spain.png", width = 680, height = 680)
  multiplot(p1, p2, p3, p4, cols =  2)
  dev.off()

  return(Spain_data)
}




############################### NOT USED IN THIS VERSION OF THE CODE ###############################

# get_crown_per_sp_spain <- function(df){
#   species <-read.table("data/SpainNFI/Species.txt",sep = "\t", header = TRUE, encoding="UTF-8")
#   species  <- species[species$CodeIFN2 != "xxx",]
#   species  <- species[species$CodeIFN2 != "",]
#   
#   species$CodeIFN2  <- as.numeric(species$CodeIFN2)
#   species$Nombre.Cientifico <- as.character(species$Nombre.Cientifico)
#   library(dplyr)
#   res <- df %>% group_by(Especie2) %>% summarise(N_CD = sum(!is.na(DIAMCOPA)),
#                                           N_CR = sum(!is.na(ALTUPV)))
#   
#   res2 <- left_join( res, species, by = c("Especie2" = "CodeIFN2"))
#   return(res2)
# }

#####################################################################################################



