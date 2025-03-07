# French NFI data from 2008 to 2019

# French_NFI_crown_data <- function(){
#  library(dplyr)
#  df <- read_tree(path = "data/FrenchNFI", years = 2008:2019)
#  species <- read.csv(file.path("data/FrenchNFI","species.csv"))
#  species <- species[!duplicated(species$code), ]
#  species$code <- as.character(species$code)

#  df <- left_join(df, species[ , c("code", "Latin_name")], by = c("espar" = "code"))
#  df <- df %>% mutate(DBH_cm = c13/pi, sp = gsub("_", " ", Latin_name)) %>% rename(HT_m = htot)
#  df[df$DBH_cm < 10 & df$HT_m >40, ]$HT_m <- NA
#  df <- df[!is.na(df$sp),]
#  df <- df %>% select(sp, DBH_cm, HT_m)

#  return(df)
# }


read_ifn_data <- function(variable, path, years,
                          addyear = TRUE, zipped = TRUE, ...) {
  if (zipped) {
    zipfiles <- list.files(path, pattern = "\\.zip$")
    zipfiles <- file.path(path, zipfiles[sapply(years, grep, x = zipfiles)])
    files <- file.path(paste0(years, "-fr"), paste0(variable, "_", years, ".csv"))
    connections <- mapply(unz, zipfiles, files, SIMPLIFY = FALSE)
  } else {
    connections <- file.path(path, years,
                             paste0(variable, "_", years, ".csv"))
  }
  data_list <- lapply(connections, utils::read.table,
                      header = TRUE, sep = ";", dec = ".",
                      encoding = "UTF-8",
                      stringsAsFactors = FALSE)
  
  if (addyear) {
    for (i in seq_along(data_list)) {
      names(data_list[[i]]) <- gsub("\\.", "", gsub("X", "", names(data_list[[i]])))
      data_list[[i]]$year = years[[i]]
    }
  }
  data.table::rbindlist(data_list, use.names = TRUE, fill = TRUE)
}


french_NFI_crown_data <- function(){

  df <- read_ifn_data(variable = "arbres_foret", path = "data/FrenchNFI", years = 2008:2019)
  species <- read.csv(file.path("data/FrenchNFI","species.csv"))
  species <- species[!duplicated(species$code), ]
  species$code <- as.character(species$code)
  
  df <- left_join(df, species[ , c("code", "Latin_name")], by = c("espar" = "code"))
  df <- df %>% mutate(DBH_cm = c13/pi, sp = gsub("_", " ", Latin_name)) %>% rename(HT_m = htot)
  df[df$DBH_cm < 10 & df$HT_m >40, ]$HT_m <- NA
  df <- df[!is.na(df$sp),]
  df <- df %>% dplyr::select(sp, DBH_cm, HT_m, idp, w)
  
  # checking for minimum DBH
  df <- df[df$DBH_cm >= 7.5,]
  
  coordinates <- read_ifn_data(variable = "placettes_foret", path = "data/FrenchNFI", years = 2008:2019)
  coordinates <- coordinates[,c("idp", "xl93", "yl93")]
  colnames(coordinates) <- c("idp", "longitude", "latitude")
  
  df <- left_join(df, coordinates, by = "idp")
  
  extract_coordinates <- df[,c("latitude", "longitude")]
  
  coordinates(extract_coordinates) <- ~longitude + latitude  
  proj4string(extract_coordinates) <- CRS("+init=epsg:2154")
  extract_coordinates <- spTransform(extract_coordinates, CRS("+init=epsg:4326"))
  coordinatesWGS <- extract_coordinates@coords
  coordinatesWGS <- as.data.frame(coordinatesWGS)
  
  df$latitude_plot <- NULL
  df$longitude_plot <- NULL 
  
  df$latitude_plot <- coordinatesWGS$latitude
  df$longitude_plot <- coordinatesWGS$longitude
  
  df$location_ID <- paste(df$idp, "FrenchNFI", sep = "_")
  df$continent <- "E_U"
  df$latitude_tree <- NA
  df$longitude_tree <- NA
  df$W <- df$w
  
  # plots to check unit within the dataset
  p1 <-  ggplot(df[!is.na(df$HT_m), ], aes(x = DBH_cm, y = HT_m)) + geom_point(alpha = 0.1) 
  png("figures/unit_check/unit_check_FrenchNFI.png", width = 340, height = 340)
  p1
  dev.off()
  
  return(df)
  
}




