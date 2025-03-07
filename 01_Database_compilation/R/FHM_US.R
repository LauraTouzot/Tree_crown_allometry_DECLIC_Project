#' Read IFN data for each year
#'
#' @description This is the base function to process the different types of FHM
#' data. The FHM data must already be downloaded and can be zipped or unzipped.
#' @export
#' @param variable Which variable to read
#' @param path Path to data folder
#' @param states Vector with states for which data are to be read
#' @param ... Further arguments to @seealso [data.table::fread()].
#' @details The variable must be either: "TREE", TODO
#' @return a data.table object

read_FHM_data <- function(variable = "TREE", path = "data/FHM", 
                          states = c("AL", "CA", "CO", "CT", "DE", "FL", 
                                     "GA", "ID", "IL", "IN", "MA", "MD", "ME", 
                                     "MI", "MN", "MO", "NC", "NH", "NJ", "NV", 
                                     "OR", "PA", "RI", "SC", "TN", "UT", "VA", 
                                     "VT", "WA", "WI", "WV", "WY"),
                          ...) {
    zipfiles <- list.files(path, pattern = "\\.zip$")
    zipfiles <- file.path(path, zipfiles[sapply(states, grep, x = zipfiles)])
    lapply(zipfiles, unzip, exdir = path)
    
    f <- function(x, variable){
        df <-  unzip(x,list = TRUE)
        x1 <- sub(".zip", "", sub("data/FHM/", "", x))
        return(grep(paste0(x1,"_", variable, "_"), df$Name,value = TRUE))
    }
    
    ls_files <- lapply(zipfiles,f, variable = variable)
    
  read.FHM_State_year <- function(file){
    res <- utils::read.csv(file,header = TRUE, sep = ",", 
                           dec = ".",quote = "",
                           stringsAsFactors = FALSE)
    names(res) <- gsub("_19..", "", names(res))
    return(res)
  }
  
  # according to Purves et al. 2007 before 1999 no height data but crown radius data 
  # after 1999 height measurement but no crown radius data. 
  # before 1999 some Tree height to be available in SITE_TREE variable TREE_HEIGHT but only 8000 obs
  
  data_list <- lapply(file.path(path, unlist(ls_files)), read.FHM_State_year)
  df <- data.table::rbindlist(data_list, use.names = TRUE, fill = TRUE)
  
  species <- read.csv("data/FHM/FHM_Species_codes.csv", sep = ";")
  df <- left_join(df, dplyr::select(species, FHM_SPECIES, COMMON_NAME, GENUS, SPECIES_NAME), by = "FHM_SPECIES")
  df <- df %>% filter(CURRENT_TREE_HISTORY %in% 1:3)
  
  df$DBH_cm <- conv_unit(df$DBH, "inch", "cm")
  df$HT_m <- conv_unit(df$TREE_HEIGHT, "ft", "m")
  df$HT_m[df$MENSURATION_YEAR == 1992] <- NA  # because weird measurement in 1992
  df$CR <- df$CROWN_RATIO/100
  df$C_diam_m <- 2* (sqrt((pi * (conv_unit(df$CROWN_DIAMETER_WIDE, "ft", "m")/2) * (conv_unit(df$CROWN_DIAMETER_90, "ft", "m")/2))/pi))
  df$C_depth_m <- df$HT_m * df$CR
  
  df <- df %>% mutate(TreeID = paste(STATE_ABBREV, COUNTY, FHM_REGION, PROJECT,  PLOT_NBR, POINT_NBR, TREE_NBR))
  df <- df %>% mutate(ID = paste(STATE_ABBREV, COUNTY, FHM_REGION, PROJECT, PLOT_NBR, "FHM", sep = "_"))
  df <- df %>% dplyr::mutate(sp = paste(GENUS, SPECIES_NAME))  
  
  # for each plot, select year with the highest amount of data
  summary_perlocation <- df %>% dplyr::select(ID, YEAR, HT_m, CR, C_diam_m) %>%
    gather(key = variable, value = value, HT_m, CR, C_diam_m) %>%
    filter(!is.na(value)) %>%
    group_by(ID, YEAR, variable) %>%
    summarize(n = n()) %>%
    spread(key = variable, value = n)
  
  summary_perlocation$nobs <- as.numeric(apply(summary_perlocation[,c(3:5)], 1, sum, na.rm = TRUE))
  
  to_keep <- summary_perlocation %>% group_by(ID) %>% 
    arrange(desc(nobs), .by_group = TRUE) %>%
    slice(stay = unique(max(nobs)),.by_group = TRUE) 
  
  to_keep <- to_keep[!duplicated(to_keep$ID),]
  to_keep$staying_IDs <- paste(to_keep$ID, to_keep$YEAR, sep = "_")
  
  df <- df %>% mutate(location_ID = paste(STATE_ABBREV, COUNTY, FHM_REGION, PROJECT, PLOT_NBR, "FHM", YEAR, sep = "_"))
  
  df <- df[df$location_ID %in% to_keep$staying_IDs,]

  
  df$latitude_plot <- NA 
  df$longitude_plot <- NA 
  df$latitude_tree <- NA 
  df$longitude_tree <- NA
  df$continent <- "N_A"
  
  return(df)
}


get_FHM_weight <- function(FHM) {

  plot_areas <- as.data.frame(matrix(nrow = 2, ncol = 2))
  colnames(plot_areas) <- c("PLOT_TYPE", "W")
  plot_areas$PLOT_TYPE <- c(1,2)
  plot_areas$W <- c(10000/(4*(conv_unit(24, "ft", "m")^2*pi)), 10000/(4*(conv_unit(6.8, "ft", "m")^2*pi)))
  FHM_data <- left_join(FHM, plot_areas, by = "PLOT_TYPE")
  FHM_data <- FHM_data %>% ungroup()
  
  # plots to check unit within the dataset
  p1 <-  ggplot(FHM_data[!is.na(FHM_data$HT_m), ], aes(x = DBH_cm, y = HT_m)) + geom_point(alpha = 0.1) 
  p2 <-  ggplot(FHM_data[!is.na(FHM_data$CR), ], aes(x = DBH_cm, y = CR)) + geom_point(alpha = 0.1) 
  p3 <-  ggplot(FHM_data[!is.na(FHM_data$C_depth_m), ], aes(x = DBH_cm, y = C_depth_m)) + geom_point(alpha = 0.1) 
  p4 <-  ggplot(FHM_data[!is.na(FHM_data$C_diam_m), ], aes(x = DBH_cm, y = C_diam_m)) + geom_point(alpha = 0.1) 
  png("figures/unit_check/unit_check_FHM.png", width = 680, height = 680)
  multiplot(p1, p2, p3, p4, cols =  2)
  dev.off()
  
  return(FHM_data)
  
  # no minimum DBH in FHM data
}



############################### NOT USED IN THIS VERSION OF THE CODE ###############################

# format_FHM_sp <- function(df){
#  res_CD <-   df %>% group_by(FHM_SPECIES) %>% summarise(CD_obs = sum(!is.na(CROWN_DIAMETER_WIDE)))
#  sp <- read.csv("data/FHM/FHM_Species_codes.csv", sep = ";")
#  res_CD <- left_join(res_CD, sp, by = "FHM_SPECIES")
#  write.csv(res_CD, "output/res_CD_FHM.csv")  
#  return(res_CD)
# }

#####################################################################################################
