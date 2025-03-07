# READ ICP data
read_ICP <- function(){
  
  df_ipm <- read.csv("data/ICP/fmd_gr_ipm.csv", sep = ";")
  df_inv <- read.csv("data/ICP/fmd_gr_inv.csv", sep = ";")
  df_pli <- read.csv("data/ICP/fmd_gr_pli.csv", sep = ";")
  species <- read.csv("data/ICP/SPECIES.csv", sep = ",")
  d_removal <- read.csv("data/ICP/d_removal_mortality_ccgr.csv", sep = ",")
  continent <- read.csv("data/ICP/COUNTRY.csv", sep = ";")

  # remove data with height smaller than crown base height
  df_ipm[df_ipm$height < df_ipm$crown_base_height & !is.na(df_ipm$crown_base_height) & !is.na(df_ipm$height), c("height", "crown_base_height" )] <- NA
  df_ipm[df_ipm$diameter == 0 & !is.na(df_ipm$diameter), ]$diameter <- NA
  
  # no minimum DBH
  df_ipm$PlotID <- paste(df_ipm$code_country, df_ipm$code_plot, sep = "_")
  df_ipm$TreeID <- paste(df_ipm$code_country, df_ipm$code_plot, df_ipm$tree_number, sep = "_")
  df_inv$PlotID <- paste(df_inv$code_country, df_inv$code_plot, sep = "_")
  
  # set a value of NA for height base_height and cronw width if bad quality codebased on removal code
  df_ipm$height[df_ipm$code_height_qc %in% c(3, 5, 9) & !is.na(df_ipm$code_height_qc)] <- NA
  df_ipm$crown_base_height[df_ipm$code_base_height_qc %in% c(3, 5, 9) & !is.na(df_ipm$code_base_height_qc)] <- NA
  df_ipm$crown_width[df_ipm$code_crown_width_qc %in% c(3, 5, 9) & !is.na(df_ipm$code_crown_width_qc)] <- NA
  # remove harvested and dead trees
  df_ipm <- df_ipm %>% filter(!(code_removal %in% c(11:19, 31:49) & is.na(code_removal)))
  #4 tree not in sample any more 11:19, 31:49
  
  df_ipm[df_ipm$code_removal %in% 11:49 & is.na(df_ipm$code_removal),  c("height", "crown_base_height", "crown_width")] <- NA
  #tree not used for crown allo 21:29
  
  df_ipm <- left_join(df_ipm, species, by = "code_tree_species")
  
  # THIS IS NOT WORKING FOR FRANCE
  df_inv$PlotID_b <- paste(df_inv$code_country, df_inv$code_plot, df_inv$survey_year, sep = "_")
  df_pli$PlotID_b <- paste(df_pli$code_country, df_pli$code_plot, df_pli$survey_year, sep = "_")
  df_ipm$PlotID_b <- paste(df_ipm$code_country, df_ipm$code_plot, df_ipm$survey_year, sep = "_")
  df_pli <- df_pli[!duplicated(df_pli$PlotID_b ) & !is.na(df_pli$longitude), ]
  
  df_inv <- df_inv[!df_inv$PlotID_b %in% names(table( df_inv$PlotID_b ))[table( df_inv$PlotID_b )>1], ]  
  df_pli$latitude_degree <- substr(df_pli$latitude, 1, 2)
  df_pli$latitude_minute <- substr(df_pli$latitude, 3, 4)
  df_pli$latitude_second <- substr(df_pli$latitude, 5, 6)
  df_pli$lon_sign <- "+"
  df_pli$lon_sign[grepl("-", df_pli$longitude, fixed = TRUE)  ] <- "-"
  
  df_pli$longitude2 <- df_pli$longitude
  df_pli$longitude <- gsub("-", "", df_pli$longitude, fixed = TRUE)
  df_pli$longitude_degree <- substr(df_pli$longitude, 1, 2)
  df_pli$longitude_minute <- substr(df_pli$longitude, 3, 4)
  df_pli$longitude_second <- substr(df_pli$longitude, 5, 6)
  # 1
  df_pli$longitude_degree[nchar(df_pli$longitude)==1] <- "00"
  df_pli$longitude_minute[nchar(df_pli$longitude)==1] <- "00"
  df_pli$longitude_second[nchar(df_pli$longitude)==1] <- substr(df_pli$longitude[nchar(df_pli$longitude)==1], 1, 1)
  # 3
  df_pli$longitude_degree[nchar(df_pli$longitude)==3] <- "00"
  df_pli$longitude_minute[nchar(df_pli$longitude)==3] <- paste0("0", substr(df_pli$longitude[nchar(df_pli$longitude)==3], 1, 1))
  df_pli$longitude_second[nchar(df_pli$longitude)==3] <- substr(df_pli$longitude[nchar(df_pli$longitude)==3], 2, 3)
  # 4
  df_pli$longitude_degree[nchar(df_pli$longitude)==4] <- "00"
  df_pli$longitude_minute[nchar(df_pli$longitude)==4] <- substr(df_pli$longitude[nchar(df_pli$longitude)==4], 1, 2)
  df_pli$longitude_second[nchar(df_pli$longitude)==4] <- substr(df_pli$longitude[nchar(df_pli$longitude)==4], 3, 4)
  # 5
  df_pli$longitude_degree[nchar(df_pli$longitude)==5] <- paste0("0", substr(df_pli$longitude[nchar(df_pli$longitude)==5], 1, 1))
  df_pli$longitude_minute[nchar(df_pli$longitude)==5] <- substr(df_pli$longitude[nchar(df_pli$longitude)==5], 2, 3)
  df_pli$longitude_second[nchar(df_pli$longitude)==5] <- substr(df_pli$longitude[nchar(df_pli$longitude)==5], 4, 5)
  df_pli[, c("latitude_degree", "latitude_minute", "latitude_second")] <- lapply(df_pli[, c("latitude_degree", "latitude_minute", "latitude_second")], as.numeric )
  df_pli[, c("longitude_degree", "longitude_minute", "longitude_second")] <- lapply(df_pli[, c("longitude_degree", "longitude_minute", "longitude_second")], as.numeric )
  df_pli$lon <- df_pli$longitude_degree +df_pli$longitude_minute/60 + df_pli$longitude_second/3600
  df_pli$lat <- df_pli$latitude_degree +df_pli$latitude_minute/60 + df_pli$latitude_second/3600
  df_pli$lon[df_pli$lon_sign=="-"] <- -1*df_pli$lon[df_pli$lon_sign=="-"] 
  

  df_ipm <- left_join(df_ipm, continent, by = "code_country")
 
  df_ipm <- left_join(df_ipm, df_pli %>% dplyr::select(PlotID_b, sample_plot_size, lon, lat), by = "PlotID_b")


  df <- df_ipm %>% rename(HT_m = height, DBH_cm = diameter, C_diam_m = crown_width, latitude_plot = lat, longitude_plot = lon) %>%
    mutate(C_depth_m = HT_m - crown_base_height,
          CR = C_depth_m/HT_m,
          data = "ICP", W = 1/sample_plot_size,year = survey_year) %>% 
   dplyr::select(sp,DBH_cm, C_diam_m, CR, C_depth_m, HT_m, data, PlotID, continent, latitude_plot, longitude_plot, W, code_country, year,PlotID_b,PlotID, country) %>%
    rename(PlotIDYear =  PlotID_b) %>%
   relocate(sp,DBH_cm, C_diam_m, CR, C_depth_m, HT_m, data, PlotID, continent, latitude_plot, longitude_plot, W) %>% ungroup()
  df$W[df$year < 2000] <- NA
  df$C_diam_m[df$C_diam_m > 27] <- NA
  df$DBH_cm[df$DBH_cm < 0.1] <- NA
  df$DBH_cm[df$DBH_cm > 300] <- NA
  
  dd <- df %>% group_by(PlotID, year) %>% summarise(W = unique(W), lon = unique(longitude_plot), lat = unique(latitude_plot))
  df <- df %>% select(sp, DBH_cm, C_diam_m, CR, C_depth_m, HT_m, PlotID, year, continent)
  
  res <- left_join(df, dd, by = c("PlotID", "year"))
  
  res$location_ID <- paste(res$PlotID, res$year, "ICP", sep = "_")
  res$latitude_plot = res$lat
  res$longitude_plot = res$lon
  res$latitude_tree <- NA
  res$longitude_tree <- NA
  
  # plots to check unit within the dataset
  p1 <-  ggplot(res[!is.na(res$HT_m), ], aes(x = DBH_cm, y = HT_m)) + geom_point(alpha = 0.1) 
  p2 <-  ggplot(res[!is.na(res$CR), ], aes(x = DBH_cm, y = CR)) + geom_point(alpha = 0.1) 
  p3 <-  ggplot(res[!is.na(res$C_depth_m), ], aes(x = DBH_cm, y = C_depth_m)) + geom_point(alpha = 0.1) 
  p4 <-  ggplot(res[!is.na(res$C_diam_m), ], aes(x = DBH_cm, y = C_diam_m)) + geom_point(alpha = 0.1) 
  png("figures/unit_check/unit_check_ICP.png", width = 680, height = 680)
  multiplot(p1, p2, p3, p4, cols =  2)
  dev.off()
  
  return(res)
  
#  res <- df %>% group_by(PlotID, year) %>% summarise(Ntree = length(DBH_cm),
#                                                    Ntree_ha = length(DBH_cm) * unique(W),
#                                                    BA = sum((DBH_cm/200)^2*pi, na.rm = TRUE) * unique(W),
#                                                    country = unique(country),
#                                                    year = unique(year),
#                                                    lon = unique(longitude_plot),
#                                                    lat = unique(latitude_plot),
#                                                    W = unique(W))
 
# ggplot(res, aes(year, BA, color = PlotID)) + geom_point()+ geom_line() +
#   facet_wrap(~country)+ theme(legend.position = "none")
# 
# ggplot(df[!duplicated(df$PlotIDYear), ], aes(longitude_plot, latitude_plot, color = factor(country), shape =factor(code_country))) + geom_point()+ 
#   theme(legend.position = "none")
# 
# ggplot(res, aes(x=factor(country), y=BA)) + 
#   geom_boxplot()
# 
# ggplot(df, aes(DBH_cm, HT_m, color = PlotID)) + geom_point()+ 
#   facet_wrap(~country)+ theme(legend.position = "none")
# ggplot(df, aes(DBH_cm, C_diam_m, color = PlotID)) + geom_point()+ 
#   facet_wrap(~country)+ theme(legend.position = "none")
# ggplot(df, aes(HT_m, C_depth_m, color = PlotID)) + geom_point()+ 
#   facet_wrap(~factor(country))+ theme(legend.position = "none")

# plot  A unique number (per country) given to the permanent plot during the selection or installation
# INV code_country	ID of the country
# INV code_plot	Observation plot number
# IPM	tree_number	Tree number
# IPM	code_tree_species	Coded tree species
# IPM	diameter	Diameter at breast height [dbh] cm
# IPM	height	Tree/shoot height (vertical distance between the highest point of the crown and the ground surface) m
# IPM	crown_base_height	Height to the crown base (vertical distance between the crown base and the ground level) m
# IPM	crown_width	Crown width (averaged) m
# IPM	code_removal	Removal/Mortality - status of tree
# IPM	code_diameter_qc	Diameter Quality code
# IPM	code_height_qc	Height quality code
# IPM	code_base_height_qc	Height to crown quality code
# IPM	code_crown_width_qc	Crown width quality code

# Description: Latitude in Degrees Minutes Seconds (WGS84)
# Remarks: E.g. ‘+505852’ = 50 Degree, 58 Minutes, 52 Seconds

# check if multiple height per tree

#df_pli$PlotID <- paste(df_pli$code_country, df_pli$code_plot, sep = "_")
#df_inv$PlotID <- paste(df_inv$code_country, df_inv$code_plot, sep = "_")

}
