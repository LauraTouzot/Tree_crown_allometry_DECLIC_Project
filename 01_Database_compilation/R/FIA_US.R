# # Read FIA data
# # Get a unique tree id df$TreeID <- paste(df$STATECD, df$UNITCD, df$COUNTYCD,df$PLOT , df$SUBP, df$TREE)
# # STATECD, INVYR, UNITCD, COUNTYCD
# 
# get_FIA_data_CR <- function(){
# Sys.setenv('R_MAX_VSIZE'=320000000000)
# 
#   df <- ls_fetch_usa(state = "all")
# 
#   df_c_h <- df$tree 
#   rm(df)
#   gc()
#   df_c_h <- df_c_h %>% dplyr::mutate(TreeID = paste(STATECD, UNITCD, COUNTYCD, PLOT, SUBP, TREE, sep = "_")) %>% 
#                        dplyr::select(PREV_TRE_CN, INVYR, STATECD, UNITCD,COUNTYCD,
#                               PLOT, SUBP, TREE, CONDID, AZIMUTH, DIST,
#                               PREVCOND, STATUSCD, SPCD, SPGRPCD, DIA, DIAHTCD,
#                               HT, HTCD, ACTUALHT, TREECLCD, CR, CCLCD, TREEGRCD,
#                               AGENTCD, DAMLOC1, DAMTYP1, DAMSEV1, DAMLOC2,
#                               DAMTYP2, DAMSEV2, DIACHECK, MORTYR, UNCRCD,
#                               CPOSCD, CLIGHTCD, CVIGORCD, CDENCD, CDIEBKCD,
#                               TRANSCD, TREEHISTCD, DIACALC, BHAGE, TOTAGE,
#                               MORTCD, HTDMP, RECONCILECD, PREVDIA, CYCLE,
#                               SUBCYCLE, TreeID, PLT_CN, TPA_UNADJ) %>% 
#                               filter(INVYR != 9999 & INVYR >= 2002 & STATUSCD == 1)
#   
#   species_code <- read.csv("data/FIA/REF_SPECIES.csv") %>%
#         dplyr::select("SPCD", "GENUS", "SPECIES")
#   
#   df_c_h <- left_join(df_c_h, species_code, by = "SPCD")
#   
#   df_c_h$DIA[df_c_h$DIAHTCD != 1 & !is.na(df_c_h$DIAHTCD)] <- NA
#   df_c_h$DBH_cm <- conv_unit(df_c_h$DIA, "inch", "cm")
#   
#   df_c_h$HT[df_c_h$HTCD != 1 & !is.na(df_c_h$HTCD)] <- NA  
#   df_c_h$HT_m_est <- conv_unit(df_c_h$HT, "feet", "m")
#   df_c_h$HT_m <- conv_unit(df_c_h$ACTUALHT, "feet", "m")
#   
#   df_c_h$CR <- df_c_h$CR/100 
#   df_c_h$C_depth_m <- df_c_h$HT * df_c_h$CR # crown depth estimated from total tree height and crown ratio
#   df_c_h <- df_c_h %>% dplyr::mutate(sp = paste(GENUS, SPECIES))
#   
#   return(df_c_h)
#   
#   ## res_HT <-   df$tree %>% group_by(SPCD) %>% summarise(HT_obs = sum(!is.na(HT)))
#   ## res_HT <- left_join(res_HT, species_code, by = "SPCD")
#   ## res_CR <-   df$tree %>% group_by(SPCD) %>% summarise(CR_obs = sum(!is.na(CR)))
#   ## res_CR <- left_join(res_CR, species_code, by = "SPCD")
#   ## write.csv(res_CR, "output/res_CR_FIA.csv")
#   ## write.csv(res_HT, "output/res_HT_FIA.csv")
#   
#   # no minimum DBH in FIA data
# }
# 
# 
# merge_FIA_tree_plot <- function(FIA_data_tree) {
#   
#   Sys.setenv('R_MAX_VSIZE'=320000000000)
#   
#   get_df_plot <- ls_fetch_usa(state = "all", what = "plot")
#   df_plot_complete <- get_df_plot$plot 
#   df_plot_clean <- df_plot_complete %>% dplyr::select(CN, LAT, LON, STATECD, DESIGNCD, UNITCD, COUNTYCD, PLOT)
#   df_plot_clean <- df_plot_clean %>% mutate(location_ID = paste(STATECD, UNITCD, COUNTYCD, PLOT, sep = "_"))
#   
#   FIA_tree_plot <- left_join(FIA_data_tree, df_plot_clean, by = c("PLT_CN" = "CN"))
#   
#   # for each plot, select year with the highest amount of data
#   summary_perlocation <- FIA_tree_plot %>% dplyr::select(location_ID, INVYR, HT_m, CR) %>%
#     gather(key = variable, value = value, HT_m, CR) %>%
#     filter(!is.na(value)) %>%
#     group_by(location_ID, INVYR, variable) %>%
#     summarize(n = n()) %>%
#     spread(key = variable, value = n)
#   
#   summary_perlocation$nobs <- as.numeric(apply(summary_perlocation[,c(3:4)], 1, sum, na.rm = TRUE))
#   
#   to_keep <- summary_perlocation %>% group_by(location_ID) %>% 
#     arrange(desc(nobs), .by_group = TRUE) %>%
#     slice(stay = unique(max(nobs)),.by_group = TRUE) 
#   
#   to_keep <- to_keep[!duplicated(to_keep$location_ID),]
#   to_keep$staying_IDs <- paste(to_keep$location_ID, to_keep$INVYR, sep = "_")
#   
#   FIA_tree_plot <- FIA_tree_plot %>% mutate(location_ID = paste(location_ID, INVYR, sep = "_"))
#   FIA_tree_plot <- FIA_tree_plot[FIA_tree_plot$location_ID %in% to_keep$staying_IDs,]
#   
#   # convert weight to get a coefficient in hectare
#   FIA_data_tree_plot <- FIA_tree_plot %>% mutate(W = 1/((1/TPA_UNADJ)*0.404686)) 
#   
#   FIA_data_tree_plot <- FIA_data_tree_plot %>% ungroup()
#   
#   # plots to check unit within the dataset
#   p1 <-  ggplot(FIA_data_tree_plot[!is.na(FIA_data_tree_plot$HT_m), ], aes(x = DBH_cm, y = HT_m)) + geom_point(alpha = 0.1) 
#   p2 <-  ggplot(FIA_data_tree_plot[!is.na(FIA_data_tree_plot$CR), ], aes(x = DBH_cm, y = CR)) + geom_point(alpha = 0.1) 
#   p3 <-  ggplot(FIA_data_tree_plot[!is.na(FIA_data_tree_plot$C_depth_m), ], aes(x = DBH_cm, y = C_depth_m)) + geom_point(alpha = 0.1) 
#   png("figures/unit_check/unit_check_FIA.png", width = 680, height = 680)
#   multiplot(p1, p2, p3, cols =  2)
#   dev.off()
#   
#   return(FIA_data_tree_plot)
# }
# 
# 
# 
# get_FIA_coord <- function(FIA_data_plot){
# 
#   coordinates <- FIA_data_plot %>% dplyr::select(TreeID, STATECD.x, DBH_cm, HT_m, C_depth_m, sp, LAT, LON, location_ID, W, CR) %>%
#                                    rename(STATECD = STATECD.x, latitude = LAT, longitude = LON)
# 
#   Pacific_Islandes <- c(15, 60, 64, 66, 68, 69, 70)
#   
#   States_data <- coordinates[!(coordinates$STATECD %in% Pacific_Islandes),]
#   States_data_NA <- States_data[is.na(States_data$latitude) & is.na(States_data$longitude),]
#   States_data <- States_data[!(is.na(States_data$latitude)) & !(is.na(States_data$longitude)),]
#   PI_data <- coordinates[coordinates$STATECD %in% Pacific_Islandes,]
#   
#   extract_coordinates_states <- States_data[,c("longitude", "latitude")]
#   coordinates(extract_coordinates_states) <- ~longitude + latitude  
#   proj4string(extract_coordinates_states) <- CRS("+init=epsg:4269")
#   extract_coordinates_states <- spTransform(extract_coordinates_states, CRS("+init=epsg:4326"))
#   coordinatesWGS <- extract_coordinates_states@coords
#   coordinatesWGS <- as.data.frame(coordinatesWGS)
#   
#   States_data$latitude <- coordinatesWGS$latitude
#   States_data$longitude <- coordinatesWGS$longitude
#   
#   clean_coordinates <- dplyr::bind_rows(PI_data, States_data, States_data_NA) # no duplicates left in location coordinates
# 
#   clean_coordinates <- clean_coordinates %>% ungroup()
#   
#   clean_coordinates$latitude_plot <- clean_coordinates$latitude
#   clean_coordinates$longitude_plot <- clean_coordinates$longitude
#   clean_coordinates$latitude_tree <- NA
#   clean_coordinates$longitude_tree <- NA
#   clean_coordinates$continent <- continent <- "N_A"
#   
#   return(clean_coordinates)
# }








############################### NOT USED IN THIS VERSION OF THE CODE ###############################

# FIA_data_HT_per_sp <- function(df){
#   library(dplyr)
#   res_HT <-   df %>% group_by(SPCD) %>% summarise(HT_obs = sum(!is.na(HT)))
#   write.csv(res_HT, "output/res_HT_FIA.csv")
#   res_HT
# }
# 
# FIA_data_CR_per_sp <- function(df){
#   library(dplyr)
#   res_CR <-   df %>% group_by(SPCD) %>% summarise(CR_obs = sum(!is.na(CR)))
#   write.csv(res_CR, "output/res_CR_FIA.csv")
#   res_CR
# }

#####################################################################################################

