# New Canada NFI

canada_NFI_crown_data <- function(){
  
  df_small_plots <- read.csv("data/CanadaNew/NFI_national_GP_rounded/all_gp_trees_approx/all_gp_stp_tree.csv")
  df_small_plots_h <- read.csv("data/CanadaNew/NFI_national_GP_rounded/all_gp_trees_approx/all_gp_stp_header.csv")
  df_large_plots <- read.csv("data/CanadaNew/NFI_national_GP_rounded/all_gp_trees_approx/all_gp_ltp_tree.csv")
  df_large_plots <- df_large_plots[df_large_plots$lgtree_status != "DS",] # removing dead trees
  df_large_plots_h <- read.csv("data/CanadaNew/NFI_national_GP_rounded/all_gp_trees_approx/all_gp_ltp_header.csv")
  df_site_info <- read.csv("data/CanadaNew/NFI_national_GP_rounded/all_gp_trees_approx/all_gp_site_info.csv")
  species <- read.csv("data/CanadaNew/NFI_genus_sp_latin.csv")
  
  species <- species %>% mutate(sp_code = paste0(genus, species, variety)) %>%
    select(sp_code, latin_name)
  
  # convert data with good utm zone
  df_site_info$lon <- NA
  df_site_info$lat <- NA
  
  for (i in unique(df_site_info$utm_zone)) {
    
    df <- df_site_info[df_site_info$utm_zone == i, ]
    coordinates(df) <- c("utm_e",  "utm_n")  
    proj4string(df) <- CRS(paste0("+proj=utm +zone=", i, " +datum=WGS84  +units=m")) # utm good zone
    df <- spTransform(df,  CRS("+init=epsg:4326"))
    df_site_info[df_site_info$utm_zone == i, ]$lon <- df$utm_e
    df_site_info[df_site_info$utm_zone == i, ]$lat <- df$utm_n
    
  }
  
  coordinates <- df_site_info[, c("nfi_plot",  "meas_date", "lon", "lat")]
  
  df_large_plots$year <- substr(df_large_plots$meas_date, 1, 4)
  df_large_plots$height[df_large_plots$meas_est_dbh %in% c("C") | df_large_plots$height <= 0] <- NA
  df_large_plots$dbh[ df_large_plots$dbh <= 0] <- NA
  df_large_plots$crown_base[ df_large_plots$crown_base <= 0] <- NA
  df_large_plots$crown_length[ df_large_plots$crown_length <= 0] <- NA

  names(df_large_plots)
  names(df_small_plots_h)
  
  df_large_plots <- left_join(df_large_plots, df_large_plots_h[, c("nfi_plot", "meas_date", "meas_plot_size")], by = c("nfi_plot", "meas_date"))
  df <- df_large_plots %>%
    dplyr::select(nfi_plot, loc_id, meas_date, tree_num, lgtree_genus, lgtree_species,  lgtree_variety,
           dbh, height, crown_base, stem_cond, year, stem_cond, meas_plot_size) %>%
    mutate(sp_code = paste0(lgtree_genus, lgtree_species, lgtree_variety),
           tree_id = paste(nfi_plot, tree_num))
  
  df_small_plots$smtree_ht[df_small_plots$smtree_measest_ht %in% c("C") | df_small_plots$smtree_ht <= 0] <- NA
  df_small_plots$smtree_dbh[ df_small_plots$smtree_dbh <= 0] <- NA
  df_small_plots$year <- substr(df_small_plots$meas_date, 1, 4)
  df_small_plots <- left_join(df_small_plots, df_small_plots_h[, c("nfi_plot", "meas_date", "meas_plot_size")], by = c("nfi_plot", "meas_date"))
  
  # stem_cond
  df2 <- df_small_plots %>%
    dplyr::select(nfi_plot, loc_id, meas_date, smtree_num, smtree_genus, smtree_species, 
           smtree_variety,smtree_dbh, smtree_ht,  stem_cond, year, meas_plot_size) %>%
    mutate(sp_code = paste0(smtree_genus, smtree_species, smtree_variety),
           tree_id = paste(nfi_plot, smtree_num))

  df <- left_join(df, species, by =c("sp_code"))
  df2 <- left_join(df2, species, by =c("sp_code"))
  
  df$crown_base[df$HT_m < df$crown_base] <- NA
 
  df <- df %>% 
    rename(DBH_cm = dbh, HT_m = height,sp = latin_name, TreeID = tree_id) %>% 
    mutate(C_depth_m = HT_m - crown_base,
           CR = C_depth_m/ HT_m) %>%
    filter(HT_m > crown_base)
  
 # removing allometry data for broken trees  
 df$HT_m[df$stem_cond == "B"] <- NA
 df$C_depth_m[df$stem_cond == "B"] <- NA
 df$CR[df$stem_cond == "B"] <- NA
 
 df <- df %>% dplyr::select(nfi_plot, meas_date, TreeID, sp, DBH_cm, HT_m, C_depth_m, CR, meas_plot_size, stem_cond, year) %>%
   mutate(CR = C_depth_m/ HT_m, W = 1/meas_plot_size)
 
 df <- df %>% dplyr::select(nfi_plot, meas_date, TreeID, sp, DBH_cm, HT_m, C_depth_m, CR, stem_cond, year, W)

 df2 <- df2 %>% 
   rename(DBH_cm = smtree_dbh, HT_m = smtree_ht,sp = latin_name, TreeID = tree_id) %>% 
   mutate(C_depth_m = NA,
          CR = NA,
          W = 1/meas_plot_size)
 df2 <- df2 %>% dplyr::select(nfi_plot, meas_date, TreeID, sp, DBH_cm, HT_m, C_depth_m, CR, stem_cond, year, W)
 
 df_all <- rbind(df, df2)
 
 Canada_complete <- left_join(df_all, coordinates, by = c("nfi_plot", "meas_date"))
 
 Canada_complete$latitude_tree <- NA
 Canada_complete$longitude_tree <- NA
 Canada_complete$longitude_plot <- Canada_complete$lon
 Canada_complete$latitude_plot <- Canada_complete$lat
 
 # for each plot, select year with the highest amount of data
 summary_perlocation <- Canada_complete %>% dplyr::select(nfi_plot, year, HT_m, CR) %>%
   gather(key = variable, value = value, HT_m, CR) %>%
   filter(!is.na(value)) %>%
   group_by(nfi_plot, year, variable) %>%
   summarize(n = n()) %>%
   spread(key = variable, value = n)
 
 summary_perlocation$nobs <- as.numeric(apply(summary_perlocation[,c(3:4)], 1, sum, na.rm = TRUE))

 to_keep <- summary_perlocation %>% group_by(nfi_plot) %>% 
   arrange(desc(nobs), .by_group = TRUE) %>%
   slice(stay = unique(max(nobs)),.by_group = TRUE) 
 
 to_keep <- to_keep[!duplicated(to_keep$nfi_plot),]
 to_keep$staying_IDs <- paste(to_keep$nfi_plot, to_keep$year, sep = "_")
 
 Canada_complete <- Canada_complete %>% mutate(location_ID = paste(nfi_plot, year, sep = "_"))
 
 Canada_complete <- Canada_complete[Canada_complete$location_ID %in% to_keep$staying_IDs,]
 
 Canada_complete$location_ID <- paste(Canada_complete$nfi_plot, "Canada", sep = "_")
 Canada_complete$continent <- "N_A"
  
 Canada_complete <- Canada_complete %>% ungroup()
 
 # plots to check unit within the dataset
 p1 <-  ggplot(Canada_complete[!is.na(Canada_complete$HT_m), ], aes(x = DBH_cm, y = HT_m)) + geom_point(alpha = 0.1) 
 p2 <-  ggplot(Canada_complete[!is.na(Canada_complete$CR), ], aes(x = DBH_cm, y = CR)) + geom_point(alpha = 0.1) 
 p3 <-  ggplot(Canada_complete[!is.na(Canada_complete$C_depth_m), ], aes(x = DBH_cm, y = C_depth_m)) + geom_point(alpha = 0.1) 
 png("figures/unit_check/unit_check_Canada.png", width = 680, height = 680)
 multiplot(p1, p2, p3, cols =  2)
 dev.off()
  
 return(Canada_complete)
  
}
