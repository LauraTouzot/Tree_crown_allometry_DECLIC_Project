# Quebec data height

quebec_NFI_crown_data <- function(){

  df <- read.csv("data/Quebec/data.quebec.csv")
  species <- read.csv("data/Quebec/Essences.csv", sep ="\t", skip = 2)
  df <- left_join(df, species[ , c("CODE", "Latin_name")], by = c("ESSENCE" = "CODE"))
  
  # Hauteur en dcm HAUTEUR.1 et HAUTEUR.R
  # DIAM mm DHPMM.1 et DHPMM.2
  # removing dead and harvested trees (ETAT.1)
    
  df <- df %>% filter(DHPMM.1> 90 & !is.na(DHPMM.1) & !(ETAT.1%in% c(12, 32, 42, 14, 26, 44, 54, 52))) %>%
               mutate(DBH_cm = DHPMM.1/10, HT_m = HAUTEUR.1/10) %>% rename(sp = Latin_name)
  df <- df %>% select(sp, DBH_cm, HT_m, ID_PEP, TreeID, YEAR.1)

  # for each plot, select year with the highest amount of data
  summary_perlocation <- df %>% dplyr::select(ID_PEP, YEAR.1, HT_m) %>%
    gather(key = variable, value = value, HT_m) %>%
    filter(!is.na(value)) %>%
    group_by(ID_PEP, YEAR.1, variable) %>%
    summarize(n = n()) %>%
    spread(key = variable, value = n)

  to_keep <- summary_perlocation %>% group_by(ID_PEP) %>% 
    arrange(desc(HT_m), .by_group = TRUE) %>%
    slice(stay = unique(max(HT_m)),.by_group = TRUE) 
  
  to_keep <- to_keep[!duplicated(to_keep$ID_PEP),]
  to_keep$staying_IDs <- paste(to_keep$ID_PEP, to_keep$YEAR.1, sep = "_")
  
  df_b <- df %>% mutate(location_ID = paste(ID_PEP, YEAR.1, sep = "_"))
  
  df_b <- df_b[df_b$location_ID %in% to_keep$staying_IDs,]
  
  coordinates <- read.csv("data/Quebec/localis.csv")
  infogen <- read.csv("data/Quebec/infogen.csv")
  infogen <- infogen %>% dplyr::select(ID_PEP, ID_PEP_MES)
  
  coordinates <- left_join(coordinates, infogen, by = "ID_PEP_MES")
  coordinates <- coordinates %>% dplyr::select(ID_PEP, ID_PEP_MES, LATITUDE, LONGITUDE)

  df_complete <- left_join(df_b, subset(coordinates, subset =! duplicated(coordinates[["ID_PEP"]])), by = "ID_PEP")
  
  df_complete$latitude_plot <- df_complete$LATITUDE
  df_complete$longitude_plot <- df_complete$LONGITUDE
  df_complete$latitude_tree <- NA
  df_complete$longitude_tree <- NA
  df_complete$location_ID <- paste(df_complete$ID_PEP, "QuebecNFI", sep = "_")
  df_complete$continent <- "N_A"
  
  df_complete <- df_complete %>% mutate(W = ifelse(DBH_cm > 31, 1/(pi * 14.19^2 *10^-4), 1/(pi * 11.28^2 *10^-4)))
  
  # plots to check unit within the dataset
  p1 <-  ggplot(df_complete[!is.na(df_complete$HT_m), ], aes(x = DBH_cm, y = HT_m)) + geom_point(alpha = 0.1) 
  png("figures/unit_check/unit_check_Quebec.png", width = 340, height = 340)
  p1
  dev.off()
                     
  return(df_complete)
}

