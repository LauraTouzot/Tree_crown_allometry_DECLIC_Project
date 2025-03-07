# Data paper

crown_data_paper <- function(){
  
  data_m <- read.csv("data/datapaper/data_m.csv")
  data_m <- data_m[!data_m$code_status %in% c(9991, 9990),] # removing dead and harvested trees
  
  data_c <- read.csv("data/datapaper/data_c.csv")
  data_p <- read.csv("data/datapaper/data_p.csv")
  data_p <- data_p[,c("plot_id", "lat", "long", "area")]
  
  species_code <- read.csv("data/datapaper/species_code.csv")
  names(species_code) <- c("code_species", "sp")

  data_m <-data_m %>% filter(!duplicated(tree_id))
  df <- left_join(data_m,data_c, by = "tree_id")
  df <- left_join(df, species_code, by = "code_species")
  names(data_p) <- c("plot_id.x", "lat", "long", "area")
  df <- left_join(df, data_p, by = "plot_id.x")
  df <- df %>% rename(HT_m = h_tot, DBH_cm = dbh)
  df$A1 <- (df$crown_r1 * df$crown_r2 * pi) / 4
  df$A2 <- (df$crown_r2 * df$crown_r3 * pi) / 4
  df$A3 <- (df$crown_r3 * df$crown_r4 * pi) / 4
  df$A4 <- (df$crown_r4 * df$crown_r1 * pi) / 4
  df$C_diam_m <- 2 * sqrt((df$A1 + df$A2 + df$A3 + df$A4) / pi)
  df <- df %>% mutate(C_height_m = rowMeans(dplyr::select(.,crown_h1, crown_h2, crown_h3, crown_h4), na.rm = TRUE),
                      C_depth_m = HT_m - C_height_m,
                      CR = C_depth_m/HT_m)
  
  # removing data on broken and leaned trees from the allometry database 
  df[df$code_status %in% 8881:8882,]$HT_m <- NA
  df[df$code_status %in% 8881:8882,]$C_diam_m <- NA
  df[df$code_status %in% 8881:8882,]$C_depth_m <- NA
  df[df$code_status %in% 8881:8882,]$CR <- NA
  
  # checking minimum DBH
  df <- df[df$DBH_cm >= 7.5,]
  
  # adding weight and continent
  df$W <- 1/df$area
  df$continent <- "E_U"
  
  df <- df %>% ungroup()
  
  # for each plot, select year with the highest amount of data
  summary_perlocation <- df %>% dplyr::select(plot_id.x, year, HT_m, CR, C_diam_m, C_depth_m) %>%
    gather(key = variable, value = value, HT_m, CR, C_diam_m, C_depth_m) %>%
    filter(!is.na(value)) %>%
    group_by(plot_id.x, year, variable) %>%
    summarize(n = n()) %>%
    spread(key = variable, value = n)
  
  summary_perlocation$nobs <- as.numeric(apply(summary_perlocation[,c(3:6)], 1, sum, na.rm = TRUE))
  
  to_keep <- summary_perlocation %>% group_by(plot_id.x) %>% 
    arrange(desc(nobs), .by_group = TRUE) %>%
    slice(stay = unique(max(nobs)),.by_group = TRUE) 
  
  to_keep <- to_keep[!duplicated(to_keep$plot_id.x),]
  to_keep$staying_IDs <- paste(to_keep$plot_id.x, to_keep$year, sep = "_")
  
  df <- df %>% mutate(location_ID = paste(plot_id.x, year, sep = "_"))
  
  df <- df[df$location_ID %in% to_keep$staying_IDs,]
  
  df$location_ID <- paste(df$location_ID, "Fuhr2017", sep = "_")
  
  # plots to check unit within the dataset
  p1 <-  ggplot(df[!is.na(df$HT_m), ], aes(x = DBH_cm, y = HT_m)) + geom_point(alpha = 0.1) 
  p2 <-  ggplot(df[!is.na(df$CR), ], aes(x = DBH_cm, y = CR)) + geom_point(alpha = 0.1) 
  p3 <-  ggplot(df[!is.na(df$C_depth_m), ], aes(x = DBH_cm, y = C_depth_m)) + geom_point(alpha = 0.1) 
  p4 <-  ggplot(df[!is.na(df$C_diam_m), ], aes(x = DBH_cm, y = C_diam_m)) + geom_point(alpha = 0.1) 
  png("figures/unit_check/unit_check_Fuhr.png", width = 680, height = 680)
  multiplot(p1, p2, p3, p4, cols =  2)
  dev.off()
  
  return(df)
}


############################### NOT USED IN THIS VERSION OF THE CODE ###############################

# crown_sp_data_paper <- function(df){
#   library(dplyr)
#   res_crown <-   df %>% group_by(sp) %>% summarise(CD_obs = sum(!is.na(C_diam_m)),
#                                                              CR_obs = sum(!is.na(C_depth_m)),
#                                                              HT_obs = sum(!is.na(HT_m)))
#   return(res_crown)
# }

#####################################################################################################

