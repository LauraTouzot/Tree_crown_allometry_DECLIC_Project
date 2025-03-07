### Last cleaning of the database to be extra sure no mistakes have been left behind
## remove observations if
# dbh = NA (or 0 even though it is not supposed to still be in the dataset)
# total tree height > crown depth
# total tree height > 1.3
# continent different than Europe or North America

last_cleaning <- function(all_crown) {
  
  # keeping only data located in Europe and North America
  all_crown <- all_crown[all_crown$continent %in% c("N_A", "E_U"),]
  
  id <- as.data.frame(c(1:dim(all_crown)[1]))
  colnames(id) <- "id"
  df <- cbind(id, all_crown)
  df <- df %>% filter(DBH_cm > 0 & !is.na(DBH_cm))
  
  # weird thing when removing HT_m = 0, that's why I'm using this uggly code
  df_b <- df %>% filter(HT_m == 0) 
  df_b$HT_m <- NA
  df <- df %>% filter(!id %in% df_b$id)
  
  df_c <- df %>% filter(C_depth_m == 0) 
  df_c$C_depth_m <- NA
  df <- df %>% filter(!id %in% df_c$id)
 
  df_d <- df %>% filter(C_diam_m == 0) 
  df_d$C_diam_m <- NA
  df <- df %>% filter(!id %in% df_d$id)

  df <- rbind(df, df_b, df_c, df_d)
  
  return(df)
  
}






