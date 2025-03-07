### Computing final and complete version of the allometry database

complete_allometry <- function(allometry_supp_variables, location_variables_checked) {
  
  df <- left_join(allometry_supp_variables, location_variables_checked, by = "location_ID")
  df <- df %>% ungroup()
  
  return(df)
  
}