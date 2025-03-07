## Comparing NFI vs. allometry databases

binding_databases <- function(NFI_data, summary_species) {

  df <- left_join(NFI_data, summary_species, by = c("checked_name" = "sp"))
  df <- df %>% dplyr::select(checked_name, continent, nplot, ntree, nplot_crown, ntree_crown, nobs_HT, nobs_Cdiam, nobs_Cdepth, nobs_CR) %>% 
               rename(sp = checked_name, nplot_NFI = nplot, ntree_NFI = ntree)
  
  return(df)
}