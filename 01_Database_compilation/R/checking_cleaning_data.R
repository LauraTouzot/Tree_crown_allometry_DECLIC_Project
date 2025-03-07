## checking allometry measures at the tree level and returning warning signals if detected

checking_tree_data <- function(all_crown_clean) {
  
  df <- all_crown_clean %>% mutate(tree_issues = NA_real_)
  df$tree_issues[df$DBH_cm > 1000 | df$HT_m > 120 | df$C_diam_m > 30] <- "TRUE"
  df$tree_issues[is.na(df$tree_issues)] <- "FALSE"
  
  return (df)

}


checking_plot_data <- function(location_variables) {
  
  df <- location_variables %>% mutate(plot_issues = NA_real_)
  df$plot_issues[df$Ntree_ha > 3000] <- "TRUE"
  df$plot_issues[is.na(df$plot_issues)] <- "FALSE"
  
  p1 <- ggplot(df, aes(x = dataset, y = Ntree_ha)) + geom_boxplot()
  p2 <- ggplot(df, aes(x = dataset, y = ba_plot)) + geom_boxplot()
  png("figures/unit_check/plot_variables_check.png", width = 800, height = 500)
  multiplot(p1, p2, cols =  1, rows = 2)
  dev.off()
  
  return (df)
  
}


### Cleaning data based on the first filter before the taxonomy check

cleaning_after_checking <- function(all_crown_checked) {
  
  df <- all_crown_checked[all_crown_checked$location_ID != "TN_127_NA_1_1_FHM_1999",] # 3 extremely high measures of crown diameter (> 30m) for small dbh (around 25 cm) in this plot 
  df$HT_m[df$location_ID == "Middle_Ural,_Nizhnie_Sergi_Usoltev" & df$DBH_cm == 12.5 & df$sp == "Pinus sibirica"] <- NA
  
  return(df)
  
}
  

