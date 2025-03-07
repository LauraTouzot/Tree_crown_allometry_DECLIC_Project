# Read data Dettmann 2018


data_dettmann <- function(file = "data/Dettmann_MacFarlane_2018/Dettmann_MacFarlane_2018.csv"){

  df <- read.csv(file)
  df <- df %>% rename(sp = SPECIES,
                      DBH_cm = DBH_CM,
                      CR = LIVE_CROWN_RATIO)
  df <- df %>% mutate(C_diam_m = sqrt(CROWN_AREA_M/pi)*2)
  
  # cleaning for minimum DBH
  df <- df[df$DBH_cm > 2,]
  
  df$latitude_plot <- NA
  df$longitude_plot <- NA
  df$latitude_tree <- NA
  df$longitude_tree <- NA
  df$location_ID <- NA
  df$continent <- "N_A"
  df$W <- NA
  
  # plots to check unit within the dataset
  p1 <-  ggplot(df[!is.na(df$C_diam_m), ], aes(x = DBH_cm, y = C_diam_m)) + geom_point(alpha = 0.1) 
  p2 <-  ggplot(df[!is.na(df$CR), ], aes(x = DBH_cm, y = CR)) + geom_point(alpha = 0.1) 
  png("figures/unit_check/unit_check_Dettmann.png", width = 680, height = 340)
  multiplot(p1, p2, cols =  2)
  dev.off()
  
  return(df)
}


############################### NOT USED IN THIS VERSION OF THE CODE ###############################

# plot_WD_Dettmann <- function(df){
#   res <- lm(log(LEAF_DRY_MASS_KG)~log(DBH_cm) + log(SPECIFIC_GRAVITY_WOOD) + log(COMPETITION_INDEX+0.0001) +log(CR) + log(BASAL_AREA_INCREMENT_5YEAR_MEAN), data = df)
#   library(effects)
#   par(mfrow = c(1,2))
#   require(gridExtra)
#   p1 <- ggplot(df, aes(x = log(SPECIFIC_GRAVITY_WOOD), y = log(df$LEAF_DRY_MASS_KG/df$DBH_cm))) + geom_point() + theme_bw()
#   p2 <- plot(effect("log(SPECIFIC_GRAVITY_WOOD)",res), main = NA)
#   ggsave("figures/Leaf_Mass_WD_Dettmann.pdf", grid.arrange(p1, p2)) 
# }  

#####################################################################################################

