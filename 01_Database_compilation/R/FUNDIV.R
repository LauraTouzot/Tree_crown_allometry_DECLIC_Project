read_FUNDIV_data <- function() {

  plots <- read.csv("data/FUNDIV_tables/FUNDIV_plots_original.csv", sep = " ")
  trees <- read.csv("data/FUNDIV_tables/FUNDIV_tree_original.csv", sep = " ")
  species <- read.csv("data/FUNDIV_tables/FunDiv_species_Nadja.csv", stringsAsFactors = FALSE, fileEncoding = "ISO-8859-15")
  harvs <- read.csv("data/FUNDIV_tables/FunDivEUROPE_plot_management.csv", stringsAsFactors = FALSE)
  
  # clean tree data and get dbh values in cm
  trees$HT_m <- trees$height1
  trees$DBH_cm <- trees$dbh1/10
  
  # combine plot and tree data to get all required variables
  data <- left_join(trees, plots, by = "plotcode")
  
  # remove Wallonia from data (dbh measured at a different height)
  data <- data[which(data$country.x != "WA"),]
  
  # remove characters at the end of the species name
  species$species[species$id ==277] <- "pubescens"
  species <- species %>% mutate(sp = paste(genus, species)) %>% dplyr::select(c(id, sp))
  species[species$id == 277, "sp"] <- "Quercus pubescens"
  species[species$id == 48, "sp"] <- "Betula"
  
  # add species within the data table
  data <- left_join(data, species, by = c("speciesid" = "id"))
  
  # remove plot with harvesting
  data <- group_by(data, plotcode) %>% mutate(N_harv = sum(treestatus_th == 3)) %>% filter(N_harv < 1)
  
  # remove plots with harvesting based on plot code
  plots_with_harv <- harvs$plotcode[harvs$management2 >0 & !is.na(harvs$management2)]
  data <- filter(data, !plotcode %in% plots_with_harv)
  
  # remove plots in the Canary Islands
  data <- filter(data, latitude > 30)
  
  # change FG to FR
  data <- data %>% mutate(country = replace(country.x, country.x =='FG', 'FR')) %>% as.data.frame()
  
  # clean data based on the minimum DBH (specific to each country)
  SW <- data %>% filter(country == "SW", DBH_cm > 10) 
  FI <- data %>% filter(country == "FI") # no minimum DBH in the Finish protocol
  FR <- data %>% filter(country == "FR", DBH_cm >= 7.5) 
  DE <- data %>% filter(country == "DE", DBH_cm > 7) 
  ES <- data %>% filter(country == "ES", DBH_cm >= 7.5)
  
  df <- bind_rows(SW, FI, FR, DE, ES)
  
  df$location_ID <- paste(df$plotcode, "FUNDIV_HT", sep = "_")
  df$W <- df$weight1
  
  df$latitude_plot <- df$latitude
  df$longitude_plot <- df$longitude
  df$latitude_tree <- NA
  df$longitude_tree <- NA
  df$continent <- "E_U"

  df <- df %>% ungroup()
  
  # plots to check unit within the dataset
  p1 <-  ggplot(df[!is.na(df$HT_m), ], aes(x = DBH_cm, y = HT_m)) + geom_point(alpha = 0.1) 
  png("figures/unit_check/unit_check_FUNDIV.png", width = 340, height = 340)
  p1
  dev.off()
  
  return(df)
  
}




############################### NOT USED IN THIS VERSION OF THE CODE ###############################

read_FUNDIV_plot_data <- function(file_plot, file_clim){
  plots <- read.csv(file_plot,
                    stringsAsFactors=FALSE)
  clims <- read.csv(file_clim,
                    stringsAsFactors=FALSE)
  clims <- clims %>% dplyr::select(-longitude, -latitude, -country )

  # merging data
  data <- left_join( plots, clims, by = "plotcode")
  return(data)
}

read_FUNDIV_tree_data <- function(file_trees, file_plot, file_clim, file_species, data_manag, remove_harv = TRUE) {
  #function to read and select data
  require(dplyr)
  trees <- read.csv(file_trees,
                    stringsAsFactors=FALSE)
  # combine Betula pendula and pubescens with the Betula genus (ids 46 and 47 with 48)
  print("read tree")
  trees$speciesid[trees$speciesid %in% c(46,47)] <- 48
  species <- read.csv(file_species,
                      stringsAsFactors=FALSE, fileEncoding = "ISO-8859-15")
  # remove characters at the end of the species name
  species$species[species$id ==277] <- "pubescens"
  species <- species %>% mutate(sp = paste(genus, species)) %>%
      dplyr::select(c(id, sp))
  species[species$id == 277, "sp"] <- "Quercus pubescens"
  species[species$id == 48, "sp"] <- "Betula"
  plots <- read.csv(file_plot,
                    stringsAsFactors=FALSE)
  # not sure if this is only me but I have encoding problem with the species csv file
  plots <- plots %>%
    dplyr::select(-country, -ba_ha1, -ba_ha2,-surveydate1, -surveydate2)
  clims <- read.csv(file_clim,
                    stringsAsFactors=FALSE)
    print("read climatic data")
  clims <- clims %>% dplyr::select(-longitude, -latitude, -country )
  # merging data
  data <- left_join(trees, species, by = c("speciesid" = "id"))
  data <- left_join(data, plots, by = "plotcode")
  data <- left_join(data, clims, by = "plotcode")

  data <-  data %>% group_by(plotcode) %>%
      mutate(BATOT_ha1 = sum(ba_ha1))

  ## funBASUP <- function(df) {
  ##     data.frame(treecode = df$treecode,
  ##                BASUP = sapply(df$dbh1,
  ##                               function(x, dd) sum(dd$ba_ha1[dd$dbh1>x]),
  ##                               dd = df))
  ## }
  ## dfBASUP <- data %>% do(funBASUP(.))
  ## data <- left_join(data, dfBASUP, by = "treecode")


  if (remove_harv){
  # remove plot with harvesting
  data <- group_by(data, plotcode) %>%
           mutate(N_harv = sum(treestatus_th == 3)) %>%
           filter(N_harv <1)
  
  # remove French plots where management has been recorded, we don't have details
  # on individually harvested trees, the management column contains a 1 for those French
  # plots in which management has been recorded
  data$management[is.na(data$management)] <- 0
  data <- filter(data, management ==0)
  
  # remove plots with harvesting based on plot code
  harvs <- read.csv(data_manag,
                    stringsAsFactors=FALSE)
  plots_with_harv <- harvs$plotcode[harvs$management2 >0 & !is.na(harvs$management2)]
  data <- filter(data, ! plotcode %in% plots_with_harv)
  }
  # remove plots in the Canary Islands
  data <- filter(data, latitude >30)
  # calculate the number per hectare from the weight
  # data$n_ha1 <-  NA
  # data$n_ha1[data$country %in% c('ES','FI','SW','WA') &
            # data$treestatus_th>1 &
            # !is.na(data$weight1)] <- 1/((data$weight1[data$country %in%
                                                    #   c('ES','FI','SW','WA') &
                                                    #   data$treestatus_th>1 &
                                                    #   !is.na(data$weight1)]*
                                        #  data$weight1[data$country %in%
                                                     #  c('ES','FI','SW','WA') &
                                                     #  data$treestatus_th>1 &
                                                     #  !is.na(data$weight1)]*
                                        #  3.14159265)/10000)
  # data$n_ha2 <-  NA
  # data$n_ha2[data$country %in% c('ES','FI','SW','WA') &
         #    data$treestatus_th<3 &
         #    !is.na(data$weight2)] <- 1/((data$weight2[data$country %in%
                                                  #     c('ES','FI','SW','WA') &
                                                  #     data$treestatus_th<3 &
                                                  #    !is.na(data$weight2)]*
                                     #     data$weight2[data$country %in%
                                                  #     c('ES','FI','SW','WA') &
                                                  #     data$treestatus_th<3 &
                                     #             #     !is.na(data$weight2)]*
                                     #     3.14159265)/10000)
 # data$n_ha1[data$country =='DE'] <- data$weight1[data$country =='DE']
 # data$n_ha2[data$country =='DE'] <- data$weight2[data$country =='DE']
  
 # weight is missing in France recompute it base on ba and ba_ha
  data$n_ha1[data$country =='FG' &
               data$treestatus_th>1]  <- data$ba_ha1[data$country =='FG' &
                                              data$treestatus_th>1]/data$ba1[data$country =='FG' &
                                                                                data$treestatus_th>1]
  data$n_ha2[data$country =='FG' &
               data$treestatus_th<3]  <- data$ba_ha2[data$country =='FG' &
                                                        data$treestatus_th<3]/data$ba2[data$country =='FG' &
                                                        data$treestatus_th<3]
  
  #compute competitors basal area
  # data <- data %>% mutate(BATOTcomp = BATOT_ha1 - ba_ha1)
  
  # CHANGE FG to FR
  data <- data %>% mutate(country = replace(country, country=='FG', 'FR')) %>%
     as.data.frame()

  data <- filter(data, !is.na(BATOT_ha1) & !is.na(surveydate2)& !is.na(surveydate1))
  data$treecode2 <- paste0("code_", seq_len(length.out = nrow(data)))
  ## remove unused variables
  if (remove_harv){
  data <- dplyr::select(data, -c(treecode, speciesid, ba1,
                          ba_ha1, ba2, bachange_ha_yr, dbh1_mod,
                          ba1_mod, bachange_ha_yr_mod, weight2,
                          biome, surveydate1, end_year, tile,
                          bio1, BATOT_ha1, N_harv, n_ha1))
  }else{
  data <- dplyr::select(data, -c(treecode, speciesid, ba1,
                          ba_ha1, ba2, bachange_ha_yr, dbh1_mod,
                          ba1_mod, bachange_ha_yr_mod, weight2,
                          biome, surveydate1, end_year, tile,
                          bio1, BATOT_ha1, n_ha1))
  }
    return(data)
}


LatLongSp <- function(spsel, df){
df %>% dplyr::filter(sp == spsel) %>%  select(plotcode, longitude, latitude)
}

FUNDIV_data_for_sdm <- function(df,sps){
ll <- lapply(sps$sp, LatLongSp, df = df)
names(ll) <- sps$sp
saveRDS(ll, file.path("output", "sdm_data.rds"))
}

select_species <- function(data, nlim_tree  = 2000, nlim_plot = 500){
  # select the species with sufficient number of individuals
  data <- filter(data, sp != " " & treestatus_th != 1)
  nplot <- group_by(data, sp) %>% summarise(n_plot= n_distinct(plotcode))
  abund <- group_by(data, sp) %>% summarise(n_tree= n())
  res <- filter(arrange(left_join(nplot, abund, "sp"), desc(n_tree)),  n_tree> nlim_tree & n_plot > nlim_plot)
  #exclude exotic and castanea
  print("species selected")
  res <- filter(res, !sp %in% c("Pinus radiata", "Eucalyptus globulus",
                                "Castanea sativa",
                                "Pseudotsuga menziesii",
                                "Robinia pseudacacia",
                                "Quercus "))
  return(res)
}


get_height_FUNDIV_sp <- function(df){
  library(dplyr)
  res <- df %>% group_by(sp) %>% summarise(N_HT = sum(!is.na(height1)))
  return(res)  
}

#####################################################################################################

