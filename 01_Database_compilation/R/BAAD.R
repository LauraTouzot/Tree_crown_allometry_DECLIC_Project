## READ BAAD data
## From code of Duursma and Falster 2016 New Phytologist https://github.com/remkoduursma/baadanalysis

extract_baad_data <- function(baad) {
  baad$data
}

extract_baad_dictionary <- function(baad) {
  baad$dictionary
}


# Convert conifer leaf area to projected leaf area
convertConiferLA <- function(baad) {

  # One sided total leaf area (half total surface area)
  # to projected area. Average of species in Barclay & Goodman 2000
  # For Pinus we use the value for lodgepole pine (see below)
  lambda1 <- c(0.873, 0.92, 0.879, 0.864, 0.839) # From Table 3 in Barclay & Goodman 2000
  ola_pla <- mean(1/lambda1)

  conv_pine <- function(x, method){
    if(method %in% c("","?","ax"))
      method <- "unknown"

    cv <- 1/0.778
    switch(method,
           a4 = x * cv,
           a5 = x/2,
           a6 = x,
           a7 = x/2,
           unknown = x
    )
  }

  conv_nonpine <- function(x, method){
    if(method %in% c("","?","ax"))
      method <- "unknown"

    switch(method,
           a4 = x * ola_pla,
           a5 = x/2,
           a6 = x,
           a7 = x/2,
           unknown = x
    )
  }

  convf <- function(x, method, species, pft){

    if(pft %in% c("DA","EA"))return(x)

    if(grepl("Pinus", species,ignore.case=TRUE))
      conv_pine(x, method)
    else
      conv_nonpine(x, method)
  }
  conv <- Vectorize(convf)

  alfmeth <- baad$methods[,c("studyName","a.lf")]
  alfmeth$method_alf <- str_extract(alfmeth$a.lf,"a[4-7]{1}")
  alfmeth$method_alf[is.na(alfmeth$method_alf)] <- ""
  baad$data <- merge(baad$data, alfmeth[,c("studyName","method_alf")],all=T)

  with(baad$data, conv(a.lf, method_alf, speciesMatched, pft))

  baad
}


BasalA_fit <- function(baad){

  # Predicted basal diameter. See R/predict_dba...R
  test <- subset(baad$data, !is.na(a.stba) & !is.na(a.stbh) & !is.na(h.t) & !is.na(h.bh) & h.t > h.bh)
  nls(d.ba ~ d.bh * h.t^(c0*h.t^c1) /(h.t - h.bh)^(c0*h.t^c1), start = list(c0=0.9, c1 = 0.7),
      data=test)
  # Fit on diameter not cross-sectional area because a much better fit is obtained this way (not sure why!),
  # and more homogeneous error variance (because large cross-sectional area outliers)
}

predictBasalA <- function(dat, baad){

  fit <- BasalA_fit(baad)

  d.ba2 <- predict(fit, dat)
  d.ba2[dat$h.bh >= dat$h.t] <- NA
  d.ba2[!is.na(dat$d.ba)] <- dat$d.ba[!is.na(dat$d.ba)]

  (pi/4)*d.ba2^2
}

# Modified function

DBH_fit <- function(baad){
  dataset <- baad$data
  # Predicted  diameter at breast height from basal diameter
  df <- dataset[!is.na(dataset$d.ba) & !is.na(dataset$d.bh)  & !is.na(dataset$h.t) & dataset$h.t > 1.3, ]
  nls(d.bh ~ a * d.ba^b, data = df, start = list(a = 0.7, b = 0.9))
}

predictDBH <- function(dat, baad){

  fit <- DBH_fit(baad)

  d.bh2 <- predict(fit, dat)
  d.bh2[is.na(dat$d.bh) ] <- NA
  d.bh2[!is.na(dat$d.bh)] <- dat$d.bh[!is.na(dat$d.bh)]
  d.bh2
}


prepare_dataset_baad <- function(baad){

  # Prepare dataset for analysis.
  # - remove non-field grown plants
  ### growingCondition: FW (field wild) / FE (field experimental) / GH (glass house)
  ### PU (plantation unmanaged) / PM (plantation managed) / GC (growth chamber) / CG (common garden)
  # - add 'Group', interaction of species and studyName (i.e. species in different datasets
  # are assumed to be independent, not entirely a fair assumption but will account for large
  # environmental/management/measurement methods differences.)

  baad <- convertConiferLA(baad)

  # use only field studies with no thinning, and get rid of deciduous gymnosperm
  dataset <- droplevels(subset(baad$data, growingCondition %in% c("FW","PM","PU")))

  # add Group
  dataset <- within(dataset, {
    Group <- paste(studyName, speciesMatched)
  })

  # add unique location ID
  dataset <- within(dataset, {
    location_ID <- paste(location, "baad", sep = "_")
  })


  # Predict basal diameter from breast-height
  dataset$a.stba2 <- predictBasalA(dataset, baad)
  dataset$d.bh2 <- predictDBH(dataset, baad)

  # no minimum DBH in BAAD data

  dataset
}


data_formatted_baad <- function(dataset){
  data_crown <- dataset[!is.na(dataset$d.bh),] # & (!is.na(dataset$d.cr) | !is.na(dataset$c.d) | !is.na(dataset$h.t)), ]
  data_crown$DBH_cm <-  data_crown$d.bh*100
  data_crown$C_diam_m <-  data_crown$d.cr
  data_crown$CR <- data_crown$c.d/data_crown$h.t
  data_crown$C_depth_m <- data_crown$c.d
  data_crown$HT_m <- data_crown$h.t
  data_crown <- data_crown %>% dplyr::rename(sp = speciesMatched)
  data_crown$latitude_plot <- data_crown$latitude
  data_crown$longitude_plot <- data_crown$longitude
  data_crown$latitude_tree <- NA
  data_crown$longitude_tree <- NA
  data_crown$W <- NA

  return(data_crown)
}


# extract continents for all locations in the baad database

extract_continents_baad = function(baad_crown) {

  states_map <- readOGR(dsn="data/ne_10m_admin_1_states_provinces/",
                        layer="ne_10m_admin_1_states_provinces")

  countriesSP <- getMap(resolution = 'low')
  ID <- as.data.frame(1:dim(baad_crown)[1])
  colnames(ID) <- "ID_unique"
  baad_crown <- cbind(ID, baad_crown)

  no_coordinates <- baad_crown[is.na(baad_crown$latitude),]
  no_coordinates$continent <- "A_S"
  no_coordinates$cont <- NA

  baad_crown <- baad_crown[!is.na(baad_crown$latitude),]

  df <- baad_crown[,c("ID_unique", "latitude", "longitude")]
  df <- df[!is.na(df$latitude),]

  coordinates <- df[,c("longitude", "latitude")]

  pointsSP = SpatialPoints(coordinates, proj4string = CRS(" +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
  indices = over(pointsSP, countriesSP)

  locations <- as.data.frame(indices$REGION)
  colnames(locations) <- "cont"
  locations <- cbind(df$ID, locations)
  colnames(locations) <- c("ID_unique", "cont")

  continent <- as.data.frame(c("E_U", "N_A", "A_S", "S_A", "A_F", "A_U"))
  cont <- as.data.frame(c("Europe", "North America", "Asia", "South America", "Africa", "Australia"))

  continents <- cbind(continent, cont)
  colnames(continents) <- c("continent", "cont")

  locations <- left_join(locations, continents, by = "cont")
  baad_crown <- left_join(locations, baad_crown, by = "ID_unique")

  missing_continents <- baad_crown[is.na(baad_crown$continent),]
  loc <- as.data.frame(unique(missing_continents$location))
  adding_continents <- as.data.frame(c("A_S", "A_U", "A_S", "E_U", "A_S"))
  missing_data <- cbind(loc, adding_continents)
  colnames(missing_data) <- c("location", "continent")
  missing_continents <- missing_continents[, -3]
  missing_data <- left_join(missing_continents, missing_data, by = "location")

  baad_crown <- baad_crown[!is.na(baad_crown$continent),]
  baad_crown_complete <- rbind(missing_data, baad_crown, no_coordinates)

  baad_crown_complete <- baad_crown_complete %>% ungroup()

  # plots to check unit within the dataset
  p1 <-  ggplot(baad_crown_complete[!is.na(baad_crown_complete$HT_m), ], aes(x = DBH_cm, y = HT_m)) + geom_point(alpha = 0.1)
  p2 <-  ggplot(baad_crown_complete[!is.na(baad_crown_complete$CR), ], aes(x = DBH_cm, y = CR)) + geom_point(alpha = 0.1)
  p3 <-  ggplot(baad_crown_complete[!is.na(baad_crown_complete$C_depth_m), ], aes(x = DBH_cm, y = C_depth_m)) + geom_point(alpha = 0.1)
  p4 <-  ggplot(baad_crown_complete[!is.na(baad_crown_complete$C_diam_m), ], aes(x = DBH_cm, y = C_diam_m)) + geom_point(alpha = 0.1)
  png("figures/unit_check/unit_check_baad.png", width = 680, height = 680)
  multiplot(p1, p2, p3, p4, cols =  2)
  dev.off()

  return(baad_crown_complete)

}



############################### NOT USED IN THIS VERSION OF THE CODE ###############################

# d.cr crown width

# c.d crown depth

# ggplot(data = dataset, aes(d.bh, d.cr, color = pft)) +geom_point(alpha = 0.3)+
#   scale_x_continuous(trans='log10') + scale_y_continuous(trans='log2')
#
# ggplot(data = dataset, aes(d.bh, d.cr, color =  growingCondition)) +geom_point()+
#   scale_x_continuous(trans='log10') + scale_y_continuous(trans='log2') + facet_grid(cols = vars(pft)) +
#   geom_vline(xintercept =  0.01)
#
# ggplot(data = dataset, aes(d.bh, d.cr, color =  pft)) +geom_point()+
#   scale_x_continuous(trans='log10') + scale_y_continuous(trans='log2') + facet_grid(cols = vars(pft)) +
#   geom_vline(xintercept =  0.05)
#
# ggplot(data = dataset, aes(d.bh, d.cr, color = pft)) +geom_point(alpha = 0.3)+
#   scale_x_continuous(trans='log10') + scale_y_continuous(trans='log2') + facet_grid(cols = vars(growingCondition), rows = vars(pft))
#
# ggplot(data = dataset, aes(h.t, d.cr, color = pft)) +geom_point(alpha = 0.3)+
#   scale_x_continuous(trans='log10') + scale_y_continuous(trans='log2')
#
# ggplot(data = dataset, aes(h.t, c.d, color = pft)) +geom_point(alpha = 0.3)+
#   scale_x_continuous(trans='log10') + scale_y_continuous(trans='log2')

#####################################################################################################


