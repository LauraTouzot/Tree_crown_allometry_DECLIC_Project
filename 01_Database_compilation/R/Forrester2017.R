#############
###  data from Forrester paper
data_forrester  <- function(file = file.path("data", "Forrester", "1-s2.0-S0378112717301238-mmc1.xlsx")){
  
  data <- read_excel(file,
                     sheet = "Table A.2", skip = 11)
  spvec <- unique(data$Species)
  spname <- c("Abies alba", "Alnus glutinosa", "Alnus incana",
              "Acer pseudoplatanus", "Betula pendula Or B. pubescens",
              "Betula pendula", "Betula pendula Or B. pubescens",
              "Betula pubescens", "Carpinus betulus",
              "Castanea sativa", "Fraxinus excelsior",
              "Fagus sylvatica", "Larix decidua", "Picea abies",
              "Populus alba", "Prunus avium", "Pinus cembra",
              "Pseudotsuga menziesii", "Pinus nigra",
              "Populus hybrids", "Pinus pinaster",
              "Prunus serotina", "Pinus sylvestris",
              "Populus tremula", "Populus hybrids",
              "Quercus ilex", "Quercus petraea",
              "Quercus robur", "Quercus petraea Or Q. robur",
              "Robinia pseudoacacia", "Sorbus aucuparia",
              "Tilia cordata", "Tilia cordata Or T. platyphyllos")
  
  names(spname) <- spvec
  data$Species2 <- spname[data$Species]
  data$Genus  <- get_genus(data$Species2)
  data<- data[data$Component == "Foliage",]
  names(data)  <- gsub(".","",gsub(",","",gsub(")", "", gsub("(","", gsub(" ", "_", names(data)),
                                                             fixed = TRUE), fixed = TRUE), fixed = TRUE), fixed = TRUE)
  
  old_names  <- names(data)
  old_names[16] <- "DBH"
  old_names[20] <- "LeafArea"
  old_names[21] <- "LeafMass"
  old_names[22] <- "LeafMassMSE"
  names(data) <- old_names
  
  return(data)
}  


data_forrester_bin  <- function(data){
  
  cut_diam  <- seq(from = 0, to =170, by = 5)
  dbh_cut <- cut_diam[-length(cut_diam)] +mean(diff(cut_diam))/2
  data$dbh_c <- dbh_cut[cut(data$DBH, breaks = cut_diam, labels = FALSE)]
  data_m <-  data %>% group_by(Genus, dbh_c) %>%
  summarise(LeafMass_m = mean(LeafMass, na.rm = TRUE),
            LeafMassMSE_m = mean(LeafMassMSE, na.rm = TRUE),
            LeafArea_m = mean(LeafArea, na.rm = TRUE))

return(data_m)
  
}


param_forrester  <- function(file = file.path("data", "Forrester", "1-s2.0-S0378112717301238-mmc1.xlsx")){

  params <- read_excel(file,
                       sheet = "Table A.4", skip = 22)
  names(params)  <-  gsub(substr(names(params)[6], 4,  4), "B",
                          names(params),
                          fixed = TRUE)
  names(params)  <-  gsub("(", "_",
                          names(params),
                          fixed = TRUE)
  names(params)  <-  gsub(")", "",
                          names(params),
                          fixed = TRUE)
  names(params)  <-  gsub(" ", "",
                          names(params),
                          fixed = TRUE)
  params$Species <- gsub("Castania", "Castanea", params$Species)
  params$Species <- gsub("Pseudotuga", "Pseudotsuga", params$Species)
  params$Genus  <- get_genus(params$Species)
  
  params_sel  <- params[params$Component == "Foliage mass" & params$Equation == 3 &
                          params$ValidationorFinal == "Final model" &
                          is.na(params$Invalidequations), ]
  params_sel  <- tidyr::separate(params_sel, "ln_B0",
                                 into = c("lnB0", "lnB0_sd"), sep = "  \\(")
  params_sel  <- tidyr::separate(params_sel, "Bforln_d",
                                 into = c("Blnd", "Blnd_sd"), sep = "  \\(")
  params_sel$lnB0 <- as.numeric( params_sel$lnB0)
  params_sel$Blnd <- as.numeric( params_sel$Blnd)
  params_sel  <-  as.data.frame(params_sel)
  
 return(params_sel)
  
}

data_SLA_Forrester <- function(file = file.path("data", "Forrester", "1-s2.0-S0378112717301238-mmc1.xlsx")){
  # read SLA Forrester
  SLAf <- read_excel(file,
                   sheet = "Table A.1", skip = 3)
  names(SLAf) <- c("Species", "Age", "SLA_m2_kg", "Reference")
  SLAf$Species <- gsub("Castania", "Castanea", SLAf$Species)
  SLAf$Species <- gsub("Pseudotuga", "Pseudotsuga", SLAf$Species)

return(SLAf)
  
}

model_forrester_refit <- function(data, data_m, param_sel){
  # fit equation

  mat_res  <- matrix(NA, nrow = length(unique(data$Genus)), ncol = 3)
  colnames(mat_res) <- c("lnB0", "Blnd", "CF")
  row.names(mat_res)  <- unique(data$Genus)
  dfp <- data.frame(Genus = NA, dbh = NA,
                    LeafMass =NA)
  mat_resb  <- matrix(NA, nrow = length(unique(data$Genus)), ncol = 3)
  colnames(mat_resb) <- c("lnB0", "Blnd", "CF")
  row.names(mat_resb)  <- unique(data$Genus)
  dfpb <- data.frame(Genus = NA, dbh = NA,
                     LeafMass =NA)
  
  dbhs  <- 1:100
  for (gs in unique(data$Genus)){
    df  <- data[data$Genus == gs,]
    if(length(unique(df$Equation_Number)) >1){
      res <- lmer(log(LeafMassMSE)~log(DBH) + (1 | Equation_Number), df)
      mat_res[gs, 1:2] <- fixef(res)
      mat_res[gs, 3]  <- mean(df$LeafMassMSE, na.rm = TRUE)/mean(exp(predict(res,
                                                                             re.form = NA)), na.rm = TRUE)
    }else{
      res <- lm(log(LeafMassMSE)~log(DBH), df)
      mat_res[gs, 1:2] <- coefficients(res)
      mat_res[gs, 3]  <- mean(df$LeafMassMSE, na.rm = TRUE)/mean(exp(predict(res
      )), na.rm = TRUE)
      
    }
    res <- lm(log(LeafMassMSE)~log(DBH), df)
    mat_resb[gs, 1:2] <- coefficients(res)
    mat_resb[gs, 3]  <- mean(df$LeafMassMSE, na.rm = TRUE)/mean(exp(predict(res
    )), na.rm = TRUE)
    
    dfp <- rbind(dfp, data.frame(Genus = gs, dbh = dbhs,
                                 LeafMass = exp( mat_res[gs, 1] +
                                                   log(dbhs)*mat_res[gs, 2])*mat_res[gs, 3]))
    dfpb <- rbind(dfpb, data.frame(Genus = gs, dbh = dbhs,
                                   LeafMass = exp( mat_resb[gs, 1] +
                                                     log(dbhs)*mat_resb[gs, 2])*mat_resb[gs, 3]))
  }
  dfp <- dfp[-1,]
  dfpb <- dfpb[-1,]
  
  dfp2 <- data.frame(Genus = NA, dbh = NA,
                     LeafMass =NA)
  mat_res2  <- matrix(NA, nrow = length(unique(data_m$Genus)), ncol = 3)
  colnames(mat_res2) <- c("lnB0", "Blnd", "CF")
  row.names(mat_res2)  <- unique(data_m$Genus)
  dbhs  <- 1:100
  for (gs in unique(data_m$Genus)){
    df  <- data_m[data_m$Genus == gs,]
    res <- lm(log(LeafMassMSE_m)~log(dbh_c), df)
    mat_res2[gs, 1:2] <- coefficients(res)
    mat_res2[gs, 3]  <- mean(df$LeafMassMSE_m, na.rm = TRUE)/mean(exp(predict(res)), na.rm = TRUE)
    dfp2 <- rbind(dfp2, data.frame(Genus = gs, dbh = dbhs,
                                   LeafMass = exp( coefficients(res)[1] +
                                                     log(dbhs)*coefficients(res)[2])*mat_res2[gs, 3]))
  }
  dfp2 <- dfp2[-1,]
  dfp$type <- "random"
  dfpb$type <- "all"
  dfp2$type <- "bins"
  dfall <- rbind(dfp, dfp2)
  
  return(dfall)
}

# look at data from Forrester paper
plot_forrester_data_refit  <- function(dfall,data,data_m){

    ggplot(data, aes(x = DBH, y = LeafMassMSE))+geom_point()+
        geom_point(data = data_m, aes(x=dbh_c, y = LeafMassMSE_m), color = "blue")+
        geom_line(data = dfall, aes(x=dbh, y = LeafMass, color = type, linetype = type))+
        xlim(c(1,100))+ ylim(c(0, 300))+ xlab("Tree diameter (cm)") + ylab("Leaf Mass (kg)")+
        facet_wrap(vars(Genus), scales = "free") + theme_bw()
    ggsave(("figures/unit_check/LeafMassRefitForrester.pdf"), width = 20, units ="cm")
  
}

# look at data from Forrester paper
plot_forrester_data_SLA_TRY  <- function(SLA, SLAf, dfall, data){

  SLA_mean  <- SLAf %>% group_by(Species) %>%
    summarise(SLA = mean(SLA_m2_kg, na.rm = TRUE))
  ### GET SLA genus   
  SLA$sp <- gsub("Robinia pseudacacia", "Robinia pseudoacacia", SLA$sp)
  vecspsla  <-  c(unique(data$Species2), "Betula")
  SLA  <-  SLA[SLA$sp %in% vecspsla, ]   
  SLA$Genus <- get_genus(SLA$sp)
  SLAb  <- left_join(SLA_mean, SLA, by = c("Species" = "sp"))
  SLAgt <- SLA %>% group_by(Genus) %>% summarise(SLAtry = mean(SLAp_mean, na.rm = TRUE))
  SLA_mean$Genus  <-  get_genus(SLA_mean$Species)  
  SLAgf <- SLA_mean %>% group_by(Genus) %>% summarise(SLAf = mean(SLA, na.rm = TRUE))
  ### multiply by SLA
  
  SLAf  <-  SLAgf$SLAf
  SLAtry  <-  SLAgt$SLAtry
  names(SLAtry)  <- SLAgt$Genus 
  names(SLAf)  <- SLAgf$Genus 
  dfall1 <- dfall2 <- dfall
  dfall1$LeafArea  <- dfall$LeafMass * SLAf[dfall$Genus]
  dfall1$SLAtype <- "Forrester"
  dfall2$LeafArea  <- dfall$LeafMass * SLAtry[dfall$Genus]
  dfall2$SLAtype <- "TRY"
  dfall_s <- rbind(dfall1, dfall2)
  
  
 ggplot(dfall_s, aes(x = dbh, y = LeafArea, color = SLAtype, linetype = type))+
    geom_line() +ylim(c(0, 1000)) + xlab("Tree diameter (cm)") +
    ylab("Leaf area (m2)") + theme_bw()+ facet_wrap(vars(Genus), scales = "free")
  ggsave(("figures/unit_check/ForresterGeorgesRefitTrySLALeafPred.pdf"))
  
}



plot_LeafMass_Equ <- function(){
  ## read params Forrester
  params <- read_excel(file.path("data", "Forrester", "1-s2.0-S0378112717301238-mmc1.xlsx"),
                       sheet = "Table A.4", skip = 22)
  names(params)  <-  gsub(substr(names(params)[6], 4,  4), "B",
                          names(params),
                          fixed = TRUE)
  names(params)  <-  gsub("(", "_",
                          names(params),
                          fixed = TRUE)
  names(params)  <-  gsub(")", "",
                          names(params),
                          fixed = TRUE)
  names(params)  <-  gsub(" ", "",
                          names(params),
                          fixed = TRUE)
  params$Species <- gsub("Castania", "Castanea", params$Species)
  params$Species <- gsub("Pseudotuga", "Pseudotsuga", params$Species)
  params$Genus  <- get_genus(params$Species)
  
  params_sel  <- params[params$Component == "Foliage mass" & params$Equation == 3 &
                          params$ValidationorFinal == "Final model" &
                          is.na(params$Invalidequations), ]
  params_sel  <- tidyr::separate(params_sel, "ln_B0",
                                 into = c("lnB0", "lnB0_sd"), sep = "  \\(")
  params_sel  <- tidyr::separate(params_sel, "Bforln_d",
                                 into = c("Blnd", "Blnd_sd"), sep = "  \\(")
  params_sel$lnB0 <- as.numeric( params_sel$lnB0)
  params_sel$Blnd <- as.numeric( params_sel$Blnd)
  params_sel  <-  as.data.frame(params_sel)
  dbh_seq  <- 1:100
  cols <- pals::glasbey(nrow(params_sel))[2:nrow(params_sel)]
  
  
  pdf(("figures/unit_check/ForresterLeafMassEqu.pdf"), width = 10)
  i  <-  1
  plot(dbh_seq, exp(params_sel[i, "lnB0"] +
                      log(dbh_seq)* params_sel[i, "Blnd"]) *
         params_sel[i, "CF"], xlab = "DBH", ylab = "Foliage mass (kg)",
       type = "l", ylim = c(0,200), lwd = 2, col = cols[1])
  dfpF <- data.frame(Genus = params_sel[i, "Genus"],
                     Species = params_sel[i, "Species"],
                     dbh = dbh_seq,
                     LeafMass =exp(params_sel[i, "lnB0"] +
                                     log(dbh_seq)* params_sel[i, "Blnd"]) *
                       params_sel[i, "CF"])
  
  for (i in 2:nrow(params_sel)){
    lines(dbh_seq, exp(params_sel[i, "lnB0"] +
                         log(dbh_seq)* params_sel[i, "Blnd"]) *
            params_sel[i, "CF"], col = cols[i], lty = i,
          lwd = 2)
    dfpF <- rbind(dfpF, data.frame(Genus = params_sel[i, "Genus"],
                                   Species = params_sel[i, "Species"],
                                   dbh = dbh_seq,
                                   LeafMass =exp(params_sel[i, "lnB0"] +
                                                   log(dbh_seq)* params_sel[i, "Blnd"]) *
                                     params_sel[i, "CF"]))
    
  }
  legend(0, 200, legend = params_sel$Species,
         col = cols, lty = 1:nrow(params_sel),
         bty = "n", lwd = 2)
  dev.off()
  
  
}

plot_LeafArea_Equ <- function(){
  ## read params Forrester
  params <- read_excel(file.path("data", "Forrester", "1-s2.0-S0378112717301238-mmc1.xlsx"),
                       sheet = "Table A.4", skip = 22)
  names(params)  <-  gsub(substr(names(params)[6], 4,  4), "B",
                          names(params),
                          fixed = TRUE)
  names(params)  <-  gsub("(", "_",
                          names(params),
                          fixed = TRUE)
  names(params)  <-  gsub(")", "",
                          names(params),
                          fixed = TRUE)
  names(params)  <-  gsub(" ", "",
                          names(params),
                          fixed = TRUE)
  params$Species <- gsub("Castania", "Castanea", params$Species)
  params$Species <- gsub("Pseudotuga", "Pseudotsuga", params$Species)
  params$Genus  <- get_genus(params$Species)
  params_selA  <- params[params$Component == "Leaf area" & params$Equation == 3 &
                           params$ValidationorFinal == "Final model" &
                           is.na(params$Invalidequations), ]
  params_selA  <- tidyr::separate(params_selA, "ln_B0",
                                  into = c("lnB0", "lnB0_sd"), sep = "  \\(")
  params_selA  <- tidyr::separate(params_selA, "Bforln_d",
                                  into = c("Blnd", "Blnd_sd"), sep = "  \\(")
  params_selA$lnB0 <- as.numeric( params_selA$lnB0)
  params_selA$Blnd <- as.numeric( params_selA$Blnd)
  params_selA  <-  as.data.frame(params_selA)
  dbh_seq  <- 1:100
  cols <- pals::glasbey(nrow(params_selA))[2:nrow(params_selA)]
  
  pdf(("figures/unit_check/ForresterLeafAreaEqu.pdf"), width = 10)
  i  <-  1
  plot(dbh_seq, exp(params_selA[i, "lnB0"] +
                      log(dbh_seq)* params_selA[i, "Blnd"]) *
         params_selA[i, "CF"], xlab = "DBH", ylab = "Leaf area (m2)",
       type = "l", ylim = c(0,1500), lwd = 2, col = cols[1])
  
  dfpAF <- data.frame(Genus = params_selA[i, "Genus"],
                      Species = params_selA[i, "Species"],
                      dbh = dbh_seq,
                      LeafAreaOF =exp(params_selA[i, "lnB0"] +
                                        log(dbh_seq)* params_selA[i, "Blnd"]) *
                        params_selA[i, "CF"])
  
  for (i in 2:nrow(params_selA)){
    lines(dbh_seq, exp(params_selA[i, "lnB0"] +
                         log(dbh_seq)* params_selA[i, "Blnd"]) *
            params_selA[i, "CF"], col = cols[i], lty = i,
          lwd = 2)
    dfpAF <- rbind(dfpAF, data.frame(Genus = params_selA[i, "Genus"],
                                     Species = params_selA[i, "Species"],
                                     dbh = dbh_seq,
                                     LeafAreaOF =exp(params_selA[i, "lnB0"] +
                                                       log(dbh_seq)* params_selA[i, "Blnd"]) *
                                       params_selA[i, "CF"]))
    
  }
  legend(0, 1500, legend = params_selA$Species,
         col = cols, lty = 1:nrow(params_selA),
         bty = "n", lwd = 2)
  dev.off()
  
}

