# get SLA from TRY public to look at petiol effect

extract_traits_public_try_light <- function(df){
  
  spvec <- names(table(df$sp)[table(df$sp)>50])
  print(length(spvec))
  TRYdata1 <- data.table::fread(file.path("data", "5374.txt"),
                                header = T, sep = "\t", dec = ".",
                                quote = "", fill = TRUE, na.strings= "")
  TRYdata2 <- data.table::fread(file.path("data", "5385.txt"),
                                header = T, sep = "\t", dec = ".", quote = "",
                                fill = TRUE, na.strings= "")
  TRYdata <- rbind(TRYdata1, TRYdata2)
  rm(TRYdata1, TRYdata2)
  # TRAITS to select
  LeafArea <- c(1, 3108:3113)
  
  SLAnp <- c(3115) # no petiol , 3117
  SLAp <- c(3116) # with petiol
  SLA <- c(3116, 3115, 3117) # with or without petiol
  df_SLA <- bind_rows(lapply(spvec, extractAll, TRYdata, SLA))
  df_SLAp <- bind_rows(lapply(spvec, extractAll, TRYdata, SLAp))
  df_SLAnp <- bind_rows(lapply(spvec, extractAll, TRYdata, SLAnp))
  df_SLA$Unit <- as.vector(na.omit(unique(TRYdata$UnitName[TRYdata$TraitID %in% SLA &
                                                             !is.na(TRYdata$TraitID)])))
  df_SLAp$Unit <- as.vector(na.omit(unique(TRYdata$UnitName[TRYdata$TraitID %in% SLAp &
                                                              !is.na(TRYdata$TraitID)])))
  df_SLAnp$Unit <- as.vector(na.omit(unique(TRYdata$UnitName[TRYdata$TraitID %in% SLAnp &
                                                               !is.na(TRYdata$TraitID)])))
  names(df_SLA) <- paste0("SLA_", names(df_SLA))
  names(df_SLAp) <- paste0("SLAp_", names(df_SLAp))
  names(df_SLAnp) <- paste0("SLAnp_", names(df_SLAnp))
  
  df_LeafArea <- bind_rows(lapply(spvec, extractAll, TRYdata, LeafArea))
  df_LeafArea$Unit <- as.vector(na.omit(unique(TRYdata$UnitName[TRYdata$TraitID %in% LeafArea &
                                                                  !is.na(TRYdata$TraitID)])))
  names(df_LeafArea) <- paste0("LeafArea_", names(df_LeafArea))
  res <- bind_cols(df_SLA, df_SLAp, df_SLAnp, df_LeafArea)
  res <- data.frame(sp = spvec, res)
  return(res)
}


plot_SLA_petiol <- function(SLA, SLAf){
  
  SLA_mean  <- SLAf %>% group_by(Species) %>%
    summarise(SLA = mean(SLA_m2_kg, na.rm = TRUE))
  SLA_f_try <- left_join(SLA_mean,SLA, by = c("Species" = "sp"))
  
  require(gridExtra)
  
  p1 <- ggplot(SLA, aes(x = SLAp_mean, y = SLAnp_mean))+
    geom_point() + xlab("SLA TRY with petiol") +
    geom_abline(slope = 1, intercept = 0)+
    ylab("SLA TRY without petiol") + theme_bw()
  p2 <- ggplot(SLA, aes(x = SLAp_mean, y = SLA_mean))+
    geom_point() + xlab("SLA TRY with petiol") +
    geom_abline(slope = 1, intercept = 0)+
    ylab("SLA TRY all") + theme_bw()
  p3 <- ggplot(SLA_f_try, aes(x = SLAp_mean, y = SLA))+
    geom_point() + xlab("SLA TRY with petiol") +
    geom_abline(slope = 1, intercept = 0)+
    ylab("SLA Forrester") + theme_bw()
  p4 <- ggplot(SLA_f_try, aes(x = SLA_mean, y = SLA))+
    geom_point() + xlab("SLA TRY all") +
    geom_abline(slope = 1, intercept = 0)+
    ylab("SLA Forrester") + theme_bw()
  
  ggsave(("figures/SLA_petiol.pdf"), arrangeGrob(p1, p2, p3, p4, ncol=2, nrow = 2))  
}
