# Test Laura

library(dplyr)  
## Attention ne pas utiliser plyr et dplyr car il y a des fonctions qui ont le même nom!!!
library(tidyr)
library(ggplot2)
library(MuMIn)
samsara <- read.csv(file = "data/out_light_compet.csv", sep = ";")
no_comp <- read.csv(file = "data/out_light_alone.csv", sep = ";")
samsara <- samsara %>%  tidyr::separate_wider_delim(id_pattern, "_",
                                                                names = c("dbh", "comp", "rep"))

no_comp <- no_comp %>% dplyr::filter(dbh_cm == 15, 
                                     bal_m2ha == 0)
res <- left_join(samsara, no_comp, 
                 by = c("id_allom", "species", "order"),
                 suffix = c("_comp", "_nocomp")) %>% 
  mutate(lci = e_comp/e_nocomp,
         rep_cd_r = cradius_m_comp/cradius_m_nocomp,
         rep_h_r = h_m_comp/h_m_nocomp,
         rep_cr_r = cratio_comp/cratio_nocomp,
         rep_h = h_m_comp - h_m_nocomp,
         rep_cd = cradius_m_comp - cradius_m_nocomp,
         rep_cr = cratio_comp - cratio_nocomp)%>%
  filter((dbh == "15"& comp == "low") | (dbh == "30"& comp == "high")) %>%
  mutate(comp2 = paste(dbh, comp, sep = "_"))
  

# a few species ti id_allom give completely crazy lci because of a crown diameter response
boxplot(lci~species, res, las = 2)
abline(a = 1, b = 0, col = "red")
table(res$lci>1.3)/length(res$lci)
# il y a 1.8 % des données avec un lci > 1.3 
boxplot(rep_cd_r~species, res, las = 2)
ggplot(res, aes(h_m_nocomp, lci, col= species))+geom_point() +
  geom_hline(yintercept = 1)+facet_grid(order~comp2, scale = "free")
ggplot(res, aes(rep_h, lci, col= species))+geom_point() +
  geom_hline(yintercept = 1)+facet_grid(order~comp2, scale = "free")
ggplot(res, aes(rep_cd, lci, col= species))+geom_point() +
  geom_hline(yintercept = 1)+facet_grid(order~comp2, scale = "free")
ggplot(res, aes(rep_cr, lci, col= species))+geom_point() +
  geom_hline(yintercept = 1)+facet_grid(order~comp2, scale = "free")
 
pb <- as.data.frame(res) %>% group_by(species, id_allom) %>% 
  dplyr::summarise(n_pb_lci = sum(lci>10), 
                   n_pb_cd = sum(rep_cd >3))

pb_sp <- as.data.frame(res) %>% group_by(species) %>% 
  dplyr::summarise(n_pb_lci = sum(lci>10), 
                   n_pb_cd = sum(rep_cd >3))

# C'est principalement lié à qque espèces mais surtout qque fit
# I propose to remove this fit arguing that an increase of more than a factor 3 of crown diameter is unrealistic 
table(res$rep_cd>3)/length(res$lci) # less than 1% of the data
res <- res %>%
  filter(rep_cd_r < 3) 

df <- res %>% dplyr::filter(dbh == 30 & comp == "high") %>% 
      dplyr::select(rep,species,order,rep_cd_r,rep_h_r,rep_cr_r) %>% 
  pivot_longer(cols = c(rep_cd_r,rep_h_r,rep_cr_r), names_to = "crown_dim", values_to = "comp_rep")

ggplot(df, aes(x=crown_dim, y=comp_rep)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(limits = quantile(df$comp_rep, c(0.1, 0.9)))+
  facet_wrap(.~ order)
names(res)
pb <- as.data.frame(res) %>% group_by(species, id_allom) %>% 
  dplyr::summarise(n_pb_lci = sum(lci>2), 
                   n_pb_cd = sum(rep_cd >2))
ggplot(res, aes(h_m_nocomp, lci, col= species))+geom_point() +
  geom_hline(yintercept = 1)+facet_grid(order~comp2, scale = "free")
ggplot(res, aes(rep_h, lci, col= species))+geom_point() +
  geom_hline(yintercept = 1)+facet_grid(order~comp2, scale = "free")
ggplot(res, aes(rep_cd, lci, col= species))+geom_point() +
  geom_hline(yintercept = 1)+facet_grid(order~comp2, scale = "free")
ggplot(res, aes(rep_cr, lci, col= species))+geom_point() +
  geom_hline(yintercept = 1)+facet_grid(order~comp2, scale = "free")


# group_by group_by(species, order, comp2) 

res_sp <- res %>% group_by(species, order, comp2) %>% dplyr::summarise(cv_lci = sd(lci)/mean(lci),
                                                         lci = mean(lci),
                                                  h_nocomp = mean(h_m_nocomp),
                                                  rep_height = mean(rep_h),
                                                  rep_diam = mean(rep_cd),
                                                  rep_ratio = mean(rep_cr)) %>%
  ungroup() %>%group_by(order) %>% 
  mutate(across(c(lci, h_nocomp,rep_height, rep_diam, rep_ratio), scale))


res_sp <- res %>% group_by(order) %>% mutate(across(c(lci, h_m_nocomp,rep_h, rep_cd, rep_cr), scale))%>% 
  group_by(species, order, comp2) %>% dplyr::summarise(cv_lci = sd(lci)/mean(lci),
                                                                       lci = mean(lci),
                                                                       h_nocomp = mean(h_m_nocomp),
                                                                       rep_height = mean(rep_h),
                                                                       rep_diam = mean(rep_cd),
                                                                       rep_ratio = mean(rep_cr)) %>%
  ungroup() %>%group_by(order) 

ggplot(res_sp, aes(h_nocomp, lci, col= species, size = 1/cv_lci))+geom_point() +
  geom_hline(yintercept = 1)+facet_grid(order~comp2, scale = "free")
ggplot(res_sp, aes(rep_height, lci, col= species, size = 1/cv_lci))+geom_point() +
  geom_hline(yintercept = 1)+facet_grid(order~comp2, scale = "free")
ggplot(res_sp, aes(rep_diam, lci, col= species, size = 1/cv_lci))+geom_point() +
  geom_hline(yintercept = 1)+facet_grid(order~comp2, scale = "free")
ggplot(res_sp, aes(rep_ratio, lci, col= species, size = 1/cv_lci))+geom_point() +
  geom_hline(yintercept = 1)+facet_grid(order~comp2, scale = "free")


# Analysis codominance
options(na.action = "na.fail")
# Angio
df_a_codo <- na.omit(res_sp[res_sp$order =="A" & res_sp$comp2 == "15_low",])
model_angio_codo <- lm(lci ~ h_nocomp + rep_height + rep_ratio + rep_diam, weights = (1/cv_lci), 
                  data = df_a_codo , na.action = na.fail) 
angio_all_models_codo <- dredge(model_angio_codo) 
best_angio_model_codo <- get.models(angio_all_models_codo, 1)[[1]]
summary_angio_best_model_codo <- as.data.frame(angio_all_models_codo)

## Storing results
mod_angio_codo <- as.data.frame(confint(best_angio_model_codo))
mod_angio_codo$mean <- summary(best_angio_model_codo)$coefficients[, "Estimate"]
mod_angio_codo$p.value <- summary(best_angio_model_codo)$coefficients[, "Pr(>|t|)"]
mod_angio_codo <- mod_angio_codo[-1,] # remove the Intercept line of results
colnames(mod_angio_codo) <- c("est.inf", "est.sup", "est", "p.value")
mod_angio_codo <- mod_angio_codo %>% dplyr::mutate(mod = NA)

# Gymno
df_g_codo <- na.omit(res_sp[res_sp$order =="G" & res_sp$comp2 == "15_low",])
model_gymno_codo <- lm(lci ~ h_nocomp + rep_height + rep_ratio + rep_diam, weights = (1/cv_lci), 
                  data = df_g_codo , na.action = na.fail) 
gymno_all_models_codo <- dredge(model_gymno_codo) 
best_gymno_model_codo <- get.models(gymno_all_models_codo, 1)[[1]]
summary_gymno_best_model_codo <- as.data.frame(gymno_all_models_codo)

## Storing results
mod_gymno_codo <- as.data.frame(confint(best_gymno_model_codo))
mod_gymno_codo$mean <- summary(best_gymno_model_codo)$coefficients[, "Estimate"]
mod_gymno_codo$p.value <- summary(best_gymno_model_codo)$coefficients[, "Pr(>|t|)"]
mod_gymno_codo <- mod_gymno_codo[-1,] # remove the Intercept line of results
colnames(mod_gymno_codo) <- c("est.inf", "est.sup", "est", "p.value")
mod_gymno_codo <- mod_gymno_codo %>% dplyr::mutate(mod = NA)

# Analysis suppresed

# Angio
df_a_supr <- na.omit(res_sp[res_sp$order =="A" & res_sp$comp2 == "30_high",])
model_angio_supr <- lm(lci ~ h_nocomp + rep_height + rep_ratio + rep_diam, weights = (1/cv_lci), 
                  data = df_a_supr , na.action = na.fail) 
angio_all_models_supr <- dredge(model_angio_supr) 
best_angio_model_supr <- get.models(angio_all_models_supr, 1)[[1]]
summary_angio_best_model_supr <- as.data.frame(angio_all_models_supr)

## Storing results
mod_angio_supr <- as.data.frame(confint(best_angio_model_supr))
mod_angio_supr$mean <- summary(best_angio_model_supr)$coefficients[, "Estimate"]
mod_angio_supr$p.value <- summary(best_angio_model_supr)$coefficients[, "Pr(>|t|)"]
mod_angio_supr <- mod_angio_supr[-1,] # remove the Intercept line of results
colnames(mod_angio_supr) <- c("est.inf", "est.sup", "est", "p.value")
mod_angio_supr <- mod_angio_supr %>% dplyr::mutate(mod = NA)


# Gymno
df_g_supr <- na.omit(res_sp[res_sp$order =="G" & res_sp$comp2 == "30_high",])
model_gymno_supr <- lm(lci ~ h_nocomp + rep_height + rep_ratio + rep_diam, weights = (1/cv_lci), 
                  data = df_g_supr , na.action = na.fail) 
gymno_all_models_supr <- dredge(model_gymno_supr) 
best_gymno_model_supr <- get.models(gymno_all_models_supr, 1)[[1]]
summary_gymno_best_model_supr <- as.data.frame(gymno_all_models_supr)

## Storing results
mod_gymno_supr <- as.data.frame(confint(best_gymno_model_supr))
mod_gymno_supr$mean <- summary(best_gymno_model_supr)$coefficients[, "Estimate"]
mod_gymno_supr$p.value <- summary(best_gymno_model_supr)$coefficients[, "Pr(>|t|)"]
mod_gymno_supr <- mod_gymno_supr[-1,] # remove the Intercept line of results
colnames(mod_gymno_supr) <- c("est.inf", "est.sup", "est", "p.value")
mod_gymno_supr <- mod_gymno_supr %>% dplyr::mutate(mod = NA)


