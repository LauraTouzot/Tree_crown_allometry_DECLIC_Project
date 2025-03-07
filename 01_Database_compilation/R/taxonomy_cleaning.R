# ### Clean allometry database based on TNRS (using cleaned files from Georges as kind of tricky to get right version of curl on Laura's computer)
# 
# check_for_taxonomy_allometry <- function(all_crown_checked_bis) {
#   
#   clean_names <- readRDS(file = "data/TNRS_check.rds")
#   df <- left_join(all_crown_checked_bis, clean_names, by = "sp")
#   
#   return(df)
#   
# }


## Clean allometry database using TNRS

check_for_taxonomy_allometry <- function(all_crown_checked_bis) {

  # cleaning awkward species name before running the next step (au cas par cas)
  all_crown_checked_bis$sp <- str_replace_all(all_crown_checked_bis$sp, c("Mezcla de peque\xf1as frondosas" = "Unknown tree",
                                                  "Otros \xe1rboles rip\xedcolas" = "Unknown tree", "_" = " ",
                                                  "--- ---" = "Unknown tree", "Unknown" = "Unknown tree",
                                                  "Tree unknown" = "Unknown tree", "other" = "Unknown tree",
                                                  "indetermined species" = "Unknown tree", "Other broadleaves" = "Unknown tree",
                                                  "Other conifers" = "coniferous", "Family Arecaceae not listed above" = "Arecaceae",
                                                  "inconnu" = "Unknown tree"))

  species_names <- as.data.frame(unique(all_crown_checked_bis$sp))
  row_number <- as.data.frame(c(1:dim(species_names)[1]))

  names_to_check <- cbind(row_number, species_names)
  colnames(names_to_check) <- c("row_number", "species")

  # using 3 different sources to compare the obtained results and increase our chances to find a match for each species name
  results_1 <- TNRS(taxonomic_names = names_to_check, sources = "wfo", mode = "resolve") # World Flora
  results_1 <- results_1 %>% dplyr::select(Name_submitted, Accepted_name)

  results_2 <- TNRS(taxonomic_names = names_to_check, sources = "usda", mode = "resolve") # American database
  results_2 <- results_2 %>% dplyr::select(Name_submitted, Accepted_name)

  results_3 <- TNRS(taxonomic_names = names_to_check, sources = "tropicos", mode = "resolve") # Tropicos
  results_3 <- results_3 %>% dplyr::select(Name_submitted, Accepted_name)

  results <- cbind(results_1, results_2$Accepted_name, results_3$Accepted_name)
  colnames(results) <- c("submitted_name", "match_wfo", "match_usda", "match_tropicos")

  a <- coalesce(results$match_wfo, results$match_usda)
  checked_name <- coalesce(a, results$match_tropicos)

  results <- cbind(results, checked_name)
  results <- results %>% dplyr::select(submitted_name, checked_name)

  allometry_names_check <- left_join(all_crown_checked_bis, results, by = c("sp" = "submitted_name"))

  allometry_names_check$checked_name <- str_replace_all(allometry_names_check$checked_name, c("Nyssa sylvatica var. biflora" = "Nyssa sylvatica",
                                                                                              "Tilia americana var. heterophylla" = "Tilia americana",
                                                                                              "Populus x canadensis" = "Populus canadensis",
                                                                                              "Abies magnifica var. shastensis" = "Abies magnifica",
                                                                                              "Abies lasiocarpa var. arizonica" = "Abies lasiocarpa",
                                                                                              "Pinus flexilis var. reflexa" = "Pinus flexilis",
                                                                                              "Populus deltoides subsp. monilifera" = "Populus deltoides",
                                                                                              "Juniperus virginiana var. silicicola" = "Juniperus virginiana",
                                                                                              "Cercocarpus montanus var. paucidentatus" = "Cercocarpus montanus",
                                                                                              "Pinus cembroides var. bicolor" = "Pinus cembroides",
                                                                                              "Taxodium distichum var. imbricatum" = "Taxodium distichum",
                                                                                              "Carya ovata var. australis" = "Carya ovata",
                                                                                              "Abies x borisii-regis" =  "Abies borisii-regis",
                                                                                              "Fagus x taurica" = "Fagus taurica",
                                                                                              "Pinus nigra subsp. laricio" = "Pinus nigra",
                                                                                              "Pinus nigra subsp. salzmannii" = "Pinus nigra",
                                                                                              "Tilia x europaea" = "Tilia x europaea",
                                                                                              "Carpinus caroliniana subsp. virginiana" = "Carpinus caroliniana",
                                                                                              "Pinus contorta var. latifolia" = "Pinus contorta",
                                                                                              "Pseudotsuga menziesii var. glauca" = "Pseudotsuga menziesii",
                                                                                              "Ulmus davidiana var. japonica" = "Ulmus davidiana",
                                                                                              "Larix gmelinii var. olgensis" = "Larix gmelinii",
                                                                                              "Alnus incana subsp. rugosa" = "Alnus incana",
                                                                                              "Acer saccharum subsp. leucoderme" = "Acer saccharum",
                                                                                              "Tilia x europaea" = "Tilia europaea",
                                                                                              "Ledum palustre subsp. groenlandicum" = "Ledum palustre",
                                                                                              "Picea x lutzii" = "Picea lutzii",
                                                                                              "Amelanchier sanguinea var. gaspensis" = "Amelanchier sanguinea",
                                                                                              "Alnus alnobetula subsp. crispa" = "Alnus alnobetula",
                                                                                              "Populus x hybrida" = "Populus hybrida",
                                                                                              "Carya glabra var. odorata" = "Carya glabra"))

  return(allometry_names_check)

}



check_for_taxonomy_NFI <- function(species_all_NFI) {

  df <- species_all_NFI %>% dplyr::select(submitted_name, nplot, ntree, continent)

  species_names <- as.data.frame(unique(df$submitted_name))
  row_number <- as.data.frame(c(1:dim(species_names)[1]))

  names_to_check <- cbind(row_number, species_names)
  colnames(names_to_check) <- c("row_number", "species")

  # using 3 different sources to compare the obtained results and increase our chances to find a match for each species name
  results_1 <- TNRS(taxonomic_names = names_to_check, sources = "wfo", mode = "resolve") # World Flora
  results_1 <- results_1 %>% dplyr::select(Name_submitted, Accepted_name)

  results_2 <- TNRS(taxonomic_names = names_to_check, sources = "usda", mode = "resolve") # American database
  results_2 <- results_2 %>% dplyr::select(Name_submitted, Accepted_name)

  results_3 <- TNRS(taxonomic_names = names_to_check, sources = "tropicos", mode = "resolve") # Tropicos
  results_3 <- results_3 %>% dplyr::select(Name_submitted, Accepted_name)

  results <- cbind(results_1, results_2$Accepted_name, results_3$Accepted_name)
  colnames(results) <- c("submitted_name", "match_wfo", "match_usda", "match_tropicos")

  a <- coalesce(results$match_wfo, results$match_usda)
  checked_name <- coalesce(a, results$match_tropicos)

  results <- cbind(results, checked_name)
  results <- results %>% dplyr::select(submitted_name, checked_name)

  NFI_names_check <- left_join(df, results, by = "submitted_name")

  NFI_names_check$checked_name <- str_replace_all(NFI_names_check$checked_name, c("Populus x canescens" = "Populus canescens",
                                                                                  "Nyssa sylvatica var. biflora" = "Nyssa sylvatica",
                                                                                  "Taxodium distichum var. imbricatum" = "Taxodium distichum",
                                                                                  "Acer saccharum subsp. nigrum" = "Acer saccharum",
                                                                                  "Abies magnifica var. shastensis" = "Abies magnifica",
                                                                                  "Acer saccharum subsp. grandidentatum" = "Acer saccharum",
                                                                                  "Carya glabra var. odorata" = "Carya glabra",
                                                                                  "Pinus cembroides var. bicolor" = "Pinus cembroides",
                                                                                  "Acer saccharum subsp. leucoderme" = "Acer saccharum"))

  return(NFI_names_check)

}

