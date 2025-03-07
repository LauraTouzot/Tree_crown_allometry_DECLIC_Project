# Landclim leaf area functions


fun_plot_LANDCLIM_LeafArea  <-  function(){
  library(XML)
  
  read_species_xml <- function(file) {
    xmlToDataFrame(file)
  }
  
  species <- read_species_xml("data/species-configuration.xml")
  
  SpecificFoliageAreaVector <- c("DECIDUOUS" = 6,
                                 "BROADLEAFEVERGREEN" = 4.2,
                                 "EVERGREENSHRUB" = 4.2,
                                 "PLASTICDECIDUOUS" = 6,
                                 "EVERGREEN" = 2.4,
                                 "HEATH" = 2.4,
                                 "HERBS" = 1,
                                 "HERBS_ANNUAL" = 1,
                                 "HERBS_SEASONAL" = 1
  )
  
  foliage_df <- data.frame("foliageType" = 1:5,
                           "FoliageWeightExponent" = c(1.43, 1.43, 1.7, 1.4, 1.56),
                           "FoliageWeightFactor" = c(0.08, 0.1, 0.06, 0.17, 0.23))
  
  leaf_df <- unique(species[species$PlantForm == "TREE", c("LeafHabit", "FoliageType")])
  Diameter <- seq(0, 100, len = 101)
  
  leafArea <- function(Diameter, SpecificFoliageArea, FoliageWeightFactor, FoliageWeightExponent) {
    SpecificFoliageArea * FoliageWeightFactor * Diameter^FoliageWeightExponent
  }
  
  leaf_df  <- leaf_df[c(5,2,4,3,1), ]
  la <- apply(leaf_df, 1,
              function(x) leafArea(Diameter,
                                   SpecificFoliageArea = SpecificFoliageAreaVector[x["LeafHabit"]],
                                   FoliageWeightFactor = foliage_df$FoliageWeightFactor[as.numeric(x["FoliageType"])],
                                   FoliageWeightExponent = foliage_df$FoliageWeightExponent[as.numeric(x["FoliageType"])]))
  
  pdf("figures/Landclim.pdf")
  matplot(la, type = "l", xlab = "Tree diameter (cm)",
          ylab = "Leaf area (m2)", lty = 1, main = "LANDCLIM Bjoern")
  
  species_groups <- apply(leaf_df, 1, function(x) species$Name[species$LeafHabit == x["LeafHabit"] & species$FoliageType == x["FoliageType"]])
  
  legend("topleft", col = 1:5, lty = 1,
         legend = sapply(species_groups, paste, collapse = " "), bty = "n",
         cex = 0.75)
  dev.off()
  
  ## plot Forclim from Bugmann PhD assuming kC2 is double side area so divided by 2 to be in same range as Landclim
  pdf("figures/Forclim.pdf")
  i  <-  1
  plot(Diameter, 0.5*12*foliage_df[i,"FoliageWeightFactor"]*Diameter^foliage_df[i, "FoliageWeightExponent"],
       xlab = "Tree diameter (cm)", ylab = "Leaf area (m2)",
       ylim = c(0, 1000), type = "l", main = "FORCLIM Bugmann PhD")
  for (i in 2:3) lines(Diameter, 0.5*12*foliage_df[i,"FoliageWeightFactor"]*Diameter^foliage_df[i, "FoliageWeightExponent"],col = i)
  
  for (i in 4:5) lines(Diameter, 0.5*6*foliage_df[i,"FoliageWeightFactor"]*Diameter^foliage_df[i, "FoliageWeightExponent"],col = i, lty = 1)
  legend("topleft", col = 1:5, lty = c(1,1,1,1, 1),
         legend = c("D1", "D2", "D3", "C4", "C5"), bty = "n",
         cex = 0.75)
  dev.off()
}
