# Tree Crown Allometry - DECLIC Project

Step 1: compiling the database
The code for compiling the database is available in the 01_Database_compilation folder. 
Data have not been added to the folder, but the links and references to the various databases, repositories or data papers are available in Supporting Information A; Table S1 (see also below). 
The _targets.R file provides the entire pipeline for compiling and cleaning the database, while the functions required to implement it are available in full in the R folder. 


Step 2: fitting and selecting allometry models
The code for fitting all considered allometry models (i.e., asymptotic, power-law and beta-regression models, with and without competition, using either BAT or BAL as competition index and whether weighting or not the protocol effects) is available in the 02_Allometry_relationships folder. 
Data has not been added to the folder but can be computed in Step 1. The database compiled in the first step is used as it is in Step 2. 

Step 3: performing all analyses
The code, data and output to perform all analyses presented in the manuscript, along with producing the figures and tables, are provided in the main folder and organized as required by the used functions. 

In addition, mean parameters and their standard deviations obtained for each studied species and crown characteristic, with and without competition (based on the selected models) are provided in files mean_parameters_nocompetition.csv and mean_parameters_competition.csv. 

Package targets is needed to run all scripts according to the pipelines available in _targets.R files.


Databases, repositories and data papers' references: 
- BAAD: a Biomass and Allometry Database for woody plants (https://github.com/dfalster/baad)
- FHM: Forest Health Monitoring (https://www.fs.usda.gov/foresthealth/protecting-forest/forest-health-monitoring/index.shtml)
- FIA: Forest Inventory and Analysis (https://www.fs.usda.gov/research/inventory/FIA)
- Spanish National Forest Inventory (https://www.miteco.gob.es/en/biodiversidad/servicios/banco-datos-naturaleza/informacion-disponible/ifn3_bbdd_descargas.htm.aspx; https://www.miteco.gob.es/en/biodiversidad/servicios/banco-datos-naturaleza/informacion-disponible/ifn2_descargas.aspx)
- FunDivEUROPE: Functional significance of forest biodiversity (Germany: https://bwi.info/Download/de/BWI-Basisdaten/ACCESS2003/; Finland and Sweden:
https://doi.org/10.5061/dryad.wm37p⟨200b⟩vmkw)
- MONTANE (https://doi.org/10.57745/GFITOT)
- GenTree (Opgenoorth, L., Dauphin, B., Benavides, R., Heer, K., Alizoti, P., Martínez-Sancho, E., ... & Cavers, S. (2021). The GenTree Platform: growth traits and tree-level environmental data in 12 European forest tree species. GigaScience)
- LegacyTreeData (https://legacytreedata.org/)
- ICP Forests (International Co-operative Programme on Air Pollution Effects on Forests) (https://icp-forests.org/documentation/Introduction/index.html)
- EuMIXFOR (https://doi.org/10.5061/dryad.8v04m)
- French National Forest Inventory (https://inventaire-forestier.ign.fr/dataifn/)
- Canada’s National Forest Inventory (https://nfi.nfis.org/en/datarequest)
- Quebec’s National Forest Inventory (inventaires.forestiers@mffp.gouv.qc.ca)
- Tallo (https://zenodo.org/record/6637599#.ZFoYT-xBw-Q)
- Sample tree biomass data for Eurasian forests (https://elar.usfeu.ru/handle/123456789/4931)
- Harvard Forest Data Archive (https://harvardforest1.fas.harvard.edu/exist/apps/datasets/showData.html?id=HF339)
- Data papers
    - Fuhr, M., Cordonnier, T., Courbaud, B., Kunstler, G., Mermin, E., Riond, C., & Tardif, P. (2017). Long‐term tree inventory data from mountain forest plots in France.
    - Evans, M. R., Moustakas, A., Carey, G., Malhi, Y., Butt, N., Benham, S., Pallett, D. & Schäfer, S. (2015). Allometry and growth of eight tree taxa in United Kingdom woodlands. Scientific Data.
    - Dettmann, G. T., & MacFarlane, D. W. (2019). Trans‐species predictors of tree leaf mass. Ecological Applications.
    - Anderson‐Teixeira, K. J., McGarvey, J. C., Muller‐Landau, H. C., Park, J. Y., Gonzalez‐Akre, E. B., Herrmann, V., ... & McShea, W. J. (2015). Size‐related scaling of tree form and function in a mixed‐age forest. Functional Ecology.
    - Dalponte, M., & Coomes, D. A. (2016). Tree‐centric mapping of forest carbon density from airborne laser scanning and hyperspectral data. Methods in Ecology and Evolution. 
