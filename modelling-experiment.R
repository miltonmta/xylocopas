file.edit("readme.R")
# citation(package = "kernlab", lib.loc = NULL)
RStudio.Version()
R.Version()

# ****  Loading functions                            ----
source("./Auxiliary_functions.R")
# source("./ensemble.R")
# source("./Multiple_ENMs.R")
# source("./predict_enfa.R")

# file.edit(c("./Auxiliary_functions.R", "./Multiple_ENMs.R"))

# ****  Loading packages                             ----
load_pak(c("tidyverse", "raster", "rgdal", "abind", "spThin", "dismo", "kernlab", "vegan", "maps", "psych", "rJava", "dendextend", "beepr", "data.table", "adehabitatHS"))

# ***************************************************************************************
## 01.  Read aogcms models                           ----
## all climatic data were imported from huberi.Rproj workspace

##### .............. current                                     
current_list <- read_current(dir = "../huberi/data/climatic_vars/current/")

current_mat <- current_list[["matrix"]] # required at Exploratoty Factor Analysis
current     <- current_list[["raster"]] # required for salving selected variables and creating "back"
beep(2)

# ***************************************************************************************
## 02.  Variable selection                           ----

### Current

fa.parallel(current_mat[ , -c(1:2)], fa = 'fa') #scree plot
current_fa <- fa(current_mat[ , -c(1:2)], nfactors = 5, rotate = 'varimax')
loadings <- loadings(current_fa)
write.table(loadings, "./data/climatic_vars/selected/varimax_loadings.txt")
beep(2)

### Selected variables 
# bio02, bio03, bio10, bio14, bio16.

### RCPs
# #. Based on a cluster analisis, we've seleceted the following AOGCMS from all RCPs:
# # 1 == "CCS4"
# # 2 == "IPSL-CMSA-LR"
# # 3 == "MIROC-ESM"


### Reading selected variables

current_select <- stack("../huberi/data/climatic_vars/selected/current/current.grd")

beep(2)

# ***************************************************************************************
## 03.  Occurrences data                             ----

###.............. Reading data
# the names of columms in the file raw_data must be: "SPEC", "LONG", "LAT".
occur_raw <- read.table("./data/occurrences/raw_data.txt", h = T)
 
###.............. Filtering occurrences in the geographical space
sp <- gsub("C[1-9]","", occur_raw$SPEC)
sp_names <- unique(sp)
occur_thinned <- NULL
for(i in 1:length(sp_names))
{
  sp <- occur_raw[occur_raw[, 1] == sp_names[i], ]
  occur <- thin(sp, thin.par = 20, reps = 1, max.files = 1, out.dir = "./data/occurrences/", locs.thinned.list.return = T)
  occur <- occur[[1]]
  spnames_vec <- rep(sp_names[i], nrow(occur))
  occur_sp <- cbind(spnames_vec, occur)
  names(occur_sp) <- names(occur_raw)
  occur_thinned <- rbind(occur_thinned, occur_sp)
}
beep(2)

write.table(occur_thinned, "./data/occurrences/occur_thinned.txt", sep = ";", row.names = FALSE)

###.............. Extrancting bio variables based on the ocurrence
for(i in 1:length(sp_names))
{
  var <- create_var(occur_thinned[occur_thinned[, 1] == sp_names[i], ], sp_names[i])
}
beep(2)

# ***************************************************************************************
## 04.  Background Sampling                          ----

###.............. Background files 
# Creating and saving the object "back" for each studied species

var_files <- list.files("./data/occurrences/", pattern = "var", full.names = TRUE)
for(i in 1:length(var_files))
{
  var_file <- read.table(var_files[i], h = T, sep = ";")
  create_back(var_file, sp_names[i])
}
beep(2)


## plotting occurrencies                             ----
sp_names
cores <-  colorRampPalette(c("white","slateblue", "orangered", "red" ))(200)
plot(current$bio02, col = cores)

brasil <- shapefile("./arquivos/poligono/Brasil.shp")
current.BR <- mask(crop(current, brasil), brasil)

plot(current.BR[[1]], col = cores) #'a funcao 'mask' corta pela forma, nao pela extensao. Crop corta pela extensao.
map(brasil, add=T)
map(add=T)
points(occur_thinned[occur_thinned[, 1] == "Xylocopa_abbreviata", ][,-1], 
       cex = 3, pch = 15, col = "blue")
points(occur_thinned[occur_thinned[, 1] == "Xylocopa_vestita", ][,-1], 
       cex = 3, pch = 17, col = "orangered")

points(occur_thinned[occur_thinned[, 1] == "Xylocopa_truxali", ][,-1], 
       cex = 3, pch = 19, col = "black")


# ***************************************************************************************
## 05.  Creating trainning-testing subsets           ----
###.....................................

#  For variation control, the occurrence will be modeled with the same random subsets. (see cross_validation loop in Multiple_ENMs)

for (i in 1:length(sp_names)) 
{
  occur <- read.table(paste0("./data/occurrences/var_", sp_names[i] ,".txt"),  sep = ";", h = T)
  back  <- read.table(paste0("./data/occurrences/back_", sp_names[i] ,".txt"), sep = ";", h = T)
  
  cross_validation <- 10
  for (j in 1:cross_validation)
  {
    sample_occur <- sample(1:nrow(occur), round(0.75 * nrow(occur), 0))
    trainning <- prepareData(x = current_select, 
                             p = occur[sample_occur,  1:2], 
                             b = back[sample_occur,  1:2]) 
    testing   <- prepareData(x = current_select, 
                             p = occur[-sample_occur, 1:2], 
                             b = back[-sample_occur, 1:2])
    
    write.table(trainning, paste0("./data/occurrences/subsets-", sp_names[i],"/trainning", j, ".txt"), sep = ";")
    write.table(testing,   paste0("./data/occurrences/subsets-", sp_names[i],"/testing",   j, ".txt"), sep = ";")
  }
}


# ***************************************************************************************

## 06.  Modelling                                    ----
###.....................................
rm(list = ls())
source("./Multiple_ENMs-gower.R")
source("./predict_enfa.R")


## ... Loading species names 
occur_thinned <- read.table("./data/occurrences/occur_thinned.txt", sep = ";", h = T) 
sp <- gsub("C[1-9]","", occur_thinned$SPEC)
sp_names <- unique(sp)
sp_names
rm(sp, occur_thinned)

for (i in 1:length(sp_names))
{
  ###.............. Running the modedelling experiment
  result <- multiple_ENMs(occurrence       = paste0("./data/occurrences/var_",  sp_names[i], ".txt"),
                          background       = paste0("./data/occurrences/back_", sp_names[i], ".txt"),
                          biovar_current   = "../huberi/data/climatic_vars/selected/current/",
                          biovar_rcp26     = "../huberi/data/climatic_vars/selected/rcp26/",
                          biovar_rcp45     = "../huberi/data/climatic_vars/selected/rcp45/",
                          biovar_rcp60     = "../huberi/data/climatic_vars/selected/rcp60/",
                          biovar_rcp85     = "../huberi/data/climatic_vars/selected/rcp85/",
                          newvar_current   = 0,
                          newvar_rcp26     = 0,
                          newvar_rcp45     = 0,
                          newvar_rcp60     = 0,
                          newvar_rcp85     = 0,
                          trainning        = paste0("./data/occurrences/subsets-", sp_names[i], "/trainning"),
                          testing          = paste0("./data/occurrences/subsets-", sp_names[i], "/testing"),
                          AOGCMs           = c(1, 2, 3),
                          Pout             = paste0("./data/outputs/XP1/Pout/xp1_", sp_names[i], "_"),
                          cross_validation = 10)
  
  ###.............. Saving evaluation data
  write.table(result[["TPR_c"]],        paste0("./data/outputs/XP1/gower_", sp_names[i], "_TPR_current.txt"), sep = "\t", row.names = F)
  write.table(result[["Threshold_c"]],  paste0("./data/outputs/XP1/gower_", sp_names[i], "_t_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["Pred_area_c"]],  paste0("./data/outputs/XP1/gower_", sp_names[i], "_d_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["Threshold_sens"]],  paste0("./data/outputs/XP1/gower_", sp_names[i], "_t_sens_current.txt"),   sep = "\t", row.names = F)
  write.table(result[["AUC"]],  paste0("./data/outputs/XP1/gower_", sp_names[i], "_AUC_current.txt"),   sep = "\t", row.names = F)
  
  rm(result)
  gc()
}
beep(8)

## 12.  Ensembles and Final Outputs                  ----
#   .....................................................................................
source("./ensemble.R")
sp_names
for (i in 1:length(sp_names))
{
  result <- ensemble(Pout           =  "./data/outputs/XP1/Pout/",
                     Alld           =  paste0("./data/outputs/XP1/", sp_names[i], "_d_current.txt"),
                     sp             =  sp_names[i],
                     AOGCMs         =  c(1, 2, 3),
                     biovar_current =  "../huberi/data/climatic_vars/selected/current/")
  
### Saving predictions
  
### Saving Ensembles
  write.table(result[["Ensemble"]],   paste0("./data/outputs/ensembles/xp1_", sp_names[i], "_ENSEMBLES.txt"), sep = "\t", row.names = F)
  
  rm(result)
  gc()
}
beep(8)

### Plotting Ensembles
source("./Auxiliary_functions.R")
PlotEnsemble(occur     = "./data/occurrences/occur_thinned.txt",
             ensb_data = "./data/outputs/ensembles/xp1_",
             sp        = sp_names,
             output    = "./data/outputs/ensembles/xp1_plot_ind_",
             brasil    = shapefile("./arquivos/poligono/Brasil.shp"))
             # cores     = colorRampPalette(c("white", "gray90", "slateblue","orangered", "red4"))(1))

# ***************************************************************************************
## 12.  Get Range                                    ----
   
#........... calculating range.
source("./get_range.R")

#...........
for(i in 1:length(sp_names))
{
  XP1 <- get.range(biovar_current = "../huberi/data/climatic_vars/selected/current",
                   sp_name        = paste0(sp_names[i], "_"),
                   pout           = "./data/outputs/XP1/Pout/xp1_",
                   thr            = "./data/outputs/XP1/",
                   cross_val      = 10,
                   AOGCMs         = 3)
  
  write.table(XP1[["pres"]], paste0("./data/outputs/range/xp1_", sp_names[i], "_current.txt"), sep = "\t", row.names = F)
  write.table(XP1[["fut26"]], paste0("./data/outputs/range/xp1_", sp_names[i], "_rcp26.txt"), sep = "\t", row.names = F)
  write.table(XP1[["fut45"]], paste0("./data/outputs/range/xp1_", sp_names[i], "_rcp45.txt"), sep = "\t", row.names = F)
  write.table(XP1[["fut60"]], paste0("./data/outputs/range/xp1_", sp_names[i], "_rcp60.txt"), sep = "\t", row.names = F)
  write.table(XP1[["fut85"]], paste0("./data/outputs/range/xp1_", sp_names[i], "_rcp85.txt"), sep = "\t", row.names = F)
}


## 13.  Preparing analysis factors                   ----
###.....................................
i <- 1
thr_no_om<- paste0("./data/outputs/XP1/", sp_names[i], "_t_current.txt")
Allthr <- read.table(thr_no_om, sep = "\t", h = T)
Allthr

# bioclim_t <- Alld[, "bioclim"]
# maxent_t  <- Alld[, "maxent"]
# SVM_t     <- Alld[, "SVM"]
# gower_t   <- Alld[, "gower"]
# 
# dStat_c <- c(bioclim_d, maxent_d, SVM_d, gower_d)

Pout <- "./data/outputs/XP1/Pout/"
sp <- sp_names[i]

### Reading predictions data
#.......
bioclim_c     <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp), 
                                                list.files(Pout, pattern = "bioclim_c"))
)))
gower_c       <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp),
                                                list.files(Pout, pattern = "gower_c"))
)))
maxent_c      <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp), 
                                                list.files(Pout, pattern = "maxent_c"))
)))
SVM_c         <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp), 
                                                list.files(Pout, pattern = "SVM_c"))
)))

#.......
bioclim_rcp26     <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp),
                                                    list.files(Pout, pattern = "bioclim_rcp26"))
)))
gower_rcp26       <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp),
                                                    list.files(Pout, pattern = "gower_rcp26"))
)))
# Enfa_rcp26       <- stack(paste0( Pout,  (intersect(list.files(Pout, pattern = sp), 
#                                                     list.files(Pout, pattern = "Enfa_rcp26"))
#                                          )))
maxent_rcp26      <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp),
                                                    list.files(Pout, pattern = "maxent_rcp26"))
)))
SVM_rcp26         <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp),
                                                    list.files(Pout, pattern = "SVM_rcp26"))
)))

#.......
bioclim_rcp45     <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp),
                                                    list.files(Pout, pattern = "bioclim_rcp45"))
)))
gower_rcp45       <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp),
                                                    list.files(Pout, pattern = "gower_rcp45"))
)))
# Enfa_rcp45       <- stack(paste0( Pout,  (intersect(list.files(Pout, pattern = sp), 
#                                                     list.files(Pout, pattern = "Enfa_rcp45"))
#                                          )))
maxent_rcp45      <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp),
                                                    list.files(Pout, pattern = "maxent_rcp45"))
)))
SVM_rcp45         <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp),
                                                    list.files(Pout, pattern = "SVM_rcp45"))
)))

#.......
bioclim_rcp60     <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp),
                                                    list.files(Pout, pattern = "bioclim_rcp60"))
)))
gower_rcp60       <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp),
                                                    list.files(Pout, pattern = "gower_rcp60"))
)))
# Enfa_rcp60        <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp), 
#                                                     list.files(Pout, pattern = "Enfa_rcp60"))
#                                          )))
maxent_rcp60      <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp),
                                                    list.files(Pout, pattern = "maxent_rcp60"))
)))
SVM_rcp60         <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp),
                                                    list.files(Pout, pattern = "SVM_rcp60"))
)))

#.......
bioclim_rcp85     <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp),
                                                    list.files(Pout, pattern = "bioclim_rcp85"))
)))
gower_rcp85       <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp),
                                                    list.files(Pout, pattern = "gower_rcp85"))
)))
# Enfa_rcp85        <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp), 
#                                                     list.files(Pout, pattern = "Enfa_rcp85"))
#                                          )))
maxent_rcp85      <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp),
                                                    list.files(Pout, pattern = "maxent_rcp85"))
)))
SVM_rcp85         <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp),
                                                    list.files(Pout, pattern = "SVM_rcp85"))
)))


### Standardizing suitabilities
all_val <- values(all_output)
all_pad <- deconstante(all_val, "standardize", 2)

# ***************************************************************************************
## 14.  Uncertainty Evaluation                       ----


dado <- values(stack(bioclim, gower, maha,SVM,GLM))
metodo <- c(bioclim.metodo, gower.metodo, maha.metodo, SVM.metodo, GLM.metodo)
tempo <- c(bioclim.tempo, gower.tempo, maha.tempo, SVM.tempo, GLM.tempo)


MSout <- NULL
for(i in 1:nrow(dado))
{
  suit.i <- data.frame(suit=dado[i,], metodo=metodo, tempo=tempo)
  anova.i <- summary.aov(lm(suit ~ tempo+metodo%in%tempo, data=suit.i))
  MSout <- rbind(MSout, anova.i[[1]][,"Mean Sq"])
} 
#fecha fo "i"


MS <- data.frame( xyFromCell(bioclim0k, 1:ncell(bioclim0k)), MSout)colnames(MS) <- c("long","lat","tempo","residuo")
MS <- rasterFromXYZ(MS)



# ****  List of improvements to the scritp           ----

# 1. Implement validation by the checkerboards method.
# 2. Implement occurrence filtering at the ambiental space and compare with the geographical space one (spThin).
# 3. Transform maps in frequencies instead of suitabilities.
# 4. Impelement multi cores for running several models simultaneously.
# 5. Reduce code by implementing subfunctions, lopps.
# 6. Rewrite the code using tidy



