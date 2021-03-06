
ensemble <- function(Pout, 
                     Alld, 
                     sp, 
                     AOGCMs, 
                     biovar_current, 
                     newvar_current)
{
  ### Reading predictions             -----
  ###......................................
  
  #.......
  bioclim_c     <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp), 
                                                  list.files(Pout, pattern = "bioclim_c"))
                                       )))
  gower_c       <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp),
                                                  list.files(Pout, pattern = "gower_c"))
                                       )))
  # Enfa_c        <- stack(paste0( Pout, (intersect(list.files(Pout, pattern = sp), 
  #                                                 list.files(Pout, pattern = "Enfa_c"))
  #                                      )))
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

  ### Partial outputs                 -----
  ###......................................
  
  bioclim_Pout_c     <- na.omit(values(bioclim_c))
  bioclim_Pout_rcp26 <- na.omit(values(bioclim_rcp26))
  bioclim_Pout_rcp45 <- na.omit(values(bioclim_rcp45))
  bioclim_Pout_rcp60 <- na.omit(values(bioclim_rcp60))
  bioclim_Pout_rcp85 <- na.omit(values(bioclim_rcp85))
  rm(bioclim_c, bioclim_rcp26, bioclim_rcp45, bioclim_rcp60, bioclim_rcp85)
  gc()
  
  gower_Pout_c       <- na.omit(values(gower_c))
  gower_Pout_rcp26   <- na.omit(values(gower_rcp26))
  gower_Pout_rcp45   <- na.omit(values(gower_rcp45))
  gower_Pout_rcp60   <- na.omit(values(gower_rcp60))
  gower_Pout_rcp85   <- na.omit(values(gower_rcp85))
  rm(gower_c, gower_rcp26, gower_rcp45, gower_rcp60, gower_rcp85)
  gc()
  # 
  # Enfa_Pout_c       <- na.omit(values(Enfa_c))
  # Enfa_Pout_rcp26   <- na.omit(values(Enfa_rcp26))
  # Enfa_Pout_rcp45   <- na.omit(values(Enfa_rcp45))
  # Enfa_Pout_rcp60   <- na.omit(values(Enfa_rcp60))
  # Enfa_Pout_rcp85   <- na.omit(values(Enfa_rcp85))
  # rm(Enfa_c, Enfa_rcp26, Enfa_rcp45, Enfa_rcp60, Enfa_rcp85)
  # gc()
  # 
  maxent_Pout_c      <- na.omit(values(maxent_c))
  maxent_Pout_rcp26  <- na.omit(values(maxent_rcp26))
  maxent_Pout_rcp45  <- na.omit(values(maxent_rcp45))
  maxent_Pout_rcp60  <- na.omit(values(maxent_rcp60))
  maxent_Pout_rcp85  <- na.omit(values(maxent_rcp85))
  rm(maxent_c, maxent_rcp26, maxent_rcp45, maxent_rcp60, maxent_rcp85)
  gc()

  SVM_Pout_c         <- na.omit(values(SVM_c))
  SVM_Pout_rcp45     <- na.omit(values(SVM_rcp45))
  SVM_Pout_rcp26     <- na.omit(values(SVM_rcp26))
  SVM_Pout_rcp60     <- na.omit(values(SVM_rcp60))
  SVM_Pout_rcp85     <- na.omit(values(SVM_rcp85))
  rm(SVM_c, SVM_rcp26, SVM_rcp45, SVM_rcp60, SVM_rcp85)
  gc()
  
  ### Saving the predictive methods   -----
  ###......................................
  
  output_current <- cbind(Bioclim = bioclim_Pout_c,
                          Gower   = gower_Pout_c,
                          # Enfa    = Enfa_Pout_c,
                          Maxent  = maxent_Pout_c,
                          SVM     = SVM_Pout_c )
  rm(bioclim_Pout_c, maxent_Pout_c, SVM_Pout_c, gower_Pout_c)
  gc()
  
  output_rcp26   <- cbind(Bioclim = bioclim_Pout_rcp26,
                          Gower   = gower_Pout_rcp26,
                          # Enfa    = Enfa_Pout_rcp26,
                          Maxent  = maxent_Pout_rcp26,
                          SVM     = SVM_Pout_rcp26 )
  rm(bioclim_Pout_rcp26, maxent_Pout_rcp26, SVM_Pout_rcp26, gower_Pout_rcp26)
  gc()

  output_rcp45   <- cbind(Bioclim = bioclim_Pout_rcp45,
                          Gower   = gower_Pout_rcp45,
                          # Enfa    = Enfa_Pout_rcp45,
                          Maxent  = maxent_Pout_rcp45,
                          SVM     = SVM_Pout_rcp45 )
  rm(bioclim_Pout_rcp45, maxent_Pout_rcp45, SVM_Pout_rcp45, gower_Pout_rcp45)
  gc()

  output_rcp60   <- cbind(Bioclim = bioclim_Pout_rcp60,
                          Gower   = gower_Pout_rcp60,
                          # Enfa    = Enfa_Pout_rcp60,
                          Maxent  = maxent_Pout_rcp60,
                          SVM     = SVM_Pout_rcp60 )
  rm(bioclim_Pout_rcp60, maxent_Pout_rcp60, SVM_Pout_rcp60, gower_Pout_rcp60)
  gc()

  output_rcp85   <- cbind(Bioclim = bioclim_Pout_rcp85,
                          Gower   = gower_Pout_rcp85,
                          # Enfa    = Enfa_Pout_rcp85,
                          Maxent  = maxent_Pout_rcp85,
                          SVM     = SVM_Pout_rcp85 )
  rm(bioclim_Pout_rcp85, maxent_Pout_rcp85, SVM_Pout_rcp85, gower_Pout_rcp85)
  gc()

  
  ### Ensembles                       ----
  ###.....................................
  
  ## ..... Standardizing predictions
  id_col_fut <- rep(1:ncol(output_current), each = 3)
  id_time    <- c(rep("c", nrow(output_current)), rep(c("rcp26", "rcp45", "rcp60", "rcp85"),
                                                      each = nrow(output_current) * length(AOGCMs)))
  
  pad_c <- pad_rcp26 <- pad_rcp45 <- pad_rcp60 <- pad_rcp85 <- NULL
  for(p in 1:ncol(output_current))
  {
    suit     <- cbind(output_current[, p],
                      output_rcp26[, which(id_col_fut == p)],
                      output_rcp45[, which(id_col_fut == p)],
                      output_rcp60[, which(id_col_fut == p)],
                      output_rcp85[, which(id_col_fut == p)])
    suit     <- as.numeric(suit)
    suit_pad <- decostand(x = suit, method = "range") # requires vegan
    
    pad_c     <- cbind(pad_c, suit_pad[which(id_time == "c"), 1])
    
    pad_rcp26 <- cbind(pad_rcp26, matrix(data = suit_pad[which(id_time == "rcp26"), 1],
                                         nrow = nrow(output_current), ncol = length(AOGCMs)))

    pad_rcp45 <- cbind(pad_rcp45, matrix(data = suit_pad[which(id_time == "rcp45"), 1],
                                         nrow = nrow(output_current), ncol = length(AOGCMs)))

    pad_rcp60 <- cbind(pad_rcp60, matrix(data = suit_pad[which(id_time == "rcp60"), 1],
                                         nrow = nrow(output_current), ncol = length(AOGCMs)))

    pad_rcp85 <- cbind(pad_rcp85, matrix(data = suit_pad[which(id_time == "rcp85"), 1],
                                         nrow = nrow(output_current), ncol = length(AOGCMs)))
  }
  rm(output_current, output_rcp26, output_rcp45, output_rcp60, output_rcp85)
  gc()
  
  ### checar e incluir os outros dois novos métodos daqui para baixo..
  pad_rcp26_ls <- list(pad_rcp26[, 1:30], pad_rcp26[, 31:60], pad_rcp26[, 61:90])
  pad_rcp45_ls <- list(pad_rcp45[, 1:30], pad_rcp45[, 31:60], pad_rcp45[, 61:90])
  pad_rcp60_ls <- list(pad_rcp60[, 1:30], pad_rcp60[, 31:60], pad_rcp60[, 61:90])
  pad_rcp85_ls <- list(pad_rcp85[, 1:30], pad_rcp85[, 31:60], pad_rcp85[, 61:90])
  
  ## ..... Calculating Ensembles
  

  Alld <- read.table(Alld, sep = "\t", h = T)
  bioclim_d <- Alld[, "bioclim"]
  maxent_d  <- Alld[, "maxent"]
  SVM_d     <- Alld[, "SVM"]
  gower_d   <- Alld[, "gower"]
  
  dStat_c <- c(bioclim_d, maxent_d, SVM_d, gower_d)
  dStat_fut <- rep(dStat_c, each = length(AOGCMs))
  
  ensemble_c     <- apply(pad_c,     1, function(x) sum(x * dStat_c)   / sum(dStat_c))
  rm(pad_c)
  gc()
  ensemble_rcp26 <- apply(pad_rcp26, 1, function(x) sum(x * dStat_fut) / sum(dStat_fut))
  rm(pad_rcp26)
  gc()
  ensemble_rcp45 <- apply(pad_rcp45, 1, function(x) sum(x * dStat_fut) / sum(dStat_fut))
  rm(pad_rcp45)
  gc()
  ensemble_rcp60 <- apply(pad_rcp60, 1, function(x) sum(x * dStat_fut) / sum(dStat_fut))
  rm(pad_rcp60)
  gc()
  ensemble_rcp85 <- apply(pad_rcp85, 1, function(x) sum(x * dStat_fut) / sum(dStat_fut))
  rm(pad_rcp85)
  gc()
  
  ensemble_rcp26_ls <- ensemble_rcp45_ls <- ensemble_rcp60_ls <- ensemble_rcp85_ls <- list()
  for (i in AOGCMs)
  {
    ensemble_rcp26_ls[[i]] <- apply(pad_rcp26_ls[[i]], 1, function(x) sum(x * dStat_c) / sum(dStat_c))
    gc()
    ensemble_rcp45_ls[[i]] <- apply(pad_rcp45_ls[[i]], 1, function(x) sum(x * dStat_c) / sum(dStat_c))
    gc()
    ensemble_rcp60_ls[[i]] <- apply(pad_rcp60_ls[[i]], 1, function(x) sum(x * dStat_c) / sum(dStat_c))
    gc()
    ensemble_rcp85_ls[[i]] <- apply(pad_rcp85_ls[[i]], 1, function(x) sum(x * dStat_c) / sum(dStat_c))
    gc()
  }
  rm(pad_rcp26_ls, pad_rcp45_ls, pad_rcp60_ls, pad_rcp85_ls)
  gc()
  
  ### Saving Ensembles with coords    ----
  ###.....................................
  
  current <- stack(list.files(biovar_current,  pattern = ".grd$", full.names = TRUE))
  
  coords <- na.omit(cbind(xyFromCell(current, 1:ncell(current)), values(current)))[,1:2]
  rm(current)
  gc()
 
  # output_current <- rasterFromXYZ(cbind(coords, output_current))
  # gc()
  # output_rcp26   <- rasterFromXYZ(cbind(coords, output_rcp26))
  # gc()
  # output_rcp45   <- rasterFromXYZ(cbind(coords, output_rcp45))
  # gc()
  # output_rcp60   <- rasterFromXYZ(cbind(coords, output_rcp60))
  # gc()
  # output_rcp85   <- rasterFromXYZ(cbind(coords, output_rcp85))
  # gc()

  FULLensemble       <- data.frame(coords,
                                   Ensemble_present = ensemble_c,
                                   Ensemble_rcp26   = ensemble_rcp26,
                                   Ensemble_rcp45   = ensemble_rcp45,
                                   Ensemble_rcp60   = ensemble_rcp60,
                                   Ensemble_rcp85   = ensemble_rcp85,
                                   Ensemble_rcp26.1 = ensemble_rcp26_ls[[1]],
                                   Ensemble_rcp26.2 = ensemble_rcp26_ls[[2]],
                                   Ensemble_rcp26.3 = ensemble_rcp26_ls[[3]],
                                   Ensemble_rcp45.1 = ensemble_rcp45_ls[[1]],
                                   Ensemble_rcp45.2 = ensemble_rcp45_ls[[2]],
                                   Ensemble_rcp45.3 = ensemble_rcp45_ls[[3]],
                                   Ensemble_rcp60.1 = ensemble_rcp60_ls[[1]],
                                   Ensemble_rcp60.2 = ensemble_rcp60_ls[[2]],
                                   Ensemble_rcp60.3 = ensemble_rcp60_ls[[3]],
                                   Ensemble_rcp85.1 = ensemble_rcp85_ls[[1]],
                                   Ensemble_rcp85.2 = ensemble_rcp85_ls[[2]],
                                   Ensemble_rcp85.3 = ensemble_rcp85_ls[[3]])
                                   
  # ### Returning Function data         ----
  # ###.....................................
  # 
  return(list(#"output_current"  = output_current,
              # "output_rcp26"    = output_rcp26,
              # "output_rcp45"    = output_rcp45,
              # "output_rcp60"    = output_rcp60,
              # "output_rcp85"    = output_rcp85,
              "Ensemble"        = FULLensemble))
}
  