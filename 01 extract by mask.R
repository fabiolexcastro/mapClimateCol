

# Fabio Castro - Llanos 
# Alliance Bioversity - CIAT 
# December 2023

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, rgeos, gtools, glue, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Functions to use --------------------------------------------------------
ext.msk <- function(sspe, prdo){
  
  # sspe <- ssps[4]; prdo <- '2050s'
  cat('To process: ', basename(sspe), prdo, '\n')
  gcms <- dir_ls(sspe, regexp = prdo) %>% dir_ls()
  
  # To extract by mask
  rstr <- map(.x = gcms, .f = function(g){
    
    try(expr = {
      
      print(g)
      file <- dir_ls(g) %>% as.character() %>% dir_ls()
      
      # Prec
      prec <- grep('prec', file, value = T) %>% as.character() %>% rast() %>% crop(., col) %>% mask(., col)
      
      # Tmin 
      tmin <- grep('tmin', file, value = T) %>% as.character() %>% rast() %>% crop(., col) %>% mask(., col)
      
      # Tmax 
      tmax <- grep('tmax', file, value = T) %>% as.character() %>% rast() %>% crop(., col) %>% mask(., col)
      
      # Make a stack and return the file
      stck <- list(prec, tmin, tmax) %>% reduce(., c)
      return(stck)
      
    })
    
  })
  
  # Now to write the rasters 
  map(.x = 1:length(rstr), .f = function(i){
    
    try(expr = {
      
      mdl <- basename(names(rstr)[i])
      print(mdl)
      rst <- rstr[[i]]
      dir <- glue('./data/tif/cmip6/30s/{basename(sspe)}')
      dir_create(dir)
      terra::writeRaster(x = rst, filename = glue('{dir}/rstr_{mdl}_{prdo}.tif'), overwrite = TRUE)
      print('Done!')
      
      
    })
    
  })
  
}
clc.ens <- function(sspe, prdo){
  
  cat('To process: ', sspe, prdo, '\n')
  fles <- dir_ls(dirs, regexp = sspe, type = 'directory') %>% dir_ls(., regexp = prdo) 
  fles <- grep(prdo, fles, value = T)
  
  # Prec 
  prec <- map(.x = 1:length(fles), .f = function(i){
    
    print(i)
    fle <- fles[i]
    rst <- rast(fle)
    rst <- rst[[grep('prec', names(rst), value = F)]]
    return(rst)
    
  }) %>% 
    reduce(., c)
  
  prec.avrg <- map(.x = 1:12, .f = function(m){
    
    print(m)
    mnt <- ifelse(m < 10, paste0('0', m), as.character(m))
    ppt <- prec[[grep(paste0('_', mnt), names(prec))]]
    ppt <- terra::app(ppt, 'mean')
    names(ppt) <- glue('prec_{mnt}')
    return(ppt)
    
  }) %>% 
    reduce(., c)
  
  # Tmin 
  tmin <- map(.x = 1:length(fles), .f = function(i){
    
    print(i)
    fle <- fles[i]
    rst <- rast(fle)
    rst <- rst[[grep('tmin', names(rst), value = F)]]
    return(rst)
    
  }) %>% 
    reduce(., c)
  
  tmin.avrg <- map(.x = 1:12, .f = function(m){
    
    print(m)
    mnt <- ifelse(m < 10, paste0('0', m), as.character(m))
    tmn <- tmin[[grep(paste0('_', mnt), names(tmin))]]
    tmn <- terra::app(tmin, 'mean')
    names(tmn) <- glue('tmin_{mnt}')
    return(tmn)
    
  }) %>% 
    reduce(., c)
  
  # Tmax
  tmax <- map(.x = 1:length(fles), .f = function(i){
    
    print(i)
    fle <- fles[i]
    rst <- rast(fle)
    rst <- rst[[grep('tmax', names(rst), value = F)]]
    return(rst)
    
  }) %>% 
    reduce(., c)
  
  tmax.avrg <- map(.x = 1:12, .f = function(m){
    
    print(m)
    mnt <- ifelse(m < 10, paste0('0', m), as.character(m))
    tmx <- tmax[[grep(paste0('_', mnt), names(tmax))]]
    tmx <- terra::app(tmax, 'mean')
    names(tmx) <- glue('tmax_{mnt}')
    return(tmx)
    
  }) %>% 
    reduce(., c)
  
  # To write the final rasters
  terra::writeRaster(x = prec.avrg, filename = glue('./data/tif/cmip6/30s/prec_ensemble_{sspe}_{prdo}.tif'), overwrite = TRUE)
  terra::writeRaster(x = tmin.avrg, filename = glue('./data/tif/cmip6/30s/tmin_ensemble_{sspe}_{prdo}.tif'), overwrite = TRUE)
  terra::writeRaster(x = tmax.avrg, filename = glue('./data/tif/cmip6/30s/tmax_ensemble_{sspe}_{prdo}.tif'), overwrite = TRUE)
  print('Done!')
  
}

# Load data ---------------------------------------------------------------
wrld <- ne_countries(scale = 50, returnclass = 'sf')
col  <- wrld[wrld$sov_a3 == 'COL',]

# Raster data
root <- '//ALLIANCEDFS.ALLIANCE.CGIAR.ORG/data_cluster19/GLOBAL/climate/Agroclimas/data/ipcc_6ar_wcl_downscaled'
ssps <- dir_ls(root)

# To extract by mask ------------------------------------------------------
ext.msk(sspe = ssps[3], prdo = '2050')
ext.msk(sspe = ssps[4], prdo = '2050')

# To calculate the ensemble -----------------------------------------------
dirs <- './data/tif/cmip6/30s'
sspe <- basename(ssps[3])
prdo <- '2050s'

clc.ens(sspe = 'ssp_370', prdo = '2050s')


