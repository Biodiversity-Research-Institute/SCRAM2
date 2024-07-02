## Script to take species output from one-state model (SCRAM 2) and place into SF 
#  (to combine with GPS model outputs into ensemble, in the case of REKN)

##Author: Holly Goyert
##Date: 05/30/2024

# #R v4.3.1

## load required packages
library(R2jags)
library(R2WinBUGS)
library(tidyverse)
library(sf)
library(raster)
library(rgeos)
library(MCMCvis) #MCMCsummary rather than summary
try(library(motus)) #to install: browseURL("https://motuswts.github.io/motus/articles/02-installing-packages.html")

## set directories
#dir.home <- file.path('P:') #Home directory (network)
#dir.root <- file.path(dir.home, 'SCRAM', 'BRI') #Root directory (network)

usr.name <- 'user.name' #hard-coded (needs manual update)!
dir.home <- file.path('C:', 'Users', usr.name, 'OneDrive - Biodiversity Research Institute', 'Documents') #Home directory (local)
dir.root <- file.path(dir.home, 'ProjectsOffshoreWind', 'Atlantic', 'USFWS', 'SCRAM', 'BRI') #Root directory (local)

dir.proj <- file.path(dir.root, 'Integrated Movement Models') #Project directory
dir.data <- file.path(dir.proj, 'data') #GIS data directory

dir.work <- dir.proj
setwd(dir.work) # setwd <- dir.work 

## REKN ####
spp <- 'REKN' #hard-coded (needs manual update)!
inc <- 24 #hard-coded (needs manual update)! refers to the increment used to fit to detections and fit to regular interval increments
out.name <- paste0(spp, '_2023') #name assigned to output files, refers to year Phase 2 started
out.nam2 <- paste0('_2021_1state_0D_var_hr',inc) #hard-coded (needs manual update)! 
#2021 refers the version of species data (release year, from Loring et al. 2021)
#1state refers to the single state movement model
#0D refers to the exclusion of a drift parameter
#var refers to the conversion of the variance-covariance matrix to independent normal variance
#n240 retains only fall migration (removes spring migratory trajectories prior to 6-28; keeps post 6-21)

## to load/read data for JAGS from 1-state model ####
N <- 240 #hard-coded (needs manual update)!

dir.out <- file.path(dir.proj, 'Motus Movement Models', spp, paste0(out.name, out.nam2)) #subfolder for output files

# dir.local <- file.path('P:') #server (or local data drive on high performance computer)
# dir.local <- file.path(dir.local, 'SCRAM', 'BRI', 'Integrated Movement Models', 'Motus Movement Models', spp, paste0(out.name, out.nam2)) 
# 
# load(file=file.path(dir.out, paste0('occ_post_', out.name, '_studyarea', '.RData')))
# 
# str(occ_post)
# #num [1:1000, 1:12, 1:512]
# #iterations, months, grid cells
# 
# ## load shapefile ####
# occup_BOEM_halfdeg_grid_aea_sf <- readRDS(file.path('P:', 'SCRAM', 'BRI', 'Movement', 'SCRAM Movement Models', 'occup_BOEM_halfdeg_grid_aea_sf_motus_gps1.RDS'))
# str(occup_BOEM_halfdeg_grid_aea_sf)
# # Classes ‘sf’ and 'data.frame':	512 obs. of  18 variables:
# #   $ id                : chr  "1" "2" "3" "4" ...
# # $ lon               : num  -66.7 -67.2 -67.7 -68.2 -68.2 ...
# # $ lat               : num  44.5 44.5 44.5 44.5 44 ...
# # $ geometry          :sfc_POLYGON of length 512; first list element: List of 1
# # ..$ : num [1:5, 1:2] 2142950 2125153 2160613 2178708 2142950 ...
# # ..- attr(*, "class")= chr [1:3] "XY" "POLYGON" "sfg"
# # $ area_sqkm         : num  2209 2209 2209 2209 2227 ...
# # $ mean_cell_width_km: num  47 47 47 47 47.2 ...
# # $ Jan               : num [1:512, 1:3, 1:1000] NA NA NA NA NA NA NA NA NA NA ...
# # $ Feb               : num [1:512, 1:3, 1:1000] NA NA NA NA NA NA NA NA NA NA ...
# # $ Mar               : num [1:512, 1:3, 1:1000] NA NA NA NA NA NA NA NA NA NA ...
# # $ Apr               : num [1:512, 1:3, 1:1000] NA NA NA NA NA NA NA NA NA NA ...
# # $ May               : num [1:512, 1:3, 1:1000] NA NA NA NA NA NA NA NA NA NA ...
# # $ Jun               : num [1:512, 1:3, 1:1000] NA NA NA NA NA NA NA NA NA NA ...
# # $ Jul               : num [1:512, 1:3, 1:1000] NA NA NA NA NA NA NA NA NA NA ...
# # $ Aug               : num [1:512, 1:3, 1:1000] NA NA NA NA NA NA NA NA NA NA ...
# # $ Sep               : num [1:512, 1:3, 1:1000] NA NA NA NA NA NA NA NA NA NA ...
# # $ Oct               : num [1:512, 1:3, 1:1000] NA NA NA NA NA NA NA NA NA NA ...
# # $ Nov               : num [1:512, 1:3, 1:1000] NA NA NA NA NA NA NA NA NA NA ...
# # $ Dec               : num [1:512, 1:3, 1:1000] NA NA NA NA NA NA NA NA NA NA ...
# # - attr(*, "sf_column")= chr "geometry"
# # - attr(*, "agr")= Factor w/ 3 levels "constant","aggregate",..: NA NA NA NA NA NA NA NA NA NA ...
# # ..- attr(*, "names")= chr [1:17] "id" "lon" "lat" "area_sqkm" ...
# 
# monthLabels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",  "Aug", "Sep", "Oct", "Nov", "Dec")
# 
# dim(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][['Jan']])
# dim(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][[1]]) #same
# #512    3 1000
# #cell,mod,iteration
# 
# dim(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][[1]]) #array
# dim(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][1]) #data.frame (first values of array dimensions & geometry)
# #512   2
# 
# sum(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][[7]][,1,1]!= #same
# st_drop_geometry(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels])[,7][,1,1], na.rm=T)
# 
# occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][7] 
# occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][,7] #should be same?
# 
# dim(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels]) #contains geometry as 13
# #512  13
# 
# dim(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][['Jan']][,1,]) #first model (Motus)
# #512 1000
# 
# dim(t(occ_post[,1,])) #Jan
# #512 1000
# 
# for (i in 1:length(monthLabels)) {
# occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][[i]][,1,] <- t(occ_post[,i,])
# }
# 
# #calculate mean across iterations and models for months with data/predictions
# m=7 #test
# dim(apply(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][[m]], c(1,2), function(x) mean(x, na.rm = TRUE)))
# #512   3
# occup_BOEM_halfdeg_grid_aea_sf_mean <- occup_BOEM_halfdeg_grid_aea_sf
# 
# for (m in which(rowSums(apply(occ_post, c(2,3), function(x) sum(x, na.rm = TRUE)), na.rm = TRUE)>0)) {
#   occup_BOEM_halfdeg_grid_aea_sf_mean[,monthLabels][[m]][,,1] <- apply(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][[m]], c(1,2), function(x) mean(x, na.rm = TRUE))
# }
# 
# #Plot months with data/predictions (need mean otherwise plots first of 1000 iterations)
# for (m in which(rowSums(apply(occ_post, c(2,3), function(x) sum(x, na.rm = TRUE)), na.rm = TRUE)>0)) {
#   plot(occup_BOEM_halfdeg_grid_aea_sf_mean[,monthLabels][m], logz = TRUE)
# }
# #comparable to P:\SCRAM\BRI\Movement\SCRAM Movement Models\winter2023results\rekn_2023_2021_1state_0D_var_n240\occ_post_rekn_2023_med_studyarea.pdf
# 
# saveRDS(occup_BOEM_halfdeg_grid_aea_sf, file.path('P:', 'SCRAM', 'BRI', 'Integrated Movement Models', 'occup_BOEM_halfdeg_grid_aea_sf_rekn_motus_gps2.RDS')) 

load(file=file.path(dir.out, paste0('occ_post_', out.name, '_studyarea', '.RData')))
dir.local <- file.path('F:') #server (or local data drive on high performance computer)
dir.local <- file.path(dir.local, 'SCRAM', 'BRI', 'Integrated Movement Models', 'Motus Movement Models', spp, paste0(out.name, out.nam2)) 
dir.create(dir.local)

## Plot all months and export shapefile for study area 
BOEM.sf <- st_read('P:/SCRAM/BRI/Integrated Movement Models/data',
                   'BOEM_halfdeg_grid_latlon_2021')
occu_med <- apply(occ_post, c(3, 2), function(x) mean(x, na.rm = TRUE))

BOEM.sf <- cbind(BOEM.sf, occu_med)
BOEM.sf <- BOEM.sf[BOEM.sf$top > 36.4914 & BOEM.sf$bottom < 42.2453,] #reduce to study area extent for mapping
#str(BOEM.sf) 

monthLabels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",  "Aug", "Sep", "Oct", "Nov", "Dec")

for (m in monthLabels[6:11]) { #which(N_bymonth>0)
  dir.create(file.path(dir.local, paste0('occ_post_', out.name, '_', m, '_med_shp')))
  try(st_write(BOEM.sf[paste0(m, '')], file.path(dir.local, paste0('occ_post_', out.name, '_', m, '_med_shp'), paste0('occ_post_', out.name, '_', m, '_med.shp')), append=F))
  jpeg(filename = file.path(dir.local, paste0('occ_post_', out.name, '_', m, '_med.jpg')))
  try(plot(BOEM.sf[paste0(m, '')], logz = TRUE))
  dev.off()
  print(summary(BOEM.sf[paste0(m, '')], na.rm=T))
}

# Deleting layer `occ_post_REKN_2023_Jun_med' using driver `ESRI Shapefile'
# Writing layer `occ_post_REKN_2023_Jun_med' to data source 
# `F:/SCRAM/BRI/Integrated Movement Models/Motus Movement Models/REKN/REKN_2023_2021_1state_0D_var_hr24/occ_post_REKN_2023_Jun_med_shp/occ_post_REKN_2023_Jun_med.shp' using driver `ESRI Shapefile'
# Writing 202 features with 1 fields and geometry type Polygon.
#       Jun               geometry  
#  Min.   : NA   POLYGON      :202  
#  1st Qu.: NA   epsg:4326    :  0  
#  Median : NA   +proj=long...:  0  
#  Mean   :NaN                      
#  3rd Qu.: NA                      
#  Max.   : NA                      
#  NA's   :202                      
# Deleting layer `occ_post_REKN_2023_Jul_med' using driver `ESRI Shapefile'
# Writing layer `occ_post_REKN_2023_Jul_med' to data source 
# `F:/SCRAM/BRI/Integrated Movement Models/Motus Movement Models/REKN/REKN_2023_2021_1state_0D_var_hr24/occ_post_REKN_2023_Jul_med_shp/occ_post_REKN_2023_Jul_med.shp' using driver `ESRI Shapefile'
# Writing 202 features with 1 fields and geometry type Polygon.
#       Jul                    geometry  
#  Min.   :0.000000   POLYGON      :202  
#  1st Qu.:0.000000   epsg:4326    :  0  
#  Median :0.000000   +proj=long...:  0  
#  Mean   :0.014028                      
#  3rd Qu.:0.003014                      
#  Max.   :0.615167                      
# Deleting layer `occ_post_REKN_2023_Aug_med' using driver `ESRI Shapefile'
# Writing layer `occ_post_REKN_2023_Aug_med' to data source 
# `F:/SCRAM/BRI/Integrated Movement Models/Motus Movement Models/REKN/REKN_2023_2021_1state_0D_var_hr24/occ_post_REKN_2023_Aug_med_shp/occ_post_REKN_2023_Aug_med.shp' using driver `ESRI Shapefile'
# Writing 202 features with 1 fields and geometry type Polygon.
#       Aug                     geometry  
#  Min.   :0.0000000   POLYGON      :202  
#  1st Qu.:0.0002353   epsg:4326    :  0  
#  Median :0.0014559   +proj=long...:  0  
#  Mean   :0.0447687                      
#  3rd Qu.:0.0060607                      
#  Max.   :2.4827132                      
# Deleting layer `occ_post_REKN_2023_Sep_med' using driver `ESRI Shapefile'
# Writing layer `occ_post_REKN_2023_Sep_med' to data source 
# `F:/SCRAM/BRI/Integrated Movement Models/Motus Movement Models/REKN/REKN_2023_2021_1state_0D_var_hr24/occ_post_REKN_2023_Sep_med_shp/occ_post_REKN_2023_Sep_med.shp' using driver `ESRI Shapefile'
# Writing 202 features with 1 fields and geometry type Polygon.
#       Sep                    geometry  
#  Min.   :0.000028   POLYGON      :202  
#  1st Qu.:0.001655   epsg:4326    :  0  
#  Median :0.005301   +proj=long...:  0  
#  Mean   :0.069670                      
#  3rd Qu.:0.018613                      
#  Max.   :4.695648                      
# Deleting layer `occ_post_REKN_2023_Oct_med' using driver `ESRI Shapefile'
# Writing layer `occ_post_REKN_2023_Oct_med' to data source 
# `F:/SCRAM/BRI/Integrated Movement Models/Motus Movement Models/REKN/REKN_2023_2021_1state_0D_var_hr24/occ_post_REKN_2023_Oct_med_shp/occ_post_REKN_2023_Oct_med.shp' using driver `ESRI Shapefile'
# Writing 202 features with 1 fields and geometry type Polygon.
#       Oct                    geometry  
#  Min.   :0.000000   POLYGON      :202  
#  1st Qu.:0.000559   epsg:4326    :  0  
#  Median :0.002697   +proj=long...:  0  
#  Mean   :0.103348                      
#  3rd Qu.:0.011072                      
#  Max.   :7.106189                      
# Deleting layer `occ_post_REKN_2023_Nov_med' using driver `ESRI Shapefile'
# Writing layer `occ_post_REKN_2023_Nov_med' to data source 
# `F:/SCRAM/BRI/Integrated Movement Models/Motus Movement Models/REKN/REKN_2023_2021_1state_0D_var_hr24/occ_post_REKN_2023_Nov_med_shp/occ_post_REKN_2023_Nov_med.shp' using driver `ESRI Shapefile'
# Writing 202 features with 1 fields and geometry type Polygon.
#       Nov                     geometry  
#  Min.   :0.0000000   POLYGON      :202  
#  1st Qu.:0.0000805   epsg:4326    :  0  
#  Median :0.0017356   +proj=long...:  0  
#  Mean   :0.0707748                      
#  3rd Qu.:0.0164914                      
#  Max.   :3.0178391

#Oct Max.   :7.10619 

## PIPL ####
spp <- 'PIPL' #hard-coded (needs manual update)!
inc <- 24 #hard-coded (needs manual update)! refers to the increment used to fit to detections and fit to regular interval increments
out.name <- paste0(spp, '_2023') #name assigned to output files, refers to year Phase 2 started
out.nam2 <- paste0('_2018_1state_0D_var_hr',inc) #hard-coded (needs manual update)! 
#2018 refers the version of species data (release year, from Loring et al. 2019)
#1state refers to the single state movement model
#0D refers to the exclusion of a drift parameter
#var refers to the conversion of the variance-covariance matrix to independent normal variance

## to load/read data for JAGS from 1-state model ####
N <- 107 #hard-coded (needs manual update)!

dir.out <- file.path(dir.proj, 'Motus Movement Models', spp, paste0(out.name, out.nam2)) #subfolder for output files

# dir.local <- file.path('P:') #server (or local data drive on high performance computer)
# dir.local <- file.path(dir.local, 'SCRAM', 'BRI', 'Integrated Movement Models', 'Motus Movement Models', spp, paste0(out.name, out.nam2)) 
# 
# load(file=file.path(dir.out, paste0('occ_post_', out.name, '_studyarea', '.RData')))
# 
# str(occ_post)
# #num [1:1000, 1:12, 1:512]
# #iterations, months, grid cells
# 
# ## load shapefile ####
# occup_BOEM_halfdeg_grid_aea_sf <- readRDS(file.path('P:', 'SCRAM', 'BRI', 'Movement', 'SCRAM Movement Models', 'occup_BOEM_halfdeg_grid_aea_sf.RDS'))
# str(occup_BOEM_halfdeg_grid_aea_sf)
# # Classes ‘sf’ and 'data.frame':	512 obs. of  18 variables:
# 
# monthLabels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",  "Aug", "Sep", "Oct", "Nov", "Dec")
# 
# dim(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][['Jan']])
# dim(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][[1]]) #same
# #512    3 1000
# #cell,mod,iteration
# 
# dim(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][[1]]) #array
# dim(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][1]) #data.frame (first values of array dimensions & geometry)
# #512   2
# 
# sum(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][[7]][,1,1]!= #same
#       st_drop_geometry(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels])[,7][,1,1], na.rm=T)
# 
# occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][7] 
# occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][,7] #should be same?
# 
# dim(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels]) #contains geometry as 13
# #512  13
# 
# dim(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][['Jan']][,1,]) #first model (Motus)
# #512 1000
# 
# dim(t(occ_post[,1,])) #Jan
# #512 1000
# 
# for (i in 1:length(monthLabels)) {
#   occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][[i]][,1,] <- t(occ_post[,i,])
# }
# 
# #calculate mean across iterations and models for months with data/predictions
# m=7 #test
# dim(apply(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][[m]], c(1,2), function(x) mean(x, na.rm = TRUE)))
# #512   3
# occup_BOEM_halfdeg_grid_aea_sf_mean <- occup_BOEM_halfdeg_grid_aea_sf
# 
# for (m in which(rowSums(apply(occ_post, c(2,3), function(x) sum(x, na.rm = TRUE)), na.rm = TRUE)>0)) {
#   occup_BOEM_halfdeg_grid_aea_sf_mean[,monthLabels][[m]][,,1] <- apply(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][[m]], c(1,2), function(x) mean(x, na.rm = TRUE))
# }
# 
# #Plot months with data/predictions (need mean otherwise plots first of 1000 iterations)
# for (m in which(rowSums(apply(occ_post, c(2,3), function(x) sum(x, na.rm = TRUE)), na.rm = TRUE)>0)) {
#   plot(occup_BOEM_halfdeg_grid_aea_sf_mean[,monthLabels][m], logz = TRUE)
# }
# #comparable to P:\SCRAM\BRI\Movement\SCRAM Movement Models\winter2023results\pipl_2023_2018_1state_0D_var\occ_post_pipl_2023_med.pdf
# 
# saveRDS(occup_BOEM_halfdeg_grid_aea_sf, file.path('P:', 'SCRAM', 'BRI', 'Integrated Movement Models', 'occup_BOEM_halfdeg_grid_aea_sf_pipl_motus2.RDS')) 

load(file=file.path(dir.out, paste0('occ_post_', out.name, '_studyarea', '.RData')))
dir.local <- file.path('F:') #server (or local data drive on high performance computer)
dir.local <- file.path(dir.local, 'SCRAM', 'BRI', 'Integrated Movement Models', 'Motus Movement Models', spp, paste0(out.name, out.nam2)) 
dir.create(dir.local)

## Plot all months and export shapefile for study area 
BOEM.sf <- st_read('P:/SCRAM/BRI/Integrated Movement Models/data',
                   'BOEM_halfdeg_grid_latlon_2021')
occu_med <- apply(occ_post, c(3, 2), function(x) mean(x, na.rm = TRUE))

BOEM.sf <- cbind(BOEM.sf, occu_med)
BOEM.sf <- BOEM.sf[BOEM.sf$top > 36.4914 & BOEM.sf$bottom < 42.2453,] #reduce to study area extent for mapping
#str(BOEM.sf) 

monthLabels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",  "Aug", "Sep", "Oct", "Nov", "Dec")

for (m in monthLabels[5:9]) { #which(N_bymonth>0)
  dir.create(file.path(dir.local, paste0('occ_post_', out.name, '_', m, '_med_shp')))
  try(st_write(BOEM.sf[paste0(m, '')], file.path(dir.local, paste0('occ_post_', out.name, '_', m, '_med_shp'), paste0('occ_post_', out.name, '_', m, '_med.shp')), append=F))
  jpeg(filename = file.path(dir.local, paste0('occ_post_', out.name, '_', m, '_med.jpg')))
  try(plot(BOEM.sf[paste0(m, '')], logz = TRUE))
  dev.off()
  print(summary(BOEM.sf[paste0(m, '')], na.rm=T))
}

# Deleting layer `occ_post_PIPL_2023_May_med' using driver `ESRI Shapefile'
# Writing layer `occ_post_PIPL_2023_May_med' to data source 
# `F:/SCRAM/BRI/Integrated Movement Models/Motus Movement Models/PIPL/PIPL_2023_2018_1state_0D_var_hr24/occ_post_PIPL_2023_May_med_shp/occ_post_PIPL_2023_May_med.shp' using driver `ESRI Shapefile'
# Writing 202 features with 1 fields and geometry type Polygon.
#       May                   geometry  
#  Min.   :0.00000   POLYGON      :202  
#  1st Qu.:0.00000   epsg:4326    :  0  
#  Median :0.00000   +proj=long...:  0  
#  Mean   :0.03449                      
#  3rd Qu.:0.00000                      
#  Max.   :5.81543                      
# Deleting layer `occ_post_PIPL_2023_Jun_med' using driver `ESRI Shapefile'
# Writing layer `occ_post_PIPL_2023_Jun_med' to data source 
# `F:/SCRAM/BRI/Integrated Movement Models/Motus Movement Models/PIPL/PIPL_2023_2018_1state_0D_var_hr24/occ_post_PIPL_2023_Jun_med_shp/occ_post_PIPL_2023_Jun_med.shp' using driver `ESRI Shapefile'
# Writing 202 features with 1 fields and geometry type Polygon.
#       Jun                   geometry  
#  Min.   :0.00000   POLYGON      :202  
#  1st Qu.:0.00000   epsg:4326    :  0  
#  Median :0.00000   +proj=long...:  0  
#  Mean   :0.05990                      
#  3rd Qu.:0.00124                      
#  Max.   :5.22259                      
# Deleting layer `occ_post_PIPL_2023_Jul_med' using driver `ESRI Shapefile'
# Writing layer `occ_post_PIPL_2023_Jul_med' to data source 
# `F:/SCRAM/BRI/Integrated Movement Models/Motus Movement Models/PIPL/PIPL_2023_2018_1state_0D_var_hr24/occ_post_PIPL_2023_Jul_med_shp/occ_post_PIPL_2023_Jul_med.shp' using driver `ESRI Shapefile'
# Writing 202 features with 1 fields and geometry type Polygon.
#       Jul                    geometry  
#  Min.   :0.000000   POLYGON      :202  
#  1st Qu.:0.000043   epsg:4326    :  0  
#  Median :0.000832   +proj=long...:  0  
#  Mean   :0.065530                      
#  3rd Qu.:0.010390                      
#  Max.   :3.414173                      
# Deleting layer `occ_post_PIPL_2023_Aug_med' using driver `ESRI Shapefile'
# Writing layer `occ_post_PIPL_2023_Aug_med' to data source 
# `F:/SCRAM/BRI/Integrated Movement Models/Motus Movement Models/PIPL/PIPL_2023_2018_1state_0D_var_hr24/occ_post_PIPL_2023_Aug_med_shp/occ_post_PIPL_2023_Aug_med.shp' using driver `ESRI Shapefile'
# Writing 202 features with 1 fields and geometry type Polygon.
#       Aug                    geometry  
#  Min.   :0.000000   POLYGON      :202  
#  1st Qu.:0.000435   epsg:4326    :  0  
#  Median :0.003935   +proj=long...:  0  
#  Mean   :0.050121                      
#  3rd Qu.:0.018391                      
#  Max.   :3.307652                      
# Writing layer `occ_post_PIPL_2023_Sep_med' to data source 
# `F:/SCRAM/BRI/Integrated Movement Models/Motus Movement Models/PIPL/PIPL_2023_2018_1state_0D_var_hr24/occ_post_PIPL_2023_Sep_med_shp/occ_post_PIPL_2023_Sep_med.shp' using driver `ESRI Shapefile'
# Writing 202 features with 1 fields and geometry type Polygon.
#       Sep               geometry  
#  Min.   : NA   POLYGON      :202  
#  1st Qu.: NA   epsg:4326    :  0  
#  Median : NA   +proj=long...:  0  
#  Mean   :NaN                      
#  3rd Qu.: NA                      
#  Max.   : NA                      
#  NA's   :202

#May Max.   :5.81543

## ROST ####
spp <- 'ROST' #hard-coded (needs manual update)!
inc <- 24 #hard-coded (needs manual update)! refers to the increment used to fit to detections and fit to regular interval increments
out.name <- paste0(spp, '_2023') #name assigned to output files, refers to year Phase 2 started
out.nam2 <- paste0('_2018_1state_0D_var_hr',inc) #hard-coded (needs manual update)! 
#2018 refers the version of species data (release year, from Loring et al. 2019)
#1state refers to the single state movement model
#0D refers to the exclusion of a drift parameter
#var refers to the conversion of the variance-covariance matrix to independent normal variance

## to load/read data for JAGS from 1-state model ####
N <- 134 #hard-coded (needs manual update)!

dir.out <- file.path(dir.proj, 'Motus Movement Models', spp, paste0(out.name, out.nam2)) #subfolder for output files

# dir.local <- file.path('P:') #server (or local data drive on high performance computer)
# dir.local <- file.path(dir.local, 'SCRAM', 'BRI', 'Integrated Movement Models', 'Motus Movement Models', spp, paste0(out.name, out.nam2)) 
# 
# load(file=file.path(dir.out, paste0('occ_post_', out.name, '_studyarea', '.RData')))
# 
# str(occ_post)
# #num [1:1000, 1:12, 1:512]
# #iterations, months, grid cells
# 
# ## load shapefile ####
# occup_BOEM_halfdeg_grid_aea_sf <- readRDS(file.path('P:', 'SCRAM', 'BRI', 'Movement', 'SCRAM Movement Models', 'occup_BOEM_halfdeg_grid_aea_sf.RDS'))
# str(occup_BOEM_halfdeg_grid_aea_sf)
# # Classes ‘sf’ and 'data.frame':	512 obs. of  18 variables:
# 
# monthLabels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",  "Aug", "Sep", "Oct", "Nov", "Dec")
# 
# dim(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][['Jan']])
# dim(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][[1]]) #same
# #512    3 1000
# #cell,mod,iteration
# 
# dim(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][[1]]) #array
# dim(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][1]) #data.frame (first values of array dimensions & geometry)
# #512   2
# 
# sum(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][[7]][,1,1]!= #same
#       st_drop_geometry(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels])[,7][,1,1], na.rm=T)
# 
# occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][7] 
# occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][,7] #should be same?
# 
# dim(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels]) #contains geometry as 13
# #512  13
# 
# dim(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][['Jan']][,1,]) #first model (Motus)
# #512 1000
# 
# dim(t(occ_post[,1,])) #Jan
# #512 1000
# 
# for (i in 1:length(monthLabels)) {
#   occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][[i]][,1,] <- t(occ_post[,i,])
# }
# 
# #calculate mean across iterations and models for months with data/predictions
# m=7 #test
# dim(apply(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][[m]], c(1,2), function(x) mean(x, na.rm = TRUE)))
# #512   3
# occup_BOEM_halfdeg_grid_aea_sf_mean <- occup_BOEM_halfdeg_grid_aea_sf
# 
# for (m in which(rowSums(apply(occ_post, c(2,3), function(x) sum(x, na.rm = TRUE)), na.rm = TRUE)>0)) {
#   occup_BOEM_halfdeg_grid_aea_sf_mean[,monthLabels][[m]][,,1] <- apply(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][[m]], c(1,2), function(x) mean(x, na.rm = TRUE))
# }
# 
# #Plot months with data/predictions (need mean otherwise plots first of 1000 iterations)
# for (m in which(rowSums(apply(occ_post, c(2,3), function(x) sum(x, na.rm = TRUE)), na.rm = TRUE)>0)) {
#   plot(occup_BOEM_halfdeg_grid_aea_sf_mean[,monthLabels][m], logz = TRUE)
# }
# #comparable to P:\SCRAM\BRI\Movement\SCRAM Movement Models\winter2023results\rost_2023_2018_1state_0D_var\occ_post_rost_2023_med.pdf
# 
# saveRDS(occup_BOEM_halfdeg_grid_aea_sf, file.path('P:', 'SCRAM', 'BRI', 'Integrated Movement Models', 'occup_BOEM_halfdeg_grid_aea_sf_rost_motus2.RDS')) 

load(file=file.path(dir.out, paste0('occ_post_', out.name, '_studyarea', '.RData')))
dir.local <- file.path('F:') #server (or local data drive on high performance computer)
dir.local <- file.path(dir.local, 'SCRAM', 'BRI', 'Integrated Movement Models', 'Motus Movement Models', spp, paste0(out.name, out.nam2)) 
dir.create(dir.local)

## Plot all months and export shapefile for study area 
BOEM.sf <- st_read('P:/SCRAM/BRI/Integrated Movement Models/data',
                   'BOEM_halfdeg_grid_latlon_2021')
occu_med <- apply(occ_post, c(3, 2), function(x) mean(x, na.rm = TRUE))

BOEM.sf <- cbind(BOEM.sf, occu_med)
BOEM.sf <- BOEM.sf[BOEM.sf$top > 36.4914 & BOEM.sf$bottom < 42.2453,] #reduce to study area extent for mapping
#str(BOEM.sf) 

monthLabels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",  "Aug", "Sep", "Oct", "Nov", "Dec")

for (m in monthLabels[6:9]) { #which(N_bymonth>0)
  dir.create(file.path(dir.local, paste0('occ_post_', out.name, '_', m, '_med_shp')))
  try(st_write(BOEM.sf[paste0(m, '')], file.path(dir.local, paste0('occ_post_', out.name, '_', m, '_med_shp'), paste0('occ_post_', out.name, '_', m, '_med.shp')), append=F))
  jpeg(filename = file.path(dir.local, paste0('occ_post_', out.name, '_', m, '_med.jpg')))
  try(plot(BOEM.sf[paste0(m, '')], logz = TRUE))
  dev.off()
  print(summary(BOEM.sf[paste0(m, '')], na.rm=T))
}

# Deleting layer `occ_post_ROST_2023_Jun_med' using driver `ESRI Shapefile'
# Writing layer `occ_post_ROST_2023_Jun_med' to data source 
# `F:/SCRAM/BRI/Integrated Movement Models/Motus Movement Models/ROST/ROST_2023_2018_1state_0D_var_hr24/occ_post_ROST_2023_Jun_med_shp/occ_post_ROST_2023_Jun_med.shp' using driver `ESRI Shapefile'
# Writing 202 features with 1 fields and geometry type Polygon.
#       Jun                   geometry  
#  Min.   :0.00000   POLYGON      :202  
#  1st Qu.:0.00000   epsg:4326    :  0  
#  Median :0.00000   +proj=long...:  0  
#  Mean   :0.05415                      
#  3rd Qu.:0.00000                      
#  Max.   :5.62099                      
# Deleting layer `occ_post_ROST_2023_Jul_med' using driver `ESRI Shapefile'
# Writing layer `occ_post_ROST_2023_Jul_med' to data source 
# `F:/SCRAM/BRI/Integrated Movement Models/Motus Movement Models/ROST/ROST_2023_2018_1state_0D_var_hr24/occ_post_ROST_2023_Jul_med_shp/occ_post_ROST_2023_Jul_med.shp' using driver `ESRI Shapefile'
# Writing 202 features with 1 fields and geometry type Polygon.
#       Jul                   geometry  
#  Min.   :0.00000   POLYGON      :202  
#  1st Qu.:0.00000   epsg:4326    :  0  
#  Median :0.00000   +proj=long...:  0  
#  Mean   :0.08522                      
#  3rd Qu.:0.00000                      
#  Max.   :9.14283                      
# Deleting layer `occ_post_ROST_2023_Aug_med' using driver `ESRI Shapefile'
# Writing layer `occ_post_ROST_2023_Aug_med' to data source 
# `F:/SCRAM/BRI/Integrated Movement Models/Motus Movement Models/ROST/ROST_2023_2018_1state_0D_var_hr24/occ_post_ROST_2023_Aug_med_shp/occ_post_ROST_2023_Aug_med.shp' using driver `ESRI Shapefile'
# Writing 202 features with 1 fields and geometry type Polygon.
#       Aug                   geometry  
#  Min.   :0.00000   POLYGON      :202  
#  1st Qu.:0.00000   epsg:4326    :  0  
#  Median :0.00000   +proj=long...:  0  
#  Mean   :0.08564                      
#  3rd Qu.:0.00000                      
#  Max.   :7.22129                      
# Deleting layer `occ_post_ROST_2023_Sep_med' using driver `ESRI Shapefile'
# Writing layer `occ_post_ROST_2023_Sep_med' to data source 
# `F:/SCRAM/BRI/Integrated Movement Models/Motus Movement Models/ROST/ROST_2023_2018_1state_0D_var_hr24/occ_post_ROST_2023_Sep_med_shp/occ_post_ROST_2023_Sep_med.shp' using driver `ESRI Shapefile'
# Writing 202 features with 1 fields and geometry type Polygon.
#       Sep                   geometry  
#  Min.   :0.00000   POLYGON      :202  
#  1st Qu.:0.00000   epsg:4326    :  0  
#  Median :0.00000   +proj=long...:  0  
#  Mean   :0.06538                      
#  3rd Qu.:0.00000                      
#  Max.   :5.53556

#Jul Max.   :9.14283

#Create study area polygon
BOEM.sf <- st_read('P:/SCRAM/BRI/Integrated Movement Models/data',
                   'BOEM_halfdeg_grid_latlon_2021')
BOEM.sf <- BOEM.sf[BOEM.sf$top > 36.4914 & BOEM.sf$bottom < 42.2453,] #keep bounding grid cells that contain those latitudes

BOEM.sf$dissolve <- 1
BOEM.sf.outline <-aggregate(
  x=BOEM.sf,
  by=list(BOEM.sf$dissolve),
  FUN=length,
  do_union = TRUE,
  simplify = TRUE,
  join = st_intersects
)

BOEM.sf.outline <- st_geometry(BOEM.sf.outline)

str(BOEM.sf.outline)
plot(BOEM.sf.outline)
dir.create(file.path('P:', 'SCRAM', 'BRI', 'Integrated Movement Models', 'Motus Movement Models', 'MotusStudyAreaShp'))
st_write(BOEM.sf.outline, file.path('P:', 'SCRAM', 'BRI', 'Integrated Movement Models', 'Motus Movement Models', 'MotusStudyAreaShp', 'MotusStudyArea.shp'), append=F)
