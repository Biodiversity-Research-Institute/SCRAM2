## Script to create shapefile of species output from crawl model (SCRAM 2) for REKN

##Author: Holly Goyert
##Date: 06/05/2024

# #R v4.3.1

## load required packages
library(R2jags)
library(R2WinBUGS)
library(tidyverse)
library(sf)
library(raster)
library(rgeos)
library(MCMCvis) #MCMCsummary rather than summary
#try(library(motus)) #to install: browseURL("https://motuswts.github.io/motus/articles/02-installing-packages.html")

## set directories
#dir.home <- file.path('P:') #Home directory (network)
#dir.root <- file.path(dir.home, 'SCRAM', 'BRI') #Root directory (network)

usr.name <- 'user.name' #hard-coded (needs manual update)!
dir.home <- file.path('C:', 'Users', usr.name, 'OneDrive - Biodiversity Research Institute', 'Documents') #Home directory (local)
dir.root <- file.path(dir.home, 'ProjectsOffshoreWind', 'Atlantic', 'USFWS', 'SCRAM', 'BRI') #Root directory (local)

dir.proj <- file.path(dir.root, 'Integrated Movement Models') #Project directory

spp <- 'REKN' #hard-coded (needs manual update)!
dir.data <- file.path(dir.proj, 'data', spp, 'GPS') #data directory

dir.work <- file.path(dir.proj, 'GPS Movement Models')
setwd(dir.work) # setwd <- dir.work 

dir.local <- file.path(dir.work, spp)
#dir.local <- file.path('P:', 'SCRAM', 'BRI', 'Integrated Movement Models', 'GPS Movement Models', spp)

## load shapefile ####
occup_BOEM_halfdeg_grid_aea_sf <- readRDS(file.path('P:', 'SCRAM', 'BRI', 'Movement', 'SCRAM Movement Models', 'occup_BOEM_halfdeg_grid_aea_sf_motus_gps1.RDS'))
str(occup_BOEM_halfdeg_grid_aea_sf)
# Classes ‘sf’ and 'data.frame':	512 obs. of  18 variables:
#   $ id                : chr  "1" "2" "3" "4" ...
# $ lon               : num  -66.7 -67.2 -67.7 -68.2 -68.2 ...
# $ lat               : num  44.5 44.5 44.5 44.5 44 ...
# $ geometry          :sfc_POLYGON of length 512; first list element: List of 1
# ..$ : num [1:5, 1:2] 2142950 2125153 2160613 2178708 2142950 ...
# ..- attr(*, "class")= chr [1:3] "XY" "POLYGON" "sfg"
# $ area_sqkm         : num  2209 2209 2209 2209 2227 ...
# $ mean_cell_width_km: num  47 47 47 47 47.2 ...
# $ Jan               : num [1:512, 1:3, 1:1000] NA NA NA NA NA NA NA NA NA NA ...
# $ Feb               : num [1:512, 1:3, 1:1000] NA NA NA NA NA NA NA NA NA NA ...
# $ Mar               : num [1:512, 1:3, 1:1000] NA NA NA NA NA NA NA NA NA NA ...
# $ Apr               : num [1:512, 1:3, 1:1000] NA NA NA NA NA NA NA NA NA NA ...
# $ May               : num [1:512, 1:3, 1:1000] NA NA NA NA NA NA NA NA NA NA ...
# $ Jun               : num [1:512, 1:3, 1:1000] NA NA NA NA NA NA NA NA NA NA ...
# $ Jul               : num [1:512, 1:3, 1:1000] NA NA NA NA NA NA NA NA NA NA ...
# $ Aug               : num [1:512, 1:3, 1:1000] NA NA NA NA NA NA NA NA NA NA ...
# $ Sep               : num [1:512, 1:3, 1:1000] NA NA NA NA NA NA NA NA NA NA ...
# $ Oct               : num [1:512, 1:3, 1:1000] NA NA NA NA NA NA NA NA NA NA ...
# $ Nov               : num [1:512, 1:3, 1:1000] NA NA NA NA NA NA NA NA NA NA ...
# $ Dec               : num [1:512, 1:3, 1:1000] NA NA NA NA NA NA NA NA NA NA ...
# - attr(*, "sf_column")= chr "geometry"
# - attr(*, "agr")= Factor w/ 3 levels "constant","aggregate",..: NA NA NA NA NA NA NA NA NA NA ...
# ..- attr(*, "names")= chr [1:17] "id" "lon" "lat" "area_sqkm" ...

monthLabels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",  "Aug", "Sep", "Oct", "Nov", "Dec")

dim(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][['Jan']])
dim(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][[1]]) #same
#512    3 1000
#cell,mod,iteration

dim(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][[1]]) #array
dim(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][1]) #data.frame (first values of array dimensions & geometry)
#512   2

sum(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][[7]][,1,1]!= #same
st_drop_geometry(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels])[,7][,1,1], na.rm=T)

occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][7]
occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][,7] #should be same?

dim(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels]) #contains geometry as 13
#512  13

dim(occup_BOEM_halfdeg_grid_aea_sf[,monthLabels][['Jan']][,2,]) #second model (GPS)
#512 1000

move_model <- array(NA, c(512, 12, 1000))
#copy crawl model to first of three models
for (m in 1:12){
  move_model[,m,] <- occup_BOEM_halfdeg_grid_aea_sf[[monthLabels[m]]][,2,]
  }

colnames(move_model) <- monthLabels

BOEM.sf <- st_read('P:/SCRAM/BRI/Integrated Movement Models/data',
                   'BOEM_halfdeg_grid_latlon_2021')
BOEM.sf$occu_med <- apply(move_model, c(1, 2), function(x) mean(x, na.rm = TRUE))
occu_med <- apply(move_model, c(1, 2), function(x) mean(x, na.rm = TRUE))
BOEM.sf <- cbind(BOEM.sf, occu_med)
str(BOEM.sf)

BOEM.sf$occu_med_mo <- apply(BOEM.sf$occu_med[,apply(BOEM.sf$occu_med, 2, sum, na.rm=TRUE)>0], 1, mean, na.rm=TRUE)
summary(BOEM.sf$occu_med_mo)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.000115 0.000503 0.043772 0.001398 7.532798
dir.create(file.path(dir.local, paste0('gps_crawl_fig')))
dir.create(file.path(dir.local, paste0('gps_crawl_fig'), paste0('gps_crawl_', 'med_shp')))
jpeg(filename = file.path(dir.local, paste0('gps_crawl_fig'), paste0('gps_crawl_', 'med_shp'), paste0('gps_crawl_', 'med.jpg')))
plot(na.omit(BOEM.sf["occu_med_mo"]), logz = TRUE)
dev.off()

st_write(BOEM.sf["occu_med_mo"], file.path(dir.local, paste0('gps_crawl_fig'), paste0('gps_crawl_', 'med_shp'), paste0('gps_crawl_', 'med.shp')), append=F)

#Plot months with data/predictions
for (month in monthLabels[6:12]) { 
  dir.create(file.path(dir.local, paste0('gps_crawl_fig'), paste0('gps_crawl_', month, '_med_shp')))
  try(st_write(BOEM.sf[paste0(month)], file.path(dir.local, paste0('gps_crawl_fig'), paste0('gps_crawl_', month, '_med_shp'), paste0('gps_crawl_', month, '_med.shp')), append=F))

  jpeg(filename = file.path(dir.local, paste0('gps_crawl_fig'), paste0('gps_crawl_', month, '_med_shp'), paste0('gps_crawl_', month, '_med.jpg')))
  try(plot(BOEM.sf[paste0(month)], logz = TRUE))
  dev.off()
  print(summary(BOEM.sf[paste0(month)], na.rm=T))
}

# Deleting layer `gps_crawl_Jun_med' using driver `ESRI Shapefile'
# Writing layer `gps_crawl_Jun_med' to data source 
# `C:/Users/holly.goyert/OneDrive - Biodiversity Research Institute/Documents/ProjectsOffshoreWind/Atlantic/USFWS/SCRAM/BRI/Integrated Movement Models/GPS Movement Models/REKN/gps_crawl_fig/gps_crawl_Jun_med_shp/gps_crawl_Jun_med.shp' using driver `ESRI Shapefile'
# Writing 512 features with 1 fields and geometry type Polygon.
#       Jun               geometry  
#  Min.   : NA   POLYGON      :512  
#  1st Qu.: NA   epsg:4326    :  0  
#  Median : NA   +proj=long...:  0  
#  Mean   :NaN                      
#  3rd Qu.: NA                      
#  Max.   : NA                      
#  NA's   :512                      
# Deleting layer `gps_crawl_Jul_med' using driver `ESRI Shapefile'
# Writing layer `gps_crawl_Jul_med' to data source 
# `C:/Users/holly.goyert/OneDrive - Biodiversity Research Institute/Documents/ProjectsOffshoreWind/Atlantic/USFWS/SCRAM/BRI/Integrated Movement Models/GPS Movement Models/REKN/gps_crawl_fig/gps_crawl_Jul_med_shp/gps_crawl_Jul_med.shp' using driver `ESRI Shapefile'
# Writing 512 features with 1 fields and geometry type Polygon.
#       Jul               geometry  
#  Min.   : NA   POLYGON      :512  
#  1st Qu.: NA   epsg:4326    :  0  
#  Median : NA   +proj=long...:  0  
#  Mean   :NaN                      
#  3rd Qu.: NA                      
#  Max.   : NA                      
#  NA's   :512                      
# Deleting layer `gps_crawl_Aug_med' using driver `ESRI Shapefile'
# Writing layer `gps_crawl_Aug_med' to data source 
# `C:/Users/holly.goyert/OneDrive - Biodiversity Research Institute/Documents/ProjectsOffshoreWind/Atlantic/USFWS/SCRAM/BRI/Integrated Movement Models/GPS Movement Models/REKN/gps_crawl_fig/gps_crawl_Aug_med_shp/gps_crawl_Aug_med.shp' using driver `ESRI Shapefile'
# Writing 512 features with 1 fields and geometry type Polygon.
#       Aug                    geometry  
#  Min.   :0.000000   POLYGON      :512  
#  1st Qu.:0.000000   epsg:4326    :  0  
#  Median :0.000392   +proj=long...:  0  
#  Mean   :0.026497                      
#  3rd Qu.:0.002023                      
#  Max.   :7.651327                      
# Deleting layer `gps_crawl_Sep_med' using driver `ESRI Shapefile'
# Writing layer `gps_crawl_Sep_med' to data source 
# `C:/Users/holly.goyert/OneDrive - Biodiversity Research Institute/Documents/ProjectsOffshoreWind/Atlantic/USFWS/SCRAM/BRI/Integrated Movement Models/GPS Movement Models/REKN/gps_crawl_fig/gps_crawl_Sep_med_shp/gps_crawl_Sep_med.shp' using driver `ESRI Shapefile'
# Writing 512 features with 1 fields and geometry type Polygon.
#       Sep                    geometry  
#  Min.   :0.000000   POLYGON      :512  
#  1st Qu.:0.000012   epsg:4326    :  0  
#  Median :0.000426   +proj=long...:  0  
#  Mean   :0.035677                      
#  3rd Qu.:0.001138                      
#  Max.   :5.301300                      
# Deleting layer `gps_crawl_Oct_med' using driver `ESRI Shapefile'
# Writing layer `gps_crawl_Oct_med' to data source 
# `C:/Users/holly.goyert/OneDrive - Biodiversity Research Institute/Documents/ProjectsOffshoreWind/Atlantic/USFWS/SCRAM/BRI/Integrated Movement Models/GPS Movement Models/REKN/gps_crawl_fig/gps_crawl_Oct_med_shp/gps_crawl_Oct_med.shp' using driver `ESRI Shapefile'
# Writing 512 features with 1 fields and geometry type Polygon.
#       Oct                     geometry  
#  Min.   : 0.000000   POLYGON      :512  
#  1st Qu.: 0.000060   epsg:4326    :  0  
#  Median : 0.000303   +proj=long...:  0  
#  Mean   : 0.058902                      
#  3rd Qu.: 0.001326                      
#  Max.   :15.023101                      
# Deleting layer `gps_crawl_Nov_med' using driver `ESRI Shapefile'
# Writing layer `gps_crawl_Nov_med' to data source 
# `C:/Users/holly.goyert/OneDrive - Biodiversity Research Institute/Documents/ProjectsOffshoreWind/Atlantic/USFWS/SCRAM/BRI/Integrated Movement Models/GPS Movement Models/REKN/gps_crawl_fig/gps_crawl_Nov_med_shp/gps_crawl_Nov_med.shp' using driver `ESRI Shapefile'
# Writing 512 features with 1 fields and geometry type Polygon.
#       Nov                    geometry  
#  Min.   : 0.00000   POLYGON      :512  
#  1st Qu.: 0.00000   epsg:4326    :  0  
#  Median : 0.00000   +proj=long...:  0  
#  Mean   : 0.05401                      
#  3rd Qu.: 0.00000                      
#  Max.   :10.88481                      
# Deleting layer `gps_crawl_Dec_med' using driver `ESRI Shapefile'
# Writing layer `gps_crawl_Dec_med' to data source 
# `C:/Users/holly.goyert/OneDrive - Biodiversity Research Institute/Documents/ProjectsOffshoreWind/Atlantic/USFWS/SCRAM/BRI/Integrated Movement Models/GPS Movement Models/REKN/gps_crawl_fig/gps_crawl_Dec_med_shp/gps_crawl_Dec_med.shp' using driver `ESRI Shapefile'
# Writing 512 features with 1 fields and geometry type Polygon.
#       Dec               geometry  
#  Min.   : NA   POLYGON      :512  
#  1st Qu.: NA   epsg:4326    :  0  
#  Median : NA   +proj=long...:  0  
#  Mean   :NaN                      
#  3rd Qu.: NA                      
#  Max.   : NA                      
#  NA's   :512

# Oct Max.   :15.023101
