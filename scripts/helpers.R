library(shiny)
library(shinyBS)
library(shinydashboard)
library(shinyFiles)
library(shinybusy)
library(dplyr)
library(tidyverse)
library(magrittr)
library(msm)
library(shinyjs)
library(shinyWidgets)
library(data.table)
library(DT)
library(zip)
library(RColorBrewer)
library(pracma)
library(rcmdcheck)
library(sp)
library(sf)
library(rgdal)
library(geosphere)
library(readr)
library(rgeos)
library(htmltools)
library(rmarkdown)
library(fs)
library(latexpdf)
library(Hmisc)
library(tinytex)
library(gmailr)
library(ipc)
library(leaflet)
library(leaflet.esri)
library(ggplot2)
library(cowplot)
library(basemaps)
library(knitr)
library(kableExtra)
# library(future)
# library(promises)
# plan(multisession)


#load the species movement model mean monthly probability data
#generate all species data
data_dir <- "data"

#add disclaimer about potentially higher collision risk estimate due to inclusion of land 
coastal_disclaimer <- "Cells that overlap land have higher collision estimates,\nas passage rate estimates include birds that are on land as well as overwater."

# for (species in c("Red_Knot", "Piping_Plover", "Roseate_Tern")){ #, "Common_Tern")){
#   # "Red_Knot_monthly_prob_BOEM_half_deg_trunc.RData"
#   data_layer <-  paste0(species, "_monthly_prob_BOEM_half_deg_trunc")
#   load(file.path(data_dir, paste0(species, "_monthly_prob_BOEM_half_deg_trunc.RData")))
#   assign(data_layer, spp_monthly_prob_BOEM_half_deg)
#   data_layer2 <-  paste0(species, "_monthly_prob_BOEM_half_deg_last")
#   load(file.path(data_dir, paste0(species, "_monthly_prob_BOEM_half_deg_last.RData")))
#   assign(data_layer2, spp_monthly_prob_BOEM_half_deg)
# }
# rm(spp_monthly_prob_BOEM_half_deg)
BOEM_lease_outlines <- sf::read_sf("data/BOEMWindLeaseOutlines_6_1_2022.shp") %>% st_transform(3857)
BOEM_planning_area_outlines <- sf::read_sf("data/BOEMWindPlanningAreas_06_01_2022.shp") %>% st_transform(3857)
states_sf <- sf::read_sf("data/statesp020.shp") %>% st_transform(3857)

# https://stackoverflow.com/questions/6177629/how-to-silence-the-output-from-this-r-package
# small function that silences cat() and print() (but not message() or warning()) and returns whatever the expression returned:
shut_up = function(expr) {
  #temp file
  f = file()
  
  #write output to that file
  sink(file = f)
  
  #evaluate expr in original environment
  y = eval(expr, envir = parent.frame())
  
  #close sink
  sink()
  
  #get rid of file
  close(f)
  
  # y
}

