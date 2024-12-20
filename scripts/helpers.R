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
library(geosphere)
library(readr)
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
library(USA.state.boundaries)
library(ggplot2)
library(cowplot)

# set.seed(11)
options(scipen=100, digits=6)
#load the species movement model mean monthly probability data
#generate all species data
data_dir <- "data"

#add disclaimer about potentially higher collision risk estimate due to inclusion of land 
coastal_disclaimer <- "Cells that overlap land have higher collision estimates,\nas passage rate estimates include birds that are on land as well as overwater."

#ATG - 010423 - updated lease area and planning area outlines from BOEM - 09 Oct 2024
BOEM_lease_outlines <- sf::read_sf("data/BOEM_Wind_Lease_Outlines_06_06_2024.shp") %>% st_transform(3857)
BOEM_planning_area_outlines <- sf::read_sf("data/BOEM_Wind_Planning_Area_Outlines_04_29_2024.shp") %>% st_transform(3857)
states_sf <- sf::read_sf("data/statesp020.shp") %>% st_transform(3857)

# buffer the state layer for intersection with corridors because the coastal outline for migration corridors is very crude and 
# sometimes lands in the water with no state assignment.
state_boundaries_ea <- st_transform(state_boundaries_wgs84, 9822) %>%  #from Package USA.state.boundaries 
  st_buffer(dist = 4000)

MotusStudyArea_sf <- read_sf("data/MotusStudyArea.shp") %>% st_transform(4326)
MotusStudyArea_sf$studyarea = "Motus study area"

MotusStudyArea_lines_sf <- st_cast(MotusStudyArea_sf, "LINESTRING")

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

#function to remove UI reference from server
#https://www.r-bloggers.com/2020/02/shiny-add-removing-modules-dynamically/
remove_shiny_inputs <- function(id, .input) {
  invisible(
    lapply(grep(id, names(.input), value = TRUE), function(i) {
      .subset2(.input, "impl")$.values$remove(i)
    })
  )
}

#function to create an east-west line of known lenght and center point a specified distance
#center_pt_lat_long is the center of the line lat long values vector
create_EW_line <- function(center_pt_lat_long, length_km){
  west_pt <- geosphere::destPointRhumb(center_pt_lat_long, b = -90, d = length_km * 1000)
  east_pt <- geosphere::destPointRhumb(center_pt_lat_long, b = 90, d = length_km * 1000)
  coords <- cbind(west_pt, east_pt)
  EW_line <- st_sfc(
    lapply(1:nrow(coords),
           function(i){
             sf::st_linestring(matrix(coords[i,],ncol=2,byrow=TRUE))
           }))
  
  EW_line <- sf::st_sf(EW_line, crs = sf::st_crs(4326))
  return(EW_line)
}

create_west_line <- function(center_pt_lat_long, length_km){
  west_pt <- geosphere::destPointRhumb(center_pt_lat_long, b = -90, d = length_km * 1000)
  coords <- cbind(west_pt, center_pt_lat_long)
  west_line <- st_sfc(
    lapply(1:nrow(coords),
           function(i){
             sf::st_linestring(matrix(coords[i,],ncol=2,byrow=TRUE))
           }))
  
  west_line <- sf::st_sf(west_line, crs = sf::st_crs(4326))
  return(west_line)
}

#function to perform linear interpolation between two variables
interp <- function(a, b, f){
  return(a * (1 - f) + (b * f))
}


weighted_mean <- function(x, y, wt_x) case_when(
  is.na(x) & is.na(y) ~ NA,
  is.na(x) & !is.na(y) ~ y,
  !is.na(x) & is.na(y) ~ x,
  !is.na(x) & !is.na(y) ~ x * wt_x + y * (1 - wt_x)
)

gen_weights <- function(x, y, wt_x) case_when(
  is.na(sum(x)) & is.na(sum(y)) ~ NA,
  is.na(sum(x)) & !is.na(sum(y)) ~ NA,
  !is.na(sum(x)) & is.na(sum(y)) ~ 1,
  !is.na(sum(x)) & !is.na(sum(y)) ~ 1 - wt_x
)


SC_states_north = c("South Carolina","North Carolina","Virginia","Maryland","Delaware","New Jersey","New York",
                    "Connecticut","Rhode Island","Massachusetts","New Hampshire","Maine","Eastern Canada")
NC_states_north = c("North Carolina","Virginia","Maryland","Delaware","New Jersey","New York","Connecticut", 
                    "Rhode Island", "Massachusetts", "New Hampshire", "Maine","Eastern Canada")
VA_states_north = c("Virginia","Maryland","Delaware","New Jersey","New York","Connecticut", "Rhode Island", 
                    "Massachusetts","New Hampshire","Maine","Eastern Canada")
MD_states_north = c("Maryland","Delaware","New Jersey","New York","Connecticut","Rhode Island","Massachusetts",
                    "New Hampshire","Maine","Eastern Canada")
DE_states_north = c("Delaware","New Jersey","New York","Connecticut","Rhode Island","Massachusetts","New Hampshire",
                    "Maine","Eastern Canada")
NJ_states_north = c("New Jersey","New York","Connecticut","Rhode Island","Massachusetts","New Hampshire","Maine","Eastern Canada")
NY_states_north = c("New York","Connecticut","Rhode Island","Massachusetts","New Hampshire","Maine","Eastern Canada")
CT_states_north = c("Connecticut","Rhode Island","Massachusetts","New Hampshire","Maine","Eastern Canada")
RI_states_north = c("Rhode Island","Massachusetts","New Hampshire","Maine","Eastern Canada")
MA_states_north = c("Massachusetts","New Hampshire","Maine","Eastern Canada")
NH_states_north = c("New Hampshire","Maine","Eastern Canada")
ME_states_north = c("Maine","Eastern Canada")
