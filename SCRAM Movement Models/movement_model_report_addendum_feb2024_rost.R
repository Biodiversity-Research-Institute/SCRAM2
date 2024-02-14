##################################################################################################################################
### Script to 1) prepare Motus data files for analyses, 2) run movement model in JAGS, 
### 3) use posterior estimates to obtain occupancy estimates for user-specified spatial units, and 
### 4) create simple plots showing mean and variance of results over space. 

### Depending on which components are run, requires: 'false_pos_filter.R', 'remove_false_pos.R', 
### 'remove_dup_bursts.R', 'serial_time.R', 'Pr_occupancy_season_trunc.R'
##################################################################################################################################

###Authors: Holly Goyert, Evan Adams (BRI), Chris Field (URI)
###Date: 02/01/2024

# ## SCRAM 2.0.0 movement models (Goyert et al. 2024), updated from:
# ## SCRAM 1.0.3 (Adams et al. 2022)
# ## Based on hDCRWS from Baldwin et al. 2018 and Jonsen et al. 2016 (bsam)
# ## hierarchical (to estimate parameters jointly across multiple individual tracking datasets)
# ## Discrete-time, continuous-space, Correlated Random Walk (DCRW for location filtering)
# ## with a Spatial measurement-error observation model (DCRWS for location filtering and behavioural state estimation)
# ## - includes random deviate from hDCRWS to force gamma[1] > gamma[2] since b=1 is migratory with high correlation between subsequent movements and b=2 is staging with low correlation
# ## - excludes Drift and relaxes assumption of covariance among error terms to independent normal variance for observed and estimated locations

# ## - predicts all species only to study area containing active stations during study period (as USFWS recommended)

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

dir.proj <- file.path(dir.root, 'Movement', 'SCRAM Movement Models') #Project directory
dir.data <- file.path(dir.proj, 'data') #GIS data directory

dir.work <- dir.proj
setwd(dir.work) # setwd <- dir.work #updated from "P:/SCRAM/BRI/Movement/SCRAM Movement Models"

out.name <- 'rost_2023' #name assigned to output files
out.nam2 <- '_2018_1state_0D_var' #hard-coded (needs manual update)! 
#2018 refers the version of species data (release year, from Loring et al. 2019)
#1state refers to the single state movement model
#0D refers to the exclusion of a drift parameter
#var refers to the conversion of the variance-covariance matrix to independent normal variance


dir.create(file.path(dir.proj, paste0(out.name, out.nam2))) #subfolder for output files
dir.out <- file.path(dir.proj, paste0(out.name, out.nam2)) #file.path(dir.proj, 'winter2023results', paste0(out.name, out.nam2))

#dir.local <- file.path('D:') #local data drive on high performance computer / server
dir.local <- file.path('F:') #local data drive on high performance computer / server
dir.local <- file.path(dir.local, 'SCRAM', 'BRI', 'Movement', 'SCRAM Movement Models') #create manually otherwise need to specify subfolder by subfolder
dir.create(file.path(dir.local, paste0(out.name, out.nam2)))

## load/read (and combine) Motus data for species of interest in movement analyses ######################
proj_14_2015_rost_final <- read.csv(file.path(dir.data,'Appendix_H-Motus_Detection_Data_ROST_2015_v20190708.csv'), header=TRUE, stringsAsFactors = FALSE)
proj_14_2016_rost_final <- read.csv(file.path(dir.data,'Appendix_H-Motus_Detection_Data_ROST_2016_v20190708.csv'), header=TRUE, stringsAsFactors = FALSE)
proj_14_2017_rost_final <- read.csv(file.path(dir.data,'Appendix_H-Motus_Detection_Data_ROST_2017_v20190708.csv'), header=TRUE, stringsAsFactors = FALSE)
allindvs <- rbind(proj_14_2015_rost_final, proj_14_2016_rost_final, proj_14_2017_rost_final)

## remove missing coordinates ####			 
allindvs <- allindvs[!is.na(allindvs$lat)&!is.na(allindvs$lon), ]
str(allindvs)

##test subset of data
# set.seed(1234)
# (test.ids <- sample(unique(allindvs$id), size=3))
# allindvs <- allindvs[which((allindvs$id %in% test.ids)), ]

length(unique(allindvs$id))

## create columns for time variables for movement model ##########################################################################
# create a new column for hour of detection with minutes represented by a decimal
hours_dec <- round(as.numeric(substr(allindvs$ts_gmt, 12, 13)) + 
                     (as.numeric(substr(allindvs$ts_gmt, 15, 16))/60), 2)
allindvs <- cbind(allindvs, hours_dec)

# create a new column for month
month <- as.numeric(substr(allindvs$ts_gmt, 6, 7))
allindvs <- cbind(allindvs, month)

# create a new column for day
day <- as.numeric(substr(allindvs$ts_gmt, 9, 10))
allindvs <- cbind(allindvs, day)

# create a new column for year
year <- as.numeric(substr(allindvs$ts_gmt, 1, 4))
allindvs<- cbind(allindvs, year)

# order the detections by individual, year, day, hour, and minute
allindvs_ordered <- allindvs[order(allindvs$id, allindvs$year, substr(allindvs$ts_gmt, 6, 11)), ] #includes month-day only, though same as using hh-mm-ss

### FILTERING ####
## pass data through filters for false positives, duplicate bursts, and stationary individuals ###################################
# one option for filtering motus data is 'filterByActivity.R' from the 'motus' package
# filterByActivity() from 'motus' requires the data in an SQL table
# using filtering scripts that run on flat data frames instead until Motus releases pre-filtered downloadable data

# add column for burst length
source('false_pos_filter.R')
### false_pos_filter() identifies "bursts" as the number of detections for an individual at the same site, on the same day, in the same hour, within 1 minute of each other
### The function takes a data frame of detections and a vector specifying an index for each row of the data frame (e.g. 1: length(x[,1]))
### This function can take a long time to run for data frames with greater than ~ 50000 rows

#trying to figure out how this is calculated
# bursts <- unique(allindvs_ordered[, c('id', 'month', 'day')])
# bursts$key <- paste0(bursts[1, 1:3])
# bursts <- allindvs_ordered %>% 
#   group_by(id, month, day) %>% 
#   summarize(N = length(id))
# allindvs_ordered$burst_length <- rep(bursts$N, each = bursts$N)

#old code that we don't have the fn for
burst_length <- false_pos_filter(allindvs_ordered)
allindvs_ordered <- cbind(allindvs_ordered, burst_length)

#I'll use runLen as burst_length for the moment
# burst_length <- allindvs_ordered$runLen
# allindvs_ordered <- cbind(allindvs_ordered, burst_length)

# remove false positives 
source('remove_false_pos.R') 
### remove_false_pos() removes rows for which the value for the number of detections within a burst
### is less than or equal to the value specified for y
### This function works with the results of false_pos_filter(), x
# allindvs_ordered_filtered <- remove_false_pos(allindvs_ordered, 3)
allindvs_ordered_filtered <- remove_false_pos(allindvs_ordered, 3)

# remove duplicate detections by defining bursts as all detections within a 24-hour period
source('remove_dup_bursts.R')
### remove_dup_bursts() removes all duplicate detections within a burst, preserving only the first detection. 
### The burst can be defined as detections that occur within the same 24 hour period or hour.
### x is a data frame of detections and y is a string specifying whether bursts should be defined by 'day' or 'hour'
allindvs_ordered_filtered_nodups <- remove_dup_bursts(allindvs_ordered_filtered, "day")

## create variables for indexing the JAGS movement model (variable terminology as in Baldwin et al. 2018) #############################################
# create a vector (Sind_obs) that references, for each individual, the position of the vector (for the observed data) that denotes its first detection
id_index <- mat.or.vec(length(allindvs_ordered_filtered_nodups[,1]), 1)
allindvs_ordered_filtered_nodups <- cbind(allindvs_ordered_filtered_nodups, id_index)
  # create an index for individual
uni_inds <- unique(allindvs_ordered_filtered_nodups$id)
Sind_obs <- unique(allindvs_ordered_filtered_nodups$id)
for(i in 1:length(uni_inds)){
  allindvs_ordered_filtered_nodups[allindvs_ordered_filtered_nodups$id==uni_inds[i], 'id_index'] <- i
  Sind_obs[i] <- min(which(allindvs_ordered_filtered_nodups$id==uni_inds[i]))
}

source('serial_time.R')
### serial_time() adds a column to the data frame for the number of days or hours since the first record in the data frame
### x is the data frame to which to add the column, and y specifies 'days' or 'hours'
allindvs_ordered_filtered_nodups <- serial_time(allindvs_ordered_filtered_nodups, 'day')

## Prepare pre-filtered data structure for JAGS model ####
str(allindvs_ordered_filtered_nodups)

season_length <- max(allindvs_ordered_filtered_nodups$days_since)

# number of individuals
(N <- length(uni_inds))

# create a vector (Xidx) that references, for each individual, the position of the vector (for the regular time steps) that denotes its first detection
Xidx <- mat.or.vec(N, 1)
# Xidx2 indexes the last day of the season, with respect to x, for each individual
Xidx2 <- mat.or.vec(N, 1)
# Xidx3 is used for posterior checks of movement model 
Xidx3 <- mat.or.vec(N, 1) #the day of the season that corresponds to each individual's last observation
Zidx <- mat.or.vec(N, 1)  #the day of the season that corresponds to each individual's first observation
for(i in 1:N){
  Xidx[i] <- min(allindvs_ordered_filtered_nodups[allindvs_ordered_filtered_nodups$id==uni_inds[i], 'days_since']) + season_length*(i-1)
  Xidx2[i] <- season_length*i
  Xidx3[i] <- max(allindvs_ordered_filtered_nodups[allindvs_ordered_filtered_nodups$id==uni_inds[i], 'days_since'])
  Zidx[i] <- min(allindvs_ordered_filtered_nodups[allindvs_ordered_filtered_nodups$id==uni_inds[i], 'days_since'])
}

# create the idx file, which indexes for each observation, which "regular" day it is related to, 
# taking into account the fact that the data for all individuals is specified as one vector
a <- 1
for(i in 1:N){
  # season length * number of individuals already recorded is added to the days since the start of the season
  a <- c(a, (allindvs_ordered_filtered_nodups[allindvs_ordered_filtered_nodups$id==uni_inds[i], 'days_since']+((season_length)*(i-1))))
}
idx <- a[-1]

# locations


y <- allindvs_ordered_filtered_nodups[, c("lat", "lon")] #note lat first coordinate, then lon
# there is an additional position added to the end of the vector to ensure that the index loops stay in bounds (as in Baldwin et al. 2018)
Xidx <- c(Xidx, (season_length*N+1))
Xidx2 <- c(Xidx2, (season_length*N+1))
# there is an additional position added to the end of the vector to ensure that the index loops stay in bounds (as in Baldwin et al. 2018)
Yidx <- c(Sind_obs, (length(y[,1])+1))[1:(N+1)]
Y <- season_length

# priors for x and b, specifying NAs for nodes that are not observed
x1 <- rep(NA, season_length*N)
x2 <- rep(NA, season_length*N)
x <- cbind(x1, x2)
b <- rep(NA, season_length*N)
for(k in 1:N){
  x[(Xidx[k]+1):Xidx2[k], 1] <- mean(y[(Yidx[k]+1):(Yidx[k+1]-1),1])
  x[(Xidx[k]+1):Xidx2[k], 2] <- mean(y[(Yidx[k]+1):(Yidx[k+1]-1),2])
  b[(Xidx[k]+1):(Xidx2[k]-1)] <- 1
}

lat_min <- min(y[,1]) - sd(y[,1])/3
lat_max <- max(y[,1]) + sd(y[,1])/3
lon_min <- min(y[,2]) - sd(y[,2])/3
lon_max <- max(y[,2]) + sd(y[,2])/3

# number of iterations for the JAGS model
# n.iter <- 200000
# n.burnin <- 100000
n.iter <- 30000 #100 #20000
n.burnin <- 30000 #10 #10000
nc <- 3 #1
n.thin <- 10 #5 #1 #reduce memory output		   

save.image(file.path(dir.out, paste0('pre_proc_', out.name, '.RData'))) #save image of data processed prior to running model
#load(file.path(dir.out, paste0('pre_proc_', out.name, '.RData')))

## run JAGS model ################################################################################################################
# setwd("~/")
# 2 state

##################################################################################################################################
### JAGS model and associated posterior checks with plots
##################################################################################################################################

# general one-state model structure modified from Baldwin et al. 2018
MOVE <- function(){
  ## prior on measurement error model: tag variance (measurement-error scalar) due to noise or attachment
  sd ~ dunif(0, 0.1) # sd ~ dunif(0, 0.001) #EA: constant gives very low chance of detections longer range than 50km with the vast majority with 22km of the tower
  # sd <- 0.1 # (converted from decimal degrees using the average size of a lat/lon grid cell). 
  tau <- 1/(sd*sd) #observed measurement error (tag noise) is the same for latitude and longitude and treated the same across animals

  # priors on process model
  # phi ~ dunif(0, 1) #2-state; prob of switching between states

  # priors on process uncertainty (multivariate normally distributed variability in movement around unobserved location state Xt): variance-covariance matrix
  # sigma_c[1, 1] ~ dexp(1) # it cannot be negative
  # sigma_c[1, 2] <- 0
  # sigma_c[2, 1] <- 0
  # # sigma_c[2, 2, 1] ~ dexp(1) #2-state
  # sigma_c[2, 2] ~ dexp(1)
  # # sigma_c[1, 2, 2] <- 0 #2-state
  # # sigma_c[2, 1, 2] <- 0 #2-state
  # # sigma_c[2, 2, 2] ~ dexp(1) #2-state

  # Rho[1, 1] <- 1.0
  # Rho[1, 2] ~ dunif(-0.8, 0.8)
  # Rho[2, 1] <- Rho[1, 2]
  # Rho[2, 2] <- 1.0
  # # Rho[1, 1, 2] <- 1.0 #2-state
  # # Rho[1, 2, 2] ~ dunif(-0.8, 0.8) #2-state
  # # Rho[2, 1, 2] <- Rho[1, 2, 2] #2-state
  # # Rho[2, 2, 2] <- 1.0 #2-state

  # for(j in 1:2){ # j is the index for lat and lon
  #   sd[j] ~ dunif(0, 0.1)
  #   tau[j] <- 1/(sd[j]*sd[j])
  #   x.sd[j] ~ dunif(0, 0.1)
  #   x.tau[j] <- 1/(x.sd[j]*x.sd[j])
  # } 
  x.sd ~ dunif(0, 0.1)
  x.tau <- 1/(x.sd*x.sd)
  # Sigma = sigma_c %*% Rho %*% sigma_c
  # Sigma[1:2, 1:2, 2] <- sigma_c[,,2] %*% Rho[,,2] %*% sigma_c[,,2] #2-state 

  # # priors on strength of autoregressive component for two states (correlation between subsequent movements)
  # gamma[1] ~ dunif(0, 1) #2-state; A high value for gamma[1] describes migratory movements: rapid persistent movement over distance in a direction with high correlation between subsequent movements. 
  # gamma[2] ~ dunif(0, 1) #2-state; A low  value for gamma[2] describes foraging or staging: or relatively slow, more random turning (area-restricted search).
  ##gamma[1] ~ dbeta(2, 1.5)
  ##dev ~ dbeta(1, 1) #from hDCRWS
  ##gamma[2] <- gamma[1] * dev #random deviate to restrict: gamma[1] > gamma[2] #from hDCRWS
  # modify the model so it has only 1 state to choose from
  gamma ~ dunif(0, 1) #1-state prior for gamma (the behavioral-state-specific correlation parameter)
  
  # # priors for drift component (one for each of 2 states, x and y dimensions, so 4 options for drift); follows a normal distribution to allow for positive (NE) or negative (SW) drift; units in dec. deg. 
  # D[1, 1] ~ dnorm(0, 0.1) #2-state; D[coordinate, b[t]]; up to a tenth of a degree additional movement in lat for migration per time step (vague prior): S(-) or N(+), where b[t] 1 = migration, b[t] 2 = staging, coordinate 1 = latitude, coordinate 2 = longitude
  # D[1, 2] ~ dnorm(0, 0.1) #2-state; D[coordinate, b[t]]; up to a tenth of a degree additional movement in lon for migration per time step (vague prior): S(-) or N(+), where b[t] 1 = migration, b[t] 2 = staging, coordinate 1 = latitude, coordinate 2 = longitude
  # D[2, 1] ~ dnorm(0, 0.1) #2-state; D[coordinate, b[t]]; up to a tenth of a degree additional movement in lon for staging   per time step (vague prior): W(-) or E(+), where b[t] 1 = migration, b[t] 2 = staging, coordinate 1 = latitude, coordinate 2 = longitude
  # D[2, 2] ~ dnorm(0, 0.1) #2-state; D[coordinate, b[t]]; up to a tenth of a degree additional movement in lat for staging   per time step (vague prior): W(-) or E(+), where b[t] 1 = migration, b[t] 2 = staging, coordinate 1 = latitude, coordinate 2 = longitude
  
  # # priors for drift component; units in dec. deg.
  # D ~ dunif(-0.5, 0) #1-state; dunif(-0.5, 0) D[1] for south; D on displace[,2] for west only
  # D[2] ~ dnorm(0, 0.1) #1-state; dunif(-0.5, 0) for west

  # k is the individual index, N is number of animals
  for(k in 1:N){
  # get the first location for each individual
    first.loc[k, 1] <- y[Yidx[k], 1]
    first.loc[k, 2] <- y[Yidx[k], 2]

    # mute measurement error components
    #logpsi[k] ~ dunif(-10, 10) #for identifiability: inflation/deflation factor for estimation error
    #psi[k] <- exp(logpsi[k]) #psi is animal-specific tag variance (measurement-error scalar) due to noise or attachment

    # b2[Xidx[k]] <- 1 #2-state; b2 is the constrained transition vector, the 1/0 outcome of switching, given a previous transition; b2[t]=b[t]-1 because staging (b[t]=2) is likely to switch, migration (b[t]=1) is not
    # b[Xidx[k]] <- 2 #2-state; assign first behavioral state for each individual (b is unobserved true state): animal in behavioral state 2 (staging) at first detection, not behavioral state 1 (migration)
    
    for(j in 1:2){ # j is the index for lat and lon
      # mute measurement error components
      #x[Xidx[k], j] ~ dt(first.loc[k, j], itau2[Xidx[k], j]*psi[k], nu[Xidx[k], j]) #Prior for first location; nu is the observation/measurement error in observed lat/lon following bivariate normal distribution
      # get x for the first location #x is the unobserved location state at regular time step t
      x[Xidx[k], j] <- first.loc[k, j] #Xidx references the position in the vector of the first detection for each individual in the regular time steps 
    ## Process equation (describes transitions between estimated unobserved locational, x, and behavioral, b, states at regular time intervals: estimates gamma, phi - not alpha, elements of Sigma)
    # get x for the second recorded location, which does not use autoregressive component (gamma)
    x[(Xidx[k]+1), j] ~ dnorm(x[Xidx[k],j], x.tau) #x is the unobserved true location #dmnorm.vcov not dmnorm as in hDCRWS
    # mute loglikelihood monitor for first two locations
    #log_lik[(Xidx[k])] <- logdensity.mnorm(x[(Xidx[k]), 1:2], x[(Xidx[k]), ], iSigma[, ])
    #log_lik[(Xidx[k]+1)] <- logdensity.mnorm(x[(Xidx[k]+1), 1:2], x[(Xidx[k]+1), ], iSigma[, ])
    }

    ### for regular time steps, t, add effect of switching to remaining behavioral states (second time step to last detection)
    for(t in (Xidx[k]+1):(Xidx2[k]-1)){ #Xidx2 indexes the last detection (of the season) for each individual; t is the regular (daily) time step between the first and last detection
      # phi2[t] <- max(phi*b2[t-1], 0.001) #2-state; phi = prob of switching between states; phi2 = constrained prob of transition; given a transition from the previous behavior, the probability of switching again cannot exceed 0.001 (once migrating, not likely to switch back to staging)
      # b2[t] ~ dbern(phi2[t]) #2-state; b2 is the constrained transition vector, the 1/0 outcome of switching, given a previous transition
      # b[t] <- b2[t] + 1 #2-state; behavioral state at timestep is equal to 1 or 2; b2[t]=b[t]-1 because staging (b[t]=2) is likely to switch, migration (b[t]=1) is not
	  
      # displacement at time t + 1, from t, includes an autoregressive component (gamma), drift (D), and a max hop distance (difference between previous and current)
      # for x and y dims (lat and lon). #displacement for each time step cannot exceed 10 (equivalent to ~46.25 km/h max speed)
      # combine x and y dims. and add displacement to the previous time step, t, to feed into dmnorm()
      displace[t, 1] <- max(-10, min((x[t, 1] - x[t-1, 1])*gamma, 10)) #displace[t, 1] = latitude; displacement is the difference between lat/lons in estimated location each day. Drift removed.
      displace[t, 2] <- max(-10, min((x[t, 2] - x[t-1, 2])*gamma, 10)) #displace[t, 2] = longitude; displacement is the difference between lat/lons in estimated location each day. Drift removed.
      # combine lat and lon and add displacement to the previous time step, t, to feed into dmnorm()
      x.mn[t, 1:2] <- x[t, 1:2] + displace[t, 1:2]
      # get x for t + 1
    for(j in 1:2){ # j is the index for lat and lon
      x[t+1, j] ~ dnorm(x.mn[t, j], x.tau) #could use tau[j] to group by lat/lon; monitored estimates of x will correspond to the estimated daily lat/lons of each individual in series
    }
    }

    ## Measurement equation
    for(i in (Yidx[k]+1):(Yidx[k+1]-1)){ #Yidx is the index of the first detection for each individual
      for(j in 1:2){ #j is lat:lon (1:2)
        # mute the component that interpolates multiple detections within a single time step
        #yhat[i, j] <- w[i]*x[idx[i], j] + (1 - w[i])*x[idx[i+1], j] #w is the proportion of the time step elapsed prior to the ith observation
        yhat[i, j] <- x[idx[i], j] #yhat is position model, feeds into movement process model
        # mute the measurement error component
        #y[i, j] ~ dt(yhat[i, j], itau2[i, j]*psi[k], nu[i, j]) #for each observed, temporally irregular location y, i is the detection, given the previous detection j; #dt(mu, tau, k) is the Student t-distribution
        y[i, j] ~ dnorm(yhat[i, j], tau) #could use tau[j] to group by lat/lon; y is the outcome of normally distributed observed longitudinal and latitudinal error tau around yhat, the estimated true location (position)
      }
    }
  }
}

  if (is.R()){
    filename <- file.path(dir.out, "MOVE.bug")} #removes comments from file
write.model(MOVE, filename)
# inits <- list(list(gamma=c(0.5, 0.5), x=x, D=matrix(c(0, 0, 0, 0), 2, 2), sigma_c=array(c(0.1, NA, NA, 0.1, 0.1, NA, NA, 0.1), c(2, 2, 2)), phi=0.1, sd=0.01, #sd=0.000001 #removes sd prior (now a constant)
                   # Rho=array(c(NA, NA, 0.01, NA, NA, NA, 0.01, NA), c(2, 2, 2))), #2-state
inits <- list(list(gamma=c(0.5), D=c(0,0), sd=0.01, x.sd=0.01), # #D=c(0) #sd=0.000001, removes sd prior (now a constant)
              list(gamma=c(0.5), D=c(0,0), sd=0.01, x.sd=0.01), # #1-state
              list(gamma=c(0.5), D=c(0,0), sd=0.01, x.sd=0.01)) # #1-state
str(jags.data <- list(Xidx=Xidx, Xidx2=Xidx2, Yidx=Yidx, y=as.matrix(y), idx=idx, N=N))

parameters <- c("gamma", "Sigma", "D", "x", 'tau', 'x.tau') #, 'sigma_c', 'Rho' #1-state removes "b", "phi", adds adds ('sd',) 'sigma_c', 'Rho'

#(.packages()) #list loaded libraries
library(dclone)
library(rjags)

stime=Sys.time() #added "try" to allow execution to continue after error - for evalTime of error
MOVE <- try(jags.parfit(cl=makePSOCKcluster(names=nc), 
                           data   = jags.data,
                           params = parameters,
                           model  = filename,
                           inits  = inits,
                           n.chains = nc,
                           n.adapt  = n.burnin,
                           n.update = 0,
                           thin     = n.thin,
                           n.iter   = n.iter)) #not including burn in
(evalTime = difftime(Sys.time(),stime)) #Time difference of 

save(MOVE, file = file.path(dir.local, paste0(out.name, out.nam2), paste0(out.name, out.nam2,'.RData'))) #save to local data drive on high performance machine

#load the movement model
#load(file.path(dir.local, paste0(out.name, out.nam2), paste0(out.name, out.nam2,'.RData'))) #load from local data drive on high performance machine

# For output from jags.parfit ####
# Check convergence on global params ###
(n.mcmc <- length(MOVE) * nrow(MOVE[[1]])) #Posterior sample size
(param.names <- dimnames(MOVE[[1]])[[2]][!grepl('b\\[', dimnames(MOVE[[1]])[[2]]) & 
                                                !grepl('x\\[', dimnames(MOVE[[1]])[[2]])])

pdf(file.path(dir.out, paste0('occ_post_', out.name, '_traceplot_global.pdf')))
plot(as.mcmc.list(lapply(MOVE, function(x) x[, param.names])))
dev.off()

MOVE.mcmc <- lapply(MOVE, function(x) as.mcmc(x))
MOVE.mcmc.global <- lapply(MOVE, function(x) as.mcmc(x[,param.names]))

fm.sum <- MCMCsummary(MOVE.mcmc.global) #MCMCvis
write.csv(fm.sum, file.path(dir.out, paste0('occ_post_', out.name, '_rhat_global.csv')))

# Post-processing ####
# create a vector of the number of individuals tagged by each month,
# to ensure number of individuals observed is divided by the correct denominator to get occupancy

##count the number of tagged individuals per month, from first month of detection to last month of detection
id_mo.mat <- as.matrix(table(unique(allindvs_ordered_filtered_nodups[,c('id','month')]))) #same as above

id_mo.dat <- data.frame(
  mo_min = sapply(split(allindvs_ordered_filtered_nodups[,c('id', 'month')], allindvs_ordered_filtered_nodups$id),
                  function(x) min(x[,'month'])),
  mo_max = sapply(split(allindvs_ordered_filtered_nodups[,c('id', 'month')], allindvs_ordered_filtered_nodups$id),
                  function(x) max(x[,'month']))
)

for (i in 1:nrow(id_mo.mat)) { #fill in any middle months with missing detections (between first and last detection)
  id_mo.mat[i,id_mo.dat[i,1]:id_mo.dat[i,2] - 
              (min(allindvs_ordered_filtered_nodups$month)-1)] <- 1
}

N_bymonth <- rep(0,12)

N_bymonth[as.numeric(dimnames(id_mo.mat)[[2]])] <- colSums(id_mo.mat)
N_bymonth #accounts for months with missing detections between first and last

## use posterior estimates from JAGS model to obtain estimates for user-defined spatial units ####################################
# convert output from JAGS to create arrays (one for lat; one for lon) that have the posterior samples (rows)
# for each time step of the movement model (columns), and for each individual (depth)
posts_lat <- array(NA, c(n.iter/n.thin*nc, season_length, N))
posts_lon <- array(NA, c(n.iter/n.thin*nc, season_length, N))

stime=Sys.time() #using jags.parfit
for(e in 1:N){
  for(i in Zidx[e]:season_length){
    posts_lat[,i,e] <- c(as.array(MOVE[,paste("x[", ((e-1)*season_length + i), "," , 1, "]", sep=""),]))
    posts_lon[,i,e] <- c(as.array(MOVE[,paste("x[", ((e-1)*season_length + i), "," , 2, "]", sep=""),]))
  }
}
(evalTime = difftime(Sys.time(),stime)) #Time difference of 2.4-5 mins

#constrain the predictions to times with detections to avoid skew from including data past where we have detections

start_x <- end_x <- c()
for(i in 1:e){

  start_x[i] <- allindvs_ordered_filtered_nodups$days_since[Yidx[i]]
  end_x[i] <- allindvs_ordered_filtered_nodups$days_since[Yidx[i + 1] - 1] - allindvs_ordered_filtered_nodups$days_since[Yidx[i]] + 1

}

#cut the predictions > the final data point

posts_lat_t <- posts_lat
posts_lon_t <- posts_lon

stime=Sys.time() 
for(i in 1:e){
  for(j in 1:season_length){
    if(j > (start_x[i] + end_x[i])){
      posts_lat_t[ , j, i] <- NA
      posts_lon_t[ , j, i] <- NA}
  }
}
(evalTime = difftime(Sys.time(),stime)) #Time difference of 5 secs

# read a file that has the extent of user-defined spatial units as rows, labeled "xmin", "ymax", "xmax", "ymin"
# can use the field calculator in QGIS to get a data file with the extent of each cell and centroid coordinates
# e.g. x_min($geometry) on polygons and $x on a .shp for centroids

# BOEM.sf <- st_read(file.path(dir.data, 'BOEM_halfdeg_grid_latlon_2021')) #file path too long
BOEM.sf <- st_read('P:/SCRAM/BRI/Movement/SCRAM Movement Models/data',
                   'BOEM_halfdeg_grid_latlon_2021')

# ## Create bounding box based on active stations during study period (as USFWS recommended)
# BOEM.sf <- BOEM.sf[BOEM.sf$top > 36.4914 & BOEM.sf$bottom < 42.2453,] #keep bounding grid cells that contain those latitudes

BOEM <- st_drop_geometry(BOEM.sf)
#BOEM <- read.csv('data/BOEM_halfdeg_grid_latlon_att.csv', header=TRUE)
# label columns
colnames(BOEM)[1] <- "xmin"
colnames(BOEM)[2] <- "ymax"
colnames(BOEM)[3] <- "xmax"
colnames(BOEM)[4] <- "ymin"

# this script makes it possible to loop through every nth value
# useful for running through a spatial input file that has too many rows to run wihtin a reasonable time
nth <- function(x,n){
  x[x%%n==0]
}
x = 1:length(BOEM[,1])
index <- nth(x,1)

## Propagate stochasticity and uncertainty through the posterior estimates of position ###########################
# extrapolation includes JAGS model uncertainty and the uncertainty from estimating the proportion of
# individuals who likely crossed into the specified extent; gives estimates by month instead of year
posts_latb <- posts_lat_t[sample(1:(n.iter/n.thin*nc), 1000), , ] #sample from multiple chains
posts_lonb <- posts_lon_t[sample(1:(n.iter/n.thin*nc), 1000), , ] #sample from multiple chains

source(file.path(dir.work, 'Pr_occupancy_season_trunc.R'))
occ_post <- array(0, c(1000, 12, length(BOEM[,1])))

stime=Sys.time() 
for(z in 1:length(index)){
  #for(z in 1:1){
  lat_input_min <- BOEM$ymin[z]
  lat_input_max <- BOEM$ymax[z]
  lon_input_min <- BOEM$xmin[z]
  lon_input_max <- BOEM$xmax[z]
  occ_post[, ,z] <- Pr_occupancy_season(posts_latb, posts_lonb, lat_input_min, lat_input_max, lon_input_min, lon_input_max, side, n.iter, N_bymonth, season_length, "time_step_bymo")
}
(evalTime = difftime(Sys.time(),stime)) #Time difference of ~25 mins.

saveRDS(occ_post, file.path(dir.out, paste0('occ_post_', out.name, '.rds'))) #this file gets converted into what SCRAM needs to assess movement uncertainty

## Map/visualize posterior results ####

#pull in a basemap
require(rnaturalearth)
require(viridis)

world <- ne_countries(continent = 'North America', returnclass = 'sf')
ocean <- ne_download(scale = 110, type = 'ocean', category = 'physical', returnclass = 'sf')

#add the crs
GCSWGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

#plot the position posteriors
x_post <- data.frame(x = c(posts_lonb[, , ]), y = c(posts_latb[, , ]))
x_post <- remove_missing(x_post)
x_post <- st_as_sf(x_post, coords = c(1, 2), crs = GCSWGS84)

bb <- st_bbox(x_post)

#intersect x_post with the BOEM shapefile to determine what the overall activity per grid cell should be
#look at the median estimates

med_lat <- apply(posts_latb, c(2,3), mean)
med_lon <- apply(posts_lonb, c(2,3), mean)

x_id <- data.frame(ID = rep(1:N, each = season_length), lat = c(med_lat[,]), lon = c(med_lon[,]))
x_id <- remove_missing(x_id)
x_id <- st_as_sf(x_id, coords = c(3, 2), crs = GCSWGS84)
x_id$Time <- 1:nrow(x_id)

plot(x_id['Time'])
dev.off()

bb <- st_bbox(x_id)

jpeg(filename = file.path(dir.out, paste0('occ_post_', out.name, '_med_ID.jpg')))
ggplot() +
  geom_sf(data = ocean)  +
  geom_sf(data = x_id, aes(col = ID)) + coord_sf(xlim = bb[c(1,3)], ylim = bb[c(2,4)], expand = T) + 
  scale_color_viridis_c(option = "magma",begin = 0.1)
dev.off()

#plot time by ID
dir.create(file.path(dir.out, paste0('occ_post_', out.name, '_med_ID')))

for (i in 1:N) {
  jpeg(filename = file.path(dir.out, paste0('occ_post_', out.name, '_med_ID'), 
                            paste0('occ_post_', out.name, '_med_ID', uni_inds[i], '.jpg')), 
       height=7, width=5, units = "in", res=300) 
  med_ID <- ggplot() +
    geom_sf(data = ocean)  +
    geom_sf(data = x_id[x_id$ID==i,], aes(col = Time)) + coord_sf(xlim = bb[c(1,3)], ylim = bb[c(2,4)], expand = T) + 
    scale_color_viridis_c(option = "magma",begin = 0.1) +
    ggtitle(uni_inds[i])
  print(med_ID) #otherwise blank
  dev.off()
}

#now let's look at the upper 95% (so extreme NE of the positions)

med_lat <- apply(posts_lat_t, c(2,3), function(x) quantile(x, probs = 0.975, na.rm = TRUE))
med_lon <- apply(posts_lon_t, c(2,3), function(x) quantile(x, probs = 0.975, na.rm = TRUE))
#takes a couple minutes each

x_id <- data.frame(ID = rep(1:N, each = season_length), lat = c(med_lat), lon = c(med_lon))
x_id <- remove_missing(x_id)
x_id <- st_as_sf(x_id, coords = c(3, 2), crs = GCSWGS84)

plot(x_id['ID'])
dev.off()

bb <- st_bbox(x_id)

jpeg(filename = file.path(dir.out, paste0('occ_post_', out.name, '_97_5_ID.jpg')))
ggplot() +
  geom_sf(data = ocean)  +
  geom_sf(data = x_id, aes(col = ID)) + coord_sf(xlim = bb[c(1,3)], ylim = bb[c(2,4)], expand = T) + 
  scale_color_viridis_c(option = "magma",begin = 0.1)
dev.off()

#plot occu_post features

hist(occ_post, breaks = length(seq(from=0,to=max(occ_post, na.rm=T),by=0.01)), freq=F)
hist(occ_post, breaks = length(seq(from=0,to=1,by=0.001)), xlim=c(0,1), freq=F)
dev.off()

jpeg(filename = file.path(dir.out, paste0('occ_post_', out.name, '_hist.jpg')))
hist(occ_post, breaks = length(seq(from=0,to=0.1,by=0.0001)), xlim=c(0,0.1), freq=F, 
     xlab="Cumulative Daily Occupancy (Truncated)", main="SCRAM 2.0") #title hard-coded (needs manual update)!
dev.off()

occu_med <- apply(occ_post, c(3, 2), function(x) mean(x, na.rm = TRUE))
occu_sd <- apply(occ_post, c(3, 2), function(x) sd(x, na.rm = TRUE))
occu_cv <- occu_sd/occu_med

BOEM.sf$occu_med <- apply(occ_post, c(3, 2), function(x) mean(x, na.rm = TRUE))
BOEM.sf$occu_sd <- apply(occ_post, c(3, 2), function(x) sd(x, na.rm = TRUE))
BOEM.sf$occu_cv <- BOEM.sf$occu_sd/BOEM.sf$occu_med

BOEM.sf <- cbind(BOEM.sf, occu_med, occu_sd, occu_cv)
#str(BOEM.sf) 

#Plot average monthly occupancy days ("cumulative daily occupancy" averaged across months) and export to shapefile for mapping 
BOEM.sf$occu_med_mo <- apply(BOEM.sf$occu_med[,apply(BOEM.sf$occu_med, 2, sum, na.rm=TRUE)>0], 1, mean, na.rm=TRUE)
summary(BOEM.sf$occu_med_mo)

dir.create(file.path(dir.out, paste0('occ_post_', out.name, '_med_shp')))
jpeg(filename = file.path(dir.out, paste0('occ_post_', out.name, '_med_shp'), paste0('occ_post_', out.name, '_med.jpg')))
plot(na.omit(BOEM.sf["occu_med_mo"]), logz = TRUE)
dev.off()

st_write(BOEM.sf["occu_med_mo"], file.path(dir.out, paste0('occ_post_', out.name, '_med_shp'), paste0('occ_post_', out.name, '_med.shp')), append=F)

#Plot all months
for (m in which(N_bymonth>0)) {
  jpeg(filename = file.path(dir.out, paste0('occ_post_', out.name, '_X', m, '_med.jpg')))
  try(plot(BOEM.sf[paste0('X', m, '')], logz = TRUE))
  dev.off()
  jpeg(filename = file.path(dir.out, paste0('occ_post_', out.name, '_X', m, '_sd.jpg')))
  try(plot(BOEM.sf[paste0('X', m, '.1')]))
  dev.off()
  jpeg(filename = file.path(dir.out, paste0('occ_post_', out.name, '_X', m, '_cv.jpg')))
  try(plot(BOEM.sf[paste0('X', m, '.2')], logz = TRUE))
  dev.off() 
}

for (m in which(N_bymonth>0)) {
  dir.create(file.path(dir.out, paste0('occ_post_', out.name, '_X', m, '_med_shp')))
  try(st_write(BOEM.sf[paste0('X', m, '')], file.path(dir.out, paste0('occ_post_', out.name, '_X', m, '_med_shp'), paste0('occ_post_', out.name, '_X', m, '_med.shp')), append=F))
}

save(occ_post, file=file.path(dir.out, paste0('occ_post_', out.name, '.RData')))
#load(file=file.path(dir.out, paste0('occ_post_', out.name, '.RData')))

## From create_baked_files_CRM.R
#after runing the species specific movement data analysis script to get occ_post:
monthLabels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",  "Aug", "Sep", "Oct", "Nov", "Dec")

dir.create(file.path(dir.out, "movements2_0"))

#Remove months with five or less tracked individuals for collision risk estimation 
#Assign NA's to grid cells outside study area
N_thresh <- 5 #hard-coded (needs manual update)!!
occ_post[,N_bymonth<=N_thresh,] <- NA #updated to soft-code: months without detections
occ_post[,,which(!(BOEM.sf$top > 36.4914 & BOEM.sf$bottom < 42.2453))] <- NA #updated to soft-code: grid cells outside study area
colnames(occ_post) <- monthLabels

save(occ_post, file=file.path(dir.out, paste0('occ_post_', out.name, '_studyarea', '.RData')))
#load(file=file.path(dir.out, paste0('occ_post_', out.name, '_studyarea', '.RData')))

for(i in 1:length(occ_post[1,1,])){
  write.csv(occ_post[,,i], file=file.path(dir.out, "movements2_0", paste0("MovementBaked_", out.name, "_", i, ".csv")), row.names=FALSE)
}

#Recalculate cumulative daily occupancy  averaged across only months meeting sample size threshold 
occu_med <- apply(occ_post, c(3, 2), function(x) mean(x, na.rm = TRUE))
occu_sd <- apply(occ_post, c(3, 2), function(x) sd(x, na.rm = TRUE))
occu_cv <- occu_sd/occu_med

BOEM.sf$occu_med <- apply(occ_post, c(3, 2), function(x) mean(x, na.rm = TRUE))
BOEM.sf$occu_sd <- apply(occ_post, c(3, 2), function(x) sd(x, na.rm = TRUE))
BOEM.sf$occu_cv <- BOEM.sf$occu_sd/BOEM.sf$occu_med

BOEM.sf <- cbind(BOEM.sf, occu_med, occu_sd, occu_cv)
#str(BOEM.sf) 

#Plot average monthly occupancy days ("cumulative daily occupancy" averaged across months) and export to shapefile for mapping 
BOEM.sf$occu_med_mo <- apply(BOEM.sf$occu_med[,apply(BOEM.sf$occu_med, 2, sum, na.rm=TRUE)>0], 1, mean, na.rm=TRUE)
summary(BOEM.sf$occu_med_mo)

BOEM.sf[!(BOEM.sf$top > 36.4914 & BOEM.sf$bottom < 42.2453),"occu_med_mo"] <- NA #reduce to study area

jpeg(filename = file.path(dir.out, paste0('occ_post_', out.name, '_med_shp'), paste0('occ_post_', out.name, '_med_N.jpg')))
plot(BOEM.sf["occu_med_mo"], logz = TRUE)
dev.off()

st_write(BOEM.sf["occu_med_mo"], file.path(dir.out, paste0('occ_post_', out.name, '_med_shp'), paste0('occ_post_', out.name, '_med_N.shp')), append=F)

BOEM.sf <- BOEM.sf[BOEM.sf$top > 36.4914 & BOEM.sf$bottom < 42.2453,] #reduce to study area extent for mapping
st_write(BOEM.sf["occu_med_mo"], file.path(dir.out, paste0('occ_post_', out.name, '_med_shp'), paste0('occ_post_', out.name, '_med_studyarea.shp')), append=F)
