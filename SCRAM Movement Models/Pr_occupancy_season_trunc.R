### Pr_occupancy_season.R() estimates the probability that an individual of a species will be 
##### within the specified extent (lat_input_min, lat_input_max, lon_input_min, lon_input_max)
### mc is the number of iterations to run to estimate the uncertainty distribution
### posts_lat and posts_lon are the formatted posterior distributions for the estimated lat and long 
##### coordinates for each time step of the movement model
### N is the number of species; season_length is the number of time steps in the movement model
### There are four versions that can be specified (method), which deal with uncertainty differently
##### 'JAGS' extrapolates with full uncertainty; 'individual' only includes uncertainty from estimating 
##### the proprotion of individuals who likely crossed into the specified extent; 'time_step' only includes 
##### uncertanity from estimating locations for each time step in the movement model; 'fast' starts with 
##### posterior means to get rapid point estimates

#EMA: update this to use daily activity estimates

Pr_occupancy_season <- function(posts_latb, posts_lonb, lat_input_min, lat_input_max, lon_input_min, lon_input_max, side, n.iter, N, season_length, method){
  
  #ptm <- proc.time()
  lamin <- (posts_latb - lat_input_min) > 0
  lamax <- (lat_input_max - posts_latb) > 0
  lomin <- (posts_lonb - lon_input_min) > 0
  lomax <- (lon_input_max - posts_lonb) > 0
  inorout <- lamin + lamax + lomin + lomax
  inorout[inorout < 4] <- 0
  inorout[inorout == 4] <- 1
  #proc.time() - ptm
  
  
  # this version of extrapolation collapses the uncertainty over monthly time steps
  # (to get whether there was a detection at least once per month), and propagates the uncertanity from individuals
  if(method=="time_step_bymo"){
    # a n.iter X # ind. matrix for whether or not each ind. was detected at least once during the season
    first_rec_day_mo <- allindvs_ordered_filtered_nodups[order(allindvs_ordered_filtered_nodups$days_since, allindvs_ordered_filtered_nodups$month, allindvs_ordered_filtered_nodups$day, allindvs_ordered_filtered_nodups$hours_dec), c("month", "day")][1,]
    last_rec_day_mo <- allindvs_ordered_filtered_nodups[order(allindvs_ordered_filtered_nodups$days_since, allindvs_ordered_filtered_nodups$month, allindvs_ordered_filtered_nodups$day, allindvs_ordered_filtered_nodups$hours_dec), c("month", "day")][length(allindvs_ordered_filtered_nodups$day),]
    calendar <- cbind(c(1:31, 1:28, 1:31, 1:30, 1:31, 1:30, 1:31, 1:31, 1:30, 1:31, 1:30, 1:31),
                      c(rep(1, 1, 31), rep(2, 1, 28), rep(3, 1, 31), rep(4, 1, 30), rep(5, 1, 31), 
                        rep(6, 1, 30), rep(7, 1, 31), rep(8, 1, 31), rep(9, 1, 30), rep(10, 1, 31), 
                        rep(11, 1, 30), rep(12, 1, 31)))
    mo <- calendar[which(calendar[,2]==first_rec_day_mo[,1]&calendar[,1]==first_rec_day_mo[,2]):which(calendar[,2]==last_rec_day_mo[,1]&calendar[,1]==last_rec_day_mo[,2]), 2]
    inorout_summed <- apply(inorout, 2, rowSums, na.rm=TRUE)
    #create an index for month (each position in the vector is a day) from the first month with a detection until the end of the year
    mo_index <- min(mo):12
    inorout_summed_mo <- mat.or.vec(length(inorout_summed[,1]), length(mo_index))
    for(i in 1:length(mo_index)){
      if(length(which(mo_index[i]==mo))>0){
        inorout_summed_mo[,i] <-  rowSums(inorout_summed[,which(mo_index[i]==mo)], na.rm=TRUE)
      }else{
        inorout_summed_mo[,i] <- 0
      }
    }
    inorout_summed_mo <- t(t(inorout_summed_mo)/N[mo_index])
    inorout_summed_mo_out <- mat.or.vec(length(inorout_summed_mo[,1]), 12)
    inorout_summed_mo_out[, unique(mo_index)] <- inorout_summed_mo
    occ_post_b <- inorout_summed_mo_out
  }
  
  
  # this version of extrapolation is the full method that propagates all uncertainty 
  # the resulting probability distribution includes the uncertainty from location estimation
  # (the proportion of the sampled birds at the location) as well uncertainty from extrapolating from the sample
  # but does not allow one to decompose these sources of uncertainty
  if(method=="JAGS"){
    inorout_summed <- apply(inorout, 3, rowSums)
    inorout_summed[inorout_summed>0] <- 1
    occ_post_b <- mat.or.vec(n.iter, 1)
    for(i in 1:n.iter){
      sum_occ <- sum(inorout_summed[i,])
      
      bin_jags <- function(){
        p ~ dunif(0, 1)
        sum_occ ~ dbinom(p, N)
      }
      
      if (is.R()){
        filename <- file.path(tempdir(), "bin_jags.bug")}
      write.model(bin_jags, filename)
      inits <- list(list(p=0.5))
      data <- list("N", "sum_occ")
      parameters <- c("p")
      bin_jags <- jags(data=data, inits=inits, parameters.to.save=parameters, filename,
                       n.chains=1, n.burnin=100000, n.iter=200000, n.thin=1, DIC=TRUE)
      bin_jags.mcmc <- as.mcmc(bin_jags)
      
      occ_post_b[i] <- sample(bin_jags$BUGSoutput$sims.array[, , "p"], 1)
    }
  }
  
  
  # this version of extrapolation collapses the uncertainty over time steps 
  # (to get whether there was a detection at least once in the season), and propagates the uncertanity from individuals
  if(method=="time_step"){
    # a n.iter X # ind. matrix for whether or not each ind. was detected at least once during the season
    inorout_summed <- apply(inorout, 3, rowSums, na.rm=TRUE)
    inorout_summed[inorout_summed>0] <- 1
    # then take the row means/# of inds. to get the prob. of an ind. being in that space at least once during the season
    # to keep uncertainty; one value for each MCMC step; this version gives the uncertainty from location estimation in full
    # (the proportion of the sampled birds at the location) but does not include uncertainty from extrapolating from the sample
    # to the population
    occ_post_b <- rowSums(inorout_summed, na.rm=TRUE)/N
    # to take the mean of MCMC steps to output a scalar
    #occ_post_b <- mean(rowSums(inorout_summed)/N)
  }
  
  # this version of extrapolation collapses over MCMC steps and propagates uncertainty from time steps 
  # the probability that an individual crosses into the specified extent at least once during the season is then used 
  # to estimate the probabily of an individual crossing into the specified extent
  # this version should give the same mean as method="time_step" if the posterior distributions are approximately symmetrical 
  if(method=="individual"){
    inorout_summed <- apply(inorout, 3, colSums)/n.iter
    # find probability by season; this line can be revised to look at other time periods, or skipped to export daily probabilities
    prnotin <- 1 - inorout_summed
    pr_byind <- 1 - apply(prnotin, 2, prod)
    # get the probability of an individual crossing into the specified extent for each day of the season
    occ_post_b <- sum(pr_byind)/N
  }
  
  # this version of extrapolation simplifies the process by collapsing over MCMC step and comparing the means that are obtained 
  # to the specified extent
  # the resulting logicals are then summed over time steps to get results in terms of occurence at least once during the season
  if(method=="fast"){
    occ_post_b <- length(which(colSums(inorout)>0))/N
  }
  occ_post_b
}
