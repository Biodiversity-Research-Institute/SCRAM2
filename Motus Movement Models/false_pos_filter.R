### false_pos_filter() identifies "bursts" as the number of detections for an individual at the same site, on the same day, in the same hour, within 1 minute of each other
### The function takes a data frame of detections and a vector specifying an index for each row of the data frame (e.g. 1: length(x[,1]))
### This function can take a long time to run for data frames with greater than ~ 50000 rows

#false_pos_filter <- function(x, input_data){
#  burst_length <- length(which(input_data$site == input_data$site[x] & 
#                                 input_data$id == input_data$id[x] & 
#                                 input_data$month == input_data$month[x] &
#                                 input_data$day == input_data$day[x] &
#                                 abs(input_data$hours_dec - input_data$hours_dec[x]) < (1/60)
#  ))
#  burst_length
#}

false_pos_filter <- function(input_data){
  burst_length <- 1
  species_id <- unique(input_data$id)
  for(i in 1:length(species_id)){
    temp <- input_data[input_data$id == species_id[i], ]
    dup_logical <- duplicated(temp[, c("month", "day", "site")])
    min_logical <- mat.or.vec(length(temp[,1]), 1)
    for(z in 1:length(temp[,1])){
      min_logical[z] <- length(which(temp$id == temp$id[z] & temp$month == temp$month[z] 
                                     & temp$day == temp$day[z] &abs(temp$hours_dec - temp$hours_dec[z]) < (1/60)))
    }
    burst_length <- c(burst_length, min_logical)
  }
  burst_length <- burst_length[-1]
  burst_length
}




