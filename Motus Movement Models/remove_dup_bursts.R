### remove_dup_bursts() removes all duplicate detections within a burst, preserving only the first detection.
### The burst can be defined as detections that occur within the same 24 hour period or hour.
### x is a data frame of detections and y is a string specifying whether bursts should be defined by 'day' or 'hour'

# create a function to remove duplicates by day or by hour
# currently this function uses for loops, but could potentially be sped up using apply functions or which(unique())
remove_dup_bursts <- function(x, y){
  if(y!='hour'&y!='day'){stop("y must be 'day' or 'hour'")}
  if(y=='day'){
    bursts <- unique(x[, c('id', 'month', 'day')])
    temp <- x[1,]
    for(i in 1:length(bursts[,1])){
      temp <- rbind(temp, x[x[, 'id'] == bursts[i, 'id']&x[, 'month'] == bursts[i, 'month']
                            &x[, 'day'] == bursts[i, 'day'], ][1,])
    }
    temp <- temp[-1,]
  }
  if(y=='hour'){
    hours_dec_floor <- floor(x$hours_dec)
    x2 <- cbind(x, hours_dec_floor)
    bursts <- unique(x2[, c('id', 'month', 'day', 'hours_dec_floor')])
    temp <- x[1,]
    for(i in 1:length(bursts[,1])){
      temp <- rbind(temp, x[x[, 'id'] == bursts[i, 'id']&x[, 'month'] == bursts[i, 'month']&x[, 'day'] == bursts[i, 'day']
                            &floor(x2[, 'hours_dec_floor']) == floor(bursts[i, 'hours_dec_floor']), ][1,])
    } 
    temp <- temp[-1,]
  }
  temp
}
