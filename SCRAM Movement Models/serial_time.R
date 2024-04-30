### serial_time() adds a column to the data frame for the number of days or hours since the first record in the data frame
### x is the data frame to which to add the column, and y specifies 'days' or 'hours'

serial_time <- function(x, y){
  if(y!='hour'&y!='day'){stop("y must be 'day' or 'hour'")}
  first_record_date <- x[order(x$month, x$day, x$hours_dec), c("month", "day", "hours_dec")][1,]
  if(y == 'day'){
    # create a matrix that indexes month and day for a non-leap year
    calendar <- cbind(c(1:31, 1:28, 1:31, 1:30, 1:31, 1:30, 1:31, 1:31, 1:30, 1:31, 1:30, 1:31),
                      c(rep(1, 1, 31), rep(2, 1, 28), rep(3, 1, 31), rep(4, 1, 30), rep(5, 1, 31), 
                        rep(6, 1, 30), rep(7, 1, 31), rep(8, 1, 31), rep(9, 1, 30), rep(10, 1, 31), 
                        rep(11, 1, 30), rep(12, 1, 31)))
    
    # days since the first detection (first detection is day 1)
    days_since <- mat.or.vec(length(x[,1]), 1)
    for(i in 1:length(x[,1])){
      days_since[i] <- which(calendar[, 2]==x[i, c("month")]&calendar[, 1]==x[i, c("day")]) -
        which(calendar[, 2]==first_record_date[1, c("month")]&calendar[, 1]==first_record_date[1, c("day")]) + 1
    }
    # correct values for "days since" that happened in the following year from the first record
    days_since[days_since<0] <- days_since[days_since<0] + 365
    x <- cbind(x, days_since)
  }
  if(y == 'hour'){
    # create a matrix that indexes month, day, and hour for a non-leap year
    daysinmo <- c(21, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    temp2 <- NA
    temp3 <- NA
    for(i in 1:12){
      temp <- mat.or.vec(24, daysinmo[i])
      for(e in 1:daysinmo[i]){
        temp[, e] <- rep(e, 1, 24)
      }
      temp2 <- c(temp2, as.vector(temp))
      temp3 <- c(temp3, rep(i, 1, length(as.vector(temp))))
    }
    temp2 <- temp2[-1]
    temp3 <- temp3[-1]
    calendar_byhour <- cbind(temp2, temp3, 1:24)
    
    # hours since the first detection (first detection is hour 1)
    hours_since <- mat.or.vec(length(x[,1]), 1)
    for(i in 1:length(x[,1])){
      hours_since[i] <- which(calendar_byhour[, 2]==x[i, "month"]&calendar_byhour[, 1]==x[i, "day"]
                              &calendar_byhour[, 3]==max(ceiling(x[i, "hours_dec"]), 1)) -
        which(calendar_byhour[, 2]==first_record_date[1, "month"]&calendar_byhour[, 1]==first_record_date[1, "day"]
              &calendar_byhour[, 3]==max(ceiling(first_record_date[1, "hours_dec"]), 1)) + 1
    }
    x <- cbind(x, hours_since)
  }
  x
}