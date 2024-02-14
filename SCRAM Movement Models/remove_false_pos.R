### remove_false_pos() removes rows for which the value for the number of detections within a burst 
### is less than or equal to the value specified for y
### This function works with the results of false_pos_filter(), x

remove_false_pos <- function(x, y){
  x[x$burst_length > y, ]
}