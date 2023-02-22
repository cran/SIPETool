
#' Time_alignment
#'
#' This function takes as input two data set containing a time vector and
#' a data vector and return the two data sets aligned. This is done
#' by reducing the dimensions of the data set with higher points. The first data set
#' is the one coming from the CPET-ESE and the second one from the SIFT-MS
#'
#' @param Cy CPET-ESE output file
#' @param sift SIFT-MS refined file
#'
#' @return A plot and the SIFT-MS data file resized for the alignment with the CPET-ESE file
#'
#' @examples
#' data(SIFT_time)
#' data(CPET_time)
#' time_filter(CPET_time, SIFT_time)
#'
#' @importFrom utils read.csv read.delim write.csv
#'
#' @export



time_filter <- function(Cy = file.choose(), sift = file.choose()) {

if (typeof(Cy) == "character" && typeof(sift) == "character") {
  Cy = read.csv(Cy, sep = ";")
  sift = read.csv(sift)
}



for (i in 1:length(Cy[[1]])) {
  time_holder <- strsplit(Cy[[1]][i], ":")
  Cy[[1]][i] <- as.numeric(time_holder[[1]][1]) + as.numeric(time_holder[[1]][2])/60
}


index_buffer <- c(1)
for (i in 2:length(Cy[[1]])) {
  for (j in (index_buffer[length(index_buffer)]+1):length(sift[[1]])) {
    if ((abs(sift[[1]][j] - as.numeric(Cy[[1]][i])) < as.numeric(Cy[[1]][i])*0.1/i)) {
      index_buffer <- append(index_buffer, j)
      break
   }
  }
}

sift = sift[index_buffer,]

plot(sift[[1]], Cy[[1]])

return(sift)
}

