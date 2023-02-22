#' Data normalizer
#'
#' This function takes as input a vector and returns it normalized between
#' a specified range
#'
#' @param dat the vector to normalize
#' @param norm_range the range used for normalization
#'
#' @return vector normalized between norm_range
#'
#' @examples
#' normalizer(c(1:10,10:1,1:10,10:1,1:10,10:1,1:10,10:1,1:10,10:1))
#'
#' @importFrom Convolutioner MA
#'
#' @export


normalizer <- function(dat, norm_range = c(0,1)) {

  if (sum(is.na(dat)) || sum(is.character(dat))) {
    dat = dat[-which(is.na(dat))]
  }

  holder <- c()

  for (i in 1:length(dat)) {
    holder[i] <- norm_range[1] + ((dat[i]- min(dat, na.rm = T))*(norm_range[2]- norm_range[1]))/(max(dat, na.rm = T)- min(dat, na.rm = T))
  }

  return(holder)

}

#' Trend finder
#'
#' This function takes as input a vector and returns
#' the trend of each column expressed as the difference
#' between two consecutive elements
#'
#' @param dat the vector to analyze
#'
#' @return vector containing the trend of the each column
#'
#' @examples
#' trend(c(1:10,10:1,1:10,10:1,1:10,10:1,1:10,10:1,1:10,10:1))
#'
#' @export

trend <- function(dat) {

  if (sum(is.na(dat)) || sum(is.character(dat))) {
    dat = dat[-which(is.na(dat))]
  }

  holder <- c()
  for (i in 1:(length(dat)-1)) {

    holder[i] <- (dat[i+1] - dat[i])
  }

  return(holder)
}

#' Sign detection
#'
#' This function takes as input a vector and returns the sign of each element
#'
#' @param dat the vector to be used
#'
#' @return vector with the signs of each element of the original matrix
#'
#' @examples
#' sign_detect(c(1:10,10:1,1:10,10:1,1:10,10:1,1:10,10:1,1:10,10:1))
#'
#' @export


sign_detect <- function(dat) {

  if (sum(is.na(dat)) || sum(is.character(dat))) {
    dat = dat[-which(is.na(dat))]
  }

  sign_holder <- c()

  for (i in 1:length(dat)) {
    if (dat[i] > 0) {
      sign_holder <- append(sign_holder, 1)
    } else if (dat[i] < 0) {
      sign_holder <- append(sign_holder, -1)
    } else if (i != 1){
      sign_holder <- append(sign_holder, sign_holder[i-1])
    } else {
      sign_holder <- append(sign_holder, 0)
    }
  }
  return(sign_holder)
}

#' Tidal finder
#'
#' This function takes as input a matrix and returns for
#' each column the end tidals depending of the threshold set. It is
#' possible to set a custom time frame for the search of the tidals.
#' Note: a minimum amount of 45 points are necessary.
#'
#' @param dat the input matrix
#' @param height_threshold the minimum height of the tidal
#' @param refine refine the dataset
#' @param time_frame_index custom time frame
#'
#' @return matrix with the tidals for each column
#'
#' @examples
#' tidal_finder(c(1:10,10:1,1:10,10:1,1:10,10:1,1:10,10:1,1:10,10:1))
#'
#' @export

tidal_finder <- function(dat, height_threshold = 0.2, refine = FALSE, time_frame_index = NA) {

  if (refine) {
    starting_index <- data_indexer(dat, time_frame_index)
    MA_buffer <- 5
  } else {
    starting_index <- c(1,length(dat))
    MA_buffer <- 45
  }
  dat_holder <- normalizer(dat[starting_index[1]:starting_index[2]])
  dat_holder <- MA(dat_holder, buffer_size = MA_buffer)
  trend_holder <- trend(dat_holder)

  tidals <- c()
  signs <- sign_detect(trend_holder)

  for (i in 2:length(trend_holder)) {
    if (signs[i] != signs[i-1] & signs[i] != 0)  {
      tidals <- append(tidals, i-1)
    }
  }

  keep_peaks <- c()

  for (i in 2:length(tidals)) {
    if (signs[tidals[i]] %in% c(0,1) && signs[tidals[i-1]] == -1) {
      if ((max(dat_holder[tidals[i-1]:tidals[i]]) - dat_holder[tidals[i-1]]) > height_threshold) {
        keep_peaks <- append(keep_peaks, c((tidals[i-1] + starting_index[1]-1), (tidals[i] + starting_index[1]-1)))
      }
    }
  }

  if (refine == FALSE) {
    negative_peaks <- c()
    for (i in 2:length(tidals)) {
      if (signs[tidals[i]] %in% c(0,-1) && signs[tidals[i-1]] == 1) {
        if ((min(dat_holder[tidals[i-1]:tidals[i]]) - dat_holder[tidals[i-1]]) < -height_threshold) {
          negative_peaks <- append(negative_peaks, c((tidals[i-1] + starting_index[1]-1), (tidals[i] + starting_index[1]-1)))
        }
      }
    }
    return(c(keep_peaks[2], negative_peaks[length(negative_peaks)-1]))
  }


  return(keep_peaks)

}

#' Tidal analyzer
#'
#' This function takes as input a csv file containing a time
#' column and data columns and returns the position of the end tidals for
#' each data column maximazing the syncronization between data. This function
#' was originally devised for the analysis of the end tidals coming from
#' exhaled breath analyzed through SIFT-MS technology
#'
#' @param setdir working directory
#' @param input_name csv file
#' @param output_name name of the output file
#' @param starting_threshold initial value for the dynamic threshold
#' @param time_frame custom data range from the time column
#' @param out_file flag for the export of a csv file
#'
#' @return csv containing the end tidals, their maximum, average, frequency, and timing
#'
#' @examples
#' data(SIFT_filtered)
#' tidal_analyzer(input_name = head(SIFT_filtered, n = 100), output_name = "out", out_file = FALSE)
#'
#' @importFrom stats sd
#' @importFrom utils read.csv read.delim write.csv
#'
#' @export

tidal_analyzer <- function(setdir = getwd(), input_name = file.choose(),
                           output_name, starting_threshold = 0.03, time_frame = NA,
                           out_file = TRUE){

  oldwd <- getwd()
  on.exit(setwd(oldwd))
  setwd(setdir)

  if (typeof(input_name) == "character") {
    input = read.csv(input_name, header = T)
  } else {
    input = input_name
  }


  tidal_holder <- c()

  outputlist <- list()
  peak_threshold <- c(rep(starting_threshold,(ncol(input)-1)))
  interrupt <- 0
  interation <- 1
  sd_threshold <- 1.1
  time_index <- NA

  if (sum(is.na(time_frame)) == F && length(time_frame) > 2) {
      return(message("Invalide time frame"))
  } else if (sum(is.na(time_frame)) == F) {
     time_index <- c(which(round(input[[1]], digits = 1) == time_frame[1])[1], which(round(input[[1]], digits = 1) == time_frame[2])[1])
  }

  while (interrupt == 0) {
    npeaks_holder <- c()
    for (i in 2:ncol(input)) {

      tidal_holder <- tidal_finder(input[[i]], height_threshold = peak_threshold[i-1], refine = TRUE, time_frame_index = time_index)

      outputdf <- data.frame(rep("", length(tidal_holder)/2))
      names(outputdf) <- "npeaks"
      outputdf$npeaks[1] = length(tidal_holder)/2
      outputdf$freq = rep("", length(tidal_holder)/2)
      outputdf$freq[[1]] = round((length(tidal_holder)/2)/(input[tidal_holder[length(tidal_holder)],1] - input[tidal_holder[1],1]))

      for (j in 1:length(tidal_holder)) {
        if (j %% 2 == 0) {
          outputdf$maximum[j/2] = round(max(input[tidal_holder[j-1]:tidal_holder[j],i]))
          outputdf$average[j/2] = round(mean(input[tidal_holder[j-1]:tidal_holder[j],i]))
          outputdf$std[j/2] = round(sd(input[tidal_holder[j-1]:tidal_holder[j],i]))
          outputdf$ti[j/2] = round(input[tidal_holder[j-1],1],2)
          outputdf$tf[j/2] = round(input[tidal_holder[j],1],2)
        }

      }
      outputlist[[i-1]] = outputdf
      npeaks_holder <- append(npeaks_holder, outputdf$npeaks[1])
    }

    if (sd(npeaks_holder) < sd_threshold) {
      interrupt <- 1
    } else {
      if (max(peak_threshold) < 0.5 && interation < 30) {
        message(paste("Standard Deviation betweek peaks:", round(sd(npeaks_holder), digits = 2)))
        for (i in 1:length(npeaks_holder)) {
          if (as.numeric(npeaks_holder[i]) > mean(as.numeric(npeaks_holder[-i]))) {
            peak_threshold[i] = peak_threshold[i] + 0.0001*sd(npeaks_holder)
          } else if (as.numeric(npeaks_holder[i]) < mean(as.numeric(npeaks_holder[-i]))){
            peak_threshold[i] = peak_threshold[i] - 0.0001*sd(npeaks_holder)
          }
        }
        interation <- interation + 1
      } else {
        message("no consistency between analytes found")
        interrupt <-1
      }
    }

  }
  if (out_file == TRUE) {
    for (i in 2:ncol(input)) {
      write.csv(outputlist[[i-1]], file = paste(paste(output_name, names(input)[i], sep = "_"),".csv", sep = ""))
    }
  }

  return(outputlist)
}


#' data indexer
#'
#' This function takes as input a vector and return the
#' data index according to the selected time frame
#'
#' @param dat input vector
#' @param time_frame_index custom data range from the time column
#'
#' @return a vector indexed according to the specified time frame
#'
#' @examples
#' data_indexer(c(1:10,10:1,1:10,10:1,1:10,10:1,1:10,10:1,1:10,10:1))
#'
#' @export

data_indexer <- function(dat, time_frame_index = NA) {

  if (sum(is.na(dat)) || sum(is.character(dat))) {
    dat = dat[-which(is.na(dat))]
  }

  if (sum(is.na(time_frame_index)) == T) {
    indexes = tidal_finder(dat, height_threshold = 0.1, refine = FALSE, time_frame_index)

    if (length(indexes) != 2) {
      message("Data has not been cleaned due to an anomaly")
      message("Manual time frame is suggested")
      indexes <- c(1, length(dat))
      return(indexes)

    } else if (indexes[1] != indexes[2]) {
      return(indexes)

    } else {
      message("Data has not been cleaned due to an anomaly")
      message("Manual time frame is suggested")
      indexes <- c(1, length(dat))
      return(indexes)
    }

  } else {
    indexes <- c(time_frame_index[1], time_frame_index[2])
    return(indexes)
  }
}




