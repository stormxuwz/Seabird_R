library(dplyr)

read_csv_file <- function(filename) {
  data <- read_csv(filename)
  return(data)
}


read_cnv_file <- function(filename) {
  con <- file(filename, "r")
  lines <- readLines(con)
  
  feature_columns <- c()
  
  n <- length(lines)
  
  counter <- 1
  for (line in lines) {
    if (startsWith(line, "# name")) {
      feature_columns <- c(feature_columns, read_cnv_metaline(line))
    }
    
    if (line == "*END*") {
      break
    }
    counter <- counter + 1
  }
  
  counter <- counter + 1
  
  data <- lines[counter: length(lines)] |> sapply(read_cnv_dataline) |> t() |> data.frame(row.names = NULL)
  names(data) <- feature_columns
  
  # add data frame rename
  
  close(con)
  return(data)
}

read_cnv_metaline <- function(line) {
  tmp <- strsplit(line, "=")[[1]][2] |> trimws()
  return(strsplit(tmp, ":")[[1]][1])
}

read_cnv_dataline <- function(line) {
  tmp <- strsplit(line, " ")[[1]]
  result <- tmp[sapply(tmp, function(x) {x != ""})] |> 
    sapply( function(x){as.double(x)}) |> 
    setNames(NULL)
  return(result)
}

rename_data <- function(df) {
  df <- rename(df, 
    Depth = depFM,
    Temperature = t090C,
    Fluorescence = flSP
  )
  return(df)
}
