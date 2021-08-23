setwd("~/Developer/Seabird_R/")

rm(list = ls())

source("thermocline.R")
source("file_parser.R")
source("preprocessing.R")
source("dcl.R")

require(dplyr)

cnv_file <- "sample.cnv"

preprocessing_config <- list(
  bad_depth_threshold = 3,
  depth_interval = 0.25,
  
  # smoothing_method config has no effect for now as all features are smoothed using hanning
  smoothing_method = list(
    Temperature = c("window", 9),
    Fluorescence = c("window", 9),
    Other = c("window", 9)
  )
)

thermocline_config <- list(
  max_error = 0.3,
  stable_gradient_relaxed = 0.25,
  stable_gradient = 0.1,
  min_trm_gradient = 0.15
)

dcl_config <- list(
  peak_minimum_magnitude = 0.3,
  peak_height = 0.2,
  peak_minimum_interval = 10,
  peak_size = 2.5
)

data <- read_cnv_file(cnv_file) %>% rename_data() %>% preprocess(preprocessing_config)

downcast <- data$downcast
cleandata <- data$filtered_data

# thermocline_result contains the thermocline features
thermocline_result <- detect_thermocline(cleandata, thermocline_config)
print(names(thermocline_result))
plot_thermocline(cleandata, thermocline_result)

# dcl_result contains the dcl features
dcl_result <- detect_dcl(cleandata, downcast, thermocline_result$lep_depth, thermocline_result$lep_depth, dcl_config)
print(names(dcl_result))
plot_dcl(dcl_result, cleandata)

