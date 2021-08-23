require(dplyr)
library(dplR)

init_filter <- function(df, depth_threshold) {
  df %>% 
    filter(Depth > depth_threshold - 0.5) %>% 
    filter(Temperature> 0) %>%
    filter(Depth < 1000)
}

separate <- function(df) {
  max_depth_idx <- which.max(df$Depth)
  downcast <- df[1:max_depth_idx, ]
  upcast <- df[(max_depth_idx + 1) : nrow(df), ]
  
  result <- list(
    downcast=downcast,
    upcast=upcast
  )
  return(result)
}

resample <- function(df, interval=0.25) {
  depth <- df$Depth
  df <- subset(df, select=-c(Depth))
  
  # create new depth grid
  new_depth <- seq(from = ceiling(min(depth)), to = max(depth), by = interval)
  
  result <- data.frame()
  
  for(d in new_depth) {
    agg_range <- depth <= (d + interval) & depth >= (d - interval)
    result <- rbind(result, df[agg_range,] %>% colMeans())
  } 
  
  names(result) <- names(df)
  result$Depth <- new_depth
  
  return(result)
}


preprocess <- function(df, config) {
  downcast <- separate(df)$downcast
  pre_data <- init_filter(df, config$bad_depth_threshold)
  
  if (nrow(pre_data) < 1) {
    return(NULL)
  }
    
  pre_data_resample <- resample(pre_data, interval=config$depth_interval)
  
  if (nrow(pre_data_resample) %% 2 > 0) {
    pre_data_resample <- pre_data_resample[1:(nrow(pre_data_resample) - 1), ]
  }
  
  result <- list(
    downcast = downcast,
    filtered_data = smooth(pre_data_resample)
  )
  
  return(result)
}

smooth <- function(df, window_len = 9) {
  feature_names <- names(subset(df, select = -c(Depth)))
  
  prev_idx <- seq(from = window_len, to=2, by = -1)
  post_idx <- seq(from = nrow(df) - 1, to = nrow(df) - window_len + 2, by=-1)
  
  n = nrow(df)
  new_df <- rbind(
    df[prev_idx,],
    df,
    df[post_idx,]
  )
  
  for(f in feature_names) {
    s <- hanning(new_df[[f]], window_len) # use hanning
    df[[f]] <- s[(length(prev_idx) + 1):(n + length(prev_idx))]
  }
  
  return(df)
}


  