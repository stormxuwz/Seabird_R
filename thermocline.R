MAX_TEMPERATURE_CHANGE <- 2.0

create_line <- function(x, method = "regression") {
  ind <- 1:length(x)
  regression_line <- lm(x ~ ind + 1) 
  return(predict(regression_line))
}


calculate_error <- function(x, x_hat) {
  return(
     (x - x_hat) %>% abs() %>% max()
  )
}

get_merge_cost <- function(left_segment, right_segment) {
  merged_segment <- c(left_segment, right_segment)
  line <- create_line(merged_segment)
  
  return(calculate_error(line, merged_segment))
}


generate_segment <- function(x, max_error) {
  n <- length(x)
  
  # initialize segment index
  segment_index <- seq(from = 1, to = n, by = 2) %>% 
    lapply(function(i) {return(c(i, i + 1))})
  
  error_list <- get_all_merge_cost(x, segment_index)
  
  while (TRUE) {
    min_error_index <- which.min(error_list)
    # merge two segments with minimum merge cost
    segment_index <- merge_next(segment_index, min_error_index)
    
    if (length(segment_index) == 3) {
      break
    }
    
    # recalculate merge cost
    error_list <- get_all_merge_cost(x, segment_index)
    
    if (min(error_list) > max_error) {
      break
    }
  }
  
  segment_index <- segment_index %>% 
    lapply(function(seg) {
      return(
        list(
          values = create_line(x[seg]),
          idx = seg
        )
      )
    })
  
  return(segment_index)
}


merge_next <- function(segment_index, index) {
  result <- list()
  i <- 1
  j <- 1
  while (i <= length(segment_index)) {
      new_segment <- segment_index[[i]]
      if (i == index) {
        new_segment <- c(new_segment, segment_index[[i + 1]])
        i <- i + 1
      }
      i <- i + 1
      result[[length(result)+1]] <- new_segment
  }
  
  return(result)
}


get_all_merge_cost <- function(x, segment_index) {
  error_list <- 1:(length(segment_index) - 1) %>% 
    sapply(
      function(i) {
        seg1 <- x[segment_index[[i]]]
        seg2 <- x[segment_index[[i + 1]]]
        return(get_merge_cost(seg1, seg2))
      }
    )
  
  return(error_list)
}

get_gradient_from_segment <- function(segment, depth_interval) {
  return((segment$values[1] - segment$values[2]) / depth_interval)
}


# function to detect thermocline
detect_thermocline <- function(df, config) {
  
  segment_list <- generate_segment(df$Temperature, config$max_error)
  gradient <- sapply(segment_list, get_gradient_from_segment, depth_interval = 0.25)
  max_gradient_index <- which.max(gradient)
  depth <- df$Depth
  
  if (max_gradient_index == 1) {
    # if the first segment has the largest gradient
    segment_depth <- depth[segment_list[[1]]$idx]
    if (max(segment_depth) - min(segment_depth) < 2) {
      # first segment is noise
      # remove segment
    }
  }
  
  result <- list(
    lep_depth = NULL,
    uhy_depth = NULL,
    trm_depth = NULL,
    lep_idx = NULL,
    uhy_idx = NULL,
    trm_idx = NULL,
    num_segments = length(segment_list),
    segment_list = segment_list
  )
  
  
  if (gradient[max_gradient_index] > config$min_trm_gradient) {
    
    segment_list_len <- length(segment_list)
    
    # find the location of thermocline
    result$trm_idx <- segment_list[[max_gradient_index]]$idx %>% mean() %>% as.integer()
    result$trm_depth <- depth[result$trm_idx]
    
    gradient_is_stable <- abs(gradient) < config$stable_gradient
    gradient_is_stable[1] <- abs(gradient[1]) < config$stable_gradient_relaxed
    gradient_is_stable[length(gradient_is_stable)] <- abs(gradient[length(gradient_is_stable)]) < config$stable_gradient_relaxed
    
    
    # find LEP and UHY
    surface_temperature <- segment_list[[1]]$values %>% head(1)
    bottom_temperature <- segment_list[[segment_list_len]]$values %>% tail(1)
    
    lep_index = NULL
    uhy_index = NULL
    
    # find LEP
    for (i in 1: (max_gradient_index - 1)) {
      if (!(gradient_is_stable[i] && abs(tail(segment_list[[i]]$values, 1) - surface_temperature) < MAX_TEMPERATURE_CHANGE )) {
        # segment that is not stable
        if (i > 1) {
          lep_index <- segment_list[[i - 1]]$idx %>% tail(1) # lep is the last point of the previous segment 
          break
        }
        
        if (i == max_gradient_index - 1) {
          lep_index <- segment_list[[i]]$idx %>% tail(1)
        }
      }
    }
    
    # find UHY
    for (i in seq(from = segment_list_len, to = max_gradient_index + 1, by = -1)) {
      if (!(gradient_is_stable[i] && abs(tail(segment_list[[i]]$values, 1) - bottom_temperature) < MAX_TEMPERATURE_CHANGE )) {
        # segment that is not stable
        if (i < segment_list_len) {
          uhy_index <- segment_list[[i + 1]]$idx %>% head(1) # lep is the first point of the previous segment 
          break
        }
        
        if (i == max_gradient_index - 1) {
          uhy_index <- segment_list[[i]]$idx %>% head(1)
        }
      }
    }
    
    if(!is.null(lep_index)) {
      result$lep_depth = depth[lep_index]
      result$lep_idx = lep_index
    }
    
    if(!is.null(uhy_index)) {
      result$uhy_depth = depth[uhy_index]
      result$uhy_idx = uhy_index
    }
  } else {
    print("minimum gradient is too low")
  }
  
  return(result)
}

plot_thermocline <- function(df, thermocline_result) {
  library(ggplot2)
  p <- ggplot(data = df) + geom_line(aes(y=-Depth, x=Temperature)) 
  p <- p + geom_hline(aes(yintercept = -thermocline_result$trm_depth))
  p <- p + geom_hline(aes(yintercept = -thermocline_result$uhy_depth))
  p <- p + geom_hline(aes(yintercept = -thermocline_result$lep_depth))
  p  
}

