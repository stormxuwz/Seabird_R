
plot_gaussian_fit <- function(y, y_hat) {
  x = 1: length(y)
  ggplot() + geom_line(aes(x=x,y=y)) + geom_line(aes(x=x,y=y_hat))
}

gaussian_function <- function(x, a, x0, sigma, y0, k) {
  return(a * exp(-(x - x0) ^ 2 / (2 * sigma ^ 2)) + y0 + k * (x - x0))
}

fit_gaussian <- function(x, y, x_mean){
  maxK <- abs(
    (y[1] - y[length(y)]) / (x[1] - x[length(x)])
  )
  
  lower_bound <- c(0, 0, -0.15 * maxK)
  upper_bound <- c(Inf, Inf, 0.15 * maxK)
  
  f <- function(p, x, y){
    d <- gaussian_function(x, p[1], x_mean, p[2], max(y) - p[1], p[3])
    return(sum((d - y) ^ 2))
  }
  
  # a, sigma, k
  a0 <- max(y) - min(y)
  sigma <- length(y) * 0.05
  k <- 0
  par <- optim(par = c(a0, sigma, k), f, x = x, y = y, upper = upper_bound, lower = lower_bound, method="L-BFGS-B")$par
  
  fit_y <- gaussian_function(x, par[1], x_mean, par[2], max(y) - par[1], par[3])
  # plot_gaussian_fit(y, fit_y)
  print(par)
  result <- list(
    popt= par,
    fit_y=fit_y
  )
  
  return(result)
}


fit_shape <- function(values, direction) {
  y <- values
  x <- 1: length(values)
  
  if (direction == "left") {
    fit_result <- fit_gaussian(x, y, length(values))
  } else {
    fit_result <- fit_gaussian(x, y, 1)
  }
  
  result <- list(
    fit_y = fit_result$fit_y,
    x_for_fit = x,
    y_for_fit  = y,
    popt = fit_result$popt
  )
 
  return(result)
}

zero_crossing <- function(x, mode) {
  n_x <- length(x)
  all_index <- which( x[2 : n_x] * x[1: (n_x - 1)] < 0)
  
  if(mode == 0) {
    return(all_index + 1)
  } else if (mode == 1) {
    return(all_index[x[all_index] > 0] + 1)
  } else if (mode == 2) {
    return(all_index[x[all_index] < 0] + 1)
  } else {
    return(NULL)
  }
}

fit_error <- function(x, xhat) {
  return(cor(x, xhat))^2
}


detect_peak <- function(values, config) {
  peak_minimum_mangitude <- config$peak_minimum_magnitude
  peak_height <- config$peak_height
  peak_minimum_interval <- config$peak_minimum_interval
  peak_size <- config$peak_size
  
  # the minimum magnitude of a peak
  threshold <- (max(values) - min(values)) * peak_minimum_mangitude + min(values)
  peak_height_threshold = (max(values) - min(values)) * peak_height
  
  # filter unqualified peaks
  x_gradient <- diff(values) # x_i - x_{i - 1}
  raw_peaks <- zero_crossing(x_gradient, 1)
  raw_peaks <- peak_raw_filters(raw_peaks, values, threshold, peak_minimum_interval)
  raw_peaks <- c(0, raw_peaks, length(values))
  
  # filter peaks
  while (TRUE) {
      peak_heights <- find_peak_height(values, raw_peaks)
      
      if (length(peak_heights) == 0) {
        break
      }
      
      minimum_peak_height_index <- which.min(peak_heights)
      if (peak_heights[minimum_peak_height_index] < peak_height_threshold) {
        raw_peaks <- c(raw_peaks[1: minimum_peak_height_index], raw_peaks[(minimum_peak_height_index + 2): length(raw_peaks)])
      } else {
        break
      }
  }
  
  boundaries <- find_boundaries(values, raw_peaks, peak_size)
  peak_features <- extract_peak_features(boundaries)
  return(peak_features)
}

extract_peak_features <- function(boundaries) {
  result <- data.frame()
  
  for(peak in boundaries) {
    one_result <- c(
      peakIndex = peak$middleNode,
      leftSigma = peak$leftShape$popt[2],
      rightSigma = peak$rightShape$popt[2],
      leftErr = peak$leftShape_err,
      rightErr = peak$rightShape_err,
      leftIndex_fit = peak$leftBoundary_fit,
      rightIndex_fit = peak$rightBoundary_fit
    ) %>% t() %>% data.frame()
    
    result <- rbind(result, one_result)
  }
  
  return(result) 
}


find_peak_height <- function(values, raw_peaks) {
  peak_heights <- c()
  
  for (i in 2:(length(raw_peaks) - 1)) {
    prev_peak <- raw_peaks[i - 1]
    curr_peak <- raw_peaks[i]
    next_peak <- raw_peaks[i + 1]
    
    left_boundary_index <- which.min(values[prev_peak:curr_peak]) + prev_peak
    right_boundary_index <- which.min(values[curr_peak : next_peak]) + curr_peak
    
    if (i == 2) {left_boundary_index = 0}
    if (i == length(raw_peaks) - 1) {right_boundary_index <- length(values)}
    
    left_values <- values[left_boundary_index : curr_peak]
    right_values <- values[curr_peak : right_boundary_index]
    
    left_peak_height <- values[curr_peak] - min(left_values)
    right_peak_height <- values[curr_peak] - min(right_values)
    
    peak_heights <- c(peak_heights, min(left_peak_height, right_peak_height))
  }
  
  return(peak_heights)
}


find_boundaries <- function(values, peaks, peak_size) {
  boundaries <- list()
  
  for(i in 2: (length(peaks) - 1)) {
    prev_peak <- peaks[i - 1]
    curr_peak <- peaks[i]
    next_peak <- peaks[i + 1]
    
    left_boundary_index <- which.min(values[prev_peak:curr_peak]) + prev_peak
    right_boundary_index <- which.min(values[curr_peak: next_peak]) + curr_peak
    
    if(i == 2) {
      left_boundary_index <- 1
    }
    
    if(i == length(peaks) - 1) {
      right_boundary_index <- length(values)
    }
    
    left_data <- values[left_boundary_index:curr_peak]
    right_data <- values[curr_peak:right_boundary_index]
    

    left_shape <- fit_shape(left_data, "left")
    left_shape_error <- fit_error(left_shape$fit_y, left_shape$y_for_fit)
    
    right_shape <- fit_shape(right_data, "right")
    right_shape_error <-fit_error(right_shape$fit_y, right_shape$y_for_fit)
    
    if(i == 1) {
      sigma <- left_shape$popt[2]
      delta <-  (peak_size * sigma) %>% ceiling() %>% as.integer()
      left_boundary_index <- max(1, curr_peak - delta)
    }
    
    if(i == length(peaks) - 1) {
      sigma <- right_shape$popt[2]
      delta <-  (peak_size * sigma) %>% ceiling() %>% as.integer()
      right_boundary_index <- min(length(values) - 1, curr_peak + delta)
    }
    
    result <- list(
      middleNode = curr_peak,
      leftShape_err = left_shape_error,
      rightShape_err = right_shape_error,
      leftShape = left_shape,
      rightShape = right_shape,
      leftBoundary_fit = left_boundary_index,
      rightBoundary_fit = right_boundary_index,
      leftSigma = left_shape,
      rightSigma = right_shape
    )
    
    boundaries[[i - 1]] <- result
  }
  
  return(boundaries)
}


peak_raw_filters <- function(raw_peaks, values, threshold, peak_minimum_interval) {
  # function to filter raw peaks
  
  n <- length(values)
  raw_peaks <- raw_peaks[values[raw_peaks] > threshold]
  
  raw_peaks_new <- c()
  
  if (length(raw_peaks) > 0) {
    raw_peaks_new <- c(raw_peaks[1])
    
    for (i in 2: length(raw_peaks)) {
      if (raw_peaks[i] - tail(raw_peaks_new, 1) > peak_minimum_interval) {
        raw_peaks_new <- c(raw_peaks_new, raw_peaks[i])
      } else {
        if (values[raw_peaks[i]] > values[tail(raw_peaks_new, 1)]) {
          raw_peaks_new[length(raw_peaks_new)] <- raw_peaks[i]
        }
      }
    }
  }
  
  return(raw_peaks_new)
}


detect_dcl <- function(df, rawdata, chl_peak_min_depth, chl_peak_upper_depth_boundary, config) {
  all_peaks <- detect_peak(df$Fluorescence, config)
  
  chl_peak_upper_depth_boundary_index <- findInterval(chl_peak_upper_depth_boundary, df$Depth)
  
  features <- list(
    allConc = sum(df$Fluorescence),
    allConc_upper = sum(df$Fluorescence[1:chl_peak_upper_depth_boundary_index]),
    peakNums = length(all_peaks)
  )
  
  if(features$peakNums == 0) {
    features$DCL_exists <- 0
    return(features)
  } else {
    chl_peak_depth <- df$Depth[all_peaks$peakIndex]
    
    # if all peaks are above minimum depth, no DCL is detected
    if(!is.null(chl_peak_min_depth) && all(chl_peak_depth < chl_peak_min_depth)) {
      features$DCL_exists <- 0
      return(features)
    } 
    
    if (is.null(chl_peak_min_depth)) {
      dcl_idx <- which.max(df$Fluorescence[all_peaks$peakIndex])
    } else {
      dcl_idx <- which.max(df$Fluorescence[all_peaks$peakIndex] * (chl_peak_depth > chl_peak_min_depth))
    }
    
    features$DCL_exists <- 1
    features$DCL_depth <- chl_peak_depth[dcl_idx]
    features$DCL_conc <- get_raw_concentration(features$DCL_depth, rawdata)
    features$DCL_bottomDepth_fit <- df$Depth[all_peaks$rightIndex_fit[dcl_idx]]
    features$DCL_upperDepth_fit <- df$Depth[all_peaks$leftIndex_fit[dcl_idx]]
    
    if (features$DCL_upperDepth_fit < chl_peak_upper_depth_boundary) {
      features$DCL_upperDepth_fit <- chl_peak_upper_depth_boundary
    }
  }
  return(features)
}

get_raw_concentration <- function(depth, rawdata) {
  tmp <- (rawdata$Depth < depth + 0.5) & (rawdata$Depth < depth - 0.5)
  return(max(rawdata$Fluorescence[tmp]))
}

plot_dcl <- function(features, cleandata) {
  ggplot(data=cleandata) + geom_point(aes(x = Fluorescence, y = -Depth)) + 
    geom_hline(aes(yintercept = -features$DCL_depth)) +
    geom_hline(aes(yintercept = -features$DCL_upperDepth_fit)) + 
    geom_hline(aes(yintercept = -features$DCL_bottomDepth_fit))
}
