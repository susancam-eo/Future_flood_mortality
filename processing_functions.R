
# title: "Functions to process input data"
# author: "Susan Campbell"
# date: "2024-06-21"


# Return periods and flood simulations ---------------------------------------

# function that calculates slope of change in return period over time
calculate_rp_slope <- function(start_rp = initial_return_period, 
                               projected_rp, time_window = 100) {
  if (!is.numeric(projected_rp)) {
    # if cell value is NA, return time series filled with NA values
    return(NA)
  } else {
    # calculate slope for linear model of change in return period over time
    slope <- (projected_rp - start_rp) / (projection_year - start_year)
    return(slope)
  }
}

# function that calculates new rp for a given year and generates binomial draw
calculate_year_specific_rp <- function(calculation_year, cell_rp_slope, 
                                       start_rp = initial_return_period) {
  if (!is.numeric(cell_rp_slope)) {
    return(NA)
  } else {
    # linear interpolation / extrapolation of return period
    year_specific_rp <- cell_rp_slope * (calculation_year - start_year) + start_rp
    # correct values less than 1 (not considering multiple floods in one year in one cell)
    year_specific_rp <- pmax(year_specific_rp, 1)
    
    flood_draw <- rbinom(1, 1, prob = (1/year_specific_rp))
    return(flood_draw)
  }
}


# Aggregate values over time series ------------------------------------------

sum_over_time <- function(time_series_raster) {
  app(time_series_raster, sum, na.rm = TRUE)
}

# TODO: fix this!
pop_wght_avg_over_time <- function(time_series_raster, pop_raster) {
  
  aggregated_results <- list()
  
  # loop through each year
  for (year in 1:nlyr(time_series_raster)) {
    year_data <- time_series_raster[[year]]
    
    # calculate mean for each cell (taking into account population)
    year_data_mean <- global(year_data, fun = "mean", na.rm = TRUE)
    year_pop_sum <- global(pop_raster, fun = "sum", na.rm = TRUE)
    
    df_year_data_mean <- as.data.frame(year_data_mean, xy = TRUE)
    df_population <- as.data.frame(year_pop_sum, xy = TRUE)
    
    # calculate population-weighted averages
    weighted_year_data <- sum(df_year_data_mean[, -c(1:2)], na.rm = TRUE) / sum(df_population[, -c(1:2)], na.rm = TRUE)
    
    aggregated_results[[year]] <- weighted_year_data
  }
  
  aggregated_results_df <- bind_rows(aggregated_results, .id = "year")
  return(aggregated_results_df)
}


# Hazard and PAF -------------------------------------------------------------

# function to calculate population attributable fraction from relative risk
calculate_paf_from_rr <- function(rr, proportion_year = 1, proportion_cell_exposed = 1) {
  paf <- ( proportion_cell_exposed * (rr - 1) ) / (proportion_cell_exposed * rr)
  paf <- paf * proportion_year
  return(paf)
}

# note: not used as of June 21st
generate_rr_draws <- function(rr_est, rr_lower, rr_upper, ndraws = 100) {
  se <- (rr_upper - rr_lower) / (2 * 1.96)
  rr_draws <- rnorm(ndraws, mean = rr_est, sd = se)
}


# Plotting maps and histograms -----------------------------------------------

## Population ----------------------------

plot_pop_map <- function(pop_data, title = "Population in grid cells") {
  # prep data
  pop_data <- as.data.frame(pop_data, xy = TRUE)
  names(pop_data) <- c("x", "y", "layer")
  
  # define the color breaks and corresponding colors
  breaks <- c(0, 1, 1000, 1E4, 1E5, 1E6, max(pop_data$layer, na.rm = TRUE))
  colors <- c("white", "#000000", "#2F2F2F", "#5D5C5C", "#8D8B8B", "#C1C0C0")
  
  # plot
  fig <- ggplot(pop_data, aes(x = x, y = y, fill = layer)) +
    geom_tile() +
    scale_fill_gradientn(
      name = "Population",
      colors = colors,
      values = scales::rescale(breaks),
      na.value = "white") +
    labs(title = title,
         x = NULL,
         y = NULL) +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  return(fig)
}

plot_pop_hist <- function(pop_data, title = "Population in grid cells") {
  # prep data
  pop_data <- as.data.frame(pop_data, xy = TRUE)
  names(pop_data) <- c("x", "y", "layer")
  
  # define custom breakpoints
  breaks <- c(0, 1, 10, 100, 1000, 10000, Inf)
  labels <- c("0", "1-10", "10-100", "100-1,000", "1,000-10,000", ">10,000")
  
  # create bins
  pop_data$bin <- cut(pop_data$layer, breaks = breaks, labels = labels, right = FALSE)
  
  # plot the histogram with log-transformed y-axis
  fig <- ggplot(pop_data, aes(x = bin)) +
    geom_histogram(stat = "count", fill = "#5D5C5C", 
                   color = "black", binwidth = 1) +
    scale_y_log10() +
    labs(title = title,
         x = "Number of people in cell",
         y = "Log of number of global cells") +
    theme_minimal()
  return(fig)
}


## Flood occurrence ----------------------------

plot_number_floods_map <- function(occurrence_layer, year, title = "Flood occurrence") {
  # prep data
  occurrence_layer <- as.data.frame(occurrence_layer[[year]], xy = TRUE)
  names(occurrence_layer) <- c("x", "y", "layer")
  occurrence_layer$layer <- as.factor(occurrence_layer$layer)
  
  # define the color breaks and corresponding colors
  # breaks <- c(0, 1, Inf)
  colors <- c("white","#6BAED6")
  
  fig <- ggplot(occurrence_layer, aes(x = x, y = y, fill = layer)) +
    geom_tile() +
    scale_fill_manual(
      name = "Floods",
      values = colors,
      drop = FALSE,
      na.value = "white"
    ) +
    labs(title = title,
         x = NULL,
         y = NULL) +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  return(fig)
}

plot_number_floods_hist <- function(df, year, title = "Flood occurrence") {
  values_df <- values(df[[year]]) %>%
    as.data.frame() %>%
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value")
  
  # define custom breakpoints
  breaks <- c(0, 1, Inf)
  labels <- c("0", "1")
  
  # create bins
  values_df$bin <- cut(values_df$value, breaks = breaks, labels = labels, right = FALSE)
  
  # plot the histogram with log-transformed y-axis
  fig <- ggplot(values_df, aes(x = bin)) +
    geom_histogram(stat = "count", fill = "#6BAED6", 
                   color = "black", binwidth = 1) +
    scale_y_log10() +
    labs(title = title,
         x = "Cells with and without floods",
         y = "Log of number of global cells") +
    theme_minimal()
  return(fig)
}

plot_occurrence_over_time_map <- function(input_raster, 
                                          title = "Sum of flood occurrence over time",
                                          fixed_binned_scales = FALSE) {
  # prep data
  input_raster <- as.data.frame(input_raster, xy = TRUE)
  names(input_raster) <- c("x", "y", "layer")
  
  # define the color breaks and corresponding colors
  breaks <- c(0, 1, 10, 30, 100)
  #colors <- c("white","#02397F", "#2458C1", "#2D83F1")
  colors <- brewer.pal(4, "Blues")
  labels <- c("0", "1-10", "10-30", "30-100", ">100")
  
  # plot version with continuous color definitions (poor legend)
  fig_continuous_freescale <- ggplot(input_raster, aes(x = x, y = y, fill = layer)) +
    geom_tile() +
    scale_fill_gradientn(
      name = "Floods",
      colors = colors,
      values = scales::rescale(breaks),
      na.value = "white") +
    labs(title = title,
         x = NULL,
         y = NULL) +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  input_raster$layer_binned <- cut(input_raster$layer, 
                                   breaks = breaks, 
                                   include.lowest = TRUE, 
                                   labels = labels[1:length(labels)-1])
  
  # plot version with binned color definitions (good legend)
  fig_binned_comparable <- ggplot(input_raster, aes(x = x, y = y, fill = layer_binned)) +
    geom_tile() +
    scale_fill_manual(
      name = "Floods",
      values = colors,
      drop = FALSE,
      na.value = "white"
    ) +
    labs(title = title,
         x = NULL,
         y = NULL) +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  if(fixed_binned_scales == TRUE) {
    return(fig_binned_comparable)
  } else {
    return(fig_continuous_freescale)
  }
}


## Flood exposures -----------------------------

plot_number_exposed_map <- function(input_raster, year, title = "People exposed to floods",
                                    fixed_binned_scales = FALSE) {
  # prep data
  input_raster <- as.data.frame(input_raster[[year]], xy = TRUE)
  names(input_raster) <- c("x", "y", "layer")
  
  # define the color breaks and corresponding colors
  # breaks <- c(0, 1, 1000, 10000, max(exposed_layer$layer, na.rm = TRUE))
  # colors <- c("white", "#390666", "#7229B3", "#AA51FA")
  breaks <- c(0, 1, 1E3, 1E4, 1E6)
  colors <- brewer.pal(5, "Purples")
  labels <- c("0", "1 - 1K", "1K - 10K", "10K - 1M", "> 1M")
  
  # plot version with continuous color definitions (poor legend)
  fig_continuous_freescale <- ggplot(input_raster, aes(x = x, y = y, fill = layer)) +
    geom_tile() +
    scale_fill_gradientn(
      name = "Exposed \npopulation",
      colors = colors,
      values = scales::rescale(breaks),
      na.value = "white") +
    labs(title = title,
         x = NULL,
         y = NULL) +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  input_raster$layer_binned <- cut(input_raster$layer, 
                                   breaks = breaks, 
                                   include.lowest = TRUE, 
                                   labels = labels[1:length(labels)-1])
  
  # plot version with binned color definitions (good legend)
  fig_binned_comparable <- ggplot(input_raster, aes(x = x, y = y, fill = layer_binned)) +
    geom_tile() +
    scale_fill_manual(
      name = "Exposed \npopulation",
      values = colors,
      drop = FALSE,
      na.value = "white"
    ) +
    labs(title = title,
         x = NULL,
         y = NULL) +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  if(fixed_binned_scales == TRUE) {
    return(fig_binned_comparable)
  } else {
    return(fig_continuous_freescale)
  }
}

plot_number_exposed_hist <- function(df, year, title = "Exposure to floods") {
  values_df <- values(df[[year]]) %>%
    as.data.frame() %>%
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value")
  
  # define custom breakpoints
  breaks <- c(0, 1, 10, 100, 1E3, 1E4, 1E6, Inf)
  labels <- c("0", "1-10", "10-100", "100-1,000", "1,000-10,000", "10,000-1,000,000", ">1,000,000")
  
  # create bins
  values_df$bin <- cut(values_df$value, breaks = breaks, labels = labels, right = FALSE)
  
  # plot the histogram with log-transformed y-axis
  fig <- ggplot(values_df, aes(x = bin)) +
    geom_histogram(stat = "count", fill = "#756BB1", 
                   color = "black", binwidth = 1) +
    scale_y_log10() +
    labs(title = title,
         x = "Number of people exposed to floods",
         y = "Log of number of global cells") +
    theme_minimal()
  return(fig)
}

plot_exposure_over_time_map <- function(input_raster, 
                                        title = "Sum of people exposed to floods over time",
                                        fixed_binned_scales = FALSE) {
  
  input_raster <- as.data.frame(input_raster, xy = TRUE)
  names(input_raster) <- c("x", "y", "layer")
  
  # define the color breaks and corresponding colors
  # breaks <- c(0, 10000, 100000, 500000, max(input_raster$layer, na.rm = TRUE))
  # colors <- c("white", "#390666", "#7229B3", "#AA51FA")
  breaks <- c(0, 1, 1E3, 1E4, 1E6)
  colors <- brewer.pal(length(breaks) - 1, "Purples")
  labels <- c("0 - 1", "1 - 1K", "1K - 10K", "10K - 1M", "> 1M")
  
  # plot version with continuous color definitions (poor legend)
  fig_continuous_freescale <- ggplot(input_raster, aes(x = x, y = y, fill = layer)) +
    geom_tile() +
    scale_fill_gradientn(
      name = "Exposed \npopulation",
      colors = colors,
      breaks = breaks,
      labels = labels,
      values = scales::rescale(breaks),
      na.value = "white") +
    labs(title = title,
         x = NULL,
         y = NULL) +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  input_raster$layer_binned <- cut(input_raster$layer, 
                                   breaks = breaks, 
                                   include.lowest = TRUE, 
                                   labels = labels[1:length(labels)-1])
  
  # plot version with binned color definitions (good legend)
  fig_binned_comparable <- ggplot(input_raster, aes(x = x, y = y, fill = layer_binned)) +
    geom_tile() +
    scale_fill_manual(
      name = "Exposed \npopulation",
      values = colors,
      drop = FALSE,
      na.value = "white"
    ) +
    labs(title = title,
         x = NULL,
         y = NULL) +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  if(fixed_binned_scales == TRUE) {
    return(fig_binned_comparable)
  } else {
    return(fig_continuous_freescale)
  }
}


## Excess deaths -----------------------------

plot_excess_deaths_map_from_raster <- function(input_raster, year, title = "Excess all-cause deaths",
                                   fixed_binned_scales = FALSE) {
  # prep data
  input_raster <- as.data.frame(input_raster[[year]], xy = TRUE)
  names(input_raster) <- c("x", "y", "layer")
  
  # define the color breaks and corresponding colors
  # breaks <- c(0, 1, 1000, 10000, max(exposed_layer$layer, na.rm = TRUE))
  # colors <- c("white", "#390666", "#7229B3", "#AA51FA")
  breaks <- c(0, 1, 1E2, 1E3, 1E4, 1E6)
  colors <- brewer.pal(6, "Reds")
  labels <- c("0", "1 - 100", "100 - 1K", "1K - 10K", "10K - 1M", "> 1M")
  
  # plot version with continuous color definitions (poor legend)
  fig_continuous_freescale <- ggplot(input_raster, aes(x = x, y = y, fill = layer)) +
    geom_tile() +
    scale_fill_gradientn(
      name = "Excess \ndeaths",
      colors = colors,
      values = scales::rescale(breaks),
      na.value = "white") +
    labs(title = title,
         x = NULL,
         y = NULL) +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  input_raster$layer_binned <- cut(input_raster$layer, 
                                   breaks = breaks, 
                                   include.lowest = TRUE, 
                                   labels = labels[1:length(labels)-1])
  
  # plot version with binned color definitions (good legend)
  fig_binned_comparable <- ggplot(input_raster, aes(x = x, y = y, fill = layer_binned)) +
    geom_tile() +
    scale_fill_manual(
      name = "Excess \ndeaths",
      values = colors,
      drop = FALSE,
      na.value = "white"
    ) +
    labs(title = title,
         x = NULL,
         y = NULL) +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  if(fixed_binned_scales == TRUE) {
    return(fig_binned_comparable)
  } else {
    return(fig_continuous_freescale)
  }
}

plot_excess_deaths_map <- function(df, year, title = "Flood-attributable all-cause deaths", 
                                   fixed_binned_scales = FALSE) {
  year <- year + start_year - 1
  input_df <- df %>% filter(year == !!year)
  
  breaks <- c(0, 1, 10, 1E2, 1E3, 1E4, 1E5)
  colors <- brewer.pal(length(breaks) - 1, "Reds")
  labels <- c("0 - 1", "1 - 10", "10 - 100", "100 - 1K", "1K - 10K", "10K - 100K", ">100K")
  
  fig_continuous_freescale <- ggplot(data = input_df) +
    geom_sf(aes(fill = flood_attributable_deaths), color = "black") +
    theme_minimal() +
    scale_fill_gradientn(
      name = "Flood \nattributable \ndeaths",
      colors = colors,
      na.value = "white") +
    labs(title = title,
         x = NULL,
         y = NULL)
  
  input_df$binned_val <- cut(input_df$flood_attributable_deaths,
                             breaks = breaks,
                             include.lowest = TRUE,
                             labels = labels[1:length(labels) - 1])
  
  fig_binned_comparable <- ggplot(data = input_df) +
    geom_sf(aes(fill = binned_val), color = "black") +
    theme_minimal() +
    scale_fill_manual(
      name = "Flood \nattributable \ndeaths",
      values = colors,
      drop = FALSE,
      na.value = "white") +
    labs(title = title,
         x = NULL,
         y = NULL)
  
  if (fixed_binned_scales == TRUE) {
    return(fig_binned_comparable)
  } else {
    return(fig_continuous_freescale)
  }
}

plot_number_deaths_hist <- function(df, year, title = "Excess deaths") {
  values_df <- values(df[[year]]) %>%
    as.data.frame() %>%
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value")
  
  # define custom breakpoints
  breaks <- c(0, 1, 10, 100, 1000, 10000, Inf)
  labels <- c("0", "1-10", "10-100", "100-1,000", "1,000-10,000", ">10,000")
  
  # create bins
  values_df$bin <- cut(values_df$value, breaks = breaks, labels = labels, right = FALSE)
  
  # plot the histogram with log-transformed y-axis
  fig <- ggplot(values_df, aes(x = bin)) +
    geom_histogram(stat = "count", fill = "#B51B13", 
                   color = "black", binwidth = 1) +
    scale_y_log10() +
    labs(title = title,
         x = "Number of excess deaths",
         y = "Log of number of global cells") +
    theme_minimal()
  return(fig)
}

plot_deaths_over_time_map_from_raster <- function(input_raster, 
                                      title = "Sum of deaths due to floods over time",
                                      fixed_binned_scales = FALSE) {
  # prep data
  input_raster <- as.data.frame(input_raster, xy = TRUE)
  names(input_raster) <- c("x", "y", "layer")
  
  # define the color breaks and corresponding colors
  breaks <- c(0, 1, 1E2, 1E3, 1E5, 1E7)
  colors <- brewer.pal(6, "Reds")
  labels <- c("0", "1 - 100", "100 - 1K", "1K - 100K", "100K - 10M", "> 10M")
  
  # plot version with continuous color definitions (poor legend)
  fig_continuous_freescale <- ggplot(input_raster, aes(x = x, y = y, fill = layer)) +
    geom_tile() +
    scale_fill_gradientn(
      name = "Excess \ndeaths",
      colors = colors,
      breaks = breaks,
      labels = labels,
      values = scales::rescale(breaks),
      na.value = "white") +
    labs(title = title,
         x = NULL,
         y = NULL) +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  input_raster$layer_binned <- cut(input_raster$layer, 
                                   breaks = breaks, 
                                   include.lowest = TRUE, 
                                   labels = labels[1:length(labels)-1])
  
  # plot version with binned color definitions (good legend)
  fig_binned_comparable <- ggplot(input_raster, aes(x = x, y = y, fill = layer_binned)) +
    geom_tile() +
    scale_fill_manual(
      name = "Excess \ndeaths",
      values = colors,
      drop = FALSE,
      na.value = "white"
    ) +
    labs(title = title,
         x = NULL,
         y = NULL) +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  if(fixed_binned_scales == TRUE) {
    return(fig_binned_comparable)
  } else {
    return(fig_continuous_freescale)
  }
}

plot_deaths_over_time_map <- function(df, title = "Sum of deaths due to floods over time", 
                                      fixed_binned_scales = FALSE) {
  df <- st_as_sf(df)
  # define the color breaks and corresponding colors
  breaks <- c(0, 1, 1E2, 1E3, 1E5, 1E7)
  colors <- brewer.pal(6, "Reds")
  labels <- c("0", "1 - 100", "100 - 1K", "1K - 100K", "100K - 10M", "> 10M")
  
  # plot version with continuous color definitions (poor legend)
  fig_continuous_freescale <- ggplot(data = df) +
    geom_sf(aes(fill = total_deaths), color = "black") +
    theme_minimal() +
    scale_fill_gradientn(
      name = "Flood \nattributable \ndeaths",
      colors = colors,
      na.value = "white") +
    labs(title = title,
         x = NULL,
         y = NULL)
  
  df$binned_val <- cut(df$total_deaths, 
                       breaks = breaks, 
                       include.lowest = TRUE, 
                       labels = labels[1:length(labels) - 1])
  
  # plot version with binned color definitions (good legend)
  fig_binned_comparable <- ggplot(data = df) +
    geom_sf(aes(fill = binned_val), color = "black") +
    theme_minimal() +
    scale_fill_manual(
      name = "Flood \nattributable \ndeaths",
      values = colors,
      drop = FALSE,
      na.value = "white") +
    labs(title = title,
         x = NULL,
         y = NULL)
  
  if (fixed_binned_scales) {
    return(fig_binned_comparable)
  } else {
    return(fig_continuous_freescale)
  }
}


## Excess death rate -----------------------------

plot_death_rate_map_from_raster <- function(input_raster, year, 
                                title = "Excess all-cause death rate") {
  # prep data
  input_raster <- as.data.frame(input_raster[[year]], xy = TRUE)
  names(input_raster) <- c("x", "y", "layer")
  
  # define the color breaks and corresponding colors
  breaks <- c(0, 0.001, 0.01, 0.1, 1)
  colors <- brewer.pal(5, "Oranges")
  #colors <- c("white", "#6E3706", "#D3721D", "#FD963C")
  
  # plot
  fig <- ggplot(input_raster, aes(x = x, y = y, fill = layer)) +
    geom_tile() +
    scale_fill_gradientn(
      name = "Excess\ndeath\nrate",
      colors = colors,
      values = scales::rescale(breaks),
      na.value = "white") +
    labs(title = title,
         x = NULL,
         y = NULL) +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  return(fig)
}

plot_death_rate_hist <- function(df, year, title = "Excess all-cause death rate") {
  values_df <- values(df[[year]]) %>%
    as.data.frame() %>%
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value")
  
  # define custom breakpoints
  breaks <- c(0, 0.001, 0.01, 0.1, 1)
  labels <- c("0-0.001", "0.001-0.01", "0.01-0.1", "0.1-1")
  
  # create bins
  values_df$bin <- cut(values_df$value, breaks = breaks, labels = labels, right = FALSE)
  
  # plot the histogram with log-transformed y-axis
  fig <- ggplot(values_df, aes(x = bin)) +
    geom_histogram(stat = "count", fill = "orange", 
                   color = "black", binwidth = 1) +
    scale_y_log10() +
    labs(title = title,
         x = "Excess death rate",
         y = "Log of number of global cells") +
    theme_minimal()
  return(fig)
}

# plot death
plot_excess_death_rate_map <- function(df, year, title = "Flood-attributable all-cause death rate", 
                                       fixed_binned_scales = FALSE) {
  year <- year + start_year - 1
  input_df <- df %>% filter(year == !!year)

  breaks <- c(0, 1E-6, 1E-5, 1E-4, 1E-3)
  colors <- brewer.pal(length(breaks) - 1, "Oranges")
  labels <- c("0 - 0.1", "0.1 - 1", "1 - 10", "10 - 100", "> 100")
  
  fig_continuous_freescale <- ggplot(data = input_df) +
    geom_sf(aes(fill = flood_attributable_death_rate), color = "black") +
    theme_minimal() +
    scale_fill_gradientn(
      name = "Flood \nattributable \ndeath rate \nper 100K",
      colors = colors,
      na.value = "white") +
    labs(title = title,
         x = NULL,
         y = NULL)
  
  input_df$binned_val <- cut(input_df$flood_attributable_death_rate,
                             breaks = breaks,
                             include.lowest = TRUE,
                             labels = labels[1:length(labels) - 1])
  
  fig_binned_comparable <- ggplot(data = input_df) +
    geom_sf(aes(fill = binned_val), color = "black") +
    theme_minimal() +
    scale_fill_manual(
      name = "Flood \nattributable \ndeath rate\n per 100K",
      values = colors,
      drop = FALSE,
      na.value = "white") +
    labs(title = title,
         x = NULL,
         y = NULL)
  
  if (fixed_binned_scales == TRUE) {
    return(fig_binned_comparable)
  } else {
    return(fig_continuous_freescale)
  }
}


# Aggregate to country and region --------------------------------------------

aggregate_by_country <- function(raster_stack, boundary_shp, 
                                 fun = sum, na.rm = TRUE, method = 'simple',
                                 first_year = start_year) {
  aggregated_list <- list()
  # for each year, aggregate from grid level to country level
  for (year in 1:nlyr(raster_stack)) {
    year_name <- first_year + year - 1
    aggregated <- extract(raster_stack[[year]], boundary_shp, 
                          fun = fun, na.rm = na.rm, method = method)
    aggregated_df <- as.data.frame(aggregated)
    names(aggregated_df) <- c("country", year_name)
    aggregated_list[[year]] <- aggregated_df
  }
  aggregated_df <- reduce(aggregated_list, full_join, by = "country")
  return(aggregated_df)
}


# Bootstrapping and processing population and mortality projections ---------------------------

# function to extrapolate population and mortality datausing ARIMA
extrapolate_future_years_with_arima <- function(df, value_col, end_year) {
  df_list <- df %>% group_split(loc_id)
  
  result_list <- df_list %>% map(function(sub_df) {
    loc_id <- unique(sub_df$loc_id)
    
    ts_data <- ts(sub_df[[value_col]], start = min(sub_df$year), frequency = 1)
    
    # apply log transformation to the time series data to avoid forecasting negative values
    log_ts_data <- log(ts_data + 1)  # +1 to handle zeros

    arima_model <- auto.arima(log_ts_data)
    
    n_ahead <- end_year - max(sub_df$year)
    forecast_values <- forecast(arima_model, h = n_ahead)
    
    # back-transform the predictions to the original scale
    predicted_values <- exp(forecast_values$mean) - 1  # Subtract 1 to undo the +1 during transformation
    # ensure that no values are below zero
    predicted_values[predicted_values < 0] <- 0
    
    future_years <- seq(max(sub_df$year) + 1, end_year)
    future_df <- data.frame(
      loc_id = loc_id,
      year = future_years,
      value = predicted_values
    )
    colnames(future_df)[3] <- value_col
    bind_rows(sub_df, future_df)
  })
  
  result_df <- bind_rows(result_list)
  return(result_df)
}

# function to extrapolate population and mortality projections using ARIMA
extend_forecasts <- function(df, ages, locs, end_year) {
  
  df <- df %>%
    filter(age_group_name == ages, location_id %in% locs) %>%
    dplyr::select(sex, location_id, location_name, year_id, val, upper, lower) %>%
    group_by(location_id, year_id) %>%
    summarise(val = sum(val, na.rm = TRUE),
              upper = sum(upper, na.rm = TRUE),
              lower = sum(lower, na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::select(location_id, year_id, val, lower, upper) %>%
    distinct()
  
  names(df) <- c("loc_id", "year", "val", "lower", "upper")
  
  # Apply ARIMA extrapolation for each column
  val_df <- extrapolate_future_years_with_arima(df, "val", end_year) %>%
    dplyr::select(c(loc_id, year, val))
  lower_df <- extrapolate_future_years_with_arima(df, "lower", end_year) %>%
    dplyr::select(c(loc_id, year, lower))
  upper_df <- extrapolate_future_years_with_arima(df, "upper", end_year) %>%
    dplyr::select(c(loc_id, year, upper))
  
  # Combine the results
  result <- val_df %>%
    left_join(lower_df, by = c("loc_id", "year")) %>%
    left_join(upper_df, by = c("loc_id", "year"))
  
  return(result)
}

# generate a draw given an estimate (50th) and bounds (2.5th and 97.5th percentile)
generate_draw <- function(est, est_lower, est_upper, ndraws = num_draws) {
  se <- (est_upper - est_lower) / (2 * 1.96)
  vals <- rnorm(ndraws, mean = est, sd = se)
}

# function to generate many draws 
execute_all_draws <- function(df, prefix) {
  
  df <- df %>%mutate(draws = pmap(list(val, lower, upper), 
                                  ~ generate_draw(..1, ..2, ..3, num_draws))) %>%
    unnest(draws) %>%
    group_by(loc_id, year) %>%
    mutate(draw_id = row_number()) %>%
    pivot_wider(names_from = draw_id, values_from = draws, names_prefix = prefix)
  
  return(df)
}


# Process flood exposures and counts --------------------------------------

# calculate the fraction of a country that has been exposed to floods for each year
aggregate_exposure_to_country <- function(exposed_stack, boundary_shp, country_pops) {
  
  exposed_list <- list()
  
  # loop through each year
  for (year_val in 1:nlyr(exposed_stack)) {
    
    # aggregate population exposed to country level
    country_exposed <- extract(exposed_stack[[year_val]], 
                               boundary_shp[, c("name", "loc_id", "geometry")], 
                               fun = sum, na.rm = TRUE, method = 'simple')
    country_exposed <- as.data.frame(country_exposed)
    names(country_exposed) <- c("ID", "exposed")
    
    country_exposed$year <- year_val + start_year - 1
    
    country_exposed <- left_join(country_id_names_dict, country_exposed, 
                                 join_by(index == ID)) %>%
      left_join(country_pops, join_by(index == ID)) %>%
      dplyr::select(-index, -code)
    
    exposed_list[[year_val]] <- country_exposed
  }
  
  country_exposed_df <- bind_rows(exposed_list)
  
  country_exposed_df <- country_exposed_df %>%
    mutate(fraction_exposed = exposed / worldpop)
  
  return(country_exposed_df)
}

# count floods in a grid cell across time/stacks
count_floods <- function(stack, start_index, end_index) {
  flood_counts <- sum(stack[[start_index:end_index]], na.rm = TRUE)
  return(flood_counts)
}

# count floods within country boundaries
extract_flood_counts <- function(flood_counts, boundaries) {
  
  id_map <- boundaries %>%
    as.data.frame() %>%
    dplyr::select(c("loc_id", "name")) %>%
    mutate(id = row_number())
  
  flood_extracted <- terra::extract(flood_counts, boundaries, cells = TRUE)
  flood_extracted <- data.frame(id = flood_extracted$ID,
                                flood_counts = flood_extracted[ , 2]) %>%
    left_join(id_map, join_by("id")) %>%
    dplyr::select(-id)
  
  return(flood_extracted)
}

# summarize the proportions of grid cells with the listed flood counts within each country
summarize_flood_counts <- function(flood_extracted, boundaries, region_dict) {
  
  flood_extracted <- flood_extracted %>%
    left_join(st_as_sf(boundaries), join_by("loc_id", "name")) %>%
    left_join(region_dict, join_by("loc_id", "name"))
  
  flood_summary <- flood_extracted %>%
    group_by(region_name) %>%
    summarize(
      `none` = sum(is.na(flood_counts)) / n(),
      `0` = sum(flood_counts == 0, na.rm = TRUE) / n(),
      `1-2` = sum(flood_counts >= 1 & flood_counts <= 2, na.rm = TRUE) / n(),
      `3-4` = sum(flood_counts >= 3 & flood_counts <= 4, na.rm = TRUE) / n(),
      `5-7` = sum(flood_counts >= 5 & flood_counts <= 7, na.rm = TRUE) / n(),
      `8-10` = sum(flood_counts >= 8 & flood_counts <= 10, na.rm = TRUE) / n(),
      `11-15` = sum(flood_counts >= 11 & flood_counts <= 15, na.rm = TRUE) / n(),
      `>15` = sum(flood_counts > 16, na.rm = TRUE) / n()
    )
  
  return(flood_summary)
}


# Calculate decade trends for countries -------------------------------------

# calculate the decade mean of population, number of floods, exposed people, 
  # deaths, and death rate, for all countries
summarize_decades <- function(yearly_totals, fun = mean, region_map = NULL, uncertainty = FALSE) {
  
  yearly_totals$decade <- (yearly_totals$year) - (yearly_totals$year %% 10)
  
  if(uncertainty == FALSE) {
    decade_summary <- yearly_totals %>%
      group_by(code, name, loc_id, decade) %>%
      summarize(
        population = fun(population, na.rm = TRUE),
        fraction_exposed = fun(fraction_exposed, na.rm = TRUE),
        #number_exposed = fun(exposed),
        flood_attributable_deaths = fun(flood_attributable_deaths, na.rm = TRUE),
        flood_attributable_death_rate = fun(flood_attributable_death_rate, na.rm = TRUE),
        .groups = 'drop')
    
  } else {
    decade_summary <- yearly_totals %>%
      group_by(name, loc_id, decade) %>%
      summarize(
        fraction_exposed = fun(fraction_exposed, na.rm = TRUE),
        #number_exposed = fun(exposed),
        flood_attributable_deaths = fun(flood_attributable_deaths, na.rm = TRUE),
        flood_attributable_deaths_lower = fun(flood_attributable_deaths_lower, na.rm = TRUE),
        flood_attributable_deaths_upper = fun(flood_attributable_deaths_upper, na.rm = TRUE),
        flood_attributable_death_rate = fun(flood_attributable_death_rate, na.rm = TRUE),
        flood_attributable_death_rate_lower = fun(flood_attributable_death_rate_lower, na.rm = TRUE),
        flood_attributable_death_rate_upper = fun(flood_attributable_death_rate_upper, na.rm = TRUE),
        .groups = 'drop')
  }
  
  if(!is.null(region_map)) {
    decade_summary <- left_join(decade_summary, region_map, by = c("name", "loc_id"))
  }
  return(decade_summary)
}

# calculate percent change in decade mean for each value 
relative_decade_change <- function(decade_summary, region_map = NULL, uncertainty = FALSE) {
  
  if(uncertainty == FALSE) {
    result <- decade_summary %>%
      group_by(code, name, loc_id) %>%
      arrange(decade) %>% 
      mutate(
        population_change = (population - first(population)) / first(population),
        #floods_change = (floods - first(floods)) / first(floods),
        fraction_exposed_change = (fraction_exposed - first(fraction_exposed)) / first(fraction_exposed),
        #number_exposed_change = (number_exposed - first(number_exposed)) / first(number_exposed),
        flood_attributable_deaths_change = (flood_attributable_deaths - first(flood_attributable_deaths)) / first(flood_attributable_deaths),
        flood_attributable_death_rate_change = (flood_attributable_death_rate - first(flood_attributable_death_rate)) / first(flood_attributable_death_rate)
      ) %>%
      ungroup()
    
  } else {
    result <- decade_summary %>%
      group_by(name, loc_id) %>%
      arrange(decade) %>% 
      mutate(
        fraction_exposed_change = (fraction_exposed - first(fraction_exposed)) / first(fraction_exposed),
        #number_exposed_change = (number_exposed - first(number_exposed)) / first(number_exposed),
        flood_attributable_deaths_change = 
          (flood_attributable_deaths - first(flood_attributable_deaths)) / first(flood_attributable_deaths),
        flood_attributable_deaths_lower_change = 
          (flood_attributable_deaths_lower - first(flood_attributable_deaths_lower)) / first(flood_attributable_deaths_lower),
        flood_attributable_deaths_upper_change = 
          (flood_attributable_deaths_upper - first(flood_attributable_deaths_upper)) / first(flood_attributable_deaths_upper),
        flood_attributable_death_rate_change = 
          (flood_attributable_death_rate - first(flood_attributable_death_rate)) / first(flood_attributable_death_rate),
        flood_attributable_death_rate_lower_change = 
          (flood_attributable_death_rate_lower - first(flood_attributable_death_rate_lower)) / first(flood_attributable_death_rate_lower),
        flood_attributable_death_rate_upper_change = 
          (flood_attributable_death_rate_upper - first(flood_attributable_death_rate_upper)) / first(flood_attributable_death_rate_upper)
      ) %>%
      ungroup()
  }

  if(!is.null(region_map)) {
    decade_summary <- left_join(decade_summary, region_map, by = "name")
  }
  return(result)
}


# Summarize trends -------------------------------------

# aggregate a given metric (take the mean), e.g. for death rate, to a higher geographical level
summarize_by_geography <- function(df, variable, summary_level, uncertainty = FALSE) {
  
  if(uncertainty == FALSE) {
    df <- df %>%
      dplyr::select(!!sym(summary_level), decade, !!sym(variable)) %>%
      group_by(!!sym(summary_level), decade) %>%
      summarize(avg_var = mean(!!sym(variable), na.rm = TRUE), .groups = 'drop') %>%
      pivot_wider(names_from = decade, values_from = avg_var, names_prefix = "decade_") %>%
      mutate(percent_change = ((decade_2120 - decade_2020) / decade_2020) * 100)
  } else {
    df <- df %>%
      dplyr::select(!!sym(summary_level), decade, !!sym(variable),
                    !!sym(paste0(variable, "_lower")), !!sym(paste0(variable, "_upper"))) %>%
      group_by(!!sym(summary_level), decade) %>%
      summarize(avg_var = mean(!!sym(variable), na.rm = TRUE),
                avg_var_lower = mean(!!sym(paste0(variable, "_lower")), na.rm = TRUE),
                avg_var_upper = mean(!!sym(paste0(variable, "_upper")), na.rm = TRUE),
                .groups = 'drop') %>%
      pivot_wider(names_from = decade, values_from = avg_var, names_prefix = "decade_")
  }
  return(df)
}

# calculate number of deaths attributable to flooding
calculate_deaths <- function(ssp_df) {
  ssp_df <- ssp_df %>%
    mutate(flood_attributable_deaths_corrected = 10 * flood_attributable_deaths) %>%
    dplyr::select(name, loc_id, decade, flood_attributable_deaths_corrected, 
                  region_name, superregion_name)
  return(ssp_df)
}


# Plotting trends at country, region, and superregion level -----------------

# makes a heatmap of a given metric (e.g., deaths) for each decade, and combines plots for three climate scenarios 
plot_heatmap_changes <- function(df_best, df_mid, df_worst, 
                                 variable, title) {
  
  breaks <- c(-1, -0.5, -0.25, 0, 0.1, 0.25, 0.5, 1)
  colors <- brewer.pal(8, "Oranges")
  labels <- c("50-100% lower", "25-50% lower", "0-25% lower", "0-10% higher",
              "10-25% higher", "25-50% higher", "50-100% higher", ">100% higher")
  
  color_scale_anchor <- df_best %>% 
    filter(decade == 2020) %>% 
    pull(variable) %>% quantile(probs = 0.95, na.rm = TRUE)
  
  make_plot <- function(df, title) {
    
    df$binned_val <- cut(df[[variable]],
                         breaks = breaks,
                         include.lowest = TRUE,
                         labels = labels[1:length(labels) - 1])
    
    ggplot(df, aes(x = decade, y = name, fill = binned_val)) +
      geom_tile() +
      # scale_fill_gradient(low = "white", high = "red", na.value = "gray80") +
      # scale_fill_viridis_c(option = "magma", limits = c(0, color_scale_anchor), na.value = "white") +
      scale_fill_manual(values = colors, drop = FALSE, 
                        na.value = "white", 
                        name = "Flood \nattributable \ndeath rate")
    labs(title = title, 
         x = NULL, 
         y = NULL) +
      theme_minimal()
  }
  
  plot_best <- make_plot(df_best, title = paste(title, "in SSP 126"))
  plot_mid <- make_plot(df_mid, title = paste(title, "in SSP 245"))
  plot_worst <- make_plot(df_worst, title = paste(title, "in SSP 585"))
  
  combined_plot <- plot_grid(plot_best, plot_mid, plot_worst, nrow = 1, ncol = 3)
  return(combined_plot)  
}

# plot maps of three decades for a given metric for all countries, for three climate scenarios (nine maps total)
plot_decade_maps_all_futures <- function(df_best, df_mid, df_worst, 
                                         geom_file, join_col, variable, title) {
  # set decades of interest
  first_decade <- 2020
  mid_decade <- 2070
  last_decade <- 2120
  
  # to make sure that color scheme is the same for all individual plots
  color_scale_anchor <- df_best %>% 
    filter(decade == last_decade) %>% 
    pull(variable) %>% quantile(probs = 0.9, na.rm = TRUE)
  
  # function to filter and plot for a specific decade, with a discrete color scale
  plot_for_decade <- function(df, decade, title_suffix) {
    
    if(variable == "flood_attributable_death_rate") {
      breaks <- c(0, 1E-7, 5E-6, 1E-6, 5E-5, 1E-5, 5E-4, 1E-4)
      colors <- brewer.pal(length(breaks) - 1, "Oranges")
      labels <- c("0 - 0.01", "0.01 - 0.05", "0.05 - 0.1", "0.1 - 0.5",
                  "0.5 - 1", "1 - 5", "5 - 10", ">10")
      legend_name <- "Flood \nattributable \ndeath rate \nper 100K"
      
    } else if(variable == "flood_attributable_deaths") {
      breaks <- c(0, 10, 100, 1E3, 1E4, 1E5)
      colors <- brewer.pal(length(breaks) - 1, "Reds")
      labels <- c("0 - 10", "10-100", "100 - 1K", "1K - 10K", "10K - 100K", "> 100K")
      legend_name <- "Flood \nattributable \ndeaths"
    }
    
    df_decade <- df %>% filter(decade == !!decade)
    df_decade <- df_decade %>%
      mutate(binned_val = cut(get(variable), 
                              breaks = breaks, 
                              labels = labels[1:length(labels) - 1], 
                              include.lowest = TRUE))
    
    merged_df <- left_join(df_decade, geom_file, by = join_col)
    
    ggplot(data = merged_df) +
      geom_sf(aes(fill = binned_val, geometry = geometry), color = "black") +
      scale_fill_manual(values = colors, drop = FALSE, 
                        na.value = "white", name = legend_name) +
      labs(title = paste(title, title_suffix), fill = variable) +
      theme_minimal()
  }
  
  # create three plots for best scenario: SSP 126
  plot_first_decade_best <- plot_for_decade(df_best, first_decade, 
                                            paste("in", first_decade, "decade, SSP 126"))
  plot_middle_decade_best <- plot_for_decade(df_best, mid_decade, 
                                             paste("in", mid_decade, "decade, SSP 126"))
  plot_last_decade_best <- plot_for_decade(df_best, last_decade, 
                                           paste("in", last_decade, "decade, SSP 126"))
  # create three plots for mid scenario: SSP 245
  plot_first_decade_mid <- plot_for_decade(df_mid, first_decade, 
                                           paste("in", first_decade, "decade, SSP 245"))
  plot_middle_decade_mid <- plot_for_decade(df_mid, mid_decade, 
                                            paste("in", mid_decade, "decade, SSP 245"))
  plot_last_decade_mid <- plot_for_decade(df_mid, last_decade, 
                                          paste("in", last_decade, "decade, SSP 245"))
  # create three plots for worst scenario: SSP 585
  plot_first_decade_worst <- plot_for_decade(df_worst, first_decade, 
                                             paste("in", first_decade, "decade, SSP 585"))
  plot_middle_decade_worst <- plot_for_decade(df_worst, mid_decade, 
                                              paste("in", mid_decade, "decade, SSP 585"))
  plot_last_decade_worst <- plot_for_decade(df_worst, last_decade, 
                                            paste("in", last_decade, "decade, SSP 585"))
  
  plot_list <- list(plot_first_decade_best, plot_middle_decade_best, plot_last_decade_best,
                    plot_first_decade_mid, plot_middle_decade_mid, plot_last_decade_mid,
                    plot_first_decade_worst, plot_middle_decade_worst, 
                    plot_last_decade_worst)
  
  # combine the plots into one figure
  combined_plot <- plot_grid(plotlist = plot_list, nrow = 3, ncol = 3)
  return(combined_plot)
}

# make lineplots of trends in metric over time for three climate scenarios (three lineplots)
plot_decade_trends_for_regions <- function(ssp_126, ssp_245, ssp_585, variable) {
  
  # ymax <- ssp_585 %>% 
  #   filter(decade == 2120) %>% 
  #   pull(gsub("region", "flood_attributable", variable)) %>% quantile(probs = 1, na.rm = TRUE)
  # ymax <- ymax + 1
  
  plot_superregion_changes <- function(df, scenario) {
    
    df$scenario <- scenario
    
    df <- df %>% 
      dplyr::select(c("scenario", "name", "loc_id", "region_name", "superregion_name", 
                      "decade", "population",
                      "flood_attributable_deaths", "flood_attributable_death_rate",
                      "flood_attributable_deaths_change", 
                      "flood_attributable_death_rate_change")) %>%
      group_by(superregion_name, decade) %>%
      summarize(
        region_deaths = mean(flood_attributable_deaths, na.rm = TRUE),
        region_death_rate = mean(flood_attributable_death_rate, na.rm = TRUE) * 1E5,
        region_death_change = mean(flood_attributable_deaths_change, na.rm = TRUE),
        region_death_rate_change = mean(flood_attributable_death_rate_change, 
                                        na.rm = TRUE)) %>%
      filter(!is.na(superregion_name))
    
    fig <- ggplot(data = df) +
      geom_line(aes(x = decade, y = get(variable), color = superregion_name), size = 1) +
      geom_point(aes(x = decade, y = get(variable), color = superregion_name), size = 3) +
      theme_minimal() +
      scale_color_brewer(
        name = "Region",
        palette = "Set1",
        na.value = "white") +
      ylim(0, 25) +
      labs(title = paste("Death rates for different regions, for", scenario),
           x = NULL,
           y = "Flood attributable death rate per 100K") +
      theme(legend.position = "bottom",
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 18),
            title = element_text(size = 20))
    
    return(fig)
  }
  
  ssp_126_fig <- plot_superregion_changes(ssp_126, "SSP 126")
  ssp_245_fig <- plot_superregion_changes(ssp_245, "SSP 245")
  ssp_585_fig <- plot_superregion_changes(ssp_585, "SSP 585")
  
  plot_list <- list(ssp_126_fig, ssp_245_fig, ssp_585_fig)
  
  # combine the plots into one figure
  combined_plot <- plot_grid(plotlist = plot_list, nrow = 1, ncol = 3)
  return(combined_plot)
}

# make lineplots of trends (with uncertainty ribbons) in metric over time for three climate scenarios (three lineplots)
plot_decade_trends_for_regions_with_uncty <- function(ssp_126, ssp_245, ssp_585, variable, location_type) {
  
  plot_superregion_changes <- function(df, scenario) {
    
    df$scenario <- scenario
    
    df <- df %>%
      group_by(.data[[location_type]], decade) %>%
      summarize(
        region_deaths = mean(flood_attributable_deaths, na.rm = TRUE),
        region_deaths_lower = mean(flood_attributable_deaths_lower, na.rm = TRUE),
        region_deaths_upper = mean(flood_attributable_deaths_upper, na.rm = TRUE),
        region_death_rate = mean(flood_attributable_death_rate, na.rm = TRUE) * 1E5,
        region_death_rate_lower = mean(flood_attributable_death_rate_lower, na.rm = TRUE) * 1E5,
        region_death_rate_upper = mean(flood_attributable_death_rate_upper, na.rm = TRUE) * 1E5,
        #region_death_change = mean(flood_attributable_deaths_change, na.rm = TRUE),
        #region_death_lower_change = mean(flood_attributable_deaths_lower_change, na.rm = TRUE),
        #region_death_upper_change = mean(flood_attributable_deaths_upper_change, na.rm = TRUE),
        #region_death_rate_change = mean(flood_attributable_death_rate_change, na.rm = TRUE),
        #region_death_rate_lower_change = mean(flood_attributable_death_rate_lower_change, na.rm = TRUE),
        #region_death_rate_upper_change = mean(flood_attributable_death_rate_upper_change, na.rm = TRUE)
        )
    
    df_filtered <- df %>% filter(!is.na(.data[[location_type]]))
    
    num_regions <- length(unique(df_filtered[[location_type]]))
    
    if(location_type == "region_name") {
      combined_palette <- c(
        brewer.pal(9, "Set1"),  # 9 colors
        brewer.pal(8, "Set2"),  # 8 colors
        brewer.pal(12, "Set3")[1:(num_regions - 17)]) 
    } else {
      combined_palette <- brewer.pal(9, "Set1")
    }
    
    fig <- ggplot(data = df_filtered) +
      geom_line(aes(x = decade, y = .data[[variable]], color = .data[[location_type]]), size = 1) +
      geom_point(aes(x = decade, y = .data[[variable]], color = .data[[location_type]]), size = 3) +
      geom_ribbon(aes(x = decade,
                      ymin = get(paste0(variable, "_lower")),
                      ymax = get(paste0(variable, "_upper")),
                      fill = get(location_type)),
                  alpha = 0.2) +
      theme_minimal() +
      scale_color_manual(values = combined_palette) +
      scale_fill_manual(values = combined_palette) +
      ylim(0, 150) +
      labs(title = paste("Death rates for different regions, for", scenario),
           x = NULL,
           y = "Flood attributable death rate per 100K") +
      theme(legend.position = "bottom",
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 18),
            title = element_text(size = 20))
    
    return(fig)
  }
  
  ssp_126_fig <- plot_superregion_changes(ssp_126, "SSP 126")
  ssp_245_fig <- plot_superregion_changes(ssp_245, "SSP 245")
  ssp_585_fig <- plot_superregion_changes(ssp_585, "SSP 585")
  
  plot_list <- list(ssp_126_fig, ssp_245_fig, ssp_585_fig)
  
  # combine the plots into one figure
  combined_plot <- plot_grid(plotlist = plot_list, nrow = 1, ncol = 3)
  return(combined_plot)
}



# Miscellaneous --------------------------------------------------------------

# read in results file after specifying scenario of interest
get_file <- function(flood_size, climate_future, result_type) {
  
  file <- paste0(dir, "results/", flood_size, climate_future,
                 "/country/country_", result_type, "_findings_",
                 flood_size, climate_future, ".csv")
  return(file)
}

# note: not used as of June 21
# function to perform division while handling NA and zero values
divide_safely <- function(deaths, pop) {
  result <- deaths / pop
  result[is.infinite(result) | is.nan(result)] <- NA
  return(result)
}

# function to sample a portion of raster cells
sample_raster <- function(raster, sample_fraction = 0.01) {
  # Get the number of cells
  n_cells <- ncell(raster)
  # Calculate the number of cells to sample
  n_sample <- ceiling(n_cells * sample_fraction)
  # Sample cell indices
  sample_indices <- sample(1:n_cells, n_sample)
  # Extract sampled cells
  sampled_values <- values(raster)[sample_indices]
  # Create a new raster with the sampled values
  sampled_raster <- raster
  values(sampled_raster) <- NA
  values(sampled_raster)[sample_indices] <- sampled_values
  return(sampled_raster)
}



