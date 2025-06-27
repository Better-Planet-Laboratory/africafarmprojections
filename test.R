library(tidyverse)
library(sf)
library(ggplot2)

# =============================================================================
# AFRICAN FARM SIZE DISTRIBUTION ANALYSIS - FIXED VERSION
# =============================================================================

# Configuration ----------------------------------------------------------------
MAX_SAMPLE_SIZE <- 5000  # Maximum samples per farm size class (computational efficiency)
YEARS_TARGET <- c(2000, 2010, 2020, 2030, 2040, 2050, 2060)

# All African ISO3 codes (54 sovereign states and territories)
AFRICAN_ISO3_CODES <- c(
  "DZA", "AGO", "BEN", "BWA", "BFA", "BDI", "CPV", "CMR", "CAF",  
  "TCD", "COM", "COD", "COG", "CIV", "DJI", "EGY", "GNQ", "ERI",  
  "SWZ", "ETH", "GAB", "GMB", "GHA", "GIN", "GNB", "KEN", "LSO",  
  "LBR", "LBY", "MDG", "MWI", "MLI", "MRT", "MUS", "MYT", "MAR",  
  "MOZ", "NAM", "NER", "NGA", "REU", "RWA", "SHN", "STP", "SEN",  
  "SYC", "SLE", "SOM", "ZAF", "SSD", "SDN", "TZA", "TGO", "TUN",  
  "UGA", "ESH", "ZMB", "ZWE"
)

# Farm size class column names
FARM_SIZE_COLS <- c("N0_1", "N1_2", "N2_5", "N5_10", "N10_20", "N20_")

# Helper Functions =============================================================

#' Fill missing admin unit values with country averages
fill_missing_admin_units <- function(data) {
  countries_with_na <- data %>%
    st_drop_geometry() %>%
    filter(if_any(all_of(FARM_SIZE_COLS), is.na)) %>%
    pull(GID_0) %>%
    unique()
  
  if(length(countries_with_na) > 0) {
    cat("Countries with NA values in admin units:", 
        paste(countries_with_na, collapse = ", "), "\n")
  }
  
  data %>%
    group_by(GID_0) %>%
    mutate(across(
      all_of(FARM_SIZE_COLS),
      ~ifelse(is.na(.), mean(., na.rm = TRUE), .)
    )) %>%
    ungroup()
}

#' Calculate proportions of farms by admin unit within countries
calculate_proportions <- function(data) {
  data %>%
    st_drop_geometry() %>%
    mutate(across(all_of(FARM_SIZE_COLS), ~replace_na(., 0))) %>%
    mutate(
      total_farms_admin = rowSums(select(., all_of(FARM_SIZE_COLS)))
    ) %>%
    group_by(GID_0) %>%
    mutate(
      total_farms_country = sum(total_farms_admin, na.rm = TRUE),
      prop_admin_of_country = ifelse(
        total_farms_country > 0, 
        total_farms_admin / total_farms_country, 
        0
      ),
      across(
        all_of(FARM_SIZE_COLS),
        list(prop = ~ifelse(total_farms_admin > 0, . / total_farms_admin, 0)),
        .names = "{.col}_prop"
      )
    ) %>%
    ungroup()
}

#' Sample farm sizes within a size class with weighted sampling
sample_farm_sizes <- function(n_farms, min_ha, max_ha, 
                              is_large_class = FALSE, 
                              max_sample = MAX_SAMPLE_SIZE) {
  
  if (is.na(n_farms) || !is.finite(n_farms) || n_farms <= 0) {
    return(numeric(0))
  }
  
  # Determine number of samples and their weight
  if (n_farms < 1) {
    num_samples <- 1
    sample_weight <- n_farms 
  } else {
    num_samples <- min(round(n_farms), max_sample)
    if (num_samples == 0) num_samples <- 1
    sample_weight <- n_farms / num_samples
  }
  
  # Generate samples
  if (is_large_class) {
    # Log-normal for large farms (20+ ha)
    sizes <- rlnorm(num_samples, meanlog = log(50), sdlog = 0.8)
    sizes <- pmax(min_ha, pmin(max_ha, sizes))
  } else {
    # Uniform for other classes
    sizes <- runif(num_samples, min = min_ha, max = max_ha)
  }
  
  attr(sizes, "weight") <- sample_weight
  sizes
}

#' Apply mean shift to match target agricultural area
mean_shift_distribution <- function(sizes, target_mean, current_mean) {
  if (length(sizes) == 0 || current_mean == 0 || 
      !is.finite(current_mean) || !is.finite(target_mean)) {
    return(sizes)
  }
  
  scaling_factor <- target_mean / current_mean
  shifted_sizes <- sizes * scaling_factor
  
  attr(shifted_sizes, "weight") <- attr(sizes, "weight")
  shifted_sizes
}

#' Rebin farm sizes into standard classes
rebin_farms <- function(sizes, weights = NULL) {
  if (length(sizes) == 0) {
    return(c(N0_1 = 0, N1_2 = 0, N2_5 = 0, N5_10 = 0, N10_20 = 0, N20_ = 0))
  }
  
  if (is.null(weights)) {
    weights <- attr(sizes, "weight")
    if (is.null(weights)) weights <- rep(1, length(sizes))
  }
  
  if (length(weights) == 1 && length(sizes) > 1) {
    weights <- rep(weights, length(sizes))
  }
  
  c(
    N0_1 = sum(weights[sizes > 0 & sizes <= 1], na.rm = TRUE),
    N1_2 = sum(weights[sizes > 1 & sizes <= 2], na.rm = TRUE),
    N2_5 = sum(weights[sizes > 2 & sizes <= 5], na.rm = TRUE),
    N5_10 = sum(weights[sizes > 5 & sizes <= 10], na.rm = TRUE),
    N10_20 = sum(weights[sizes > 10 & sizes <= 20], na.rm = TRUE),
    N20_ = sum(weights[sizes > 20], na.rm = TRUE)
  )
}

#' Preprocess fsfix data to handle NA NAME_1 values
preprocess_fsfix <- function(fsfix_data) {
  fsfix_name1_present <- fsfix_data %>% filter(!is.na(NAME_1))
  fsfix_name1_missing <- fsfix_data %>% filter(is.na(NAME_1))
  
  if (nrow(fsfix_name1_missing) > 0) {
    cat("Found", nrow(fsfix_name1_missing), 
        "rows with NA NAME_1. Aggregating by country...\n")
    
    fsfix_name1_aggregated <- fsfix_name1_missing %>%
      group_by(GID_0, NAME_0) %>%
      summarise(
        across(all_of(FARM_SIZE_COLS), sum, na.rm = TRUE),
        geometry = st_union(geometry),
        .groups = "drop"
      ) %>%
      mutate(NAME_1 = NAME_0)
    
    return(bind_rows(fsfix_name1_present, fsfix_name1_aggregated))
  }
  
  return(fsfix_name1_present)
}

#' Calculate continental average distribution
get_continental_distribution <- function(fsfix_props) {
  fsfix_props %>%
    filter(total_farms_country > 0) %>%
    select(GID_0, all_of(FARM_SIZE_COLS), total_farms_country) %>%
    group_by(GID_0) %>%
    summarise(
      across(all_of(FARM_SIZE_COLS), \(x) sum(x, na.rm = TRUE)),
      country_total = first(total_farms_country),
      .groups = "drop"
    ) %>%
    summarise(
      N0_1_prop = mean(N0_1 / country_total, na.rm = TRUE),
      N1_2_prop = mean(N1_2 / country_total, na.rm = TRUE),
      N2_5_prop = mean(N2_5 / country_total, na.rm = TRUE),
      N5_10_prop = mean(N5_10 / country_total, na.rm = TRUE),
      N10_20_prop = mean(N10_20 / country_total, na.rm = TRUE),
      N20__prop = mean(N20_ / country_total, na.rm = TRUE)
    )
}

#' Process countries with both fsfix and SSP2T data (Case A)
process_case_a <- function(fsfix_props, year_proj, countries_both, 
                           ethiopia_growth, year) {
  fsfix_props %>%
    filter(GID_0 %in% countries_both) %>%
    left_join(year_proj, by = c("GID_0" = "ISO3_CODE")) %>%
    mutate(
      Year = year,
      growth_if_missing = approx(
        ethiopia_growth$Year, 
        ethiopia_growth$growth_factor, 
        year, 
        rule = 2
      )$y,
      final_total_farms = ifelse(
        is.na(Median.smooth), 
        total_farms_country * growth_if_missing, 
        Median.smooth
      ),
      treatment_flag = ifelse(
        is.na(Median.smooth),
        "Case A/C Hybrid: fsfix data exists, but SSP2T missing for this year. Used ETH growth.",
        "Case A: Data from both fsfix and SSP2T"
      ),
      across(
        all_of(FARM_SIZE_COLS),
        ~ ifelse(total_farms_country > 0, .x / total_farms_country * final_total_farms, 0),
        .names = "{.col}_corrected"
      )
    )
}

#' Process countries only in SSP2T but with fsfix geometry (Case B with geometry)
process_case_b_with_geometry <- function(year_proj, countries_fsfix_no_data_ssp2t, 
                                         continental_dist, year, fsfix_props) {
  
  # Get the fsfix structure for these countries (NAME_0, NAME_1, geometry info)
  fsfix_structure <- fsfix_props %>%
    filter(GID_0 %in% countries_fsfix_no_data_ssp2t) %>%
    select(GID_0, NAME_0, NAME_1)
  
  year_proj %>%
    filter(ISO3_CODE %in% countries_fsfix_no_data_ssp2t) %>%
    cross_join(continental_dist) %>%
    mutate(
      Year = year,
      GID_0 = ISO3_CODE,
      N0_1_corrected = N0_1_prop * Median.smooth,
      N1_2_corrected = N1_2_prop * Median.smooth,
      N2_5_corrected = N2_5_prop * Median.smooth,
      N5_10_corrected = N5_10_prop * Median.smooth,
      N10_20_corrected = N10_20_prop * Median.smooth,
      N20__corrected = N20__prop * Median.smooth,
      treatment_flag = "Case B: SSP2T only - continental average disaggregation"
    ) %>%
    # Join with fsfix structure to get the correct NAME_0, NAME_1 combinations
    left_join(fsfix_structure, by = "GID_0")
}

#' Process countries only in fsfix (Case C)
process_case_c <- function(fsfix_props, countries_only_fsfix, 
                           ethiopia_growth, year) {
  # Get growth factor for the year
  growth <- 1
  if (year %in% ethiopia_growth$Year) {
    growth <- ethiopia_growth$growth_factor[ethiopia_growth$Year == year]
  } else if (year > max(ethiopia_growth$Year)) {
    growth <- ethiopia_growth$growth_factor[ethiopia_growth$Year == max(ethiopia_growth$Year)]
  } else if (year < min(ethiopia_growth$Year)) {
    growth <- ethiopia_growth$growth_factor[ethiopia_growth$Year == min(ethiopia_growth$Year)]
  } else {
    growth <- approx(ethiopia_growth$Year, ethiopia_growth$growth_factor, year)$y
  }
  
  fsfix_props %>%
    filter(GID_0 %in% countries_only_fsfix) %>%
    mutate(
      Year = year,
      across(
        all_of(FARM_SIZE_COLS),
        ~ .x * growth,
        .names = "{.col}_corrected"
      ),
      Median.smooth = total_farms_country * growth,
      Agarea.kha = NA,
      treatment_flag = "Case C: fsfix only - Ethiopia growth trend applied"
    )
}

#' Apply mean shift and rebin for a country
process_country_mean_shift <- function(country_data, country_name, year) {
  # Skip if no area data available
  if (all(is.na(country_data$Agarea.kha)) || 
      all(is.na(country_data$Median.smooth)) || 
      all(country_data$Median.smooth == 0)) {
    return(country_data %>%
             mutate(across(ends_with("_corrected"), ~., 
                           .names = "{str_remove(.col, '_corrected')}")))
  }
  
  # Calculate target mean farm size
  A_m <- unique(country_data$Agarea.kha) * 1000 / unique(country_data$Median.smooth)
  
  if (!is.finite(A_m)) {
    warning(paste("Invalid target mean farm size for", country_name, 
                  "in year", year, ". Skipping mean shift."))
    return(country_data %>%
             mutate(across(ends_with("_corrected"), ~., 
                           .names = "{str_remove(.col, '_corrected')}")))
  }
  
  # Generate farm size distributions
  all_sizes <- list()
  all_sizes_df <- list()
  
  for (i in 1:nrow(country_data)) {
    admin_sizes <- list(
      sample_farm_sizes(country_data$N0_1_corrected[i], 0.01, 1),
      sample_farm_sizes(country_data$N1_2_corrected[i], 1, 2),
      sample_farm_sizes(country_data$N2_5_corrected[i], 2, 5),
      sample_farm_sizes(country_data$N5_10_corrected[i], 5, 10),
      sample_farm_sizes(country_data$N10_20_corrected[i], 10, 20),
      sample_farm_sizes(country_data$N20__corrected[i], 20, 1000, is_large_class = TRUE)
    )
    all_sizes[[i]] <- admin_sizes
    
    if (length(unlist(admin_sizes)) > 0) {
      dist_df <- tibble(
        Year = year,
        GID_0 = country_name,
        NAME_0 = country_data$NAME_0[i],
        NAME_1 = country_data$NAME_1[i],
        farm_size = unlist(admin_sizes),
        weight = unlist(lapply(admin_sizes, function(s) {
          if(length(s) > 0) rep(attr(s, "weight"), length(s)) else numeric(0)
        }))
      )
      all_sizes_df[[i]] <- dist_df
    }
  }
  
  # Store distribution data
  if (length(all_sizes_df) > 0) {
    assign(paste0("dist_", country_name, "_", year), 
           bind_rows(all_sizes_df), 
           envir = .GlobalEnv)
  }
  
  # Calculate observed mean
  all_sizes_combined <- unlist(all_sizes)
  weights_combined <- unlist(lapply(all_sizes, function(x) {
    unlist(lapply(x, function(s) {
      if(length(s) > 0) rep(attr(s, "weight"), length(s)) else numeric(0)
    }))
  }))
  
  A_o <- if (length(all_sizes_combined) > 0 && sum(weights_combined, na.rm = TRUE) > 0) {
    weighted.mean(all_sizes_combined, weights_combined, na.rm = TRUE)
  } else {
    A_m
  }
  
  # Apply mean shift and rebin
  shifted_data <- country_data
  for (i in 1:nrow(country_data)) {
    if (length(all_sizes[[i]]) > 0 && A_o > 0) {
      shifted_sizes <- lapply(all_sizes[[i]], function(sizes) {
        if (length(sizes) > 0) mean_shift_distribution(sizes, A_m, A_o) else sizes
      })
      
      all_shifted_combined <- unlist(shifted_sizes)
      weights_shifted <- unlist(lapply(shifted_sizes, function(s) {
        if(length(s) > 0) rep(attr(s, "weight"), length(s)) else numeric(0)
      }))
      
      rebinned <- rebin_farms(all_shifted_combined, weights = weights_shifted)
      
      # Preserve total farm count
      original_total <- sum(country_data[i, paste0(FARM_SIZE_COLS, "_corrected")], na.rm = TRUE)
      current_total <- sum(rebinned, na.rm = TRUE)
      
      if (abs(original_total) < .Machine$double.eps^0.5) {
        rebinned[] <- 0
      } else if (abs(current_total) < .Machine$double.eps^0.5) {
        rebinned <- as.numeric(country_data[i, paste0(FARM_SIZE_COLS, "_corrected")])
        names(rebinned) <- names(rebin_farms(numeric(0)))
      } else {
        rebinned <- rebinned * (original_total / current_total)
      }
      
      shifted_data[i, FARM_SIZE_COLS] <- as.list(rebinned)
    } else {
      shifted_data[i, FARM_SIZE_COLS] <- 0
    }
  }
  
  shifted_data
}

#' Create validation plots
create_validation_plots <- function(final_output, farm_projections, 
                                    shifted_distributions) {
  # Plot 1: Total farms validation
  validation_data <- final_output %>%
    st_drop_geometry() %>%
    filter(!is.na(Median.smooth)) %>%
    group_by(Year, GID_0) %>%
    summarise(
      total_output = sum(N0_1 + N1_2 + N2_5 + N5_10 + N10_20 + N20_, na.rm = TRUE),
      expected_total = first(Median.smooth),
      .groups = "drop"
    )
  
  p1 <- ggplot(validation_data, aes(x = expected_total, y = total_output)) +
    geom_jitter(aes(color = as.factor(Year)), alpha = 0.6, 
                width = 0.1, height = 0.1) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma) +
    labs(
      title = "Validation: SSP2T Projections vs Output Totals",
      subtitle = "Points should fall on the red line",
      x = "SSP2T Expected Total (log scale)",
      y = "Output Total (log scale)",
      color = "Year"
    ) +
    theme_minimal()
  
  ggsave("figs/validation_totals_check.png", p1, width = 10, height = 8, dpi = 300)
  
  # Plot 2: Map comparison
  map_data <- final_output %>%
    filter(Year %in% c(2000, 2050)) %>%
    mutate(
      total_farms = N0_1 + N1_2 + N2_5 + N5_10 + N10_20 + N20_,
      farms_category = cut(
        total_farms,
        breaks = c(0, 1000, 5000, 10000, 50000, 100000, 500000, Inf),
        labels = c("<1k", "1-5k", "5-10k", "10-50k", "50-100k", "100-500k", ">500k"),
        include.lowest = TRUE
      )
    )
  
  p2 <- ggplot(map_data) +
    geom_sf(aes(fill = farms_category), color = "white", size = 0.1) +
    facet_wrap(~Year, ncol = 2) +
    scale_fill_viridis_d(name = "Total Farms", option = "plasma", na.value = "grey90") +
    labs(
      title = "Farm Numbers Across Africa: 2000 vs 2050",
      subtitle = "Total farms per administrative unit"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom", axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  ggsave("figs/africa_farms_2000_2050.png", p2, width = 14, height = 8, dpi = 300)
  
  # Plot 3: Farm size distribution evolution
  distribution_plot_data <- shifted_distributions %>%
    filter(Year %in% c(2000, 2020, 2050)) %>%
    mutate(Year = as.factor(Year))
  
  p3 <- ggplot(distribution_plot_data, 
               aes(x = farm_size_shifted, weight = weight, 
                   color = Year, fill = Year)) +
    geom_density(alpha = 0.3, adjust = 1.5) +
    scale_x_log10(
      breaks = c(0.1, 0.5, 1, 2, 5, 10, 20, 50, 100, 500),
      labels = c("0.1", "0.5", "1", "2", "5", "10", "20", "50", "100", "500")
    ) +
    scale_color_viridis_d(option = "plasma", end = 0.8) +
    scale_fill_viridis_d(option = "plasma", end = 0.8) +
    labs(
      title = "Evolution of Farm Size Distribution in Africa",
      subtitle = "Kernel density estimates for all countries",
      x = "Farm Size (hectares, log scale)",
      y = "Density"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave("figs/africa_farm_size_distribution_evolution.png", p3, 
         width = 10, height = 6, dpi = 300)
  
  # Plot 4: Agricultural area validation
  validation_area_data <- shifted_distributions %>%
    select(Year, GID_0, farm_size_shifted, weight) %>%
    left_join(
      farm_projections %>%
        distinct(Year, ISO3_CODE, Median.smooth, Agarea.kha), 
      by = c("Year" = "Year", "GID_0" = "ISO3_CODE")
    ) %>%
    filter(!is.na(Agarea.kha) & !is.na(Median.smooth) & Median.smooth > 0) %>%
    group_by(Year, GID_0) %>%
    summarise(
      output_mean_farm_size = weighted.mean(farm_size_shifted, weight, na.rm = TRUE),
      expected_total_farms = first(Median.smooth),
      original_ag_area_kha = first(Agarea.kha),
      .groups = "drop"
    ) %>%
    mutate(
      calculated_output_area_ha = output_mean_farm_size * expected_total_farms,
      expected_area_ha = original_ag_area_kha * 1000
    ) %>%
    filter(is.finite(calculated_output_area_ha) & 
             is.finite(expected_area_ha) & 
             expected_area_ha > 0)
  
  p4 <- ggplot(validation_area_data, 
               aes(x = expected_area_ha, y = calculated_output_area_ha)) +
    geom_jitter(aes(color = as.factor(Year)), alpha = 0.6, 
                width = 0.1, height = 0.1) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma) +
    labs(
      title = "Validation: Agricultural Area",
      subtitle = "Points should fall on the red 1:1 line",
      x = "SSP2T Projected Agricultural Area (hectares, log scale)",
      y = "Output Implied Agricultural Area (hectares, log scale)",
      color = "Year"
    ) +
    theme_minimal()
  
  ggsave("figs/validation_areas_check.png", p4, width = 10, height = 8, dpi = 300)
}
# Main Analysis ================================================================

main <- function() {
  cat("Starting African Farm Size Distribution Analysis\n")
  cat("==============================================\n\n")
  
  # Load data ------------------------------------------------------------------
  cat("Loading data...\n")
  fsfix_original <- st_read("input/fsfix.shp", quiet = TRUE)
  ssp2t <- readRDS("input/ssp2t.rds")
  farm_projections <- ssp2t[[4]] %>%
    filter(ISO3_CODE %in% AFRICAN_ISO3_CODES)
  
  # Preprocess fsfix data ------------------------------------------------------
  cat("Preprocessing fsfix data...\n")
  fsfix <- preprocess_fsfix(fsfix_original)
  
  # Process fsfix data ---------------------------------------------------------
  cat("Processing fsfix data...\n")
  fsfix_filled <- fill_missing_admin_units(fsfix)
  fsfix_props <- calculate_proportions(fsfix_filled)
  
  # FIXED: Identify country categories more carefully -------------------------
  # Countries with actual farm data (non-zero totals)
  countries_with_farm_data <- fsfix_props %>%
    filter(total_farms_country > 0) %>%
    pull(GID_0) %>%
    unique()
  
  # Countries that appear in fsfix (regardless of farm data)
  countries_in_fsfix <- unique(fsfix_props$GID_0)
  
  # Countries in SSP2T
  countries_ssp2t <- unique(farm_projections$ISO3_CODE)
  
  # Case A: Countries with actual farm data in fsfix AND in SSP2T
  countries_both <- intersect(countries_with_farm_data, countries_ssp2t)
  
  # Case B: Countries in fsfix with geometry but no farm data, also in SSP2T
  # These get SSP2T data with continental average disaggregation + fsfix geometry
  countries_fsfix_no_data_but_ssp2t <- intersect(
    setdiff(countries_in_fsfix, countries_with_farm_data), 
    countries_ssp2t
  )
  
  # Case C: Countries with farm data in fsfix but not in SSP2T
  countries_only_fsfix <- setdiff(countries_with_farm_data, countries_ssp2t)
  
  cat("\nCountry categories:\n")
  cat("- Both datasets (Case A):", length(countries_both), "countries\n")
  cat("- fsfix geometry but no farm data, has SSP2T (Case B):", length(countries_fsfix_no_data_but_ssp2t), "countries\n")
  cat("- Only fsfix (Case C):", length(countries_only_fsfix), "countries\n\n")
  
  # Get continental distribution and Ethiopia growth ---------------------------
  continental_dist <- get_continental_distribution(fsfix_props)
  
  ethiopia_growth <- farm_projections %>%
    filter(ISO3_CODE == "ETH", Year >= min(YEARS_TARGET)) %>%
    arrange(Year) %>%
    mutate(growth_factor = Median.smooth / first(Median.smooth)) %>%
    select(Year, growth_factor)
  
  # Process each year ----------------------------------------------------------
  corrected_numbers_list <- list()
  final_results_list <- list()
  full_distributions_list <- list()
  
  for (year in YEARS_TARGET) {
    cat("\nProcessing year", year, "...\n")
    
    year_proj <- farm_projections %>%
      filter(Year == year) %>%
      select(ISO3_CODE, Median.smooth, Agarea.kha)
    
    # Process different country cases
    corrected_both <- process_case_a(fsfix_props, year_proj, countries_both, 
                                     ethiopia_growth, year)
    
    # Case B: Countries with fsfix geometry but no farm data - treat as Case B with geometry
    corrected_b_fsfix_geo_ssp2t <- process_case_b_with_geometry(year_proj, countries_fsfix_no_data_but_ssp2t, 
                                                                continental_dist, year, fsfix_props)
    
    # Case C: fsfix only
    corrected_fsfix_only <- process_case_c(fsfix_props, countries_only_fsfix, 
                                           ethiopia_growth, year)
    
    # Combine results
    year_corrected <- bind_rows(
      corrected_both, 
      corrected_b_fsfix_geo_ssp2t,
      corrected_fsfix_only
    )
    corrected_numbers_list[[as.character(year)]] <- year_corrected
    
    # Apply mean shift
    cat("  Applying mean shift...\n")
    year_distributions <- list()
    
    year_results <- year_corrected %>%
      group_by(GID_0) %>%
      group_modify(function(country_data, key) {
        country_name <- key$GID_0
        result <- process_country_mean_shift(country_data, country_name, year)
        
        # Extract stored distributions if any
        dist_name <- paste0("dist_", country_name, "_", year)
        if (exists(dist_name, envir = .GlobalEnv)) {
          year_distributions[[country_name]] <<- get(dist_name, envir = .GlobalEnv)
          rm(list = dist_name, envir = .GlobalEnv)
        }
        
        result
      }) %>%
      ungroup()
    
    if (length(year_distributions) > 0) {
      full_distributions_list[[as.character(year)]] <- bind_rows(year_distributions)
    }
    
    final_results_list[[as.character(year)]] <- year_results
  }
  
  # Combine all results --------------------------------------------------------
  cat("\nCombining results...\n")
  corrected_numbers_all <- bind_rows(corrected_numbers_list)
  final_results <- bind_rows(final_results_list)
  full_distributions <- bind_rows(full_distributions_list)
  
  # Save distributions ---------------------------------------------------------
  cat("Saving distribution data...\n")
  saveRDS(full_distributions, "output/full_farm_distributions_all_countries.rds")
  
  # Create shifted distributions
  shifted_distributions_list <- list()
  
  for (year in names(full_distributions_list)) {
    year_dist <- full_distributions_list[[year]]
    year_results <- final_results_list[[year]]
    
    country_means <- year_results %>%
      filter(!is.na(Agarea.kha) & !is.na(Median.smooth) & Median.smooth > 0) %>%
      group_by(GID_0) %>%
      summarise(A_m = first(Agarea.kha) * 1000 / first(Median.smooth), 
                .groups = "drop")
    
    shifted_dist <- year_dist %>%
      left_join(country_means, by = "GID_0") %>%
      group_by(GID_0) %>%
      mutate(
        A_o = weighted.mean(farm_size, weight, na.rm = TRUE),
        scaling_factor = ifelse(!is.na(A_m) & is.finite(A_m) & A_o > 0, 
                                A_m / A_o, 1),
        farm_size_shifted = farm_size * scaling_factor
      ) %>%
      ungroup()
    
    shifted_distributions_list[[year]] <- shifted_dist
  }
  
  shifted_distributions <- bind_rows(shifted_distributions_list)
  saveRDS(shifted_distributions, "output/full_farm_distributions_shifted.rds")
  
  # FIXED: Create spatial output without duplicates ---------------------------
  cat("Creating spatial output...\n")
  
  # Start with all fsfix geometries
  base_spatial <- fsfix %>%
    select(GID_0, NAME_0, NAME_1, geometry) %>%
    crossing(Year = YEARS_TARGET)
  
  # Join with final results - this will automatically handle all cases
  # since final_results contains all processed countries with the right flags
  final_output <- base_spatial %>%
    left_join(
      final_results %>%
        select(Year, GID_0, NAME_0, NAME_1, all_of(FARM_SIZE_COLS), 
               treatment_flag, Median.smooth, Agarea.kha),
      by = c("Year", "GID_0", "NAME_0", "NAME_1")
    ) %>%
    # Handle countries that have no farm data at all
    mutate(
      across(all_of(FARM_SIZE_COLS), ~replace_na(., 0)),
      treatment_flag = case_when(
        !is.na(treatment_flag) ~ treatment_flag,
        GID_0 %in% countries_in_fsfix ~ "No farm data available",
        TRUE ~ "Not in analysis"
      )
    ) %>%
    st_as_sf()
  
  # Add countries that are in SSP2T but not in fsfix (Case B without geometry)
  # These need synthetic geometry
  countries_needing_synthetic <- countries_b_no_geometry
  
  if (length(countries_needing_synthetic) > 0) {
    synthetic_results <- final_results %>%
      filter(GID_0 %in% countries_needing_synthetic) %>%
      mutate(geometry = st_sfc(st_point(c(0, 0)), crs = st_crs(fsfix))) %>%
      st_as_sf() %>%
      select(GID_0, NAME_0, NAME_1, Year, all_of(FARM_SIZE_COLS), 
             treatment_flag, Median.smooth, Agarea.kha, geometry)
    
    # Combine with main results
    final_output <- bind_rows(final_output, synthetic_results)
  }
  # Save outputs ---------------------------------------------------------------
  cat("Saving final outputs...\n")
  cat("Saving final outputs...\n")
  
  # Only save geometries that are proper polygons
  final_output_polygons <- final_output %>% 
    filter(st_geometry_type(geometry) %in% c("POLYGON", "MULTIPOLYGON"))
  
  st_write(final_output_polygons, "output/fsfix_corrected_2000_2060.shp", 
           delete_dsn = TRUE, quiet = TRUE)
  saveRDS(final_output, "output/fsfix_corrected_2000_2060.rds")
  
  # Create summary -------------------------------------------------------------
  summary_stats <- final_output %>%
    st_drop_geometry() %>%
    group_by(Year, treatment_flag) %>%
    summarise(
      n_countries = n_distinct(GID_0),
      n_admin_units = n(),
      total_farms = sum(N0_1 + N1_2 + N2_5 + N5_10 + N10_20 + N20_, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("\nAnalysis complete! Summary by year and treatment:\n")
  print(summary_stats)
  
  # Check for duplicates
  duplicate_check <- final_output %>%
    st_drop_geometry() %>%
    group_by(Year, GID_0, NAME_1) %>%
    summarise(n_rows = n(), .groups = "drop") %>%
    filter(n_rows > 1)
  
  if (nrow(duplicate_check) > 0) {
    cat("\nWARNING: Found duplicates for these country-admin combinations:\n")
    print(duplicate_check)
  } else {
    cat("\nNo duplicates found - issue resolved!\n")
  }
  
  # Create validation plots ----------------------------------------------------
  cat("\nCreating validation plots...\n")
  create_validation_plots(final_output, farm_projections, shifted_distributions)
  
  cat("\nAll processing complete!\n")
}

# Run analysis
main()

