# This script produces statistical summaries and graphical results
# from the output of ba_rand_test.R in the different regions and subregions

library(tidyverse)

cell_km2 <- 0.215 # MODIS raster cell area in km2

# Codes (used in file names) and names (to display in plots) for regions and subregions
reg_codes <- c("northam", "eurasia", "westna", "eastna", 
               "scand", "eurus", "wsib", "esib")
reg_names <- c("North America", "Eurasia", "West North Am.", "East North Am.",
               "Scandinavia", "Eur. Russia", "West Siberia", "East Siberia")


# Calculate summary statistics from randomization test --------------------

calc_stats <- function(reg_code) {
    rand_out <- readRDS(paste0("res_rand_multi_", reg_code, ".rds"))
    obs_out <- readRDS(paste0("res_obs_", reg_code, ".rds"))    
    
    # Combine the distribution of # of years burned by cell from the
    # 200 randomizations (identified by "sim" ID column) into one data frame
    # and add 0s (when no cell with that burn count) with complete function
    burn_counts <- map_dfr(rand_out, ~ as.data.frame(.$tab), .id = "sim") %>%
        mutate(burn_counts = as.integer(as.character(burn_counts))) %>%
        complete(sim, burn_counts, fill = list(Freq = 0))
    
    # Calculate the mean and 95% interval of cell frequencies 
    #  for each value of burn_counts (# of years with fire) across simulations
    burn_stats <- group_by(burn_counts, burn_counts) %>%
        summarize(mean = mean(Freq), lo = quantile(Freq, 0.025),
                  hi = quantile(Freq, 0.975))
    
    # Combine with distribution from original data
    burn_obs <- as.data.frame(obs_out$tab) %>%
        mutate(burn_counts = as.integer(as.character(burn_counts))) %>%
        rename(obs = Freq)
    burn_stats <- full_join(burn_stats, burn_obs)
    # Replace NAs with 0s (when a burn_counts value is absent from simulated or observed data)
    burn_stats <- replace_na(burn_stats, list(mean = 0, lo = 0, hi = 0, obs = 0))
    
    # Calculate the distribution of time between fires for each randomization output
    # and combine in one data frame
    ret_counts <- map_dfr(rand_out, ~ as.data.frame(table(.$ret$dt)), .id = "sim") %>%
        mutate(Var1 = as.integer(as.character(Var1))) %>%
        complete(sim, Var1, fill = list(Freq = 0))
    # Similar to above, get the mean and 95% interval for the cell frequencies for
    # each value of dt (years between fires) across simulations, then
    # combine with observed values in original data
    ret_stats <- rename(ret_counts, dt = Var1) %>%
        group_by(dt) %>%
        summarize(mean = mean(Freq), lo = quantile(Freq, 0.025),
                  hi = quantile(Freq, 0.975)) 
    ret_obs <- as.data.frame(table(obs_out$ret$dt)) %>%
        mutate(Var1 = as.integer(as.character(Var1))) %>%
        rename(dt = Var1, obs = Freq)
    ret_stats <- full_join(ret_stats, ret_obs)
    ret_stats <- replace_na(ret_stats, list(mean = 0, lo = 0, hi = 0, obs = 0))
    
    # Return output as a list and save to disk
    stats_out <- lst(burn_stats, ret_stats)
    saveRDS(stats_out, paste0("res_stats_", reg_code, ".rds"))
    stats_out
}

# Apply function above to all subregions and combine results into list

#res <- map(reg_codes, calc_stats) %>%
#    setNames(reg_names)

# or get from disk
res <- map(paste0("res_stats_", reg_codes, ".rds"), readRDS) %>%
    setNames(reg_names)


# Years with fire ---------------------------------------------------------

# The number of cells with 0 fires in burn_stats table includes cells that are
# as located in water or outside boreal biomes, need to count them from mask
# and subtract numbers from table

mask_dir <- "data/cell_masks"

count_na <- function(reg_code) {
    rast <- raster(file.path(mask_dir, paste0("cells_na025_", reg_code, ".tif")))
    cellStats(rast == 0, sum)
}

reg_na_counts <- map_dbl(reg_codes, count_na)

for (i in seq_along(res)) {
    for (j in c("mean", "lo", "hi", "obs")) {
        res[[i]]$burn_stats[1, j] <- res[[i]]$burn_stats[1, j] - reg_na_counts[i]
    }
}

# Combine all the burn_stats tables from all regions into one data frame
burn_stats <- map_df(res, "burn_stats", .id = "region")

# Sum counts for cells with 4+ fires into same category
burn_stats2 <- burn_stats %>%
    mutate(burn_counts = ifelse(burn_counts > 4, 4, burn_counts)) %>%
    group_by(region, burn_counts) %>%
    summarize_all(.funs = "sum")
# Manually add row of 0s for 4+ fires category in Scandinavia
burn_stats2 <- rbind(burn_stats2, data.frame(region = "Scandinavia", burn_counts = 4,
                                             mean = 0, lo = 0, hi = 0, obs = 0))
# Convert numbers of cells to area in km2
burn_stats2 <- mutate(burn_stats2, mean = mean * cell_km2,
                      lo = lo * cell_km2, hi = hi * cell_km2, obs = obs * cell_km2)

# Pivot data in burn_stats2 to put observed and simulated (lo/mean/hi) statistics
# side to side (for results presentation, redundant values need to be removed in some columns)
burn_tab <- pivot_longer(burn_stats2, cols = c("hi", "mean", "lo"),
                         names_to = "stat", values_to = "area") %>%
    ungroup() %>%
    nest_by(region, burn_counts, obs) %>%
    pivot_wider(names_from = "burn_counts", values_from = c("obs", "data")) %>%
    unnest()

# Bias in total fire area (relative difference between mean of simulations and observation)
# The is negative and due to fires being "pushed out" of study area by random translation
group_by(burn_stats, region) %>%
    summarize(bias = sum(burn_counts * mean) / sum(burn_counts * obs) - 1)


# Return times ------------------------------------------------------------

# Combine all the ret_stats tables from all regions into one data frame,
# convert cell counts to areas in km2
ret_stats <- map_df(res, "ret_stats", .id = "region")
ret_stats2 <- mutate(ret_stats, mean = mean * cell_km2,
                     lo = lo * cell_km2, hi = hi * cell_km2, obs = obs * cell_km2)
ret_stats2$region <- factor(ret_stats2$region, levels = reg_names)

# Produce graph of simulated vs. observed distribution of years between fires by region
ggplot(ret_stats2, aes(x = dt, y = mean)) +
    geom_pointrange(aes(ymin = lo, ymax = hi), fatten = 2) +
    geom_point(aes(y = obs), color = "red") +
    geom_line(aes(y = obs), color = "red") +
    labs(x = "Time between fires", y = "Area (sq. km)") +
    facet_wrap(~ region, ncol = 2, scale = "free_y") +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))


# Map of study area -------------------------------------------------------

library(raster)
library(stars)
library(sf)
library(spData)
data(world)

na_mask <- read_stars(file.path(mask_dir, "cells_na025_northam.tif"), proxy = TRUE)
eu_mask <- read_stars(file.path(mask_dir, "cells_na025_eurasia.tif"), proxy = TRUE)
bbox = st_bbox(na_mask)

ggplot(world) +
    labs(x = "", y = "") +
    geom_stars(data = na_mask, downsample = 10) +
    geom_stars(data = eu_mask, downsample = 10) +
    geom_sf(fill = NA) +
    coord_sf(crs = st_crs(na_mask), ylim = c(5000000, 8000000)) +
    scale_fill_gradient(low = "white", high = "darkgreen") +
    theme_minimal() +
    theme(legend.position = "none")
