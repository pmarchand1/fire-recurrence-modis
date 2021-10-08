# This script performs 200 randomizations (random translations) 
# of the annual burned area maps in a given region and reports certain statistics 
# (see Output) for the randomized maps as well as the original data.
#
# Usage: This script is designed to be called from the command line 
#        (e.g. on a remote computing server) with two arguments representing
#        the region (northam and eurasia) and the subregion (could be the same as region)
# e.g. "Rscript ba_rand_test.R northam westna" for Western North America
#
# Input: The region/subregion masks (e.g. cells_na025_westna.tif) and the masked
#        annual burned area maps (e.g. northam2001_mask.tif) need to be present
#        in the same folder as the script is run from.
#
# Output: For each of the 200 randomized burned area annual maps as well as the 
#         observed (unrandomized) maps, a list with two elements:
#  "tab", a table with the number of raster cells that burned 0, 1, 2, ... times
#  "ret", a data frame of all cells that burned more than once,
#         with coordinates (x,y) of the cell, year of 2nd burn and 
#         dt, the number of years since last burn

library(dplyr)
library(tidyr)
library(stringr)
library(parallel)
library(raster)
# This is used to manually set the storage folder for temp files on the computing cluster
rasterOptions(tmpdir = "/scratch/marchanp") 

# Get command-line arguments
argv <- commandArgs(trailingOnly=TRUE)
region <- argv[1]
subregion <- argv[2]

# This function takes a raster stack containing the annual burned area maps
# and computes the "tab" and "ret" outputs described at the top of the file
calc_return_time <- function(rstack) {
    # Get the number of years with burn for each cell
    burn_counts <- getValues(calc(rstack, sum, na.rm = TRUE))
    # Extract cell indices that burned on multiple years; if none, just return 
    # the table of number of years with burn, and a NULL output for "ret"
    mult_burn <- which(burn_counts >= 2)
    if (length(mult_burn) == 0) {
        return(list(tab = table(burn_counts), ret = NULL))
    }
    # The following section extracts the (x,y) coordinates and annual burn (0/1)
    # values for each cell that burns more than once, and pivots the result 
    # into a 4-column data frame: x, y, year, burn
    cell_xy <- xyFromCell(rstack, mult_burn)
    cell_vals <- raster::extract(rstack, cell_xy)
    cell_vals <- cbind(cell_xy, cell_vals)
    cell_vals <- pivot_longer(as.data.frame(cell_vals), cols = -c(x, y),
                              names_to = "year", values_to = "burn")
    cell_vals <- mutate(cell_vals, year = as.integer(str_extract(year, "[[:digit:]]+")))
    # Retain the rows with burn presence for each cell and calculate
    # the intervals dt between consecutive burns in each cell
    cell_ret <- filter(cell_vals, burn == 1) %>%
        dplyr::select(-burn) %>%
        group_by(x, y) %>%
        mutate(dt = year - lag(year)) %>%
        filter(!is.na(dt))
    
    list(tab = table(burn_counts), ret = cell_ret)
}

# This function performs a random translation of a raster layer rast, 
# i.e. translating by a random number of cells between nmin and nmax 
#      with a random sign (+ or -) in both x and y
# Note: Because of the shift some cells on one edge will "fall off" 
#       and there will be NA cells on the opposite edge.
rand_shift <- function(rast, nmin, nmax) {
    ext <- extent(rast)
    reso <- res(rast)[1]
    dcells <- round(runif(2, nmin, nmax)) * sign(runif(2, -1, 1))
    extent(rast) <- ext + reso * rep(dcells, each = 2)
    extend(crop(rast, ext), ext)
}

# This function combines the two previous functions to perform a random translation
# of the annual burned area maps in filelist (independent translation for each year)
# and apply the calc_return_time function to the randomized raster stack
rand_return_time <- function(filelist, nmin, nmax, rmask) {
    lapply(filelist, raster) %>%
        lapply(crop, rmask) %>%
        lapply(rand_shift, nmin, nmax) %>%
        stack() %>%
        mask(mask = rmask, maskvalue = 0) %>%
        calc_return_time()
}

# Get the annual burned area maps for the region and the subregion mask
region_files <- dir(pattern = region)
mask025 <- raster(paste0("cells_na025_", subregion, ".tif"))

# Apply the rand_return_time function 200 times with a translation range of
#  20 to 100 cells (~ 10 to 50 km given the MODIS resolution).
# mc.cores can be changed depending of the number of parallel processors to use
res_rand <- mclapply(1:200, 
    function(i) rand_return_time(region_files, 20, 100, mask025), 
    mc.cores = 8)

# Apply calc_return_time for the observed (original, not randomized) burned area maps
res_obs <- lapply(region_files, raster) %>%
    lapply(crop, mask025) %>%
    stack() %>%
    mask(mask = mask025, maskvalue = 0) %>%
    calc_return_time()

# Save separately the list of 200 randomized outputs and the original data output
saveRDS(res_rand, paste0("res_rand_multi_", subregion, ".rds"))
saveRDS(res_obs, paste0("res_obs_", subregion, ".rds"))
