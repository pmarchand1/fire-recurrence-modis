# This script prepares the data files needed for the burned area
# randomization test in ba_rand_test.R, from the original burned area
# and land use data products available from MODIS.

library(sf)
library(MODIS)
library(raster)
library(rgdal)
library(tidyverse)

# Aggregate in .tif file by tile and year ---------------------------------

# The input data (in folder given as ba_orig_dir) is the MODIS monthly 
# burned area product (MCD64A1 v006). The folder should contain one subfolder
# by tile (ex.: h06v03) that contains all the monthly .hdf files associated
# with that tile.

ba_orig_dir <- "data/ba_original"
ba_annual_dir <- "data/ba_annual" # Directory to save annual aggregates
tiles <- dir(ba_orig_dir) # Get the list of tiles present in data folder
years <- 2001:2019 # Years included in analysis

# This function first gets the 12 filenames for the monthly burned area data
# for a given tile and year (by looking for the AYYYY pattern in the filenames,
# where YYYY is the year). It then loads the 12 burn date rasters and puts 
# them into a stack. Each raster contains a value between 1 and 366 indicating the
# day where there was a fire in that cell, or a value <= 0 if no fire.
# Therefore, if the max value across the 12 layers is >0, there was a fire
# in that cell for that year. The result data_sum (a layer of 0/1 for 
# absence/presence of fire) is saved as a GeoTIFF.
agg_tile_year <- function(tile, year) {
    filenames <- dir(file.path(ba_dir, tile), paste0("A", year), 
                     full.names = TRUE)
    # This gets the list of scientific datsets (with getSds) then picks the path
    # for the first dataset of the "SDS4gdal" element, corresponding to burn date.
    layernames <- map(filenames, getSds) %>%
        map(list("SDS4gdal", 1))
    data_stack <- map(layernames, readGDAL, as.is = TRUE) %>%
        map(raster) %>%
        stack()
    data_sum <- max(data_stack, na.rm = TRUE) > 0
    writeRaster(data_sum, file.path(ba_annual_dir, paste0(tile, "_", year, ".tif")))
}

# This part applies the function above to all tiles and years
df <- expand.grid(tile = tiles, year = years)
pwalk(df, agg_tile_year)


# Find boreal subset of tiles ------------------------------------------------

# This function counts the total number of cells with fire
# across multiple files given in filenames.
sum_ba <- function(filenames) {
    sum(map_dbl(filenames, ~ cellStats(raster(.), sum)))
}    

# This part takes all annual burned area files and 
# extracts the tile and year from the filename as part of a dataframe
all_files <- dir(ba_annual_dir, pattern = "tif")
file_df <- data.frame(filename = all_files) %>%
    separate(filename, c("tile", "year"), sep = "_", remove = FALSE) %>%
    mutate(filename = file.path(ba_annual_dir, filename))
# Nest the data frame by tile and apply the 
# sum_ba function above to the filenames for each tile.
file_df <- nest_by(file_df, tile)
file_df$sum_ba <- map_dbl(file_df$data, ~ sum_ba(.$filename))

# Load the MODIS sinusoidal grid and take a large section including the boreal region
modis <- st_read("data/modis_grid/modis_sinusoidal_grid_world.shp") %>%
    filter(h >= 8, h <= 28, v <= 5)
# Extract the taiga (6) and tundra (11) biomes from the WWF map
# "Terrestrial Ecoregions of the World".
biomes <- st_read("data/wwf_terr_ecos_oRn")
biomes <- filter(biomes, BIOME %in% c(6, 11)) %>%
    st_transform(st_crs(modis)) %>% # transform to MODIS coordinate system
    # st_crop removes the area outsdie the MODIS boreal section
    st_crop(modis) 
# Plot boreal biomes over MODIS grid
ggplot() + 
    geom_sf(data = modis) + 
    geom_sf(data = biomes, aes(fill = BIOME), color = NA) + 
    geom_sf_text(data = modis, aes(label = paste0("h", h, "v", v)))

# Based on sum of cells with fire in each tile and the boreal map above,
# here are all MODIS tiles with fire located in the boreal region
# North America: v02, h10-13; v03, h11-14; v04, h12-13 
# Eurasia: v2, h18-25; v03, h19-26

# If we take out tiles where only a small part of the tile is in 
# the boreal region, we get the following "core" tiles:
# North America: v02, h11-13; v03, h11-14
# Eurasia: v02, h18-24; v03, h19-24

northam_tiles <- c("h11v02", "h12v02", "h13v02", 
                   "h11v03", "h12v03", "h13v03", "h14v03")
eurasia_tiles <- c("h18v02", "h19v02", "h20v02", "h21v02", "h22v02", "h23v02", "h24v02", 
                   "h19v03", "h20v03", "h21v03", "h22v03", "h23v03", "h24v03")


# Create masks to exclude water and non-boreal biomes --------------------------

landuse_dir <- "data/landuse" # Directory for MODIS landuse data
mask_dir <- "data/cell_masks" # Directory to save resulting masks

# This function creates a the land and biome mask for region "name" comprised of
# list of tiles in "tiles", based on MODIS land cover maps (MCD12Q1) and WWF biomes map
# The result is a layer where cells in water or outside target biomes are NA
create_landbiome_mask <- function(name, tiles) {
    lu_files <- dir(landuse_dir, pattern = paste(tiles, collapse = "|"),
                    full.names = TRUE)
    # Get layer 13 (land/water map, where water = 1, land = 2) for each tile
    land_water <- map(lu_files, getSds) %>%
        map(~ readGDAL(.$SDS4gdal[13], as.is = TRUE)) %>%
        map(raster)
    
    # Merge all tiles in one raster with the given filename
    land_water$filename <- file.path(landuse_dir, paste0(name, "_landwater.tif"))
    do.call(merge, land_water)
    land_water <- raster(land_water$filename)
    
    # Transform biomes to same coordinate system as land/water mask, 
    # then crop to region covered by mask
    biomes_sub <- st_transform(biomes, st_crs(land_water)) %>%
        st_crop(st_bbox(land_water))
    
    land_biome <- mask(land_water, biomes_sub) # set all cells not in biomes to NA
    land_biome[land_biome == 1] <- NA # set all water cells to NA
    
    writeRaster(land_biome, file.path(mask_dir, paste0(name, "_landbiome_mask.tif")))
}

create_landbiome_mask("northam", northam_tiles)
create_landbiome_mask("eurasia", eurasia_tiles)


# Combine core maps -------------------------------------------------------

# Directory to save merged fire maps
ba_merged_dir <- "data/ba_merged"

merge_by_year <- function(year, tiles, name) {
    files <- file.path(ba_annual_dir, paste0(tiles, "_", year, ".tif"))
    rast_list <- map(files, raster)
    rast_list$filename <- file.path(ba_merged_dir, paste0(name, year, ".tif"))
    do.call(merge, rast_list)
}

walk(2001:2019, merge_by_year, tiles = northam_tiles, name = "northam")
walk(2001:2019, merge_by_year, tiles = eurasia_tiles, name = "eurasia")


# Mask core maps ----------------------------------------------------------

# Function to apply land mask to burn area raster for a given region (name) and year
# setting raster cells in water or outside boreal biomes to the value NA
mask_by_year <- function(year, name, mask) {
    rast_in <- raster(file.path(ba_merged_dir, paste0(name, year, ".tif")))
    mask(rast_in, mask, maskvalue = NA, 
         filename = file.path(ba_merged_dir, paste0(name, year, "_mask.tif")))
}

# Apply function above to each year for both North America and Eurasia
land_biome_na <- raster(file.path(mask_dir, "northam_landbiome_mask.tif"))
walk(2001:2019, mask_by_year, "northam", land_biome_na)

land_biome_eu <- raster(file.path(mask_dir, "eurasia_landbiome_mask.tif"))
walk(2001:2019, mask_by_year, "eurasia", land_biome_eu)


# Masks based on fraction of neighbours in water or non-boreal biome ---------

calc_frac_cell_masks <- function(region, landbiome_mask) { # region: "northam" or "eurasia"
    # Calculates fraction of cells in water or non-boreal (NA values in 
    # landbiome_mask) for a 201x201 square centered on each cell.
    na100 <- focal(is.na(landbiome_mask), w = matrix(rep(1/40401, 40401), nrow = 201),
                   filename = file.path(mask_dir, paste0("mean_na_100cells_", region, ".tif")))
    
    # mask025 is a new mask = 0 if the cell is in water, in a non-boreal biome or
    # has >25% of neighbours in water or non-boreal biomes, = 1 otherwise
    # This will be use to later exclude cells from the analysis if they are
    # near many cells outside of region of interest
    mask025 <- na100 < 0.25 & !is.na(landbiome_mask)
    mask025[is.na(mask025)] <- 0
    writeRaster(mask025, file.path(mask_dir, paste0("cells_na025_", region, ".tif")))
}

calc_frac_cell_masks("northam", land_biome_na)
calc_frac_cell_masks("eurasia", land_biome_eu)


# Define subregion masks --------------------------------------------------

# Data frame with the region name, subregion name, and min/max longitude values
regions_df <- read.csv("regions.csv")

# This function takes a region mask created in the previous section and applies a new mask
# to only have "1" values within subregion defined by longitude range (xmin, xmax)
create_subregion_mask <- function(region, subregion, xmin, xmax) {
    mask025 <- raster(file.path(mask_dir, paste0("cells_na025_", region, ".tif")))
    # Create a polygon for subregion boundaries
    trim_poly <- st_bbox(c(xmin = xmin, xmax = xmax, ymin = 40, ymax = 80), crs = 4326) %>%
        st_as_sfc() %>%
        st_transform(st_crs(mask025)) %>%
        st_as_sf()
    # Set areas of mask025 outside trim_poly to 0, 
    # then trim raster to remove bordering rows and columns with just 0s
    mask025 <- mask(mask025, trim_poly, updatevalue = 0) %>%
        trim(values = 0)
    writeRaster(mask025, file.path(mask_dir, paste0("cells_na025_", subregion, ".tif")))
}

# Apply function above to all subregions
pwalk(regions_df, create_subregion_mask)
