# Boreal fire recurrence study - a randomization approach using MODIS data

The code in this repository was used to study the effect of previous forest fires in limiting future fire occurrence locally, and how that effect varies across the boreal zone. It is based on a randomization approach: a random and independent spatial shift is applied to the fire maps for each year, and summary statistics of fire recurrence (distribution of number of fire years per cell, years between fires in the same cell) are calculated for multiple randomizations as well as the original data. If previous fires limit future fire activity, we expect the original data to show less cells with many fires compared with the randomized versions, and longer time periods between successive fires in the same cell. 

## File list

- `ba_prepare_inputs.R`: This script includes all the steps to produce the input maps for the main analysis script, based on publicly available datasets (500-meter resolution burned area maps and landuse maps from MODIS).

- `ba_rand_test.R`: This script performs the randomization of fire maps and outputs the distribution of years with fire and data on recurring fires in the same cell, for each randomized dataset and the original one.

- `ba_process_outputs`: This script calculates summary statistics and produces analysis graphs based on the resutls of `ba_rand_test.R`.

- `regions.csv`: This file contains the longitudinal boundaries for the different subregions (of the North American and Eurasian boreal regions) to which the analysis is applied.

## Required input data

- MODIS burned area product ([MCD64A1](https://lpdaac.usgs.gov/products/mcd64a1v006/));
- MODIS land cover product ([MCD12Q1](https://lpdaac.usgs.gov/products/mcd12q1v006/)), which is used to filter out water cells from the analysis;
- WWF Terrestrial Ecoregions of the World map ([shapefile here](https://cmerwebmap.cr.usgs.gov/catalog/item/get/508fece8e4b0a1b43c29ca22?files.sort=name&files.order=asc&files.metadataFirst=false)), which is used to filter out non-boreal biomes from the analysis.

## Design choices

- The spatial scale of the random shift must be large enough so that most fires would leave their original footprint, but small enough so that most fires would remain in an area that is similar in terms of fire conditions (climate, vegetation). In other words, we suppose that local feedbacks due to fire history operate at a smaller scale than climate and ecosystem drivers of fire occurrence, and thus aim to choose a randomization scale that "breaks" the former pattern while maintaining the latter. We chose a uniform shift of 20 to 100 cells (~ 10 km to 50 km) with a random direction (+ or -) in both $x$ and $y$ spatial coordinates.

- The spatial randomization process aims to move fires to similar areas nearby, but this assumption applies less near the edge of the study area (for cells bordering lots of water or non-boreal biomes). To reduce these edge effects, we calculate the fraction of water or non-boreal cells in a 201x201 square centered on each cell and eliminate from the analysis cells where this fraction is >25%. That is, fire in those cells can move within the study area, but those cells are not counted in the reported statistics. Despite this adjustment, the total burned area in randomized data is less than the observed one, as some fires from the study area get shifted to water or non-boreal areas and are "lost" from the analysis. Further work could be needed to determine which level of edge exclusion (instead of 25%) would minimize this bias.
