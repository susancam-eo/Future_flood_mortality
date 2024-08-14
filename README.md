# Future_flood_mortality
In short: generates estimates of future floods for countries and regions across the globe and produces estimates of number people exposed, deaths, and death rates attributable to flooding under different climate futures

There are four scripts:
1. processing_functions.R. Contains most of the functions that are used multiple times in the other scripts
2. generate_floods_exposures.Rmd. Processes input .tif files of flood return period projections and simulates flood occurrence and calculates number of people exposed to those floods, for each of the next n years, for ~25*25 km grid cells across the globe
3. process_deaths.Rmd. Aggregates from grid floods and exposures to country-level exposures and calculates flood-attributable deaths and death rates for countries, regions, and super-regions
4. compilation_and_figures. For each flood type, produces figures for results of deaths and death rates and compares results under different climate futures
5. figures_by_ssp. For each climate scenario, produces figures that aggregate deaths and death rates across flood types

The scripts require several inputs:
a. A .tif file with flood return period values for a future date for as many flood types and climate scenarios as you would like. Built to work for 10-, 30-, and 100-year floods for SSP 126, SSP 245, and SSP 585
b. A .tif file with gridded world population data for the most recent year
c. A .shp file defining countries/regions of interest
d. A .csv file with projections of population counts and 95% uncertainty intervals for all countries/regions of interest
e. A .csv file with projections of death counts and 95% uncertainty intervals for all countries/regions of interest

To produce correct results, you will need to:
a. Install R and dependencies and create the expected results file structure
b. Specify several settings, including the time window of interest, the initial return period, relative risk values, etc.
c. Check and modify code to ensure that all steps function as expected
