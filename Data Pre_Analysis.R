#******************************************************************************#
#* This script was used to analyze changes in bird counts in China before, during and after the pandemic 
#* Bird counts are compared before and during the pandemic, and during and after the pandemic, respectively. 
#* It consists of the following three main parts:
#*     1. eBird data cleaning
#*     2. GLMM
#*     3. Data visualization
#* Note: The following parameters are set to compare before and during the pandemic
#******************************************************************************#

# Libraries --------------------------------------------------------------------
library(tidyverse) # for general operation
library(auk) # package for manipulating eBird data
library(lubridate) # for dealing with dates
library(dggridR) # for creating the hexagonal grid for the spatial subsample
library(sf) # for various gis functions
library(foreach)
library(doParallel)
library(dplyr)
select <- dplyr::select 

# Directories ------------------------------------------------------------------
# Data directory for project datasets
data_dir <- "D:/RWork/COVID_19/data"

# Subset full data -------------------------------------------------------------

# file paths to store the project ebd and sed files
f_ebd <- file.path(data_dir, "CN_ebd_2020-09.txt")
f_sed <- file.path(data_dir, "CN_sed_2020-09.txt")

# columns to retain from original eBird data:
keep_cols <- c("sampling_event_identifier", "observer_id", "group_identifier", 
               "category", "scientific_name", "subspecies_scientific_name", 
               "observation_count", 
               "latitude", "longitude", 
               "country_code", "state_code", "county_code", 
               "locality_id", "locality_type", 
               "protocol_type", "protocol_code", 
               "observation_date", "time_observations_started", 
               "duration_minutes", "effort_distance_km", "effort_area_ha", 
               "number_observers", 
               "all_species_reported")

# Run the extraction (only if the project file doesn't already exist)
if (!file.exists(f_ebd)) {
  auk_ebd("D:/RWork/eBird_Work/data_CN/ebird/ebd_CN_relApr-2023.txt", 
          "D:/RWork/eBird_Work/data_CN/ebird/ebd_sampling_relApr-2023.txt") %>%
    # Only the CN:
    auk_country(country = c("CN")) %>% 
    # Mar-May, 2017-2020
    auk_year(year = 2017:2020) %>% 
    auk_date(date = c("*-03-01", "*-05-31")) %>% 
    # < 5 h, < 5 km, traveling and stationary counts
    auk_protocol(protocol = c("Stationary", "Traveling")) %>% 
    auk_distance(distance = c(0, 5)) %>% 
    auk_duration(duration = c(0, 5 * 60)) %>% 
    # Only complete checklists:
    auk_complete() %>% 
    # Run the filter:
    auk_filter(f_ebd, file_sampling = f_sed, keep = keep_cols)
}

# Taxonomy data ----------------------------------------------------------------
tax <- dplyr::select(ebird_taxonomy, scientific_name, species_code)
glimpse(tax)
head(tax)
# write_csv(tax_select, path = "data/taxonomy.csv")
# write_csv(tax_select, path = "data/taxonomy_selected.csv")

# Load filtered EBD and SED data -----------------------------------------------
ebd <- read_ebd(f_ebd, unique = FALSE)
sed <- read_sampling(f_sed, unique = FALSE)

# add "county code" ------
# Assign a county code to each record, since the Chinese ebird data does not have this attribute.
library(raster)
library(dplyr)

# ebd
raster_data <- raster("D:/RWork/eBird_Work/data_CN/china_boundary/city_big_code_sin_toRaster.tif") # tif file of county code in China
point_coords <- ebd %>%
  dplyr::select(longitude, latitude) 
point_sp <- SpatialPoints(point_coords, proj4string=CRS("+proj=longlat +datum=WGS84")) 
point_sp_transformed <- spTransform(point_sp, CRS(proj4string(raster_data))) 
point_cells <- cellFromXY(raster_data, point_sp_transformed) 
point_values <- raster_data[point_cells] 
ebd <- cbind(ebd, city_code = point_values) 
ebd_county <- ebd %>% 
  filter(!is.na(city_code)) 
glimpse(ebd_county)
# write.csv(ebd, "ebd_city_code_removeNA.csv", row.names = FALSE) 

# sed
raster_data <- raster("D:/RWork/eBird_Work/data_CN/china_boundary/city_big_code_sin_toRaster.tif") 
point_coords <- sed %>%
  dplyr::select(longitude, latitude) 
point_sp <- SpatialPoints(point_coords, proj4string=CRS("+proj=longlat +datum=WGS84")) 
point_sp_transformed <- spTransform(point_sp, CRS(proj4string(raster_data))) 
point_cells <- cellFromXY(raster_data, point_sp_transformed) 
point_values <- raster_data[point_cells] 
sed <- cbind(sed, city_code = point_values) 
sed_county <- sed %>% 
  filter(!is.na(city_code)) 
glimpse(sed_county)
# write.csv(sed, "sed_city_code_removeNA.csv", row.names = FALSE) 

# Observer experience ------
# Function to check if a user is 'new' in 2020:
# (Users are considered 'new' if they did not submit any checklists from the UK during March-June of 2017-2019)
check_status <- function(x) {
  yr <- lubridate::year(x)
  if (any(yr < 2020)) {
    if (any(yr == 2020)) {
      "continuing"
    } else {
      "old"
    }
  } else {
    return("new")
  }
}

# identify observer status
obs_status <- sed_county %>% 
  group_by(observer_id) %>% 
  summarise(observer_status = check_status(observation_date)) %>% 
  ungroup() %>% 
  mutate(experienced_obs = (observer_status != "new"))
obs_status

# identify checklist status (in relation to whether observers are new)
checklist_status <- sed_county %>% 
  mutate(checklist_id = coalesce(group_identifier, 
                                 sampling_event_identifier)) %>% 
  inner_join(obs_status, by = "observer_id") %>% 
  group_by(checklist_id) %>% 
  summarise(experienced_obs = any(experienced_obs)) %>% 
  ungroup()
checklist_status
# status_counts <- table(obs_status$observer_status) 
# status_counts

# Remove duplicate checklists (those shared among multiple users) and any without experienced observers
ebd_unique <- auk_unique(ebd_county)
ebd_unique
sed_unique <- auk_unique(sed_county, checklists_only = TRUE) %>% 
  inner_join(checklist_status, by = "checklist_id") 
  # inner_join(checklist_status, by = "checklist_id") %>%
  # filter(experienced_obs)
sed_unique

# zero fill--------
zf <- auk_zerofill(ebd_unique, sed_unique, collapse = TRUE)
zf
glimpse(zf)

# convert scientific names to codes------
zf <- zf %>% 
  inner_join(tax, by = "scientific_name") %>% 
  select(checklist_id, species_code, species_observed, observation_count, city_code)
zf

# prevalence------
prevalence <- checklist_status %>% 
  filter(experienced_obs) %>%  
  semi_join(zf, ., by = "checklist_id") %>% 
  group_by(city_code, species_code) %>% 
  summarise(n_checklists = n(),
            n_pos_obs = sum(species_observed),
            prevalence = mean(species_observed)) %>% 
  ungroup() %>% 
  mutate(county_code = city_code) %>% 
  select(county_code, species_code,
         n_checklists, n_pos_obs, prevalence)
prevalence

# prevalence_total-------
prevalence_total <- prevalence %>% 
  group_by(species_code) %>% 
  summarise(n_checklists = sum(n_checklists),
            n_pos_obs = sum(n_pos_obs)) %>% 
  ungroup() %>% 
  mutate(county_code = "total",  
         prevalence = n_pos_obs / n_checklists) %>% 
  select(any_of(names(prevalence)))
prevalence_total

prevalence <- prevalence %>% 
  mutate(county_code = as.character(county_code)) 
bind_rows(prevalence_total, prevalence) %>% 
  write_csv("data/prevalence.csv")

# taxonomy
ebird_taxonomy %>% 
  filter(category == "species") %>% 
  select(species_code, common_name, scientific_name, family) %>% 
  inner_join(prevalence_total %>% select(-county_code), by = "species_code") %>% 
  write_csv("data/ebird-taxonomy.csv")


# convert x to -1, drop zeros ------------
zf_drop0 <- zf %>% 
  mutate(observation_count = if_else(observation_count == "X", -1L, 
                                     as.integer(observation_count))) %>% 
  filter(species_observed) %>% 
  dplyr::select(checklist_id, species_code, observation_count)
zf_drop0
glimpse(zf_drop0)
# write.csv(zf_drop0, "data/observations.csv", row.names = FALSE)

# Subsetting checklists ---------------------------------------------------

# Spatial functions
library(dggridR) # Creating Hexagonal grid
library(sf) # for gis functions

# Remove any checklists with
sed_subset <- sed_unique %>% 
  # Removing a couple of grouped checklists that evaded the first round
  distinct(checklist_id, .keep_all = T) %>% 
  # give zero distance to stationary counts
  mutate(effort_distance_km = if_else(protocol_type == "Stationary",0, 
                                      effort_distance_km)) %>%
  # Additional filters
  filter(number_observers <= 10,
         effort_distance_km <= 5, 
         experienced_obs == T) %>%
  filter(locality_type == "P" | locality_type == "H") %>%
  # Add year column
  mutate(observation_year = lubridate::year(observation_date)) %>%
  # "Pandemic" column indicating Pre (FALSE) or During (TRUE)
  mutate(pandemic = (observation_year == 2020))
glimpse(sed_subset)

# Spatial subsample -------------------------------------------------------

# Creates a hexagonal grid and samples an equal number of checklists (up to a pre-defined max) from each cell.
# First, generate hexagonal grid with ~ 3 km between cells
dggs <- dgconstruct(spacing = 3)

# For each checklist, get the gridcell ID, and add that to the data
checklists_f <- sed_subset %>% 
  mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum)
checklists_f

# Function to subsample. This randomly selects an equal number of checklists from pre- and post-pandemic, 
# up to a maximum threshold in each county (thresh). It selects them in pairs from each gridcell in the county, 
# looping through gridcells until the threshold is reached (or until every gridcell in the county no longer has paired pre- and post- checklists to choose from.
SpatialSubsample <- function(dat, thresh) {
  # Object for full selection:
  select <- NULL
  
  # Run selection algorithm
  for (i in 1:length(unique(dat$city_code))) {
    
    # Selection from county
    cty_select <- NULL
    
    # Current county
    cty <- unique(dat$city_code)[i]
    
    # Grid cells covered
    cells <- dat %>% filter(city_code == cty) %>% distinct(cell)
    
    # Split data into pre and post pandemic
    pre <- dat %>% filter(city_code == cty,
                          pandemic == FALSE)
    post <- dat %>% filter(city_code == cty,
                           pandemic == TRUE)
    
    # Step through cells to sample
    while(length(cty_select) < thresh &
          intersect(
            (pre %>% group_by(cell, .drop = FALSE) %>% count())$cell,
            (post %>% group_by(cell, .drop = FALSE) %>% count())$cell
          ) %>% length() > 0) {
      
      # Loop through each cell
      for (j in 1:nrow(cells)) {
        
        # Current cell
        cl <- cells$cell[j]
        
        # Sample
        if(any(pre$cell == cl) &
           any(post$cell == cl) &
           length(cty_select) < thresh){
          # Sample one from each
          s1 <- sample((pre %>% filter(cell == cl))$checklist_id, size = 1)
          s2 <- sample((post %>% filter(cell == cl))$checklist_id, size = 1)
          
          # Add to county selection
          cty_select <- c(cty_select, s1, s2)
          
          # Remove from "supply"
          pre <- pre[-which(pre$checklist_id == s1),]
          post <- post[-which(post$checklist_id == s2),]
        }
      }
    }
    
    # Add to full selection
    select <- c(select, cty_select)
  }
  
  return(select)
  
}

# Set the random seed, so that this next step is reproducible:
set.seed(22)
# tt <- proc.time() # For run time, later
selection <- SpatialSubsample(dat = checklists_f, thresh = 1000)
# run_time <- (proc.time()[3] - tt[3]) / 60 # Run time in minutes
glimpse(selection)

checklists_sample <- checklists_f %>% filter(checklist_id %in% selection)
checklists_sample
glimpse(checklists_sample)

checklists_sample$county_code <- as.character(checklists_sample$city_code)
checklists_sample
# write.csv(checklists_sample, "data/checklists_1000.csv", row.names = FALSE)

# checklist data, subset columns, get distinct locations------------------------
locations_sub <- checklists_sample %>% 
  mutate(county_code = city_code) %>% 
  dplyr::select(latitude, longitude, country_code, state_code, county_code, locality_id, locality_type) %>% 
  distinct(locality_id, .keep_all = T)
glimpse(locations_sub)

# Merge objects
locations_sub <- merge(locations_sub, location_count_sub, by = "locality_id")
glimpse(locations_sub)
# Save list
write_csv(locations_sub, path = "data/locations_subset_1000.csv") # for adding 'Distance to airport', 'Distance to road' in ArcMAP

# extract land cover info of each checklist ------------------------------------
library(sf)
library(raster)

#
GAIA <- raster("data/GAIA_city26_Reclassify_sin.tif")
#
locations_sub_sf <- st_read("data/locations_subset_1000_sin_buffer100.shp") %>% 
  # st_geometry() %>% 
  st_transform(crs = projection(GAIA))
locations_sub_sf
#
GAIA_checklists <- exact_extract(GAIA, locations_sub_sf, progress = TRUE) %>% 
  map_dfr(~tibble(urban = mean(.$value, na.rm = TRUE))) %>% 
  bind_cols(locations_sub_sf)
GAIA_checklists
glimpse(GAIA_checklists)
#
# write_csv(GAIA_checklists, path = "data/GAIA_checklists.csv")
# I convert 0-0.05 to 1 (rural), 0.05-1 to 0 (urban) in ArcMAP


# Extract data for focal species -----------------------------------------------

# setting the threshold:
# drop observations from counties with lower than this prevalence
prevalence_thresh <- 0.05

# drop counties with less than this number of sampled checklists
check_thresh <- 100 

# input data
# prevalence
prevalence <- read.csv("data/prevalence.csv") 
glimpse(prevalence)

prevalence_select <- prevalence %>% 
  filter(county_code != "total",
         prevalence > prevalence_thresh)
glimpse(prevalence_select)

# Get list of counties to keep 
keep_cty <- read.csv("data/checklists_1000.csv") %>% 
  count(county_code) %>% 
  filter(n > 0)
glimpse(keep_cty)

# checklist data
checklists <- read.csv("data/checklists_1000.csv")
glimpse(checklists)
checklists$county_code <- as.character(checklists$county_code) 
glimpse(checklists)

checklists <- checklists %>% 
  # drop counties below threshold
  semi_join(prevalence, by = "county_code") %>%
  # drop counties with poor coverage
  filter(county_code %in% keep_cty$county_code) 
glimpse(checklists)

# counts
observations <- read.csv("data/observations.csv") %>%
  # filter(species_code == )
  mutate(species_observed = TRUE)
glimpse(observations)

# zero fill
zf2 <- left_join(checklists, observations, by = "checklist_id") %>% 
  # deal with X/-1 and fill zeros
  mutate(species_observed = coalesce(species_observed, FALSE),
         # observation_count = coalesce(observation, 0L),
         observation_count = if_else(observation_count == -1, NA_integer_,
                                     observation_count))
glimpse(zf2)


# # Function to Obtain and zero-fill a single species data product:
SpeciesPrep_1000 <- function(species) {
  # define species
  taxonomy <- here("data", "ebird-taxonomy.csv") %>%
    read_csv() %>%
    filter(common_name == species)

  # input data
  # prevalence
  prevalence <- here("data", "prevalence.csv") %>%
    read_csv() %>%
    filter(species_code == taxonomy$species_code,
           county_code != "total",
           prevalence > prevalence_thresh)

  # Get list of counties to keep
  keep_cty <- here("data", "checklists_1000_v2.rds") %>%
    read_rds() %>%
    count(county_code) %>%
    filter(n >= check_thresh)

  # checklist data
  checklists <- here("data", "checklists_1000_v2.rds") %>%
    read_rds() %>%
    # drop counties below threshold
    semi_join(prevalence, by = "county_code") %>%
    # drop counties with poor coverage
    filter(county_code %in% keep_cty$county_code) 

  # counts
  observations <- here("data", "observations.rds") %>%
    read_rds() %>%
    # subset to species
    filter(species_code == taxonomy$species_code) %>%
    mutate(species_observed = TRUE)

  # zero fill
  zf <- left_join(checklists, observations, by = "checklist_id") %>%
    # deal with X/-1 and fill zeros
    mutate(species_observed = coalesce(species_observed, FALSE),
           observation_count = coalesce(observation_count, 0L),
           observation_count = if_else(observation_count == -1, NA_integer_,
                                       observation_count),
           species_code = coalesce(species_code, taxonomy$species_code))

  return(zf)
}


# Model data prep --------------------------------------------------------------

library(tidyverse)
library(lubridate)

# Directory for model output:
mod_dir <- "D:/RWork/COVID_19/model_output"

# Load checklist metadata 
checklists <- read.csv("data/checklists_1000.csv")
glimpse(checklists)

# Time variables 
# Day of Year
checklists <- checklists %>% mutate(day_of_year = lubridate::yday(observation_date))
# Month
checklists <- checklists %>% mutate(month = lubridate::month(observation_date, label = FALSE))
glimpse(checklists)

# Load covariates
# Location-specific data (extracted from GIS)
library(readxl)
loc_cov_data <- read_excel("data/locations_subset_1000_sin_buffer100_Air_Road_Land.xls")
glimpse(loc_cov_data)
loc_cov_data

# confirmed cases in each county
cases_count <- read.csv("D:/RWork/COVID_19/data/cases_count.csv")
cases_count$cityCode <- as.double(cases_count$cityCode)
cases_count <- cases_count %>% 
  mutate(county_cod = cityCode) %>% 
  select(county_cod, confirmed_sum, confirmed_mean)
glimpse(cases_count)
head(cases_count)
#
loc_cov_data <- loc_cov_data %>% 
  inner_join(cases_count, by = "county_cod")
glimpse(loc_cov_data)


# Distance to Road/Airport/Stringency
# Add centered, logged distance to major road and distance airport
# Negative = closer than average to road/airport,
# Positive = further than average to road/airport
checklists2 <- loc_cov_data %>% 
  mutate(delta_traffic = confirmed_mean*-1) %>% 
  mutate(c_delta_traffic = scale(delta_traffic, scale = FALSE)) %>% 
  mutate(air = scale(log(Dis_Air), scale = FALSE)) %>% 
  mutate(road = scale(log(Dis_Road), scale = FALSE)) %>% 
  select(locality_id, c_delta_traffic, air, road) %>% 
  left_join(checklists, ., by = "locality_id")

# Developed/undeveloped, based on 50m buffer of landcover
# Rounded and subtracted from 1 so that:
# 0 = mostly developed (i.e. 'urban')
# 1 = mostly undeveloped (i.e. 'rural')
checklists3 <- loc_cov_data %>% 
  select(locality_id, land) %>% # I convert 0-0.05 to 1 (rural), 0.05-1 to 0 (urban) in ArcMAP
  left_join(checklists2, ., by = "locality_id")
glimpse(checklists3)


# Prep Function 

checklists3_select <- checklists3 %>% 
  select(checklist_id, day_of_year, c_delta_traffic, air, road, land)
glimpse(checklists3_select)

zf3 <- zf2 %>% 
  filter(!is.na(observation_count)) %>% 
  # Convert Pandemic variable
  # -1 = Pre(2017-19); 0 = Post(2020)
  mutate(pandemic = as.numeric(pandemic) - 1) 
  # Add other data
glimpse(zf3)

# join the checklists and obsercations
zf4 <- left_join(zf3, checklists3_select, by = "checklist_id")
glimpse(zf4)

# Remove any NA landcover checklists (should only be from one bad location)
zf4_remove <- zf4 %>% 
  filter(!is.na(land))
glimpse(zf4_remove)

# Return only the data needed for the model
dat <- zf4_remove %>% 
  mutate(sp = species_code) %>% 
  mutate(y = observation_count) %>% 
  mutate(cty = county_code) %>% 
  mutate(delta_traf = c_delta_traffic) %>% 
  mutate(dist = effort_distance_km) %>% 
  mutate(dur = duration_minutes) %>% 
  select(sp, checklist_id, day_of_year,
         y, cty, delta_traf,
         pandemic, dist, dur, road, air, land)
head(dat)
glimpse(dat)

#
delta_traf <- dat$delta_traf
road <- dat$road
air <- dat$air
delta_traf <- apply(delta_traf, 1, unlist)
road <- apply(road, 1, unlist)
air <- apply(air, 1, unlist)

dat_remove <- dat[, !colnames(dat) %in% c("delta_traf","road","air")]
dat_full <- cbind(dat_remove, delta_traf, road, air)
head(dat_full)
glimpse(dat_full)
# write_csv(dat_full, path = "CN_dat_full_delta_traf_centered.csv")



# species filtering
species_counts <- dat_full %>% 
  group_by(sp) %>% 
  summarise(count = n()) %>% 
  ungroup()
species_counts
# write_csv(species_counts, path = "species_counts.csv")

species_select <- species_counts %>% 
  filter(count >= 100) %>% 
  pull(sp)
species_select

#
dat_full_filter <- dat_full %>% 
  filter(sp %in% species_select)
head(dat_full_filter)
glimpse(dat_full_filter)
# write_csv(dat_full_filter, path = "CN_dat_full_filter_delta_traf_centered.csv")

# Running Model -----------------------------------------------------------

setwd("D:/RWork/COVID_19")
getwd()


# libraries
library(dplyr)
library(parallel)
num_cores <- detectCores()
cat("Number of CPU cores:", num_cores, "\n")

# Set directories 
out_dir <- "D:/RWork/COVID_19/model_output_V4"

# species
species_unique <- read.csv("CN_species_family_order_unique.csv")
species_unique

species_list <- as.list(species_unique$species_code)
species_list
glimpse(species_list)

# species loop
for (species in species_list) {
  # Start error sink 
  errors <- file(paste("D:/RWork/COVID_19/model_output/errors_", species, ".txt", sep = ""), open="wt")
  sink(errors, type="message")
  # Stan settings
  options(mc.cores = 5)
  rstan::rstan_options(auto_write = TRUE)
  # Data 
  dat <- read.csv("CN_dat_full_filter.csv") %>%
    filter(sp == species)
  tt <- proc.time() # For run time, later
  # Run Model 
  fit <- rstanarm::stan_glmer(y ~ (1|cty) + dist + dur +
                                pandemic + 
                                air + air:pandemic +
                                road + road:pandemic +
                                land + land:pandemic,
                              data = dat,
                              family = rstanarm::neg_binomial_2, 
                              chains = 5, 
                              warmup = 300, 
                              iter = 1300, 
                              thin = 1, 
                              cores = 5, 
                              seed = 121212,
                              refresh = 1) 
  run_time <- (proc.time()[3] - tt[3]) / 60 
  # Save Output:
  saveRDS(fit, file = file.path(out_dir, paste("fit_", species, ".rds", sep = "")))
  saveRDS(run_time, file = file.path(out_dir, paste("time_", species, ".rds", sep = "")))
  # Finish error sink 
  sink(type="message")
  close(errors)
  #
  print(species)
  
}


# Interpreting results ---------------------------------------------------------
setwd("D:/RWork/COVID_19")
getwd()
# Libraries 
library(tidyverse)
library(MCMCvis)
library(bayesplot)

# Set directories 
mod_dir <- "D:/RWork/COVID_19/model_output_V4"



# Species List 
sp_list <- read.csv("CN_species_family_order_unique.csv") %>% 
  mutate(band_code = species_code) %>% 
  dplyr::select(common_name, scientific_name, band_code)
sp_list

sp_list <- sp_list %>% 
  filter(band_code != "martit2")
sp_list

# Parameter posteriors to summarize
pars <- c("dist",
          "dur",
          "pandemic",
          "delta_traf",
          "air",
          "road",
          "land",
          "pandemic:delta_traf",
          "pandemic:air",
          "pandemic:road",
          "pandemic:land",
          "reciprocal_dispersion")

# Summary stats desired for each posterior
coln <- c("mean",
          "sd",
          "lower",
          "median",
          "upper",
          "Rhat",
          "ESS")

# Create object to store results:
res_name_list <- vector(length = length(coln)*length(pars) + 3)
res_name_list[1:3] <- c("Species_common", "Species_sci", "Band_code")
for(i in 1:length(pars)) {
  for (j in 1:length(coln)) {
    res_name_list[3 + length(coln) * (i - 1) + j] <- paste0(pars[i],"-",coln[j])
  }
}
res <- matrix(nrow = nrow(sp_list), ncol = length(res_name_list))
dimnames(res) <- list(
  NULL,
  res_name_list
)
res[,1:3] <- dplyr::select(sp_list, common_name, scientific_name, band_code) %>% as.matrix()
res

# Fill each row with appropriate species results
for (i in 1:nrow(res)) {
  species <- res[i,"Band_code"]
  
  # Extract species results
  fit <- readRDS(file = file.path(mod_dir, paste("fit_", species, ".rds", sep = "")))
  
  # Fill row with summary stats, rounded
  res[i,4:ncol(res)] <- MCMCsummary(fit, params = pars, HPD = FALSE) %>%
    as.matrix() %>% t() %>% as.vector() %>% signif(digits = 4)
  #
  print(i)
}

post_summary <- as.data.frame(res)

# # Save as csv
# write_csv(post_summary, path = file.path(mod_dir, "post_summary_V4.csv"))
# # reload
# post_summary <- read_csv(file = file.path(mod_dir, "post_summary_V4.csv"))
# post_summary

# 
sp_list_family <- read.csv("CN_species_family_order_unique.csv") %>% 
  mutate(Band_code = species_code) %>% 
  dplyr::select(common_name, family, Band_code) 
  # filter(Band_code != "martit2")
sp_list_family

post_summary2 <- inner_join(post_summary, sp_list_family, by = "Band_code")
glimpse(post_summary2)


# Function to plot all interaction CIs and medians, ordered by family and common name---------

PostSumBW <- function(trm) {
  # Extract data
  dat <- post_summary2 %>%
    dplyr::select(Species_common, family, Band_code,
                  paste("pandemic:", trm, "-median", sep = ""),
                  paste("pandemic:", trm, "-upper", sep = ""),
                  paste("pandemic:", trm, "-lower", sep = "")) %>%
    rename(common_name = paste("Species_common"), 
           family = paste("family"),
           band_code = Band_code,
           median = paste("pandemic:", trm, "-median", sep = ""),
           upper = paste("pandemic:", trm, "-upper", sep = ""),
           lower = paste("pandemic:", trm, "-lower", sep = "")) %>%
    # arrange(median)   
    arrange(desc(family), desc(common_name))   
  
  # Conversion to the correct sign (due to reversing axis direction for all variables other than 'overlap') # ignore this
  if (trm == "overlap") {
    dat <- dat
  } else {
    dat <- dat %>%
      mutate(median_old = median,
             lower_old = lower,
             upper_old = upper,
             median = median_old * -1,
             lower = upper_old * -1,
             upper = lower_old * -1) %>%
      # arrange(median) 
      arrange(desc(family), desc(common_name))
  }
  
  # Plotting info:
  y <- -(1:nrow(dat))
  xlims <- c(min(dat$lower), max(dat$upper))
  
  # Margin parameters
  par(mar = c(0.1, 10, 3.1, 0.1), 
      mgp = c(2, 0.5, 0)) 
  
  # Blank plot
  plot(0,
       xaxt = 'n',
       yaxt = 'n',
       bty = 'n',
       pch = '',
       ylab = '', xlab = '',
       xlim = xlims, ylim = c(-1,-length(y))
  )
  
  # X axis
  axis(side = 3, tck = -0.01, cex.axis = 0.75)
  mtext(text = paste("pandemic:", trm, " interaction", sep = ""), side = 3, line = 2)
  
  # Species Guidelines
  for (i in 1:nrow(dat)) {
    lines(x = xlims,
          y = c(y[i], y[i]),
          lty = 3,
          col = c("#bdbdbd"))
  }
  
  # Zero Line
  abline(v = 0, lty = 2)
  
  # Data
  for (i in 1:nrow(dat)) {
    points(x = dat$median[i], y = y[i],
           pch = 16, cex = 0.75,
           col = if((dat$upper[i] > 0 & dat$lower[i] > 0) |
                    (dat$upper[i] < 0 & dat$lower[i] < 0)) {
             c("#005a32")} else {
               c("#74c476")
             })
  }
  
  
  for (i in 1:nrow(dat)) {
    lines(x = c(dat$lower[i], dat$upper[i]),
          y = c(y[i], y[i]),
          lwd = 1.5,
          col = if((dat$upper[i] > 0 & dat$lower[i] > 0) |
                   (dat$upper[i] < 0 & dat$lower[i] < 0)) {
            c("#005a32")} else {
              c("#74c476")
            })
  }
  
  # Species names
  par(mgp = c(0, 0, 0))
  for (i in 1:nrow(dat)) {
    axis(labels = paste(dat$common_name[i], "(", dat$family[i], ")", sep = ""),
         at = y[i],
         side = 2, las = 2, tick = FALSE,
         cex.axis = 0.6,
         col = if((dat$upper[i] > 0 & dat$lower[i] > 0) |
                  (dat$upper[i] < 0 & dat$lower[i] < 0)) {
           c("#005a32")} else {
             c("#74c476")
           },
         font = if((dat$upper[i] > 0 & dat$lower[i] > 0) |
                   (dat$upper[i] < 0 & dat$lower[i] < 0)) {
           2} else {1}
    )
  }
  
}


# plot and export

# stringency
png("D:/RWork/COVID_19/Figure_PostSum_Before_During/delta_traf_order.png", width=6, height=10, res=600, units="in") 
PostSumBW(trm = "delta_traf")
dev.off()

# air
png("D:/RWork/COVID_19/Figure_PostSum_Before_During/air_order.png", width=6, height=10, res=600, units="in") 
PostSumBW(trm = "air")
dev.off()

# road
png("D:/RWork/COVID_19/Figure_PostSum_Before_During/road_order.png", width=6, height=10, res=600, units="in") 
PostSumBW(trm = "road")
dev.off()

# land
png("D:/RWork/COVID_19/Figure_PostSum_Before_During/land_order.png", width=6, height=10, res=600, units="in") 
PostSumBW(trm = "land")
dev.off()

# Marginal effects/Interaction Plots -------------------------------------------

### Function to make the species plots (all panels)
SpeciesIntPlot <- function(species) {
  # Load results
  fit <- readRDS(file = file.path(mod_dir, paste("fit_", species, ".rds", sep = "")))
  
  # Panel arrangement 
  par(mfrow = c(3,2), 
      mar = c(5, 5, 1, 1), 
      mgp = c(2, 0.3, 0), 
      las = 1) 
  
  # Custom functions for floor and ceiling 
  floor_dec <- function(x, level = 1) round(x - 5 * 10^(-level - 1), level) 
  ceiling_dec <- function(x, level = 1) round(x + 5  *10^(-level - 1), level) 
  
  # Function to predict and plot each variable
  IntPlot <- function(trm) {
    
    # Decide whether to show the plot:
    if ( # if the posterior does not occur all above or all below zero (i.e. overlaps zero)
      !(((post_summary %>% filter(Band_code == species) %>%
          dplyr::select(paste("pandemic:", trm, "-lower", sep = "")) %>%
          pull > 0) &
         (post_summary %>% filter(Band_code == species) %>%
          dplyr::select(paste("pandemic:", trm, "-upper", sep = "")) %>%
          pull > 0)) |
        ((post_summary %>% filter(Band_code == species) %>%
          dplyr::select(paste("pandemic:", trm, "-lower", sep = "")) %>%
          pull < 0) &
         (post_summary %>% filter(Band_code == species) %>%
          dplyr::select(paste("pandemic:", trm, "-upper", sep = "")) %>%
          pull < 0)))
    ) {
      # Plot for 'non-significant' variable
      plot(0, 
           xaxt = 'n', yaxt = 'n', bty = 'n',
           pch = '',
           ylab = '', xlab = '',
           xlim = c(-1,1), ylim = c(-1,1))
      
      text(x = 0, y = 0,
           labels = paste("posterior distribution 95% CI \n for interaction with\n*",
                          trm, "* includes zero", sep = ""),
           cex = 1)
      
    } else { # if the posterior did not overlap zero
      
      # Set limits and resolution for ggpredict (auto limits are not great)
      if (trm == "land") {
        mn <- 0
        mx <- 1
      } else {
        mn <- fit$x[,trm] %>% min() %>% floor_dec(level = 1)
        mx <- fit$x[,trm] %>% max() %>% ceiling_dec(level = 1)
      }
      # tms <- c(paste(trm, " [", mn, ":", mx, ", by = 1]", sep = ""), "pandemic") #--
      if (trm == "delta_traf") {
        tms <- c(paste(trm, " [", mn, ":", mx, ", by = 10]", sep = ""), "pandemic")
      } else {
        if (trm == "land") {
          tms <- c(paste(trm, " [", mn, ":", mx, ", by = 1]", sep = ""), "pandemic")
        } else {tms <- c(paste(trm, " [", mn, ":", mx, ", by = 0.1]", sep = ""), "pandemic")}
      }
      
      
      # Run predict function:
      temp <- ggeffects::ggpredict(fit, tms)
      
      temp <- as.data.frame(temp)
      
      ### Values to plot
      # Median bird values
      y_plot <- temp %>% filter(group == 0) %>% dplyr::select(predicted) %>%
        as.data.frame() %>% rename(post = predicted)
      y_plot <- temp %>% filter(group == -1) %>% dplyr::select(predicted) %>%
        as.data.frame() %>% rename(pre = predicted) %>% add_column(y_plot)
      
      # CI (upper)
      y_up_plot <- temp %>% filter(group == 0) %>% dplyr::select(conf.high) %>%
        as.data.frame() %>% rename(post = conf.high)
      y_up_plot <- temp %>% filter(group == -1) %>% dplyr::select(conf.high) %>%
        as.data.frame() %>% rename(pre = conf.high) %>% add_column(y_up_plot)
      
      # CI (lower)
      y_low_plot <- temp %>% filter(group == 0) %>% dplyr::select(conf.low) %>%
        as.data.frame() %>% rename(post = conf.low)
      y_low_plot <- temp %>% filter(group == -1) %>% dplyr::select(conf.low) %>%
        as.data.frame() %>% rename(pre = conf.low) %>% add_column(y_low_plot)
      
      # Axes limits and values
      ylims <- c(min(y_low_plot) - min(y_low_plot) * 0.1,
                 max(y_up_plot) + max(y_up_plot) * 0.1)
      
      # Predictor Values
      x <- temp %>% filter(group == 0) %>% dplyr::select(x) %>%
        as.data.frame() %>% rename(post = x)
      x <- temp %>% filter(group == -1) %>% dplyr::select(x) %>%
        as.data.frame() %>% rename(pre = x) %>% add_column(x)
      
      # X values for plotting
      if (trm == "overlap") {
        x_plot <- x
      } else {
        x_plot <- -x
      }
      xlims <- c(
        min(x_plot) %>% floor,
        max(x_plot) %>% ceiling
      )
      
      # Other X axis inputs:
      if (trm == "overlap") {
        x_lab_pos <- c(-1, -0.5, 0) 
        x_lab_val <- c(0, 50, 100) 
        x_lab_txt <- "Percent Overlap" 
      } else {
        if (trm == "delta_traf") {
          
          x_lab_val <- seq(from = -1000, to = 0, by = 100)
          x_lab_pos <- (x_lab_val - -748.1376) * -1
          
          x_lab_txt <- "Stringency of lockdown"
        } else {
          if (trm == "land") {
            x_lab_pos <- c(-1, 0)
            x_lab_val <- c("Rural", "Urban")
            x_lab_txt <- "Landcover"
          } else {
            if (trm == "air") {
              
              x_lab_val <- c(1, 10, 100, 500)
              x_lab_pos <- (c(1, 10, 100, 500) * 1000) %>% log() %>% -10.26713 * -1
              x_lab_txt <- "Distance to Airport (km)"
            } else {
              
              x_lab_val <- c(0.001, 0.01, 0.1, 1, 10, 100)
              x_lab_pos <- (c(0.001, 0.01, 0.1, 1, 10, 100) * 1000) %>% log() %>% -5.491703 * -1
              x_lab_txt <- "Distance to Road (km)"
              
            }
          }
        }
      }
      
      # plot 
      if (trm == "land") {
        ### Draw panel
        # Draw blank plot with specific axes
        plot(0,
             xaxt = 'n',
             pch = '',
             ylab = '', xlab = '',
             xlim = xlims, ylim = ylims,
             cex.axis = 1, tck = -0.01)
        title(ylab = "Expected Count", line = 2, cex.lab = 1)
        title(xlab = x_lab_txt, line = 2, cex.lab = 1)
        axis(side = 1, 
             at = x_lab_pos, 
             labels = x_lab_val, 
             cex.axis = 1,
             tck = -0.01, 
             mgp = c(0.9, 0, 0),
             padj = 0.3
        )
        
        # Add pre-CI---
        segments(x_plot$pre, y_low_plot$pre,
                 x_plot$pre, y_up_plot$pre,
                 col = adjustcolor("#01665e", alpha.f = 0.8),
                 lwd = 3)
        # Add pre-point
        points(x = x_plot$pre,
               y = y_plot$pre,
               type = "p",
               pch = 16,
               col = adjustcolor("#01665e", alpha.f = 0.8),
               cex = 3)
        # Add pre-line
        lines(x = x_plot$pre,
              y = y_plot$pre,
              lty = 2,
              lwd = 3, 
              col = adjustcolor("#01665e", alpha.f = 0.8))
        
        # Add post-CI---
        segments(x_plot$post+0.01, y_low_plot$post,
                 x_plot$post+0.01, y_up_plot$post,
                 col = adjustcolor("#8c510a", alpha.f = 0.8),
                 lwd = 3)
        # Add post-point
        points(x = x_plot$post+0.01,
               y = y_plot$post+0.01,
               type = "p",
               pch = 16,
               col = adjustcolor("#8c510a", alpha.f = 0.8),
               cex = 3)
        # Add post-line
        lines(x = x_plot$post+0.01,
              y = y_plot$post+0.01,
              lty = 2,
              lwd = 3, 
              col = adjustcolor("#8c510a", alpha.f = 0.8))
      } else {
        ### Draw panel
        # Draw blank plot with specific axes
        plot(0,
             xaxt = 'n',
             pch = '',
             ylab = '', xlab = '',
             xlim = xlims, ylim = ylims,
             cex.axis = 1, tck = -0.01)
        title(ylab = "Expected Count", line = 2, cex.lab = 1)
        title(xlab = x_lab_txt, line = 2, cex.lab = 1)
        axis(side = 1,
             at = x_lab_pos,
             labels = x_lab_val,
             cex.axis = 1, 
             tck = -0.01, 
             mgp = c(0.9, 0, 0),
             padj = 0.3
        )
        
        # Add pre-CI 
        polygon(x = c(x_plot$pre, rev(x_plot$pre)),
                y = c(y_low_plot$pre, rev(y_up_plot$pre)),
                border = NA, col = adjustcolor("#01665e", alpha.f = 0.1))
        # Add pre-line
        points(x = x_plot$pre,
               y = y_plot$pre,
               type = "l",
               col = adjustcolor("#01665e", alpha.f = 0.8),
               lwd = 3)
        
        
        # Add post-CI 
        polygon(x = c(x_plot$post, rev(x_plot$post)),
                y = c(y_low_plot$post, rev(y_up_plot$post)),
                border = NA, col = adjustcolor("#8c510a", alpha.f = 0.1))
        # Add post-line
        points(x = x_plot$post,
               y = y_plot$post,
               type = "l",
               col = adjustcolor("#8c510a", alpha.f = 0.8),
               lwd = 3)
      }
      
    }
  }
  
  
  # Loop through all five variables
  vars <- c(
    "delta_traf",
    "land",
    "air",
    "road")
  
  for (i in 1:length(vars)) {
    IntPlot(trm = vars[i])
  }
  
  
  # Plot name and legend 
  plot(0,
       xaxt = 'n', yaxt = 'n', bty = 'n', 
       pch = '',
       ylab = '', xlab = '',
       xlim = c(-2,2), ylim = c(-1.5,1.5))
  
  
  text(x = 0, y = 1,
       labels = sp_list %>% filter(band_code == species) %>%
         dplyr::select(common_name) %>% paste(),
       cex = 1.5 
  )
  text(x = 0, y = 0.5, font = 3, 
       labels = sp_list %>% filter(band_code == species) %>%
         dplyr::select(scientific_name) %>% paste(),
       cex = 1) 
  
  legend(x = 0, y = -0.5, ncol = 1, cex = 1.5,
         legend = c("Before pandemic (2017-2019)", "During pandemic (2020)"), # 
         col = c("#01665eCC","#8c510aCC"),  
         lwd = 3, 
         xjust = 0.5, yjust = 0.5) 
  
}

# # for test
# png("D:/RWork/COVID_19/Figure_MarginalEffect_Before_During/plot.png", width=8, height=11, res=600, units="in") 
# SpeciesIntPlot(species = "bkhgul")
# dev.off()

sp_list
species_l <- sp_list$band_code
species_l
for (species in species_l) {
  png(paste0("D:/RWork/COVID_19/Figure_MarginalEffect_Before_During/", species, ".png"), width=8, height=11, res=600, units="in")
  SpeciesIntPlot(species = species)
  dev.off()
  print(species)
}

