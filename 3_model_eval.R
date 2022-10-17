#script purpose: evaluate UK-PLS model performances

##################################################################################################
# SETUP
##################################################################################################

# Clear workspace of all objects and unload all extra (non-base) packages
rm(list = ls(all = TRUE))
if (!is.null(sessionInfo()$otherPkgs)) {
  res <- suppressWarnings(
    lapply(paste('package:', names(sessionInfo()$otherPkgs), sep=""),
           detach, character.only=TRUE, unload=TRUE, force=TRUE))
}

pacman::p_load(tidyverse,
               parallel, #mclapply; detectCores()
               sf,  
               units # set_units()
)    

set.seed(1)

##################################################################################################
# LOAD DATA
##################################################################################################
# mapping variables
project_crs <- 4326  #lat/long
m_crs <- 32148
 
# uk predictions
predictions <- readRDS(file.path("Output", "UK Predictions", "all_predictions.rda")) %>% 
  #gather("reference", "estimate", contains("estimate")) %>%
  pivot_longer(contains("estimate"), names_to = "reference", values_to = "estimate") %>%
  # some campaigns don't have "estimtes" for test set locations
  drop_na(estimate)

# simulation details
sims <- readRDS(file.path("Output", "annual_training_set.rda")) %>%  
  distinct(location, route, visits, campaign, design, version, spatial_temporal)


#location lat/long 
loc_lat_long <- readRDS(file.path("Output", "location_lat_long.rda")) %>%
  st_as_sf(coords = c('longitude', 'latitude'), crs=project_crs, remove = F) %>%
  st_transform(m_crs)  

# monitoring area shp 
monitoring_region <- readRDS(file.path("..", "..", "1. Our Campaign", "Our Campaign R", "Data", "Output", "GIS", "monitoring_land_shp.rda")) %>%
  st_transform(m_crs)

# spatial clusters
#clusters0 <- readRDS(file.path("Output", "spatial_cluster.rda"))

##################################################################################################
# update datasets

monitoring_area <- st_area(monitoring_region) %>% set_units("km2") %>% drop_units()

##################################################################################################
# CV STATS FUNCTION
##################################################################################################
 
# Fn returns RMSE and MSE-based R2 for a given dataset
validation_stats <- function(dt, prediction, reference){

  # MSE of predictions
  MSE_pred <- mean((dt[[reference]] - dt[[prediction]])^2)
  # MSE of observations (for R2 denominator)
  MSE_obs <- mean((dt[[reference]] - mean(dt[[reference]]))^2)
  
  RMSE = sqrt(MSE_pred)
  MSE_based_R2 = max(1 - MSE_pred/MSE_obs, 0)
  #reg_based_R2 = cor(dt[[reference]], dt[[prediction]], method = "pearson")^2
  
  result <- distinct(dt, campaign, design, version, variable, out_of_sample, reference) %>%
    mutate(
      no_sites = nrow(dt),
      RMSE = RMSE,
      MSE_based_R2 = MSE_based_R2,
      #reg_based_R2 = reg_based_R2
    )
  
  return(result)
  
}

saveRDS(validation_stats, file.path("Output", "validation_stats_fn.rda"))
##################################################################################################
# don't do traditional assessment for spatial clustering - distance analysis
message("calculating performance statistics")

model_perf <- mclapply(group_split(predictions,
                                   campaign, design, version, variable, out_of_sample, reference), 
                       mc.cores = 5,
                       validation_stats, prediction = "prediction", reference = "estimate") %>%
  bind_rows()
 

##################################################################################################
# SAVE DATA
##################################################################################################
select(model_perf , -no_sites) %>%
  saveRDS(., file.path("Output", "model_eval.rda"))


message("done with 3_model_eval.R")

