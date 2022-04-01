#script purpose: evaluate UK-PLS model performances

##################################################################################################
# to run in terminal:
# 1. change directory:  cd ~/OneDrive\ -\ UW/Documents/Post\ Doc/Study\ Projects/ACT\ TRAP\ MM/ACT\ HEI\ Supp/act_hei_aim1a
# 2. run script:        Rscript 3_model_eval.R 

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

# ## for future.apply::future_replicate()
# plan(multisession, workers = 6)

##################################################################################################
# LOAD DATA
##################################################################################################
# mapping variables
project_crs <- 4326  #lat/long
m_crs <- 32148
 
# uk predictions
predictions <- readRDS(file.path("Output", "UK Predictions", "all_predictions.rda"#, #"predictions.rda"
                                 )) %>% 
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
 
# dt = group_split(predictions, campaign, design, version, variable, out_of_sample, reference)[[3]]
# prediction = "prediction"
# reference = "estimate"

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

##################################################################################################

model_perf <- mclapply(group_split(predictions, campaign, design, version, variable, out_of_sample, reference), 
                       mc.cores = 5,
                       validation_stats, prediction = "prediction", reference = "estimate") %>%
  bind_rows()
 
##################################################################################################
# CALCULATE GEOGRAPHIC DISTANCES
##################################################################################################

# Additionally, for the spatial sims
# Fewer random sites sims: RMSE vs Monitoring density (monitors/km2); avg distance between training monitors; avg distance between training and test (“cohort”) monitors 

# fn calculates the distance between 2 sf point datasets & calculates a distance (mean/min/max) summary statistic between two datasets
## ? idea: it's the closer locations that influence kriging in UK most 

# pt_dist <- function(locs1, locs2, summary_stat = "min", 
#                     #dataset with location info
#                     loc_lat_long. = loc_lat_long
#                     ) {
#   #make datasets into sf objects
#   locs1 <- filter(loc_lat_long., location %in% locs1)
#   locs2 <- filter(loc_lat_long., location %in% locs2)
#   
# 
#   dist_matrix <- st_distance(locs1,locs2) %>%  
#     #make numeric
#     drop_units() %>%
#     #replace all "0" distances (for same sites - e.g. training-training pt distances) w/ NA
#     na_if(0)
#   
#   #calc distance summary statistic (e.g. mean/min/max) of each dataset
#   result <- apply(dist_matrix, 1, summary_stat, na.rm=T) %>% 
#     #overall avg for both datasets
#     mean()
#   
#   return(result)
#   
#   }
# 
# ##################################################################################################
# # calc [mean] distances between training & prediction locations
# 
# train_test_dist <- mclapply(group_split(predictions, campaign, design, version, out_of_sample), mc.cores = 5, function(x) {
# 
#   # prediction locations
#   pred_locs <- unique(x$location)
#   
#   #training locations used to build models for each sim
#   train_locs <- sims %>% 
#     filter(
#       campaign == first(x$campaign),
#       design == first(x$design),
#       version == first(x$version),
#       ) %>%  
#     distinct(location) %>% pull()
#   
#   # calc mean minimum distance between training & prediction locations, for any given campaign
#   ## 1) min distance between each training-test pair; 2) average of all of those
#   dist = pt_dist(train_locs, pred_locs)
#   
#   x %>%
#     distinct(campaign, design, version, spatial_temporal, out_of_sample) %>%
#     mutate(
#       training_sites = length(train_locs),
#       prediction_sites = length(pred_locs),
#       dist = dist)
#   }) %>%
#   
#   bind_rows() %>%
#   mutate(dist = dist/1000) %>%
#   rename(mean_min_train_to_pred_dist_km = dist) %>%
#   #calc monitoring density: 1 monitor every X km2
#   mutate(training_monitors_per_100km2 = training_sites/monitoring_area*100
#   ) 


##################################################################################################
# COMBINE MODEL EVAL & GEOGRAPHIC DISTANCES
##################################################################################################

# train_test_dist, model_perf
# note that only temporal sims have diff estimates for gs_estimate & campaign_estimate
#model_eval <- left_join(train_test_dist, select(model_perf , -no_sites))

##################################################################################################
# SAVE DATA
##################################################################################################
select(model_perf , -no_sites) %>%
  saveRDS(., file.path("Output", "model_eval.rda"))

print("done with 3_model_eval.R")





##################################################################################################
# TEST
##################################################################################################
v1 <- c("RMSE", "MSE_based_R2")

lapply(v1, function(x){
  # x=v1[2]
  #R2
  df <- model_perf %>%
    filter(grepl("Spatial", out_of_sample)) %>%
    pivot_longer(cols = c("RMSE", "MSE_based_R2"), names_to = "model_eval", values_to = "value") %>%  
    
    filter(model_eval %in% x) %>%
    mutate(version = factor(version, levels = readRDS(file.path("Output", "cluster_levels.rda")))  )  
    
    if(x=="RMSE") {
      p <- df %>%
        ggplot(., aes(x=version, y=value, col=variable, group=variable)) +
        facet_wrap(model_eval~design+variable, scales = "free") + 
        geom_point() + 
        geom_line()  
      }
     
  if(x=="MSE_based_R2") {
    p <- df %>%
      ggplot(., aes(x=version, y=value, col=variable, group=variable)) +
      facet_wrap(model_eval~design, scales = "free") + 
      geom_point() + 
      geom_line() 
  }
  print(p)
  
  }) %>%
  ggpubr::ggarrange(plotlist = ., common.legend = T, legend = "bottom",
                    ncol = 1, heights = c(6,4)
                    ) %>%
  ggpubr::annotate_figure(top = 
                            #"Similar distance bins, but a different number of sites per bin"
                            "Similar number of sites per bin but different distance bin ranges"
                          )
  
  
  
