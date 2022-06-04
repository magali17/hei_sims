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
# don't do traditional assessment for spatial clustering - distance analysis
model_perf <- mclapply(group_split(#filter(predictions, !grepl("Spatial", out_of_sample)),
                                    predictions,
                                   campaign, design, version, variable, out_of_sample, reference), 
                       mc.cores = 5,
                       validation_stats, prediction = "prediction", reference = "estimate") %>%
  bind_rows()
 

##################################################################################################
# SAVE DATA
##################################################################################################
select(model_perf , -no_sites) %>%
  saveRDS(., file.path("Output", "model_eval.rda"))


print("done with 3_model_eval.R")

##################################################################################################
# TEST - SPATIAL CLUSTER VALIDATION APPROACH FOR PREDICTION ERROR BY MONITOR-PREDICTION DISTANCE
##################################################################################################
# TEST- different dataset for distance validation if want to look at exact prediction distance
cluster_perf <- predictions %>%
  filter(grepl("Spatial", out_of_sample)) %>%
  mutate(
    absolute_error = abs(prediction-estimate), # ~RMSE for 1 pt only
    square_error = absolute_error^2, #~MSE
  )

saveRDS(cluster_perf, file.path("Output", "model_eval_cluster_dist.rda"))


v1 <- c("absolute_error", "square_error")

lapply(v1, function(x){
  # x=v1[2]
  df <- cluster_perf %>%
    filter(grepl("Spatial", out_of_sample)) %>%
    pivot_longer(cols = c("absolute_error", "square_error"), names_to = "model_eval", values_to = "value") %>%  
    
    filter(model_eval %in% x) %>%
    mutate(version = as.numeric(version))  
    
      p <- df %>%
        ggplot(., aes(x=version, y=value, col=variable, group=variable)) +
        facet_grid(variable~design, scales = "free", switch = "both") +
        geom_hline(yintercept = 0, linetype=2, alpha=0.5) +
        #geom_smooth(method = "lm") +
        geom_smooth() +
        #geom_point(alpha=0.3) + 
        labs(title = x,
             x = "Distance Between Each Prediction Site & the Closest Monitor"
             )
      
      print(p)
  
  })  %>%
  ggpubr::ggarrange(plotlist = ., common.legend = T, legend = "bottom",
                    ncol = 2 
                    ) %>%
  ggpubr::annotate_figure(top = "Out-of-spatial-cluster prediction error by monitor-prediction distance"
                          )
  
  
############################################################
# TEST: COMPARE VALIDATION APPROACHES: RANDOM VS SPATIAL

# --> NOTE: need to set scales="free" for RMSE

model_perf %>% 
  filter(grepl("Spatial|CV", out_of_sample),
         grepl("full|distance", design),
         reference == "gs_estimate"
         ) %>%  
  mutate(
    design = ifelse(design == "full", "Random", 
                    ifelse(grepl("PCA", design), "Spatially Clustered\nby PCA distance", "Spatially Clustered\nby Geographic\nDistance (m)" ))
    ) %>%  
  ggplot(aes(x="", y=RMSE, col=design)) +
  facet_wrap(~variable, scales="free", 
             nrow = 1) +
  geom_point() + 
  labs(col = "Cross-Validation\nApproach",
       x="")





