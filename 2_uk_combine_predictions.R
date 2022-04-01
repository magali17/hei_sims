##################################################################################################
# setup
##################################################################################################

# Clear workspace of all objects and unload all extra (non-base) packages
rm(list = ls(all = TRUE))
if (!is.null(sessionInfo()$otherPkgs)) {
  res <- suppressWarnings(
    lapply(paste('package:', names(sessionInfo()$otherPkgs), sep=""),
           detach, character.only=TRUE, unload=TRUE, force=TRUE))
}

pacman::p_load(tidyverse, sf)    

set.seed(1)

run_buffered_loo <- TRUE

##################################################################################################
# LOAD DATA & PREDICTIONS
##################################################################################################
# estimates
annual <- readRDS(file.path("Output", "annual_training_set2.rda"))
annual_test_set <- readRDS(file.path("Output", "annual_test_set2.rda"))

#predictions
test_set_predictions <- readRDS(file.path("Output", "UK Predictions", "test_set_predictions.rda"))
cv_predictions <- readRDS(file.path("Output", "UK Predictions", "cv_predictions.rda"))
cluster_predictions <- readRDS(file.path("Output", "UK Predictions", "cluster_predictions.rda"))

if (run_buffered_loo == TRUE) {
  buff_loo_predictions <- readRDS(file.path("Output", "UK Predictions", "buffered_loo_predictions.rda"))
  }

##################################################################################################
# COMBINE PREDICTIONS
##################################################################################################

# add campaign-specific and gold-standard estimates to predictions 
# in UK, can't predict at same places, so can't compare estimates & predictions w/o doing CV

#GS estimates for 278 sites
annual_gs_estimates <- annual %>% st_drop_geometry() %>%
  filter(grepl("full", design)) %>%
  distinct(location, variable, value) %>%
  rename(gs_estimate = value)

#GS estimates for 31 sites
test_gs_estimates <- annual_test_set %>% st_drop_geometry() %>%
  distinct(location, variable, value) %>%
  rename(gs_estimate = value)

#GS estimates for all 309 sites
gs_estimates <- bind_rows(annual_gs_estimates, test_gs_estimates)

# estimates from specific campaign simultaions (n=278 sites)
## note that these only change w/ temporal sims where fewer visits/site are used to estimate annual averages
campaign_estimates <- annual %>% st_drop_geometry() %>%
  distinct(location, design, version, campaign, variable, value) %>% 
  rename(campaign_estimate = value)


# combine predictions & estimates
predictions0 <- rbind(test_set_predictions, cv_predictions) %>%
  #spatial clustered CV
  rbind(cluster_predictions) 

if (run_buffered_loo ==TRUE) {predictions0 <- rbind(predictions0, buff_loo_predictions)}

predictions <- predictions0 %>%
  #left join b/c locations w/ predictions may be fewer than the 309 sites if dno't do 10FCV
  left_join(gs_estimates) %>%
  left_join(campaign_estimates) %>%
  #put back on native scale before evaluating
  mutate_at(vars(contains("estimate"), prediction), ~exp(.)) 


##################################################################################################
# SAVE DATA
##################################################################################################
print("saving predictions")
saveRDS(predictions, file.path("Output", "UK Predictions", "all_predictions.rda"))

# if (run_buffered_loo ==TRUE) {
#   saveRDS(predictions, file.path("Output", "UK Predictions", "all_predictions_with_buffered_loo.rda"))
# }

print("done")