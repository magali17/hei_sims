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

pacman::p_load(tidyverse,
               parallel, #mclapply; detectCores()
               pls, gstat, sf # UK-PLS MODEL
               )    

#load the prediction workspace
load(file.path("Output", "uk_workspace.rdata"))

set.seed(1)

##################################################################################################
# STANDARD CROSS-VALIDATION
##################################################################################################

## each df in the list is a simulation; running 10-FCV within each simulation to get CV predictions
print(paste0("random ", k, " FCV"))

cv_predictions0 <- mclapply(group_split(annual, spatial_temporal, design, version, campaign, variable), 
                                  mc.cores = use_cores, 
                                  FUN = do_cv, fold_name = "random_fold") %>%
  bind_rows() 

cv_predictions <- cv_predictions0 %>% 
  select(all_of(common_vars)) %>%
  mutate(out_of_sample = "CV") %>%
  
  #note, the few sites & visits design has a handful of sites w/o predicitons - b/c the no. of sites is small(25)?
  ### --> WHY???
  drop_na(prediction)

##################################################################################################
# SAVE PREDICTIONS
##################################################################################################

saveRDS(cv_predictions, file.path("Output", "UK Predictions", "cv_predictions.rda"))

print("done")