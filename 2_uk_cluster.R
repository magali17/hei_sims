##################################################################################################
# setup
##################################################################################################

# Clear workspace of all objects and unload all extra (non-base) packages
# rm(list = ls(all = TRUE))
# if (!is.null(sessionInfo()$otherPkgs)) {
#   res <- suppressWarnings(
#     lapply(paste('package:', names(sessionInfo()$otherPkgs), sep=""),
#            detach, character.only=TRUE, unload=TRUE, force=TRUE))
# }

pacman::p_load(tidyverse,
               parallel, #mclapply; detectCores()
               pls, gstat, sf,
               units# UK-PLS MODEL
               )    

#load the prediction workspace
load(file.path("Output", "uk_workspace.rdata"))

set.seed(1)

# number of bins to split the geographic & PCA distance versions into
distance_bins <- 3

##################################################################################################
# SPATIAL CLUSTERS (GEOGRAPHIC DISTANCE)
##################################################################################################
# K-means clustering with 5 clusters of sizes 54, 51, 67, 45, 61
temp0 <- assign_cluster(cluster_covariates = c("lambert_x", "lambert_y")) %>%
  mutate(min_train_to_pred_dist = NA)


# geographic distance between a prediction site and the nearest training site
for(i in 1:nrow(temp0)) {
  # i=1
  # prediction location
  pred_loc <- temp0$location[i] 
  pred_cluster <- temp0$cluster[i]
  
  # sites in a different spatial cluster
  train_locs <- filter(temp0, cluster != pred_cluster) %>%
    pull(location)
  
  # minimum distance between prediction & training locations
  temp0$min_train_to_pred_dist[i] <- round(pt_dist(pred_loc, train_locs)/1000, 1)
}  

saveRDS(temp0, file.path("Output", "spatial_cluster.rda"))

cluster_df <- annual %>%
  filter(grepl("full", design)) %>%
  left_join(temp0) %>%
  mutate(
    version = cut_number(min_train_to_pred_dist, n = distance_bins, dig.lab=1),
    #version = cut(min_train_to_pred_dist, breaks = distance_bins),
    design = "geographic distance 2 (m)",
    out_of_sample = "Spatial CV",
    spatial_temporal = "spatial"
  ) %>%
  st_as_sf()


# # check the location of the clusters
# cluster_df %>%
#   filter(variable == first(variable)) %>%
#   ggplot(aes(x=longitude, y=latitude, col=factor(cluster))) +
#   geom_point()

print(paste0("spatial/clustered ", k, " FCV"))

# predict out-of-cluster
cluster_predictions <- mclapply(group_split(cluster_df, variable), FUN = do_cv, fold_name = "cluster", mc.cores = use_cores
) %>%
  bind_rows() %>% 
  select(all_of(common_vars))  


##################################################################################################
# PCA CLUSTERS (PCA DISTANCE) 
##################################################################################################
#only need covariates here
geocov0 <- annual %>% 
  st_drop_geometry() %>%
  filter(variable == first(variable),
         grepl("full", design))

pca_model <- prcomp(geocov0[cov_names], center = T, scale. = T, rank. = 3)
# pca_s <- summary(pca_model)
# variance explained by each component 
pca_var <- c(0.3755, 0.0959,  0.05857) #sum = 0.52997

# pacman::p_load(ggbiplot)  #devtools
# ggbiplot(pca_model, labels = geocov$location)

#############################################
pca_sf <- geocov0 %>%
  #pull component score
  mutate(pca1 = pca_model$x[,"PC1"],
         pca2 = pca_model$x[,"PC2"],
         pca3 = pca_model$x[,"PC3"],
  ) %>%
  # calculate 3D distance matrix
  st_as_sf(coords = c("pca1", "pca2", "pca3"), remove=F) 

# assign clusters based on PC1-3 (these summarize geocovariates; we could also use "cov_names" here, in principal)
pca_temp0 <- assign_cluster(dt = pca_sf, cluster_covariates = c("pca1", "pca2", "pca3")) %>%
  mutate(min_train_to_pred_dist = NA)

# PCA distance between a prediction site and the nearest training site
for(i in 1:nrow(pca_temp0)) {
  # i=1
  # prediction location
  pred_loc <- pca_temp0$location[i] 
  pred_cluster <- pca_temp0$cluster[i]
  
  # sites in a different spatial cluster
  train_locs <- filter(pca_temp0, cluster != pred_cluster) %>% pull(location)
  
  # minimum distance between prediction & training locations
  pca_temp0$min_train_to_pred_dist[i] <- pt_dist(loc_lat_long. = pca_sf, pred_loc, train_locs, drop_units.=F)
  
}  

saveRDS(pca_temp0, file.path("Output", "pca_cluster.rda"))


pca_df <- annual %>%
  filter(grepl("full", design)) %>%
  left_join(pca_temp0) %>%
  mutate(
    version = cut_number(min_train_to_pred_dist, n = distance_bins, dig.lab=1),
    #version = cut(min_train_to_pred_dist, breaks = distance_bins),
    design = "PCA distance 2",
    out_of_sample = "PCA Spatial CV",
    spatial_temporal = "spatial"
  ) %>%
  st_as_sf()

# # # check the location of the clusters
# pca_df %>%
#   filter(variable == first(variable)) %>%
#   ggplot(aes(x=longitude, y=latitude, col=factor(cluster))) +
#   geom_point()

print(paste0("PCA spatial ", k, " FCV"))

# predict out-of-cluster
pca_cluster_predictions <- mclapply(group_split(pca_df, variable), FUN = do_cv, fold_name = "cluster", 
                                    mc.cores = use_cores
                                    ) %>%
  bind_rows() %>%
  select(all_of(common_vars))  

##################################################################################################
# SAVE PREDICTIONS
##################################################################################################
rbind(cluster_predictions, pca_cluster_predictions) %>%
  mutate(out_of_sample = "Spatial Clusters") %>%
  
  saveRDS(., file.path("Output", "UK Predictions", "cluster_predictions.rda"))

#save version levels
saveRDS(c(levels(cluster_df$version), levels(pca_df$version)), file.path("Output", "cluster_levels.rda"))


print("done")