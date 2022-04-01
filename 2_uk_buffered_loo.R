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

pacman::p_load(tidyverse, kableExtra,
               parallel, #mclapply; detectCores()
               pls, gstat, sf # UK-PLS MODEL
               )    

#load the prediction workspace
load(file.path("Output", "uk_workspace.rdata"))

set.seed(1)


##################################################################################################
# GEOGRAPHIC DISTANCE - BUFFERED LOO
##################################################################################################

print("geographic distance")

buffers <- c(200, #closest set of monitors are ~263 m apart
             2e3, 4e3, 8e3,12e3, 16e3
             )

#place to save buffered LOO predictions
buff_loo <- data.frame()

for(b in buffers) {
  
  temp <- mclapply(group_split(filter(annual, grepl("full", design)), variable),
                   mc.cores = use_cores, 
                   FUN = buffered_loo, buff = b) %>%
    bind_rows() 
  
  buff_loo <- rbind(buff_loo, temp)
}


buff_loo_predictions <- select(buff_loo, all_of(common_vars))  

##################################################################################################
# fn returns the proportion of sites kept in each buffer.

no_total_sites <- length(unique(annual$location))

sites_used_fn <- function(dt, no_total_sites. = no_total_sites, combine_n_pct = TRUE) {
  
  df <- dt %>%
    #all variables have same counts
    filter(variable == first(variable)) %>%  
    group_by(version) %>%
    summarize(
      n = n(),
      min = min(no_modeling_sites),
      median = round(median(no_modeling_sites)),
      max = max(no_modeling_sites),
      
      min_pt = min/no_total_sites*100,
      median_pt =median/no_total_sites*100,
      max_pt =max/no_total_sites*100,
      
    ) 
  
  #clean up for table output
  if(combine_n_pct == TRUE) {
    df <- df %>%
      mutate_if(is.numeric, ~round(., 1)) %>%
      mutate(
        min = paste0(min, " (", min_pt, "%)"),
        median = paste0(median, " (", median_pt, "%)"),
        max = paste0(max, " (", max_pt, "%)"),
      ) %>%
      select(-contains("_pt" ), -n)
  }
  
  return(df)
  
}


##################################################################################################
# probably want to keep at least 50% of the data in each fold

# sites_used_fn(dt=buff_loo) %>%
#   kable(caption = "Distribution of the number of sites included in each buffered leave-one-out validation model using geographic distance (N = 278 total sites)",
#         format.args = list(big.mark=","), col.names = c("Buffer (m)", "Min", "Median", "Max") ) %>%
#   kable_styling()

#save for summarizing later
sites_used_fn(dt=buff_loo, combine_n_pct = F) %>%
  saveRDS(file.path("Output", "SI", "sites_used_geo_dist.rda"))


##################################################################################################
# COVARAITE DISTANCE - BUFFERED LOO
##################################################################################################
# geocovariate blocks from divisions in principal component ranges

# PCA
# 1. prep dataset 

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
pca_3d_distances <- geocov0 %>%
  #pull component score
  mutate(pca1 = pca_model$x[,"PC1"],
         pca2 = pca_model$x[,"PC2"],
         pca3 = pca_model$x[,"PC3"],
  ) %>%
  #select(location, contains("pca")) %>%
  # calculate 3D distance matrix
  st_as_sf(coords = c("pca1", "pca2", "pca3")) %>% 
  st_distance() %>% 
  
  #as_tibble(.) %>%
  as.data.frame() %>% 
  
  rename_all(., ~geocov0$location) %>%
  rownames_to_column(var = "location") %>% mutate(location = geocov0$location) %>%
  pivot_longer(geocov0$location, names_to = "location2", values_to = "distance") %>% 
  filter(location != location2) 

#############################################

# range of combined PCA scores
# range(pca_3d_distances$distance) #0.01536253 43.93710358

##################################################################################################
# # example plot 
# pct_val <- 10
# 
# s1 <- geocov$pca1[1]
# buf1 <- data.frame(buff =  c(s1 + one_pct*pct_val, s1 - one_pct*pct_val))
# 
# ggplot(geocov, aes(x=pca1)) + geom_density() + 
#   geom_vline( aes(xintercept = s1, col = "prediction site")) +
#   geom_vline(data=buf1, aes(xintercept = buff, col = paste0(pct_val, "% buffer"))) +
#   labs(x = "PCA PC1 Score",
#        title = paste0("Distribution of PCA PC1 Scores \nwith a ", pct_val, "% Buffer LOO Example"),
#        col = ""
#        )
# 
# ggsave(file.path(image_path, "SI", "pca_dist_loo.png"), width = 5, height = 5)


############################################################################################

print("covariate distance")

# PCA distances
# # check that the selected buffers will have similar training sites in each buffer as the distance buffers 
# pca_dist <- 10
# pca_3d_distances %>% filter(distance >pca_dist) %>% group_by(location) %>% summarize(n = n()) %>% summarize(median(n))

buffers = c(0.5, 1, 2, 5, 7, 9)

#place to save buffered LOO predictions
pca_dist_loo <- data.frame()

for(b in buffers) {
  temp <- mclapply(group_split(filter(annual, grepl("full", design)), variable),
                   mc.cores = 3, 
                   FUN = pca_dist, buff = b) %>%
    bind_rows() %>% as.data.frame()
  
  pca_dist_loo <- rbind(pca_dist_loo, temp)
}


# # TEST - Check proportion of sites in each buffer
# pca_dist_loo %>%
#   filter(variable == first(variable)) %>%
#   group_by(version) %>%
#   summarize(min = min(no_modeling_sites),
#             median=median(no_modeling_sites),
#             max = max(no_modeling_sites)
#   )



pca_dist_loo_predictions <- select(pca_dist_loo, all_of(common_vars))

##################################################################################################
# check the proportion of sites kept in each buffer.
# probably want to keep at least 50% of the data in each fold

# sites_used_fn(dt=pca_dist_loo) %>% 
#   kable(caption = "Distribution of the number of sites included in each buffered leave-one-out validation model using PCA distance (N = 278 total sites)",
#         format.args = list(big.mark=","), col.names = c("Buffer", "Min", "Median", "Max") ) %>%
#   kable_styling()

#save for summarizing later
sites_used_fn(dt=pca_dist_loo, combine_n_pct = F) %>%
  saveRDS(file.path("Output", "SI", "sites_used_pca_dist.rda"))
  

##################################################################################################
# SAVE PREDICTIONS
##################################################################################################
rbind(buff_loo_predictions, pca_dist_loo_predictions) %>%
  mutate(out_of_sample = "Buffered LOO")  %>%
  saveRDS(., file.path("Output", "UK Predictions", "buffered_loo_predictions.rda"))

print("done")
