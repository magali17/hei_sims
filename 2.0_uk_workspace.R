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
               kableExtra, 
               parallel, #mclapply; detectCores()
               #future.apply, #future_replicate()
               pls, 
               gstat, #variogram()
               sf, #for spatial data; st_
               units
)    

set.seed(1) #19

# ## for future.apply::future_replicate()
# plan(multisession, workers = 6)
use_cores <- 5

image_path <- file.path("..", "Manuscript", "Images")


##################################################################################################
# LOAD DATA
##################################################################################################
# mapping variables
# save coordinate systems as variables
project_crs <- 4326  #lat/long
m_crs <- 32148
# Lambert Conic projection (meters)
lambert_proj <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

##################################################################################################
# mm annual estimates
annual <- readRDS(file.path("Output", "annual_training_set.rda")) %>%
  #add covariates
  left_join(readRDS(file.path("Output", "mm_cov_train_set_hei.rda"))) %>%
  #convert to sf
  st_as_sf(coords = c('longitude', 'latitude'), crs=project_crs, remove = F) %>%
  st_transform(m_crs) %>%
  # log transform before modeling
  mutate_at(vars(ma200_ir_bc1:pnc_noscreen), ~log(.)) %>%
  #log-transform UFPs since these are right skewed. all other variables look fine.
  #mutate(pnc_noscreen = log(pnc_noscreen)) %>%
  gather("variable", "value", ma200_ir_bc1:pnc_noscreen) %>%
  #pivot_longer(ma200_ir_bc1:pnc_noscreen,  names_to = "variable", values_to = "value") %>%
  relocate(variable, value, .before = longitude)

## same for test set
annual_test_set <- readRDS(file.path("Output", "annual_test_set.rda")) %>%
  #add covariates
  left_join(readRDS(file.path("Output", "mm_cov_test_set.rda"))) %>%
  #convert to sf
  st_as_sf(coords = c('longitude', 'latitude'), crs=project_crs, remove = F) %>%
  st_transform(m_crs) %>%
  # log transform before modeling
  mutate(value = log(value))

saveRDS(annual_test_set, file.path("Output", "annual_test_set2.rda"))

# ##################################################################################################
# # GRID COVARIATES (for maps)
# ##################################################################################################
# 
# grid_cov0 <- readRDS(file.path("..", "..", "1. Our Campaign", "Our Campaign R", "Data", "Original", "Geocovariates", "dr0311_grid_covars.rda"))
# 
# # load a spatial file of the original monitoring area to assess spatial extrapolation later
# monitoring_area <- readRDS(file.path("..","..", "1. Our Campaign", "Our Campaign R","Data", "Output", "GIS", "monitoring_land_zero_water_shp.rda"))  
# 
# ###########################################################################################
# # GENERATE NEW COVARIATES FOR THE DATASET
# 
# # created some new proximity variables  
# # log transform land proximity variables (e.g., distance to roadways)
# 
# combine_a23_ll <- function(df) {
#   #find buffers for a2-3 length variables
#   buffers <- str_subset(names(df), "ll_a[2:3]") %>% str_extract("s[0:9].*")
#   
#   #for each buffer, calculate sum of a2+a3 length
#   for (i in seq_along(buffers)) {
#     old_vars <- paste0(c("ll_a2_", "ll_a3_"), buffers[i])
#     new_var <- paste0("ll_a23_", buffers[i])
#     
#     df[new_var] <- apply(df[old_vars], 1, sum)
#   }
#   return(df)
# }
# 
# generate_new_vars <- function(df) {
#   # for the NO2 covariate, use the average levels from several available years
#   no2_behr_vars <- c("no2_behr_2005","no2_behr_2006", "no2_behr_2007")
#   
#   df <- df %>%
#     rowwise() %>%
#     mutate(m_to_a123 = min(m_to_a1, m_to_a2, m_to_a3),
#            m_to_a23 = min(m_to_a2, m_to_a3),
#            no2_behr = mean(!!as.symbol(no2_behr_vars))
#     ) %>%
#     ungroup() %>%
#     #make min distance 1 m before log transforming
#     mutate_at(vars(starts_with("m_to_")), ~ifelse(.==0, 1, .) %>% log(.)) %>%
#     rename_at(vars(starts_with("m_to_")), ~gsub("m_to_", "log_m_to_", .)) %>%
#     # calculate sum of a2 and a3 roads in each buffer
#     combine_a23_ll()
# }
# 
# ###########################################################################################
# grid_cov <- generate_new_vars(grid_cov0) 
# 
# ###########################################################################################
# # ADD LOCATION INDICATORS
# 
# # add indicators of whether or not the prediction locations are in the monitoring area
# grid_cov$in_monitoring_area <- suppressMessages(
#   grid_cov %>%
#     st_as_sf(coords = c('longitude', 'latitude'), crs= project_crs) %>%
#     st_intersects(., monitoring_area, sparse = F) %>%
#     apply(., 1, any)
# )
# 
# #convert to sf
# grid_cov <- grid_cov %>%
#   filter(in_monitoring_area==TRUE) %>%
#   st_as_sf(coords = c('longitude', 'latitude'), crs=project_crs, remove = F) %>%
#   st_transform(m_crs)

##################################################################################################
# COMMON VARIABLES
##################################################################################################

# for modeling 
cov_names <- st_drop_geometry(annual) %>% ungroup() %>%
  select(log_m_to_a1:last_col()) %>% names() # 224 covarites

pls_comp_n <- 2

#k-folds for CV
k <- 5 #10

var_names <- unique(annual$variable)

##################################################################################################
# SETUP
##################################################################################################

# # everything look fairly normally distributed after UFPs are log-transformed
# annual %>%
#   filter(grepl("full", design)) %>%
#   ggplot(aes(x=value)) +
#   facet_wrap(~variable, scales="free") +
#   geom_histogram() +
#   labs(title = "Distribution of site annual average concentrations")







##################################################################################################
# CREATE RANDOM VALIDATION FOLDS FOR EACH VARIABLE-DESIGN-VERSION
##################################################################################################

random_fold <- function(df, k.=k) {
  #make sure temporal sims receive same fold designation
  set.seed(2)
  
  result <- df %>% st_drop_geometry() %>%
    distinct(location, spatial_temporal, design, version, campaign) %>%
    mutate(random_fold = sample(1:k.,size = nrow(.), replace = T ))
  
  return(result)
}

############################################################################################################

random_fold_df <- lapply(group_split(annual, spatial_temporal, design, version, campaign),
                         function(x) random_fold(x, k.=k)) %>%
  bind_rows()

#join fold to annual 
annual <- suppressMessages(left_join(annual, random_fold_df)) %>%
  select(random_fold, everything())

saveRDS(annual, file.path("Output", "annual_training_set2.rda"))

#################################################################################################
# FN TO CALCULATE GEOGRAPHIC DISTANCES
##################################################################################################
# fn calculates the distance between 2 sf point datasets & calculates a distance (mean/min/max) summary statistic between two datasets

loc_lat_long <- annual %>%
  filter(grepl("full", design),
         variable == first(variable)) %>%
  distinct(location, longitude, latitude)


pt_dist <- function(locs1, locs2, summary_stat = "min", 
                    #dataset with location info
                    loc_lat_long. = loc_lat_long,
                    drop_units. = TRUE
) {
  #make datasets into sf objects
  locs1 <- filter(loc_lat_long., location %in% locs1)
  locs2 <- filter(loc_lat_long., location %in% locs2)
  
  
  dist_matrix <- st_distance(locs1,locs2) #%>% 
  # drop_units() %>%
  # #replace all "0" distances (for same sites - e.g. training-training pt distances) w/ NA
  # na_if(0)  
  
  #make numeric
  if(drop_units. == TRUE){ 
    dist_matrix <- drop_units(dist_matrix) %>%
      #replace all "0" distances (for same sites - e.g. training-training pt distances) w/ NA
      na_if(0)  
    
  }
  
  #calc distance summary statistic (e.g. mean/min/max) of each dataset
  result <- apply(dist_matrix, 1, summary_stat, na.rm=T) %>% 
    #overall avg, if applicable
    mean()
  
  return(result)
  
}


##################################################################################################
# CLUSTERING FN
##################################################################################################
# fn returns locations assigned to clusters that are based on specific cluster covariates (i.e., lat/long, or geocovariates)

# dt = annual
# cluster_covariates = cov_names

assign_cluster <- function(dt = annual, cluster_covariates, no_centers=k) {
  
  temp <- dt %>% st_drop_geometry() %>% ungroup() %>%
    #variables used for clustering
    distinct(location, !!!syms(cluster_covariates)) %>%
    mutate_at(cluster_covariates, ~scale(.))
  
  set.seed(3)
  clusters <- kmeans(select(temp, all_of(cluster_covariates)), centers = no_centers, algorithm = "Lloyd", iter.max = 30)
  
  ## resulting clusters each have 45-67 sites (range = 22)
  # clusters$size %>% range() %>% diff()
  temp <- temp %>%
    select(location#-cluster_covariates
    ) %>%
    mutate(cluster = clusters$cluster)
  
  return(temp)
}

############################################################################################################
# UK-PLS FUNCTION 
############################################################################################################
# fn returns UK-PLS predictions. inputs are two spatial objects (simple features). fn automatically transforms these to a lambert projection

# modeling_data = group_split(filter(annual, variable == "pnc_noscreen"), spatial_temporal, design, version, campaign)[[1]]
# new_data = filter(annual_test_set, variable == "pnc_noscreen")
# cov_names. = cov_names  #covariates to be used in modeling
# pls_comp_n. = pls_comp_n
# fn_result = "models"
# var_choice = "Exp"
#####
uk_pls <- function(modeling_data, # data for fitting pls-uk models
                   new_data, #prediction locations
                   cov_names. = cov_names,  #covariates to be used in modeling
                   pls_comp_n. = pls_comp_n, #pls components to use
                   fn_result = "predictions", #can be: "predictions" or "models"; return the model predictions or the model fit information
                   # optional: selects the best variogram fit, unless otherwise stated
                   var_choice = "" #Exp
) {
  
  #lambert projection for UK model
  lambert_proj <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
  ############################################################################################################
  # fit PLS model to estimate fewer components from geocovariates
  
  pls_model <- plsr(as.formula(paste('value ~', paste(cov_names., collapse = "+"))),
                    data = modeling_data, ncomp = pls_comp_n., scale=T, center=T)
  
  #extract compoent scores for UK
  modeling_data_scores <- predict(pls_model, type = "scores") %>% data.frame() %>% #head()
    # add location & value info
    cbind(data.frame(select(modeling_data, -all_of(cov_names.)))) %>%
    #convert back to sf. geom is dropped otherwise
    st_as_sf()
  
  new_data_scores <- predict(pls_model, type = "scores", newdata = new_data) %>% data.frame() %>%  
    # add location & value info
    cbind(data.frame(select(new_data, -all_of(cov_names.)))) %>%
    #convert back to sf. geom is dropped otherwise
    st_as_sf()
  
  ############################################################################################################
  # fit UK models & predict at new locations
  
  # UK formula 
  uk_formula <- as.formula(paste("value ~", paste0("Comp.", 1:pls_comp_n., collapse = "+")))
  
  # estimate the variogram model: fit a variogram model, offering to the function several different model options (exponential, spherical, and Matern):
  # using lambert coordinates b/c vertical/horizontal units represent the same jump
  # the default distance in gstat is 1/3 of the maximum distance (use cutoff option to change this)
  v_uk <- variogram(uk_formula, st_transform(modeling_data_scores, lambert_proj) )
  
  #select the best fitting variogram
  #m_uk <- fit.variogram(v_uk, vgm(c("Exp", "Sph", "Mat")) )
  
  if(var_choice == "") {
    m_uk <- fit.variogram(v_uk, vgm(c("Exp", "Sph", "Mat")))
  } else {
    # or select a specific variogram if it's not left blank
    m_uk <- fit.variogram(v_uk, vgm(var_choice))
  }
  
  #make sure Exp/Sph range estimate is at least 0 when little/no correlation in the data 
  m_uk$range[2] <- max(m_uk$range[2], 1)
  
  # fit UK to the modeling data and predict at the new data locations
  uk_model <- krige(formula = uk_formula, st_transform(modeling_data_scores, lambert_proj), 
                    newdata = st_transform(new_data_scores, lambert_proj), 
                    model = m_uk)
  
  #save predictions
  predictions <- select(new_data, -all_of(cov_names.)) %>%
    mutate(prediction = uk_model$var1.pred)
  
  #return(result)
  # return the desired output: either the predictions or the modeling specifications
  if(fn_result == "predictions") {return(predictions)}
  if(fn_result == "models") {
    result = list(
      variable = first(modeling_data$variable),
      pls_model = pls_model, 
      variogram_model = m_uk
    )
    return(result)
  }
  
}

##################################################################################################
# CV function
##################################################################################################
# function returns cross-valited predictions for a given dataset

do_cv <- function (x, fold_name) {
  
  #code to make sure this fn works even if folds don't have all numbers in a sequence (e.g., if too few sites)
  k = sort(unique(x[[fold_name]]))
  
  df <- data.frame()
  
  for(f in k) {
    # f=k[4]
    modeling_data0 = filter(x, !!as.symbol(fold_name) != f)
    new_data0 = filter(x, !!as.symbol(fold_name) == f)
    
    temp <- uk_pls(modeling_data = modeling_data0, new_data = new_data0) %>% st_drop_geometry() 
    df <- rbind(df, temp)
  }
  
  return(df)
}


##################################################################################################
# BUFFERED LOO - GEOGRAPHIC DISTANCE FN
##################################################################################################

# fn returns UK-PLS predictions from buffered leave-one-out validation  
## dt is df for a specific pollutant

buffered_loo <- function(dt, buff) {
  #place to save predictions
  df <- data.frame()
  
  #for each site (row)...
  for (i in 1:nrow(dt) ) {
    #i=1 
    prediction_site <- dt[i,]
    
    #only model with keep sites outside the buffer
    keep_sites <- !st_intersects(st_buffer(prediction_site, buff), dt, sparse = F)[1,]
    modeling_sites <- dt[keep_sites, ]
    
    temp <- uk_pls(modeling_data = modeling_sites, new_data = prediction_site) %>% st_drop_geometry() %>%
      #add number of sites used in model - don't want too little data used
      mutate(no_modeling_sites = sum(keep_sites))
    
    df <- rbind(df, temp) 
  }
  
  df <- df %>%
    mutate(
      spatial_temporal = "spatial",
      design = "geographic distance (m)",
      version = buff
    )
  
  return(df)
  
}
  

##################################################################################################
# BUFFERED LOO - PCA DISTANCE FN
##################################################################################################


# fn returns CV predictions at locations with X increment differences in the geocovariate space
pca_dist <- function(x, buff, distances = pca_3d_distances) {
  
  #place to save predictions
  df <- data.frame()
  
  #sites with similar PCA loadings are test site
  for(i in 1:nrow(x) ) {
    # i=1
    # 1 prediction site w/ AP data & covariates
    pred_site <- x[i,]
    
    # sites with PCA distances beyond the buffer
    training_sites <- distances %>%
      # sites beyond a pca distance from the prediction site
      filter(location == pred_site$location, 
             distance > buff) %>%
      select(-location) %>%
      select(location = location2) %>%
      # join training site observations & covariates
      left_join(x)
    
    temp <- uk_pls(modeling_data = training_sites, new_data = pred_site) %>% st_drop_geometry() %>%
      #add number of sites used in model - don't want too little data used
      mutate(no_modeling_sites = nrow(training_sites))
    
    df <- rbind(df, temp)
  }
  
  df <- df %>%
    mutate(
      spatial_temporal = "spatial",
      design = "PCA covariate distance",
      version = buff  
    )
  
  return(df)
  
}

##################################################################################################
# COMMON VARIABLES
##################################################################################################

common_vars <- c("location", "route", "visits", "campaign", "design", "version", "spatial_temporal", "variable", "prediction")


##################################################################################################
# SAVE WORKSPACE
##################################################################################################

save.image(file.path("Output", "uk_workspace.rdata"))

print("done")
