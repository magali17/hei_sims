# Clear workspace of all objects and unload all extra (non-base) packages
rm(list = ls(all = TRUE))
if (!is.null(sessionInfo()$otherPkgs)) {
  res <- suppressWarnings(
    lapply(paste('package:', names(sessionInfo()$otherPkgs), sep=""),
           detach, character.only=TRUE, unload=TRUE, force=TRUE))
}



pacman::p_load(tidyverse,
               parallel, #mclapply; detectCores()
               #future.apply, #future_replicate()
               pls, 
               gstat, #variogram()
               sf#, #for spatial data; st_
               #sp #has muse datasedt       # ? nEED?
)    

set.seed(1)

# ## for future.apply::future_replicate()
# plan(multisession, workers = 6)
use_cores <- 5

image_path <- file.path("..", "Manuscript", "Images")

##################################################################################################
# LOAD DATA
##################################################################################################
# mapping variables
project_crs <- 4326  #lat/long
m_crs <- 32148

# save coordinate systems as variables
# WGS84 latitude-longitude
#latlong_proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
# Lambert Conic projection (meters)
lambert_proj <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

##################################################################################################
# mm annual estimates
annual <- readRDS(file.path("Output", "annual_training_set.rda")) %>%
  #add covariates
  left_join(readRDS(file.path("Output", "mm_cov_train_set.rda"))) %>%
  #convert to sf
  st_as_sf(coords = c('longitude', 'latitude'), crs=project_crs, remove = F) %>%
  st_transform(m_crs) %>%
  # log transform before modeling
  mutate_at(vars(ma200_ir_bc1:pnc_noscreen), ~log(.)) %>%
  #log-transform UFPs since these are right skewed. all other variables look fine.
  #mutate(pnc_noscreen = log(pnc_noscreen)) %>%
  gather("variable", "value", ma200_ir_bc1:pnc_noscreen) %>%
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

##################################################################################################
# COMMON VARIABLES
##################################################################################################

# for modeling 
cov_names <- st_drop_geometry(annual) %>% ungroup() %>%
  select(log_m_to_a1:last_col()) %>% names() # 224 covarites

pls_comp_n <- 2

#k-folds for CV
k=10

var_names <- unique(annual$variable)

 
##################################################################################################
# CREATE RANDOM VALIDATION FOLDS FOR EACH VARIABLE-DESIGN-VERSION
##################################################################################################

random_fold <- function(df, k.=k) {
  #make sure temporal sims receive same fold designation
  set.seed(2)
  
  result <- df %>% st_drop_geometry() %>%
    distinct(location,
             spatial_temporal, design, version, campaign
    ) %>%
    mutate(random_fold = sample(1:k.,size = nrow(.), replace = T, ))
  
  return(result)
}

############################################################################################################

random_fold_df <- lapply(group_split(annual, spatial_temporal, design, version, campaign),
                         function(x) random_fold(x, k.=k)) %>%
  bind_rows()

#join clusters to annual 
annual <- suppressMessages(left_join(annual, random_fold_df)) %>%
  select(#cluster, 
    random_fold, everything())


############################################################################################################
# UK-PLS FUNCTION 
############################################################################################################
# fn returns UK-PLS predictions. inputs are two spatial objects (simple features). fn automatically transforms these to a lambert projection


# modeling_data = x[[1]]
# new_data = filter(annual_test_set, variable == var_names[i])

uk_pls <- function(modeling_data, # data for fitting pls-uk models
                   new_data, #prediction locations
                   cov_names. = cov_names,  #covariates to be used in modeling
                   pls_comp_n. = pls_comp_n #pls components to use
) {
  
  #lambert projection for UK model
  lambert_proj <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
  ############################################################################################################
  # fit PLS model to estimate fewer components from geocovariates
  
  print(paste(first(modeling_data$version), first(modeling_data$variable), "campaign ", first(modeling_data$campaign) #, "test set ", c(1:10)[!c(1:10) %in% unique(modeling_data$cluster)]
  ))
  
  pls_model <- plsr(as.formula(paste('value ~', paste(cov_names., collapse = "+"))),
                    data = modeling_data, ncomp = pls_comp_n., scale=T, center=T)
  
  #extract compoent scores for UK
  modeling_data_scores <- predict(pls_model, type = "scores") %>% data.frame() %>% #head()
    # add location & value info
    cbind(data.frame(select(modeling_data, -cov_names.))) %>%
    #convert back to sf. geom is dropped otherwise
    st_as_sf()
  
  new_data_scores <- predict(pls_model, type = "scores", newdata = new_data) %>% data.frame() %>%  
    # add location & value info
    cbind(data.frame(select(new_data, -cov_names.))) %>%
    #convert back to sf. geom is dropped otherwise
    st_as_sf()
  
  ############################################################################################################
  # fit UK models & predict at new locations
  
  # UK formula 
  uk_formula <- as.formula(paste("value ~", paste0("Comp.", 1:pls_comp_n., collapse = "+")))
  
  # estimate the variogram model: fit a variogram model, offering to the function several different model options (exponential, spherical, and Matern):
  # using lambert coordinates b/c vertical/horizontal units represent the same jump
  # the default distance in gstat is 1/3 of the maximum distance (use cutoff option to change this)
  
  v.uk <- variogram(uk_formula, st_transform(modeling_data_scores, lambert_proj) )
  m.uk <- fit.variogram(v.uk, vgm(c("Exp", "Sph", "Mat")) )
  #make sure Exp/Sph range estimate is at least 0 when little/no correlation in the data 
  m.uk$range[2] <- max(m.uk$range[2], 1)
  
  #plot(v.uk, m.uk)
  
  # fit UK to the modeling data and predict at the new data locations
  uk_model <- krige(formula = uk_formula, st_transform(modeling_data_scores, lambert_proj), 
                    newdata = st_transform(new_data_scores, lambert_proj), 
                    model = m.uk)
  
  #save predictions
  result <- select(new_data, -cov_names.) %>%
    mutate(prediction = uk_model$var1.pred)
  
  return(result)
  
}

##################################################################################################
# CV function
##################################################################################################
# function returns cross-valited predictions for a given dataset

# x = filter(annual, grepl("full", design), grepl("ma200", variable))
# fold_name = "random_fold"

do_cv <- function (x, fold_name) {
  
  #code to make sure this fn works even if folds don't have all numbers in a sequence (e.g., if too few sites)
  #k = length(unique(x[[fold_name]]))
  k = sort(unique(x[[fold_name]]))
  
  df <- data.frame()
  
  for(f in seq_along(k)) {
    #f=1
    #modeling_data0 = filter(x, !!as.symbol(fold_name) != f)
    modeling_data0 = filter(x, !!as.symbol(fold_name) != k[f])
    #new_data0 = filter(x, !!as.symbol(fold_name) == f)
    new_data0 = filter(x, !!as.symbol(fold_name) == k[f])
    
    temp <- uk_pls(modeling_data = modeling_data0, new_data = new_data0) %>% st_drop_geometry() 
    df <- rbind(df, temp)
  }
  
  return(df)
}

##################################################################################################
# PREDICT
##################################################################################################

common_vars <- c(#"cluster", 
  "location", "route", "visits", "campaign", "design", "version", "spatial_temporal", "variable", #"value",
  "prediction"
)

# # 


##################################################################################################
# GEOCOVARAITE EXTRAPOLLATION - BUFFERED LOO
##################################################################################################
# geocovariate blocks from divisions in principal component ranges

#o	Could weight 2nd dimension by variability: e.g.,	38/total vs 9.7/total
##################################################################################################
# PCA
# prep datasedt 

#only need covariates here
geocov <- annual %>% 
  st_drop_geometry() %>%
  filter(variable == first(variable),
         grepl("full", design))

pca_model <- prcomp(geocov[cov_names], center = T, scale. = T, rank. = 2)
# pca_s <- summary(pca_model)
# prop of variance explained: PC1 = 0.3718, PC2 = 0.09

# pacman::p_load(ggbiplot) 
# ggbiplot(pca_model, labels = dt$location)

#PCA unit equal to a 10% different in geocovariate space
one_pct <- diff(range(pca_model$x[,1]))/100


#pull component score
geocov <- geocov %>%
  mutate(pca1 = pca_model$x[,1]) %>%
  select(location, pca1)

#new dt with PCA comp scores
annual_pca <- annual %>%
  filter(grepl("full", design)) %>%
  left_join(geocov) 

# x =group_split(annual_pca, variable)[[5]] #UFPs
# p1 <- fmsb::percentile(x$pca1)
# grep(pattern = 10, x = p1) # ~ -10

# ggplot(data=x, aes(x=pca1)) + geom_density()

##################################################################################################
# fn returns CV predictions at locations with 10% increment differences in the geocovariate space

# % jump increment
# buff = c(10,30,50,70) [3]
# x = group_split(annual_pca, variable)[[5]] #UFP

pca_dist <- function(x, buff, one_pct. = one_pct) {
  pca_jump = buff*one_pct.
  
  #place to save predictions
  df <- data.frame()
  
  #sites with similar PCA loadings are test site
  for(i in 1:nrow(x) ) {
    # i=2
    # sites with similar PCA loadings
    pred_site <- x[i,]
    # sites with different PCA loadings
    training_sites <- filter(x, pca1 >= pred_site$pca1 + pca_jump | pca1 <= pred_site$pca1 - pca_jump)
     
    
    temp <- uk_pls(modeling_data = training_sites, new_data = pred_site) %>% st_drop_geometry() 
    df <- rbind(df, temp)
  }
  
  df <- df %>%
    mutate(
      spatial_temporal = "spatial",
      design = "covariate distance",
      version = paste0(buff, "%")
    )
  
  return(df)
  
}

############################################################################################
buffers = c(10, #20, 
            30) 


#place to save buffered LOO predictions
pca_dist_loo <- data.frame()

for(b in buffers) {
  # b = buffers[3]
  # x = group_split(annual_pca, variable)[[5]] #UFP
  temp <- mclapply(group_split(annual_pca, variable),
                   mc.cores = use_cores, 
                   FUN = pca_dist, buff = b) %>%
    bind_rows() %>% as.data.frame()
  
  pca_dist_loo <- rbind(pca_dist_loo, temp)
}

pca_dist_loo_predictions <- pca_dist_loo %>% 
  select(common_vars) %>%
  mutate(out_of_sample = "Buffered LOO")  


############################################################################################
print("done with test.R")