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
# UPLOAD DATA
##################################################################################################

grid_cov0 <- readRDS(file.path("..", "..", "1. Our Campaign", "Our Campaign R", "Data", "Original", "Geocovariates", "dr0311_grid_covars.rda"))

# load a spatial file of the original monitoring area to assess spatial extrapolation later
monitoring_area <- readRDS(file.path("..","..", "1. Our Campaign", "Our Campaign R","Data", "Output", "GIS", "monitoring_land_zero_water_shp.rda"))  

##################################################################################################
# PREP THE GRID
##################################################################################################
# GENERATE NEW COVARIATES FOR THE DATASET

# created some new proximity variables  
# log transform land proximity variables (e.g., distance to roadways)

combine_a23_ll <- function(df) {
  #find buffers for a2-3 length variables
  buffers <- str_subset(names(df), "ll_a[2:3]") %>% str_extract("s[0:9].*")
  
  #for each buffer, calculate sum of a2+a3 length
  for (i in seq_along(buffers)) {
    old_vars <- paste0(c("ll_a2_", "ll_a3_"), buffers[i])
    new_var <- paste0("ll_a23_", buffers[i])
    
    df[new_var] <- apply(df[old_vars], 1, sum)
  }
  return(df)
}

generate_new_vars <- function(df) {
  # for the NO2 covariate, use the average levels from several available years
  no2_behr_vars <- c("no2_behr_2005","no2_behr_2006", "no2_behr_2007")
  
  df <- df %>%
    rowwise() %>%
    mutate(m_to_a123 = min(m_to_a1, m_to_a2, m_to_a3),
           m_to_a23 = min(m_to_a2, m_to_a3),
           no2_behr = mean(!!as.symbol(no2_behr_vars))
    ) %>%
    ungroup() %>%
    #make min distance 1 m before log transforming
    mutate_at(vars(starts_with("m_to_")), ~ifelse(.==0, 1, .) %>% log(.)) %>%
    rename_at(vars(starts_with("m_to_")), ~gsub("m_to_", "log_m_to_", .)) %>%
    # calculate sum of a2 and a3 roads in each buffer
    combine_a23_ll()
}

###########################################################################################
grid_cov <- generate_new_vars(grid_cov0) 

###########################################################################################
# ADD LOCATION INDICATORS

# add indicators of whether or not the prediction locations are in the monitoring area
grid_cov$in_monitoring_area <- suppressMessages(
  grid_cov %>%
    st_as_sf(coords = c('longitude', 'latitude'), crs= project_crs) %>%
    st_intersects(., monitoring_area, sparse = F) %>%
    apply(., 1, any)
)

#convert to sf
grid_cov <- grid_cov %>%
  filter(in_monitoring_area==TRUE) %>%
  st_as_sf(coords = c('longitude', 'latitude'), crs=project_crs, remove = F) %>%
  st_transform(m_crs)


##################################################################################################
# GRID PREDICTIONS
##################################################################################################

# for computational efficiency, select only a few designs to map later
annual_fewer <- annual %>%
  filter(design == "full" |
           (design == "fewer total stops" & version %in% c("4_visits 250_sites", "20_visits 150_sites"#, "20_visits 200_sites"
           )) |
           (design == "balanced seasons" & version %in% c("1", "2")) |
           (design == "fewer days" & version %in% c("weekday")) |
           (design == "fewer hours" & version %in% c("business", "rush"))
  )


print("grid predictions")

# x = group_split(annual_fewer, variable, spatial_temporal, design, version, campaign)[[1]]

grid_predictions0 <- lapply(group_split(annual_fewer, variable, spatial_temporal, design, version, campaign), 
                            #mc.cores = use_cores,
                            function(x) {
                              df = uk_pls(modeling_data = x, new_data = grid_cov) %>%
                                #fn has binding issues later if don't drop geom 
                                st_drop_geometry() %>%
                                #add info to new dataset about the prediction model
                                mutate(
                                  spatial_temporal = first(x$spatial_temporal), 
                                  design = first(x$design), 
                                  version  = first(x$version), 
                                  campaign = first(x$campaign),
                                  variable = first(x$variable)
                                )  
                            }) %>%
  bind_rows() 

# common_vars
grid_predictions <- grid_predictions0 %>% 
  select(location_id, longitude, latitude, spatial_temporal, design, version, campaign, variable, prediction) %>%
  mutate(prediction = exp(prediction))



# # Check that things look OK
# grid_predictions %>%
#   group_by(location_id, longitude, latitude, spatial_temporal, design, version, variable) %>%
#   summarize(prediction_median = median(prediction)) %>%
#   ggplot(aes(x=longitude, y=latitude, col=prediction_median)) +
#   facet_wrap(~variable + design+version) +
#   geom_point()


##################################################################################################
# GRID PREDICTIONS
##################################################################################################

saveRDS(grid_predictions, file.path("Output", "UK Predictions", "grid_predictions.rda"))
