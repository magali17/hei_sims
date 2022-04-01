# to run scripts all scripts (1_act_annaulR - 3_model_eval.R scripts) in terminal:
# 1. change directory: cd ~/Documents/Post\ Doc/Study\ Projects/ACT\ TRAP\ MM/ACT\ HEI\ Supp/act_hei_aim1a
# 2. run bash script: bash run_scripts.bash

# to run individual script in terminal:
# 2. enter in terminal: Rscript 1_act_annual.R 

# examples of parallel processing: https://jennybc.github.io/purrr-tutorial/ls12_different-sized-samples.html 

##################################################################################################
# **Purpose** 
#   
#   * We want to better understand how sampling design may impact air pollution predictions models. We use data collected in the ACT TRAP mobile monitoring campaign to conduct sampling simulations and answer questions related to:   
#   - visits per site    
# - stop sampling density 
# - the size of the monitoring area relative to the cohort locations    
# 
# * We noted that prediction models in the California simulations did not always perfom well. Conducting simulations with the ACT TRAP campaign data rather than CA regulatory data should be more revealing since our models perform well


##################################################################################################
# SETUP
##################################################################################################
tictoc::tic()

# Clear workspace of all objects and unload all extra (non-base) packages
rm(list = ls(all = TRUE))
if (!is.null(sessionInfo()$otherPkgs)) {
  res <- suppressWarnings(
    lapply(paste('package:', names(sessionInfo()$otherPkgs), sep=""),
      detach, character.only=TRUE, unload=TRUE, force=TRUE))
}

pacman::p_load(tidyverse,
               parallel, #mclapply; detectCores()
               future.apply, #future_replicate()
               lubridate # %within%
               )    

set.seed(1)

## for future.apply::future_replicate()  
# availableCores() #8 
plan(multisession, workers = 6)
 

##################################################################################################
# LOAD DATA
##################################################################################################
#stop data with temporal variables included
stops_all <- readRDS(file.path("..", "..", "1. Our Campaign", "Our Campaign R", "Data", "Output", "stop_data_win_medians.rda")) %>%
  #drop duplicate UFP instruments
  filter(!variable %in% c("pmdisc_number", "pnc_screen", "ns_total_conc")) %>%
  #don't need these
  select(-c(contains(c("pollutant", "range", "ufp_instrument")))) %>%
  # add route
  mutate(route = as.numeric(gsub(".*R0", "", runname)),
         tow2 = ifelse(day %in% c("Sat", "Sun"), "weekend", "weekday")
         ) 

# only use non-test sites for simulations
non_test_sites <- readRDS(file.path("Output", "mm_cov_train_set_hei.rda")) %>%
  distinct(location) %>% pull()

stops <- filter(stops_all, location %in% non_test_sites)

##########

# stops_all %>% 
#   filter(!location %in% non_test_sites) %>% 
#   
#   group_by(variable, location) %>% 
#   # no. observations per site & pollutant
#   summarize(n=n()) %>% 
#   # distribution
#   summarize(min = min(n), mean=mean(n), max=max(n))




##################################################################################################
# COMMON VARIABLES
##################################################################################################
# number of simulations
sim_n <- 30

rush_hours <- c(7:10, 15:18)
business_hours <- c(9:17)


unique_seasons <- unique(stops$season) %>% as.character()
variable_names <- unique(stops$variable)

# number of samples for fewer hours, seasons, reduced balance

fewer_hrs_seasons_n <- 12

##################################################################################################
# ONLY KEEP STOPS W/ ALL POLLUTANTS
##################################################################################################

keep_times0 <- stops_all %>%
  distinct(time, variable, value) %>%
  pivot_wider(names_from = "variable", values_from = "value") %>% 
  mutate(original_stops = n()) %>%
  drop_na() %>%
  mutate(
    remaining_stops = n(),
    prop_remaining_stops = remaining_stops/original_stops
    ) 

# 9,047 total original stops (309 site x ~26 visits/site)
# after dropping stops w/o all 5 pollutant measures, we have: 8,137 in the training & 810 (see later code) in the test set
# keep_times0

keep_times <- keep_times0 %>% #7345 stops left (90.12%)
  distinct(time) %>% pull()

# range(keep_times) # "2019-03-08 13:51:35 PST" "2020-03-17 23:12:29 PDT"

stops <- filter(stops, time %in% keep_times)

#stops_w <- spread(stops, variable, value)
stops_w <- pivot_wider(data = stops, names_from = "variable",values_from =  "value")

saveRDS(stops_w, file.path("Output", "stops_used.rda"))



### --> repeat for the test set, for consistency






##################################################################################################
# TRUE ANNUAL AVERAGE
##################################################################################################
true_annual <- stops_w %>%
  group_by(location, route) %>%
  mutate(visits = n()) %>%
  summarize_at(vars(variable_names, visits), ~mean(.)) %>%
  mutate(
    campaign = 1,
    design = "full",
    version = "all training data"
  ) %>%
  ungroup()
  
   
# save est set annual averages for validation later
## this also uses the winsorized stop data & the same date ranges
annual_test_set <- stops_all %>%
  filter(!location %in% non_test_sites,
         #same date range as training data; only keep times where all pollutants have values
         time %in% keep_times
         #time >= min(as.Date(keep_times)) & time <= max(as.Date(keep_times))
         ) %>%
  group_by(variable, location) %>%
  summarize(value = mean(value),
            visits = n(),
            campaign = 1,
            design = "test set",
            version = "test set"
  )

saveRDS(annual_test_set, file.path("Output", "annual_test_set.rda"))

##################################################################################################
# SAMPLING DESIGNS
##################################################################################################
# FEWER VISITS

##################################################################################################
## functions 

### sampling functions
s_tow2_sample_fn <- function(df, wkday_visits = 2, wkend_visits = 1) {
  
  no_visits <- ifelse(first(df$tow2) == "weekday", wkday_visits, wkend_visits)
  df0 <- slice_sample(df, n  = no_visits)
  
  return(df0)
  }

# fn takes one random sample from each list item, according to FUN, and calculates a 'value' average
#list should be in wide format (variable/pollutant names) so that same [temporal] samples are collected across pollutants
                             #list   #sample fn    #descriptor
one_sample_avg <- function(my_list,# = stops_w_list, 
                             my_sampling_fn # =s_tow2_sample_fn  
                             ) {
  result <- suppressWarnings( 
    mclapply(my_list, mc.cores = 6, FUN = my_sampling_fn) %>%
    #unlist results 
    bind_rows() %>%
    group_by(location, route) %>%
    mutate(visits = n()) %>%
    #calculate annual average
    summarize_at(all_of(c(variable_names, "visits")), ~mean(.))
  )
    
    return(result)
} 

##################################################################################################
# fewer days

# x = group_split(stops_w, location, tow2)[[2]]

days <- future_replicate(n = sim_n,
                         simplify = F,
                        expr = mclapply(group_split(stops_w, location, tow2), mc.cores = 6, 
                                        FUN = function(x) {slice_sample(x, n= fewer_hrs_seasons_n, replace=T) }) %>%
                          #unlist results 
                          bind_rows() %>%
                          group_by(location, route, tow2) %>%
                          mutate(visits = n(),
                                 version = first(tow2),
                                 design = "fewer days",
                                 ) %>%
                            group_by(location, route, design, version, visits) %>%
                          #calculate annual average
                          summarize_at(all_of(variable_names), ~mean(.)) ) %>%
  bind_rows() %>%
  group_by(location, version) %>%
  mutate(campaign = row_number()) %>%
  #for merging with other sims
  select(names(true_annual)) %>%
  as.data.frame()
 

##################################################################################################
#2. fewer seasons. keep the number of samples the same. 6 is approximately the # of samples/site/season (i.e., max for season ==1)

## notice that a few sites (e.g., MS0138 in winter) have < 6 samples/season, so we'll sample w/ replacement to keep the numebr of samples the same

season_n <- c(1:4)  

season_times <- data.frame()

for (i in seq_along(season_n)) {

  temp <- future_replicate(n = sim_n,  
                           simplify = F,
                                    expr = one_sample_avg(my_list = group_split(stops_w, #variable, 
                                                                                location), 
                                                          my_sampling_fn = function(x)  {
                                                            #? keep to make sure locations pick same season combos across simulations ?? or will lapply do this anyways?
                                                            #set.seed(1)
                                                            seasons_to_sample <- sample(c("spring", "summer", "fall", "winter"), size = season_n[i], replace = F)
                                                            
                                                            df <- filter(x, season %in% seasons_to_sample) %>%
                                                              group_by(season) %>%
                                                              # a few sites have < 6 samples/season, so use replacement=True here
                                                              slice_sample(n = fewer_hrs_seasons_n/season_n[i], replace=T)
                                                            }
                                                          )
                           ) %>%
    #unlist
    bind_rows() %>%
    #add simulation number & label
    group_by(location) %>%
    mutate(campaign =  row_number(),
           design = "balanced seasons",
           version = paste0(season_n[i])
           )%>%
    ungroup()
  
  season_times <- rbind(season_times, temp)

}


##################################################################################################
## randomly select fewer samples for the campaign (easy to simulate but less impractical?)

visit_n <- seq(4,24,4) 

# random_fewer <- list()
# 
# for(i in seq_along(visit_n)) {
#   #i=1
#   temp <- future_replicate(n = sim_n, simplify = F,
#                                    expr = one_sample_avg(my_list = group_split(stops_w, location),
#                                                          my_sampling_fn = function(x) slice_sample(x, n  = visit_n[i]))
#   ) %>%
#     #unlist
#     bind_rows() %>%
#     #add simulation number & label
#     group_by(location, visits) %>%
#     mutate(campaign =  row_number(),
#            design = "fewer visits",
#            version = paste0(visit_n[i], "")
#            )
#   
#   random_fewer[[i]] <- temp
#   
#   }
# 
# random_fewer <- bind_rows(random_fewer) %>% ungroup()

##################################################################################################

rh_bh <- list(sort(unique(c(rush_hours, business_hours))),
              business_hours, rush_hours
              )

names(rh_bh) <- c("business & rush", "business", "rush")

rh_bh_df <- data.frame()

for(i in seq_along(rh_bh)) {
  #i=1
  temp <- future_replicate(n = sim_n,
                           simplify = F,
                           expr = one_sample_avg(my_list =group_split(stops_w, location), 
                           #mc.cores = 4, 
                           my_sampling_fn = function(x) {
                             x %>% filter(tow2 == "weekday",
                                          hour %in% rh_bh[[i]]) %>%
                               # NOTE: using sampling w/ replacement to ensure x samples/site 
                               slice_sample(n=fewer_hrs_seasons_n, replace=T)
                             }
                           )
                           ) %>%
    bind_rows() %>%
    group_by(location) %>%
    mutate(
      campaign = row_number(),
      version = names(rh_bh)[i],
      design = "fewer hours"
    ) %>%
    as.data.frame()
  
  rh_bh_df <- rbind(rh_bh_df, temp)
  
  }


##################################################################################################
# combine TEMPORAL simulation results
temporal_sims <- rbind(
  true_annual,
  days,
  season_times,
  #random_fewer,
  rh_bh_df
  ) %>%
  mutate(spatial_temporal = ifelse(grepl("full", design), "gold standard", "temporal"))


##################################################################################################
## SPATIAL simulations: fewer stops (lower density) 
##################################################################################################
### fewer number of stops/sites 

site_n <- c(25, seq(50, 250, 50))

# fewer_sites <- future_replicate(n = sim_n, simplify = F, 
#                                 expr = mclapply(site_n, mc.cores = 6, FUN = function(x) {
#                                   slice_sample(true_annual, n  = x) %>%
#                                     mutate(version = paste0(x, ""),
#                                            design = "fewer sites",
#                                            )}) ) %>%
#     bind_rows() %>%
#     mutate(campaign = rep(1:sim_n, each=sum(site_n)))
      
##################################################################################################
# ## smaller area/fewer routes
# ### many sources: 1 (urban center, many participants), 6 (airport, Tacoma/suburban/industrial/port)   
# 
# routes_list <- list(c(1:2), c(1:4), c(1:6)) 
# 
# fewer_routes <- mclapply(routes_list, mc.cores = 6, FUN = function(x) {
#   filter(true_annual, route %in% x) %>%
#     mutate(version = paste(range(x), collapse = "-"),
#            design = "fewer routes",
#            campaign = 1)} ) %>%
#   bind_rows()  

##################################################################################################
# # combine spatial simulations
# spatial_sims <- rbind(fewer_sites#, 
#                       #fewer_routes
#                       ) %>%
#   mutate(spatial_temporal = "spatial")



##################################################################################################
# FEWER SITES & VISITS (TOTAL STOPS)
##################################################################################################
# site_n <- c(25, seq(50, 250, 50))
# visit_n <- seq(4,24,4) 

site_n2 <- append(site_n, length(unique(stops_w$location)))
visit_n2 <- append(visit_n, round(mean(true_annual$visits))) 

site_visit_df <- data.frame()

for(v in visit_n2) {
  # v = visit_n2[1]
  temp <- replicate(n = sim_n, simplify = F,
                         expr = mclapply(site_n2, mc.cores = 5, function(x) {
                                                 #sample sites
                                                 sample_sites <- sample(unique(stops_w$location), size = x, replace = F)
                                                 
                                                 #sample visits
                                                 filter(stops_w, location %in% sample_sites) %>%
                                                   group_by(location, route) %>%
                                                   slice_sample(n = v) %>% 
                                                   mutate(visits = n()) %>%
                                                   #calculate annual average
                                                   summarize_at(all_of(c(variable_names, "visits")), ~mean(.)) %>%
                                                   mutate(
                                                     version = paste0(v, "_visits ", x, "_sites")
                                                   )
                                                 })
                         ) %>%
  #unlist
  bind_rows() %>%
    
  ungroup() %>%
  mutate(
    #campaign = rep(1:sim_n, each=s),
    campaign = rep(1:sim_n, each=nrow(.)/sim_n),
    design = "fewer total stops",
    spatial_temporal = "spatial temporal"
  ) 
  
  site_visit_df <- rbind(site_visit_df, temp)
  
  }

site_visit_df <- select(site_visit_df, names(temporal_sims))

# # how many total stops are there in the 278 sites x 26 visits since some sites have a max od 21-25 visits? # 7080
# site_visit_df  %>% 
#   filter(version == "26_visits 278_sites") %>% 
#   distinct(location, visits) %>% 
#   summarize(max_total_stops = sum(visits))

##################################################################################################
# combine spatial and temporal simulations

annual_training_set <- rbind(temporal_sims, #spatial_sims,
                             site_visit_df
                             )

##################################################################################################
# SAVE DATA
##################################################################################################
saveRDS(annual_training_set, file.path("Output", "annual_training_set.rda"))

print("done with 1_act_annual.R")

tictoc::toc()


