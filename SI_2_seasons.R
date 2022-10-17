##################################################################################################
# setup
##################################################################################################

#load the prediction workspace
load(file.path("Output", "uk_workspace.rdata"))

pacman::p_load(tidyverse,
               #parallel, 
               future.apply, 
               #mclapply,#; detectCores()
               pls, gstat, sf # UK-PLS MODEL
)
set.seed(1)

image_path <- file.path("..", "Manuscript", "Images")

sim_n <- 30
fewer_hrs_seasons_n <- 12
variable_names <- c("ma200_ir_bc1", "co2_umol_mol", "no2", "pm2.5_ug_m3", "pnc_noscreen")

stops_w <- read_rds(file.path("Output", "stops_used.rda"))

##################################################################################################
# functions
##################################################################################################

# fn takes one random sample from each list item, according to FUN, and calculates a 'value' average
#list should be in wide format (variable/pollutant names) so that same [temporal] samples are collected across pollutants
#list   #sample fn    #descriptor
one_sample_avg <- function(my_list,# = stops_w_list,
                           my_sampling_fn # =s_tow2_sample_fn
) {
  result <- suppressWarnings(
    lapply(my_list, #mc.cores = 6, 
             FUN = my_sampling_fn) %>%
      #unlist results
      bind_rows() %>%
      group_by(location) %>%
      mutate(visits = n()) %>%
      #calculate annual average
      summarize_at(all_of(c(variable_names, "visits")), ~mean(.))
  )

  return(result)
}

validation_stats <- readRDS(file.path("Output", "validation_stats_fn.rda"))

label_pollutant <- readRDS(file.path("Output", "lable_pollutant_fn.rda"))
##################################################################################################
# design sampling
##################################################################################################

##################################################################################################
# SI - look at each two-season pair
message("sampling sites")

two_seasons <- list(
  "spring", "summer", "fall", "winter",
  c("spring", "summer"),
  c("spring", "fall"),
  c("spring", "winter"),
  c("summer", "fall"),
  c("summer", "winter"),
  c("fall", "winter"))

season_lvls <- lapply(two_seasons, function(x) paste(x, collapse = ", ")) %>%
  unlist()

two_seasons_df <- data.frame()

for (i in two_seasons) {
  
  temp <- future_replicate(n = sim_n,
                           simplify = F,
                           expr = one_sample_avg(my_list = group_split(stops_w, location),
                                                 my_sampling_fn = function(x)  {

                                                   df <- filter(x, season %in% i) %>%
                                                     group_by(season) %>%
                                                     # a few sites have < 6 samples/season, so use replacement=True here
                                                     slice_sample(n = fewer_hrs_seasons_n/length(i), replace=T)
                                                 })) %>%
    #unlist
    bind_rows() %>%
    #add simulation number & label
    group_by(location) %>%
    mutate(campaign =  row_number(),
           design = ifelse(length(i)==1, "1 season",
                           ifelse(length(i)==2, "2 seasons", NA)),
           version = paste0(i, collapse = ", ")
    )%>%
    ungroup()

  two_seasons_df <- rbind(two_seasons_df, temp)

}

# add geocovariates to annual estimates & reformat data
two_seasons_df2 <- two_seasons_df %>%
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


# true estimates 
annual <- readRDS(file.path("Output", "annual_training_set2.rda"))
annual_gs_estimates <- annual %>% st_drop_geometry() %>%
  filter(grepl("full", design)) %>%
  distinct(location, variable, value) %>%
  rename(gs_estimate = value)
##################################################################################################
# STANDARD CROSS-VALIDATION
##################################################################################################
message("running models")

folds <- two_seasons_df2 %>%
  distinct(location, design, version, campaign) %>%
  mutate(random_fold = sample(1:k,size = nrow(.), replace = T ))

two_seasons_df2 <- left_join(two_seasons_df2, st_drop_geometry(folds))

## each df in the list is a simulation; running 10-FCV within each simulation to get CV predictions

cv_predictions0 <- lapply(group_split(two_seasons_df2, design, version, campaign, variable), 
                           # mc.cores = use_cores, 
                            FUN = do_cv, fold_name = "random_fold") %>%
  bind_rows() 

cv_predictions <- cv_predictions0 %>%
  select(-var1.var) %>%
  mutate(out_of_sample = "CV",
         #reference = ""
         ) %>%
  rename(campaign_estimate = value) %>%
  # add GS estimates
  left_join(annual_gs_estimates) %>%
  #put back on native scale before evaluating
  mutate_at(vars(contains("estimate"), prediction), ~exp(.)) %>%
  pivot_longer(contains("estimate"), names_to = "reference", values_to = "estimate") %>%
  # some campaigns don't have "estimtes" for test set locations
  drop_na(estimate)

##################################################################################################
# SAVE PREDICTIONS
##################################################################################################
saveRDS(cv_predictions, file.path("Output", "UK Predictions", "specific_season_cv_predictions.rda"))

##################################################################################################
# don't do traditional assessment for spatial clustering - distance analysis
message("calculating performance statistics")

model_perf <- lapply(group_split(cv_predictions, campaign, design, version, variable, out_of_sample, reference), 
                       #mc.cores = 5,
                       validation_stats, prediction = "prediction", reference = "estimate") %>%
  bind_rows()


##################################################################################################
# SAVE DATA
##################################################################################################
select(model_perf , -no_sites) %>%
  saveRDS(., file.path("Output", "specific_season_model_eval.rda"))


message("done with SI_2_seasons.R")



##################################################################################################
# PLOT
##################################################################################################
# ggplot settings
theme_set(theme_bw())
theme_update(legend.position = "bottom",
             # WHY do I have to do this all of a sudden?? plots have a black background otherwise
             plot.background = element_rect(fill = "white", color = NA)
)

model_perf0 <- readRDS(file.path("Output", "specific_season_model_eval.rda"))

model_perf1 <- model_perf0 %>%
  filter(reference == "gs_estimate") %>%
  pivot_longer(cols = c("RMSE", "MSE_based_R2"), names_to = "validation_stat" ) %>%
  mutate(validation_stat = ifelse(validation_stat=="MSE_based_R2", "R2",
                                  ifelse(validation_stat=="RMSE", "RMSE", NA))) %>%
  label_pollutant() %>%
  mutate(
    version = factor(version, levels = season_lvls)
  )

model_perf1 %>%
  filter(validation_stat != "RMSE") %>%
  ggplot(aes(x=version, y=value, fill=variable)) +
  facet_grid(validation_stat~design,  switch = "y", scales = "free") +
  geom_boxplot() +
  labs(x="Sampling Seasons", y="", fill="Pollutant")

ggsave(file.path(image_path, "SI", "r2_by_season.png"), width = 12, height = 5)





model_perf1 %>%
  #filter(validation_stat != "RMSE") %>%
  group_by(variable, design, version, validation_stat) %>%
  summarize(
    Min = min(value),
    Q25 = quantile(value, 0.25),
    Median = median(value),
    Q75 = quantile(value, 0.75),
    Max = max(value)
  ) %>%
  mutate_if(is.numeric, ~round(.,2)) %>%
  select(variable, design, version, Median, validation_stat) %>%
  pivot_wider(names_from = variable, values_from = Median) %>%
  arrange(validation_stat, design, version)
