# functions are taken from nox_v2.Rmd

#######################################################################################################################
# functions to prep for PLS  
#######################################################################################################################

## e.g., drop covariates that don't vary enough, etc.

log_transform_distance <- function(all.data, lowerBound=10, removeOrig=FALSE)
{
  distance.vars <- grep("^m_to_", colnames(all.data))
  new.varnames <- c()
  for (i in distance.vars)
  {
    newcol.index <- 1 + ncol(all.data)
    all.data[, newcol.index] <- log( sapply( all.data[, i], function(x) { max(lowerBound, x) } ) )
    colnames(all.data)[newcol.index] <- paste('log_', colnames(all.data)[i], sep='')
  }
  if (removeOrig) all.data <- all.data[, -distance.vars]
  return (all.data)
}


fail_most_common_value <- function(mon.data, vars.all, thres=0.2)
{
  thres <- dim(mon.data[,vars.all])[1]*thres
  
  fail <- apply( mon.data[,vars.all], 2, function(x) {
    tmp <- split(x,x)
    most.common.value <- tmp[[which.max(sapply(tmp, length))]][1]
    sum(x != most.common.value, na.rm=T) } ) < thres
  fail <- names(mon.data[,vars.all])[fail]
  return(fail)
}


fail_low_landuse <- function(mon.data, vars.all, lowValue=10) 
{
  lu.vars <- grep("^rlu_|^lu", vars.all, value=T)
  fail <- sapply(lu.vars, function(x) return (max(mon.data[, x]) < lowValue))
  fail <- names(fail)[grep("TRUE", fail)]
  return (fail)
}



#######################################################################################################################
# fn returns CV PLS predictions for EACH Campaign  
#######################################################################################################################

n_comp <- 2

cv_pls_p <- function(dt,
                     y_val,
                     x_predictors,
                     #CV folds
                     k = 10,
                     # PLS component no.
                     pls_comp = n_comp,
                     label = "" 
) {  
  
  
  dt <- dt %>% rename(y_val = y_val) %>%
    drop_na() %>%
    # place to save predictions
    mutate(cv_p = NA) %>%
    # remove grouping??
    as.data.frame()
  
  # place to save cv predictions
  dt2 <- data.frame()

  #for(p in seq_along(unique_pollutants)) {
    #p=3
    #dt_temp_pollutant <-  dt #%>% filter(Parameter.Name == unique_pollutants[p])
    
    #want the same test set sites across diff designs for any given campaign (note: test sets do change across campaigns, but in the same way across designs. could move this seed into the next for() loop to use same training/test sets every single time)
    set.seed(1)  
    
    # loop through each campaign (e.g., n=30)
    for (i in seq_len(max(dt$Campaign))) {
      #i=1
      
      #1 Campaign at a time
      dt_temp <- dt %>% #dt_temp_pollutant
        filter(Campaign == i)
      
      #create folds for test/train set
      set = sample(c(1:k), size = nrow(dt_temp), replace = T)
      
      # loop through each fold (e.g., k=10)
      for(f in seq_len(k)) {
        #f=1
        
        train_grp <- set != f
        
        dt_temp_train <- dt_temp %>% filter(train_grp)  
        dt_temp_test <- dt_temp %>% filter(!train_grp)   
        
        #fit PLS to training data
        pls_train <- plsr(y_val ~.,
                          data=dt_temp_train[,c("y_val", x_predictors)], 
                          ncomp = pls_comp,
                          scale=T)
        
        #save CV predictions for each test set
        dt_temp$cv_p[!train_grp] <- predict(object = pls_train,
                                            newdata = dt_temp_test
                                            ) %>%
          as.data.frame() %>%
          #make into a vector
          pull()
      }
      
      dt2 <- rbind(dt2, dt_temp) 
      
    }
    
  #}
  
  #change back to original name
  names(dt2)[names(dt2) == "y_val"] <- y_val
  
  #option to relabel prediction variable
  if(label != "") {
    names(dt2)[names(dt2) == "cv_p"] <- label
    }
  
  return(dt2)
  
}



#######################################################################################################################
# fn returns design-sampling-type-verson names (data in long format) as separate columns 
#######################################################################################################################

#str_extract(cv_p_names, "[a-z]+_hr")

#### --> EDIT 

separate_design_type_v <- function(dt, 
                                   var = "Design" #,  designs
) {
  
  dt <- dt %>%
    rename(var = var)
  
  result <- dt %>%
    mutate(
      Samples_per_site = factor(gsub("n", "", str_extract(var, "n[0-9]+"))),
      Version = factor(gsub("v", "", str_extract(var, "v[0-9]"))),
      Season_duration = factor(str_extract(var, "[0-9]wk")),
      Type = ifelse(grepl("true", var), "Long-Term", "Short-Term"),
      Hours = factor(str_to_title(gsub("_hr", "", str_extract(var, "[a-z]+_hr")))),
       
       
    )
  
  #change ver name back
  names(result) [grepl("var", names(result))] <- var
  
  return(result)
}



#######################################################################################################################
# facet_wrap_equal() and facet_grid_equal() functions act like facet_wrap() and facet_grid() in ggplot but it sets the axes ranges (min/max) of each facet to same scale so that the 1-1 line is always down the middle :D !!
#######################################################################################################################
# code source: https://fishandwhistle.net/post/2018/modifying-facet-scales-in-ggplot2/ 

FacetEqualWrap <- ggproto(
  "FacetEqualWrap", FacetWrap,
  
  train_scales = function(self, x_scales, y_scales, layout, data, params) {
    
    # doesn't make sense if there is not an x *and* y scale
    if (is.null(x_scales) || is.null(x_scales)) {
      stop("X and Y scales required for facet_equal_wrap")
    }
    
    # regular training of scales
    ggproto_parent(FacetWrap, self)$train_scales(x_scales, y_scales, layout, data, params)
    
    # switched training of scales (x and y and y on x)
    for (layer_data in data) {
      match_id <- match(layer_data$PANEL, layout$PANEL)
      
      x_vars <- intersect(x_scales[[1]]$aesthetics, names(layer_data))
      y_vars <- intersect(y_scales[[1]]$aesthetics, names(layer_data))
      
      SCALE_X <- layout$SCALE_X[match_id]
      ggplot2:::scale_apply(layer_data, y_vars, "train", SCALE_X, x_scales)
      
      SCALE_Y <- layout$SCALE_Y[match_id]
      ggplot2:::scale_apply(layer_data, x_vars, "train", SCALE_Y, y_scales)
    }
    
  }
)

facet_wrap_equal <- function(...) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  ggproto(NULL, FacetEqualWrap,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}


#same as above but for facet_grid()
FacetEqualGrid <- ggproto(
  "FacetEqualGrid", FacetGrid,
  
  train_scales = function(self, x_scales, y_scales, layout, data, params) {
    
    # doesn't make sense if there is not an x *and* y scale
    if (is.null(x_scales) || is.null(x_scales)) {
      stop("X and Y scales required for facet_equal_wrap")
    }
    
    # regular training of scales
    ggproto_parent(FacetGrid, self)$train_scales(x_scales, y_scales, layout, data, params)
    
    # switched training of scales (x and y and y on x)
    for (layer_data in data) {
      match_id <- match(layer_data$PANEL, layout$PANEL)
      
      x_vars <- intersect(x_scales[[1]]$aesthetics, names(layer_data))
      y_vars <- intersect(y_scales[[1]]$aesthetics, names(layer_data))
      
      SCALE_X <- layout$SCALE_X[match_id]
      ggplot2:::scale_apply(layer_data, y_vars, "train", SCALE_X, x_scales)
      
      SCALE_Y <- layout$SCALE_Y[match_id]
      ggplot2:::scale_apply(layer_data, x_vars, "train", SCALE_Y, y_scales)
    }
    
  }
)

facet_grid_equal <- function(...) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_grid(...)
  
  ggproto(NULL, FacetEqualGrid,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}
 


#######################################################################################################################
# fns calculate basic model performance statistics
#######################################################################################################################


#returns mse_fn
mse_fn <- function(obs, pred){
  mean((obs - pred)^2)
}

rmse_fn <- function(obs, pred){
  sqrt(mean((obs - pred)^2))
}

#returns mse_fn-based R2
r2_mse_based <- function(obs, pred) {
  mse_fn.est <- mse_fn(obs, pred)
  r2 <- 1- mse_fn.est/mean((obs - mean(obs))^2)
  max(0, r2)
}  




#######################################################################################################################
# fn adds sesason to a given dataset with a date variable. Uses typical equinox/solstice dates
#######################################################################################################################


add_season <- function(dt, .date_var) {
  
  pacman::p_load(lubridate)
  
  # dt <- aqs_daily
  # .date_var <- "Date.Local"
  
  winter <- "-12-21" #usually winter starts on 21st, sometimes on 22nd 
  spring <- "-03-20"
  summer <- "-06-21" #usually summer starts on 21st, sometimes on 22nd 
  fall <- "-09-23" #usually fall starts on 22nd, sometimes on 23nd. Using 23rd for 2019 mobile monitoring campaign 
  
  dt <- dt %>%
    rename(date_var = .date_var) %>%
    #make sure variable is in date format
    mutate(date_var = as.Date(date_var),
           season = factor(ifelse((date_var >= ymd(paste0((year(date_var)-1), winter)) & date_var < ymd(paste0(year(date_var), spring))) |
                                    date_var >= ymd(paste0(year(date_var), winter)), "winter",
                                  ifelse(date_var >= ymd(paste0(year(date_var), spring)) &
                                           date_var < ymd(paste0(year(date_var), summer)), "spring",
                                         ifelse(date_var >= ymd(paste0(year(date_var), summer)) &
                                                  date_var < ymd(paste0(year(date_var), fall)), "summer", 
                                                ifelse( date_var >= ymd(paste0(year(date_var), fall)) &
                                                          date_var < ymd(paste0(year(date_var), winter)), "fall", 
                                                        NA)))), 
                           levels = c("spring", "summer", "fall", "winter"))
    )
  
  #change time variable back to what it was originally
  names(dt)[names(dt) %in% "date_var"] <- .date_var
  
  return(dt)
  
}


######################################################################################################################
# * function for a standard boxplot with different whisker definitions to avoid plotting extreme/outlier points
#function takes df and returns summary statistics for plotting alternative boxplots with quantiles: 10, 25, 50, 75 and 90. this reduces the plotting of outliers, which are typically problematic when dealign with large datasets. 
######################################################################################################################

alt_boxplot <- function(df, var, min_q=0.025, max_q=0.975){
  df <- df %>%
    rename(var = var) %>%
    
    #calculate quantiles
    summarize(
      N = n(),
      Min = min(var),
      Qmin = quantile(var, min_q),
      Q25 = quantile(var, 0.25),
      Q50 = quantile(var, 0.50),
      Q75 = quantile(var, 0.75),
      Qmax = quantile(var, max_q),
      Max = max(var)
    )
  
  names(df)[names(df)==var] <- var
  
  return(df) 
  
}


################################################################################################
# fn returns coordinates for a different transformation. it convertes the dataset into a spatial object, calculates coordinates for a diff refernce system, converts these to a df, and attaches these to the original coordinates
# tutorial: https://ryanpeek.org/2017-08-03-converting-XY-data-with-sf-package/ 
################################################################################################
add_crs <- function(
  dt,
  original_crs, original_coords,
  new_crs, new_coord_names = c("long", "lat") 
) {
  
  library(sf)
  
  #convert flat file to spatial file, give it the original CRS
  dt_sp <- st_as_sf(dt, coords = original_coords, crs = original_crs)
  
  #convert to different CRS, and save the coordinates
  new_coords <- st_transform(dt_sp, crs = new_crs) %>%
    st_coordinates() %>%
    as.data.frame() %>%
    rename(new_x = X, 
           new_y = Y) %>%
    mutate(new_crs = new_crs)  
  
  #rename new columns
  names(new_coords) <- c(new_coord_names, paste0(c(new_coord_names, "crs"), collapse = "_"))
  
  # add new coords to original dt
  dt2 <- cbind(dt, new_coords)
  
  
  return(dt2)
  
}

