
#Function to filter date between 01/05 and 31/10
filter_date <- function(df){
  if ("Date" %in% names(df)) {df <- df %>% rename(date = Date)}
  res <- df %>%
    filter(month(date) >=6, month(date) <= 10)
  return(res)
}

### Function to transform reconstruct_history() to the right format
# One row per station_year, one column by day, NA if not working, 0 if working
# June 1st -> October 31st
formatting_history <- function (df) {
  fichier_problemes <- df %>%
    tidyr::pivot_longer(
      cols = -date,
      names_to = "station",
      values_to = "val") %>%
    mutate(
      annee = year(date),
      jour_julien = yday(date)-yday(as.Date(paste0(annee,"-06-01")))+1,
      val = ifelse(val == 1, 0, NA)  # 1 = no issue → 0, 0 = issue → NA
    )  %>%
    mutate(jour_julien=factor(jour_julien,levels=1:153)) %>% # first day june 1st, 153rd day october 31st
    mutate(station_annee = paste0(station, "_", annee)) %>%
    group_by(station_annee, jour_julien) %>%
    summarise(val = first(val), .groups = "drop") %>%  
    tidyr::pivot_wider(
      names_from = jour_julien,
      values_from = val,
      values_fill = NA   # missing days → NA
    ) %>%
    select(
      station_annee,
      all_of(
        colnames(.)[setdiff(seq_along(colnames(.)), 1)]  # indices without station_year
        [order(as.numeric(colnames(.)[-1]))]             
      )
    ) %>%
    as.data.frame()
  
  rownames(fichier_problemes) <- fichier_problemes$station_annee
  fichier_problemes$station_annee <- NULL
  
  return(fichier_problemes)
}

####################################################################
############## Calculation of occupancy probabilities ############## 

### 1. Creating occupancy files 

# Function to filter independent events based on species-specific treshold
filter_temporal_independence <- function(df,treshold){
  
  df <- df %>%
    mutate(date_time = as.POSIXct(paste(as.character(date), hour), 
                                  format = "%Y-%m-%d %H:%M:%S"))
  df <- df %>%
    arrange(station, date_time) %>%
    group_by(station) %>%
    filter(row_number() == 1 | difftime(date_time, dplyr::lag(date_time), units = "mins") >= treshold) %>%
    ungroup()
  
  return(df)
}


# Function which creates first occupancy file for each species. One row per station_annee, one column per day, not taking history into account
fichier_detection <- function (species_clean_filtered){ 
  animal <- species_clean_filtered %>% 
    mutate(
      annee = year(date),
      mois = month(date),
      jour_julien = yday(date)-yday(as.Date(paste0(annee,"-06-01")))+1,
      station=factor(station,levels=unique(species_clean_filtered$station))) %>%
    mutate(jour_julien=factor(jour_julien, levels = 1:153),
           station_annee = paste0(station,"_",annee))
  
  b<-as.data.frame(table(animal$station_annee,animal$jour_julien)) %>%
    tidyr::pivot_wider(names_from = Var2,
                       values_from = Freq,
                       values_fill = 0) %>%
    rename(station_annee=Var1) %>%
    as.data.frame()
  rownames(b) <- b$station_annee
  b$station_annee <- NULL
  return(b)
} 
# Function which adds camera history (days with malfunction) to the occupancy file. One column per day, adding cameras histories.
merging_detection_problemes <- function(species_clean_filtered,camera_history) {
  animal <- fichier_detection(species_clean_filtered)
  stopifnot(all(colnames(animal) == colnames(camera_history)))
  
  for (row in rownames(camera_history)){
    if (row %in% rownames(animal)) {
      camera_history[row,]<-camera_history[row,]+animal[row,]
    }
  }
  
  return(camera_history)
}

# Function which aggregate occupancy file according to the survey temporal window. Return the occupancy file and the sampling effort file (nb of working days for each station_year, and each survey)
file_occupancy <- function(df_occupancy_day,nsurvey=5){
  df1_occ <- fichier_occupancy(df_occupancy_day,nsurvey=nsurvey)$occ
  df1_eff <- fichier_occupancy(df_occupancy_day,nsurvey=nsurvey)$effort
  
  df_long_occ <- df1_occ %>%
    tibble::rownames_to_column("station") %>%
    tidyr::pivot_longer(
      cols = starts_with("survey_"),
      names_to = "survey",
      values_to = "value") %>%
    mutate(
      survey_num = as.integer(gsub("survey_", "", survey)),
      block = ceiling(survey_num / 6)) %>%
    mutate(
      survey_in_block = paste0("survey_", ((survey_num - 1) %% 6) + 1)) %>%
    dplyr::select(station, block, survey_in_block, value) %>%
    tidyr::pivot_wider(
      names_from = survey_in_block,
      values_from = value) %>%
    arrange(station, block) %>%
    mutate(station=paste0(station,"_",block)) %>%
    dplyr::select(-block) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("station")
  
  df_long_eff <- df1_eff %>%
    tibble::rownames_to_column("station") %>%
    tidyr::pivot_longer(
      cols = starts_with("survey_"),
      names_to = "survey",
      values_to = "value") %>%
    mutate(
      survey_num = as.integer(gsub("survey_", "", survey)),
      block = ceiling(survey_num / 5)) %>%
    mutate(
      survey_in_block = paste0("survey_", ((survey_num - 1) %% 6) + 1)) %>%
    dplyr::select(station, block, survey_in_block, value) %>%
    tidyr::pivot_wider(
      names_from = survey_in_block,
      values_from = value) %>%
    arrange(station, block) %>%
    mutate(station=paste0(station,"_",block)) %>%
    dplyr::select(-block) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("station")
  
  return(list(occ=df_long_occ,eff=df_long_eff))
}

# Function that runs the model with spOccupancy
model_psi_site_monthly_spocc <- function(file_occupancy_monthly,
                                         n.samples = 5000,
                                         n.burn = 2000,
                                         n.chains = 3,
                                         n.thin = 5){
  df_occ <- file_occupancy_monthly$occ
  df_occ[df_occ>1] <- 1
  site_cov <- sub("_[^_]*$", "", rownames(df_occ))
  
  library(spOccupancy)
  y <- df_occ
  site_cov_df <- data.frame(site_cov = site_cov)
  site_cov_df$site_cov <- factor(site_cov_df$site_cov)
  form_occ <- ~ site_cov
  form_det <- ~ 1
  
  all_na_sites <- apply(y, 1, function(x) all(is.na(x)))
  y2 <- y[!all_na_sites, , drop = FALSE]
  site_cov_df2 <- site_cov_df[!all_na_sites, , drop = FALSE]
  
  mod_sp <- PGOcc(
    occ.formula = ~ site_cov,
    det.formula = ~ 1,
    data = list(
      y = y2,
      occ.covs = site_cov_df2),
    n.samples = n.samples,
    n.burn = n.burn,
    n.chains = n.chains,
    n.thin = n.thin,
    verbose = TRUE)
  
  return(mod_sp)
}
# Function that evaluate the occupancy probability at each site from the spOccupancy model
evaluate_psi_site <- function(mod_monthly){
  
  beta_samp <- mod_monthly$beta.samples
  beta_summary <- data.frame(
    param = colnames(beta_samp),
    mean = apply(beta_samp, 2, mean),
    lcl  = apply(beta_samp, 2, quantile, probs = 0.025),
    ucl  = apply(beta_samp, 2, quantile, probs = 0.975)) %>%
    mutate(psi_mean=plogis(mean),
           psi_upp=plogis(ucl),
           psi_low=plogis(lcl)) %>%
    rename(station_year=param)
  beta_summary$station_year[beta_summary$station_year == "(Intercept)"] <- "site_covHAB_M1_2023"
  
  res <- beta_summary %>%
    select(station_year,psi_mean,psi_upp,psi_low) %>%
    mutate(station_year = sub("^.*site_cov", "", station_year))
  
  return(res)
}

# Functions that run the Structural equation models
## Modèle M1 : elevation, et wolf/humans => herbivores
model_SEM_M1 <- function(psi_wolf_site_monthly,psi_chamois_site_monthly,
                         psi_reddeer_site_monthly,psi_roedeer_site_monthly,
                         psi_wildboar_site_monthly,
                         camera_information){
  
  big_df_sem <- psi_wolf_site_monthly %>% select(station_year, psi_wolf = psi_mean) %>%
    left_join(psi_chamois_site_monthly %>% select(station_year, psi_chamois = psi_mean), by = "station_year") %>%
    left_join(psi_roedeer_site_monthly %>% select(station_year, psi_roedeer = psi_mean), by = "station_year") %>%
    left_join(psi_reddeer_site_monthly %>% select(station_year, psi_reddeer = psi_mean), by = "station_year") %>%
    left_join(psi_wildboar_site_monthly %>% select(station_year, psi_wildboar = psi_mean), by = "station_year") %>%
    mutate(station=str_remove(station_year, "_\\d{4}$")) %>%
    left_join(camera_information, by = c("station")) %>%
    mutate(elevation_std = scale(elevation), strava_std = scale(strava)) %>%
    mutate(
      psi_wolf = qlogis(psi_wolf),
      psi_chamois = qlogis(psi_chamois),
      psi_roedeer = qlogis(psi_roedeer),
      psi_reddeer = qlogis(psi_reddeer),
      psi_wildboar = qlogis(psi_wildboar))
  
  library(piecewiseSEM)
  SEM_M1 <- psem(
    lm(psi_wolf ~ elevation_std+strava_std,data=big_df_sem),
    lm(psi_chamois ~ elevation_std+strava_std+psi_wolf,data=big_df_sem),
    lm(psi_reddeer ~ elevation_std+strava_std+psi_wolf,data=big_df_sem),
    lm(psi_roedeer ~ elevation_std+strava_std+psi_wolf,data=big_df_sem),
    lm(psi_wildboar ~ elevation_std+strava_std+psi_wolf,data=big_df_sem))
  
  sem_summary <- summary(SEM_M1)
  return(sem_summary)
}

## Modèle M2 : environmental features cov => herbivores and wolf, et wolf/humans => herbivores
model_SEM_M2 <- function(psi_wolf_site_monthly,psi_chamois_site_monthly,
                         psi_reddeer_site_monthly,psi_roedeer_site_monthly,
                         psi_wildboar_site_monthly,
                         camera_information){
  
  big_df_sem <- psi_wolf_site_monthly %>% select(station_year, psi_wolf = psi_mean) %>%
    left_join(psi_chamois_site_monthly %>% select(station_year, psi_chamois = psi_mean), by = "station_year") %>%
    left_join(psi_roedeer_site_monthly %>% select(station_year, psi_roedeer = psi_mean), by = "station_year") %>%
    left_join(psi_reddeer_site_monthly %>% select(station_year, psi_reddeer = psi_mean), by = "station_year") %>%
    left_join(psi_wildboar_site_monthly %>% select(station_year, psi_wildboar = psi_mean), by = "station_year") %>%
    mutate(station=str_remove(station_year, "_\\d{4}$")) %>%
    left_join(camera_information, by = c("station")) %>%
    mutate(elevation_std = scale(elevation), strava_std = scale(strava),
           forest_std = scale(forest),slope_std = scale(slope),
           northness_std = scale(northness), distance_roads_std = scale(distance_roads),
           distance_trails_std = scale(distance_trails),distance_buildings_std = scale(distance_buildings),
           livestock_area_std = scale(livestock_area)) %>%
    mutate(
      psi_wolf = qlogis(psi_wolf),
      psi_chamois = qlogis(psi_chamois),
      psi_roedeer = qlogis(psi_roedeer),
      psi_reddeer = qlogis(psi_reddeer),
      psi_wildboar = qlogis(psi_wildboar))
  
  library(piecewiseSEM)
  SEM_M2 <- psem(
    lm(psi_wolf ~ forest_std+distance_roads_std+distance_buildings_std+livestock_area_std+strava_std,data=big_df_sem),
    lm(psi_chamois ~ forest_std+livestock_area_std+distance_roads_std+distance_trails_std+
         psi_wolf+strava_std,data=big_df_sem),
    lm(psi_reddeer ~ distance_roads_std+distance_trails_std+distance_buildings_std+slope_std+forest_std+northness_std+livestock_area_std+
         strava_std+psi_wolf,data=big_df_sem),
    lm(psi_roedeer ~ strava_std+forest_std+livestock_area_std+psi_wolf,data=big_df_sem),
    lm(psi_wildboar ~ strava_std+forest_std+slope_std+distance_roads_std+distance_trails_std+northness_std+livestock_area_std+psi_wolf,data=big_df_sem))
  
  sem_summary <- summary(SEM_M2)
  return(sem_summary)
}

## Modèle M3_correl : env cov => ongulés et loup, et loup/humains => ongulés et ongulés<=>ongulés par corrélations
model_SEM_M3_correl <- function(psi_wolf_site_monthly,psi_chamois_site_monthly,
                                psi_reddeer_site_monthly,psi_roedeer_site_monthly,
                                psi_wildboar_site_monthly,
                                camera_information){
  
  big_df_sem <- psi_wolf_site_monthly %>% select(station_year, psi_wolf = psi_mean) %>%
    left_join(psi_chamois_site_monthly %>% select(station_year, psi_chamois = psi_mean), by = "station_year") %>%
    left_join(psi_roedeer_site_monthly %>% select(station_year, psi_roedeer = psi_mean), by = "station_year") %>%
    left_join(psi_reddeer_site_monthly %>% select(station_year, psi_reddeer = psi_mean), by = "station_year") %>%
    left_join(psi_wildboar_site_monthly %>% select(station_year, psi_wildboar = psi_mean), by = "station_year") %>%
    mutate(station=str_remove(station_year, "_\\d{4}$")) %>%
    left_join(camera_information, by = c("station")) %>%
    mutate(elevation_std = scale(elevation), strava_std = scale(strava),
           forest_std = scale(forest),slope_std = scale(slope),
           northness_std = scale(northness), distance_roads_std = scale(distance_roads),
           distance_trails_std = scale(distance_trails),distance_buildings_std = scale(distance_buildings),
           livestock_area_std = scale(livestock_area)) %>%
    mutate(
      psi_wolf = qlogis(psi_wolf),
      psi_chamois = qlogis(psi_chamois),
      psi_roedeer = qlogis(psi_roedeer),
      psi_reddeer = qlogis(psi_reddeer),
      psi_wildboar = qlogis(psi_wildboar))
  
  library(piecewiseSEM)
  SEM_M3 <- psem(
    lm(psi_wolf ~ forest_std+distance_roads_std+distance_buildings_std+livestock_area_std+strava_std,data=big_df_sem),
    lm(psi_chamois ~ forest_std+livestock_area_std+distance_roads_std+distance_trails_std+
         psi_wolf+strava_std,data=big_df_sem),
    lm(psi_reddeer ~ distance_roads_std+distance_trails_std+distance_buildings_std+slope_std+forest_std+northness_std+livestock_area_std+
         strava_std+psi_wolf,data=big_df_sem),
    lm(psi_roedeer ~ strava_std+forest_std+livestock_area_std+psi_wolf,data=big_df_sem),
    lm(psi_wildboar ~ strava_std+forest_std+slope_std+northness_std+livestock_area_std+psi_wolf+distance_roads_std+distance_trails_std,data=big_df_sem),
    
    psi_chamois %~~% psi_reddeer,
    psi_chamois %~~% psi_roedeer,
    psi_chamois %~~% psi_wildboar,
    psi_reddeer %~~% psi_roedeer,
    psi_reddeer %~~% psi_wildboar,
    psi_roedeer %~~% psi_wildboar)
  
  sem_summary <- summary(SEM_M3)
  return(sem_summary)
}
  
