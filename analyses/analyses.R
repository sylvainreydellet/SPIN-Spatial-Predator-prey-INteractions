
### Verifying data file exists #####
# Download the files at dryad repository: https://doi.org/10.5061/dryad.h44j0zq28, and store them in the 'data' folder.

if(!file.exists("data/list_information.RData")){
  print("Download the files at dryad repository: https://doi.org/10.5061/dryad.h44j0zq28, and verify that the files are stored in the 'data' folder")
}

### Loading packages #####

packages.in <- c("dplyr","lubridate","stringr","spOccupancy","piecewiseSEM")

for (pkg in packages.in) {if(!(pkg %in% rownames(installed.packages()))) install.packages(pkg)}
for (pkg in packages.in) library(pkg, character.only = TRUE)


############## Data processing ############## 
source('analyses/functions.R')
load("data/list_information.RData")

problemes_23_24 <- list_information$problem_information
problemes_23_24[is.na(problemes_23_24)] <- 0
camera_history_23_24 <- formatting_history(filter_date(problemes_23_24))

data_bauges_cleaned <- list_information$observation
data_bauges_filtered_period <- filter_date(data_bauges_cleaned)

############## Calculation of occupancy probabilities ############## 

### 1. Creating occupancy files 

# Create species-specific dataset
reddeer_df <- data_bauges_filtered_period %>% filter(prediction=="red deer")
chamois_df <- data_bauges_filtered_period %>% filter(prediction=="chamois")
roedeer_df <- data_bauges_filtered_period %>% filter(prediction=="roe deer")
wildboar_df <- data_bauges_filtered_period %>% filter(prediction=="wild boar")
wolf_df <- data_bauges_filtered_period %>% filter(prediction=="wolf")

# Filter independent events (tresholds based on Vanderlocht et al. 2026 results)
reddeer_df_filtered <- filter_temporal_independence(reddeer_df,treshold = 21)
chamois_df_filtered <- filter_temporal_independence(chamois_df,treshold = 16)
roedeer_df_filtered <- filter_temporal_independence(roedeer_df,treshold = 20)
wildboar_df_filtered <- filter_temporal_independence(wildboar_df,treshold = 20)
wolf_df_filtered <- filter_temporal_independence(wolf_df, treshold = 20)

# Create occupancy file : station_year in rows, day in columns, taking CT issues into account
reddeer_occupancy_day <- merging_detection_problemes(reddeer_df_filtered,camera_history_23_24)
chamois_occupancy_day <- merging_detection_problemes(chamois_df_filtered,camera_history_23_24)
roedeer_occupancy_day <- merging_detection_problemes(roedeer_df_filtered,camera_history_23_24)
wildboar_occupancy_day <- merging_detection_problemes(wildboar_df_filtered,camera_history_23_24)
wolf_occupancy_day <- merging_detection_problemes(wolf_df_filtered,camera_history_23_24)

# Create occupancy files : 5-day survey lengths. On line per station x year x month
reddeer_occupancy_monthly_5day <- file_occupancy(reddeer_occupancy_day,nsurvey = 5)
chamois_occupancy_monthly_5day <- file_occupancy(chamois_occupancy_day,nsurvey = 5)
roedeer_occupancy_monthly_5day <- file_occupancy(roedeer_occupancy_day,nsurvey = 5)
wildboar_occupancy_monthly_5day <- file_occupancy(wildboar_occupancy_day,nsurvey = 5)
wolf_occupancy_monthly_5day <- file_occupancy(wolf_occupancy_day,nsurvey = 5)

# Run spOccupancy models with no covariates
modele_psi_site_monthly_wolf <- model_psi_site_monthly_spocc(wolf_occupancy_monthly_5day,n.samples = 10000)
modele_psi_site_monthly_chamois <- model_psi_site_monthly_spocc(chamois_occupancy_monthly_5day,n.samples = 10000)
modele_psi_site_monthly_reddeer <- model_psi_site_monthly_spocc(reddeer_occupancy_monthly_5day,n.samples = 10000)
modele_psi_site_monthly_roedeer <- model_psi_site_monthly_spocc(roedeer_occupancy_monthly_5day,n.samples = 10000)
modele_psi_site_monthly_wildboar <- model_psi_site_monthly_spocc(wildboar_occupancy_monthly_5day,n.samples = 10000)

# Extract the values of occupancy probabilities
psi_wolf_site_monthly <- evaluate_psi_site(modele_psi_site_monthly_wolf)
psi_chamois_site_monthly <- evaluate_psi_site(modele_psi_site_monthly_chamois)
psi_reddeer_site_monthly <- evaluate_psi_site(modele_psi_site_monthly_reddeer)
psi_roedeer_site_monthly <- evaluate_psi_site(modele_psi_site_monthly_roedeer)
psi_wildboar_site_monthly <- evaluate_psi_site(modele_psi_site_monthly_wildboar)

### Run the structural equation models
# M1: only elevation as environmental covariate
model_M1 <- model_SEM_M1(psi_wolf_site_monthly,psi_chamois_site_monthly,
                         psi_reddeer_site_monthly,psi_roedeer_site_monthly,
                         psi_wildboar_site_monthly,camera_information)

# M2: whole set of physical habitat features and passive non-lethal human disturbances
model_M2 <- model_SEM_M2(psi_wolf_site_monthly,psi_chamois_site_monthly,
                         psi_reddeer_site_monthly,psi_roedeer_site_monthly,
                         psi_wildboar_site_monthly,camera_information)

# M3: whole set of physical habitat features and passive non-lethal human disturbances and intra-guild herbivores interactions
model_M3 <- model_SEM_M3_correl(psi_wolf_site_monthly,psi_chamois_site_monthly,
                               psi_reddeer_site_monthly,psi_roedeer_site_monthly,
                               psi_wildboar_site_monthly,camera_information)
