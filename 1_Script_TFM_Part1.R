
rm(list = ls()) # Clear all
# SETUP ####

# the following command sets the working directory to the folder where this
# script is located (similar to RStudio menu "Session - Set Working Directory -
# To Source File Location"):
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Check the directory :
getwd()

# LOAD PACKAGES ####

library(openxlsx) 
library(terra)  
library(geodata) 
library(sdmpredictors)  
library(fuzzySim)  
library(pastclim)
library(gridExtra)
library(ncdf4)
library(raster)
library(lattice)
library(ggplot2)
library(gam)  
library(maxnet)  
library(randomForest)  
library(gbm)  
library(embarcadero) 
library(beanplot)
library(ENMTools)
library(tidyverse)
library(magick)
library(plyr)
library(rworldmap)
library(modEvA) 
library(blockCV)  
library(landscapemetrics)
library(ncdf4)
library(fields)
library(maps)
library(ggplot2)
library(sf)
library(car)
library(pROC)
library(geosphere)


# Specify where you want to save the outputs // Pon aquí tu dirección

my_path<- " "



#_________________________________________________________________________
################# OPEN DATA AND SELECT SPECIES
#

#  "Raw.data_filtered.xlsx" should be in the same folder where
#   you have the script file

Presence_data<-read.xlsx("Raw.data_filtered.xlsx", rowNames=FALSE, 
          colNames=TRUE, sheet="Fauna")
head(Presence_data)

# for each species just write exactly the same name you have in Raw.data_filtered.xlsx below, in "species"
species <- "Homo_neand" # specify the species. 
# "Homo_neand" por "Homo_sapiens" en todo el código.
db <- subset(Presence_data,Homo_neand==1) #specify the species 

head(db) 
nrow(db) 

# Create as new data frame only with data of interest

Presence_data<- data.frame(db$`Site/Level`)
Presence_data$x <- db$x
Presence_data$y<- db$y
Presence_data$Age <- db$Age
Presence_data$SD <- db$SD
Presence_data$presence <- 1
Presence_data

# Create a new folder for each species. Here will be saved all outcomes obtained.
pres_folder <- paste0(species)
pres_folder
dir.create(pres_folder)  # '../' goes up one level from the current working directory, so this creates the 'species' folder 


# map the species occurrence data:
# Define the extent for plotting
mywindow <- c(-15, 60, 33, 75)  # xmin, xmax, ymin, ymax
#x11()
# Extract the subset of the world map data
world_map_subset <- map_data("world") %>%
  filter(long >= mywindow[1] & long <= mywindow[2] & lat >= mywindow[3] & lat <= mywindow[4] &
           region != "Morocco" & region != "Algeria" & region != "Tunisia")

# Plot the subset of the world map
ggplot() +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "black") +
  geom_point(data = Presence_data, aes(x = x, y = y), color = "red")+
  theme_classic()




#######___________________________________________________________________________

# CHECK ALL PREDICTOR VARIABLES ####

get_vars_for_dataset(dataset="Krapp2021", details = TRUE)

# if this is the first time you run the code, run the following line chunk:
 set_data_path("C:/Users/vidalma/Desktop/KRAPP_2021") # directory where .nc files will be stored

# El siguiente chunk (donwload_dataset(c("bio01", "bio04", etc.))) solo lo tienes que correr la primera vez que corres el código. Tardará un rato.
download_dataset(c("bio01","bio04", "bio05", "bio06","bio07", "bio08", "bio09", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",
                   "bio17", "bio18", "bio19", "npp", "lai", "altitude", "rugosity"),data="Krapp2021")

# Take into account that these NetCDF files are directly dowloaded with a 0.5 spatial resolution. Each file was downscaled with
# bilinear interpolation to 0.25 resolution with the OriginLab software. For details, see the main text.

# Select all variables
bio_variables = c("bio01","bio04", "bio05", "bio06", "bio07", "bio08", "bio09", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",
                  "bio17", "bio18", "bio19", "npp", "lai", "altitude", "rugosity")


Presence_data$time_bp<- Presence_data$Age * -1 # in pastclim, Age  should be negative, so multiply it by * -1 
head(Presence_data)
Presence_data_p<-Presence_data

# GET VARIABLES IN ALL PRESENCE DATA
# THIS STEP IS ONLY NECESSARY THE FIRST TIME YOU RUN THE CODE

# first, change name of x and y columns because some packages require "longitude" and "latitude" to be specifically stated:

colnames(Presence_data)[2] = "longitude"
colnames(Presence_data)[3] = "latitude"

# get the values of each variable for each location with the presence of the selected species:

Variables<- location_slice( x = Presence_data, bio_variables = c("bio01","bio04", "bio05", "bio06", "bio08", "bio09", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",
                                 "bio17", "bio18", "bio19",  "npp",  "lai", "altitude", "rugosity"), dataset = "Krapp2021", nn_interpol = TRUE, buffer=TRUE)

head(Presence_data)
head(Variables)


# Replace NA values in Variables with Presence_data values
Variables <- Variables %>%
  mutate(
    db..Site.Level.= coalesce(db..Site.Level., Presence_data$db..Site.Level),
    longitude = coalesce(longitude, Presence_data$longitude),
    latitude = coalesce(latitude, Presence_data$latitude),
    Age = coalesce(Age, Presence_data$Age),
    SD = coalesce(SD, Presence_data$SD),
    presence = coalesce(presence, Presence_data$presence),
    time_bp = coalesce(time_bp, Presence_data$time_bp)
    # Agrega más columnas aquí si es necesario
  )

# Aproxima
Variables$time_bp_slice <- round(Variables$time_bp, -3)
# Ver los primeros registros del DataFrame actualizado
head(Variables)


# save outputs in your species' folder

write.xlsx(Variables, paste0( species, "/Variables_Krapp2021_", species, ".xlsx"), colNames = TRUE, rowNames = FALSE, detectDates = FALSE)


# I you already have run the code and you have this file ("Variables"), just load it:

Presence_data <-read.xlsx(paste0( species, "/Variables_Krapp2021_", species, ".xlsx"), rowNames=FALSE, 
                          colNames=TRUE, sheet="Sheet 1")
head(Presence_data)
nrow(Presence_data)

Presence_data<- Presence_data[complete.cases(Presence_data), ]


#_________________________________________________________________
# REPEAT THE SAME PROCEDURE, BUT FILTERING PRESENCE POINTS:
#
#___________________________________________________________
# Sample bias in species occurrence data has long been a recognized issue in SDM. 
# However, environmental filtering of observation data can improve model predictions 
# by reducing redundancy in environmental (e.g. climatic) hyper-space (Varela et al. 2014). 
# The following command removes duplicates and thin the points: only one occurrence per pixel and time-step


 # _________
 #
 # BUFFER
 #
 #__________
 # Before obtaining climate variables from points, we set an area enclosing all the fossil 
 # localities where the species was present, and then we create a buffer
 # around the polygon with a radius equal to 10% of the maximum distance between actual 
 # fossil occurrences:
#x11()

# Write the directory where you have the "gadm36_adm0_r5_pk.rds" fil
countries <- readRDS(".../gadm36_adm0_r5_pk.rds") 


 pres_points <- vect(Presence_data, geom = c("longitude", "latitude"), keepgeom = TRUE, crs = "epsg:4326")  
 buff_dist <- max(distance(pres_points) * 0.1 , na.rm = TRUE)  # can take time if there are many points!
 buff_dist
 pres_buff <- aggregate(buffer(pres_points, width = buff_dist))
 plot(pres_buff)
 plot(pres_points, col = "blue", add = TRUE)
 plot(countries, border = "grey", add = TRUE)
 

 
# Get the spatial resolution of one layer for thinning:

 climate_layer <- region_slice(
   time_bp = -40000,
   bio_variables = bio_variables,
   dataset = "Krapp2021"
 )
 
# plot(climate_layer)

 plot(climate_layer[[1]])
 layers_cut <- crop(climate_layer[[1]], pres_buff, mask = TRUE)
  # plot(layers_cut)
 # Create empty dataset where outcomes will be saved
  Outcomes_v_thinned <- data.frame()

  # The following command thin presence point to avoid overrepresentation of points. There will be only 1 presence/pixel/age
  for (step in unique(Presence_data$time_bp_slice)){
  dat <- subset(Presence_data, time_bp_slice== step) # select only points with same chronology
  if (nrow(dat)> 0) 
  {
  dat <- gridRecords(rst = layers_cut[[1]], pres.coords = dat[ , c("longitude", "latitude")], plot = TRUE) # thin, for details see ?gridRecords
  dat<- subset(dat, presence==1)
  nrow(dat)
  colnames(dat)[2] = "longitude"
  colnames(dat)[3] = "latitude"
  dat$time_bp <- step
  dat<- dat[ , !(names(dat) %in% "bio01")] # the raster with the variable bios01 was used for thinning, now we remove this column from dat
  Variables_thinned<- location_slice( x = dat, bio_variables = c("bio01","bio04", "bio05", "bio06", "bio07", "bio08", "bio09", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",
                                                       "bio17", "bio18", "bio19", "npp", "lai", "altitude", "rugosity"), dataset = "Krapp2021", nn_interpol = TRUE, buffer=TRUE)
 head(Variables_thinned)
 head(dat)
   Variables_thinned <- Variables_thinned %>%
    mutate(
      presence= coalesce(presence, dat$presence),
      longitude = coalesce(longitude, dat$longitude),
      latitude = coalesce(latitude, dat$latitude),
      cells = coalesce(cells, dat$cells),
      time_bp = coalesce(time_bp, dat$time_bp),
      time_bp_slice = coalesce(time_bp_slice, dat$time_bp_slice)
   
    )
  
  
  Outcomes_v_thinned <- rbind(Outcomes_v_thinned, Variables_thinned)
  }
}


#PLOT RAW PRESENCE POINTS AND FILTERED PRESENCES
plot(countries, ext = mywindow)
plot(countries, xlim = mywindow[1:2], ylim = mywindow[3:4])
points(Presence_data[ , c("longitude", "latitude")], col = "blue")
points(Outcomes_v_thinned[ , c("longitude", "latitude")], col = "red")
nrow(Presence_data)
nrow(Outcomes_v_thinned)

Outcomes_v_thinned


# save outputs in your species' folder
write.xlsx(Outcomes_v_thinned, paste0( species, "/Outcomes_v_thinned_Krapp2021_", species, ".xlsx"))


# I you already have run the code and you have this file ("Outcomes_v_thinned"), just load it:

Presence_data_thinned <-read.xlsx(paste0( species, "/Outcomes_v_thinned_Krapp2021_", species, ".xlsx"), rowNames=FALSE, 
                          colNames=TRUE, sheet="Sheet 1")

head(Presence_data_thinned)


# CYCLE OVER 50 * n observations TIME STEPS AND SAMPLE THE ENVIRONMENTAL CONDITIONS ####
# THIS STEP IS ONLY NECESSARY THE FIRST TIME YOU RUN THE CODE

# First, we create an empty list:
climate_list<-list()

for (step in unique(Presence_data_thinned$time_bp)){
  # extract climate for the world

  step_climate <- region_slice(
    time_bp = step,
    bio_variables = bio_variables,
    dataset = "Krapp2021"
  )
  
  # subset calibration area:
  crs(step_climate) <- crs(pres_buff)
  step_climate <- terra::mask(step_climate, pres_buff)# select calibration area (press_buff) 
  # Select number of observations per time step / chronology:
  n <- subset(Presence_data_thinned, time_bp== step) # select only points with same chronology
  nrow(n) 
  # Select pesudoabsences:
  step_values <- terra::spatSample(step_climate, 50 * nrow(n), na.rm=TRUE, cells=TRUE, xy=TRUE) 
  step_values$time_bp<-step # Time step
  climate_list[[as.character(step)]]<-step_values
}

# combine all into a single matrix and include "presence", "x" (longitude) and "y" (latitude) variables
baseline_climate <- do.call(rbind, climate_list)
baseline_climate$presence<- 0
baseline_climate$longitude<- baseline_climate$x
baseline_climate$latitude<- baseline_climate$y
nrow(baseline_climate)
# save output
write.xlsx(baseline_climate, paste0( species, "/Baseline_allvars", species, ".xlsx"))

# merge datasets

head(baseline_climate)
head(Presence_data_thinned)

db<- rbind.fill(Presence_data_thinned, baseline_climate)

# remove locations without climate (because of the regridding a few observations may have been moved to sea cells or could be covered in ice)
#db <- db_climate[complete.cases(db_climate), ]
#save output

colnames(db)

# Now, lets eliminate the empty or unnecesary columns for the next steps:
db <- db[, !(colnames(db) %in% c("cells", "time_bp_slice", "cell", "x", "y"))]
colnames(db)



write.xlsx(db, paste0( species, "/ALL", species, ".xlsx"))


# I you already have run the code and you have this file ("ALL+species"), just load it:

db<- read.xlsx(paste0( species, "/ALL", species, ".xlsx"), rowNames = FALSE, colNames = TRUE, sheet="Sheet 1")
head(db)
nrow(db)

head(db)
colnames(db)

# PLOT PSEUDOABSENCE AND PRESENCE DATA: ###
 plot(countries, ext = mywindow)
 # Plot the 'countries' data within the specified extent
 plot(countries, xlim = mywindow[1:2], ylim = mywindow[3:4])
 points(subset(db, db$presence == 1, select = c("longitude", "latitude")), pch = 20, cex = 0.7, col="red")
 points(subset(db, db$presence == 0, select = c("longitude", "latitude")), pch = 20, cex = 0.3, col="blue")

#head(db)


# SPECIES DISTRIBUTION MODELS ####
 
 # First, load all variables for each 1 k years from  between 55 and 27 k yrs BP 
 # (These files are obtained from the file "Get_Env_F.R")

 
 path_env_f <- " " # Directory where you have the files
 
 Climate_55K_BP <- readRDS(paste0( path_env_f, "/Climate_55.rds"))
 Climate_54K_BP <- readRDS(paste0( path_env_f, "/Climate_54.rds"))
 Climate_53K_BP <- readRDS(paste0( path_env_f, "/Climate_53.rds"))
 Climate_52K_BP <- readRDS(paste0( path_env_f, "/Climate_52.rds"))
 Climate_51K_BP <- readRDS(paste0( path_env_f, "/Climate_51.rds"))
 Climate_50K_BP <- readRDS(paste0( path_env_f, "/Climate_50.rds"))
 Climate_49K_BP <- readRDS(paste0( path_env_f, "/Climate_49.rds"))
 Climate_48K_BP <- readRDS(paste0( path_env_f, "/Climate_48.rds"))
 Climate_47K_BP <- readRDS(paste0( path_env_f, "/Climate_47.rds"))
 Climate_46K_BP <- readRDS(paste0( path_env_f, "/Climate_46.rds"))
 Climate_45K_BP <- readRDS(paste0( path_env_f, "/Climate_45.rds"))
 Climate_44K_BP <- readRDS(paste0( path_env_f, "/Climate_44.rds"))
 Climate_43K_BP <- readRDS(paste0( path_env_f, "/Climate_43.rds"))
 Climate_42K_BP <- readRDS(paste0( path_env_f, "/Climate_42.rds"))
 Climate_41K_BP <- readRDS(paste0( path_env_f, "/Climate_41.rds"))
 Climate_40K_BP <- readRDS(paste0( path_env_f, "/Climate_40.rds"))
 Climate_39K_BP <- readRDS(paste0( path_env_f, "/Climate_39.rds"))
 Climate_38K_BP <- readRDS(paste0( path_env_f, "/Climate_38.rds"))
 Climate_37K_BP <- readRDS(paste0( path_env_f, "/Climate_37.rds"))
 Climate_36K_BP <- readRDS(paste0( path_env_f, "/Climate_36.rds"))
 Climate_35K_BP <- readRDS(paste0( path_env_f, "/Climate_35.rds"))
 Climate_34K_BP <- readRDS(paste0( path_env_f, "/Climate_34.rds"))
 Climate_33K_BP <- readRDS(paste0( path_env_f, "/Climate_33.rds"))
 Climate_32K_BP <- readRDS(paste0( path_env_f, "/Climate_32.rds"))
 Climate_31K_BP <- readRDS(paste0( path_env_f, "/Climate_31.rds"))
 Climate_30K_BP <- readRDS(paste0( path_env_f, "/Climate_30.rds"))
 Climate_29K_BP <- readRDS(paste0( path_env_f, "/Climate_29.rds"))
 Climate_28K_BP <- readRDS(paste0( path_env_f, "/Climate_28.rds"))
 Climate_27K_BP <- readRDS(paste0( path_env_f, "/Climate_27.rds"))
 

 # COMPUTE SPECIES DISTRIBUTION MODELS ###
 
 ######################### Select the predictive variables
 
 # define the modelling columns:

 names(db)
 dat <- db
 spc_col <- "presence"  # species presence/absence column is named "presence" in the dataset
 var_cols <-  names(db[,5:25]) 
 var_cols
 

  
  ##
  # RMSE (Root Mean Sqaure Error)
  
  calculate_rmse <- function(model) {
    predicted_values <- predict(model)
    observed_values <- dat$presence
    rmse <- sqrt(mean((observed_values - predicted_values)^2))
    return(rmse)
  }
  
  # VIF (Variance Inflation Factor)
  
  calculate_vif <- function(model) {
    vif_values <- tryCatch({
      vif(model)
    }, error = function(e) {
      return(6)  # 6 in case of error
    })
    return(vif_values)
  }
  
  
  bio_variables = c("bio01","bio04", "bio05", "bio06", "bio07", "bio08", "bio09", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",
                    "bio17", "bio18", "bio19", "npp", "lai", "altitude", "rugosity")

  all_predictors <- bio_variables  # 
  # Check highly correlated variables (r^2 > 0.8)
  head(dat)
  head(dat)[,5:25]
  correlation_matrix <- cor(dat[,5:25])

  # Function to select variables with correlation coeffficients < 0.8
  select_variables <- function(matrix) {
    selected_variables <- character(0)
    remaining_variables <- colnames(matrix)
    
    while (length(remaining_variables) > 0) {
      reference_variable <- remaining_variables[1]
      selected_variables <- c(selected_variables, reference_variable)
      
      # Filter variables with r < 0.8 
      filtered_variables <- remaining_variables[-1][abs(matrix[reference_variable, remaining_variables[-1]]) < 0.8]
      
      # Filtered variables
      remaining_variables <- filtered_variables
    }
    
    return(selected_variables)
  }
  
  # Get list of filtered variables
  selected_vars <- select_variables(correlation_matrix)
  
  
  # Print variables
  cat("Selected variables:", selected_vars, "\n")
  
  
  
  
  
  
  
  ## GENERALIZED LINEAR MODEL (GLM) ####

 form_glm <- as.formula(paste(spc_col, "~", paste(selected_vars, collapse = "+")))
 form_glm
 mod_glm <- glm(formula = form_glm, family = binomial, data = dat)
 summary(mod_glm)
 
 # Initialize the variables for the best model
 best_AUC_value <- -Inf
 best_rmse_value <- Inf
 best_model <- NULL
 best_predictors <- NULL
 counter<-0
 
 dev.off()
 
 for (i in 2:length(selected_vars)) {
   combinations <- combn(selected_vars, i)
   
   for (j in 1:ncol(combinations)) {
     selected_predictors <- combinations[, j]
    form_glm <- as.formula(paste(spc_col, "~", paste(selected_predictors, collapse = "+")))
     mod_glm <- glm(formula = form_glm, family = binomial, data = dat)
     summary(mod_glm)
     head(dat)

     
     p_values <-summary(mod_glm)$coefficients[, "Pr(>|z|)"]
     AUC_value<- AUC(mod_glm)$AUC
     vif_values <- calculate_vif(mod_glm)
     rmse_value<-  calculate_rmse(mod_glm)
     
         # In any VIFs > 5 = Vif_exclusion = 1
         if (any(vif_values > 5) ) {
           Vif_exclusion <- 1
         } else {
           Vif_exclusion <- 0
         }
         
         # Check whether this is the best model so far
         if (Vif_exclusion == 0 && AUC_value > best_AUC_value || (AUC_value == best_AUC_value && rmse_value < best_rmse_value)) {
  
             best_AUC_value <- AUC_value
             print(best_AUC_value)
             best_rmse_value <- rmse_value
             print(best_rmse_value)
             best_model_glm <- mod_glm
             best_predictors_glm <- selected_predictors
             print(best_predictors_glm)
           }
         }
       }
     
 
 
 # Best model found after the loop
 cat("Final best AUC: ", best_AUC_value, "\n")
 cat("Final best RMSE: ", best_rmse_value, "\n")
 cat("Final best predictors: ", best_predictors_glm, "\n")
 print(best_predictors_glm)

 # Variables used for each species
 # Write here the selected variables. For instance: selected_vars_glm <- c( "bio01" ,  "bio15" ,    "rugosity")
 selected_vars_glm <- c( ) # 
 
 
  # Run the definitive model
 form_glm <- as.formula(paste(spc_col, "~", paste(selected_vars_glm, collapse = "+")))
 form_glm
 mod_glm <- glm(formula = form_glm, family = binomial, data = dat)
 summary(mod_glm)
 AUC(mod_glm)$AUC
 
 
 ## GENERALIZED ADDITIVE MODEL (GAM) ####

 form_gam <- as.formula(paste(spc_col, "~", paste0("s(",  selected_vars, ")", collapse = "+")))  # GAM with smoothing splines ('s')
 form_gam
 
 mod_gam <- gam(formula = form_gam, family = binomial, data = dat)
 summary(mod_gam)

 # Initialize the variables for the best model
 best_AUC_value <- -Inf
 best_rmse_value <- Inf
 best_model <- NULL
 best_predictors <- NULL
 counter<-0
 dev.off()
 
 for (i in 2:length(selected_vars)) {
   combinations <- combn(selected_vars, i)
   
   for (j in 1:ncol(combinations)) {
     selected_predictors <- combinations[, j]
     form_gam <- as.formula(paste(spc_col, "~", paste0("s(",  selected_predictors, ")", collapse = "+")))  
     mod_gam <- gam(formula = form_gam, family = binomial, data = dat)
     summary(mod_gam)
     head(dat)
     
     p_values <-summary(mod_gam)$coefficients[, "Pr(>F)"]
     AUC_value<- AUC(mod_gam)$AUC
     vif_values <- calculate_vif(mod_gam)
     rmse_value<-  calculate_rmse(mod_gam)
     
     # In any VIFs > 5 = Vif_exclusion = 1
     if (any(vif_values > 5) ) {
       Vif_exclusion <- 1
     } else {
       Vif_exclusion <- 0
     }
     
     # Check whether this is the best model so far
     if (Vif_exclusion == 0 && AUC_value > best_AUC_value || (AUC_value == best_AUC_value && rmse_value < best_rmse_value)) {
       # Also check that all variables have p-values < 0.05
         best_AUC_value <- AUC_value
         print(best_AUC_value)
         best_rmse_value <- rmse_value
         print(best_rmse_value)
         best_model_gam <- mod_gam
         best_predictors_gam <- selected_predictors
         print(best_predictors_gam)
       }
     }
   }
 
 

 cat("Final best AUC: ", best_AUC_value, "\n")
 cat("Final best RMSE: ", best_rmse_value, "\n")
 cat("Final best predictors: ", best_predictors_gam, "\n")
 print(best_predictors_gam)
 

 # Variables used for each species. The same as before, write the selected variables:
 selected_vars_gam <- c() 

 
 # Run the definitive model
 form_gam <- as.formula(paste(spc_col, "~", paste0("s(",  selected_vars_gam, ")", collapse = "+")))  # GAM with smoothing splines ('s')
 form_gam
 
 mod_gam <- gam(formula = form_gam, family = binomial, data = dat)
 summary(mod_gam)
 
 
 
 ## MAXIMUM ENTROPY (MAXENT) ####
 

 mod_mxt <- maxnet(p = dat[ , spc_col], data = dat[ , selected_vars], f = maxnet.formula(dat[ , spc_col], dat[ , selected_vars])) 
 summary (mod_mxt)

# VIF cannot be computed from MAXENT in the loop, SO WE FIRST TEST THE VIF AND SELECT THE VARIABLES:
 
 corSelect(data = dat, var.cols =   selected_vars, cor.thres = 0.8, select = "VIF")
 option1<- c( "bio04",    "bio08",    "bio12",    "bio15",    "npp" ,     "lai" ,     "altitude", "rugosity") # Cambia aquí las variables de option1 por las que te salgan al correr la línea anterior
 corSelect(data = dat, var.cols =   option1, cor.thres = 0.8, select = "VIF")
 selected_vars<- option1
   
   
 selected_vars



# Initialize the variables for the best model
best_AUC_value <- -Inf
best_rmse_value <- Inf
best_model <- NULL
best_predictors <- NULL
counter<-0

for (i in 2:length(selected_vars)) {
  combinations <- combn(selected_vars, i)
  cat("Combinación ", counter, ": ", combinations, "\n")
  
  for (j in 1:ncol(combinations)) {
    selected_predictors <- combinations[, j]
    mod_mxt <- maxnet(p = dat[ , spc_col], data = dat[ , selected_predictors], f = maxnet.formula(dat[ , spc_col], dat[ , selected_predictors])) 
    summary(mod_mxt)
    # Get predictions from the model
    predicted_scores <- predict(mod_mxt, newdata = dat[ , selected_predictors])
    
    # Compute AUC
    roc_obj <- roc(dat[ , spc_col], predicted_scores)
    AUC_value <- auc(roc_obj)
    plot(AUC_value, main=paste0( species, " :", AUC_value), col="red")
    rmse_value<-  sqrt(mean((predicted_scores - dat[, spc_col])^2))
    

    # Check whether this is the best model so far
    if (AUC_value > best_AUC_value || (AUC_value == best_AUC_value && rmse_value < best_rmse_value)) {
      # Also check that all variables have p-values < 0.05
        best_AUC_value <- AUC_value
        print(best_AUC_value)
        best_rmse_value <- rmse_value
        print(best_rmse_value)
        best_model_mxt<- mod_mxt
        best_predictors_mxt <- selected_predictors
        print(best_predictors_mxt)
      }
    }
  }


cat("Final best AUC: ", best_AUC_value, "\n")
cat("Final best RMSE: ", best_rmse_value, "\n")
cat("Final best predictors: ", best_predictors_mxt, "\n")
print(best_predictors_mxt)
selected_vars

corSelect(data = dat, var.cols =  best_predictors_mxt, cor.thres = 0.8, select = "VIF")


# Check VIF: corSelect(data = dat, var.cols =  bio_variables, cor.thres = 0.8, select = "VIF")

# Select the definitive model. 
selected_vars_maxt <- c( ) 


mod_mxt <- maxnet(p = dat[ , spc_col], data = dat[ , selected_vars_maxt], f = maxnet.formula(dat[ , spc_col], dat[ , selected_vars_maxt])) 
summary (mod_mxt)

# Get model predictions
predicted_scores <- predict(mod_mxt, newdata = dat[ , selected_vars])

# Compute AUC
roc_obj <- roc(dat[ , spc_col], predicted_scores)
auc <- auc(roc_obj)
auc


# I you changed selected_vars because of the VIF, now select again the original selected_vars:

selected_vars <- select_variables(correlation_matrix)

# Print variables
cat("Selected variables:", selected_vars, "\n")


 ## BAYESIAN ADDITIVE REGRESSION TREES (BART) ####
 

 # Now, we select minimal subset of relevant variables with BART. BART runs the model with different variable combinations
 # and select those with the lowest RMSE.
  set.seed(123)
  n_iter <- 50 # if it takes too long, reduce the number
  varselect_bart <- variable.step(x.data = dat[ , selected_vars], y.data = dat[ , spc_col], n.trees = 10, iter= n_iter)
  varselect_bart

  #varselect_bart<- selected_vars
 #Write selected variables 
  varselect_bart <- c( ) # H neand
  varselect_bart <- c( ) # H sapiens
 
  
  
  # We check whether the RMSE is reduced meaningully. If so, we use "varlselect_bart" in the mod_bart_final (below),
 # If the RMSE remains identical or similar, we keep using "selected_vars"
#  varselect_bart

 mod_bart <- bart(x.train = dat[ , varselect_bart], y.train = dat[ , spc_col], keeptrees = TRUE)
 summary(mod_bart)
 
 # if you want to use this BART model in future R sessions, you need to explicitly
 # ask for the full information to be included when you next save the model object (see "Saving" section in ?bart):
 invisible(mod_bart$fit$state)
 
 summary(mod_bart)
 
 
 # Get AUC 
 match_auc <- regmatches(capture.output(summary(mod_bart)), regexpr("AUC = [0-9.]+", capture.output(summary(mod_bart))))

 auc_value <- as.numeric(gsub("AUC = ", "", match_auc)) # Guarda también el resultado en un Excel.
 
 
 # Initialize the variables for the best model
 best_AUC_value <- -Inf
 best_model <- NULL
 best_predictors <- NULL
 counter<-0
 
 for (i in 2:length(varselect_bart)) {
   combinations <- combn(varselect_bart, i)
   cat("Combinación ", counter, ": ", combinations, "\n")
   
   for (j in 1:ncol(combinations)) {
     selected_predictors <- combinations[, j]
     mod_bart <- bart(x.train = dat[ , selected_predictors], y.train = dat[ , spc_col], keeptrees = TRUE)
     summary(mod_bart)
     

     # Compute AUC
     

     match_auc <- regmatches(capture.output(summary(mod_bart)), regexpr("AUC = [0-9.]+", capture.output(summary(mod_bart))))
     AUC_value <- as.numeric(gsub("AUC = ", "", match_auc))
     
     # Check whether this is the best model so far
     if (AUC_value > best_AUC_value ) {
    
       best_AUC_value <- AUC_value
       print(best_AUC_value)

       best_model_bart<- mod_bart
       best_predictors_bart <- selected_predictors
       print(best_predictors_bart)
     }
   }
 }
 
 summary(best_model_bart)

  # Check VIF:  corSelect(data = dat, var.cols =  bio_variables, cor.thres = 0.8, select = "VIF")
 
 
 # Calculate RMSE
 true_outcomes <- dat[, spc_col]
 predictions <- colMeans(stats::predict(mod_bart, dat))

 rmse <- sqrt(mean((predictions - true_outcomes)^2))
 rmse
 
 
 
 # SAVE ALL MODELS
 
 models <- list(mod_glm ,mod_gam, mod_mxt, mod_bart)
 names(models) <- c("glm", "gam", "mxt", "bart")
 pred_folder <- paste0( species, "/predictions")
 dir.create(pred_folder, recursive = TRUE)
 saveRDS(models, file = paste0(pred_folder, "/models.rds"))
 
 
 # IF YOU HAVE PREVIOUSLY RUN THE CODE AND SAVED THE GAM, MXT AND BART MODELS, JUST LOAD IT:
 
 # Load models:
setwd(" ") # Write directory
 models <- readRDS(paste0(species, "/predictions/models.rds"))
summary(models)


summary(models$glm)
summary(models$gam)
summary(models$mxt)
summary(models$bart)


 # VALIDATION ####

  
  # Cross block validation
  
  # DIVIDE STUDY AREA INTO SPATIAL BLOCKS ###
  
  # first, convert 'dat' to a spatial points object with its CRS (mind that you
  # need to ascertain the correct CRS to which your coordinates refer! it is 4326
  # in the GBIF data):
  dat <- db
  nrow(dat)
  names(dat)
  head(dat)
  dat_points <- vect(dat, geom = c("longitude", "latitude"), crs = "epsg:4326")
  
  # the size (or range) of the spatial blocks is sometimes defined as the range of
  # spatial autocorrelation in the variables (which you can compute with function
  # 'spatialAutoRange'), but this is now discouraged (see Wadoux et al. 2021,
  # https://doi.org/10.1016/j.ecolmodel.2021.109692) here we'll use a size we
  # consider reasonable for our species and study area, e.g. 500-km blocks:
  layers_cut <- crop(climate_layer[[1]], mywindow, mask = TRUE)
  plot(layers_cut)
  blocks <- spatialBlock(speciesData = as(dat_points, "Spatial"), 
                         rasterLayer = raster(layers_cut[[1]]), theRange =  500000, # 500000 / 1000000
                         k = 5, seed = 2023)  # you can also use an optional additional argument, species = spc_col; see ?spatialBlock
  
  # add the spatial block ID to the dataset:
  dat$foldID <- blocks$foldID


  #dat <- head(dat, - 1)   
  
  # map blocks and presences:
  dev.off()
  plot(vect(blocks$blocks), "folds", col = hcl.colors(5, "geyser"), main = species)
  plot(countries, lwd = 1, add = TRUE)
  plot(subset(dat_points, dat_points$presence == 1), col = "blue", cex = 0.7, add = TRUE)
#x11()
  #save 700-460
  
  # GET PREDICTIONS FOR CROSS-VALIDATION ###
  
  # (i.e., build models again, leaving out each fold in turn)
  
  folds <- sort(unique(dat$foldID))
  
  names(dat)
  spc_col <- "presence"  
  
  # Select the predictive variables for each model

  bart_selected_vars<- c("bio01", "bio04",    "bio12" ,   "bio15",  "lai" ,     "altitude", "rugosity") # H neand. Check it
  mxt_selected_vars<- c( "bio04",    "bio08" ,   "bio12" ,   "bio15",    "npp"  ,    "altitude" ,"rugosity") # H neand. Check it
  
  
  for (f in folds) {
    cat("modelling outside fold", f, "...\n")  # inform of progress
    dat_train <- subset(dat, foldID != f)
    
    mod_glm_fold <- glm(formula = models$glm$formula, family = binomial, data = dat_train)
    dat[ , paste0("glm_fold", f, "_p")] <- predict(mod_glm_fold, newdata = dat, type = "response")
    
    mod_gam_fold <- gam(formula = models$gam$formula, family = binomial, data = dat_train)
    dat[ , paste0("gam_fold", f, "_p")] <- predict(mod_gam_fold, dat, type = "response")
 

    mod_mxt_fold <- maxnet(p = dat_train[ , spc_col], data = dat_train[ , mxt_selected_vars], f = maxnet.formula(dat_train[ , spc_col], dat_train[ , mxt_selected_vars]))  # you can add to 'maxnet.formula' e.g. classes="lq", to use only linear ('l') and quadratic ('q') features
    dat[ , paste0("mxt_fold", f, "_p")] <- predict(mod_mxt_fold, dat, clamp = FALSE, type = "cloglog")
    

    mod_bart_fold <- bart(y.train = dat_train[ , spc_col], x.train = dat_train[ , bart_selected_vars], keeptrees = TRUE, verbose = FALSE)
    dat[ , paste0("bart_fold", f, "_p")] <- colMeans(stats::predict(mod_bart_fold, dat))
    
    
    gc()  # cleanup memory before next loop iteration
  }  # end for f
  
  # see the new predictions added to the data frame:
  head(dat)
  
  
  # EVALUATE EACH MODEL ON ITS VALIDATION FOLD ###
  
  fold_cols <- grep("_fold", names(dat))
  names(dat)[fold_cols]
  
  # choose e.g. 3 metrics with complementary information (discrimination, classification, calibration):
  metrics <- c("AUC", "TSS", "MCS")
  
  # create an empty table to receive the cross-validation results:
  crossval <- as.data.frame(matrix(nrow = length(folds), ncol = length(metrics) * length(names(models))))
  colnames(crossval) <- c(outer(1, metrics, FUN = paste, sep = "_"))
  crossval  # for now it's only filled with NAs
  
  par(mfrow = c(5, 2), mar = c(2, 2, 1.1, 1))
  for (m in names(models))  for (f in folds) {
    fold_name <- paste0("fold", f)
    fold_col <- names(dat)[grep(paste0(m, "_fold", f), names(dat))]
    fold_dat <- subset(dat, foldID == f)
    crossval[f, paste(m, "AUC", sep = "_")] <- AUC(obs = fold_dat[ , spc_col], pred = fold_dat[ , fold_col], simplif = TRUE, plot = TRUE, main = paste(m, "AUC"))
    crossval[f, paste(m, "TSS", sep = "_")] <- threshMeasures(obs = fold_dat[ , spc_col], pred = fold_dat[ , fold_col], thresh = "preval", measures = "TSS", simplif = TRUE, standardize = FALSE, main = paste(m, "TSS"))
    crossval[f, paste(m, "MCS", sep = "_")] <- MillerCalib(obs = fold_dat[ , spc_col], pred = fold_dat[ , fold_col], main = paste(m, "Miller line"))$slope
  }
  
  #Export 550-600
  
  # get more information with boxplots of the cross-validation metrics:

  # remove empty columns
  empty_columns <- sapply(crossval, function(x) all(is.na(x) | x == ""))
  dev.off()
  e<- crossval[, !empty_columns]
  e
  boxplot(e, col = 1:length(metrics), each = length(names(models)), las = 2, main = species)
  abline(h = 1, col = "darkred", lty = 1, lwd = 2)  # remember Miller calibration slope (MCS) should ideally be close to 1 (not bigger = better)
  
  species
  # 500-470
  
  # save outputs
  write.csv(crossval, paste0("F:/SPD_NP_v2/Species/", species, "/performance.csv"))





