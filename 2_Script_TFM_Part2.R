


rm(list = ls()) 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

my_path <- " " # Your directory

species<- "Homo.sapiens" # Specify Homo.sapiens or Homo_neand



 # MAKE PREDICTIONS AND GET THE RESULTS FOR EACH TIME STEP OF THE MIS3 ####
  
  # the next command creates a folder for output files:
  # setwd("F:/SPD_NP_v2_final/Species")
  dir.create(paste0(species, "/Layers") )  # '../' goes up one level from the current working directory, so this creates the 'outputs' folder just outside the 'scripts' folder
  # Now two additional folders are created, one for the .TIF files and other one for the .JPEG files
  dir.create(paste0(species, "/Layers/TIF") )
  dir.create(paste0(species, "/Layers/JPEG") ) 
  
# Modyfy my_patn to include your repository:

  my_path<-"F:/SPD_NP_v2_final/Species"

path_TIF <- paste0(my_path, "/", species, "/Layers/TIF")
path_JPEG <- paste0(my_path, "/", species, "/Layers/JPEG") 
 
 
 # Now, for visual representation purposes, we will include the coast line
  coast_file <- read.csv("C:/Users/vidalma/Desktop/coastline__120m.csv")
 
 # Convert coordinates from WKT to geometric objects
  library(sf)
  library(terra)
 coast_geom <- st_as_sfc(coast_file$the_geom, crs = 4326)
 
 # Create object sf
 coast_sf <- st_sf(coast_geom, data = coast_file)
 path_env_f <- "C:/Users/vidalma/Desktop/SPD_NP_v2/Clim_basemaps/BILINEAR_INTERPOLATION"
 f <- list.files(path=path_env_f, pattern='rds$', full.names=TRUE) # create a list with all raster files

 library(RColorBrewer)
 mi_paleta <- colorRampPalette(brewer.pal(11, "RdYlBu"))
 clrs <- rev(mi_paleta(100))

 
 dev.off()
 
 
 # LOAD:
 
 species<- "Homo.sapiens"
 setwd("F:/SPD_NP_v2_final/Species/") # replace with your directory
 models <- readRDS(paste0(species, "/predictions/models.rds"))
 summary(models)
 
 summary(models$glm)
 summary(models$gam)
 summary(models$mxt)
 summary(models$bart)
 
 
 # before running the next chunk, check you have the performance_alldata.xlx file with the AUC values
 #of each model to compute the weighted mean value in the ensemble projection.
 
 for (i in 1:length(f)) {
   raster <- readRDS(f[[i]])
  
   glm_p <- predict(stack(raster), models$glm, type = "response")
   glm_fav <- Fav(pred = glm_p, sample.preval = prevalence(model = models$glm))
   plot(glm_fav, col = clrs, range = c(0, 1), main = "GLM")
   glm_fav <- as(glm_fav, "SpatRaster")
   
   gam_p <- predict(stack(raster), models$gam, type = "response")
   gam_fav <- Fav(pred = gam_p, sample.preval = prevalence(model = models$gam))
   plot(gam_fav, col = clrs, range = c(0, 1), main = "GAM")
   gam_fav <- as(gam_fav, "SpatRaster")
   
   mxt_p <- rast(predict(stack(raster), models$mxt, type = "cloglog"))  # raster format conversions needed because 'predict' for MX is not working correctly with 'terra' raster format
   plot(mxt_p, col = clrs, range = c(0, 1), main = "MXT")
   plot(coastsCoarse, add=T, col="black")
   mxt_fav <- as(mxt_p, "SpatRaster")
   
   bart_p <- rast(predict(models$bart, x.layers = stack(raster)))
   bart_fav <- Fav(pred = bart_p, sample.preval = prevalence(model = models$bart))
   plot(bart_fav,  range = c(0, 1),col=clrs, main = paste0("", (26 + i), "bart k yrs BP"))
   bart_fav <- as(bart_fav, "SpatRaster")
 
   # Ensemble during the loop
   # read csv with cross-validation outcomes:
   
   ens_df<- read.csv(paste0(my_path, "/", species, "/performance_alldata.csv"), sep = ";")

   AUC_GLM <- mean(ens_df$glm_AUC)
   AUC_GAM <- mean(ens_df$gam_AUC)
   AUC_MXT <- mean(ens_df$mxt_AUC)
   AUC_BART <- mean(ens_df$bart_AUC)
   
   AUC_GLM
   AUC_GAM 
   AUC_MXT 
   AUC_BART
   
   # Ensemble "response"
   glm_p <- as(glm_p, "SpatRaster")
   gam_p <- as(gam_p, "SpatRaster")
   
   b <- brick(x = c(glm_p, gam_p, mxt_p, bart_p))

   ensemble_r <- weighted.mean(b, w= c(AUC_GLM, AUC_GAM, AUC_MXT, AUC_BART))

   plot(ensemble_r, col = clrs, range = c(0, 1), main = "ENSEMBLE r")
   
   # Ensemble "favourability"
   glm_fav <- as(glm_fav, "SpatRaster")
   gam_fav <- as(gam_fav, "SpatRaster")
   
    b <- brick(x = c(glm_fav, gam_fav, mxt_p, bart_fav))
   
    ensemble_f <- weighted.mean(b, w= c(AUC_GLM, AUC_GAM, AUC_MXT, AUC_BART))
  #  ensemble_f <- weighted.mean(b, w= c(AUC_GAM, AUC_MXT, AUC_BART))
   #ensemble_f <- weighted.mean(b, w= c(AUC_GAM, AUC_BART))
  # ensemble_f <- weighted.mean(b, w= c(AUC_GAM, AUC_MXT))
   plot(ensemble_f, col = clrs, range = c(0, 1), main = "ENSEMBLE f")
   

   
  ## End ensemble
   
   # Save rasters as tif
   writeRaster(ensemble_r,paste0(path_TIF, "/", (26+ i), "r.tif"),options=c('TFW=YES'))
   writeRaster(ensemble_f,paste0(path_TIF, "/", (26+ i), "f.tif"),options=c('TFW=YES'))
   
   # Save rasters as jpeg to make a video
   jpeg(paste0(path_JPEG, "/", (40 - i), "p.jpeg"),  units="in", width=9, height=7, res=300)
   

   dev.off()
   
  
 }
 

 
 #### NICHE BREADTH ####
 
 
 library(ENMTools)
 library(geosphere)
 my_path<-" " # write your directory
 path_TIF <- paste0(my_path, "/", species, "/Layers/TIF")
 
 
 
 # Get the .tif files
 f <- list.files(path_TIF,
                 pattern = "\\.tif$",
                 recursive = FALSE,
                 all.files = FALSE,
                 full.names = TRUE)
 # Select the favourability files
 f <- f[grep("f\\.tif$", f)]
 
 # Get the numbers
 file_numbers <- as.integer(gsub("[^0-9]", "", basename(f)))
 
 # Order the numbers
 order_index <- order(file_numbers, decreasing = FALSE)
 
 # Order the f list accordingly 
 f_sorted <- f[order_index]
 
 # Get and Store data
 resultados <- list()
 
 
 
 # First, we create a new data frame to save the outputs obtained from the niche breadth computations
 Time_k_BP <- NA
 Output <- data.frame(Time_k_BP)
 Output$Experiment <- NA
 Output$NichBrd_1_f <- NA
 Output$NichBrd_2_f<- NA
 Output$Favourable_area<- NA
 Output$Favourable_area_abs<- NA
 Output$Area_km<- NA
 Output$Species <- species
 
 # Now run the loop
 for (i in 1:length(f_sorted)) {
   archivo_raster <- f_sorted[i]
   
   # Carga el raster desde el archivo
   raster <- raster::raster(archivo_raster)
   plot(raster)
   
   ensemble_f <- raster
   # ensemble_f <- weighted.mean(b, w= c(AUC_GAM, AUC_MXT, AUC_BART))
   # plot(ensemble_f, range = c(0, 1), main = "ENSEMBLE f")
   ensemble_f<- as(ensemble_f, "SpatRaster")
   NicheBreadth_f<- raster.breadth(ensemble_f,verbose=FALSE) # Measures the spatial heterogeneity of the distribution of suitability scores
   
   
   Time <- as.integer(gsub("[^0-9]", "", basename(archivo_raster)))
   bloque <- 29
   bloque_actual <- ((i - 1) %/% bloque) + 1
   Experiment <- bloque_actual
   NichBrd_1_f <- NicheBreadth_f$B1
   NichBrd_2_f<- NicheBreadth_f$B2
   raster <- ensemble_f >= 0.5 # Pixels with favourability values =>0.5 = presence.
   plot(raster)
   absence <- sum(raster[] == 0, na.rm = TRUE)
   Favourable_area<- sum(raster[] == 1, na.rm = TRUE) / sum(raster[] == 0, na.rm = TRUE)
   Favourable_area_abs<- sum(raster[] == 1, na.rm = TRUE)
   Species <- species
   # area in km
   r <- raster(raster)
   polys <- rasterToPolygons(r, fun=function(x){x==1}, dissolve=TRUE)
   Area_km <- sum(areaPolygon(polys)) / 1E6  # en km2
   
   
   #Save outputs
   Output[nrow(Output) + 1,] <- c(Time, Experiment, NichBrd_1_f, NichBrd_2_f, 
                                  Favourable_area, Favourable_area_abs, Area_km, Species)
   
 }
 
 
 print(Output)
 # save data
 library(xlsx)
 library(ggpubr)
 getwd()
 write.xlsx(Output, "Output.xlsx")
 
 # PLOT ANICHE BREADTH
 
 # Write your directory (...)
 dataframe_sensitivity <- read.xlsx("...Output_Sensitivity.xlsx", rowNames=FALSE,
                                    colNames=TRUE, sheet = "Hoja1")
 head(dataframe_sensitivity)
 
 
 Favourable_area <-ggplot(data = dataframe_sensitivity, aes(x = Time_k_BP, y = Area_km / 1000000, fill=Species, group = interaction(Time_k_BP, Species))) +
   geom_boxplot(position = position_dodge(width = 0.8), color = "black", width = 1, outlier.size = 1, outlier.color = "grey") +
   theme_classic() +
   #ylim(1, 3) +
   xlab("Age (kya BP)") +
   ylab(expression(paste("Area "("10" ^{6}, paste("km"^{-2})))))+
   labs(fill = "Species") +
   ggtitle(expression(italic("Homo sapiens")))+
   scale_x_reverse()+
   scale_y_continuous(limits = c(1, 3), breaks = seq(0.5, 3, by = 0.5))+
   theme(plot.title = element_text(hjust = 0.5))
 
 Favourable_area
 
 
 # Levins' B metric of niche breadth
 Niche_breadth_1_Plot<- ggplot(dataframe_sensitivity, aes(x = Species, y = NichBrd_1_f, fill=Species)) +
   geom_boxplot() +
   labs(x = "Species", y = "Levin's B1", title = "Niche breadth")+
   theme_classic()+
   stat_compare_means()
 Niche_breadth_1_Plot
 
 Niche_breadth_2_Plot<- ggplot(dataframe_sensitivity, aes(x = Species, y = NichBrd_2_f, fill=Species)) +
   geom_boxplot() +
   labs(x = "Species", y = "Levin's B2", title = "Niche breadth")+
   theme_classic()+
   stat_compare_means()
 Niche_breadth_2_Plot
 
 grid.arrange(Niche_breadth_1_Plot, Niche_breadth_2_Plot, ncol=2)
 
 # Summary statistics
 
 Neand_df<- subset(dataframe_sensitivity, Species=="Homo_neand")
 head(Neand_df)
 Neand_55_50_df<- subset(Neand_df, Time_k_BP >= 50 & Time_k_BP <= 55)
 head(Neand_55_50_df)
 mean(Neand_55_50_df$Area_km)
 sd(Neand_55_50_df$Area_km)
 
 Neand_50_45_df<- subset(Neand_df, Time_k_BP >= 45 & Time_k_BP <= 50)
 head(Neand_50_45_df)
 mean(Neand_50_45_df$Area_km)
 sd(Neand_50_45_df$Area_km)
 
 Neand_45_40_df<- subset(Neand_df, Time_k_BP >= 40 & Time_k_BP <= 45)
 head(Neand_45_40_df)
 mean(Neand_45_40_df$Area_km)
 sd(Neand_45_40_df$Area_km)
 
 Neand_40_35_df<- subset(Neand_df, Time_k_BP >= 35 & Time_k_BP <= 40)
 head(Neand_40_35_df)
 mean(Neand_40_35_df$Area_km)
 sd(Neand_40_35_df$Area_km)
 
 Neand_35_30_df<- subset(Neand_df, Time_k_BP >= 30 & Time_k_BP <= 35)
 head(Neand_35_30_df)
 mean(Neand_35_30_df$Area_km)
 sd(Neand_35_30_df$Area_km)
 
 
 Sap_df<- subset(dataframe_sensitivity, Species=="Homo_sapiens")
 head(Sap_df)
 Sap_55_50_df<- subset(Sap_df, Time_k_BP >= 50 & Time_k_BP <= 55)
 head(Sap_55_50_df)
 mean(Sap_55_50_df$Area_km)
 sd(Sap_55_50_df$Area_km)
 
 Sap_50_45_df<- subset(Sap_df, Time_k_BP >= 45 & Time_k_BP <= 50)
 head(Sap_50_45_df)
 mean(Sap_50_45_df$Area_km)
 sd(Sap_50_45_df$Area_km)
 
 Sap_45_40_df<- subset(Sap_df, Time_k_BP >= 40 & Time_k_BP <= 45)
 head(Sap_45_40_df)
 mean(Sap_45_40_df$Area_km)
 sd(Sap_45_40_df$Area_km)
 
 Sap_40_35_df<- subset(Sap_df, Time_k_BP >= 35 & Time_k_BP <= 40)
 head(Sap_40_35_df)
 mean(Sap_40_35_df$Area_km)
 sd(Sap_40_35_df$Area_km)
 
 Sap_35_30_df<- subset(Sap_df, Time_k_BP >= 30 & Time_k_BP <= 35)
 head(Sap_35_30_df)
 mean(Sap_35_30_df$Area_km)
 sd(Sap_35_30_df$Area_km)
 
 
 
 
 
 
## MAPS WITH THE HABITAT FAVOURABILITY DIFFERENCE BETWEEN SAPIENS AND NEANDERTHALS

library(RColorBrewer)
library(raster)
mi_paleta <- colorRampPalette(brewer.pal(11, "RdYlBu"))
clrs <- rev(mi_paleta(100))

getwd()

species <- "Homo.sapiens"

# Change my_path if necessary
my_path<- paste0(my_path, "/", species, "/Layers/TIF")
layers <- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)


Sapiens_55_50 <- c(
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/50f.tif",
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/51f.tif",
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/52f.tif",
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/53f.tif",
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/54f.tif",
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/55f.tif"
)
Sapiens_55_50 <- mean(stack(Sapiens_55_50))
plot(Sapiens_55_50)

Sapiens_50_45 <- c(
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/45f.tif",
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/46f.tif",
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/47f.tif",
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/48f.tif",
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/49f.tif"
)
Sapiens_50_45 <- mean(stack(Sapiens_50_45))

Sapiens_45_40 <- c(
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/40f.tif",
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/41f.tif",
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/42f.tif",
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/43f.tif",
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/44f.tif"
)
Sapiens_45_40 <- mean(stack(Sapiens_45_40))

Sapiens_40_35 <- c(
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/35f.tif",
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/36f.tif",
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/37f.tif",
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/38f.tif",
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/39f.tif"
)
Sapiens_40_35 <- mean(stack(Sapiens_40_35))

Sapiens_35_30 <- c(
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/30f.tif",
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/31f.tif",
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/32f.tif",
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/33f.tif",
  "F:/SPD_NP_v2_final/Species/Homo.sapiens/Layers/TIF/34f.tif"
)
Sapiens_35_30 <- mean(stack(Sapiens_35_30))



species <- "Homo_neand"

# Change my_path if necessary
my_path<- paste0(my_path, "/", species, "/Layers/TIF")
layers <- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)


Neanderthals_55_50 <- c(
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/50f.tif",
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/51f.tif",
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/52f.tif",
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/53f.tif",
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/54f.tif",
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/55f.tif"
)
Neanderthals_55_50 <- mean(stack(Neanderthals_55_50))


Neanderthals_50_45 <- c(
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/45f.tif",
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/46f.tif",
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/47f.tif",
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/48f.tif",
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/49f.tif"
)
Neanderthals_50_45 <- mean(stack(Neanderthals_50_45))

Neanderthals_45_40 <- c(
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/40f.tif",
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/41f.tif",
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/42f.tif",
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/43f.tif",
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/44f.tif"
)
Neanderthals_45_40 <- mean(stack(Neanderthals_45_40))

Neanderthals_40_35 <- c(
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/35f.tif",
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/36f.tif",
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/37f.tif",
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/38f.tif",
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/39f.tif"
)
Neanderthals_40_35 <- mean(stack(Neanderthals_40_35))

Neanderthals_35_30 <- c(
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/30f.tif",
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/31f.tif",
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/32f.tif",
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/33f.tif",
  "F:/SPD_NP_v2_final/Species/Homo_neand/Layers/TIF/34f.tif"
)
Neanderthals_35_30 <- mean(stack(Neanderthals_35_30))



### DIFFERENCE BETWEEN SAPIENS AND NEANDERTHALS
# Define the extent for plotting
library(dplyr)
library(gridExtra)
library(ggplot2)
mywindow <- c(-15, 60, 33, 75)  # xmin, xmax, ymin, ymax
#x11()
# Extract the subset of the world map data
world_map_subset <- map_data("world") %>%
  filter(long >= mywindow[1] & long <= mywindow[2] & lat >= mywindow[3] & lat <= mywindow[4] &
           region != "Morocco" & region != "Algeria" & region != "Tunisia")

### DIFFERENCE BETWEEN SAPIENS AND NEANDERTHALS

 library(rasterVis)
library(ggplot2)
library(RColorBrewer)    
library(openxlsx)
getwd()
Transitional_industries_df<-read.xlsx("F:/TFM_Noelia/Script/Raw.data_filtered_fixed.xlsx", rowNames=FALSE, 
                         colNames=TRUE, sheet="Fauna")
head(Transitional_industries_df)

# Select only sites transitional industries in each time-period:

##################### Between 55 and 50 kya

Chatelperronian_df_55_50<- subset(Transitional_industries_df, Chatelperronian == "1"& Age <= 55000 & Age >= 50000)

Uluzzian_df_55_50<- subset(Transitional_industries_df, Uluzzian == "1" & Age <= 55000 & Age >= 50000)

Neronian_df_55_50<- subset(Transitional_industries_df, Neronian == "1" & Age <= 55000 & Age >= 50000)

Bohunician_df_55_50<- subset(Transitional_industries_df, Bohunician == "1" & Age <= 55000 & Age >= 50000)

Szeletian_df_55_50<- subset(Transitional_industries_df, Szeletian == "1"& Age <= 55000 & Age >= 50000)

LRJ_df_55_50<- subset(Transitional_industries_df, LRJ == "1"& Age <= 55000 & Age >= 50000)

##################### Between 50 and 45 kya

Chatelperronian_df_50_45<- subset(Transitional_industries_df, Chatelperronian == "1"& Age <= 50000 & Age >= 45000)

Uluzzian_df_50_45<- subset(Transitional_industries_df, Uluzzian == "1" & Age <= 50000 & Age >= 45000)

Neronian_df_50_45<- subset(Transitional_industries_df, Neronian == "1" & Age <= 50000 & Age >= 45000)

Bohunician_df_50_45<- subset(Transitional_industries_df, Bohunician == "1" & Age <= 50000 & Age >= 45000)

Szeletian_df_50_45<- subset(Transitional_industries_df, Szeletian == "1"& Age <= 50000 & Age >= 45000)

LRJ_df_50_45<- subset(Transitional_industries_df, LRJ == "1"& Age <= 50000 & Age >= 45000)

##################### Between 45 and 40 kya

Chatelperronian_df_45_40<- subset(Transitional_industries_df, Chatelperronian == "1"& Age <= 45000 & Age >= 40000)

Uluzzian_df_45_40<- subset(Transitional_industries_df, Uluzzian == "1" & Age <= 45000 & Age >= 40000)

Neronian_df_45_40<- subset(Transitional_industries_df, Neronian == "1" & Age <= 45000 & Age >= 40000)

Bohunician_df_45_40<- subset(Transitional_industries_df, Bohunician == "1" & Age <= 45000 & Age >= 40000)

Szeletian_df_45_40<- subset(Transitional_industries_df, Szeletian == "1"& Age <= 45000 & Age >= 40000)

LRJ_df_45_40<- subset(Transitional_industries_df, LRJ == "1"& Age <= 45000 & Age >= 40000)

##################### Between 40 and 35 kya

Chatelperronian_df_40_35<- subset(Transitional_industries_df, Chatelperronian == "1"& Age <= 40000 & Age >= 35000)

Uluzzian_df_40_35<- subset(Transitional_industries_df, Uluzzian == "1" & Age <= 40000 & Age >= 35000)

Neronian_df_40_35<- subset(Transitional_industries_df, Neronian == "1" & Age <= 40000 & Age >= 35000)

Bohunician_df_40_35<- subset(Transitional_industries_df, Bohunician == "1" & Age <= 40000 & Age >= 35000)

Szeletian_df_40_35<- subset(Transitional_industries_df, Szeletian == "1"& Age <= 40000 & Age >= 35000)

LRJ_df_40_35<- subset(Transitional_industries_df, LRJ == "1"& Age <= 40000 & Age >= 35000)

##################### Between 35 and 30 kya

Chatelperronian_df_35_30<- subset(Transitional_industries_df, Chatelperronian == "1"& Age <= 35000 & Age >= 30000)

Uluzzian_df_35_30<- subset(Transitional_industries_df, Uluzzian == "1" & Age <= 35000 & Age >= 30000)

Neronian_df_35_30<- subset(Transitional_industries_df, Neronian == "1" & Age <= 35000 & Age >= 30000)

Bohunician_df_35_30<- subset(Transitional_industries_df, Bohunician == "1" & Age <= 35000 & Age >= 30000)

Szeletian_df_35_30<- subset(Transitional_industries_df, Szeletian == "1"& Age <= 35000 & Age >= 30000)

LRJ_df_35_30<- subset(Transitional_industries_df, LRJ == "1"& Age <= 35000 & Age >= 30000)

##################### NOW PLOT

# Colors
mi_paleta <- colorRampPalette(brewer.pal(11, "RdYlBu"))
                                                 
# par(mfrow = c(2, 3))

Dif_55_50 <- Neanderthals_55_50 - Sapiens_55_50
#writeRaster(Dif_40_35, "F:/TFM_Noelia/Resultados/Dif_55_50")
raster_stack_dif <- stack(Dif_55_50)
raster_df_s_dif <- as.data.frame(raster_stack_dif, xy = TRUE)

Dif_55_50<-ggplot(raster_df_s_dif, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), limits = c(-0.5, 0.4), na.value = "white") +
  theme_classic() +
  ylab("Latitude (º)")+
  xlab("Longitude (º)")+
  ggtitle("Mean favourability difference between 55-50 kya")+
  coord_cartesian(xlim = c(-10, 30), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", size=0.1)+
  geom_point(data = Chatelperronian_df_55_50, aes(x = x, y = y), size = 3, color="black", inherit.aes = FALSE)+
  geom_point(data = Uluzzian_df_55_50, aes(x = x, y = y), size = 3, color="brown", inherit.aes = FALSE)+
  geom_point(data = Neronian_df_55_50, aes(x = x, y = y), size = 3, color="grey", inherit.aes = FALSE)+
  geom_point(data = Bohunician_df_55_50, aes(x = x, y = y), size = 3, color="purple", inherit.aes = FALSE)+
  geom_point(data = Szeletian_df_55_50, aes(x = x, y = y), size = 3, color="pink", inherit.aes = FALSE)+
  geom_point(data = LRJ_df_55_50, aes(x = x, y = y), size = 3, color="green", inherit.aes = FALSE)

Dif_55_50





Dif_50_45 <-  Neanderthals_50_45 -Sapiens_50_45

raster_stack_dif <- stack(Dif_50_45)
raster_df_s_dif <- as.data.frame(raster_stack_dif, xy = TRUE)

Dif_50_45<-ggplot(raster_df_s_dif, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), limits = c(-0.5, 0.4), na.value = "white") +
  theme_classic() +
  ylab("Latitude (º)")+
  xlab("Longitude (º)")+
  ggtitle("Mean favourability difference between 50-45 kya")+
  coord_cartesian(xlim = c(-10, 30), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", size=0.1)+
  geom_point(data = Chatelperronian_df_50_45, aes(x = x, y = y), size = 3, color="black", inherit.aes = FALSE)+
  geom_point(data = Uluzzian_df_50_45, aes(x = x, y = y), size = 3, color="brown", inherit.aes = FALSE)+
  geom_point(data = Neronian_df_50_45, aes(x = x, y = y), size = 3, color="grey", inherit.aes = FALSE)+
  geom_point(data = Bohunician_df_50_45, aes(x = x, y = y), size = 3, color="purple", inherit.aes = FALSE)+
  geom_point(data = Szeletian_df_50_45, aes(x = x, y = y), size = 3, color="pink", inherit.aes = FALSE)+
  geom_point(data = LRJ_df_50_45, aes(x = x, y = y), size = 3, color="green", inherit.aes = FALSE)


Dif_50_45



Dif_45_40 <-  Neanderthals_45_40 -Sapiens_45_40

raster_stack_dif <- stack(Dif_45_40)
raster_df_s_dif <- as.data.frame(raster_stack_dif, xy = TRUE)

Dif_45_40<-ggplot(raster_df_s_dif, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), limits = c(-0.5, 0.4), na.value = "white") +
  theme_classic() +
  ylab("Latitude (º)")+
  xlab("Longitude (º)")+
  ggtitle("Mean favourability difference between 45-40 kya")+
  coord_cartesian(xlim = c(-10, 30), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", size=0.1)+
  geom_point(data = Chatelperronian_df_45_40, aes(x = x, y = y), size = 3, color="black", inherit.aes = FALSE)+
  geom_point(data = Uluzzian_df_45_40, aes(x = x, y = y), size = 3, color="brown", inherit.aes = FALSE)+
  geom_point(data = Neronian_df_45_40, aes(x = x, y = y), size = 3, color="grey", inherit.aes = FALSE)+
  geom_point(data = Bohunician_df_45_40, aes(x = x, y = y), size = 3, color="purple", inherit.aes = FALSE)+
  geom_point(data = Szeletian_df_45_40, aes(x = x, y = y), size = 3, color="pink", inherit.aes = FALSE)+
  geom_point(data = LRJ_df_45_40, aes(x = x, y = y), size = 3, color="green", inherit.aes = FALSE)


Dif_45_40


Dif_40_35 <-  Neanderthals_40_35 -Sapiens_40_35

raster_stack_dif <- stack(Dif_40_35)
raster_df_s_dif <- as.data.frame(raster_stack_dif, xy = TRUE)

Dif_40_35<-ggplot(raster_df_s_dif, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)),limits = c(-0.5, 0.4), na.value = "white") +
  theme_classic() +
  ylab("Latitude (º)")+
  xlab("Longitude (º)")+
  ggtitle("Mean favourability difference between 40-35 kya")+
  coord_cartesian(xlim = c(-10, 30), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", size=0.1)+
  geom_point(data = Chatelperronian_df_40_35, aes(x = x, y = y), size = 3, color="black", inherit.aes = FALSE)+
  geom_point(data = Uluzzian_df_40_35, aes(x = x, y = y), size = 3, color="brown", inherit.aes = FALSE)+
  geom_point(data = Neronian_df_40_35, aes(x = x, y = y), size = 3, color="grey", inherit.aes = FALSE)+
  geom_point(data = Bohunician_df_40_35, aes(x = x, y = y), size = 3, color="purple", inherit.aes = FALSE)+
  geom_point(data = Szeletian_df_40_35, aes(x = x, y = y), size = 3, color="pink", inherit.aes = FALSE)+
  geom_point(data = LRJ_df_40_35, aes(x = x, y = y), size = 3, color="green", inherit.aes = FALSE)


Dif_40_35



Dif_35_30 <-  Neanderthals_35_30 -Sapiens_35_30
#writeRaster(Dif_40_35, "F:/TFM_Noelia/Resultados/Dif_40_35")
raster_stack_dif <- stack(Dif_35_30)
raster_df_s_dif <- as.data.frame(raster_stack_dif, xy = TRUE)

Dif_35_30<-ggplot(raster_df_s_dif, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), limits = c(-0.5, 0.4), na.value = "white") +
  theme_classic() +
  ylab("Latitude (º)")+
  xlab("Longitude (º)")+
  ggtitle("Mean favourability difference between 35-30 kya")+
  coord_cartesian(xlim = c(-10, 30), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", size=0.1)+
  geom_point(data = Chatelperronian_df_35_30, aes(x = x, y = y), size = 3, color="black", inherit.aes = FALSE)+
  geom_point(data = Uluzzian_df_35_30, aes(x = x, y = y), size = 3, color="brown", inherit.aes = FALSE)+
  geom_point(data = Neronian_df_35_30, aes(x = x, y = y), size = 3, color="grey", inherit.aes = FALSE)+
  geom_point(data = Bohunician_df_35_30, aes(x = x, y = y), size = 3, color="purple", inherit.aes = FALSE)+
  geom_point(data = Szeletian_df_35_30, aes(x = x, y = y), size = 3, color="pink", inherit.aes = FALSE)+
  geom_point(data = LRJ_df_35_30, aes(x = x, y = y), size = 3, color="green", inherit.aes = FALSE)


Dif_35_30


grid.arrange(Dif_55_50, Dif_50_45, Dif_45_40, Dif_40_35, Dif_35_30)



library(grDevices)

# Generar el gráfico
Dif_50_45 <- ggplot(...) + ...

# Extraer los datos subyacentes del gráfico
plot_data <- ggplot_build(Dif_50_45)$data[[1]]

# Obtener los colores de los puntos
point_colors <- plot_data$colour

# Mapear los colores a tus categorías
categories <- character(length(point_colors))

for (i in seq_along(point_colors)) {
  if (point_colors[i] == "black") {
    categories[i] <- "Chatelperronian"
  } else if (point_colors[i] == "brown") {
    categories[i] <- "Uluzzian"
  } else if (point_colors[i] == "grey") {
    categories[i] <- "Neronian"
  } else if (point_colors[i] == "purple") {
    categories[i] <- "Bohunician"
  } else if (point_colors[i] == "pink") {
    categories[i] <- "Szeletian"
  } else if (point_colors[i] == "green") {
    categories[i] <- "LRJ"
  } else {
    categories[i] <- "Otra"
  }
}


plot_data$categoria <- categories

print(length(categories))
print(length(plot_data$colour))

getwd()

# Replace "..." with your directory

df <- read.xlsx("...Df_transitional.xlsx", rowNames=FALSE,
                                   colNames=TRUE, sheet = "Hoja1")
head(df)

ggplot(df, aes(x = Culture, y = Dif, fill = Culture)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.1), size = 3) +  # Agregar puntos con jitter
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Añadir línea horizontal en 0
  labs(x = "Culture", y = "Diferencia (%)", title = "") +
  theme_minimal()


ggplot(df, aes(x = Culture, y = Dif, fill = Culture)) +
  geom_jitter(position = position_jitter(width = 0.2), size = 3, alpha = 0.6) +  # Gráfico de dispersión con jitter
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Añadir línea horizontal en 0
  labs(x = "Culture", y = "Dif", title = "Diferencia por culturas") +
  theme_minimal()

ggplot(df, aes(x = Culture, y = Dif, fill = Culture, color = Age)) +
 # geom_boxplot() +
  geom_violin()+
  geom_point(position = position_jitter(width = 0.1), size = 3) +  # Agregar puntos con jitter
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Añadir línea horizontal en 0
  labs(x = "Culture", y = "Dif", title = "Diferencia por culturas") +
  theme_minimal() +
  scale_color_manual(values = c("35_30" = "blue", "40_35" = "green", "45_40" = "orange", "50_45" = "purple", "55_50" = "red"))  # Especificar la paleta de colores


