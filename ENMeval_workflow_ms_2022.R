####Adapted from ENMeval 2.0.0 Vignette by Jamie M. Kass, Robert Muscarella, Gonzalo E. Pinilla-Buitrago, 
#and Peter J. Galante
#https://jamiemkass.github.io/ENMeval/articles/ENMeval-2.0.0-vignette.html#references-and-resources

options(java.parameters = "-Xmx200g") #need to run this before loading packages to avoid low Java memory error
library(ENMeval) #Need ENMeval version >= 2.0
library(raster)
library(dplyr)

####### USAGE #############
occ <- read.csv("occurrence_coordinates.csv")

# Removing occurrences to avoid pseudoreplication.
occ <- occ[!duplicated(occ),]

#AICc (Akaike Information Criterion for small sample size) to compare models

#USE CLIPPED WORLDCLIM BIOVARIABLES RASTERS, (600mb is too large for a single raster .asc file)
r1 <- raster("bio_01_clip.asc")
r2 <- raster("bio_02_clip.asc")
r3 <- raster("bio_03_clip.asc")
r4 <- raster("bio_04_clip.asc")
r5 <- raster("bio_05_clip.asc")
r6 <- raster("bio_06_clip.asc")
r7 <- raster("bio_07_clip.asc")
r8 <- raster("bio_08_clip.asc")
r9 <- raster("bio_09_clip.asc")
r10 <- raster("bio_10_clip.asc")
r11 <- raster("bio_11_clip.asc")
r12 <- raster("bio_12_clip.asc")
r13 <- raster("bio_13_clip.asc")
r14 <- raster("bio_14_clip.asc")
r15 <- raster("bio_15_clip.asc")
r16 <- raster("bio_16_clip.asc")
r17 <- raster("bio_17_clip.asc")
r18 <- raster("bio_18_clip.asc")
r19 <- raster("bio_19_clip.asc")
env <- stack(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19)

#use maxent.jar algorithm to get variable contribution to model stats
enmeval <- ENMevaluate(occ, env, bg.coords = NULL, occ.grp = NULL, 
            bg.grp = NULL, RMvalues = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4), 
            fc = c("L", "LQ", "LQP", "H", "T", "LQH", "LQHP", "LQHPT"),
            categoricals = NULL, n.bg = 20000, method = "block", 
            algorithm = 'maxent.jar', overlap = FALSE, 
            bin.output = FALSE, clamp = TRUE, rasterPreds = NULL, 
            parallel = FALSE, numCores = 1, progbar = TRUE, 
            updateProgress = FALSE)

saveRDS(enmeval, "enmeval.RDS")

#####COMMAND-LINE R TO CREATE CSV OF WORLDCLIM VARIABLE VALUES FOR ALL EV OCCURRENCES AND IDENTIFY CORRELATED PREDICTORS ########:
library(sf)
#env #raster stack
evpres <- st_read("ev_pres_AK_layers/POINT.shp") #Worldclim shapefile

#create a csv with worldclim data for each occurrence
ev_wc_spatial <- extract(x = env, y = as(evpres, "Spatial"), df = TRUE ) #doesn't work on laptop
write.csv(ev_wc_spatial, "ev_pres_AK_wc_all_occs_08262021.csv")

library(psych)
env_csv <- read.csv("ev_pres_AK_wc_all_occs_08262021.csv")
rs <- corr.test(env_csv[,1:])
rs

#Use ENMTools to examine variable collinearity and then re-run ENMeval with subset of variables
library (ENMTools)
raster.cor.matrix(env, method = "pearson")


#################Exploring ENMeval results ###################################
library(dismo)
library(sp)
library(tidyverse)

enmeval <- readRDS("enmeval.RDS")
enmeval
str(enmeval, max.level=2)

enmeval@algorithm
enmeval@models
enmeval@partition.method
enmeval@results %>% head()
enmeval@models %>% str(max.level = 1)
enmeval@results
write.csv(enmeval@results, file = "enmeval_results.csv")

#investigate AICc, AUC_test and orMTP to evaluate model performance
eval.plot(enmeval@results, value = "delta.AICc", variance = NULL, legend = TRUE,
          legend.position = "topright")
eval.plot(enmeval@results, value = "AICc", variance = NULL, legend = TRUE,
          legend.position = "topright")
eval.plot(enmeval@results, value = "avg.test.AUC", variance = NULL, legend = TRUE,
          legend.position = "topright")
eval.plot(enmeval@results, value = "avg.diff.AUC", variance = NULL, legend = TRUE,
          legend.position = "topright")
eval.plot(enmeval@results, value = "avg.test.orMTP", variance = NULL, legend = TRUE,
          legend.position = "topright")