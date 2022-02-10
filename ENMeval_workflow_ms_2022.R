####Adapted from ENMeval 2.0.0 Vignette by Jamie M. Kass, Robert Muscarella, Gonzalo E. Pinilla-Buitrago, 
#and Peter J. Galante
#https://jamiemkass.github.io/ENMeval/articles/ENMeval-2.0.0-vignette.html#references-and-resources

options(java.parameters = "-Xmx200g") #need to run this before loading packages to avoid low Java memory error
library(ENMeval) #Need ENMeval version >= 2.0
library(raster)
library(dplyr)

# Set a random seed in order to be able to reproduce this analysis.
set.seed(48)

####### USAGE #############
occs <- read.csv("occurrence_coordinates.csv")

# Removing occurrences to avoid pseudoreplication.
occs <- occs[!duplicated(occs),]

#AICc (Akaike Information Criterion for small sample size) to compare models

#USE CLIPPED WORLDCLIM BIOVARIABLES RASTERS, (600 MB is too large for a single raster .asc file)
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

#use maxent.jar algorithm to get variable contribution to model statistics
enmeval <- ENMevaluate(occs, env, bg.coords = NULL, occ.grp = NULL, 
            bg.grp = NULL, RMvalues = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4), 
            fc = c("L", "LQ", "LQP", "H", "T", "LQH", "LQHP", "LQHPT"),
            categoricals = NULL, n.bg = 20000, method = "block", 
            algorithm = 'maxent.jar', overlap = FALSE, 
            bin.output = FALSE, clamp = TRUE, rasterPreds = NULL, 
            parallel = FALSE, numCores = 1, progbar = TRUE, 
            updateProgress = FALSE)

saveRDS(enmeval, "enmeval.RDS")


#################Exploring ENMeval results ###################################
library(tidyverse)

enmeval <- readRDS("enmeval.RDS")
enmeval
str(enmeval, max.level=2)
eval.results(enmeval) %>% head()
eval.bg(enmeval) %>% head()
eval.occs.grp(enmeval) %>% str()

#Write .csvs of results and variable.importance for enmeval object
write.csv(enmeval@results, file = "enmeval_v2_results.csv")
#enmeval@variable.importance
write.csv(enmeval@variable.importance, file = "enmeval_v2_results_var_import.csv")

#####CREATE CSV OF WORLDCLIM VARIABLE VALUES FOR ALL OCCURRENCES AND CREATE VARIABLE SUBSET########:
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

#Examine variable collinearity and re-run ENMeval with subset of variables: variables contributing most to 
  #model are retained and highly correlated variables are removed 
library (ENMTools)
raster.cor.matrix(env, method = "pearson")

#### Visualizing tuning results #####################

#to plot more than one statistic at once with ggplot facetting:
evalplot.stats(e = enmeval, stats = c("or.mtp", "auc.val"), color = "fc", x.var = "rm")
#to jitter the positions of overlapping points:
evalplot.stats(e = enmeval, stats = c("or.mtp", "auc.val"), color = "fc", x.var = "rm", 
               dodge = 0.5)
#to switch which variables are on the x-axis and which are symbolized by color.
# (ENMeval only accepts two variables for plotting)
evalplot.stats(e = enmeval, stats = c("or.mtp", "auc.val"), color = "rm", x.var = "fc", 
               error.bars = FALSE)

########## Model selection and overall results ##############
res <- eval.results(enmeval)
# Select the model with delta AICc equal to 0, or the one with the lowest AICc score.
opt.aicc <- res %>% filter(delta.AICc == 0)
opt.aicc

#cross-validation sequential method to select models with the lowest average test omission rate, and if needed, 
#to break ties using highest average validation AUC (Radosavljevic & Anderson 2014, Kass et al. 2020).
opt.seq <- res %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))
opt.seq

# Select the model predictions for our optimal model (based on preferred method) and plot 
pred.seq <- eval.predictions(enmeval)[[opt.aicc$tune.args]]
plot(pred.seq)

#plot the binned background points below occurrence points on top to 
# visualize locations of training data
points(eval.bg(enmeval), pch = 3, col = eval.bg.grp(enmeval), cex = 0.5)
points(eval.occs(enmeval), pch = 21, col = "white", bg = eval.occs.grp(enmeval))


########## METADATA ###############################
#create a rangeModelMetadata object (described by Merow et al. (2019)) to describe standards of metadata and ENM rating systems
rmm <- eval.rmm(enmeval)
#describe model selection approach
rmm$model$selectionRules <- "lowest delta AICc"
#define optimal model settings
rmm$model$finalModelSettings <- "rm.0.5_fc.LQHPT"
rmm$prediction$continuous$minVal <- cellStats(pred.seq, min)
rmm$prediction$continuous$maxVal <- cellStats(pred.seq, max)
rmm$prediction$continuous$units <- "suitability (cloglog transformation)"
#save the metadata to a CSV file
rangeModelMetadata::rmmToCSV(rmm, "rmm_evlg2022.csv")
