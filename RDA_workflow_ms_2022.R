#Adapted from the "Detecting multilocus adaptation using Redundancy Analysis (RDA)" vignette by Forester et al., 2018
  #https://popgen.nescent.org/2018-03-27_RDA_GEA.html

#Genetic data:
#a matrix of genotypic data; sampling alleles in an individual or allele frequencies in a population

#Spatial modelling of the distribution of samples is currently done very efficiently using Moran Eigenvector 
#Maps (MEM). These allow you to model small-scale to large-scale spatial patterns in your sampling sites,
#and uses this matrix as a set of explanatory variables in a RDA.

###########################################################################

#For population-based data, you can input the genomic data as allele frequencies within demes

# install.packages(c("psych","vegan"), dependencies=TRUE)

# Load packages
# -------------
library(psych)    # Used to investigate correlations among predictors
library(vegan)    # Used to run RDA

par(mfrow = c(1, 1))

Gev_gen <- read.csv("genotypic_data_matrix.csv")

Gev_gen$X <- as.character(Gev_gen$X) #pop names as character

#Make sure no missing data! RDA must have complete data frames
sum(is.na(Gev_gen)) #0

#Read in and screen envr predictors:
  #pop and ecotype columns need to have numeric categorical values
env <- read.csv("environmental_predictor_data.csv")

env$pop <- as.factor(env$pop) #use numbers to designate pops and ecotypes
env$ecotype <- as.factor(env$ecotype)

#Confirm that genotypes and environmental data are in the same order
#identical(rownames(gen.imp), env[,1]) 

#heat map
cor.plot(env[,5:26], colors = TRUE, show.legend = TRUE)

#with heat map using raw and adjusted probabilities of correlations
rs <- corr.test(env[,5:26])
rs

corPlot(r=rs$r,numbers=TRUE,pval=rs$p,main="Correlations scaled by probability values") 

#Show the upper and lower confidence intervals and remove highly correlated (Predictor variables with R2 > 0.75 and < -0.75 
  #(Schweizer et al., 2016; Forester et al., 2018) predictors)
cor.plot.upperLowerCi(R=rs,numbers=TRUE) 

pred <- subset(env, select=-c(worldclim.1, worldclim.2, worldclim.4, worldclim.5, worldclim.6, worldclim.7, worldclim.8,
worldclim.11, worldclim.12, worldclim.13, worldclim.16, worldclim.17, worldclim.18, worldclim.19, elev)) 
 
#investigate the distribution of factor levels in land_cover
#table(pred$land_cover) 
table(pred$ecotype) #only pop and ecotype columns are factors

#check reduced predictor set:
pred <- pred[,5:11]

colnames(pred) <- c("iso","tdq","twq","prd","prs", "MEM1", "MEM2")
pairs.panels(pred, scale=T)

################### RUN the RDA ###########################
#Note: if your predictors include any factors, you’ll need to write out the formula in the rda call 

ev.17pop.rda <- rda(Gev_gen[,-c(1)] ~ ., data=pred, scale=T) #no character columns (popnames, etc)
ev.17pop.rda

RsquareAdj(ev.17pop.rda)

summary(eigenvals(ev.17pop.rda, model = "constrained")) 

screeplot(ev.17pop.rda)

#check RDA model (full model and each contstrained axis) for significance using formal tests 
  #(F-stats(Legendre et al, 2010))

signif.full <- anova.cca(ev.17pop.rda, parallel=getOption("mc.cores")) # default is permutation=999
signif.full

#next check each constrained axis for significance

signif.axis <- anova.cca(ev.17pop.rda, by="axis", parallel=getOption("mc.cores"))
signif.axis

#Check variance inflation factors (VIF) of predictor variables used in model
 #If variables have inflation factors > 20 (Ter Braak, 1988; Pienitz et al., 1995), remove variables and re-run RDA
  
vif.cca(ev.17pop.rda)

################# PLOT RDA ############################
par(mar=c(5, 7, 1, 2) + 0.1)
plot(ev.17pop.rda, scaling=3) #default plots axes 1 and 2

## ECOTYPE LEVEL PLOT######
levels(env$ecotype) <- c("EC", "South", "North")
eco <- env$ecotype
bg <- c("#ff7f00","#1f78b4","#ffff33") # 3 colors for the ecotypes  

par(mfrow = c(1, 1))
plot(ev.17pop.rda, type="n", scaling=3)
points(ev.17pop.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(ev.17pop.rda, display="sites", pch=21, cex=1.7, col="gray32", scaling=3, bg=bg[eco]) # Ev ecotypes
text(ev.17pop.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("bottomright", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

### POP LEVEL PLOT ###################
levels(env$pop) <- c("EC", "NC", "VM", "CC", "EL", "NN", "GO", "CF", "ST", "TB", "CH", "AT", "TL", "AN",
                     "SG", "CP", "PB")
evpop <- env$pop
bg <- c("#ff7f00","#1f78b4","#ffff33",'#a6cee3', '#6a3d9a', '#e31a1c','#33a02c','#fb9a99','magenta', 
        'lightgoldenrod2', 'maroon', 'cyan2', 'greenyellow', 'grey47', 'black', 'coral4', 'palegreen1')
    #colors for 17 pops

##below plot script has larger filled squares to denote pops in legend
par(mfrow = c(1, 1))
plot(ev.17pop.rda, type="n", scaling=3)
points(ev.17pop.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(ev.17pop.rda, display="sites", pch=21, cex=1.7, col="gray32", scaling=3, bg=bg[evpop]) # the plants by pop
text(ev.17pop.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("topleft", legend=evpop, bty="n", col="gray32", cex=0.5, 
       fill=bg, pt.bg=bg)

###### ID CANDIDATE SNPS INVOLVED IN LOCAL ADAPTATION ######
#The SNP loadings are stored as species in the RDA object. SNP loadings are extracted only from significant
  #constrained axis/axes
  
load.rda <- scores(ev.17pop.rda, choices=c(1), display="species")  

#look at histograms of the loadings on each RDA axis, to check for (relatively) normal distributions. 

hist(load.rda[,1], main="Loadings on RDA1")
 
#ID SNPs in tail of distribution using a 2.0 standard deviation cutoff to id as many candidates as possible 
  #(even those under weak selection) and accept the risk of false positives 

#outliers function, where x is the vector of loadings and z is the number of standard deviations to use
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

#apply the sd cutoff to each significant constrained axis:

cand1 <- outliers(load.rda[,1],2)

ncand <- length(cand1)
ncand 

#organize results by making one data frame with the axis, SNP name, loading, 
#& correlation with each predictor:

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))

colnames(cand1) <- c("axis","snp","loading")

cand <- rbind(cand1)

cand$snp <- as.character(cand$snp)

#add in the correlations of each candidate SNP with the environmental predictors

foo <- matrix(nrow=(ncand), ncol=7)
colnames(foo) <- c("iso","tdq", "twq", "prd","prs", "MEM1","MEM2")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- Gev_gen[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)


#### INVESTIGATE THE CANDIDATES ########
#look for duplicate candidates across axes if more than one axis is significant

#length(cand$snp[duplicated(cand$snp)])
#foo <- cbind(cand$axis, duplicated(cand$snp)) 
#table(foo[foo[,1]==1,2])
#table(foo[foo[,1]==2,2])
#cand <- cand[!duplicated(cand$snp),] # remove duplicate detections

###Now look at which of the predictors each candidate SNP is most strongly correlated with:

for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,11] <- names(which.max(abs(bar[4:10]))) # gives the variable [if have 7 variables, 4:10 excludes axis, snp, loading cols]
  cand[i,12] <- max(abs(bar[4:10]))              # gives the correlation
}

colnames(cand)[11] <- "predictor"
colnames(cand)[12] <- "correlation"

table(cand$predictor)
write.csv(cand, file = "cand_pred_corr.csv") 

#Note that, in some cases, correlations may be strong for multiple variables, so should explore SNP
   #correlations to each predictor

#### PLOT THE SNPs ##########
#RDAs with SNPs colored based on predictor variable most correlated with

sel <- cand$snp
env <- cand$predictor

env[env=="iso"] <- '#a6cee3'
env[env=="MEM1"] <- '#1f78b4'
env[env=="tdq"] <- '#33a02c'
env[env=="twq"] <- 'deeppink1'
env[env=="prd"] <- '#fb9a99'
env[env=="prs"] <- '#b2df8a'
env[env=="MEM2"] <- 'yellow1'

#color by predictor:

col.pred <- rownames(ev.17pop.rda$CCA$v) # pull the SNP names

for (i in 1:length(sel)) {           # color code candidate SNPs
  foo2 <- match(sel[i],col.pred)
  col.pred[foo2] <- env[i]
}

col.pred[grep("X", col.pred)] <- '#f1eef6' # non-candidate SNPs (SNPs start with X, not chr (as in tut.))
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")

bg <- c('#a6cee3','#1f78b4','#33a02c','deeppink1','#fb9a99','#b2df8a', 'yellow1') 

plot(ev.17pop.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(ev.17pop.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
points(ev.17pop.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(ev.17pop.rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("topright", legend=c("iso", "tdq", "twq", "prd","prs","MEM1", "MEM2"), 
       bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

