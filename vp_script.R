# Install packages and load data ####
library(tidyverse)
library(vegan)
library(ade4)
library(adespatial)
library(spdep)
library(betapart)
library(spdep)
library(corrplot)
library(SoDA) 
library(rgdal)

options(scipen = 999)


# Bring in the data ####
spec <- read.csv("data/species.csv")
head(spec)

env <- read.csv("data/env_coords.csv")
head(env)
# Selecting species ####
source("code/biostats.r") # attach the code from biostats

occur <- foa.plots(spec)
rare <- which(occur[,2]<5)
common <- which(occur[,2]>95)

red_spec <- spec[,-c(rare,common)]

# Transform data ####
spe.hel <- decostand(spec, "hellinger") 
red.spe.hel <- decostand(red_spec, "hellinger")
# Preparing environmental data ####

# Drop site, date, spatial and texture proportions
red.env <- env[,-c(1:4,16,17,18)] # 

# Separate the numerical data
env.numeric <- as.data.frame(select(red.env, -PERIOD, -TEXTURE, 
																		-OWNERSHIP, -DISTURB, -INTENSITY))

# Separate the categorical and ordinal data
env.factors <- select(red.env, PERIOD, TEXTURE, OWNERSHIP, DISTURB, INTENSITY)

# Convert intensity into ordinal 
env.factors$INTENSITY <- factor(red.env$INTENSITY, 
														order = TRUE, levels = c("NONE", "LOW", "MED", "HIGH"))

# Check for collinearity ####
# Correlations among continuous variables
red.env.cont.corr <- as.data.frame(scale(env.numeric)) %>% 
								 		cor(method = c("kendall"))

red.env.cont.corr.df <- as.data.frame(as.table(red.env.cont.corr))
red.env.cont.corr.filt <- red.env.cont.corr.df %>%  
															arrange(desc(Freq)) %>% 
															filter(Freq>0.7)

corrplot(red.env.cont.corr, tl.cex = 0.5)

numeric.drop <- c("Morus.rubra", "Cercis.canadensis", "Ilex.decidua",
														 "S", "Asimina.triloba", "Carya.pallida")

# Correlations among categorical and continuous variables

# DIST
red.env.factor.corr <- red.env %>% 
							select(-PERIOD, -TEXTURE, -OWNERSHIP, -INTENSITY) %>% 
							select(DISTURB, everything())

DIST <- red.env.factor.corr[,1]
red.env.factor.corr <- red.env.factor.corr[,-1]


rsq <- vector()
pseudo.r <- vector()
for(i in 1:84){
	mod <- lm(formula = red.env.factor.corr[,i] ~ DIST, 
					data = red.env.factor.corr)
	rsq <- summary(mod)$r.squared
	pseudo.r[i] <- sqrt(rsq)
}

pseudo.r <- as.data.frame(as.table(pseudo.r)) 
pseudo.r$Var1 <- as.vector(dimnames(red.env.factor.corr)[[2]]) 
# No variables correlated with disturbance type (DIST) > 0.7

# PERIOD
red.env.factor.corr <- red.env %>% 
	select(-DISTURB, -TEXTURE, -OWNERSHIP, -INTENSITY) %>% 
	select(PERIOD, everything())

PERIOD <- red.env.factor.corr[,1]
red.env.factor.corr <- red.env.factor.corr[,-1]
rsq <- vector()
pseudo.r <- vector()
for(i in 1:84){
	mod <- lm(formula = red.env.factor.corr[,i] ~ PERIOD, 
					data = red.env.factor.corr)
	rsq <-  summary(mod)$r.squared
	pseudo.r[i] <- sqrt(rsq)
}

pseudo.r <- as.data.frame(as.table(pseudo.r)) 
pseudo.r$Var1 <- as.vector(dimnames(red.env.factor.corr)[[2]]) 
# soilNA correlated with sampling period >0.7

# TEXTURE
red.env.factor.corr <- red.env %>% 
	select(-PERIOD, -DISTURB, -OWNERSHIP, -INTENSITY) %>% 
	select(TEXTURE, everything())

TEXTURE <- red.env.factor.corr[,1]
red.env.factor.corr <- red.env.factor.corr[,-1]
rsq <- vector()
pseudo.r <- vector()
for(i in 1:84){
	mod <- lm(formula = red.env.factor.corr[,i] ~ TEXTURE, 
					data = red.env.factor.corr)
	rsq <-  summary(mod)$r.squared
	pseudo.r[i] <- sqrt(rsq)
}

pseudo.r <- as.data.frame(as.table(pseudo.r)) 
pseudo.r$Var1 <- as.vector(dimnames(red.env.factor.corr)[[2]]) 
# Quercus.pagoda correlated with soil texture >0.7

# OWNERSHIP
red.env.factor.corr <- red.env %>% 
	select(-PERIOD, -DISTURB, -TEXTURE, -INTENSITY) %>% 
	select(OWNERSHIP, everything())

OWNERSHIP <- red.env.factor.corr[,1]
red.env.factor.corr <- red.env.factor.corr[,-1]
rsq <- vector()
pseudo.r <- vector()
for(i in 1:84){
	mod <- lm(formula = red.env.factor.corr[,i] ~ OWNERSHIP, 
					data = red.env.factor.corr)
	rsq <-  summary(mod)$r.squared
	pseudo.r[i] <- sqrt(rsq)
}

pseudo.r <- as.data.frame(as.table(pseudo.r)) 
pseudo.r$Var1 <- as.vector(dimnames(red.env.factor.corr)[[2]]) 
# DEVHI correlated with land ownership  >0.7

# INTENSITY
red.env.factor.corr <- red.env %>% 
	select(-PERIOD, -DISTURB, -TEXTURE, -OWNERSHIP) %>% 
	select(INTENSITY, everything())

INTENSITY <- red.env.factor.corr[,1]
red.env.factor.corr <- red.env.factor.corr[,-1]
rsq <- vector()
pseudo.r <- vector()
for(i in 1:84){
	mod <- lm(formula = red.env.factor.corr[,i] ~ INTENSITY, 
					data = red.env.factor.corr)
	rsq <-  summary(mod)$r.squared
	pseudo.r[i] <- sqrt(rsq)
}

pseudo.r <- as.data.frame(as.table(pseudo.r)) 
pseudo.r$Var1 <- as.vector(dimnames(red.env.factor.corr)[[2]]) # nothing
# No variables correlated with disturbance intensity >0.7
# Creating the final environmental matrix

# List of species correlated with categorical variables
categorical.drop <- c("soilNA", "Quercus.pagoda", "DEVHI")

# Combined list of correlated variables to drop
drop.list <- c(numeric.drop, categorical.drop)
# Drop the correlated variables
'%!in%' <- function(x,y){ 
	!('%in%'(x,y)) # make a function to do the opposite of %in%
}

# Bring the data back together
env.temp <- cbind(env.factors, scale(env.numeric))

# Drop the collinear variables
env.final <- env.temp[which(names(env.temp) %!in% drop.list)]	
# Variable selection ####
# Full species matrix and environmental variables = 13% explained
spe.rda <- rda(spe.hel ~ ., env.final)
full.step.forward <- ordiR2step(rda(spe.hel ~ 1, data=env.final), 
scope=formula(spe.rda), R2scope = F, direction="forward", pstep=1000)

# Reduced species matrix and environmental variables = 16% explained
red.spe.rda <- rda(red.spe.hel ~ ., env.final)
red.step.forward <- ordiR2step(rda(red.spe.hel ~ 1, data=env.final), 
scope=formula(red.spe.rda), R2scope = F, direction="forward", pstep=1000)

# Create the PCNM variables ####
spatial.tmp <- as.data.frame(env[,2:3]) # extract spatial variables for later analysis
spatial.tmp <- rename(spatial.tmp, longitude = east_x, latitude = north_y)
spatial <- geoXY(spatial.tmp$longitude, spatial.tmp$latitude)

# Is there a linear trend in the mite data?
anova(rda(red.spe.hel, spatial)) # Result: significant trend
# Computation of linearly detrended mite data
red.spe.hel.det <- resid(lm(as.matrix(red.spe.hel) ~ ., data = as.data.frame(spatial)))

## Step 1. Construct the matrix of dbMEM variables
dbmem.tmp <- dbmem(spatial, silent = FALSE)
dbmem <- as.data.frame(dbmem.tmp)

# Truncation distance used above:
thr <- give.thresh(dist(spatial))

# Display and count the eigenvalues
attributes(dbmem.tmp)$values 
length(attributes(dbmem.tmp)$values) # 18 eigenvalues

## Step 2. Run the global dbMEM analysis on the detrended
## Hellinger-transformed reduced species data
dbmem.rda <- rda(red.spe.hel.det ~., dbmem)
anova(dbmem.rda) # model is significant


## Step 3. Since the R-square is significant, compute the adjusted
## R2 and run a forward selection of the dbmem variables 
R2a <- RsquareAdj(dbmem.rda)$adj.r.squared
dbmem.fwd <- forward.sel(red.spe.hel, as.matrix(dbmem),
												 adjR2thresh = R2a)

nb.sig.dbmem <- nrow(dbmem.fwd) # 3 signif. dbMEM

# Identity of the significant dbMEM in increasing order
dbmem.sign <- sort(dbmem.fwd[ ,2])

# Write the significant dbMEM to a new object 
dbmem.red <- dbmem[ ,c(dbmem.sign)]
# Visualizing the spatial descriptors ####

## Step 1. Construct the matrix of dbMEM variables
dbmem.viz.tmp <- dbmem(spatial, silent = FALSE)
dbmem.viz <- as.data.frame(dbmem.tmp)

# Truncation distance used above:
thr <- give.thresh(dist(spatial))

# Display and count the eigenvalues
attributes(dbmem.viz.tmp)$values 
length(attributes(dbmem.viz.tmp)$values) 

## Step 2. Run the global dbMEM analysis on the detrended
## Hellinger-transformed reduced species data
dbmem.viz.rda <- rda(red.spe.hel ~., dbmem.viz)
anova(dbmem.viz.rda) # model is significant


## Step 3. Since the R-square is significant, compute the adjusted
## R2 and run a forward selection of the dbmem variables 
R2a.viz <- RsquareAdj(dbmem.viz.rda)$adj.r.squared
dbmem.viz.fwd <- forward.sel(red.spe.hel, as.matrix(dbmem.viz),
												 adjR2thresh = R2a.viz)

nb.sig.dbmem.viz <- nrow(dbmem.viz.fwd) # 3 signif. dbMEM

# Identity of the significant dbMEM in increasing order
dbmem.viz.sign <- sort(dbmem.viz.fwd[ ,2])

# Write the significant dbMEM to a new object 
dbmem.viz.red <- dbmem.viz[ ,c(dbmem.viz.sign)]

## Step 4. New dbMEM analysis with 8 significant dbMEM variables ## Adjusted R-square after forward selection: R2adj = 0.2418
dbmem.viz.rda2 <- rda(red.spe.hel ~ ., data = dbmem.viz.red) 
fwd.viz.R2a <- RsquareAdj(dbmem.viz.rda2)$adj.r.squared
anova(dbmem.viz.rda2)

axes.test <- anova(dbmem.viz.rda2, by = "axis") 

## Step 5. Plot the significant canonical axes
dbmem.viz.rda2.axes <- 
	scores(dbmem.viz.rda2, 
				 choices = c(1:nb.sig.dbmem.viz), 
				 display = "lc", 
				 scaling = 1) 

source("code/sr.value.R") # attach the code from biostats

par(mfrow = c(1,nb.sig.dbmem.viz))
for(i in 1:nb.sig.dbmem.viz){
	sr.value(spatial, dbmem.viz.rda2.axes[ ,i],
	sub = paste("RDA",i), csub = 2)
}

# Variance partitioning ####
spatiotemporal <- cbind(dbmem.red, env$PERIOD)
landscape <- select(env, OWNERSHIP, WATER.150, DECID, DEVOP, CROP, GRASS)
overstory <- select(env, CANOPY, Pinus.taeda, Quercus.alba)
soil.ground <- select(env, TEXTURE, MG, PH, ROCK)

red.spe.part <- varpart(red.spe.hel, spatiotemporal, landscape,
												overstory, soil.ground)

dev.off()
plot(red.spe.part, digits = 2, Xnames = c('Space & time', 'Landscape', 
		 'Overstory', 'Soil & ground'), id.size = 0.75, bg = 2:5)
# Run anova on fractions ####
env.final.xy <- cbind(env.final, dbmem.red)

rda.all <- rda(red.spe.hel ~MEM2+MEM4+MEM5+MEM15+PERIOD+
							 	OWNERSHIP+WATER.150+DECID+DEVOP+CROP+GRASS+
							 	CANOPY+Pinus.taeda+Quercus.alba+TEXTURE+MG+
							 	PH+ROCK, data = env.final.xy)
anova(rda.all)

rda.soil.ground <- rda(red.spe.hel ~ TEXTURE + 
											 	MG + PH + ROCK, data = env.final.xy)
anova(rda.soil.ground)

rda.overstory <- rda(red.spe.hel ~ CANOPY + Quercus.alba + 
					 					Pinus.taeda, data = env.final.xy)
anova(rda.overstory)

rda.landscape <- rda(red.spe.hel ~ DECID + DEVOP + CROP + 
								WATER.150 + GRASS + OWNERSHIP, data = env.final.xy)
anova(rda.landscape)

rda.space.time <- rda(red.spe.hel ~ PERIOD + MEM2 + MEM4 +
												MEM5 + MEM15, data = env.final.xy)
anova(rda.space.time)







