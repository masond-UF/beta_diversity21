# Understory VP—11 Sept 2021
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
options(warn=-1)


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
# Make groups to use for forward selection
time <- env.final %>% select(PERIOD)

soil.ground <- env.final %>% select(TEXTURE,OM,P,K,CA,MG,ZN,PH,LICHEN,
																		LITTER,MOSS,ROCK,ROOT,SOIL,STREAM,
																		TRAMPLING,WOOD,WRACKI,CRUSTACEAN)
landscape <- env.final %>% select(OWNERSHIP,DISTURB,INTENSITY,ELEV,
																	WATER,WATER.20,WATER.150,WATER.500,
																	TRANSPORT,GRASS,MIXFOR,DEVOP,DECID,
																	WETE,PAST,OPENH20,EVRGRN,SHRSCR,
																	BARREN,DEVLO,WETW,CROP,DEVMED)
overstory <- env.final %>% select(CANOPY,Acer.rubrum,Acer.saccharinum,
																	Carpinus.caroliniana,Carya.carolinae.septentrionalis,
																	Carya.glabra,Carya.ovata,Carya.sp,
																	Carya.tomentosa,Cornus.florida,
																	Fagus.grandifolia,Fraxinus.americana,
																	Fraxinus.pennsylvanica,Fraxinus.sp,
																	Hamamelis.virginiana,Liquidambar.styraciflua,
																	Liriodendron.tulipfera,Magnolia.sp,
																	Nyssa.sylvatica,Pinus.taeda,Platanus.occidentalis,
																	Prunus.americana,Prunus.serotina,Quercus.alba,
																	Quercus.coccinea,Quercus.falcata,Quercus.lyrata,
																	Quercus.marilandica,Quercus.nigra,Quercus.phellos,
																	Quercus.stellata,Quercus.velutina,Taxodium.distichum,
																	Ulmus.alata,Ulmus.americana,Ulmus.rubra,Ulmus.sp)

length(time)+length(overstory)+length(soil.ground)+length(landscape)

# soil forward selection
soil.rda <- rda(red.spe.hel ~ ., soil.ground)
soil.step.forward <- ordiR2step(rda(red.spe.hel ~ 1, data=soil.ground), 
scope=formula(soil.rda), R2scope = F, direction="forward", pstep=1000)
# TEXTURE + PH + ZN + CA + ROCK 

# landscape forward selection
landscape.rda <- rda(red.spe.hel ~ ., landscape)
landscape.step.forward <- ordiR2step(rda(red.spe.hel ~ 1, data=landscape), 
scope=formula(landscape.rda), R2scope = F, direction="forward", pstep=1000)
# OWNERSHIP + EVRGRN + WETE + DEVOP + GRASS + WATER.150

# overstory forward selection
overstory.rda <- rda(red.spe.hel ~ ., overstory)
overstory.step.forward <- ordiR2step(rda(red.spe.hel ~ 1, data=overstory), 
scope=formula(overstory.rda), R2scope = F, direction="forward", pstep=1000)
# CANOPY + Pinus.taeda + Quercus.alba + Quercus.stellata + Ulmus.alata
# + Liriodendron.tulipfera + Carya.glabra + Quercus.falcata 

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
space <- dbmem[ ,c(dbmem.sign)]
# Soil Partial RDAs ####
TEXTURE.rda <- rda(red.spe.hel ~ TEXTURE +Condition(PERIOD,OM,P,K,CA,MG,ZN,PH,LICHEN,
																		LITTER,MOSS,ROCK,ROOT,SOIL,STREAM,
																		TRAMPLING,WOOD,WRACKI,CRUSTACEAN,
																		OWNERSHIP,DISTURB,INTENSITY,ELEV,
																		WATER,WATER.20,WATER.150,WATER.500,
																		TRANSPORT,GRASS,MIXFOR,DEVOP,DECID,
																		WETE,PAST,OPENH20,EVRGRN,SHRSCR,
																		BARREN,DEVLO,WETW,CROP,DEVMED,
																		CANOPY,Acer.rubrum,Acer.saccharinum,
																		Carpinus.caroliniana,Carya.carolinae.septentrionalis,
																		Carya.glabra,Carya.ovata,Carya.sp,
																		Carya.tomentosa,Cornus.florida,
																		Fagus.grandifolia,Fraxinus.americana,
																		Fraxinus.pennsylvanica,Fraxinus.sp,
																		Hamamelis.virginiana,Liquidambar.styraciflua,
																		Liriodendron.tulipfera,Magnolia.sp,
																		Nyssa.sylvatica,Pinus.taeda,Platanus.occidentalis,
																		Prunus.americana,Prunus.serotina,Quercus.alba,
																		Quercus.coccinea,Quercus.falcata,Quercus.lyrata,
																		Quercus.marilandica,Quercus.nigra,Quercus.phellos,
																		Quercus.stellata,Quercus.velutina,Taxodium.distichum,
																		Ulmus.alata,Ulmus.americana,Ulmus.rubra,Ulmus.sp) , data=env.final)
anova(TEXTURE.rda)
