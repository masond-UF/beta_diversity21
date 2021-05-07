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
red.env.numeric <- as.data.frame(select(red.env, -PERIOD, -TEXTURE, 
																		-OWNERSHIP, -DISTURB, -INTENSITY))

# Separate the categorical and ordinal data
red.env.factors <- select(red.env, PERIOD, TEXTURE, OWNERSHIP, DISTURB, INTENSITY)

# Convert intensity into ordinal 
red.env.factors$INTENSITY <- factor(red.env$INTENSITY, 
														order = TRUE, levels = c("NONE", "LOW", "MED", "HIGH"))

# Check for collinearity ####
# Correlations among continuous variables
red.env.cont.corr <- as.data.frame(scale(env_numeric)) %>% 
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
env.temp <- cbind(red.env.factors, scale(red.env.numeric))

# Drop the collinear variables
env.final <- env.temp[which(names(env.temp) %!in% drop.list)]	

