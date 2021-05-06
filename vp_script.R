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

# Separate the categorical and ordinal data
env_factors <- select(red.env, PERIOD, TEXTURE, OWNERSHIP, DISTURB, INTENSITY)

# Separate the numerical data
env_numeric <- select(red.env, -PERIOD, -TEXTURE, -OWNERSHIP, -DISTURB, -INTENSITY)

# Scale the data
env_numeric <- as.data.frame(scale(env_numeric))

# Bring them back together 
env.red.scaled <- cbind(env_factors, env_numeric)

# Rearrange the data so types (e.g., soil, landscape, overstory) are grouped
env.red.scaled <- select(env.red.scaled, PERIOD, OWNERSHIP, DISTURB, INTENSITY, TEXTURE, everything())

# Convert intensity into ordinal 
env.red.scaled$INTENSITY <- factor(env.red.scaled$INTENSITY, 
														order = TRUE, levels = c("NONE", "LOW", "MED", "HIGH"))
# Check for collinearity ####