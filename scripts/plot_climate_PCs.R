### R script for plotting climatic PC values extracted at sampled sites versus latitude.
### By Ivan Prates and Anna Penna, June 2018.
### Smithsonian National Museum of Natural History, Washington DC, USA.

# Packages we'll need:
#install.packages("magrittr")
#install.packages("tidyr")
#install.packages("plyr")
#install.packages("dplyr")
#install.packages("psych")

library(LEA)
library(magrittr)
library(tidyr)
library(plyr)
library(dplyr)
library(psych)

setwd("~/Dropbox (Smithsonian)/ivan_lab/2018_Anolis_GEA/2018-03/climate_PCA")

# import environmental data
all.envdata = read.csv("~/Dropbox (Smithsonian)/ivan_lab/2018_Anolis_GEA/2018-03/climate_PCA/bioclim_samples_combined.csv", header = TRUE)
rownames(all.envdata) = all.envdata$ID_number # change row names to specimen IDs

# extract columns corresponding to environmental variables
variables = names(all.envdata[c(13:31)])
envdata = all.envdata[c(13:31)]

## Select variables to keep

# Select more easily interpretable, variable in geographic space, analogous between temp. and precip.
#sel.envdata = envdata[,variables[c(1:19)]] # all

# Select more easily interpretable, variable in geographic space, analogous between temp. and precip.
#sel.envdata = envdata[,variables[c(1,5,6,10,11,12,13,14,16,17)]] # all

# Select as above, but no cold/dry related, only hot/wet
#sel.envdata = envdata[,variables[c(1,5,6,10,11,12,13,14,16,17)]] # v10
sel.envdata = envdata[,variables[c(1,5,10,11,12,13,16)]] # v7
#sel.envdata = envdata[,variables[c(1,5,10,12,13,16)]] # v6
v = 7

# test for correlations between selected variables
corr.env = cov2cor(cov(sel.envdata))
abs(corr.env) >= 0.7

# implement pca, all variables
pca = prcomp(x = sel.envdata, retx = TRUE, center = TRUE, scale. = TRUE)
pca
pca$x
summary(pca)
t(pca$rotation)*pca$sdev # to check the correlations among the  get the correlations among the principal components' scores and the initial data (correlation matrix) scores and the initial data (correlation matrix)

# Implementing PCA with varimax rotation and recovering rotated scores
#pca_rotated = psych::principal(sel.envdata, rotate="varimax", nfactors = 2, scores=TRUE)
#print(pca_rotated$scores)

# Plot PC1 versus latitude
plot(pca$x[,1], all.envdata$latitude, xlab = "PC1", ylab="Latitude", cex.lab=2, cex.axis = 1.5, col="white")
points(pca$x[90:150,1], all.envdata$latitude[90:150], pch=22, cex=3, col="black", bg="grey50") # punctatus Amz, no gen
points(pca$x[188:217,1], all.envdata$latitude[188:217], pch=22, cex=3, col="black", bg="grey50") # punctatus Amz, yes gen
points(pca$x[188:217,1], all.envdata$latitude[188:217], pch=21, cex=0.75, col="black", bg="black") # punctatus Amz, yes gen
points(pca$x[151:187,1], all.envdata$latitude[151:187], pch=22, cex=3, col="black", bg="grey90") # punctatus AF, no gen
points(pca$x[218:233,1], all.envdata$latitude[218:233], pch=22, cex=3, col="black", bg="grey90") # punctatus AF, yes gen
points(pca$x[218:233,1], all.envdata$latitude[218:233], pch=21, cex=0.75, col="black", bg="black") # punctatus AF, yes gen

points(pca$x[1:54,1], all.envdata$latitude[1:54], pch=21, cex=3, col="black", bg="grey50") # ortonii AMz, no gen
points(pca$x[67:84,1], all.envdata$latitude[67:84], pch=21, cex=3, col="black", bg="grey50") # ortonii AMz, yes gen
points(pca$x[67:84,1], all.envdata$latitude[67:84], pch=21, cex=0.75, col="black", bg="black") # ortonii AMz, yes gen
points(pca$x[55:66,1], all.envdata$latitude[55:66], pch=21, cex=3, col="black", bg="grey90") # ortonii AF, no gen
points(pca$x[85:89,1], all.envdata$latitude[85:89], pch=21, cex=3, col="black", bg="grey90") # ortonii AF, yes gen
points(pca$x[85:89,1], all.envdata$latitude[85:89], pch=21, cex=0.75, col="black", bg="black") # ortonii AF, yes gen

# Plot rotated PC1 versus for latitude
#plot(pca_rotated$scores[,1], all.envdata$latitude, xlab="Rotated PC1", ylab="Latitude", pch=20, cex=.8, col="white")
#points(pca_rotated$scores[74:151,1], all.envdata$latitude[74:151], pch=24, cex=2, col="gray30", bg="gray50") # punctatus Amz
#points(pca_rotated$scores[152:172,1], all.envdata$latitude[152:172], pch=24, cex=2, col="black", bg="#C9283E") # punctatus AF
#points(pca_rotated$scores[1:67,1], all.envdata$latitude[1:67], pch=21, cex=2, col="gray30", bg="gray40") # ortonii AMz
#points(pca_rotated$scores[68:73,1], all.envdata$latitude[68:73], pch=21, cex=2, col="black", bg="#4A33E8") # ortonii AF

# Plot rotated PC2 versus for latitude
#plot(pca_rotated$scores[,2], all.envdata$latitude, xlab="Rotated PC2", ylab="Latitude", pch=20, cex=.8, col="white")
#points(pca_rotated$scores[74:151,2], all.envdata$latitude[74:151], pch=24, cex=2, col="gray30", bg="gray50") # punctatus Amz
#points(pca_rotated$scores[152:172,2], all.envdata$latitude[152:172], pch=24, cex=2, col="black", bg="#C9283E") # punctatus AF
#points(pca_rotated$scores[1:67,2], all.envdata$latitude[1:67], pch=21, cex=2, col="gray30", bg="gray40") # ortonii AMz
#points(pca_rotated$scores[68:73,2], all.envdata$latitude[68:73], pch=21, cex=2, col="black", bg="#4A33E8") # ortonii AF