### R script to define gamma-distributed priors for the theta and tau parameters in G-PhoCS analyses.
### By Ivan Prates, June 2018.
### Smithsonian National Museum of Natural History, Washington DC, USA.

# We'll use the viridis color palletes for plotting:
#install.packages("viridis")
library(viridis)

# Given a mutation rate of:
u = 2.42*10^-9 # from Prates et al. 2016 Molecular Ecology (average mutation rate among all genes of the three lizard species)

# 1. For tau:

# The tree height in number of generations is: 
gen = 1620000 # from Prates et al. 2016 Molecular Ecology (T4 of A. punctatus)
gen = 1120000 # from Prates et al. 2016 Molecular Ecology (T4 of A. ortonii)

# The expected tree height (in coalescent units) is:
r = u * gen 
r # 0.004 for punctatus, 0.0027 for ortonii

# 2. For theta:
N = 1860000 # Between 1.86 million and 4.53 million from Prates et al. 2016 Molecular Ecology
N = 4530000
theta = 4*N*u
theta # From 0.018 to 0.044

# Then, implement gamma parameters as:
shape = c(1,2,5,10)
rate = c(10,20,50)

# Plot:
par(mfrow = c(length(rate), 1)) # number of plot rows = number of rate values
for (j in 1:length(rate)) {
for (i in 1:length(shape)) { # shapes go together, so in inside loop
x.max = qgamma(0.999, shape=shape[i], rate=rate[j])
x = seq(from=0, to=x.max, by=x.max/1000)
dens = dgamma(x, shape=shape[i], rate=rate[j])
palette = viridis(n = length(shape))
plot(x, dens, type='l', ylim = c(0,30), xlim = c(0.003,0.5), col=palette[i], lwd = 2, xlab = "theta")
title(paste("beta =", rate[j])) # titles for each plot
legend(x = 0.35, y = 30 - 3*i, paste("alpha =", shape[i]), text.col=palette[i], bg = "transparent", bty = "n") # settings and position of legend
par(new = T) # plot in same graph
}
par(new = F) # to stop plotting on top
}
par(new = F) # to stop plotting on top

# Done!
