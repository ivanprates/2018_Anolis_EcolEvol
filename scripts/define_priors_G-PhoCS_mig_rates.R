### R script to define gamma-distributed priors for the migration rate parameters in G-PhoCS analyses.
### By Ivan Prates, June 2018.
### Smithsonian National Museum of Natural History, Washington DC, USA.

# We'll use the viridis color palletes for plotting:
#install.packages("viridis")
library(viridis)

# Given a mutation rate of:
u = 2.42*10^-9 # from Prates et al. 2016 Molecular Ecology (average mutation rate among all genes of the three lizard species)

# G-PhoCS estimates m, given by:
# m = M/u

# To estimate M: 
M = m*u
M # (proportion of ind. in pop. 1 that arose from migration from pop. 2 per generation)

# m for different M
M = 1/100
M = 1/1000
M = 1/1000000
m = M/u
m # from 413 to 4.13M

# For migration, implement gamma parameters as:
shape = c(1)
rate = c(0.0000002)
mean = shape/rate
mean

scale = 1/rate
scale

variance = shape/rate*rate
variance

# Plot:
par(mfrow = c(length(rate), 1)) # number of plot rows = number of rate values
for (j in 1:length(rate)) {
  for (i in 1:length(shape)) { # shapes go together, so in inside loop
    x.max = qgamma(0.999, shape=shape[i], rate=rate[j])
    x = seq(from=0, to=x.max, by=x.max/1000)
    dens = dgamma(x, shape=shape[i], rate=rate[j])
    palette = viridis(n = length(shape))
    plot(x, dens, type='l', xlim = c(0.0,40000000), col=palette[i], lwd = 2, xlab = "migration")
    title(paste("beta =", rate[j])) # titles for each plot
    legend(x = 10, y = 10 - 1*i, paste("alpha =", shape[i]), text.col=palette[i], bg = "transparent", bty = "n") # settings and position of legend
    par(new = T) # plot in same graph
  }
  par(new = F) # to stop plotting on top
}
par(new = F) # to stop plotting on top

# Done!
