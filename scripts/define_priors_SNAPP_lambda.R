# Setting gamma distributions for SNAPP priors

#install.packages("viridis")
library(viridis)

# Setting Yule diversification process (lambda) prior
ntips = 3
rate = 213.675213675 # This is the value we need
height = 0.0039
length = 0.00936

# Then, implement gamma parameters as:
shape = 2
scale = 100

# Plot:
x.max = qgamma(0.999, shape = 2.0, scale = 100.0)
x = seq(from=0, to=x.max, by=x.max/1000)
dens = dgamma(x, shape = 2.0, scale = 100.0)
plot(x, dens, type='l')

median(dgamma(x, shape = 2.0, scale = 100.0))

mean_gam = shape*scale
mean_gam

var_gam = shape*(scale^2)
var_gam
