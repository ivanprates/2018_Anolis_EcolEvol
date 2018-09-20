### This script's goal is to run PCA analyses based on genetic (using the LEA package) and environmental data.
### By Ivan Prates and Anna Penna, June 2018.
### Smithsonian National Museum of Natural History, Washington DC, USA.

## Part 1: Getting ready

# Load packages
library(ggplot2)
library(LEA)
library(magrittr)
library(plyr)
library(viridis)
library(cowplot)

## Running PCA on the 

# Set target species
#sp = "ortonii"
sp = "punctatus"

# Set all versus candidate SNPs
snps = "all"
#snps = "candidate"

# Settings for each species
if (sp == "ortonii") {
  n = 23
  myPalette = c("#440154", "#fde725") # Color palette for plots
} else {
  n = 46
  myPalette = c("#fde725", "#35b779", "#440154")
}

# Set minimum read length
t = 70

# Set maximum number of SNPs per locus
s = 10

# Pick general path, depending on which computer I'm working at:
path = "~/Dropbox (Smithsonian)/ivan_lab/2018_Anolis_GEA/2018-03/"

# Settings to import each SNP dataset
if (snps == "all") {
  gendata = read.table(paste0(path, "data/VCFtools_LFMM_", sp, "_t", t, "_s", s, "_n", n, "/", sp, "_t", t, "_s", s, "_n", n, "_MAF10_LFMM.012"), sep = "\t", row.names = 1, header = FALSE) # read genotype file from VCFtools
  individuals = read.table(paste0(path, "data/VCFtools_LFMM_", sp, "_t", t, "_s", s, "_n", n, "/", sp, "_t", t, "_s", s, "_n", n, "_MAF10_LFMM.012.indv"), sep = "\t", header = FALSE) # read genotype file from VCFtools
  row.names(gendata) = individuals$V1
  } else {
  gendata = read.table(paste0(path, "data/", sp, "_t", t, "_s", s, "_n", n, "_candidate_SNPs/", sp, "_t", t, "_s", s, "_n", n, "_candidate_SNPs.012"), sep = "\t", row.names = 1) # read genotype file, candidate loci only
}

# Check dimensions of the genetic data
dim(gendata)

# Create a directory for the outputs of LEA analyses
dir.create(paste0(path, "LEA/PCA_", sp, "_t", t, "_s", s, "_n", n, "_", snps, "_SNPs"))

# Set this new directory as the working directory
setwd(paste0(path, "LEA/PCA_", sp, "_t", t, "_s", s, "_n", n,  "_", snps, "_SNPs"))

# Write down data in LFMM format (LEA saves it all outside R, instead of as objects)
write.lfmm(gendata, paste0(sp, "_genotypes.lfmm"))

##  Part 2: Performing PCA on the genetic data

# Create a pcaProject object: pc
pc = pca(paste0(sp,"_genotypes.lfmm"), center = TRUE, scale = FALSE)

# Display information on analysis
show(pc)

# Summarize analysis
summary(pc)

# Plot eigenvalues
#plot(pc, lwd=5, col="blue", cex = .7, xlab=("Factors"), ylab="Eigenvalues")

# Plot standard deviations
#plot(pc$sdev)

## Perfom Tracy-Widom tests for all eigenvalues

# Create file: genotypes.tracyWidom - tracy-widom test information, in the directory genotypes.pca/
#tw <- tracy.widom(pc)

# Plot the percentage of variance explained by each component
#plot(tw$percentage)

# Show p-values for the Tracy-Widom tests
#tw$pvalues

# Plotting genetic PCA axes
pcadata = pc$projections # This is what we'll plot (the pc axes)
pcadata = as.data.frame(pcadata)
pcadata$ID = rownames(gendata) # Adding names of samples based on "gendata" object

# Assigning samples to populations
if (sp == "ortonii") {
  pcadata$pop = c("pop1","pop1","pop1","pop1","pop1","pop1","pop1","pop1","pop1","pop1","pop1","pop1","pop1","pop2","pop2","pop2","pop1","pop1","pop1","pop2","pop2","pop1","pop1")
  } else {
  pcadata$pop = c("pop2","pop2","pop1","pop1","pop1","pop3","pop3","pop3","pop3","pop3","pop1","pop1","pop2","pop1","pop2","pop1","pop2","pop1","pop2","pop2","pop2","pop2","pop2","pop2","pop3","pop3","pop3","pop3","pop3","pop1","pop2","pop2","pop3","pop1","pop2","pop2","pop2","pop3","pop3","pop2","pop2","pop3","pop3","pop3","pop1","pop1")
}

# Plot with ggplot: genetic PC1 vs. genetic PC2
gen_gen_plot = pcadata %>% ggplot() + geom_point(aes(x = V1, y = V2, fill = pop), size = 8, shape = 21) +
  #geom_text(aes(x = V1, y = V2, label = ID)) + # adds individual labels
  scale_fill_manual(values = myPalette[c(3:1)], na.value = "grey80") +
  theme_light() +
  labs(y = "Genetic PC2", x = "Genetic PC1", title = paste0("Anolis ", sp, ", ", snps, " SNPs")) +
  theme(plot.title = element_text(size = 30, hjust = 0.5),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 22),
        panel.grid = element_blank())
gen_gen_plot

# Save plot with cowplot
save_plot(filename = paste0(sp, "_gen_PC1_vs_gen_PC2_", snps, "_SNPs.pdf"), plot = gen_gen_plot, base_width = 15, base_height = 10)

# Part 3: Performing PCA on the environmental data

## Importing and selecting environmental data

# Import environmental data
# We're using 46 samples of A. punctatus and 23 of A. ortonii
all.envdata = read.csv(paste0(path, "data/Bioclim/Bioclim_01-19_46punc.csv"), header = TRUE)
rownames(all.envdata) = all.envdata$ID_number # change row names to specimen IDs

# Extract columns corresponding to the 19 Bioclim environmental variables
variables = names(all.envdata[c(8:26)])
envdata = as.matrix(all.envdata[all.envdata$species == paste0("Anolis_", sp), variables[c(1:19)]]) # select 'rows, columns' corresponding to a given 'species, climate-variable' combination

# We are using PCs extracted from the following selected Bioclim variables:
sel.envdata = envdata[,variables[c(1,5,10,11,12,13,16)]] # v7
n.var = 7 # Number of climatic variables; this is just to name output files

# Test for correlations between selected variables
corr.env = cov2cor(cov(sel.envdata))
abs(corr.env) >= 0.7

# Implement PCA based on all variables
pca = prcomp(x = sel.envdata, retx = TRUE, center = TRUE, scale. = TRUE)
pca
pca$x
summary(pca)
t(pca$rotation)*pca$sdev # to check the correlations among the principal components' scores and the initial data (correlation matrix) scores and the initial data (correlation matrix)

# Choose PC variables to use
envdata.used = pca$x[,1] # PC1 from temperature and precipitation PCA
envdata.name = "pc1" # This is just to name output files
pcadata$pc1 = envdata.used

# Plot with ggplot: environmental PC1 vs. genetic PC1
gen_env_plot = pcadata %>% ggplot() + geom_point(aes(x = pc1, y = V1, fill = pop), size = 8, shape = 21) +
  scale_fill_manual(values = myPalette[c(3:1)], guide = "none", na.value = "grey80") +
  theme_light() +
  labs(y = "Genetic PC1", x = "Environmental PC1", title = paste0("Anolis ", sp, ", ", snps, " SNPs")) +
  theme(plot.title = element_text(size = 30, hjust = 0.5),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 22),
        panel.grid = element_blank())

# Save plot with cowplot
save_plot(filename = paste0(sp, "_gen_PC1_vs_env_PC1_", snps, "_SNPs.pdf"), plot = gen_env_plot, base_width = 15, base_height = 10)

# Done!
