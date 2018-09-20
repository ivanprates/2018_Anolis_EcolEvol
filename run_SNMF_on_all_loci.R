### This script's goal is to run genetic clustering analyses (using SNMF) solely for candidate SNPs recovered by LFMM.
### It will also make bar plots based on the ancestry coefficients.
### By Ivan Prates and Anna Penna, June 2018.
### Smithsonian National Museum of Natural History, Washington DC, USA.

## Part 1: Getting ready

# Load packages
library(cowplot)
library(dplyr)
library(ggplot2)
library(LEA)
library(magrittr)
library(plyr)
library(tidyr)

# Set target species
#sp = "ortonii"
sp = "punctatus"

# Settings for each species (this is very specific to my dataset, e.g., how I named my individual samples)
if (sp == "ortonii") {
  n = 23 # number of sampled individuals
  sp.short = "orto_" # That's how I named species during GBS library sequencing
  myPalette = c("#440154", "#fde725") # Color palette for plots
  bestK = 2 # best K, defined based on results (added here just to make final plots)
  a = 100 # SNMF regularization parameter (alpha), defined based on results (added here just to make final plots)
  } else {
  n = 46 # number of sampled individuals
  sp.short = "punc_" # That's how I named species during GBS library sequencing
  myPalette = c("#fde725", "#35b779", "#440154") # Color palette for plots
  bestK = 3 # best K, defined based on results (added here just to make final plots)
  a = 100 # SNMF regularization parameter (alpha), defined based on results (added here just to make final plots)
  }

# Set minimum read length 
# This is a parameter in ipyrad. I included parameter values in the name of ipyrad output files - this is why this is here)
t = 70

# Set maximum number of SNPs per locus
# This is a parameter in ipyrad. I included parameter values in the name of ipyrad output files - this is why this is here)
s = 10

# Pick general path, depending on which computer I'm working at:
path = "~/Dropbox (Smithsonian)/ivan_lab/2018_Anolis_GEA/2018-03/"

# Create a directory to place resulting SNMF bar plots
dir.create(paste0(path, "LEA/plots_SNMF_all_SNPs"))

# Create a directory for the outputs of SNMF analyses. Notice that I name directories by indicating parameter settings!
dir.create(paste0(path, "LEA/SNMF_", sp, "_t", t, "_s", s, "_n", n, "_a", a, "_all_SNPs"))

# Set this new directory as the working directory
setwd(paste0(path, "LEA/SNMF_", sp, "_t", t, "_s", s, "_n", n, "_a", a, "_all_SNPs"))

# Read genetic data
# Missing sites should be coded as 9! Don't forget to replace -1 to 9 if using vcf tools.
gendata = read.table(paste0(path, "data/VCFtools_SNMF_", sp, "_t", t, "_s", s, "_n", n, "/", sp, "_t", t, "_s", s, "_n", n, "_1SNP-locus.012"), sep = "\t", row.names = 1) # read genotype file from vcftools
dim(gendata)

# Write down data in LFMM format (LEA saves it all outside R, instead of as objects)
write.lfmm(gendata, paste0(sp, "_genotypes.lfmm"))

## Part 2: Running genetic clustering analyses

# Run sNMF
project.snmf = snmf(paste0(sp, "_genotypes.lfmm"), K = 2:5, alpha = a, entropy = TRUE, repetitions = 100, project = "new", CPU = 12)

# plot cross-entropy criterion of all runs of the project
png(filename = paste0(path, "LEA/plots_SNMF_all_SNPs/", sp, "_all_SNPs_entropy_plot_t", t, "_s", s, "_n", n, "_a", a, ".png"))
plot(project.snmf, cex = 1.2, col = "blue", pch = 19)
dev.off() # saving plot as .png image

# get the cross-entropy of the 10 runs for K = best K
ce = cross.entropy(project.snmf, K = bestK)
ce

# select the run with the lowest cross-entropy for K = best K
bestrun = which.min(ce)
bestrun

# Entropy value corresponding to best run for best K
e = round(ce[bestrun], digits = 4) # rounding to four digits
e

## Part 3: Plotting ancestry coefficients as a bar plot based on the Q matrix estimated by SNMF

# get the Q matrix, for the best run and K
qmatrix = Q(project.snmf, K = bestK, run = bestrun)
qmatrix = as.data.frame(qmatrix)
colnames(qmatrix) = c(paste0(rep("cluster", bestK), 1:bestK)) # change column names for cluster name. Will repeat "cluster" K times, then paste0 with second element from 1 to K times

# Adding ID names for plotting
individuals = read.table(file = paste0(path, "data/VCFtools_SNMF_", sp, "_t", t, "_s", s, "_n", n, "/", sp, "_t", t, "_s", s, "_n", n, "_1SNP-locus.012.indv")) # read list of sample IDs from vcftools
individuals$V1 %<>% gsub(pattern = sp.short, replacement = "", .) # removes "punc_" or "orto_" from names. "." denotes object (individuals) when using "%<>%"
qmatrix$ID = individuals$V1 # change row names to sample IDs

# Arrange order or samples in plot
sample_info = read.csv(file = paste0(path, "data/plot_order_", sp, "_n", n, ".csv"), header = TRUE)
qmatrix$plot_order = c(sample_info$plot_order)
qmatrix_ord = plyr::arrange(qmatrix, plot_order)

qmatrix = as_tibble(qmatrix) # dplyr likes the tibble format
qmatrix_ord$ID = factor(qmatrix_ord$ID, levels = unique(qmatrix_ord$ID)) # transforming ID into factor to keep order in plot
qmatrix_ord = as.data.frame(qmatrix_ord)

# Save qmatrix
write.csv(qmatrix, file = paste0("qmatrix_", sp, ".csv"), quote = FALSE, row.names = FALSE)

# Load qmatrix (if changes were made out of R)
qmatrix = read.csv(file = paste0("qmatrix_", sp, ".csv"), header = TRUE)

# "Melt" dataframe using the gather function, as required by ggplot
qmatrix_melt = gather(qmatrix_ord, key = cluster, value = coeff, 1:bestK)

# Plot with ggplot
SNMF_plot = qmatrix_melt %>% ggplot(aes(x= ID)) +
  geom_bar(aes(y = coeff, fill = cluster), stat = "identity", position = "fill") +
  scale_fill_manual("Population", values = myPalette[c(1:bestK)], labels = c("Atlantic Forest", "Amazonia 1", "Amazonia 2")) + 
  #scale_fill_manual("Population", values = myPalette[c(1:bestK)], labels = c("Amazonia", "Atlantic Forest")) + 
  scale_color_manual(values = myPalette[c(1:bestK)], guide = "none") +
  theme_minimal(base_size = 16) +
  labs(y = "Ancestrality coefficient", x = "") +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5, size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0))

# Check plot
SNMF_plot

# save plot
save_plot(filename = paste0(path, "LEA/plots_SNMF_all_SNPs/", sp, "_all_SNPs_ancestry_plot_t", t, "_s", s, "_n", n, "_e", e, "_a", a, ".pdf"), plot = SNMF_plot, base_width = 15, base_height = 10)

# Save our progress
save.image(file = paste0(sp, ".RData"))

# Done!
