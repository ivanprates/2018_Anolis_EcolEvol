### R script for environmental association analyses using LFMM (LEA package).
### By Ivan Prates and Anna Penna, June 2018.
### Smithsonian National Museum of Natural History, Washington DC, USA.

# To get the last version of LEA
#install.packages("devtools")
#devtools::install_github("bcm-uga/LEA")

# Other packages we'll need:
#install.packages("dplyr")
#install.packages("magrittr")
#install.packages("plyr")
#install.packages("tidyr")

# Load packages
library(LEA)
library(dplyr)
library(magrittr)
library(plyr)
library(tidyr)

# Pick general path, depending on which computer I'm working from:
path = "~/Dropbox (Smithsonian)/ivan_lab/2018_Anolis_GEA/2018-03/"

# Set target species
#sp = "ortonii"
sp = "punctatus"

# Set minimum read length 
# This is a parameter in ipyrad. I included parameter values in the name of ipyrad output files - this is why this is here)
t = 70

# Set maximum number of SNPs per locus
# This is a parameter in ipyrad. I included parameter values in the name of ipyrad output files - this is why this is here)
s = 10

# Settings for each species (this is very specific to my dataset, e.g., how I named my individual samples)
if (sp == "ortonii") {
  sp.short = "orto_" # That's how I named species during GBS library sequencing
  n = 23
  bestK = 2 # From sNMF clustering analyses
  K = 2 # Best K for controlling false-discovery rates based on genomic inflation factor (see below)
} else {
  sp.short = "punc_" # That's how I named species during GBS library sequencing
  n = 46
  bestK = 3 # From sNMF clustering analyses
  K = 5 # Best K for controlling false-discovery rates based on genomic inflation factor (see below)
}
  
# Define LFMM run parameters
CPU = 12
repetitions = 10
iterations = 50000
burnin = 25000

## PART 1: Importing and selecting environmental data
  
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

## PART 2: Getting the data ready and a few checkpoints

# Create a directory for the outputs of LEA analyses. Notice that I name directories by indicating parameter settings!
dir.create(paste0(path, "LEA/LFMM_", sp, "_t", t, "_s", s, "_n", n, "_", envdata.name, "_v", n.var, "_k", K))

# Set this new directory as the working directory
setwd(paste0(path, "LEA/LFMM_", sp, "_t", t, "_s", s, "_n", n, "_", envdata.name, "_v", n.var, "_k", K))

# Read genetic data
gendata = read.table(paste0(path, "data/VCFtools_LFMM_", sp, "_t", t, "_s", s, "_n", n, "/", sp, "_t", t, "_s", s, "_n", n, "_MAF10_LFMM.012"), sep = "\t", row.names = 1) # read genotype file from VCFtools
dim(gendata)

## Checkpoint: Do genetic and environmental data (and their order) match?
# Change names of rows in the genotype file for their corresponding sample IDs
individuals = read.table(file = paste0(path, "data/VCFtools_LFMM_", sp, "_t", t, "_s", s, "_n", n, "/", sp, "_t", t, "_s", s, "_n", n, "_MAF10_LFMM.012.indv")) # read list of sample IDs from vcftools
individuals$V1 %<>% gsub(pattern = sp.short, replacement = "", .) # removes "punc_" or "orto_" from names. "." denotes object (individuals) when using "%<>%"
rownames(gendata) = individuals$V1 # change row names to sample IDs
rownames(gendata) == rownames(envdata) # checking whether genetic and environmental data (and their order) match; should all be "TRUE"

# Create list of loci that will serve as a SNP map when retrieving sequences to blast later
loci = read.table(file = paste0(path, "data/VCFtools_LFMM_", sp, "_t", t, "_s", s, "_n", n, "/", sp, "_t", t, "_s", s, "_n", n, "_MAF10_LFMM.012.pos")) # read list of locus IDs from VCFtools
names(loci) = c("locus.name", "SNP.position") # add column names
loci$locus.name %<>% gsub(pattern = "locus_", replacement = "", .) # removes "locus_" from loci names. "." denotes object (loci) when using "%<>%"
loci$SNP.number = 1:nrow(loci) # adding a column with consecutive numbers to the loci object
loci # list loci
length(unique(loci$locus.name)) # How many loci?

## Checkpoint: Do SNPs correspond to the loci that we think they do?
# Change names of columns in the genotype file for their corresponding locus IDs
colnames(gendata) = loci$locus.name # change column names to locus IDs
length(gendata) == length(loci$locus.name) # checking whether number of loci is the same in .pos file and gendata object; should all be "TRUE"

# Write down data in LFMM format (LEA saves it all outside R, instead of as R objects)
write.lfmm(gendata, "genotypes.lfmm")
write.env(envdata.used, "gradients.env")

## Imputation of missing data (not needed if using already imputed data or having no missing sites)

# Run sNMF under the previously inferred best k, because imputation in LEA is based on ancestry coefficients
project.missing = snmf("genotypes.lfmm", K = bestK, entropy = TRUE, repetitions = 20, project = "new")
best = which.min(cross.entropy(project.missing, K = bestK)) # best run from replicates

# Impute missing data based on ancestry coefficients from best sNMF run under best k
impute(project.missing, "genotypes.lfmm", method = "mode", K = bestK, run = best)

# Renaming imputed file to change "." for "_", which for some reason was causing R Studio to abort!!
# This may not be a problem for people not using Linux or using other versions of R and RStudio
file.rename("genotypes.lfmm_imputed.lfmm", "genotypes_lfmm_imputed.lfmm")

## PART 3: Running genome-environment association analyses using LFMM (LEA package)

# Run LFMM; file name has to be the same as that of the imputed data
project = lfmm("genotypes_lfmm_imputed.lfmm", "gradients.env", K = K, repetitions = repetitions, iterations = iterations, burnin = burnin,
               project = "new", CPU = CPU, missing.data = FALSE)

# Save our progress
save.image(file = paste0(sp, ".RData"))

## PART 4: Compute Z scores and adjusted p-values

# Set false discovery rate (FDR) level
q = 0.00001

# Estimating z-scores
zs = z.scores(project, K = K) # Estimate z-scores
zs.median = apply(zs, MARGIN = 1, median) # Estimate median z-score across repetitions

# The median z-scores must be recalibrated before computing p-values, which is done using a rescaling factor, lambda
# Two ways of setting lambda:

# WAY 1: If setting lambda as the Genomic Inflation Factor (GIF), it's value should be close or slightly lower than 1.
# Setting lambda as GIF:
lambda = median(zs.median^2)/qchisq(0.5, df = 1) 

# Yet, it is also important to check the histogram of corrected p-values: It should look flat with a peak close to zero.

# WAY 2: If setting lambda empirically to improve the null distribution:
lambda = 0.45 # Best for punctatus with k = 5, v = 7, PC1
#lambda = 1.5 # Best for ortonii with k = 2, v = 7, PC1

# Estimate corrected p-values (method described in original LEA paper; different from the tutorial, which uses LEA functions, but same answer)
p.values = pchisq(zs.median^2/lambda, df = 1, lower = FALSE) # compute adjusted p-values from the median z-scores

# Set working directory
setwd(paste0(path, "LEA/LFMM_", sp, "_t", t, "_s", s, "_n", n, "_", envdata.name, "_v", n.var, "_k", K))

# Plot histogram of adjusted p-values
pdf(file = paste0(sp, "_k", K, "_L=", lambda, ".pdf")) # Make pdf to plot
hist(p.values, col = "red", main = paste0("Corrected p-values for Anolis ", sp), breaks = 10)
dev.off() # closing plot

## PART 5: Control of false discoveries: Correct for multiple testing using the Benjamini-Hochberg algorithm

# Find candidate SNPs based on false-discovery rate (FDR)
L = dim(gendata)[2] # L = Total number of SNPs
w = which(sort(p.values) < q * (1:L)/L) # the Benjamini-Hochberg algorithm
candidates.list = order(p.values)[w] # candidate SNPs sorted by adjusted p-value

str(candidates.list) # how many candidate SNPs did we get?

# Save our progress
save.image(file = paste0(sp, ".RData"))

## PART 5: Extract the sequences of loci that harbor the candidate SNPs flagged by LFMM
# This part is very specific to ipyrad and VCFtools outputs and the way I organized my files!

# Using loci object as a SNP map, extract those loci hosting candidate SNPs flagged by LFMM, keeping only unique entries for each locus
flagged.list = unique(loci[loci$SNP.number %in% candidates.list, ]$locus.name) # using candidate SNPs from LFMM as a mask to extract their corresponding loci
flagged.list # display flagged loci

# Read sequences of all loci (".1seq-loci" file generated using bash scripts to extract one sequence per locus from Ipyrad's loci output file)
all.loci = read.csv(paste0(path, "data/", sp, "_t", t, "_s", s, "_n", n, "_outfiles/", sp, "_t", t, "_s", s, "_n", n, ".1seq-loci"), header = FALSE)
names(all.loci) = c("locus.name", "sequence") # add column names
all.loci$locus.name %<>% gsub(pattern = "locus_", replacement = "", .) # change locus names, removing "locus_"

# Obtain sequences for each unique flagged loci and transform them in fasta format
sequences = all.loci[all.loci$locus.name %in% flagged.list, ] # using candidate SNPs from LFMM as a mask to extract their corresponding loci
if (dim(sequences)[[1]] !=0) {
    sequences %<>% mutate(., "locus.name" = paste0(">locus", locus.name)) # mutate rows to name loci following fasta format
  } else {
    print(paste("Warning: No loci associated with environmental variable")) # in case no loci is associated with predictor variable
  }

# Remove gaps from sequences for blasting
sequences$sequence %<>% gsub(pattern = "-", replacement = "", .) # removes gaps from sequences. "." denotes object when using "%<>%"

# Save list of loci while transforming into a fasta file
write.table(sequences, file = paste0(sp, "_PC1_L", lambda, "_a", q, ".fasta"), row.names = FALSE, quote = FALSE, sep = "\n", eol = " \n", col.names = FALSE)

# Extract SNPs from candidate loci
candidate.SNPs = gendata[ , colnames(gendata) %in% flagged.list]

save.image(file = paste0(sp, ".RData")) # Save our progress

# Create a directory to save matrix of candidate SNPs
dir.create(paste0(path, "data/", sp, "_t", t, "_s", s, "_n", n, "_candidate_SNPs"))

# Set this new directory as the working directory
setwd(paste0(path, "data/", sp, "_t", t, "_s", s, "_n", n, "_candidate_SNPs"))

# Save matrix of candidate SNPs in .geno format
write.table(candidate.SNPs, file = paste0(path, "data/", sp, "_t", t, "_s", s, "_n", n, "_candidate_SNPs/", sp, "_t", t, "_s", s, "_n", n, "_candidate_SNPs.012"), sep = "\t", quote = FALSE)

# Exclude candidate SNPs from the all-SNP dataset
neutral_SNPs = gendata[ , colnames(gendata)[colnames(gendata) %in% flagged.list == F]]

# Create a directory to save matrix of candidate SNPs
dir.create(paste0(path, "data/", sp, "_t", t, "_s", s, "_n", n, "_neutral_SNPs"))

# Set this new directory as the working directory
setwd(paste0(path, "data/", sp, "_t", t, "_s", s, "_n", n, "_neutral_SNPs"))

# Save matrix of candidate SNPs in .geno format
write.table(neutral_SNPs, file = paste0(path, "data/", sp, "_t", t, "_s", s, "_n", n, "_neutral_SNPs/", sp, "_t", t, "_s", s, "_n", n, "_neutral_SNPs.012"), sep = "\t", quote = FALSE)

save.image(file = paste0(sp, ".RData")) # Save our progress

# Done!