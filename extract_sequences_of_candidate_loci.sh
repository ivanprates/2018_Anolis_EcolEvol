## Code by Ivan Prates with the help of Kevin Mulder
## February 2018
## Smithsonian National Museum of Natural History, Washington DC, USA

# Extracting one sequence per locus from .loci file

# opening a loop to repeat over species:
#for s in ortonii_t70_s10_n23 # list target taxa and number of sampled individuals
for s in punctatus_t70_s10_n47 # list target taxa and number of sampled individuals
do

# create a folder and change directory to that new folder
mkdir /Users/dequeirk/"Dropbox (Smithsonian)"/ivan_lab/2018_Anolis_GEA/2018-03/data/VCFtools_LFMM_${s} 
cd /Users/dequeirk/"Dropbox (Smithsonian)"/ivan_lab/2018_Anolis_GEA/2018-03/data/VCFtools_LFMM_${s}

# Filter by MAF = 0.10 using VCF Tools
vcftools --vcf /Users/dequeirk/"Dropbox (Smithsonian)"/ivan_lab/2018_Anolis_GEA/2018-03/data/${s}_outfiles/${s}.vcf --recode --out ${s}_MAF10_LFMM --min-alleles 2 --max-alleles 2 --maf 0.10

# Now saving in 012 format
vcftools --vcf ${s}_MAF10_LFMM.recode.vcf --out ${s}_MAF10_LFMM --012

# Now replace -1 with 9 to indicate missing sites
cat ${s}_MAF10_LFMM.012 | sed -e 's/-1/9/g' ${s}_MAF10_LFMM.012 > ${s}_MAF10_LFMM_tmp.012 # replace and save in temporary file
mv ${s}_MAF10_LFMM_tmp.012 ${s}_MAF10_LFMM.012 # rename temp file as final 012 file

# Extract one sequence from each locus from alignments, save in a different file
cat /Users/dequeirk/"Dropbox (Smithsonian)"/ivan_lab/2018_Anolis_GEA/2018-03/data/${s}_outfiles/${s}.loci | # cat: Do everything on this file
grep -B 1 "/" | # select a given line and the one above
sed 's/|/locus_/' | # replace | by locus_
sed 's/^\/.*loc/loc/' |  # replace all the junk from beginning of line (^) until "loc" with nothing
gsed 's/punc.*\s//' | # replace everything from "punc" to first space with nothing
paste - - - | #  
sed 's/|//' | # get rid of second |
awk ' { t = $1; $1 = $2; $2 = t; print; } ' | # swap first and second column
cut -f 1,2 -d ' ' | 
tr ' ' ',' > /Users/dequeirk/"Dropbox (Smithsonian)"/ivan_lab/2018_Anolis_GEA/2018-03/data/${s}_outfiles/${s}.1seq-loci # extract first and second column, replace space by comma, and save it

# Same thing, one-liner:
# cat /Users/dequeirk/"Dropbox (Smithsonian)"/ivan_lab/2018_Anolis_GEA/2018-03/data/punctatus_t90_s10_n46_outfiles/punctatus_t90_s10_n46.loci | grep -B 1 "/" | sed 's/|/locus_/' | sed 's/^\/.*loc/loc/' | gsed 's/punc.*\s//' | paste - - - | sed 's/|//' | awk ' { t = $1; $1 = $2; $2 = t; print; } ' | cut -f 1,2 -d ' ' | tr ' ' ',' > /Users/dequeirk/"Dropbox (Smithsonian)"/ivan_lab/2018_Anolis_GEA/2018-03/data/punctatus_t90_s10_n46_outfiles/punctatus_t90_s10_n46.1seq-loci

# Without saving:
# cat /Users/dequeirk/"Dropbox (Smithsonian)"/ivan_lab/2018_Anolis_GEA/2018-03/data/punctatus_t90_s10_n46_outfiles/punctatus_t90_s10_n46.loci | grep -B 1 "/" | sed 's/|/locus_/' | sed 's/^\/.*loc/loc/' | gsed 's/punc.*\s//' | paste - - - | sed 's/|//' | awk ' { t = $1; $1 = $2; $2 = t; print; } ' | cut -f 1,2 -d ' ' | tr ' ' ','

# close the loop:
done

