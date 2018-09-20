### This script's goal is to filter the SNP data by MAP using VCFTools and get them ready for sNMF analyses.
### By Ivan Prates with the help of Kevin P. Mulder, June 2018.
### Smithsonian National Museum of Natural History, Washington DC, USA.

# Filtering SNPs by MAF and then extracting one SNP per locus for SNMF analyses
# From Ipyrad vcf format output file

# opening a loop to repeat over species:
for s in ortonii_t70_s10_n23 punctatus_t70_s10_n46 # list target taxa and assembly settings
do

# create a folder and change directory to that new folder
mkdir /Users/dequeirk/"Dropbox (Smithsonian)"/ivan_lab/2018_Anolis_GEA/2018-03/data/VCFtools_SNMF_${s} 
cd /Users/dequeirk/"Dropbox (Smithsonian)"/ivan_lab/2018_Anolis_GEA/2018-03/data/VCFtools_SNMF_${s}

# Filter by MAF = 0.10 using VCF Tools
vcftools --vcf /Users/dequeirk/"Dropbox (Smithsonian)"/ivan_lab/2018_Anolis_GEA/2018-03/data/${s}_outfiles/${s}.vcf --recode --out ${s}_filtered --min-alleles 2 --max-alleles 2 --maf 0.10

# Now extracting one random SNP per locus
cat ${s}_filtered.recode.vcf | grep "#" > ${s}_1SNP-locus.vcf # extract the headers (whose line starts with #), save in different file
cat ${s}_filtered.recode.vcf | grep -v "#" | sort -u -k1,1 >> ${s}_1SNP-locus.vcf # sort by 'chromosome' name and extract unique, concatenate to that file created

# Now saving in 012 format
vcftools --vcf ${s}_1SNP-locus.vcf --out ${s}_1SNP-locus --012

# Now replace -1 with 9 to indicate missing sites
cat ${s}_1SNP-locus.012 | sed -e 's/-1/9/g' ${s}_1SNP-locus.012 > ${s}_1SNP-locus_tmp.012 # replace and save in temporary file
mv ${s}_1SNP-locus_tmp.012 ${s}_1SNP-locus.012 # rename temp file as final 012 file

# close the loop:
done