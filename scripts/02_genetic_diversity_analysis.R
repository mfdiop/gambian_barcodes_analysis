
# https://www.researchgate.net/publication/248565353_Estimating_allelic_richness_and_its_diversity
# https://introtogenomics.readthedocs.io/en/latest/2021_diversityTutorial.html
# https://bookdown.org/hhwagner1/LandGenCourse_book/WE_3.html

# Phylogenetic analysis
# 
# http://www.phytools.org/eqg/Exercise_3.2/
# https://pavelmatos.files.wordpress.com/2019/10/comparativephylogenetics_r_tutorial-1.pdf
# http://www.phytools.org/Cordoba2017/ex/2/Intro-to-phylogenies.html
# https://lukejharmon.github.io/ilhabela/instruction/2015/07/02/introduction-phylogenies-in-R/
# https://justinbagley.rbind.io/2011/11/22/r-functions-for-working-with-phylogenetic-trees-in-packages-ape-geiger-and-caper-part-i/
# https://rpubs.com/abanban93/558901

library(tidyverse)
require(PopGenReport)

source("scripts/functions.R")

df <- read_delim("data/barcode_gambia.tsv")

output_dir <- "genetic_diversity"
if(!dir.exists(output_dir)) dir.create(output_dir)

data <- df %>% dplyr::filter(Location != "Nkakat Eyamba" & Location != "Ngayen Sanjal")

# Convert alleles to lower case
allele <- sapply(data[, -c(1:5)], tolower)
rownames(allele) <- data %>% pull(sample_internal_id)

# Convert matrix to DNAbin object
allele <- ape::as.DNAbin(allele, )

# Get location of sample ids
pop <- data %>% pull(Location)

# Convert DNAbin object to genind
obj.genind <- adegenet::DNAbin2genind(allele, pop=pop, exp.char=c("a","t","g","c"))

# Compute allelic richness across loci
locus.ar <- popgenreport(obj.genind, mk.allel.rich=TRUE, mk.pdf=FALSE) # allel.rich(obj.genind)

# Save Allelic Richness results
write_rds(locus.ar, file.path(output_dir, "allelic_richness.rds"))

ar.matrix <- locus.ar$allel.rich[[1]]$all.richness

# Example: Assuming 'barcode_data' is your data frame
snp_frequencies <- calculate_snp_frequencies(data[,-c(1:5)])

# Save Allele frequencies results
write_rds(snp_frequencies, file.path(output_dir, "allele_frequencies.rds"))
