
# -----------------------------------
# Transmission Dynamics Analysis
# -----------------------------------
# Cluster Analysis: Identify genetic clusters of Plasmodium falciparum 
# using methods like network analysis to infer transmission links.
# Create Network from genetic similarity indices (IBS, IBD)

# To get coordinates of gambian cities
# https://www.longitude-latitude-maps.com/country/77,Gambia/

# Install and load StAMPP if not already installed
# if (!requireNamespace(c("StAMPP", "poppr", "dartR"), quietly = TRUE))
#    install.packages(c("StAMPP", "poppr", "dartR"))

library(tidyverse)
library(StAMPP)
library(poppr)
library(dartR)

source("scripts/functions.R")

df <- read_delim("data/barcode_gambia.tsv")

data <- df %>% 
   dplyr::filter(Location != "Nkakat Eyamba" & Location != "Ngayen Sanjal") %>% 
   arrange(Year, Location)

output_dir <- "spatial_analysis"
if(!dir.exists(output_dir)) dir.create(output_dir)

# ----------------------------
# GENETIC SIMILARITY INDICES
# ----------------------------

# ----------------------
# 1. Identity-by-state
# ----------------------

# ibs_matrix <- calculateIBS_v2(data[, -c(1:5)])  # Take alleles across SNPs positions per sample
ibs_matrix <- calculateIBS2(data[, -c(1:5)])  # Take alleles across SNPs positions per sample
rownames(ibs_matrix) = colnames(ibs_matrix) <- data %>% 
   pull(sample_internal_id)

# Save IBS matrix
write.table(ibs_matrix, file.path(output_dir, "gambia_ibs.tsv"), 
          col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')

# Estimate the proportion of similarity between the samples
# Measure the running time using microbenchmark
# result <- microbenchmark(
#    ibs_matrix1 <- compute_pairwise_IBS_v1(pf.gambian.barcodes[1:100, ]), # Take sequence barcode
#    ibs_matrix2 <- compute_pairwise_IBS_v2(pf.gambian.barcodes[1:100, ]), # Take sequence barcode
#    
#    times = 100
# )
# 
# print(result)

# Perform hierarchical clustering
hc <- hclust(as.dist(1 - ibs_matrix), method = 'complete')  # Using 1 - similarity to convert to distance

# Plot a heatmap with clustering
pheatmap(gen_sim_matrix, clustering_distance_rows = as.dist(1 - ibs_matrix),
         clustering_distance_cols = as.dist(1 - ibs_matrix),
         clustering_method = 'complete')

# ------------------------------------
# 2. Jaccard Similarity Coefficient
# ------------------------------------
jaccar_matrix <- compute_jaccard_index(data[, -c(1:5)], report_progress = TRUE) 

colnames(jaccar_matrix) = rownames(jaccar_matrix) <- data %>% pull(sample_internal_id)

# Save JACCARD matrix
write.table(jaccar_matrix, file.path(output_dir, "jaccard_matrix.tsv"), 
          col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')

# -----------------------------------------------
# 3. Sørensen-Dice Coefficient (Dice Similarity)

# https://www.ars.usda.gov/ARSUserFiles/50620500/Publications/KJL/sim_coef.pdf
# -----------------------------------------------
dice_matrix <- compute_dice_index(data[, -c(1:5)], report_progress = TRUE) 
colnames(dice_matrix) = rownames(dice_matrix) <- data %>% pull(sample_internal_id)

# Save DICE matrix
write.table(dice_matrix, file.path(output_dir, "dice_matrix.tsv"), 
            col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')


my_hclust <- hclust(dist(dice_matrix), method = "complete")

# load package
library(dendextend)

as.dendrogram(my_hclust) %>%
   plot(horiz = TRUE)

# ------------------------
# 4. FST (Fixation Index) 
# ------------------------

# We will use the StAMPP
# This function calculates pairwise Fst values along with confidence intervals and p-values
# between populations according to the method proposed by
# Wright(1949) and updated by Weir and Cockerham (1984)

# Convert alleles to lower case
test <- sapply(data[, -c(1:5)], tolower)
rownames(test) <- data %>% 
   pull(sample_internal_id)

# Get location of sample ids
pop <- data %>% pull(Location)

# Convert matrix to DNAbin object
test <- as.DNAbin(test)

# Convert DNAbin object to genind
obj.genind <- DNAbin2genind(test, pop = pop, exp.char = c("a","t","g","c"))

# Convert genind to genlight object
obj.genlight <- dartR::gi2gl(obj.genind)

# Estimate Fst between population memberships
fst_result <- StAMPP::stamppFst(obj.genlight, 100, 95, 1)
print(fst_result$Fsts)

write_rds(fst_result, file.path(output_dir, "pairwise_pop_fst.rds"))

# ============ YEARLY BASES ========
# Get location of sample ids
years <- data %>% pull(Year)

# Convert DNAbin object to genind
obj.genind <- DNAbin2genind(test, pop = years, exp.char = c("a","t","g","c"))

# Convert genind to genlight object
obj.genlight <- dartR::gi2gl(obj.genind)

# Estimate Fst between population memberships
fst_result <- StAMPP::stamppFst(obj.genlight, 100, 95, 1)
print(fst_result$Fsts)

write_rds(fst_result, file.path(output_dir, "pairwise_years_fst.rds"))

# ------------------------------------
# 5. Nei's Genetic Distance (Nei 1972)
# ------------------------------------

# Assuming 'genind_obj' is a genind object from the adegenet package containing your SNP data
nei_dist <- poppr::nei.dist(obj.genind)
nei_dist <- as.matrix(nei_dist)
print(nei_dist)


mmod::Gst_Nei(obj.genind)
mmod::D_Jost(obj.genind)
mmod::diff_stats(obj.genind)
mmod::pairwise_D(obj.genind)

# ------------------------------------------
# 6. AMOVA (Analysis of Molecular Variance)
# ------------------------------------------
# https://grunwaldlab.github.io/Population_Genetics_in_R/AMOVA.html
# https://parasitesandvectors.biomedcentral.com/articles/10.1186/s13071-017-2013-z
# https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2019.00184/full

# Assuming 'genind_obj' as before and 'strata' to define hierarchical levels
my_strata <- data %>% dplyr::select(Location, Year)
names(my_strata) <- c("Pop", "subPop")

my_strata <- my_strata %>% 
   group_by(subPop,Pop) %>% 
   mutate(label = group_indices(., subPop,Pop)) %>% 
   ungroup() %>% dplyr::select(-2) %>% 
   dplyr::rename(subPop = label)

my_strata <- my_strata %>% ungroup() %>% dplyr::select(-2) %>% 
   dplyr::rename(subPop = label) %>% as.data.frame()

strata(obj.genind) <- my_strata

amova_results <- poppr::poppr.amova(obj.genind, hier = ~Pop/subPop)
print(amova_results)

amova.test <- randtest(amova.result) # Test for significance
plot(amova.test)
amova.test

# ---------------------------
# 7. Shared Allele Distance
# ---------------------------
# https://dyerlab.github.io/applied_population_genetics/genetic-distances.html

# Assuming 'gen_data' is a matrix with individuals as rows and loci as columns
calculate_shared_allele_distance <- function(gen_data) {
   num_individuals <- nrow(gen_data)
   distance_matrix <- matrix(NA, nrow = num_individuals, ncol = num_individuals)
   diag(distance_matrix) <- 0
   
   for (i in 1:(num_individuals-1)) {
      for (j in (i+1):num_individuals) {
         shared_alleles <- sum(gen_data[i, ] == gen_data[j, ], na.rm = TRUE)
         total_alleles <- sum(!(gen_data[i, ] %in% c("N", "X")) & !(gen_data[j, ] %in% c("N", "X")))
         distance_matrix[i, j] <- 1 - (shared_alleles / total_alleles)
         distance_matrix[j, i] <- distance_matrix[i, j]  # Symmetric
      }
   }
   
   return(distance_matrix)
}

shared_allele_distance <- calculate_shared_allele_distance(data[, -c(1:5)])
print(shared_allele_distance)







