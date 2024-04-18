
# https://rpubs.com/kkeenan02/divMigrate
# https://onlinelibrary.wiley.com/doi/10.1002/ece3.2096
# https://www.rdocumentation.org/packages/diveRsity/versions/1.9.89

# --------------------
# Temporal Dynamics: 
# --------------------
# Analyse how genetic variation changes over time 
# to understand the spread and dynamics of transmission.

library(tidyverse)

raw_barcode <- readxl::read_xlsx("data/FinalgambianDataset.xlsx") %>% 
   dplyr::rename(sample_internal_id = `Sample Internal ID`)

df <- read_delim("data/barcode_gambia.tsv")

data <- df %>% dplyr::filter(Location != "Nkakat Eyamba" & Location != "Ngayen Sanjal")

output_dir <- "temporal_analysis"
if(!dir.exists(output_dir)) dir.create(output_dir)

# 1. Temporal FST (Fixation Index)
# Temporal FST compares genetic differentiation between populations at different time points. 
# It measures how genetic variance is distributed over time, providing insights into the population 
# structure's temporal dynamics. A high FST value indicates significant genetic differentiation, 
# which could result from selection pressure, population bottlenecks, or founder effects.

time_points <- data %>% 
   pull(Year) %>% 
   unique() %>% 
   sort()

Fst <- tibble(rows = NULL)

for (t in time_points) {
   
   subset <- data %>% 
      dplyr::filter(Year == t)
   
   # Get number of populations
   Nbpop <- subset %>% pull(Location) %>% unique %>% length()
   
   if(Nbpop > 1){
      # Extract population IDs
      pop <- subset %>% pull(Location)
      
      # Convert alleles to lower case
      alleles <- sapply(subset[, -c(1:5)], tolower)
      
      # Rename matrix rownames
      rownames(alleles) <- subset %>% 
         pull(sample_internal_id)
      
      # Convert matrix to DNAbin object
      alleles <- ape::as.DNAbin(alleles)
      
      # Convert DNAbin object to genind
      obj.genind <- adegenet::DNAbin2genind(alleles, pop = pop, exp.char = c("a","t","g","c"))
      
      # Convert genind to genlight object
      obj.genlight <- dartR::gi2gl(obj.genind)
      
      # Estimate Fst between population memberships
      fst_result <- StAMPP::stamppFst(obj.genlight, 100, 95, 1)
      
      fst.df <- dplyr::as_tibble(fst_result$Fsts, rownames = "p1")
      
      Fst <- rbind(Fst, fst.df %>% pivot_longer(!p1, names_to = "p2", values_to = "fst") %>% add_column(Year = t))
   }
}

Fst <- Fst %>% 
   dplyr::filter(!is.na(fst)) %>% 
   dplyr::mutate(fst = ifelse(fst<0, 0, fst),
                 fst = round(fst, 3))

# Save Fst results
write.table(Fst, file.path(output_dir, "Fst_over_time.tsv"), col.names = TRUE, 
            row.names = FALSE, quote = FALSE, sep = '\t')

# install diveRsity from github
# devtools::install_github("kkeenan02/diveRsity")

library(diveRsity)

# load Test_data
data(Test_data)

div_res <- divMigrate(infile = Test_data, stat = "d_jost")

# ---------------------------------
# 2. Haplotype Diversity Over Time
# ---------------------------------
# Haplotype diversity measures the uniqueness of the genetic variants present in a population. 
# Tracking changes in haplotype diversity over time can indicate the emergence or disappearance of strains, 
# reflecting how genetic diversity is influenced by transmission dynamics, selection pressures, 
# and population size changes.

hap.matrix <- tibble(rows = NULL)

for (t in time_points) {
   
   subset <- data %>% 
      dplyr::filter(Year == t)
   
   # Get number of populations
   Nbpop <- subset %>% pull(Location) %>% unique()
   
   for (p in 1:length(Nbpop)) {
      
      pop.data <- subset %>% 
         dplyr::filter(Location == Nbpop[p])
      
      # Extract population IDs
      pop <- pop.data %>% pull(Location)
      
      # Number of samples
      nbsapls <- nrow(pop.data)
      
      # Convert alleles to lower case
      alleles <- sapply(pop.data[, -c(1:5)], tolower)
      
      # Rename matrix rownames
      rownames(alleles) <- pop.data %>% 
         pull(sample_internal_id)
      
      # Convert matrix to DNAbin object
      alleles <- ape::as.DNAbin(alleles)
      
      # Assuming 'haplotypes' is a DNA sequence object in R
      hap_div <- round(pegas::hap.div(alleles), 3)
      
      # print(hap_div)
      
      hap.matrix <- rbind(hap.matrix, tibble(Year = t, Location = Nbpop[p], nb_samples = nbsapls, hap.div = hap_div))
   }
   
}

# Save Haplotype diversity results
write.table(hap.matrix, file.path(output_dir, "hap_diversity_over_time.tsv"), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')

# -------------------
# 3. Skyline Plots
# -------------------
# Skyline plots visualize changes in effective population size over time, based on genetic data. 
# They can infer demographic history and identify periods of population expansion or contraction, 
# which may correlate with known transmission events or disease outbreaks.
# https://taming-the-beast.org/tutorials/Skyline-plots/
   

# ---------------------------
# 4. Phylodynamic Analysis
# ---------------------------
# Phylodynamics integrates phylogenetic methods with epidemiological models to study how infectious diseases spread 
# and evolve. Using molecular sequences, phylodynamic analysis can estimate the effective reproduction number (R), 
# selection pressures, and the rate of geographic spread over time. 
# Tools like BEAST (Bayesian Evolutionary Analysis Sampling Trees) can perform these analyses, 
# generating estimates of evolutionary rates and population dynamics.
# 
# -----------------------------
# 5. Molecular Clock Analysis
# -----------------------------
# Molecular clock analysis estimates the rate of mutation over time, providing insights into the evolutionary timeline 
# of pathogens. By calibrating a molecular clock, researchers can estimate the time to the most recent common ancestor 
# (TMRCA) for a set of sequences, helping to reconstruct the temporal aspects of transmission chains.

# ------------------------------
# 4. Drug Resistance Analysis
# ------------------------------
# Association with Drug Resistance: Since the SNPs are linked to drug resistance, 
# analyse the prevalence of these resistance-associated SNPs 
# in different locations and over time.

# https://wellcomeopenresearch.org/articles/7-45

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7192906/

drug_resistance <- data %>%
   dplyr::select(1) %>% 
   dplyr::left_join(., raw_barcode, by = "sample_internal_id") %>% 
   dplyr::select("sample_internal_id", "Year", "Location", "McCOIL", 14:21) %>% 
   dplyr::rename(P23.BP = `P23:BP`)

resistance_genes <- names(drug_resistance)[-c(1:4)]

# Estimate the prevalence of mutations for each gene
prevalence_results <- lapply(resistance_genes, function(gene) {
   gene_data <- drug_resistance[[gene]]
   extract_and_count_mutations(gene_data, gene)
})

# Names the list elements with the gene names for easy identification
names(prevalence_results) <- resistance_genes

# Printing results for demonstration (you might want to format or export these results)
print(prevalence_results)

write_rds(prevalence_results, file = file.path(output_dir, "overall_prevalence.rds"))

# --------------------------
# Estimate the prevalence of 
# mutations over time
# --------------------------
prevalence_results_over_time <- list()

for (gene in resistance_genes) {
   prevalence_results_over_time[[gene]] <- expand_and_extract_mutations_time(drug_resistance, gene)
}

# Save list
write_rds(prevalence_results_over_time, 
          file = file.path(output_dir, "prevalence_over_time.rds"))

# for (d in seq_along(drug_resistance_genes)) {
#    
#    # Number of samples
#    nsmpl <- nrow(drug_resistance)
#    
#    # Extract drug-resistance haplotypes
#    haplotypes <- drug_resistance %>% 
#       dplyr::pull(d) %>% 
#       str_split(., " ") %>% 
#       unlist() %>% 
#       unique()
#    
#    haplotypes <- gsub("\\[|\\]", " ", haplotypes)
# 
#    index <- which(nchar(haplotypes)>5)
#    
#    if(is_empty(index)){
#       
#       prevalence <- drug_resistance %>%
#          # dplyr::select(all_of(gene)) %>% 
#          separate_longer_delim(d, delim = " ") %>% 
#          dplyr::group_by(!!.[[d]]) %>% 
#          dplyr::summarise(N = n(), Prevalence = round((N/nsmpl)*100, 2)) %>% 
#          dplyr::rename(hap = `<chr>`) %>% 
#          dplyr::mutate(genes = d, hap = recode(hap, "-" = "None")) %>% 
#          ungroup() %>% relocate(genes, .before = hap)
#    }
#    else{
#    }
#  
# }
# 

# -----------------------------------
# 5. Super-Spreader Identification
# -----------------------------------
# Transmission Network Analysis: Construct and analyse transmission networks to identify 
# potential super-spreaders based on the number and distribution of linked cases.
library(igraph)

# Assuming 'genetic_similarity_matrix' is a matrix of genetic similarity indices between samples
g <- graph.adjacency(genetic_similarity_matrix, mode = "undirected", weighted = TRUE)
V(g)$name <- sample_ids  # Assigning sample IDs to nodes

# Calculate degree centrality
degree_cent <- degree(g, mode = "all")

# Identify potential super-spreaders as nodes with degree centrality above a certain threshold
threshold <- mean(degree_cent) + 2*sd(degree_cent)  # Example threshold
super_spreader_ids <- names(degree_cent[degree_cent > threshold])

print(super_spreader_ids)


