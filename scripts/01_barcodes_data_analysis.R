
library(tidyverse)
library(dlookr)
library(flextable)
library(microbenchmark) # Measuring running time

source("scripts/functions.R")

gambian.barcodes <- readxl::read_excel("data/FinalgambianDataset.xlsx", na = c("", "-")) %>% 
   rename(sample_internal_id = `Sample Internal ID`, sample_external_id = `Sample External ID`) %>%
   select(-c("sample_external_id", "Country", "Study", "Chloroquine", "Sulfadoxine", "Mefloquine", "Pyrimethamine"))

# # Explore the data
# gambian.barcodes %>% diagnose() %>% flextable()
# gambian.barcodes %>% diagnose_outlier() %>% flextable()
# 
# plot_outlier(gambian.barcodes)
# 
# # Generate EDA report
# gambian.barcodes %>% diagnose_web_report(output_format = "html")
# 
# # Number of monoclonal samples per year
# counts <- gambian.barcodes %>% 
#     filter(Species == "Pf" & McCOIL == 1) %>% # !is.na(McCOIL)
#     group_by(Year, McCOIL) %>% 
#      count()

# # Extract Plasmodium falciparum data
# pf.gambian.barcodes <- gambian.barcodes %>% 
#     filter(Species == "Pf" & !is.na(McCOIL)) %>% 
#     select(c(1, 3,4,7, "McCOIL", "Barcode"))
# 
# # Estimate the proportion of missing genotypes
# # for each sample
# missing <- prop.missing(pf.gambian.barcodes)

# Another method for missing data
# on samples and loci
data <- gambian.barcodes %>% 
   filter(Species == "Pf" & !is.na(McCOIL)) %>% 
   select(-c(4, 7:50))

# Extract only loci from the data
df <- data[, -c(1:5)]

# Estimate missing genotypes on samples
missing <- calculate_missing_genotypes(df, 1) # Use 1 if you want missing data on samples

# Samples with less than 20% missing data
index <- which(missing < 1)

# Remove samples with high proportion of missing data
data <- data[index,]

# Extract only loci from the data
df <- data[, -c(1:5)]

# Estimate missing genotypes on loci
missing <- calculate_missing_genotypes(df, 2)

# Loci with less than 20% missing data
index <- which(missing < 20)

# Remove samples with high proportion of missing data
df <- df[, index]

data <- bind_cols(data[, c(1:5)], df)

write.table(data, "data/barcode_gambia.tsv", col.names = TRUE, 
            row.names = FALSE, quote = FALSE, sep = "\t")
