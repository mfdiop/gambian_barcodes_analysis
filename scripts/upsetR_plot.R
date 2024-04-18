
# https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html
# https://krassowski.github.io/complex-upset/articles/Examples_R.html


# if (!requireNamespace("UpSetR", quietly = TRUE)) install.packages("UpSetR")
library(UpSetR)
library(tidyverse)

color_bar <- c('#00ff00', '#33FFFF', 'brown', "#0099FF", "magenta", "black", 
               "orange", "#1f9463", "#ff0000", "#000099", "yellow", "skyblue")
                        
# Function
generate_combinations <- function(haplotype) {
   # Identify segments split by brackets
   segments <- strsplit(haplotype, "\\[|\\]", perl = TRUE)[[1]]
   
   # Initialize a list to store options (constant or variable)
   options_list <- list()
   
   for (i in 1:length(segments)) {
      segment <- segments[i]
      if (i %% 2 == 0) {  # Even-indexed segments are variable parts
         options_list[[length(options_list) + 1]] <- strsplit(segment, "/")[[1]]
      } else {  # Odd-indexed segments are constant parts
         # Split constant segments into individual characters to maintain the structure
         options_list[[length(options_list) + 1]] <- segment
      }
   }
   
   # Flatten the list to prepare for combination
   options_list <- lapply(options_list, function(x) if (length(x) == 0) NA else x)
   
   # Generate all possible combinations
   combinations <- expand.grid(options_list, stringsAsFactors = FALSE)
   
   # Concatenate parts to form haplotypes
   combined_haplotypes <- apply(combinations, 1, function(x) paste0(x, collapse = ""))
   
   return(combined_haplotypes)
   # # Filter out incomplete haplotypes if any
   # valid_haplotypes <- combined_haplotypes[nchar(combined_haplotypes) == 5]
   # 
   # return(valid_haplotypes)
}

# End of function

# Example dataframe
df <- readxl::read_xlsx("data/FinalgambianDataset.xlsx") %>% 
   dplyr::rename(SampleID = `Sample Internal ID`, P23.BP = `P23:BP`) %>% 
   select(1,3,4,14:21) %>% 
   dplyr::filter(Location != "Nkakat Eyamba" & Location != "Ngayen Sanjal")

haplotype_upsetR_plot <- function(hap_data, gene_name, n){
   
   # Extract gene
   gene <- hap_data %>% 
      group_by({{gene_name}}, Year) %>%
      filter(n() >= 6) %>%
      # count(Year) %>%
      # arrange(Year)
      ungroup() %>% 
      mutate(gene_name = if_else(grepl("-----", {{gene_name}}), "None", {{gene_name}})) %>%
      filter(!grepl("-", {{gene_name}})) %>% 
      dplyr::select(SampleID, gene_name)
   
   # Extract drug-resistance haplotype
   
   my_matrix <- gene %>% 
      # Applying the function and expanding the dataframe
      rowwise() %>%
      mutate(Combinations = list(generate_combinations(gene_name))) %>%
      unnest(Combinations) %>%
      select(-gene_name) %>%
      rename(gene_name = Combinations) %>% 
      pivot_longer(-SampleID, names_to = "mutation", values_to = "Haplotype") %>% 
      mutate(Presence = 1) %>%
      pivot_wider(names_from = Haplotype, values_from = Presence, values_fill = list(Presence = 0)) %>%
      select(-mutation) %>%
      distinct()
   
   my_matrix <- my_matrix[!duplicated(my_matrix$SampleID), ]
   
   # Prepare the data for UpSetR (remove the sample identifier column)
   data_for_upsetr <- my_matrix[,-1]
   
   UpSetR::upset(as.data.frame(data_for_upsetr), 
         sets = colnames(data_for_upsetr),
         keep.order = TRUE, 
         mainbar.y.label = "Haplotype Intersections", sets.x.label = "% of haplotypes",
         main.bar.color = rev(rainbow(n)),
         sets.bar.color = "#000000")
   
}

pdf(file = "temporal_analysis/pfcrt_upset.pdf") # or other device
haplotype_upsetR_plot(df, PfCRT, 3)
dev.off()


pdf(file = "temporal_analysis/dhfr_upset.pdf")
haplotype_upsetR_plot(df, PfDHFR, 17)
dev.off()

pdf(file = "temporal_analysis/dhfr_upset.pdf")
haplotype_upsetR_plot(df, PfDHFR)
dev.off()

pdf(file = "temporal_analysis/dhps_upset.pdf")
haplotype_upsetR_plot(df, PfDHPS)
dev.off()

pdf(file = "temporal_analysis/pfmdr1_upset.pdf")
haplotype_upsetR_plot(df, PfMDR1)
dev.off()

haplotypes <- colnames(data_for_upsetr)

ComplexUpset::upset(as.data.frame(data_for_upsetr), haplotypes, min_size=5, 
                    themes=list(default = theme_minimal() +
                                   theme(legend.position = "none",
                                         axis.title.x = element_blank(),
                                         axis.title.y = element_text(face = 'bold', color = "#000000", size = 12),
                                         axis.text.x = element_blank(),
                                         axis.text.y = element_text(face = 'bold', color = "#000000", size = 10),
                                         axis.ticks.x = element_blank())),
                    name='Haplotypes', width_ratio=0.2,
                    base_annotations = list(
                       'Intersection size'=(
                          intersection_size(
                             text_mapping=aes(
                                label=paste0(round(!!get_size_mode('exclusive_intersection')/nrow(data_for_upsetr) * 100, 1), '%'),
                                size = 10), 
                             text_colors=c(on_background='brown', on_bar='yellow')
                          ) 
                          + ylab('Intersection %')
                          + scale_y_continuous(
                             labels=scales::percent_format(scale=100 / nrow(data_for_upsetr)))),
                       # ,breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100) / 100 * nrow(data_for_upsetr)
                       'Intersection ratio'=intersection_ratio()),
                    set_sizes=(upset_set_size(filter_intersections=FALSE) + 
                                  theme(axis.title = element_blank()))
)
