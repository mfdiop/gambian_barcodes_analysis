

library(patchwork)

prevalence <- readRDS("temporal_analysis/prevalence_over_time.rds")

plot_common_haplotypes <- function(hap_data, gene_name){
   
   hap_data <- hap_data %>%
      group_by(Haplotypes) %>%
      filter(n_distinct(Year) >= 3) %>%
      ungroup() 
   
   # Select the top 2 most common haplotypes based on average frequency
   top_haplotypes <- hap_data %>%
      group_by(Year, Haplotypes) %>%
      summarise(AvgFrequency = mean(Prevalence), .groups = 'drop') %>%
      ungroup() %>%
      filter(!grepl("-", Haplotypes) & AvgFrequency > 1)
   
   # Filter the original dataframe to keep only the top haplotypes
   # Plotting
   p <- ggplot(top_haplotypes, aes(x = Year, y = AvgFrequency, color = Haplotypes, group = Haplotypes)) +
      geom_line(linewidth = 1) + 
      geom_point(size = 2) +
      theme_minimal() +
      labs(x = "Year", y = "Frequency (%)", color = gene_name) +
      scale_color_brewer(palette = "Set1") +
      scale_x_continuous(breaks = 2016:2022, limits = c(2016, 2022)) + # Customize x-axis
      theme(axis.line = element_line(colour = "#000000", linewidth = 2),
            axis.title = element_text(color = "#000000", face = "bold", size = 11),
            axis.text.x = element_text(color = "#000000", face = "bold", angle = 90),
            axis.text.y = element_text(color = "#000000", face = "bold", size = 9),
            legend.title = element_text(color = "#000000", face = "bold", size = 11, hjust = .5))
   
   return(p)
}


# Apply the function to each element of the list, passing the element name (gene name) as an argument
plots <- mapply(plot_common_haplotypes, prevalence, names(prevalence), SIMPLIFY = FALSE)

# Display the plots
# invisible(lapply(plots, print))

combined_plot <- reduce(plots, `+`) + 
   plot_layout(ncol = 3) # Arrange plots in a single column

# Print the combined plot
combined_plot

ggsave("temporal_analysis/haplotypes_over_time.pdf", 
       width = 390, height = 320, units = "mm", dpi = 600)

#  -----------------------------------------------------
prevalence <- readRDS("temporal_analysis/overall_prevalence.rds")

plot_common_haplotypes <- function(hap_data, gene_name){
   
   # Select the top 2 most common haplotypes based on average frequency
   top_haplotypes <- hap_data %>%
      group_by(Haplotypes) %>%
      summarise(AvgFrequency = mean(Prevalence), .groups = 'drop') %>%
      ungroup() %>%
      filter(!grepl("-", Haplotypes) & AvgFrequency > 1)
   
   # Filter the original dataframe to keep only the top haplotypes
   # Plotting
   p <- ggplot(top_haplotypes, aes(x = Year, y = AvgFrequency, color = Haplotypes, group = Haplotypes)) +
      geom_line(linewidth = 1) + 
      geom_point(size = 2) +
      theme_minimal() +
      labs(x = "Year", y = "Frequency (%)", color = gene_name) +
      scale_color_brewer(palette = "Set1") +
      scale_x_continuous(breaks = 2016:2022, limits = c(2016, 2022)) + # Customize x-axis
      theme(axis.line = element_line(colour = "#000000", linewidth = 2),
            axis.title = element_text(color = "#000000", face = "bold", size = 11),
            axis.text.x = element_text(color = "#000000", face = "bold", angle = 90),
            axis.text.y = element_text(color = "#000000", face = "bold", size = 9),
            legend.title = element_text(color = "#000000", face = "bold", size = 11, hjust = .5))
   
   return(p)
}


# Apply the function to each element of the list, passing the element name (gene name) as an argument
plots <- mapply(plot_common_haplotypes, prevalence, names(prevalence), SIMPLIFY = FALSE)

# Display the plots
# invisible(lapply(plots, print))

combined_plot <- reduce(plots, `+`) + 
   plot_layout(ncol = 3) # Arrange plots in a single column

# Print the combined plot
combined_plot

ggsave("temporal_analysis/haplotypes_over_time.pdf", 
       width = 390, height = 320, units = "mm", dpi = 600)