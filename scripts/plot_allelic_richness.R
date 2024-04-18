

# https://cran.r-project.org/web/packages/tmap/vignettes/tmap-getstarted.html

# https://journals.plos.org/plosone/article/figures?id=10.1371/journal.pone.0029845

# https://www.researchgate.net/publication/357581462_Genetic_association_of_TMPRSS2_rs2070788_polymorphism_with_COVID-19_case_fatality_rate_among_Indian_populations/figures?lo=1

# https://www.researchgate.net/publication/5656026_Polymorphisms_of_TNF-enhancer_and_gene_for_FcgammaRIIa_correlate_with_the_severity_of_falciparum_malaria_in_the_ethnically_diverse_Indian_population/figures?lo=1

# https://stackoverflow.com/questions/67001602/can-you-create-multiple-polygons-in-r-from-a-dataframe-containing-the-vertices


library(tidyverse)
library(sf)
library(tmap)

# Only run this chunk if you are loading a local shapefile
shape.gmb <- st_read( "gis_coordinates/geoBoundaries-GMB-ADM3-all/geoBoundaries-GMB-ADM3_simplified.shp")

# Load coordinates
coordinates <- readxl::read_xlsx("coordinates.xlsx", sheet = 2)

# Load Allelic richness data
al.richness.data <- readRDS("genetic_diversity/allelic_richness.rds")

ar.matrix <- al.richness.data$allel.rich[[1]]$all.richness %>% as.matrix()

image.plot(ar.matrix)
contour(ar.matrix, add = TRUE)

# Extract mean population allelic richness
mean.ar <- al.richness.data$allel.rich[[1]]$mean.richness

df <- data.frame(location = names(mean.ar), allelic_richness = mean.ar) %>% 
   left_join(., coordinates, by = c("location" = "Locations")) %>% 
   drop_na()

df_sf <- st_as_sf(df[-nrow(df),], coords = c("longitude", "latitude"), crs = st_crs(shape.gmb), na.fail = FALSE, agr = "constant")

# Plotting
ggplot() +
   geom_sf(data = shape.gmb, fill = "white", color="#000000") + # Plot the base map
   geom_sf(data = df_sf, aes(color = allelic_richness), size = 5) + # Plot allelic richness
   scale_color_gradient(low = "yellow", high = "brown") + # Color gradient for allelic richness
   labs(title = "Allelic Richness Distribution in the Gambia", 
        color = "Allelic Richness") +
   ggrepel::geom_label_repel(data = df, aes(label = location, x = longitude, y = latitude), 
                             color = 'black', size = 3, box.padding = unit(1, "lines"), 
                             segment.color = '#132B43',  max.overlaps = Inf) +
   theme_void() + 
   theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
         legend.title = element_text(face = "bold", size = 12))

ggsave("genetic_diversity/allelic_richness.pdf",
       width = 190, height = 150, dpi = 600, units = "mm")

joined_sf <- st_join(shape.gmb, df_sf, join = st_nearest_feature)

ggplot(joined_sf) +
   geom_sf(aes(fill = allelic_richness)) +
   scale_fill_gradient(low = "yellow", high = "brown") +
   labs(title = "Allelic Richness Distribution in the Gambia", color = "Allelic Richness") +
   theme_void() + 
   theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
         legend.title = element_text(face = "bold", size = 12))

ggsave("genetic_diversity/allelic_richness_v1.pdf",
       width = 190, height = 150, dpi = 600, units = "mm")




