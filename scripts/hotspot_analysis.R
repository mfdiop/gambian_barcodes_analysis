
# https://cran.r-project.org/web/packages/sfhotspot/vignettes/introduction.html

# https://rpubs.com/quarcs-lab/spatial-autocorrelation

# https://rpubs.com/heatherleeleary/hotspot_getisOrd_tut
# https://michaelminn.net/tutorials/r-point-analysis/
# https://www.youtube.com/watch?v=OnMNZwJywjs
# https://www.youtube.com/watch?v=6ONCKKRyVgI

# Install and load necessary packages
if (!requireNamespace("spdep", quietly = TRUE)) install.packages("spdep")

library(sf)
library(sfdep)
library(spdep)
library(tidyverse)
library(patchwork)
library(pheatmap)

metadata <- read_delim("data/barcode_gambia.tsv") %>% 
  dplyr::filter(Location != "Nkakat Eyamba" & Location != "Ngayen Sanjal") %>% 
  select(sample_internal_id,  Year, Location)

# data munging
mtdt_x <- metadata %>%
  dplyr::rename(p1 = sample_internal_id)

mtdt_y <- metadata %>%
  dplyr::rename(p2 = sample_internal_id)

# Only run this chunk if you are loading a local shapefile
shape.gmb <- st_read( "gis_coordinates/geoBoundaries-GMB-ADM3-all/geoBoundaries-GMB-ADM3_simplified.shp")

# Load coordinates
coordinates <- readxl::read_xlsx("coordinates.xlsx", sheet = 2)
# x <- shape.gmb %>% left_join(., coordinates, by = c("shapeName" = "Locations"))

ibs_matrix <- read_delim("spatial_analysis/gambia_ibs.tsv")

# Perform hierarchical clustering
# hc <- hclust(as.dist(1 - ibs_matrix), method = 'complete')  # Using 1 - similarity to convert to distance
# 
# # Plot a heatmap with clustering
# pheatmap::pheatmap(ibs_matrix, 
#                    clustering_distance_rows = as.dist(1 - ibs_matrix),
#                    clustering_distance_cols = as.dist(1 - ibs_matrix),
#                    clustering_method = 'complete')

ibs <- broom::tidy(as.dist(ibs_matrix)) %>%  
  magrittr::set_colnames(c("p1", "p2", "ibs"))

ibs <- dplyr::left_join(ibs, mtdt_x, by = "p1") %>%
  dplyr::left_join(., mtdt_y, by = "p2") %>%
  rename_with(., ~gsub(".x", "_p1", .x, fixed = T)) %>%
  rename_with(., ~gsub(".y", "_p2", .x, fixed = T))

df <- ibs %>% 
  mutate(population = ifelse(Location_p1 == Location_p2, Location_p1, NA)) %>%
  group_by(population) %>% 
  summarise(N = n(), 
            mean_ibs = mean(ibs),
            median_ibs = median(ibs),
            sdt = sd(ibs),
            .groups = "drop") %>% 
  arrange(desc(mean_ibs)) %>%
  left_join(., coordinates, by = c("population" = "Locations")) %>% 
  drop_na()

# -------------------------
# Hotspot Identification: 
# -------------------------

# Spatial analysis to identify hotspots of high transmission rates or clusters 
# of genetically similar infections involves using statistical methods and 
# indices that can detect areas with significantly higher incidences of disease or 
# genetic similarity than would be expected by chance. 
# This analysis can be crucial for targeting interventions 
# and understanding the spatial dynamics of infectious diseases. 
# Here are some commonly used indices and methods for spatial hotspot identification:

# 1. Getis-Ord Gi (Hotspot Analysis)*
# The Getis-Ord Gi* statistic identifies hotspots by comparing local mean values to the overall mean.
# It determines the degree of clustering for high values (hotspots) or low values (cold spots) in spatial data.
# High Gi* values indicate that a location and its neighbors have higher values than expected,
# highlighting potential hotspots of transmission or genetic clusters.
# 
# The Getis Ord Gi* statistic is what is used in "Hot spot analysis."
# It looks at how local spatial variance in a single variable, as mentioned by u/Drewddit.
# It was advanced by Artur Getis and Keith Ord. MGWR is built on Geographically weighted regression,
# which is built on OLS regression. GWR adds a spatial weighted matrix into a
# normal regression equation, as opposed to assuming that there is not local variance.
# Multi-scale and Geographically and Temporally Weighted Regression both change that
# matrix to weigh variables differently in an analysis.
# GWR was advanced primarily by Alexander Stewart Fotheringham.
# 
# While MGWR is as easy as pushing a button, I want to caution that GWR is an advanced statistic 
# with ALOT of assumptions you need to consider. MGWR is obviously more advanced. 
# If you are only doing county level analysis, you probably could stick with GWR. 
# I'd recommend ALOT of reading, as you need a strong stats background to really understand 
# these and get into spatial statistics.

df_sf <- st_as_sf(df, coords = c("longitude", "latitude"), 
                  crs = st_crs(shape.gmb), na.fail = FALSE, agr = "constant")

# We consider there is a spatial clustering pattern
# Assuming 'your_data' is an sf object with coordinates and a value column 'values'
coords <- st_centroid(st_geometry(df_sf), of_largest_polygon=TRUE)

values <- df_sf$mean_ibs

# Neighbors and weights
nb <- knn2nb(knearneigh(coords, k = 4), sym = TRUE)
lw <- nb2listw(nb, style = "W")

# Check if there is some level of clustering using Getis-Ord Gi*
spdep::globalG.test(values, lw, zero.policy = TRUE) 

# Compute Local Getis-Ord Gi*
gi_star_perm <- localG_perm(values, lw, zero.policy = TRUE, nsim=499)

# Compute Local Getis-Ord Gi*
gi_star <- localG(values, lw, zero.policy = TRUE)
attributes(gi_star)$internals

# Calculate the p-values for gi_star
pvals <- pnorm(q = 2*(abs(gi_star)), lower.tail = FALSE)

# Adjust the p-values bonferroni
pvals_bon <- spdep::p.adjustSP(p = pvals, nb = nb, method = "bonferroni")

# Adjust the p-values fdr
pvals_fdr <- spdep::p.adjustSP(p = pvals, nb = nb, method = "fdr")

# Plot hotspots
# plot(st_geometry(df_sf), col = ifelse(gi_star > quantile(gi_star, 0.95), "red", "white"))

# data_gi <- cbind(df_sf, attributes(gi_star)$internals) %>% 
#   rename("gi" = "as.vector.gi_star.")

df_sf$gi <- gi_star

# Plotting
p1 <- ggplot() +
  geom_sf(data = shape.gmb, fill = "white", color = '#000000', linewidth = 0.5) + # Plot the base map
  geom_sf(data = data_gi, aes(color = mean_ibs), size = 10, lwd = 0.15) + #, show.legend = FALSE
  scale_color_gradient(low = "pink", high = "darkred") +
  labs(title = "Mean IBS Values By Village", color = "Mean IBS") +
  ggrepel::geom_label_repel(data = df, aes(label = population, x = longitude, y = latitude), 
                            color = '#000000', size = 4, box.padding = unit(1, "lines"), 
                            segment.color = '#132B43',  max.overlaps = Inf) +
  theme_void() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 12))

p2 <- ggplot() +
  geom_sf(data = shape.gmb, fill = "white", color = '#000000', linewidth = 0.5) + # Plot the base map
  geom_sf(data = data_gi, aes(color = gi), size = 10, lwd = 0.15) + #, show.legend = FALSE
  scale_color_gradient(low = "pink", high = "darkred") +
  labs(title = "Getis Ord Gi* statistic showing Hotspot Villages", color = "Getis-Ord Gi*") +
  ggrepel::geom_label_repel(data = df, aes(label = population, x = longitude, y = latitude), 
                            color = '#000000', size = 4, box.padding = unit(1, "lines"), 
                            segment.color = '#132B43',  max.overlaps = Inf) +
  theme_void() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 12))

p1/p2

# 2. Moran's I
# Moran's I is a measure of spatial autocorrelation that evaluates whether the pattern expressed is clustered, 
# dispersed, or random across the study area. It compares the value of a variable at one location with the value 
# of the same variable at nearby locations. High positive Moran's I values indicate clustering of similar values 
# (e.g., high transmission rates or genetic similarity).

gen_samples_sp <- as(df_sf, "Spatial")

# Define neighbors (e.g., using k-nearest neighbors)
coords <- coordinates(gen_samples_sp)
nb <- knn2nb(knearneigh(coords, k = 4), sym = TRUE)

# Calculate Moran's I for genetic similarity
lw <- nb2listw(nb, style = "W")
moran <- moran.test(df_sf$mean_ibs, lw)
print(moran)

# Plotting Moran's I scatterplot
moran.plot(df_sf$mean_ibs, lw)

# Compute local Moran
local <- localmoran(x = df_sf$mean_ibs, listw = lw)

# binds results to our polygon shapefile
moran.map <- cbind(df_sf, local)

ggplot() +
  geom_sf(data = shape.gmb, fill = "white", color = '#000000', linewidth = 0.5) + # Plot the base map
  geom_sf(data = moran.map, aes(color = Ii), size = 10, lwd = 0.15) + #, show.legend = FALSE
  scale_color_gradient(low = "pink", high = "darkred") +
  labs(title = "Local spatial autocorrelation", color = "local moran statistic") +
  ggrepel::geom_label_repel(data = df, aes(label = population, x = longitude, y = latitude), 
                            color = '#000000', size = 4, box.padding = unit(1, "lines"), 
                            segment.color = '#132B43',  max.overlaps = Inf) +
  theme_void() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 12))


# 3. Geary's C
# Geary's C is another measure of spatial autocorrelation, similar to Moran's I, but it is more sensitive to differences 
# between neighboring locations. It focuses on the dissimilarities between adjacent areas. 
# Values significantly lower than 1 indicate spatial clustering of similar values.
# 
# 
# 4. Kulldorff’s Spatial Scan Statistic
# Used in SaTScan software, this method scans the study area for clusters by moving a circular window across 
# the area and evaluating the likelihood of observing the number of cases within the window by chance. 
# It's widely used for detecting disease outbreaks and can be adapted 
# for identifying clusters of genetically similar infections.
# 
# Performing Spatial Analysis in R
# For spatial analysis in R, you can use packages like sp, sf (for handling spatial data), 
# spdep (for spatial dependence, including Moran's I and Geary's C), and spatialEco (for Getis-Ord Gi*).
# Assuming you have a spatial object `disease_data` with transmission rates or genetic similarity index



# Check if the clustering is significant
global_g_test(df_sf$mean_ibs, nb, lw)
