
# if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")

# https://biostatsquid.com/step-by-step-heatmap-tutorial-with-pheatmap/
# https://yunranchen.github.io/intro-net-r/igraph.html
# https://rpubs.com/writetosamadalvi/CommunityDetection

library(igraph)
library(tidyverse)
library(RColorBrewer)

set.seed(27011988)

data <- readxl::read_xlsx("data/FinalgambianDataset.xlsx") %>% 
   dplyr::filter(Location != "Nkakat Eyamba" & Location != "Ngayen Sanjal") %>% 
   dplyr::rename(sampleID = `Sample Internal ID`)

genetic_similarity <- read_delim("spatial_analysis/gambia_ibs.tsv") %>% 
   as.matrix

samples <- colnames(genetic_similarity)
rownames(genetic_similarity) <- samples

# Order data according to matrix
data <- data[match(samples, data$sampleID),]

# ======= Annotations ======= 
populations <- data %>% 
   pull(Location) %>% 
   unique() %>% 
   sort()

# Create a data frame for column annotation
ann_df <- data.frame(Year = data$Year, study_sites = data$Location)
row.names(ann_df) <- samples

# location_df <- data.frame(study_sites = data$Location)
# row.names(location_df) <- samples

ann_colors <- list(
   study_sites = c(Bansang = "#7F00FF", Basse = "#d72613", Brikama = "#077b8a", Brikamaba = "#5c3c92", 
                   Fajikunda = "#801818", Fatoto = "#FFE500", Ijede = "#12a4d9", Kafuta = "#00FFFF",
                   Kerewan = "#ff6c40", Pirang = "#000075", Sabi = "#c4a35a", Soma = "#12e761", 
                   'Sotuma Sere' = "#001AFF", Yerobawol = "#FF004D"),
   Year = brewer.pal(length(unique(data$Year)), 'PiYG')
)

# ============================ Base heatmap =======================
heat_plot <- pheatmap::pheatmap(genetic_similarity, 
                      col = rev(brewer.pal(10, 'RdYlGn')), # choose a color scale for your data
                      cluster_rows = TRUE, cluster_cols = TRUE, # set to FALSE if you want to remove the dendograms
                      clustering_distance_cols = 'euclidean',
                      clustering_distance_rows = 'euclidean',
                      clustering_method = 'ward.D',
                      # annotation_row = location_df, # row (location) annotations
                      annotation_col = ann_df, # column (sample) annotations
                      annotation_colors = ann_colors, # colors for your annotations
                      annotation_names_row = FALSE, 
                      annotation_names_col = FALSE,
                      border_color = "black", # default is grey60
                      number_color = "black",
                      fontsize_row = 10,     # row label font size
                      fontsize_col = 7,      # column label font size 
                      angle_col = 45, # sample names at an angle
                      legend_breaks = c(0.3, 0.7, 0.9), # legend customization
                      legend_labels = c("Low", "Medium", "High"), # legend customization
                      show_colnames = FALSE, show_rownames = FALSE, # displaying column and row names
                      main = "Pairwise IBS heatmap") # a title for our heatmap


# ============================ Network plot =======================
# Create a graph from the adjacency matrix
g <- graph.adjacency(genetic_similarity, mode = "undirected", weighted = TRUE, diag = FALSE)

# Set names and other attributes if necessary
V(g)$name <- samples # Assuming row names are sample IDs
V(g)$location <- data$Location

# Generate colors based on Location:
colours <- c("#7F00FF", "#d72613", "#077b8a", "#5c3c92", "#801818", "#FFE500", "#12a4d9",
             "#00FFFF", "#ff6c40", "#000075", "#c4a35a", "#12e761", "#001AFF", "#FF004D")
             
names(colours) <- populations

V(g)$colours <- colours[V(g)$location]

# ===== Calculate network properties =========
# degree_distribution <- degree(g)
# cluster_coefficient <- transitivity(g, type = "local")
# 
# # ========= Visualize the network ==============
# plot(g, vertex.size = 5, edge.width = E(g)$weight, vertex.label = NA, layout = layout_with_fr(g),
#      vertex.color = as.factor(V(g)$location), edge.curved = .9)

# ========================================================
# Retain edges with a similarity above a certain threshold
threshold <- 0.8 # Adjust based on your similarity measure and analysis needs
g1 <- delete_edges(g, E(g)[E(g)$weight < threshold])

# ===== Calculate network properties =========
degree_distribution <- degree(g1)
cluster_coefficient <- transitivity(g1, type = "local")

# ===== Get different layouts =========
layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 

# Remove layouts that do not apply to our graph.
layouts <- layouts[grepl("components|dh|nicely|fr|graphopt|kk|mds", layouts)]

# ========= Visualize the network ==============

par(mfrow=c(3,3), mar=c(1,1,1,1))
for (layout in layouts) {
   
   print(layout)
   
   l <- do.call(layout, list(g1)) 
   
   plot(g1, vertex.size = 5, edge.width = E(g1)$weight, vertex.label = NA, main = layout, 
        layout = l, vertex.color = as.factor(V(g1)$location))
}


legend("bottom", legend = populations, pch = 21, col = colours, pt.bg = colours,
       pt.cex = 2, cex = .7, bty = "n", ncol = 3, horiz  =  FALSE)

# ===================================================
# Community detection algorithms can help
# identify clusters within the network 
# that represent closely related groups of samples

# https://people.duke.edu/~jmoody77/snh/2021/CommunitiesSNH2021.nb.html
# ===================================================
# clusterlouvain <- cluster_louvain(g)
# plot(g, vertex.color = rainbow(20, alpha=0.6)[clusterlouvain$membership])

clusterlouvain <- cluster_louvain(g1)

plot(clusterlouvain, g1, vertex.color = colours[clusterlouvain$membership],
     vertex.size = 5, edge.width = E(g1)$weight, vertex.label = NA, 
     edge.arrow.size = .2, layout = layout_components(g1))


colors <- rainbow(max(clusterlouvain$membership))
plot(clusterlouvain, g1, vertex.color = colors,
     vertex.size = 5, edge.width = E(g1)$weight, vertex.label = NA, 
     edge.arrow.size = .2, layout = layout_components(g1))

#Step 1: find communities on karate
k1 = cluster_walktrap(g1)
k2 = cluster_infomap(g1)
# k3 = cluster_edge_betweenness(g1)

colors <- rainbow(max(k1$membership))
plot(g1, vertex.color = colors[k1$membership],
     vertex.size = 5, edge.width = E(g1)$weight, vertex.label = NA, 
     edge.arrow.size = .2, layout = layout_with_fr(g1))

colors <- rainbow(max(k2$membership))
plot(g1, vertex.color = colors[k2$membership],
     vertex.size = 5, edge.width = E(g1)$weight, vertex.label = NA, 
     edge.arrow.size = .2, layout = layout_with_sugiyama(g1))


genetic_similarity <- read.table("spatial_analysis/jaccard_matrix.tsv", header = TRUE) %>% 
   as.matrix

similarity_network <- function(matrix, plot = FALSE, layout = NULL, threshold = NULL){
   samples <- colnames(genetic_similarity)
   
   barcode.data <- readxl::read_xlsx("../../01_data/FinalgambianDataset.xlsx") %>% 
      dplyr::filter(Location != "Nkakat Eyamba" & Location != "Ngayen Sanjal") %>% 
      dplyr::rename(sampleID = `Sample Internal ID`)
   
   # Order data according to matrix
   barcode.data <- data[match(samples, data$SampleID),]
   
   # Create network
   g <- graph.adjacency(genetic_similarity, mode = "undirected", weighted = TRUE, diag = FALSE)
   
   # Set names and other attributes if necessary
   nbre_sites <- length(unique(barcode.data$Location))
   V(g)$name <- colnames(genetic_similarity) # Assuming row names are sample IDs
   V(g)$location <- barcode.data$Location
   V(g)$colours <- rainbow(nbre_sites)
   
   if(plot){
      # ========= Visualize the network ==============
      plot(g, vertex.size = 5, edge.width = E(g)$weight, vertex.label = NA,
           layout = layout, vertex.color = as.factor(V(g)$location))
      
      legend("topright", legend = unique(barcode.data$Location), pch = 21, 
             col = colours, pt.bg = colours,
             pt.cex = 1.2, cex = .5, bty = "n", ncol = 1, horiz  =  FALSE)
   }
   
   if(!is.null(threshold)){
      thres <- threshold # Adjust based on your similarity measure and analysis needs
      g1 <- delete_edges(g, E(g)[E(g)$weight < thres])
      
      # ========= Visualize the network ==============
      plot(g1, vertex.size = 5, edge.width = E(g1)$weight, vertex.label = NA,
           layout = layout, vertex.color = as.factor(V(g1)$location))
      
      legend("topright", legend = unique(barcode.data$Location),
             pch = 21, col = colours, pt.bg = colours,
             pt.cex = 1.2, cex = .5, bty = "n", ncol = 1, horiz  =  FALSE)
   }
   
   clusterlouvain <- cluster_louvain(g1)
   
   # Assign a color to each cluster
   # Here, we generate a color palette with a distinct color for each cluster
   colors <- rainbow(max(clusterlouvain$membership))
   
   # Convert membership IDs to labels (e.g., Community 1, Community 2, ...)
   # This step is optional and depends on whether you want custom names for communities
   community_labels <- paste("Community", clusterlouvain$membership)
   
   # Assign the labels to the graph vertices
   V(g1)$label <- community_labels
   
   # Plot the graph with vertex colors based on community and labels showing community names
   plot(g1, vertex.color = colors[clusterlouvain$membership],
        vertex.size = 7, edge.width = E(g1)$weight, 
        vertex.label = NA, vertex.label.color = "black",  # Adjust label color as needed
        vertex.label.cex = 0.8)
   
   legend("topright", legend = unique(community_labels), 
          pch = 21, col = colors, pt.bg = colors,
          pt.cex = 1.5, cex = .7, bty = "n", ncol = 1, horiz  =  FALSE)
   
}


