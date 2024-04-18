
# if (report_progress) {
#     setTxtProgressBar(pbar, i)
# }

# Function to estimate the proportion of missing data for each individual
prop.missing <- function(data){
   
   # Barcode sequence length
   seq.length <- nchar(data$Barcode) %>% 
      unique()
    missing <- str_count(data$Barcode, pattern = "X")
    
    prop.missing <- 100* (missing/seq.length)
    return(prop.missing)
}

# Estimate missing data for each locus (2) across samples and for each individual (1)
# when using alleles at each SNPs positions data
# Set margin to 1 (samples) or 2 (loci)
calculate_missing_genotypes <- function(barcode_data, margin) {
    missing_genotypes <- apply(df, margin, function(locus) {
        # str_count(toString(data[1, -c(1:5)]), pattern = "X|N")
        seq.length <- length(locus)
        missing <- str_count(toString(as.character(locus)), pattern = "X")
        return((missing/seq.length)*100)
    })
    return(as.vector(missing_genotypes))
}


# Function to estimate Identity by State (IBS) between two sequences
estimate_IBS <- function(sequence1, sequence2) {
    # Ensure sequences are of the same length
    if (nchar(sequence1) != nchar(sequence2)) {
        stop("Sequences must be of the same length")
    }
    
    # Convert sequences to vectors of individual alleles
    alleles1 <- str_split_1(sequence1, pattern = "") # strsplit(sequence1, "")[[1]]
    alleles2 <- str_split_1(sequence2, pattern = "") # strsplit(sequence2, "")[[1]]
    
    # Identify positions with valid data in both sequences
    # valid_positions <- which(!(alleles1 %in% c("X", "N") | alleles2 %in% c("X", "N")))
    valid_positions <- which(!(alleles1 %in% "X" | alleles2 %in% "X"))
    
    if (length(valid_positions) == 0) return(NA)  # Return NA if no valid positions to compare
    
    # Subset alleles to only valid positions
    valid_alleles1 <- alleles1[valid_positions]
    valid_alleles2 <- alleles2[valid_positions]
    
    # Calculate matches and mismatches
    matches <- sum(valid_alleles1 == valid_alleles2)
    total_valid <- length(valid_positions)
    
    # Calculate IBS score
    ibs_score <- matches / total_valid
    
    return(ibs_score)
}

# Function to compute pairwise IBS for all individuals (WARNING: Computationally Intensive)
compute_pairwise_IBS_v1 <- function(barcode_data) {
    num_individuals <- nrow(barcode_data)
    # Initialize a matrix to store IBS scores
    ibs_matrix <- matrix(NA, nrow = num_individuals, ncol = num_individuals)
    
    for (i in 1:(num_individuals-1)) {
        for (j in (i+1):num_individuals) {
            seq1 <- barcode_data$Barcode[i]
            seq2 <- barcode_data$Barcode[j]
            
            ibs_matrix[i, j] <- estimate_IBS(seq1, seq2)
            ibs_matrix[j, i] <- ibs_matrix[i, j]  # Symmetric matrix
        }
    }
    return(ibs_matrix)
}

compute_pairwise_IBS_v2 <- function(barcode_data) {
    # Number of samples (individuals)
    nsamp <- nrow(barcode_data)
    
    # Initialize a matrix to store IBS scores
    similarity.matrix <- matrix(NA, ncol = nsamp, nrow = nsamp)
    
    rownames(similarity.matrix) = colnames(similarity.matrix) <- barcode_data %>%
        pull(sample_internal_id)
    
    for (i in 1:nsamp) {
        j <- i + 1
        while (j <= nsamp) {
            sequence1 <- barcode_data$Barcode[i]
            sequence2 <- barcode_data$Barcode[j]
            
            similarity.matrix[j, i] <- estimate_IBS(sequence1, sequence2)
            similarity.matrix[i, j] <- similarity.matrix[j, i]  # Symmetric matrix
            j = j+1
        }
    }
    
    return(similarity.matrix)
}


calculateIBS_v1 <- function(df)
{
   require(data.table)
   require(tictoc)
   ibs_matrix <- data.frame(matrix(NA, nrow = nrow(df), ncol = nrow(df)))
   
   tic()
   for (i in 1:nrow(df))
   {
      message("********* Processing sample row [", i, "] ******")
      
      ibs_matrix[i,i] <- 1
      sample1 <- as.matrix(df[i,])
      k <- i+1
      while((i<k) & (k <= nrow(df)))
      {
         missing <- 0
         s1_vs_s2 <- numeric(length(sample1))
         sample2 <- as.matrix(df[k,])
         for (j in 1:length(sample1))
         {
            #Split alts
            scol1 <- sample1[j]
            scol2 <- sample2[j]
            if (scol1 %in% c("X", "N") | scol2 %in% c("X", "N") )
               missing <- missing+1
            else if (scol1 != scol2)
               s1_vs_s2[j] <- 0
            else
               s1_vs_s2[j] <- 1
         }
         ibs <- sum(s1_vs_s2)/(ncol(df)-missing)
         ibs_matrix[i,k] <- ibs
         ibs_matrix[k,i] <- ibs
         k <- k+1
      }
      
   }
   toc()
   return(ibs_matrix)
}

calculateIBS_v2 <- function(df)
{
   require(data.table)
   
   ibs_matrix <- data.frame(matrix(NA, nrow = nrow(df), ncol = nrow(df)))
   
   start_time = Sys.time()
   for (i in 1:nrow(df))
   {
      message("********* Processing sample row [", i, "] ******")
      
      ibs_matrix[i,i] <- 1
      sample1 <- as.matrix(df[i,])
      k <- i+1
      while((i<k) & (k <= nrow(df)))
      {
         missing <- 0
         s1_vs_s2 <- numeric(length(sample1))
         sample2 <- as.matrix(df[k,])
         for (j in 1:length(sample1))
         {
            # scol1 <- sample1[j]
            # scol2 <- sample2[j]
            # if (scol1 %in% "X" | scol2 %in% "X") missing <- missing+1
            # else if ((scol1 == "N" & scol2 %in% c("A", "T", "C", "G")) |
            #          (scol2 == "N" & scol1 %in% c("A", "T", "C", "G"))) s1_vs_s2[j] <- 0.5
            # else if (scol1 != scol2) s1_vs_s2[j] <- 0
            # else s1_vs_s2[j] <- 1
            
            pair <- c(sample1[j], sample2[j])
            if("X" %in% pair) missing <- missing+1  # Skip missing data
            else if(all(pair == "N")) s1_vs_s2[j] <- 0.5  # Both are N
            else if("N" %in% pair) s1_vs_s2[j] <- 0.5  # One is N
            else if(pair[1] == pair[2]) s1_vs_s2[j] <- 1  # Identical alleles
            # if(pair[1] != pair[2]) s1_vs_s2[j] <- 0
            else s1_vs_s2[j] <- 0
            
         }
         
         ibs <- round(sum(s1_vs_s2)/(ncol(df)-missing), 2)
         ibs_matrix[i,k] <- ibs
         ibs_matrix[k,i] <- ibs
         k <- k+1
      }
      
   }
   
   end_time = Sys.time()
   
   message("********* Processing time = ", (end_time - start_time))
   return(ibs_matrix)
}

calculateIBS2 <- function(input) {
   
   # Convert input to matrix and assign to 'IBS_Matrix'
   file = as.matrix(input)
   n = nrow(file)
   
   # Create an empty matrix to store IBS values
   ibs_matrix = matrix(1, nrow=n, ncol=n)
   rownames(ibs_matrix) = rownames(file)
   colnames(ibs_matrix) = rownames(file)
   
   start_time = Sys.time()
   for (i in 1:(n-1)) {
      # Display progress message
      message(paste("******* Processing sample row [", i , "] *******"))
      
      # Compare sample i with subsequent samples
      for (k in (i+1):n) {
         # Extract genotype data for sample i and k
         sample1 = file[i,]
         sample2 = file[k,]
         
         # Calculate the number of missing values and specific alleles ('X') for each pair of samples
         missing = sum(sample1 == "X" | sample2 == "X")
         
         # Calculate the proportion of matching alleles between sample1 and sample2
         s1_vs_s2 <- ifelse((sample1 == sample2 & sample1 != "X" & sample2 != "X"), 1,
                            ifelse((sample1 == "N" & sample2 %in% c("A", "C", "G", "T")) |
                                      (sample2 == "N" & sample1 %in% c("A", "C", "G", "T")), 0.5, 0)
         )
         
         # Calculate IBS value for the pair
         ibs = round(sum(s1_vs_s2, na.rm=TRUE) / (ncol(file) - missing), 2)
         
         # Store IBS value in the matrix for both positions
         ibs_matrix[i,k] <- ibs
         ibs_matrix[k,i] <- ibs
      }
   }
   
   end_time = Sys.time()
   
   message("********* Processing time = ", (end_time - start_time))
   return(ibs_matrix)
}



calculate_allelic_richness <- function(barcode_data) {
   allelic_richness <- apply(barcode_data, 2, function(locus) {
      alleles <- unlist(strsplit(as.character(locus), ""))
      unique_alleles <- unique(alleles[!alleles %in% c("X", "N")]) # Exclude missing data
      return(length(unique_alleles))
   })
   return(round(sum(allelic_richness)/ncol(barcode_data), 2))
}


# Compute SNP frequencies
calculate_snp_frequencies <- function(barcode_data) {
   snp_frequencies <- lapply(1:ncol(barcode_data), function(locus_index) {
      alleles <- barcode_data[[locus_index]]
      # alleles <- unlist(strsplit(as.character(locus), "\t"))
      valid_alleles <- alleles[!alleles %in% c("X", "N")] # Exclude missing data
      allele_counts <- table(valid_alleles)
      total_alleles <- sum(allele_counts)
      frequencies <- allele_counts / total_alleles
      return(frequencies)
   })
   names(snp_frequencies) <- colnames(barcode_data)
   return(snp_frequencies)
}


# --------------------------
# JACCARD SIMILARITY INDEX
# --------------------------
calculate_jaccard_alleles <- function(sample1, sample2) {

   shared = 0
   total = 0
   
   for (i in 1:length(sample1)) {
      alleles1 = as.character(sample1[i])
      alleles2 = as.character(sample2[i])
      
      if (!(alleles1  %in% c("N", "X")) & !(alleles2 %in% c("N", "X"))) {
         total = total + 1
         if (length(intersect(alleles1, alleles2)) > 0) {
            shared = shared + 1
         }
      }
   }
   
   jaccard_similarity = round(shared / total, 2)
   return(jaccard_similarity)
}

# Example usage with two samples from a genetic_data dataframe

compute_jaccard_index <- function(barcode_data, report_progress = TRUE) {
   library(progress)
   library(tcltk)
   
   num_individuals <- nrow(barcode_data)
   
   if (report_progress) {
      pbar <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                             max = num_individuals, # Maximum value of the progress bar
                             style = 3,    # Progress bar style (also available style = 1 and style = 2)
                             width = 100,   # Progress bar width. Defaults to getOption("width")
                             char = "=")   # Character used to create the bar
   }
   
   # Initialize a matrix to store Jaccard scores
   jaccar_matrix <- matrix(NA, nrow = num_individuals, ncol = num_individuals)
   diag(jaccar_matrix) <- 1
   
   for (i in 1:(num_individuals-1)) {
      for (j in (i+1):num_individuals) {
         sample1 <- data[i, ]
         sample2 <- data[j, ]
         
         jaccard_similarity = calculate_jaccard_alleles(sample1, sample2)
         # print(jaccard_similarity)
         
         jaccar_matrix[i, j] <- jaccard_similarity
         jaccar_matrix[j, i] <- jaccar_matrix[i, j]  # Symmetric matrix
      }
      
      if (report_progress) {
         setTxtProgressBar(pbar, i+1)
      }
   }
   return(jaccar_matrix)
}


# -----------------------------
#
# -----------------------------
calculate_dice_alleles <- function(sample1, sample2) {
   shared = 0
   total = 0
   
   for (i in 1:length(sample1)) {
      alleles1 = as.character(sample1[i])
      alleles2 = as.character(sample2[i])
      
      if (!(alleles1 %in% c("N", "X")) & !(alleles2 %in% c("N", "X"))) {
         total = total + length(unique(c(alleles1, alleles2)))
         shared = shared + length(intersect(alleles1, alleles2)) * 2
      }
   }
   
   dice_similarity = round(shared / total, 3)
   return(dice_similarity)
}

compute_dice_index <- function(barcode_data, report_progress = TRUE) {
   num_individuals <- nrow(barcode_data)
   
   if (report_progress) {
      pbar <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                             max = num_individuals, # Maximum value of the progress bar
                             style = 3,    # Progress bar style (also available style = 1 and style = 2)
                             width = 100,   # Progress bar width. Defaults to getOption("width")
                             char = "=")   # Character used to create the bar
   }
   
   # Initialize a matrix to store Jaccard scores
   dice_matrix <- matrix(NA, nrow = num_individuals, ncol = num_individuals)
   diag(dice_matrix) <- 1
   
   for (i in 1:(num_individuals-1)) {
      for (j in (i+1):num_individuals) {
         sample1 <- data[i, ]
         sample2 <- data[j, ]
         
         dice_similarity = calculate_dice_alleles(sample1, sample2)
         # print(dice_similarity)
         
         dice_matrix[i, j] <- dice_similarity
         dice_matrix[j, i] <- dice_matrix[i, j]  # Symmetric matrix
      }
      
      if (report_progress) {
         setTxtProgressBar(pbar, i+1)
      }
      
   }
   return(dice_matrix)
}

# -----------------------------------------
# Prevalence of drug resistance mutations
# ------------------------------------------


# Define the function to generate haplotypes from a single string
generate_haplotypes <- function(haplotype) {
   # Example: haplotype="CV[M/I][N/E][K/T]"
   # Define the input string
   # input_string <- "CV M/I N/E K/T"
   input_string <- stringr::str_squish(gsub("\\[|\\]", " ", haplotype))
   
   # Split the string into elements based on spaces
   elements <- unlist(strsplit(input_string, " "))
   
   # Initialize a list to hold options for each element
   options_list <- list()
   
   # Process each element to extract options
   for (i in 1:length(elements)) {
      element <- elements[i]
      # Check if the element contains a '/', indicating multiple options
      if (grepl("/", element)) {
         # Split the element into its options and add to the list
         options_list[[i]] <- unlist(strsplit(element, "/"))
      } else {
         # Element does not have options, add it directly
         options_list[[i]] <- element
      }
   }
   
   # Use expand.grid to generate all possible combinations of the options
   combinations_df <- expand.grid(options_list)
   
   # Combine the elements of each combination into a single string
   haplotypes <- apply(combinations_df, 1, paste, collapse = "")
   
   # Print the haplotypes
   return(haplotypes)
}


# Updated function to extract and count mutations or allele combinations
extract_and_count_mutations <- function(column_data, gene) {
   # Initialize a vector to hold all mutations and combinations
   all_mutations <- c()
   
   for (haplotype in column_data) {
      # Check for combination format and expand if necessary
      if (grepl("\\[[A-Z]/[A-Z]\\]", haplotype)) {
         expanded <- generate_haplotypes(haplotype)
         all_mutations <- c(all_mutations, expanded)
      } else {
         # Split entries with multiple mutations and add them to the list
         mutations <- unlist(strsplit(haplotype, " "))
         all_mutations <- c(all_mutations, mutations)
      }
   }
   
   # Count occurrences of each mutation
   mutation_counts <- table(all_mutations)
   
   # Calculate prevalence as a percentage
   total_samples <- length(column_data)
   mutation_prevalence <- prop.table(mutation_counts) * 100
   
   return(tibble(Drug = gene, Haplotypes = names(mutation_counts),
                 "Counts" = as.vector(mutation_counts), 
                 "Prevalence" = as.vector(mutation_prevalence)))
}

# -------------------------------------------
# Estimate prevalence of mutations over time
# -------------------------------------------
# Function to expand combination haplotypes and extract mutations (adapted for time analysis)
expand_and_extract_mutations_time <- function(barcode, gene) {
   
   barcode %>%
      group_by(Year) %>%
      do(extract_count_mutations_over_time(., gene)) %>% 
      ungroup()
}


# Updated function to extract and count mutations or allele combinations
extract_count_mutations_over_time <- function(barcode, gene) {
   
   column_data <- barcode[[gene]]
   
   # Initialize a vector to hold all mutations and combinations
   all_mutations <- c()
   
   for (haplotype in column_data) {
      # Check for combination format and expand if necessary
      if (grepl("\\[[A-Z]/[A-Z]\\]", haplotype)) {
         expanded <- generate_haplotypes(haplotype)
         all_mutations <- c(all_mutations, expanded)
      } else {
         # Split entries with multiple mutations and add them to the list
         mutations <- unlist(strsplit(haplotype, " "))
         all_mutations <- c(all_mutations, mutations)
      }
   }
   
   # Count occurrences of each mutation
   mutation_counts <- table(all_mutations)
   
   # Calculate prevalence as a percentage
   total_samples <- length(column_data)
   mutation_prevalence <- prop.table(mutation_counts) * 100
   
   return(tibble(Haplotypes = names(mutation_counts),
                 "Counts" = as.vector(mutation_counts), 
                 "Prevalence" = as.vector(mutation_prevalence)))
}
