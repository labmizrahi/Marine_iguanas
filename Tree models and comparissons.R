# Load necessary libraries
library(ape)         # For reading trees and calculating distances
library(dplyr)       # For data manipulation
library(ggplot2)     # For plotting
library(pheatmap)    # For heatmaps

# --- SETUP DIRECTORY AND EXPORTS ---
run_date <- ""
out_dir <- paste0(run_date, "_diversification_analysis")
if(!dir.exists(out_dir)) dir.create(out_dir)

# Step 1: Load the phylogenetic tree
tree_file <- "path/to/your/tree_file.newick"
tree <- read.tree(tree_file)

# Step 2: Calculate pairwise distances between sequences in the tree
distances <- cophenetic(tree)  # Generates a distance matrix from the tree

# Convert the distance matrix to a data frame
distances_df <- as.data.frame(distances)

# Ensure all data is numeric
distances_df[] <- lapply(distances_df, as.numeric)

# Step 3: Read target sequence names
target_seq_file <- "path/to/your/target_sequences.csv"
target_seq_names <- read.csv(target_seq_file, header = FALSE, stringsAsFactors = FALSE)
target_seq_names <- unlist(target_seq_names)

# Filter target sequences present in the distance data
target_seq_names <- target_seq_names[target_seq_names %in% colnames(distances_df)]

# Identify non-target sequence names
non_target_seq_names <- setdiff(colnames(distances_df), target_seq_names)

# Initialize data frame to store results
target_data <- data.frame(
  sequence = target_seq_names,
  min_distance = NA,
  closest_sequence_header = NA,
  n_neighbours = NA,
  stringsAsFactors = FALSE
)

# Step 4: Calculate minimum distances and number of target neighbours
for (i in seq_along(target_seq_names)) {
  seq_name <- target_seq_names[i]
  
  # Get distances from current target sequence to all sequences, explicitly keeping rownames
  seq_distances <- distances_df %>% 
    tibble::rownames_to_column("row_names") %>%
    dplyr::select(row_names, all_of(seq_name))
  
  same_seqs <- seq_distances %>% filter(.data[[seq_name]] < 0.005) 
  same_sequence_number <- same_seqs %>% nrow()
  
  if (same_sequence_number > 1){
    not_endemic_count <- not_endemic_count + 1
  }
  
  # Exclude distances to itself while keeping row names
  seq_distances <- seq_distances %>% filter(.data[[seq_name]] > 0)
  
  # Ensure we only consider non-target sequences present in seq_distances
  valid_non_target <- intersect(seq_distances$row_names, non_target_seq_names)
  
  if (length(valid_non_target) > 0) {
    # Minimum distance to non-target sequences
    valid_subset <- seq_distances %>% filter(row_names %in% valid_non_target)
    min_dist_non_target <- min(valid_subset[[seq_name]], na.rm = TRUE)
    
    # Find the header of the closest sequence
    closest_non_target_header <- valid_subset$row_names[which(valid_subset[[seq_name]] == min_dist_non_target)][1]
  } else {
    min_dist_non_target <- NA
    closest_non_target_header <- NA
  }
  
  # Number of target neighbours closer than min_dist_non_target
  target_subset <- seq_distances %>% filter(row_names %in% target_seq_names)
  n_neighbours <- sum(target_subset[[seq_name]] < min_dist_non_target, na.rm = TRUE)
  
  # Store results
  target_data$min_distance[i] <- min_dist_non_target
  target_data$closest_sequence_header[i] <- closest_non_target_header
  target_data$n_neighbours[i] <- n_neighbours
}

write.csv(target_data, file.path(out_dir, "observed_data.csv"), row.names = FALSE)

# Step 5: Plot observed data
p1 <- ggplot(target_data, aes(x = min_distance, y = n_neighbours)) +
  geom_point(color = "#AEC6CF") + # Pastel blue
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, color = "#FFB347") + # Pastel orange
  theme_classic() +
  xlab("Minimum Distance to Non-Target Sequences") +
  ylab("Number of Closer Target Sequences")

# Save plot using cairo_pdf so text is editable in Inkscape
ggsave(file.path(out_dir, "observed_plot.pdf"), plot = p1, width = 6, height = 4, device = cairo_pdf)
ggsave(file.path(out_dir, "observed_plot.png"), plot = p1, width = 6, height = 4, dpi = 300)

# Step 7: Generate random distribution for comparison
set.seed(123)  # For reproducibility
random_sample_size <- nrow(target_data)
min_dist_random <- numeric(random_sample_size)
n_neighbours_random <- numeric(random_sample_size)

for (i in 1:random_sample_size) {
  # Randomly select a sequence
  random_seq <- sample(colnames(distances_df), 1)
  
  # Get distances from the random sequence to all sequences
  seq_distances <- distances_df[, random_seq, drop = FALSE]
  
  # Exclude distances to itself
  seq_distances <- seq_distances[seq_distances > 0, , drop = FALSE]
  
  # Sample a random minimum distance from observed min distances
  min_dist <- sample(target_data$min_distance, 1)
  min_dist_random[i] <- min_dist
  
  # Number of neighbours closer than min_dist
  n_neighbours <- sum(seq_distances < min_dist, na.rm = TRUE)
  n_neighbours_random[i] <- n_neighbours
}

# Step 8: Build a null model by repeating random sampling
iterations <- 1000
null_model_data <- data.frame(min_distance = numeric(), n_neighbours = numeric())

for (iter in 1:iterations) {
  for (i in 1:random_sample_size) {
    random_seq <- sample(colnames(distances_df), 1)
    seq_distances <- distances_df[, random_seq, drop = FALSE]
    seq_distances <- seq_distances[seq_distances > 0, , drop = FALSE]
    min_dist <- sample(target_data$min_distance, 1)
    n_neighbours <- sum(seq_distances < min_dist, na.rm = TRUE)
    
    # Append to null model data
    null_model_data <- rbind(null_model_data, data.frame(min_distance = min_dist, n_neighbours = n_neighbours))
  }
}
write.csv(null_model_data, file.path(out_dir, "null_model_data.csv"), row.names = FALSE)

# Step 11: Fit a model to the null distribution
null_model <- lm(min_distance ~ poly(n_neighbours, 2), data = null_model_data)

# part 2:

# Step 1: Extract Coefficients 
coefficients <- coef(null_model) # Changed from non_linear_model_tree to null_model

# Step 3: Define a Function for Predictions
predict_y <- function(x_values, coefficients) {
  y_values <- coefficients[1] + coefficients[2] * x_values + coefficients[3] * x_values^2
  return(y_values)
}

x_subset <- target_data$min_distance 

# Step 4: Predict y Values for x_subset
exp_n_neighbours <- predict_y(x_subset, coefficients)
delta_n_neighbours <- target_data$n_neighbours - exp_n_neighbours

exp_obs_n_neighbours <- cbind(target_data$n_neighbours, exp_n_neighbours, delta_n_neighbours) %>% as.data.frame()

colnames(exp_obs_n_neighbours)[1] <- "n_neighbours"

p2 <- ggplot(data = exp_obs_n_neighbours, aes(x = n_neighbours, y = exp_n_neighbours)) +
  geom_point(color = "#AEC6CF") +
  geom_smooth(method = "lm", formula = y ~ poly(x, degree = 2), se = FALSE, color = "#FFB347") + 
  theme_classic()

ggsave(file.path(out_dir, "expected_vs_observed.pdf"), plot = p2, width = 6, height = 4, device = cairo_pdf)
ggsave(file.path(out_dir, "expected_vs_observed.png"), plot = p2, width = 6, height = 4, dpi = 300)

# adding it to the original table and plotting correlations to check for all
exp_obs_n_neighbours_merged <- cbind(exp_obs_n_neighbours, target_data)
write.csv(exp_obs_n_neighbours_merged, file.path(out_dir, "final_merged_data.csv"), row.names = FALSE)

