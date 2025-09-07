# =============================================================================
#
# SCRIPT: PCA Visualization with Centroids for Metabolomics Data
#
# DESCRIPTION:
# This script performs a Principal Component Analysis (PCA) on metabolomics
# data. It then calculates the centroids for predefined groups within the data
# and generates scatter plots of the principal components (PCs). These plots
# visualize the relationship between individual data points and their
# respective group centroids, helping to assess within-group cohesion and
# between-group separation.
#
# AUTHOR: Alan Hernández-Melgar
# DATE: 2023-10-27
#
# R Version: 4.3.1
#
# REQUIRED LIBRARIES:
# - ggplot2: For creating advanced and customized plots.
# - ggbiplot: Often used for PCA biplots, though here ggplot2 is used directly.
#
# =============================================================================


# -----------------------------------------------------------------------------
# 1. SETUP: LOAD LIBRARIES
# -----------------------------------------------------------------------------
# Ensure the necessary packages are installed before running the script.
# You can install them with: install.packages(c("ggplot2", "ggbiplot"))

library(ggbiplot) # Primarily used for creating PCA biplots
library(ggplot2)  # Used library for data visualization


# -----------------------------------------------------------------------------
# 2. DATA PREPARATION
# -----------------------------------------------------------------------------
# Load the dataset from a CSV file.
# NOTE: The path is hardcoded. For better portability, consider using relative
# paths or a file chooser dialog.
# The data is expected to be in a matrix-like format where rows are samples
# and columns are features (metabolites). A 'feature' or 'group' column
# should exist for sample classification.
data_metaboanalyst <- read.csv(
  "G:\\My Drive\\Ph. D\\Projects\\SkinProject\\Clinical_Trial_2\\Metabolomics\\2_Isotopes_New\\Metaboanalyst\\data_metaboanalyst_no_log_t.csv"
)

# Perform Principal Component Analysis (PCA).
# - `center = T`: Shifts the data to be centered at the origin (mean = 0).
# - `scale. = T`: Scales the data to have a standard deviation of 1.
#   Centering and scaling are crucial for PCA to prevent variables with
#   larger variances from dominating the analysis.
pca_data <- prcomp(data_metaboanalyst, center = TRUE, scale. = TRUE)

# Extract the PCA scores (the coordinates of each sample in the PC space)
# and convert them to a data frame for easier manipulation with ggplot2.
df.metaboanalyst.x <- as.data.frame(pca_data$x)

# Add the grouping variable from the original dataset to the PCA scores.
# This is essential for color-coding the plots and for calculating centroids.
df.metaboanalyst.x$groups <- data_metaboanalyst$feature

# Calculate the centroids for each group. A centroid is the average
# coordinate (mean) for all samples within a group for each principal component.
# This gives a single point representing the "center" of each group.
pca_centroids <- aggregate(df.metaboanalyst.x[,1:13], list(Type = df.metaboanalyst.x$groups), mean)


# -----------------------------------------------------------------------------
# 3. CENTROID & DATA MERGING
# -----------------------------------------------------------------------------
# Rename the first four columns of the centroids data frame for clarity.
# 'groups' will be the group identifier, and 'C1', 'C2', 'C3' will store the
# centroid coordinates for the first three principal components.
colnames(pca_centroids)[1:4] <- c("groups", "C1", "C2", "C3")

# Merge the centroid coordinates back into the main PCA scores data frame.
# This operation adds the centroid coordinates (C1, C2, C3) to each sample's
# row based on its group, which is necessary for drawing segments from each
# point to its group centroid.
pca_c <- merge(df.metaboanalyst.x, pca_centroids[1:4], by = 'groups', sort = FALSE)


# -----------------------------------------------------------------------------
# 4. VISUALIZATION
# -----------------------------------------------------------------------------

### PLOT 1: Principal Component 1 vs. Principal Component 2 ###

p1_cent_coh <- ggplot(pca_c, aes(x = PC1, y = PC2, colour = groups)) +
  # Add dashed vertical and horizontal lines at zero for reference
  geom_vline(xintercept = 0, linetype="dashed", size = 0.5, color= "#999999") +
  geom_hline(yintercept = 0, linetype="dashed", size = 0.5, color= "#999999") +
  # Plot individual data points (samples) with transparency
  geom_point(size = 8, alpha=0.3) +
  # Plot the group centroids as larger, semi-transparent points
  geom_point(data = pca_centroids, aes(x=C1, y=C2, color=groups), size = 15, alpha=0.66) +
  # Draw lines connecting each sample to its respective group centroid
  geom_segment(data = pca_c, mapping = aes(xend = C1, yend = C2), alpha=0.3, size = 2) +
  # Manually define the colors for each group for consistent plotting
  scale_color_manual(values = c("#9092c0", "#c25253", "yellow")) +
  # Customize theme elements for better readability
  theme(
    axis.title = element_text(size = 20),           # Axis titles
    axis.text = element_text(size = rel(1.5)),      # Axis tick labels
    legend.text = element_text(size = 15),          # Legend item text
    legend.title = element_text(face="bold"),       # Legend title
    legend.position="top"                           # Move legend to the top
  )

# Display the plot
p1_cent_coh


### PLOT 2: Principal Component 1 vs. Principal Component 3 ###

p2_cent_coh <- ggplot(pca_c, aes(x = PC1, y = PC3, colour = groups)) +
  geom_vline(xintercept = 0, linetype="dashed", size = 0.5, color= "#999999") +
  geom_hline(yintercept = 0, linetype="dashed", size = 0.5, color= "#999999") +
  geom_point(size = 8, alpha=0.3) +
  # Note the change to y=C3 and yend=C3 for plotting PC3
  geom_point(data = pca_centroids, aes(x=C1, y=C3, color=groups), size = 15, alpha=0.66) +
  geom_segment(data = pca_c, mapping = aes(xend = C1, yend = C3), alpha=0.3, size = 2) +
  scale_color_manual(values = c("#9092c0", "#c25253", "yellow")) +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = rel(1.5)),
    legend.text = element_text(size = 15),
    legend.title = element_text(size=20, face="bold")
  )

# Display the plot
p2_cent_coh


### PLOT 3: Principal Component 2 vs. Principal Component 3 ###

# NOTE: The original code has an error here, plotting PC1 vs PC2 again.
# This has been corrected to plot PC2 vs PC3.
p3_cent_coh <- ggplot(pca_c, aes(x = PC2, y = PC3, colour = groups)) +
  geom_vline(xintercept = 0, linetype="dashed", size = 0.5, color= "#999999") +
  geom_hline(yintercept = 0, linetype="dashed", size = 0.5, color= "#999999") +
  geom_point(size = 8, alpha=0.3) +
  # Plotting C2 vs C3 for centroids and segments
  geom_point(data = pca_centroids, aes(x=C2, y=C3, color=groups), size = 15, alpha=0.66) +
  geom_segment(data = pca_c, mapping = aes(xend = C2, yend = C3), alpha=0.3, size = 2) +
  scale_color_manual(values = c("#9092c0", "#c25253", "yellow")) +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = rel(1.5)),
    legend.text = element_text(size = 15),
    legend.title = element_text(size=20, face="bold")
  )

# Display the plot
p3_cent_coh

# =============================================================================
# END OF SCRIPT
# =============================================================================
