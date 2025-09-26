# Export Seurat data for scCODA analysis
library(Seurat)
library(dplyr)

# Assuming your Seurat object is called 'PD'

# 1. CREATE SAMPLE IDENTIFIERS
# Create unique sample identifiers
PD@meta.data$sample_id <- paste(PD@meta.data$Condition, 
                                PD@meta.data$Genotype, 
                                PD@meta.data$Timepoint, 
                                PD@meta.data$Treatment, sep = "_")

# 2. CREATE CELL COUNT MATRIX
# Count cells per sample and cell type
cell_counts <- PD@meta.data %>%
  group_by(sample_id, microglia_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = microglia_type, values_from = count, values_fill = 0)

# Convert to matrix with samples as rows
cell_counts_matrix <- as.data.frame(cell_counts)
rownames(cell_counts_matrix) <- cell_counts_matrix$sample_id
cell_counts_matrix <- cell_counts_matrix[, -1]  # Remove sample_id column

# 3. CREATE SAMPLE METADATA
library(dplyr)  # if not loaded yet

sample_metadata <- PD@meta.data %>%
  as.data.frame() %>%
  dplyr::select(sample_id, Condition, Genotype, Timepoint, Treatment) %>%
  distinct() %>%
  arrange(sample_id)
# Ensure order matches
sample_metadata <- sample_metadata[match(rownames(cell_counts_matrix), sample_metadata$sample_id), ]

# 4. EXPORT DATA
# Save as CSV files for Python
write.csv(cell_counts_matrix, "microglia_cell_counts.csv", row.names = TRUE)
write.csv(sample_metadata, "sample_metadata.csv", row.names = FALSE)

# 5. ALTERNATIVE: SAVE AS H5AD FOR SCANPY
# If you prefer to work with h5ad format
library(SeuratDisk)

# Convert to h5ad (requires SeuratDisk)
SaveH5Seurat(PD, filename = "microglia_data.h5seurat")
Convert("microglia_data.h5seurat", dest = "h5ad")

print("Data exported successfully!")
print(paste("Cell count matrix shape:", nrow(cell_counts_matrix), "x", ncol(cell_counts_matrix)))
print(paste("Number of samples:", nrow(sample_metadata)))
print(paste("Cell types:", paste(colnames(cell_counts_matrix), collapse = ", ")))
print("\nSample metadata:")
print(sample_metadata)
print("\nCell counts preview:")
print(head(cell_counts_matrix))



library(SeuratDisk)

# Convert to h5ad (requires SeuratDisk)
DefaultAssay(PD) <- "RNA"
SaveH5Seurat(PD, filename = "microglia_data.h5seurat", verbose = TRUE)
Convert("microglia_data.h5seurat", dest = "h5ad")

library(zellkonverter)

# Convert Seurat object to SingleCellExperiment
sce <- as.SingleCellExperiment(PD)
writeH5AD(sce, "microglia_data.h5ad")


kj


DotPlot(PD, feature = "CD83")
DotPlot(PD, feature = "CD83", group.by = "Genotype")
DotPlot(PD, feature = "CD83", group.by = "Timepoint")
DotPlot(PD, feature = "CD83", group.by = "Condition")
DotPlot(PD, feature = "CD83", group.by = "Treatment")

FeaturePlot(PD, feature = "CD83")

DefaultAssay(PD) = "SCT"



#publication ready prop plots
library(ggplot2)
library(dplyr)

# Define exotic microglia color palette
microglia_colors_exotic <- c(
  "Antigen-Presenting Microglia"      = "#E63946",  # Flamingo red – urgent immune alert
  "Disease-Associated Microglia"      = "#6A4C93",  # Mystic violet – neurodegenerative vibes
  "Homeostatic Microglia"             = "#2A9D8F",  # Deep jade green – peaceful balance
  "Interferon-Responsive Microglia"   = "#4CC9F0",  # Electric cyan – antiviral spark
  "Proliferating Microglia"           = "#F4A261"   # Spicy mango – hot mitotic energy
)

# Extract metadata
meta <- as.data.frame(PD@meta.data)

# Plot 1: Stacked barplot of celltype proportions by Timepoint
tab1 <- table(meta$Timepoint, meta$microglia_type)
prop1 <- prop.table(tab1, margin = 1) * 100
plot_data1 <- as.data.frame(prop1)
colnames(plot_data1) <- c("Timepoint", "microglia_type", "Percent")

p1 <- ggplot(plot_data1, aes(x = Timepoint, y = Percent, fill = microglia_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = microglia_colors_exotic) +
  labs(title = "Microglia subtype % by Timepoint", 
       y = "Percentage", 
       x = "Timepoint",
       fill = "Microglia Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Plot 2: Stacked barplot of celltype proportions by Genotype
tab2 <- table(meta$Genotype, meta$microglia_type)
prop2 <- prop.table(tab2, margin = 1) * 100
plot_data2 <- as.data.frame(prop2)
colnames(plot_data2) <- c("Genotype", "microglia_type", "Percent")

p2 <- ggplot(plot_data2, aes(x = Genotype, y = Percent, fill = microglia_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = microglia_colors_exotic) +
  labs(title = "Microglia subtype % by Genotype", 
       y = "Percentage", 
       x = "Genotype",
       fill = "Microglia Type") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Plot 3: Stacked barplot of celltype proportions by Condition
tab3 <- table(meta$Condition, meta$microglia_type)
prop3 <- prop.table(tab3, margin = 1) * 100
plot_data3 <- as.data.frame(prop3)
colnames(plot_data3) <- c("Condition", "microglia_type", "Percent")

p3 <- ggplot(plot_data3, aes(x = Condition, y = Percent, fill = microglia_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = microglia_colors_exotic) +
  labs(title = "Microglia subtype % by Condition", 
       y = "Percentage", 
       x = "Condition",
       fill = "Microglia Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Plot 4: Stacked barplot of celltype proportions by Treatment 
tab4 <- table(meta$Treatment, meta$microglia_type)
prop4 <- prop.table(tab4, margin = 1) * 100
plot_data4 <- as.data.frame(prop4)
colnames(plot_data4) <- c("Treatment", "microglia_type", "Percent")

p4 <- ggplot(plot_data4, aes(x = Treatment, y = Percent, fill = microglia_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = microglia_colors_exotic) +
  labs(title = "Microglia subtype % by Treatment", 
       y = "Percentage", 
       x = "Treatment",
       fill = "Microglia Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Display all plots
print(p1)
print(p2)
print(p3)
print(p4)

# Optional: Create a combined plot using patchwork
# library(patchwork)
# combined_plot <- (p1 + p2) / (p3 + p4)
# print(combined_plot)

# Save plots with higher resolution
ggsave("/data/outputs/cellranger_outputs/Analysis_by_Mo_Multi/Final_plots/proportions_microglia_by_timepoint.pdf", 
       p1, width = 10, height = 6, dpi = 300)
ggsave("/data/outputs/cellranger_outputs/Analysis_by_Mo_Multi/Final_plots/proportions_microglia_by_genotype.pdf", 
       p2, width = 10, height = 6, dpi = 300)
ggsave("/data/outputs/cellranger_outputs/Analysis_by_Mo_Multi/Final_plots/proportions_microglia_by_condition.pdf", 
       p3, width = 10, height = 6, dpi = 300)
ggsave("/data/outputs/cellranger_outputs/Analysis_by_Mo_Multi/Final_plots/proportions_microglia_by_treatment.pdf", 
       p4, width = 10, height = 6, dpi = 300)

# Print color mapping for reference
cat("Color mapping used:\n")
print(microglia_colors_exotic)




#publication ready area plots

# Load required packages
# Required packages
library(dplyr)
library(tidyr)
library(areaplot)



# Your custom color palette
# Your exotic microglia color palette
microglia_colors_exotic <- c(
  "Antigen-Presenting Microglia"      = "#E63946",  # Flamingo red – urgent immune alert
  "Disease-Associated Microglia"      = "#6A4C93",  # Mystic violet – neurodegenerative vibes
  "Homeostatic Microglia"             = "#2A9D8F",  # Deep jade green – peaceful balance
  "Interferon-Responsive Microglia"   = "#4CC9F0",  # Electric cyan – antiviral spark
  "Proliferating Microglia"           = "#F4A261"   # Spicy mango – hot mitotic energy
)

# Apply to Seurat UMAP
library(Seurat)
DimPlot(PD, 
        group.by = "microglia_type",
        cols = microglia_colors_exotic,
        pt.size = 0.5,
        label = F,
        label.size = 4,
        repel = TRUE) +
  theme_void() +
  theme(legend.position = "right",
        legend.text = element_text(size = 12))

# Store colors in Seurat object for future use
PD@misc$exotic_colors <- microglia_colors_exotic

# Preview your exotic color palette
scales::show_col(microglia_colors_exotic, labels = TRUE)

print("Your exotic microglia colors:")
print(microglia_colors_exotic)



library(dplyr)
library(tidyr)
library(areaplot)

# Define exotic microglia color palette
microglia_colors_exotic <- c(
  "Antigen-Presenting Microglia"      = "#E63946",  # Flamingo red – urgent immune alert
  "Disease-Associated Microglia"      = "#6A4C93",  # Mystic violet – neurodegenerative vibes
  "Homeostatic Microglia"             = "#2A9D8F",  # Deep jade green – peaceful balance
  "Interferon-Responsive Microglia"   = "#4CC9F0",  # Electric cyan – antiviral spark
  "Proliferating Microglia"           = "#F4A261"   # Spicy mango – hot mitotic energy
)

# 1. Prepare metadata and sample_id
meta <- as.data.frame(PD@meta.data)
meta$sample_id <- paste(meta$Condition, meta$Genotype, meta$Timepoint, meta$Treatment, sep = "_")

# ==================== PLOT 1: BY TREATMENT ====================

# 2. Summarize: counts per Treatment × microglia_type
df <- meta %>%
  group_by(Treatment, microglia_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Treatment) %>%
  mutate(freq = n / sum(n)) %>%  # proportional
  ungroup() %>%  # Important: ungroup before complete()
  dplyr::select(Treatment, microglia_type, freq)

# Optional: Check if all treatments have all cell types (fill missing with 0)
df <- df %>%
  complete(Treatment, microglia_type, fill = list(freq = 0))

# 3. Reshape to wide format: rows = Treatment, columns = microglia_type
df_wide <- df %>%
  pivot_wider(names_from = microglia_type, values_from = freq, values_fill = 0) %>%
  arrange(Treatment)

# 4. Prepare data for areaplot
# Convert Treatment to a factor to preserve order
treatment_levels <- unique(meta$Treatment)  # or specify custom order
x <- factor(df_wide$Treatment, levels = treatment_levels)

# Ensure y is numeric matrix/data.frame
y <- df_wide %>%
  dplyr::select(-Treatment) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.data.frame()

# Apply exotic colors to match cell types
cols <- microglia_colors_exotic[colnames(y)]

# Create the area plot
par(mar = c(5, 4, 4, 18), xpd = T)
 # Increase right margin for legend
areaplot(x = x,
         y = y,
         prop = TRUE,
         col = cols,
         border = "white",
         lwd = 1,
         lty = 1,
         legend = FALSE,  # Don't draw legend automatically
         main = "Microglia Subtype Proportions by Treatment",
         xlab = "Treatment",
         ylab = "Proportion",
         cex.main = 1.2,
         cex.lab = 1.1)
abline(h = 0, col = "white", lwd = 4)  # match background

# Add legend outside plot area
legend(x = par("usr")[2] + 0.05 * diff(par("usr")[1:2]), 
       y = par("usr")[4], 
       legend = colnames(y),
       fill = cols,
       border = "white",
       bty = "n",  # No box around legend
       cex = 0.9,
       xjust = 0,
       yjust = 1)

# ==================== PLOT 2: BY TIMEPOINT ====================

# 2. Summarize: counts per Timepoint × microglia_type
df2 <- meta %>%
  group_by(Timepoint, microglia_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Timepoint) %>%
  mutate(freq = n / sum(n)) %>%  # proportional
  ungroup() %>%  # Important: ungroup before complete()
  dplyr::select(Timepoint, microglia_type, freq)

# Optional: Check if all timepoints have all cell types (fill missing with 0)
df2 <- df2 %>%
  complete(Timepoint, microglia_type, fill = list(freq = 0))

# 3. Reshape to wide format: rows = Timepoint, columns = microglia_type
df_wide2 <- df2 %>%
  pivot_wider(names_from = microglia_type, values_from = freq, values_fill = 0) %>%
  arrange(Timepoint)

# 4. Prepare data for areaplot
# Convert Timepoint to a factor to preserve order
timepoint_levels <- unique(meta$Timepoint)  # or specify custom order
x2 <- factor(df_wide2$Timepoint, levels = timepoint_levels)

# Ensure y is numeric matrix/data.frame
y2 <- df_wide2 %>%
  dplyr::select(-Timepoint) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.data.frame()

# Apply exotic colors to match cell types
cols2 <- microglia_colors_exotic[colnames(y2)]

# Create the area plot
par(mar = c(5, 4, 4, 18), xpd = TRUE)  # Increase right margin for legend
areaplot(x = x2,
         y = y2,
         prop = TRUE,
         col = cols2,
         border = "white",
         lwd = 1,
         lty = 1,
         legend = FALSE,  # Don't draw legend automatically
         main = "Microglia Subtype Proportions by Timepoint",
         xlab = "Timepoint",
         ylab = "Proportion",
         cex.main = 1.2,
         cex.lab = 1.1)
abline(h = 0, col = "white", lwd = 4)  # match background

# Add legend outside plot area
legend(x = par("usr")[2] + 0.05 * diff(par("usr")[1:2]), 
       y = par("usr")[4], 
       legend = colnames(y2),
       fill = cols2,
       border = "white",
       bty = "n",  # No box around legend
       cex = 0.9,
       xjust = 0,
       yjust = 1)

# Print color mapping for reference
cat("Color mapping used:\n")
print(microglia_colors_exotic)