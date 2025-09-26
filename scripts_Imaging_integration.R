# Seurat Object and Image Data Integration Workflow
# ================================================

library(Seurat)
library(dplyr)
library(readr)
library(ggplot2)

# Load your data (already done)
# PD <- readRDS("Processed_harmonized_artifacts_removed_annotated_by_Dolan_paper.rds")
# image_data <- your combined image data

# 1. EXPLORE THE DATA STRUCTURE
# -----------------------------
print("=== SEURAT OBJECT OVERVIEW ===")
print(PD)

print("\n=== SAMPLE DISTRIBUTION IN SEURAT ===")
print(table(PD$orig.ident))  # or whatever column contains sample info

print("\n=== IMAGE DATA OVERVIEW ===")
print(head(image_data))
print("\nUnique SampleIDs in image data:")
print(unique(image_data$SampleID))

print("\nImage data dimensions:")
print(dim(image_data))

# 2. IDENTIFY SAMPLE MAPPING AND CREATE MAPPING FUNCTION
# ------------------------------------------------------
# Check sample identifiers in both datasets
seurat_samples <- unique(PD$orig.ident)
image_samples <- unique(image_data$SampleID)

print("\n=== SAMPLE MAPPING CHECK ===")
print("Samples in Seurat:")
print(seurat_samples)
print("\nSamples in Image data:")
print(image_samples)

# Create mapping function to standardize sample names
standardize_sample_name <- function(sample_name) {
  # Convert Seurat format to Image format
  sample_name %>%
    str_replace("GBA1-mut-", "GBA1_") %>%  # Replace GBA1-mut- with GBA1_
    str_replace("KOLF-", "KOLF_") %>%      # Replace KOLF- with KOLF_
    str_replace_all("-", "_")              # Replace remaining hyphens with underscores
}

# Create standardized sample IDs for Seurat data
seurat_samples_standardized <- sapply(seurat_samples, standardize_sample_name)
names(seurat_samples_standardized) <- seurat_samples

print("\n=== STANDARDIZED MAPPING ===")
print("Seurat -> Image format mapping:")
for(i in 1:length(seurat_samples_standardized)) {
  cat(names(seurat_samples_standardized)[i], " -> ", seurat_samples_standardized[i], "\n")
}

# Find overlapping samples after standardization
overlapping_samples <- intersect(seurat_samples_standardized, image_samples)
print(paste("\nOverlapping samples after standardization:", length(overlapping_samples)))
print(overlapping_samples)

# Add standardized sample ID to Seurat metadata
PD$SampleID_standardized <- standardize_sample_name(PD$orig.ident)

# 3. PREPARE IMAGE DATA FOR INTEGRATION
# ------------------------------------
# Your image data appears to have individual cell/nucleus measurements
# We need to aggregate by sample to get sample-level metrics

print("\n=== IMAGE DATA STRUCTURE ===")
print("Image data columns:")
print(colnames(image_data))
print("\nSample distribution in image data:")
print(table(image_data$SampleID))

# Aggregate image data by sample
image_summary <- image_data %>%
  group_by(SampleID) %>%
  summarise(
    # Count metrics
    n_nuclei = n(),
    
    # Speck-related metrics
    total_specks = sum(Speck_Count, na.rm = TRUE),
    mean_specks_per_nucleus = mean(Speck_Count, na.rm = TRUE),
    median_specks_per_nucleus = median(Speck_Count, na.rm = TRUE),
    prop_nuclei_with_specks = mean(Speck_Count > 0, na.rm = TRUE),
    
    # Speck area metrics
    mean_total_speck_area = mean(Total_Speck_Area, na.rm = TRUE),
    median_total_speck_area = median(Total_Speck_Area, na.rm = TRUE),
    
    # Speck intensity metrics
    mean_speck_intensity = mean(Mean_Speck_Intensity[Speck_Count > 0], na.rm = TRUE),
    median_speck_intensity = median(Mean_Speck_Intensity[Speck_Count > 0], na.rm = TRUE),
    
    # IBA1 intensity metrics (microglial marker)
    mean_IBA1_intensity = mean(IBA1_Intensity, na.rm = TRUE),
    median_IBA1_intensity = median(IBA1_Intensity, na.rm = TRUE),
    
    # MKI67 metrics (proliferation marker)
    mean_MKI67_intensity = mean(MKI67_Intensity, na.rm = TRUE),
    median_MKI67_intensity = median(MKI67_Intensity, na.rm = TRUE),
    prop_MKI67_positive = mean(MKI67_Positive, na.rm = TRUE),
    
    .groups = "drop"
  )

print("\n=== AGGREGATED IMAGE DATA ===")
print(head(image_summary))

# 4. CREATE SAMPLE-LEVEL METADATA
# ------------------------------
# Extract sample-level metadata from Seurat object
sample_metadata <- PD@meta.data %>%
  group_by(orig.ident, SampleID_standardized) %>%
  summarise(
    n_cells = n(),
    # Add other sample-level aggregations from scRNA-seq data
    mean_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
    mean_nCount_RNA = mean(nCount_RNA, na.rm = TRUE),
    mean_percent_mt = mean(percent.mt, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(SampleID = SampleID_standardized)

# 5. MERGE SAMPLE METADATA WITH IMAGE DATA
# ----------------------------------------
integrated_sample_data <- sample_metadata %>%
  left_join(image_summary, by = "SampleID")

print("\n=== INTEGRATED SAMPLE-LEVEL DATA ===")
print(integrated_sample_data)

# Check for missing data
missing_image_data <- integrated_sample_data %>%
  filter(is.na(n_nuclei))
if(nrow(missing_image_data) > 0) {
  print("\nSamples missing image data:")
  print(missing_image_data$orig.ident)
} else {
  print("\nAll samples have matching image data!")
}

# 6. ADD IMAGE METRICS TO SEURAT METADATA
# ---------------------------------------
# Create a mapping dataframe for adding to cell-level metadata
image_metrics_for_cells <- integrated_sample_data %>%
  select(-c(orig.ident, n_cells, mean_nFeature_RNA, mean_nCount_RNA, mean_percent_mt))

# Add image metrics to each cell's metadata
# Add image metrics to each cell's metadata
# First check what columns we have
print("Columns in PD@meta.data:")
print(colnames(PD@meta.data))
print("Columns in image_metrics_for_cells:")
print(colnames(image_metrics_for_cells))

# Remove any duplicate columns from image_metrics_for_cells that might already exist in metadata
image_cols_to_add <- setdiff(colnames(image_metrics_for_cells), colnames(PD@meta.data))
image_metrics_clean <- image_metrics_for_cells %>%
  select(SampleID, all_of(image_cols_to_add))

# Add image metrics to each cell's metadata
cell_metadata_with_images <- PD@meta.data %>%
  left_join(image_metrics_clean, by = c("SampleID_standardized" = "SampleID"))

# Update Seurat object metadata
PD@meta.data <- cell_metadata_with_images

print("\n=== UPDATED SEURAT METADATA COLUMNS ===")
print(colnames(PD@meta.data))

# 7. QUALITY CONTROL AND VISUALIZATION
# ------------------------------------
# Visualize the integration
print("\n=== CREATING INTEGRATION PLOTS ===")

# Plot 1: Sample-level correlation between cell count and nuclei count
p1 <- integrated_sample_data %>%
  ggplot(aes(x = n_cells, y = n_nuclei)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "scRNA-seq Cells vs Image Nuclei Count",
       x = "Number of scRNA-seq Cells",
       y = "Number of Nuclei (Image)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

print(p1)

# Plot 2: Speck metrics across conditions
# Properly parse the sample names
plot_data <- integrated_sample_data %>%
  mutate(
    Genotype = case_when(
      str_detect(orig.ident, "GBA1-mut") ~ "GBA1-mut",
      str_detect(orig.ident, "KOLF") ~ "KOLF",
      TRUE ~ "Unknown"
    ),
    Time = case_when(
      str_detect(orig.ident, "24h") ~ "24h",
      str_detect(orig.ident, "7d") ~ "7d",
      TRUE ~ "Unknown"
    ),
    Treatment = case_when(
      str_detect(orig.ident, "aSyn") ~ "aSyn",
      str_detect(orig.ident, "ctrl") ~ "ctrl", 
      str_detect(orig.ident, "debris") ~ "debris",
      str_detect(orig.ident, "LPS") ~ "LPS",
      TRUE ~ "Unknown"
    )
  )

p2 <- plot_data %>%
  ggplot(aes(x = Treatment, y = mean_specks_per_nucleus, fill = Genotype)) +
  geom_bar(stat = "identity", position = "dodge") +  # <- bar plot
  facet_wrap(~Time) +
  labs(title = "Mean Specks per Nucleus by Treatment",
       y = "Mean Specks per Nucleus",
       x = "Treatment") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p2)

# Plot 3: UMAP colored by image metrics
# First check which image metrics were successfully added
print("Available image metrics in Seurat metadata:")
image_metric_cols <- colnames(PD@meta.data)[grepl("mean_|median_|prop_|total_", colnames(PD@meta.data))]
print(image_metric_cols)

# Check for non-NA values
non_na_metrics <- sapply(image_metric_cols, function(x) sum(!is.na(PD@meta.data[[x]])))
print("Non-NA counts for image metrics:")
print(non_na_metrics)

# Use a metric that actually has values
if("mean_specks_per_nucleus" %in% colnames(PD@meta.data) && sum(!is.na(PD@meta.data$mean_specks_per_nucleus)) > 0) {
  p3 <- FeaturePlot(PD, features = "mean_specks_per_nucleus", 
                    cols = c("lightgrey", "red")) +
    labs(title = "UMAP: Mean Specks per Nucleus (Sample-level)")
} else if("mean_IBA1_intensity" %in% colnames(PD@meta.data) && sum(!is.na(PD@meta.data$mean_IBA1_intensity)) > 0) {
  p3 <- FeaturePlot(PD, features = "mean_IBA1_intensity", 
                    cols = c("lightgrey", "blue")) +
    labs(title = "UMAP: Mean IBA1 Intensity (Sample-level)")
} else {
  # If no image metrics available, use a basic cell count plot
  p3 <- DimPlot(PD, group.by = "orig.ident") +
    labs(title = "UMAP: Samples") +
    theme(legend.position = "right")
}

print(p3)

# Plot 4: IBA1 intensity correlation with microglial markers
p4_bar <- plot_data %>%
  ggplot(aes(x = Treatment, y = mean_IBA1_intensity, fill = Genotype)) +
  geom_bar(stat = "identity", position = "dodge") +  # <- bar plot
  facet_wrap(~Time) +
  labs(title = "IBA1 Intensity by Treatment",
       y = "Mean IBA1 Intensity",
       x = "Treatment") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        legend.position = "top")

print(p4_bar)

# 8. SAVE INTEGRATED DATA
# ----------------------
# Save updated Seurat object
saveRDS(PD, "PD_with_image_integration.rds")

# Save sample-level integrated data
write_csv(integrated_sample_data, "integrated_sample_metadata.csv")

print("\n=== INTEGRATION COMPLETE ===")
print("Updated Seurat object saved as: PD_with_image_integration.rds")
print("Sample-level data saved as: integrated_sample_metadata.csv")

# 9. VERIFICATION STEPS
# ---------------------
print("\n=== VERIFICATION ===")
print("Samples with both scRNA-seq and image data:")
complete_samples <- integrated_sample_data %>%
  filter(!is.na(n_nuclei)) %>%
  pull(orig.ident)
print(complete_samples)
print(paste("Complete samples:", length(complete_samples), "out of", length(seurat_samples)))

# Check data distribution
print("\nKey image metrics summary:")
print("Specks per nucleus:")
print(summary(integrated_sample_data$mean_specks_per_nucleus))
print("\nIBA1 intensity:")
print(summary(integrated_sample_data$mean_IBA1_intensity))
print("\nMKI67 positivity rate:")
print(summary(integrated_sample_data$prop_MKI67_positive))



DimPlot(PD)

FeatureScatter(PD, feature1 = "MKI67", feature2 = "prop_MKI67_positive")

table(PD@meta.data$prop_MKI67_positive)
















# scRNA-seq and Image Data Correlation Analysis
# =============================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(corrplot)
library(viridis)
library(patchwork)
library(stringr)

# 1. DEFINE GENES OF INTEREST
# ===========================
Homeostatic_Microglia = c("CX3CR1", "P2RY12", "TMEM119", "SLC2A5", "HEXB", "CSF1R", "GPR34", "MEF2A")

Disease_Associated_Microglia = c("APOE", "TREM2", "GPNMB", "LPL", "ABCA1", "CD9", "CTSD", "CTSB", "SPP1", "CLEC7A")

Antigen_Presenting_Microglia = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "CD74", "CIITA", "B2M")

Interferon_Responsive_Microglia = c("IFIT1", "IFIT2", "IFIT3", "MX1", "MX2", "ISG15", "STAT1", "IRF7", "OAS1", "OAS2")

Proliferating_Microglia = c("MKI67", "TOP2A", "PCNA", "STMN1", "HIST1H4C", "HIST1H1B")


# Combine all genes
all_genes_of_interest <- c(Homeostatic_Microglia, Disease_Associated_Microglia, Antigen_Presenting_Microglia, 
                           Interferon_Responsive_Microglia, Proliferating_Microglia)

# Check which genes are available in your dataset
available_genes <- intersect(all_genes_of_interest, rownames(PD))
print("=== AVAILABLE GENES OF INTEREST ===")
print(available_genes)
print(paste("Found", length(available_genes), "out of", length(all_genes_of_interest), "genes"))

# 2. CALCULATE SAMPLE-LEVEL GENE EXPRESSION
# =========================================
# Get average expression per sample for genes of interest
sample_expression <- AverageExpression(PD, features = available_genes, 
                                       group.by = "orig.ident", 
                                       assays = "SCT")$SCT

# Convert to dataframe and add sample info
sample_expression_df <- sample_expression %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "orig.ident", values_to = "Expression") %>%
  pivot_wider(names_from = Gene, values_from = Expression)

print("=== SAMPLE-LEVEL EXPRESSION DATA ===")
print(head(sample_expression_df))

# 3. MERGE WITH IMAGE DATA
# ========================
# Get the integrated sample data from previous analysis
correlation_data <- integrated_sample_data %>%
  left_join(sample_expression_df, by = "orig.ident")

print("=== CORRELATION DATA DIMENSIONS ===")
print(dim(correlation_data))
print(colnames(correlation_data))

# 4. CORRELATION ANALYSIS - IMAGE METRICS vs GENE EXPRESSION
# ==========================================================

# Define image metrics of interest
image_metrics <- c("mean_specks_per_nucleus", "prop_nuclei_with_specks", 
                   "mean_IBA1_intensity", "prop_MKI67_positive", 
                   "mean_total_speck_area", "mean_speck_intensity")

# Function to calculate correlations
calculate_correlations <- function(data, genes, metrics) {
  correlations <- list()
  p_values <- list()
  
  for(gene in genes) {
    if(gene %in% colnames(data)) {
      gene_correlations <- numeric()
      gene_p_values <- numeric()
      
      for(metric in metrics) {
        if(metric %in% colnames(data)) {
          cor_test <- cor.test(data[[gene]], data[[metric]], method = "pearson")
          gene_correlations[metric] <- cor_test$estimate
          gene_p_values[metric] <- cor_test$p.value
        }
      }
      correlations[[gene]] <- gene_correlations
      p_values[[gene]] <- gene_p_values
    }
  }
  
  # Convert to matrices
  cor_matrix <- do.call(rbind, correlations)
  p_matrix <- do.call(rbind, p_values)
  
  return(list(correlations = cor_matrix, p_values = p_matrix))
}

# Calculate correlations
cor_results <- calculate_correlations(correlation_data, available_genes, image_metrics)

print("=== CORRELATION MATRIX ===")
print(cor_results$correlations)

# 5. VISUALIZATION - CORRELATION HEATMAP
# ======================================
# Create significance mask
sig_mask <- cor_results$p_values < 0.05


# 6. TOP CORRELATIONS ANALYSIS
# ============================
# Find strongest correlations
cor_long <- cor_results$correlations %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Image_Metric", values_to = "Correlation") %>%
  left_join(
    cor_results$p_values %>%
      as.data.frame() %>%
      rownames_to_column("Gene") %>%
      pivot_longer(-Gene, names_to = "Image_Metric", values_to = "P_value"),
    by = c("Gene", "Image_Metric")
  ) %>%
  filter(!is.na(Correlation)) %>%
  mutate(Significant = P_value < 0.05,
         Abs_Correlation = abs(Correlation))

# Top positive and negative correlations
top_correlations <- cor_long %>%
  filter(Significant) %>%
  arrange(desc(Abs_Correlation)) %>%
  head(20)

print("=== TOP SIGNIFICANT CORRELATIONS ===")
print(top_correlations)

# 7. SPECIFIC GENE-METRIC SCATTER PLOTS
# =====================================
# Function to create scatter plots for top correlations
create_scatter_plots <- function(data, gene, metric, title_suffix = "") {
  p <- ggplot(data, aes_string(x = gene, y = metric)) +
    geom_point(size = 3, alpha = 0.7, color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "red", alpha = 0.3) +
    labs(title = paste(gene, "vs", metric, title_suffix),
         x = paste(gene, "Expression (log2)"),
         y = gsub("_", " ", str_to_title(metric))) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12))
  
  # Add correlation info
  cor_test <- cor.test(data[[gene]], data[[metric]])
  p <- p + annotate("text", x = Inf, y = Inf, 
                    label = paste("r =", round(cor_test$estimate, 3), 
                                  "\np =", round(cor_test$p.value, 4)),
                    hjust = 1.1, vjust = 1.1, size = 3)
  return(p)
}

# Create plots for top 6 correlations
plot_list <- list()
for(i in 1:min(6, nrow(top_correlations))) {
  gene <- top_correlations$Gene[i]
  metric <- top_correlations$Image_Metric[i]
  
  if(gene %in% colnames(correlation_data) & metric %in% colnames(correlation_data)) {
    plot_list[[i]] <- create_scatter_plots(correlation_data, gene, metric)
  }
}

# Combine plots
combined_plots <- wrap_plots(plot_list, ncol = 2)
ggsave("top_correlations_scatter.png", combined_plots, width = 12, height = 18)

# 8. CONDITION-SPECIFIC ANALYSIS
# ==============================
# Add condition information
correlation_data_annotated <- correlation_data %>%
  separate(orig.ident, into = c("Genotype", "Timepoint", "Treatment"), 
           sep = "-", remove = FALSE) %>%
  mutate(Genotype = str_replace(Genotype, "mut", "GBA1-mut"))

# Analyze correlations by condition
condition_correlations <- function(data, condition_col) {
  conditions <- unique(data[[condition_col]])
  condition_results <- list()
  
  for(cond in conditions) {
    cond_data <- data[data[[condition_col]] == cond, ]
    if(nrow(cond_data) >= 3) {  # Need at least 3 points for correlation
      cond_results[[cond]] <- calculate_correlations(cond_data, available_genes, image_metrics)
    }
  }
  return(condition_results)
}

# Analyze by genotype
genotype_cors <- condition_correlations(correlation_data_annotated, "Genotype")

print("=== GENOTYPE-SPECIFIC CORRELATIONS ===")
for(genotype in names(genotype_cors)) {
  cat("\n", genotype, ":\n")
  if(!is.null(genotype_cors[[genotype]]$correlations)) {
    print(head(genotype_cors[[genotype]]$correlations))
  }
}


# 10. SAVE RESULTS
# ===============
# Save correlation results
write.csv(cor_long, "gene_image_correlations.csv", row.names = FALSE)
write.csv(correlation_data, "integrated_correlation_data.csv", row.names = FALSE)
write.csv(top_correlations, "top_correlations.csv", row.names = FALSE)


library(corrplot)

# Create significance mask
sig_matrix <- cor_results$p_values < 0.05

library(pheatmap)

# Open a device first
png("correlation_heatmap.png", width = 1200, height = 1500, res = 150)

# Plot the heatmap
pheatmap(cor_results$correlations,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         fontsize_number = 2,
         main = "Gene Expression vs Image Metrics Correlations")

# Close the device to write the file
dev.off()



# Correlation Analysis: Imaging Features vs Cell Type Proportions
library(corrplot)
library(ggplot2)
library(reshape2)
library(Hmisc)

print("=== CORRELATION ANALYSIS APPROACH ===")

# 1. Calculate correlations with p-values
correlation_results <- data.frame(
  CellType = character(),
  ImagingFeature = character(), 
  Correlation = numeric(),
  P_value = numeric(),
  Significance = character(),
  stringsAsFactors = FALSE
)

# Get cell type proportions and imaging features
celltype_data <- merged_data[, celltype_cols]
imaging_data <- merged_data[, colnames(X)]

# Calculate correlations
for(celltype in celltype_cols) {
  for(feature in colnames(imaging_data)) {
    
    # Remove samples with missing data
    complete_idx <- !is.na(celltype_data[[celltype]]) & !is.na(imaging_data[[feature]])
    
    if(sum(complete_idx) >= 3) {  # Need at least 3 points for correlation
      cor_test <- cor.test(celltype_data[[celltype]][complete_idx], 
                           imaging_data[[feature]][complete_idx])
      
      # Determine significance
      if(cor_test$p.value < 0.001) {
        sig_level <- "***"
      } else if(cor_test$p.value < 0.01) {
        sig_level <- "**" 
      } else if(cor_test$p.value < 0.05) {
        sig_level <- "*"
      } else {
        sig_level <- "ns"
      }
      
      correlation_results <- rbind(correlation_results, data.frame(
        CellType = celltype,
        ImagingFeature = feature,
        Correlation = cor_test$estimate,
        P_value = cor_test$p.value,
        Significance = sig_level,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# 2. Create correlation matrix
cor_matrix <- matrix(NA, nrow = length(celltype_cols), ncol = ncol(imaging_data))
rownames(cor_matrix) <- celltype_cols
colnames(cor_matrix) <- colnames(imaging_data)
p_matrix <- cor_matrix

for(i in 1:nrow(correlation_results)) {
  row_idx <- which(rownames(cor_matrix) == correlation_results$CellType[i])
  col_idx <- which(colnames(cor_matrix) == correlation_results$ImagingFeature[i])
  cor_matrix[row_idx, col_idx] <- correlation_results$Correlation[i]
  p_matrix[row_idx, col_idx] <- correlation_results$P_value[i]
}

# 3. Visualize correlations
png("correlation_heatmap_celltype.png", width = 10, height = 8, units = "in", res = 300)
corrplot(cor_matrix, 
         method = "color",
         type = "full",
         order = "hclust",
         tl.cex = 0.8,
         tl.col = "black",
         cl.cex = 0.8,
         addCoef.col = "black",
         number.cex = 0.7,
         title = "Correlations: Cell Types vs Imaging Features",
         mar = c(0,0,1,0))
dev.off()

# 4. Summary table of significant correlations
significant_cors <- correlation_results[correlation_results$P_value < 0.05, ]
significant_cors <- significant_cors[order(abs(significant_cors$Correlation), decreasing = TRUE), ]

print("=== SIGNIFICANT CORRELATIONS (p < 0.05) ===")
print(significant_cors)

# 5. Create scatter plots for top correlations
top_correlations <- head(significant_cors, 6)

plot_list <- list()
for(i in 1:min(6, nrow(top_correlations))) {
  celltype <- top_correlations$CellType[i]
  feature <- top_correlations$ImagingFeature[i]
  cor_val <- round(top_correlations$Correlation[i], 3)
  p_val <- top_correlations$P_value[i]
  
  plot_data <- data.frame(
    Imaging = imaging_data[[feature]],
    CellType = celltype_data[[celltype]]
  )
  
  p <- ggplot(plot_data, aes(x = Imaging, y = CellType)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    labs(title = paste(celltype, "vs", feature),
         subtitle = paste("r =", cor_val, ", p =", round(p_val, 4)),
         x = feature,
         y = paste(celltype, "Proportion")) +
    theme_minimal()
  
  plot_list[[i]] <- p
}

# Save scatter plots
if(length(plot_list) > 0) {
  combined_scatter <- wrap_plots(plot_list, ncol = 3)
  ggsave("correlation_scatter_plots.png", combined_scatter, 
         width = 15, height = 10, dpi = 300)
}

# 6. Effect size interpretation
correlation_results$Effect_Size <- ifelse(abs(correlation_results$Correlation) >= 0.7, "Large",
                                          ifelse(abs(correlation_results$Correlation) >= 0.5, "Medium",
                                                 ifelse(abs(correlation_results$Correlation) >= 0.3, "Small", "Negligible")))

effect_summary <- table(correlation_results$Significance, correlation_results$Effect_Size)
print("=== EFFECT SIZE SUMMARY ===")
print(effect_summary)

# 7. Key findings summary
print("\n=== KEY FINDINGS ===")
n_significant <- sum(correlation_results$P_value < 0.05, na.rm = TRUE)
n_total <- nrow(correlation_results)
cat("Significant correlations:", n_significant, "out of", n_total, "\n")

strongest_positive <- correlation_results[which.max(correlation_results$Correlation), ]
strongest_negative <- correlation_results[which.min(correlation_results$Correlation), ]

cat("Strongest positive correlation:", strongest_positive$CellType, "~", 
    strongest_positive$ImagingFeature, "(r =", round(strongest_positive$Correlation, 3), ")\n")
cat("Strongest negative correlation:", strongest_negative$CellType, "~", 
    strongest_negative$ImagingFeature, "(r =", round(strongest_negative$Correlation, 3), ")\n")

# Save results
write.csv(correlation_results, "celltype_imaging_correlations.csv", row.names = FALSE)
write.csv(significant_cors, "significant_correlations.csv", row.names = FALSE)

print("\n=== ANALYSIS COMPLETE ===")
print("This approach is appropriate for small sample sizes and provides:")
print("- Statistical significance testing")
print("- Effect size estimation") 
print("- No overfitting issues")
print("- Clear biological interpretations")




# Correlation Analysis: Imaging Features vs Metadata Variables

library(corrplot)
library(ggplot2)
library(reshape2)
library(patchwork)
library(dplyr)

print("=== METADATA CORRELATION ANALYSIS ===")

# 1. PREPARE METADATA VARIABLES
# ==============================

# Derive metadata from orig.ident
integrated_sample_data <- integrated_sample_data %>%
  mutate(
    Genotype = case_when(
      grepl("GBA1", orig.ident) ~ "GBA1-mut",
      grepl("KOLF", orig.ident) ~ "KOLF"
    ),
    Timepoint = case_when(
      grepl("24h", orig.ident) ~ "24h",
      grepl("7d", orig.ident) ~ "7d"
    ),
    Treatment = case_when(
      grepl("LPS", orig.ident) ~ "LPS",
      grepl("aSyn", orig.ident) ~ "aSyn",
      grepl("ctrl", orig.ident) ~ "ctrl",
      grepl("debris", orig.ident) ~ "debris"
    )
  )

# Check available columns
print("Available metadata columns:")
print(colnames(integrated_sample_data))

# Define imaging features for X matrix
imaging_features <- c(
  "mean_nFeature_RNA", "mean_nCount_RNA", "mean_percent_mt",
  "n_nuclei", "total_specks", "mean_specks_per_nucleus",
  "median_specks_per_nucleus", "prop_nuclei_with_specks",
  "mean_total_speck_area", "median_total_speck_area",
  "mean_speck_intensity", "median_speck_intensity",
  "mean_IBA1_intensity", "median_IBA1_intensity",
  "mean_MKI67_intensity", "median_MKI67_intensity",
  "prop_MKI67_positive"
)
X <- as.matrix(integrated_sample_data[, imaging_features])
print("X dimensions:")
print(dim(X))

# Extract metadata of interest
metadata_vars <- c("Genotype", "Timepoint", "Treatment")
available_metadata <- intersect(metadata_vars, colnames(integrated_sample_data))
print(paste("Available metadata variables:", paste(available_metadata, collapse = ", ")))

# Prepare analysis data
metadata_analysis_data <- integrated_sample_data %>%
  select(orig.ident, all_of(available_metadata), all_of(colnames(X))) %>%
  na.omit()
print("Rows in metadata_analysis_data:")
print(nrow(metadata_analysis_data))
print("Missing values in metadata_analysis_data:")
print(colSums(is.na(metadata_analysis_data)))

# 2. CONVERT CATEGORICAL TO NUMERIC
# ================================

convert_categorical <- function(data, var_name) {
  if(var_name %in% colnames(data)) {
    if(is.factor(data[[var_name]]) || is.character(data[[var_name]])) {
      unique_vals <- unique(data[[var_name]])
      if(length(unique_vals) > 0) {
        print(paste("Converting", var_name, "levels:", paste(unique_vals, collapse = ", ")))
        data[[paste0(var_name, "_numeric")]] <- as.numeric(as.factor(data[[var_name]]))
      } else {
        warning(paste("No unique values for", var_name, "- skipping conversion"))
      }
    } else {
      print(paste(var_name, "is already numeric or not suitable for conversion"))
    }
  } else {
    warning(paste(var_name, "not found in data"))
  }
  return(data)
}

# Convert categorical variables
for(var in available_metadata) {
  metadata_analysis_data <- convert_categorical(metadata_analysis_data, var)
}

numeric_metadata_cols <- grep("_numeric$", colnames(metadata_analysis_data), value = TRUE)
print(paste("Numeric metadata columns:", paste(numeric_metadata_cols, collapse = ", ")))

# 3. CORRELATION ANALYSIS
# ======================

perform_metadata_correlation <- function(metadata_col, data) {
  correlation_results <- data.frame(
    MetadataVariable = character(),
    ImagingFeature = character(),
    Correlation = numeric(),
    P_value = numeric(),
    Significance = character(),
    stringsAsFactors = FALSE
  )
  
  imaging_features <- colnames(X)
  print(paste("Processing metadata:", metadata_col))
  
  for(feature in imaging_features) {
    complete_idx <- !is.na(data[[metadata_col]]) & !is.na(data[[feature]])
    print(paste("Feature:", feature, "Valid samples:", sum(complete_idx)))
    
    if(sum(complete_idx) >= 3) {
      if(is.numeric(data[[metadata_col]]) && is.numeric(data[[feature]])) {
        cor_test <- tryCatch(
          {
            cor.test(data[[metadata_col]][complete_idx], 
                     data[[feature]][complete_idx], method = "pearson")
          },
          error = function(e) {
            warning(paste("cor.test failed for", metadata_col, "and", feature, ":", e$message))
            return(NULL)
          }
        )
        
        if(!is.null(cor_test)) {
          sig_level <- if(cor_test$p.value < 0.001) "***" else 
            if(cor_test$p.value < 0.01) "**" else 
              if(cor_test$p.value < 0.05) "*" else "ns"
          
          correlation_results <- rbind(correlation_results, data.frame(
            MetadataVariable = metadata_col,
            ImagingFeature = feature,
            Correlation = as.numeric(cor_test$estimate),
            P_value = as.numeric(cor_test$p.value),
            Significance = sig_level,
            stringsAsFactors = FALSE
          ))
        }
      } else {
        warning(paste("Non-numeric data for", metadata_col, "(type:", typeof(data[[metadata_col]]), ") or", feature, "(type:", typeof(data[[feature]]), ")"))
      }
    } else {
      warning(paste("Insufficient data for", metadata_col, "and", feature, ":", sum(complete_idx), "samples"))
    }
  }
  
  print("Correlation results:")
  print(correlation_results)
  return(correlation_results)
}

# 4. ANALYZE EACH METADATA VARIABLE
# ================================

all_metadata_results <- data.frame()
correlation_matrices <- list()

analysis_vars <- numeric_metadata_cols  # Use only numeric metadata
for(meta_var in analysis_vars) {
  if(meta_var %in% colnames(metadata_analysis_data)) {
    cat("\n=== Analyzing", meta_var, "===\n")
    
    results <- perform_metadata_correlation(meta_var, metadata_analysis_data)
    
    if(nrow(results) > 0) {
      results$MetadataVariable <- meta_var
      all_metadata_results <- rbind(all_metadata_results, results)
      
      cor_vec <- results$Correlation
      names(cor_vec) <- results$ImagingFeature
      correlation_matrices[[meta_var]] <- cor_vec
      
      significant_results <- results[results$P_value < 0.05, ]
      if(nrow(significant_results) > 0) {
        cat("Significant correlations for", meta_var, ":\n")
        print(significant_results[order(abs(significant_results$Correlation), decreasing = TRUE), ])
      } else {
        cat("No significant correlations found for", meta_var, "\n")
      }
    } else {
      cat("No correlation results for", meta_var, "\n")
    }
  } else {
    warning(paste(meta_var, "not found in metadata_analysis_data"))
  }
}

# 5. VISUALIZATION
# ===============

if(length(correlation_matrices) > 0) {
  imaging_features <- colnames(X)
  combined_cor_matrix <- matrix(NA, 
                                nrow = length(analysis_vars), 
                                ncol = length(imaging_features))
  rownames(combined_cor_matrix) <- analysis_vars
  colnames(combined_cor_matrix) <- imaging_features
  
  for(i in 1:length(analysis_vars)) {
    meta_var <- analysis_vars[i]
    if(meta_var %in% names(correlation_matrices)) {
      cor_vec <- correlation_matrices[[meta_var]]
      combined_cor_matrix[i, names(cor_vec)] <- cor_vec
    }
  }
  
  png("metadata_correlation_heatmap.png", width = 12, height = 6, units = "in", res = 300)
  corrplot(combined_cor_matrix,
           method = "color",
           type = "full",
           tl.cex = 0.8,
           tl.col = "black",
           cl.cex = 0.8,
           addCoef.col = "black",
           number.cex = 0.7,
           title = "Correlations: Metadata vs Imaging Features",
           mar = c(0,0,2,0))
  dev.off()
} else {
  cat("No correlation matrices available for heatmap.\n")
}

# 6. SCATTER PLOTS FOR SIGNIFICANT CORRELATIONS
# ============================================

significant_metadata_cors <- all_metadata_results[all_metadata_results$P_value < 0.05, ]

print("\n=== ALL SIGNIFICANT METADATA CORRELATIONS ===")
if(nrow(all_metadata_results) == 0) {
  cat("No correlations computed. Check input data or metadata variables.\n")
} else if(nrow(significant_metadata_cors) == 0 || !("Correlation" %in% colnames(significant_metadata_cors))) {
  cat("No significant correlations found or Correlation column missing.\n")
} else {
  significant_metadata_cors$Correlation <- as.numeric(significant_metadata_cors$Correlation)
  significant_metadata_cors <- significant_metadata_cors[!is.na(significant_metadata_cors$Correlation), ]
  significant_metadata_cors <- significant_metadata_cors[order(abs(significant_metadata_cors$Correlation), decreasing = TRUE), ]
  print(significant_metadata_cors)
  
  top_metadata_cors <- head(significant_metadata_cors, 9)
  
  plot_list <- list()
  for(i in 1:nrow(top_metadata_cors)) {
    meta_var <- top_metadata_cors$MetadataVariable[i]
    feature <- top_metadata_cors$ImagingFeature[i]
    cor_val <- round(top_metadata_cors$Correlation[i], 3)
    p_val <- top_metadata_cors$P_value[i]
    
    plot_data <- data.frame(
      Imaging = metadata_analysis_data[[feature]],
      Metadata = metadata_analysis_data[[meta_var]]
    )
    
    p <- ggplot(plot_data, aes(x = Imaging, y = Metadata)) +
      geom_point(size = 2, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, color = "red") +
      labs(title = paste(meta_var, "vs", feature),
           subtitle = paste("r =", cor_val, ", p =", round(p_val, 4)),
           x = feature,
           y = meta_var) +
      theme_minimal() +
      theme(plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8))
    
    plot_list[[i]] <- p
  }
  
  if(length(plot_list) > 0) {
    combined_metadata_scatter <- wrap_plots(plot_list, ncol = 3)
    ggsave("metadata_correlation_scatter_plots.png", combined_metadata_scatter, 
           width = 15, height = 12, dpi = 300)
  } else {
    cat("No scatter plots generated: No significant correlations.\n")
  }
}

# 7. SUMMARY STATISTICS
# ====================

print("\n=== METADATA CORRELATION SUMMARY ===")
if(nrow(all_metadata_results) > 0) {
  summary_stats <- all_metadata_results %>%
    group_by(MetadataVariable) %>%
    summarise(
      Total_Correlations = n(),
      Significant_Correlations = sum(P_value < 0.05, na.rm = TRUE),
      Mean_Abs_Correlation = mean(abs(Correlation), na.rm = TRUE),
      Max_Abs_Correlation = max(abs(Correlation), na.rm = TRUE),
      .groups = "drop"
    )
  print(summary_stats)
} else {
  cat("No correlation results available for summary.\n")
}

# 8. SAVE RESULTS
# ===============

write.csv(all_metadata_results, "metadata_imaging_correlations.csv", row.names = FALSE)
if(nrow(significant_metadata_cors) > 0 && "Correlation" %in% colnames(significant_metadata_cors)) {
  write.csv(significant_metadata_cors, "significant_metadata_correlations.csv", row.names = FALSE)
}
if(exists("summary_stats")) {
  write.csv(summary_stats, "metadata_correlation_summary.csv", row.names = FALSE)
}

print("\n=== METADATA ANALYSIS COMPLETE ===")
print("Files saved:")
print("- metadata_correlation_heatmap.png")
if(nrow(significant_metadata_cors) > 0 && "Correlation" %in% colnames(significant_metadata_cors)) {
  print("- metadata_correlation_scatter_plots.png")
  print("- significant_metadata_correlations.csv")
}
if(exists("summary_stats")) {
  print("- metadata_correlation_summary.csv")
}
print("- metadata_imaging_correlations.csv")

if(nrow(significant_metadata_cors) > 0 && "Correlation" %in% colnames(significant_metadata_cors)) {
  cat("\nKey findings:\n")
  cat("Variables with strongest imaging associations:\n")
  top_vars <- summary_stats[order(-summary_stats$Max_Abs_Correlation), ]
  print(head(top_vars, 3))
} else {
  cat("\nNo significant correlations found between metadata and imaging features.\n")
}

