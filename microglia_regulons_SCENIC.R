# ============================================================================
# SCENIC ANALYSIS WORKFLOW
# ============================================================================

# ----------------------------------------------------------------------------
# MODULE 1: LOAD LIBRARIES
# ----------------------------------------------------------------------------

library(Seurat)
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(GENIE3)
library(pheatmap)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(dplyr)


# ----------------------------------------------------------------------------
# MODULE 2: PREPARE EXPRESSION MATRIX
# ----------------------------------------------------------------------------

# Extract expression matrix from Seurat object
exprMat <- as.matrix(GetAssayData(seu_ds, slot = "data"))

# Filter genes (expressed in at least 1% of cells)
min_cells <- ncol(exprMat) * 0.01
genes_keep <- rowSums(exprMat > 0) >= min_cells
exprMat_filtered <- exprMat[genes_keep, ]

print(paste("Genes kept:", nrow(exprMat_filtered)))


# ----------------------------------------------------------------------------
# MODULE 3: PREPARE TRANSCRIPTION FACTORS
# ----------------------------------------------------------------------------

# Download TF list
download.file(
  "https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/hs_hgnc_tfs.txt",
  "hs_tfs.txt"
)

# Load TF list
tfs <- read.table("hs_tfs.txt", header = FALSE)$V1

# Keep only TFs present in your data
tfs_in_data <- intersect(tfs, rownames(exprMat_filtered))
print(paste("TFs found in data:", length(tfs_in_data)))


# ----------------------------------------------------------------------------
# MODULE 4: GENE REGULATORY NETWORK INFERENCE (GENIE3)
# ----------------------------------------------------------------------------

# Run GENIE3 (may take 30 minutes to several hours)
set.seed(123)
weightMat <- GENIE3(exprMat_filtered, 
                    regulators = tfs_in_data,
                    nCores = 4)

# Save the weight matrix
saveRDS(weightMat, "genie3_weights.rds")

# Get regulatory links sorted by importance
linkList <- getLinkList(weightMat, threshold = 0.001)

# View top interactions
head(linkList, 20)

# Save results
write.csv(linkList, "regulatory_network.csv", row.names = FALSE)


# ----------------------------------------------------------------------------
# MODULE 5: BUILD REGULONS
# ----------------------------------------------------------------------------

# Initialize regulons list
regulons <- list()

# For each TF, get its top target genes
for(tf in tfs_in_data) {
  # Get targets for this TF
  tf_links <- linkList[linkList$regulatoryGene == tf, ]
  
  # Keep top 50 targets (or use a weight threshold)
  top_targets <- head(tf_links[order(-tf_links$weight), ], 50)
  
  if(nrow(top_targets) > 0) {
    regulons[[tf]] <- top_targets$targetGene
  }
}

# Check regulon sizes
regulon_sizes <- sapply(regulons, length)
print(head(sort(regulon_sizes, decreasing = TRUE), 20))

# Save regulons
saveRDS(regulons, "regulons.rds")


# ----------------------------------------------------------------------------
# MODULE 6: CALCULATE REGULON ACTIVITY (AUCell)
# ----------------------------------------------------------------------------

# Build gene sets for AUCell
geneSets <- regulons

# Calculate cell rankings
cells_rankings <- AUCell_buildRankings(exprMat_filtered, 
                                       nCores = 4,
                                       plotStats = FALSE)

# Calculate AUC scores
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)

# Get AUC matrix
auc_matrix <- getAUC(cells_AUC)

# View results
dim(auc_matrix)
auc_matrix[1:5, 1:5]


# ----------------------------------------------------------------------------
# MODULE 7: ADD SCENIC RESULTS TO SEURAT OBJECT
# ----------------------------------------------------------------------------

# Add AUC scores as a new assay
seu_ds[["SCENIC"]] <- CreateAssayObject(data = auc_matrix)


# ----------------------------------------------------------------------------
# MODULE 8: VISUALIZE REGULONS ON UMAP
# ----------------------------------------------------------------------------

# Set default assay
DefaultAssay(seu_ds) <- "SCENIC"

# Find most variable regulons
regulon_variance <- apply(auc_matrix, 1, var)
top_regulons <- names(sort(regulon_variance, decreasing = TRUE)[1:9])

# Plot top 9 variable regulons
FeaturePlot(seu_ds, 
            features = top_regulons,
            ncol = 3)

# UMAP with specific regulon (example)
FeaturePlot(seu_ds, 
            features = "MAFB",  # replace with your TF of interest
            min.cutoff = "q10",
            max.cutoff = "q90")


# ----------------------------------------------------------------------------
# MODULE 9: BASIC HEATMAP - TOP VARIABLE REGULONS
# ----------------------------------------------------------------------------

# Get top 20 variable regulons
regulon_variance <- apply(auc_matrix, 1, var)
top_regulons_20 <- names(sort(regulon_variance, decreasing = TRUE)[1:20])

# Subset AUC matrix
auc_subset <- auc_matrix[top_regulons_20, ]

# Basic heatmap
pheatmap(auc_subset,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         show_colnames = FALSE,
         color = viridis(100),
         main = "Top 20 Variable Regulons")


# ----------------------------------------------------------------------------
# MODULE 10: HEATMAP WITH CELL TYPE ANNOTATIONS
# ----------------------------------------------------------------------------

# Get cell type annotations from Seurat
cell_types <- seu_ds$microglia_type  # replace with your actual column name
names(cell_types) <- colnames(seu_ds)

# Create annotation dataframe
annotation_col <- data.frame(
  CellType = cell_types
)
rownames(annotation_col) <- names(cell_types)

# Heatmap with cell type annotations
pheatmap(auc_subset,
         scale = "row",
         annotation_col = annotation_col,
         show_colnames = FALSE,
         clustering_distance_rows = "correlation",
         color = viridis(100),
         main = "Regulon Activity by Cell Type")


# ----------------------------------------------------------------------------
# MODULE 11: AVERAGE REGULON ACTIVITY PER CELL TYPE
# ----------------------------------------------------------------------------

# Calculate mean AUC per cell type
auc_by_celltype <- matrix(nrow = nrow(auc_matrix), 
                          ncol = length(unique(cell_types)))
rownames(auc_by_celltype) <- rownames(auc_matrix)
colnames(auc_by_celltype) <- unique(cell_types)

for(ct in unique(cell_types)) {
  cells_in_ct <- names(cell_types)[cell_types == ct]
  cells_in_ct <- intersect(cells_in_ct, colnames(auc_matrix))
  auc_by_celltype[, ct] <- rowMeans(auc_matrix[, cells_in_ct, drop=FALSE])
}

# Get top 30 variable regulons
top_var_30 <- names(sort(apply(auc_by_celltype, 1, var), decreasing = TRUE)[1:30])

# Plot average regulon activity
pheatmap(auc_by_celltype[top_var_30, ],
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Average Regulon Activity per Cell Type",
         fontsize_row = 8)


# ----------------------------------------------------------------------------
# MODULE 12: COMPLEXHEATMAP VISUALIZATION
# ----------------------------------------------------------------------------

# Set color scale
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Create cell type colors
celltype_colors <- rainbow(length(unique(cell_types)))
names(celltype_colors) <- unique(cell_types)

# Create top annotation
ha_top <- HeatmapAnnotation(
  CellType = cell_types,
  col = list(CellType = celltype_colors),
  show_legend = TRUE
)

# Create ComplexHeatmap
Heatmap(auc_subset,
        name = "AUC",
        col = col_fun,
        top_annotation = ha_top,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 8),
        column_title = "Top 20 Variable Regulons",
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "pearson")


# ----------------------------------------------------------------------------
# MODULE 13: FIND MARKER REGULONS PER CELL TYPE
# ----------------------------------------------------------------------------

# Set identity to your cell types
DefaultAssay(seu_ds) <- "SCENIC"
Idents(seu_ds) <- "microglia_type"  # replace with your cell type column

# Find marker regulons for each cell type
regulon_markers <- FindAllMarkers(seu_ds, 
                                  only.pos = TRUE,
                                  min.pct = 0.25)

# Get top 5 regulons per cell type
top5 <- regulon_markers %>% 
  group_by(cluster) %>% 
  top_n(5, avg_log2FC)

print(top5)


# ----------------------------------------------------------------------------
# MODULE 14: PLOT CELL TYPE-SPECIFIC MARKER REGULONS
# ----------------------------------------------------------------------------

# Get unique top regulons per cell type
top_regulons_per_ct <- regulon_markers %>%
  group_by(cluster) %>%
  top_n(5, avg_log2FC) %>%
  pull(gene) %>%
  unique()

# Subset average AUC matrix
auc_markers <- auc_by_celltype[top_regulons_per_ct, ]

# Plot marker regulons
pheatmap(auc_markers,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Cell Type-Specific Marker Regulons",
         fontsize_row = 7,
         angle_col = 45)


# ----------------------------------------------------------------------------
# MODULE 15: REGULON-REGULON CORRELATION ANALYSIS
# ----------------------------------------------------------------------------

# Calculate correlation between top regulons
regulon_cor <- cor(t(auc_matrix[top_regulons_20, ]))

# Plot correlation heatmap
pheatmap(regulon_cor,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Regulon-Regulon Correlation",
         fontsize = 8)


# ----------------------------------------------------------------------------
# MODULE 16: BINARY REGULON ACTIVITY HEATMAP
# ----------------------------------------------------------------------------

# Binarize AUC scores (active vs inactive)
auc_binary <- auc_subset
thresholds <- apply(auc_subset, 1, function(x) quantile(x, 0.75))

for(i in 1:nrow(auc_binary)) {
  auc_binary[i, ] <- ifelse(auc_subset[i, ] > thresholds[i], 1, 0)
}

# Plot binary heatmap
pheatmap(auc_binary,
         annotation_col = annotation_col,
         show_colnames = FALSE,
         color = c("pink", "darkred"),
         legend_breaks = c(0, 1),
         legend_labels = c("Inactive", "Active"),
         main = "Binary Regulon Activity",
         cluster_rows = TRUE,
         cluster_cols = TRUE)