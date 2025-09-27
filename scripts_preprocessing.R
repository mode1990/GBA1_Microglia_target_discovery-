###############################################################################
# Title:     scRNA-seq Preprocessing, Doublet Removal, Integration & Clustering
# Author:    DZNE Bioinformatics - Mo Dehestani
# Date:      2025-06-19
# Project:   Single-cell transcriptomics (GEX) from 10x Genomics data
# Purpose:   Load, QC, remove doublets, normalize, integrate, and cluster 
#            single-cell RNA-seq data across 4 conditions (c1–c4)
# Input:     10X-formatted directories per sample (filtered_feature_bc_matrix/)
# Output:    Final integrated Seurat object: pbmc_final.rds
# Dependencies:
#    - Seurat >= 4.0
#    - scDblFinder
#    - harmony
#    - Matrix, dplyr, SingleCellExperiment
###############################################################################

# 1. Load Required Libraries
library(Seurat) #v5
library(dplyr)
library(Matrix)
library(scDblFinder)
library(SingleCellExperiment)
library(harmony)
library(ggplot2)

library(Seurat)

library(Seurat)

library(Seurat)

GBA1_mut_24h_LPS <- CreateSeuratObject(
  counts = Read10X("/data/cellranger_multi_outputs/GBA1-mut-24h/outs/per_sample_outs/GBA1-mut-24h-LPS/count/sample_filtered_feature_bc_matrix"),
  min.cells = 3, min.features = 200, project = "GBA1-mut-24h-LPS"
)

GBA1_mut_24h_aSyn <- CreateSeuratObject(
  counts = Read10X("/data/cellranger_multi_outputs/GBA1-mut-24h/outs/per_sample_outs/GBA1-mut-24h-aSyn/count/sample_filtered_feature_bc_matrix"),
  min.cells = 3, min.features = 200, project = "GBA1-mut-24h-aSyn"
)

GBA1_mut_24h_ctrl <- CreateSeuratObject(
  counts = Read10X("/data/cellranger_multi_outputs/GBA1-mut-24h/outs/per_sample_outs/GBA1-mut-24h-ctrl/count/sample_filtered_feature_bc_matrix"),
  min.cells = 3, min.features = 200, project = "GBA1-mut-24h-ctrl"
)

GBA1_mut_24h_debris <- CreateSeuratObject(
  counts = Read10X("/data/cellranger_multi_outputs/GBA1-mut-24h/outs/per_sample_outs/GBA1-mut-24h-debris/count/sample_filtered_feature_bc_matrix"),
  min.cells = 3, min.features = 200, project = "GBA1-mut-24h-debris"
)

GBA1_mut_7d_LPS <- CreateSeuratObject(
  counts = Read10X("/data/cellranger_multi_outputs/GBA1-mut-7d/outs/per_sample_outs/GBA1-mut-7d-LPS/count/sample_filtered_feature_bc_matrix"),
  min.cells = 3, min.features = 200, project = "GBA1-mut-7d-LPS"
)

GBA1_mut_7d_aSyn <- CreateSeuratObject(
  counts = Read10X("/data/cellranger_multi_outputs/GBA1-mut-7d/outs/per_sample_outs/GBA1-mut-7d-aSyn/count/sample_filtered_feature_bc_matrix"),
  min.cells = 3, min.features = 200, project = "GBA1-mut-7d-aSyn"
)

GBA1_mut_7d_ctrl <- CreateSeuratObject(
  counts = Read10X("/data/cellranger_multi_outputs/GBA1-mut-7d/outs/per_sample_outs/GBA1-mut-7d-ctrl/count/sample_filtered_feature_bc_matrix"),
  min.cells = 3, min.features = 200, project = "GBA1-mut-7d-ctrl"
)

GBA1_mut_7d_debris <- CreateSeuratObject(
  counts = Read10X("/data/cellranger_multi_outputs/GBA1-mut-7d/outs/per_sample_outs/GBA1-mut-7d-debris/count/sample_filtered_feature_bc_matrix"),
  min.cells = 3, min.features = 200, project = "GBA1-mut-7d-debris"
)

KOLF_24h_LPS <- CreateSeuratObject(
  counts = Read10X("/data/cellranger_multi_outputs/KOLF-24h/outs/per_sample_outs/KOLF-24h-LPS/count/sample_filtered_feature_bc_matrix"),
  min.cells = 3, min.features = 200, project = "KOLF-24h-LPS"
)

KOLF_24h_aSyn <- CreateSeuratObject(
  counts = Read10X("/data/cellranger_multi_outputs/KOLF-24h/outs/per_sample_outs/KOLF-24h-aSyn/count/sample_filtered_feature_bc_matrix"),
  min.cells = 3, min.features = 200, project = "KOLF-24h-aSyn"
)

KOLF_24h_ctrl <- CreateSeuratObject(
  counts = Read10X("/data/cellranger_multi_outputs/KOLF-24h/outs/per_sample_outs/KOLF-24h-ctrl/count/sample_filtered_feature_bc_matrix"),
  min.cells = 3, min.features = 200, project = "KOLF-24h-ctrl"
)

KOLF_24h_debris <- CreateSeuratObject(
  counts = Read10X("/data/cellranger_multi_outputs/KOLF-24h/outs/per_sample_outs/KOLF-24h-debris/count/sample_filtered_feature_bc_matrix"),
  min.cells = 3, min.features = 200, project = "KOLF-24h-debris"
)

KOLF_7d_LPS <- CreateSeuratObject(
  counts = Read10X("/data/cellranger_multi_outputs/KOLF-7d/outs/per_sample_outs/KOLF-7d-LPS/count/sample_filtered_feature_bc_matrix"),
  min.cells = 3, min.features = 200, project = "KOLF-7d-LPS"
)

KOLF_7d_aSyn <- CreateSeuratObject(
  counts = Read10X("/data/cellranger_multi_outputs/KOLF-7d/outs/per_sample_outs/KOLF-7d-aSyn/count/sample_filtered_feature_bc_matrix"),
  min.cells = 3, min.features = 200, project = "KOLF-7d-aSyn"
)

KOLF_7d_ctrl <- CreateSeuratObject(
  counts = Read10X("/data/cellranger_multi_outputs/KOLF-7d/outs/per_sample_outs/KOLF-7d-ctrl/count/sample_filtered_feature_bc_matrix"),
  min.cells = 3, min.features = 200, project = "KOLF-7d-ctrl"
)

KOLF_7d_debris <- CreateSeuratObject(
  counts = Read10X("/data/cellranger_multi_outputs/KOLF-7d/outs/per_sample_outs/KOLF-7d-debris/count/sample_filtered_feature_bc_matrix"),
  min.cells = 3, min.features = 200, project = "KOLF-7d-debris"
)

# Vector of your Seurat object variable names as strings
sample_ids <- c(
  "GBA1_mut_24h_LPS", "GBA1_mut_24h_aSyn", "GBA1_mut_24h_ctrl", "GBA1_mut_24h_debris",
  "GBA1_mut_7d_LPS", "GBA1_mut_7d_aSyn", "GBA1_mut_7d_ctrl", "GBA1_mut_7d_debris",
  "KOLF_24h_LPS", "KOLF_24h_aSyn", "KOLF_24h_ctrl", "KOLF_24h_debris",
  "KOLF_7d_LPS", "KOLF_7d_aSyn", "KOLF_7d_ctrl", "KOLF_7d_debris"
)

for (id in sample_ids) {
  seurat_obj <- get(id)  # Fetch object by name
  
  # Add QC metrics
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL")
  
  # Pre-QC plot
  png(paste0(id, "_preQC.png"), width = 1600, height = 600, res = 150)
  print(VlnPlot(seurat_obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), 
                ncol = 4, pt.size = 0.1, log = TRUE))
  dev.off()
  
  # QC filtering
  seurat_obj <- subset(
    seurat_obj,
    subset = nCount_RNA < 25000 & nCount_RNA > 1000 & percent.mt < 10 & percent.ribo > 5
  )
  
  # Post-QC plot
  png(paste0(id, "_postQC.png"), width = 1600, height = 600, res = 150)
  print(VlnPlot(seurat_obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), 
                ncol = 4, pt.size = 0.1, log = TRUE))
  dev.off()
  
  # Add condition metadata
  seurat_obj$Condition <- id
  
  # Save to disk
  saveRDS(seurat_obj, paste0(id, ".rds"))
  
  # Reassign object back to the environment with the same name
  assign(id, seurat_obj)
}


#3 Load and annotate Seurat objects
seurat_objects <- lapply(sample_ids, function(id) {
  obj <- readRDS(paste0(id, ".rds"))
  obj$Condition <- id
  
  parts <- strsplit(id, "_")[[1]]
  
  # Example: parts = c("GBA1", "mut", "24h", "LPS")
  obj$Genotype <- parts[1]           # "GBA1" or "KOLF"
  obj$Timepoint <- parts[3]          # "24h" or "7d"
  obj$Treatment <- parts[4]          # "LPS", "aSyn", "ctrl", "debris"
  
  return(obj)
})

# Name list entries by sample ID
names(seurat_objects) <- sample_ids

#4 doublet removal and merge 
library(SingleCellExperiment)
library(scDblFinder)

# Initialize filtered list
filtered_objects <- list()

# Loop over each Seurat object to run scDblFinder
for (id in sample_ids) {
  message("Processing ", id)
  
  # Convert to SCE
  sce <- as.SingleCellExperiment(seurat_objects[[id]])
  
  # Run scDblFinder
  sce <- scDblFinder(sce)
  
  # Save results back to Seurat object
  seurat_objects[[id]]$scDblFinder_class <- sce$scDblFinder.class
  
  # Keep singlets only
  filtered_objects[[id]] <- subset(seurat_objects[[id]], subset = scDblFinder_class == "singlet")
}

# ----- Doublet % stats and barplot -----
doublet_stats <- sapply(sample_ids, function(id) {
  sum(seurat_objects[[id]]$scDblFinder_class == "doublet") / ncol(seurat_objects[[id]]) * 100
})

# Plot
par(mar = c(5, 10, 4, 2))  # Bottom, Left, Top, Right
barplot(doublet_stats,
        main = "% Doublets per Sample",
        xlab = "Percentage",
        las = 1,                 # horizontal y-axis labels
        col = "tomato",
        horiz = TRUE)

# ----- Merge filtered Seurat objects -----
combined_filtered <- merge(
  filtered_objects[[1]],
  y = filtered_objects[-1],
  add.cell.ids = names(filtered_objects),
  project = "Merged_GEX"
)

# Save to disk
saveRDS(combined_filtered, "/data/outputs/cellranger_outputs/Analysis_by_Mo_Multi/all_merged_QCied_DblRemoved.rds")

# 5. Normalize, Integrate, and PCA
PD_list <- SplitObject(readRDS("all_merged_QCied_DblRemoved.rds"), split.by = "Condition")

library(future)
options(future.globals.maxSize = 2 * 1024^3)

for (i in 1:length(PD_list)) {
  PD_list[[i]] <- SCTransform(
    PD_list[[i]],
    vars.to.regress = c("percent.mt", "percent.ribo"),
    return.only.var.genes = FALSE,   
    verbose = TRUE
  )
}


# Merge 
# 1. Check the actual features in each object
lapply(PD_list, function(x) {
  cat("Total features:", nrow(x), "\n")
  cat("SCT features:", nrow(x[["SCT"]]), "\n")
  cat("First few SCT features:", head(rownames(x[["SCT"]])), "\n\n")
})

# 2. Find common features across ALL assays (not just SCT)
common_features_all <- Reduce(intersect, lapply(PD_list, rownames))
cat("Common features across all objects:", length(common_features_all), "\n")

# 3. Subset ALL objects to common features BEFORE merging
for (i in seq_along(PD_list)) {
  PD_list[[i]] <- subset(PD_list[[i]], features = common_features_all)
  cat(names(PD_list)[i], "- Features after subset:", nrow(PD_list[[i]]), "\n")
}
# 1. Merge objects
# After ensuring consistent features, try merging
PD <- merge(PD_list[[1]], y = PD_list[2:length(PD_list)], 
            add.cell.ids = names(PD_list),
            project = "GEX")
# 2. Run PCA
sct_features <- unique(unlist(lapply(PD_list, VariableFeatures)))
VariableFeatures(PD) <- sct_features
PD <- RunPCA(PD, npcs = 30, assay = "SCT", features = sct_features)

# 3. Use IntegrateLayers with correct syntax
PD <- IntegrateLayers(
  object = PD,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  group.by = "Condition"  
)

# 7. UMAP, Clustering, and Visualization
PD <- RunUMAP(PD, reduction = "harmony", dims = 1:30)
PD <- FindNeighbors(PD, reduction = "harmony", graph.name = "test", dims = 1:30)
PD <- FindClusters(PD, graph.name = "test", resolution = c(0.1, 0.2, 0.5, 0.8))

library(ggplot2)
DimPlot(PD, group.by = "test_res.0.1", label = TRUE) 
DimPlot(PD, group.by = "test_res.0.2", label = TRUE)  
DimPlot(PD, group.by = "test_res.0.5", label = TRUE)  
DimPlot(PD, group.by = "test_res.0.8", label = TRUE)  
DimPlot(PD, group.by = "Condition")
DimPlot(PD, group.by = "Timepoint")
DimPlot(PD, group.by = "Genotype")

# Save final object
saveRDS(PD, "/data/outputs/cellranger_outputs/Analysis_by_Mo_Multi/Processed_harmonized_unannotated.rds")

#clustering stability analysis
# Randomly sample cells to make silhouette calculation possible
set.seed(123)
library(cluster)

n_sample <- min(10000, ncol(PD))  # Sample max 10k cells
sample_cells <- sample(colnames(PD), n_sample)
dist_matrix_sub <- dist(Embeddings(PD, "pca")[sample_cells, 1:20])
sil_scores <- silhouette(as.numeric(PD$test_res.0.5[sample_cells]), dist_matrix_sub)
sil_summary <- aggregate(sil_scores[, "sil_width"], 
                         by = list(cluster = sil_scores[, "cluster"]), 
                         FUN = mean)
print(sil_summary)
library(factoextra)
fviz_silhouette(sil_scores) +
  theme_minimal() +
  labs(title = "Silhouette Analysis - Cluster Quality")

#Annotation by curated lead markers 
# iPSC-microglia specific markers at different differentiation stages from Dolan paper

Homeostatic_Microglia = c("CX3CR1", "P2RY12", "TMEM119", "SLC2A5", "HEXB", "CSF1R", "GPR34", "MEF2A")

Disease_Associated_Microglia = c("APOE", "TREM2", "GPNMB", "LPL", "ABCA1", "CD9", "CTSD", "CTSB", "SPP1", "CLEC7A")

Antigen_Presenting_Microglia = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "CD74", "CIITA", "B2M")

Interferon_Responsive_Microglia = c("IFIT1", "IFIT2", "IFIT3", "MX1", "MX2", "ISG15", "STAT1", "IRF7", "OAS1", "OAS2")

Proliferating_Microglia = c("MKI67", "TOP2A", "PCNA", "STMN1", "HIST1H4C", "HIST1H1B")

Transitional_Activated_Microglia = c("FOS", "JUN", "EGR1", "NR4A1", "DUSP1", "ZFP36")

DimPlot(PD_filtered, group.by = "test_res.0.5", label=T)
table(PD$test_res.0.5)
# Create comprehensive feature plot
FeaturePlot(PD, features = Homeostatic_Microglia)
FeaturePlot(PD, features = Disease_Associated_Microglia)
FeaturePlot(PD, features = Antigen_Presenting_Microglia)
FeaturePlot(PD, features = Interferon_Responsive_Microglia)
FeaturePlot(PD, features = Proliferating_Microglia)
FeaturePlot(PD, features = Transitional_Activated_Microglia)



# retrospective QC metrics comparison
VlnPlot(PD, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), 
        group.by = "test_res.0.5")

#remove artifacts 
# Define clusters to remove
Idents(PD) <- "test_res.0.5"
clusters_to_remove <- c(8, 13, 14, 15)  # Poor silhouette + QC issues
PD_filtered <- subset(PD, idents = clusters_to_remove, invert = TRUE)
DimPlot(PD_filtered, group.by = "test_res.0.5", label=T)

#now lets annotate  
cluster_labels <- c(
  "0" = "Homeostatic Microglia", 
  "1" = "Homeostatic Microglia", 
  "6" = "Homeostatic Microglia",
  
  "2" = "Disease-Associated Microglia", 
  "12" = "Disease-Associated Microglia", 
  
  "4" = "Proliferating Microglia", 
  "5" = "Proliferating Microglia", 
  "7" = "Proliferating Microglia", 
  "9" = "Proliferating Microglia", 
  
  "3" = "Antigen-Presenting Microglia",
  
  "10" = "Interferon-Responsive Microglia",
  "11" = "Interferon-Responsive Microglia"
  
)

Idents(PD_filtered)<-"test_res.0.5"
levels(Idents(PD_filtered))
# Get cluster identities
cluster_ids <- as.character(Idents(PD_filtered))
# Map subtype labels, and remove names
microglia_labels <- unname(cluster_labels[cluster_ids])
# Assign directly (by position, not by name)
PD_filtered$microglia_type <- microglia_labels
table(PD_filtered$microglia_type, useNA = "always")
DimPlot(PD_filtered, group.by = "microglia_type") +
  ggtitle("Annotated Microglia Subtypes")

#confirm labels 
DotPlot(PD_filtered, group.by = "microglia_type", features = Homeostatic_Microglia) + RotatedAxis()
DotPlot(PD_filtered, group.by = "microglia_type", features = Disease_Associated_Microglia) + RotatedAxis()
DotPlot(PD_filtered, group.by = "microglia_type", features = Antigen_Presenting_Microglia) + RotatedAxis()
DotPlot(PD_filtered, group.by = "microglia_type", features = Interferon_Responsive_Microglia) + RotatedAxis()
DotPlot(PD_filtered, group.by = "microglia_type", features = Proliferating_Microglia) + RotatedAxis()


#save
saveRDS(PD_filtered, "/data/outputs/cellranger_outputs/Analysis_by_Mo_Multi/Processed_harmonized_artifacts_removed_annotated_by_Dolan_paper.rds")



# Find markers for each cluster compared to all other cells
Idents(PD_filtered) <- "microglia_type"
DefaultAssay(PD_filtered) <- "RNA"
PD_filtered <- NormalizeData(PD_filtered)
PD_filtered <- ScaleData(PD_filtered)
PD_filtered <- JoinLayers(PD_filtered)
markers <- FindAllMarkers(PD_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(markers)
significant_markers <- markers %>% 
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::arrange(desc(avg_log2FC))

write.csv(significant_markers, file = "/data/outputs/cellranger_outputs/Analysis_by_Mo_Multi/significant_cluster_markers_minpct25_logfc25.csv", row.names = FALSE)

significant_markers <- fread("significant_cluster_markers_minpct25_logfc25.csv")
# View markers per cluster
table(significant_markers$cluster)

library(dplyr)
top_markers <- significant_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

Idents(PD)="microglia_type"
PD_downsampled <- subset(PD, downsample = 5000)

pdf("microglia_heatmap.pdf", width = 10, height = 8)
DefaultAssay(PD_downsampled) <- "SCT"
print(DoHeatmap(PD_downsampled, features = top_markers$gene))
dev.off()



home_markers <- c(
  "P2RY12",
  "TMEM119",
  "CX3CR1",
  "TREM2",
  "FCRLS",
  "SALL1",
  "SIGLEC-H",
  "OLFM3",
  "P2RY13",
  "TGFBR1"
)
DotPlot(PD_filtered, features = home_markers)+RotatedAxis()
table(PD@meta.data$Timepoint)
table(PD@meta.data$Genotype)
table(PD@meta.data$Condition)
table(PD@meta.data$Treatment)


#metadata visualization 

# Extract metadata
meta <- as.data.frame(PD@meta.data)

# Plot 1: Stacked barplot of celltype proportions by Timepoint
tab1 <- table(meta$Timepoint, meta$microglia_type)
prop1 <- prop.table(tab1, margin = 1) * 100
plot_data1 <- as.data.frame(prop1)
colnames(plot_data1) <- c("Timepoint", "microglia_type", "Percent")

p1 <- ggplot(plot_data1, aes(x = Timepoint, y = Percent, fill = microglia_type)) +
  geom_bar(stat = "identity") +
  labs(title = "Microglia subtype % by Timepoint", y = "Percentage", x = "Timepoint") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 2: Stacked barplot of celltype proportions by Genotype
tab2 <- table(meta$Genotype, meta$microglia_type)
prop2 <- prop.table(tab2, margin = 1) * 100
plot_data2 <- as.data.frame(prop2)
colnames(plot_data2) <- c("Genotype", "microglia_type", "Percent")

p2 <- ggplot(plot_data2, aes(x = Genotype, y = Percent, fill = microglia_type)) +
  geom_bar(stat = "identity") +
  labs(title = "Microglia subtype % by Genotype", y = "Percentage", x = "Genotype") 

# Plot 3: Stacked barplot of celltype proportions by Condition
tab3 <- table(meta$Condition, meta$microglia_type)
prop3 <- prop.table(tab3, margin = 1) * 100
plot_data3 <- as.data.frame(prop3)
colnames(plot_data3) <- c("Condition", "microglia_type", "Percent")

p3 <- ggplot(plot_data3, aes(x = Condition, y = Percent, fill = microglia_type)) +
  geom_bar(stat = "identity") +
  labs(title = "Microglia subtype % by Sample", y = "Percentage", x = "Condition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 4: Stacked barplot of celltype proportions by Treatment 
tab4 <- table(meta$Treatment, meta$microglia_type)
prop4 <- prop.table(tab4, margin = 1) * 100
plot_data4 <- as.data.frame(prop4)
colnames(plot_data4) <- c("Treatment", "microglia_type", "Percent")

p4 <- ggplot(plot_data4, aes(x = Treatment, y = Percent, fill = microglia_type)) +
  geom_bar(stat = "identity") +
  labs(title = "Microglia subtype % by Treatment", y = "Percentage", x = "Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display all plots
print(p1)
print(p2)
print(p3)
print(p4)



# Optional: Save plots
ggsave("/data/outputs/cellranger_outputs/Analysis_by_Mo_Multi/proportions_microglia_by_timepoint.png", p1, width = 8, height = 6)
ggsave("/data/outputs/cellranger_outputs/Analysis_by_Mo_Multi/proportions_microglia_by_genotype.png", p2, width = 8, height = 6)
ggsave("/data/outputs/cellranger_outputs/Analysis_by_Mo_Multi/proportions_microglia_by_sample.png", p3, width = 8, height = 6)
ggsave("/data/outputs/cellranger_outputs/Analysis_by_Mo_Multi/proportions_microglia_by_treatment.png", p4, width = 8, height = 6)


#proportions statistics
library(speckle)
df <- data.frame(
  cluster = PD_filtered$microglia_type,         # or seurat_clusters
  sample = PD_filtered$Condition,
  geno = PD_filtered$Genotype,
  time = PD_filtered$Timepoint,
  treatment = PD_filtered$Treatment
  
)

res_geno <- propeller(cluster = df$cluster,sample = df$sample,group = df$geno,transform = "asin")
res_trt <- propeller(cluster = df$cluster,sample = df$sample,group = df$treatment,transform = "asin")
res_tp <- propeller(cluster = df$cluster,sample = df$sample,group = df$time,transform = "asin")

head(res_geno)
head(res_trt)
head(res_tp)


#pseudotime

# slinghsot
# Convert to SCE
sce <- as.SingleCellExperiment(PD_filtered)
# Add clustering if needed (e.g., microglia_type or seurat_clusters)
sce$cluster <- PD_filtered$microglia_type  
library(slingshot)
# Use the UMAP embedding and cluster labels
sce <- slingshot(sce, clusterLabels = 'cluster', reducedDim = 'UMAP', start.clus = 'Homeostatic Microglia')
# View pseudotime
head(slingPseudotime(sce))
library(scater)
plotUMAP(sce, colour_by = 'cluster')  # original clusters
plotUMAP(sce, colour_by = 'slingPseudotime_1')  # pseudotime trajectory
plot(reducedDims(sce)$UMAP, col = viridis::viridis(100)[cut(slingPseudotime(sce)[,1], breaks = 100)], pch = 16)
lines(SlingshotDataSet(sce), lwd = 2)

#monocle
# Assuming PD_filtered is your Seurat object
data <- as.matrix(PD_filtered@assays$RNA@layers$counts)
rownames(data) <- rownames(PD_filtered@assays$RNA)  
pd <- PD_filtered@meta.data
fd <- data.frame(gene_short_name = rownames(data))
rownames(fd) <- rownames(data)
# Create Monocle CellDataSet
monocle_cds <- newCellDataSet(data,
                              phenoData = new("AnnotatedDataFrame", data = pd),
                              featureData = new("AnnotatedDataFrame", data = fd),
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
disp_table <- dispersionTable(monocle_cds)
ordering_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1.5 * dispersion_fit)$gene_id
#as monocle2 does not scale up lets downsample randomly 
set.seed(42)
sub_cells <- sample(colnames(monocle_cds), 10000)
monocle_sub <- monocle_cds[, sub_cells]
monocle_sub <- reduceDimension(monocle_sub, max_components = 2, method = "DDRTree")
monocle_sub <- orderCells(monocle_sub)
plot_cell_trajectory(monocle_sub, color_by = "Pseudotime")
plot_cell_trajectory(monocle_sub, color_by = "microglia_type")  # replace with your cluster column name
plot_cell_trajectory(monocle_sub, color_by = "State")

table(pData(monocle_sub)$State, pData(monocle_sub)$microglia_type)

monocle_sub <- orderCells(monocle_sub, root_state = "4")  

plot_cell_trajectory(monocle_sub, color_by = "Pseudotime")
plot_cell_trajectory(monocle_sub, color_by = "State")
plot_cell_trajectory(monocle_sub, color_by = "Genotype")
plot_cell_trajectory(monocle_sub, color_by = "Timepoint")
plot_cell_trajectory(monocle_sub, color_by = "Treatment")


# Test for genes that vary over pseudotime
deg_pseudotime <- differentialGeneTest(monocle_sub,
                                       fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                       cores = 4)

# Filter top genes
sig_genes <- row.names(subset(deg_pseudotime, qval < 0.01))

# Plot pseudotime heatmap
plot_pseudotime_heatmap(monocle_sub[sig_genes, ],
                        num_clusters = 4,
                        show_rownames = FALSE)


deg_genotype <- differentialGeneTest(monocle_sub,
                                     fullModelFormulaStr = "~Genotype + Pseudotime",
                                     cores = 4)

# View top DE genes
head(deg_genotype[order(deg_genotype$qval), ])
sig_genes <- row.names(subset(deg_genotype, qval < 0.01))
plot_pseudotime_heatmap(monocle_sub[sig_genes, ],
                        num_clusters = 4,
                        show_rownames = FALSE)








#Module scores based on gene sets from Dolan et al paper 
# Module score for microglia identity, maturity and artifactual activation

# 1. MICROGLIA IDENTITY GENES (from reference 21)
microglia_identity_genes <- c(
  "P2RY12", "GPR34", "C1QA", "CX3CR1", "CABLES1", 
  "BHLHE41", "TREM2", "OLFML3", "PROS1", "APOE", 
  "SLCO2B1", "SLC7A8", "PPARD", "CRYBB1"
)

# 2. MICROGLIA MATURITY GENES (top 30 from reference 15)
microglia_maturity_genes <- c(
  "SPP1", "CD74", "ACTB", "C3", "FTL", "FOS", "CSF1R", 
  "B2M", "C1QC", "C1QB", "PSAP", "A2M", "ITM2B", 
  "LAPTM5", "CTSB", "P2RY12", "C1QA", "SLCO2B1", 
  "RGS1", "APOE", "CCL4L2", "RNASET2", "NEAT1", 
  "CX3CR1", "DUSP1", "SAT1", "ZFP36", "CD81", 
  "HLA-B", "HLA-DRA"
)

# 3. ARTIFACTUAL ACTIVATION GENES (ExAMs - from reference 22)
# These are genes upregulated during single-cell isolation and cell handling
artifactual_activation_genes <- c(
  "RGS1", "HIST2H2AA1", "HIST1H4I", "NFKBIZ", "KLF2", 
  "JUNB", "DUSP1", "CCL3", "HSPA1A", "HSP90AA1", 
  "FOS", "HSPA1B", "JUN", "JUND", "NFKBID", "GEM", 
  "CCL4", "IER5", "TXNIP", "HIST1H2BC", "ZFP36", 
  "HIST1H1C", "EGR1", "ATF3", "RHOB"
)

# Print gene sets for verification
cat("Microglia Identity Genes (n =", length(microglia_identity_genes), "):\n")
cat(paste(microglia_identity_genes, collapse = ", "), "\n\n")

cat("Microglia Maturity Genes (n =", length(microglia_maturity_genes), "):\n")
cat(paste(microglia_maturity_genes, collapse = ", "), "\n\n")

cat("Artifactual Activation Genes (n =", length(artifactual_activation_genes), "):\n")
cat(paste(artifactual_activation_genes, collapse = ", "), "\n\n")

# Load the gene sets above

# Add module scores using Seurat's AddModuleScore function
# Control size of 25 as specified in the paper
PD <- AddModuleScore(
  PD_filtered, 
  features = list(microglia_identity_genes), 
  name = "Microglia_Identity_Score",
  ctrl = 25
)

PD <- AddModuleScore(
  PD, 
  features = list(microglia_maturity_genes), 
  name = "Microglia_Maturity_Score",
  ctrl = 25
)

PD <- AddModuleScore(
  PD, 
  features = list(artifactual_activation_genes), 
  name = "Artifactual_Activation_Score",
  ctrl = 25
)

# Note: AddModuleScore automatically adds "1" to the end of the name
# So your scores will be: Microglia_Identity_Score1, Microglia_Maturity_Score1, etc.

library(patchwork)

# Feature plots of the module scores
p1 <- FeaturePlot(PD, features = "Microglia_Identity_Score1") + 
  ggtitle("Microglia Identity Score")
p2 <- FeaturePlot(PD, features = "Microglia_Maturity_Score1") + 
  ggtitle("Microglia Maturity Score")  
p3 <- FeaturePlot(PD, features = "Artifactual_Activation_Score1") + 
  ggtitle("Artifactual Activation Score")

p1 | p2 | p3

# Violin plots by cluster
VlnPlot(PD, features = c("Microglia_Identity_Score1", "Microglia_Maturity_Score1", "Artifactual_Activation_Score1"), 
        group.by = "microglia_type", ncol = 3)

# Ridge plots for distribution
library(ggridges)
RidgePlot(PD_filtered, features = "Microglia_Identity_Score1", group.by = "microglia_type")
RidgePlot(PD_filtered, features = "Microglia_Maturity_Score1", group.by = "microglia_type")

# Get average scores per cluster
cluster_scores <- PD@meta.data %>%
  group_by(test_res.0.2) %>%
  summarise(
    Identity_Score = mean(Microglia_Identity_Score1),
    Maturity_Score = mean(Microglia_Maturity_Score1),
    Artifactual_Score = mean(Artifactual_Activation_Score1),
    .groups = 'drop'
  )

print(cluster_scores)

# Create annotation based on scores (adjust thresholds based on your data)
PD$Paper_Based_Annotation <- case_when(
  PD$Microglia_Identity_Score1 > 0.2 & PD$Microglia_Maturity_Score1 > 0.2 & PD$Artifactual_Activation_Score1 < 0.1 ~ "Mature_Microglia",
  PD$Microglia_Identity_Score1 > 0.1 & PD$Microglia_Maturity_Score1 < 0.1 ~ "Immature_Microglia",
  PD$Artifactual_Activation_Score1 > 0.2 ~ "Artifactual_Activation",
  PD$Microglia_Identity_Score1 < 0 ~ "Non_Microglia",
  TRUE ~ "Intermediate"
)

# Visualize the annotation
DimPlot(PD, group.by = "Paper_Based_Annotation", label = TRUE)

# Check which genes from the paper are present in your data
genes_present_identity <- microglia_identity_genes[microglia_identity_genes %in% rownames(PD)]
genes_present_maturity <- microglia_maturity_genes[microglia_maturity_genes %in% rownames(PD)]
genes_present_artifactual <- artifactual_activation_genes[artifactual_activation_genes %in% rownames(PD)]

cat("Identity genes present:", length(genes_present_identity), "/", length(microglia_identity_genes), "\n")
cat("Maturity genes present:", length(genes_present_maturity), "/", length(microglia_maturity_genes), "\n")
cat("Artifactual genes present:", length(genes_present_artifactual), "/", length(artifactual_activation_genes), "\n")

# Missing genes
cat("\nMissing Identity genes:", setdiff(microglia_identity_genes, genes_present_identity), "\n")
cat("Missing Maturity genes:", setdiff(microglia_maturity_genes, genes_present_maturity), "\n")

# Summary dot plot of key genes from each module
key_genes <- c("P2RY12", "CX3CR1", "C1QA", "TREM2",  # Identity
               "SPP1", "CD74", "CSF1R", "C1QB",      # Maturity  
               "FOS", "JUN", "EGR1", "CCL3")         # Artifactual

DotPlot(PD, features = key_genes, group.by = "test_res.0.2") + 
  RotatedAxis() + 
  ggtitle("Key Genes from Paper Module Scores")

# Correlation plot between scores
library(corrplot)
score_matrix <- PD@meta.data[, c("Microglia_Identity_Score1", "Microglia_Maturity_Score1", "Artifactual_Activation_Score1")]
cor_matrix <- cor(score_matrix, use = "complete.obs")
corrplot(cor_matrix, method = "color", addCoef.col = "black")
