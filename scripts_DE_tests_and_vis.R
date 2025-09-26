# ==================== CHECK NORMALIZATION STATUS ====================

# 1. Check what assays are available
cat("Available assays:\n")
print(names(PD@assays))

# 2. Check current default assay
cat("\nCurrent default assay:", DefaultAssay(PD), "\n")

# ==================== CHECK NORMALIZATION STATUS (Seurat v5) ====================

# 1. Check what assays are available
cat("Available assays:\n")
print(names(PD@assays))

# 2. Check current default assay
cat("\nCurrent default assay:", DefaultAssay(PD), "\n")

# 3. Check Seurat version
cat("Seurat version:", packageVersion("Seurat"), "\n")

# 4. Check if RNA assay has normalized data (Seurat v5 compatible)
if("RNA" %in% names(PD@assays)) {
  cat("\n=== RNA ASSAY ===\n")
  
  # Use GetAssayData for v5 compatibility
  tryCatch({
    # Check normalized data
    data_matrix <- GetAssayData(PD, assay = "RNA", slot = "data")
    
    if(nrow(data_matrix) > 0 && ncol(data_matrix) > 0) {
      # Sample some values to check
      sample_indices <- min(100, nrow(data_matrix))
      sample_cols <- min(10, ncol(data_matrix))
      sample_values <- as.numeric(data_matrix[1:sample_indices, 1:sample_cols])
      sample_values <- sample_values[sample_values > 0 & !is.na(sample_values)]  # Remove zeros and NAs
      
      if(length(sample_values) > 0) {
        cat("RNA data slot - sample values (first non-zero):\n")
        print(head(sample_values, 10))
        cat("Range of non-zero values:", range(sample_values), "\n")
        
        # Log-normalized data typically ranges from ~0.5 to ~8
        if(max(sample_values, na.rm = TRUE) < 15) {
          cat("✓ RNA data appears to be log-normalized (values < 15)\n")
        } else {
          cat("✗ RNA data may NOT be log-normalized (values > 15)\n")
        }
      } else {
        cat("All sampled values are zero\n")
      }
    } else {
      cat("✗ RNA data slot is empty\n")
    }
    
    # Check raw counts
    counts_matrix <- GetAssayData(PD, assay = "RNA", slot = "counts")
    if(nrow(counts_matrix) > 0) {
      sample_counts <- as.numeric(counts_matrix[1:min(100, nrow(counts_matrix)), 1:min(10, ncol(counts_matrix))])
      sample_counts <- sample_counts[sample_counts > 0 & !is.na(sample_counts)]
      if(length(sample_counts) > 0) {
        cat("Raw counts range:", range(sample_counts), "\n")
      }
    }
    
    # Check if scale.data exists
    tryCatch({
      scale_matrix <- GetAssayData(PD, assay = "RNA", slot = "scale.data")
      if(nrow(scale_matrix) > 0) {
        cat("✓ RNA scale.data exists (scaled/centered)\n")
      } else {
        cat("✗ RNA scale.data is empty\n")
      }
    }, error = function(e) {
      cat("✗ RNA scale.data is empty or not accessible\n")
    })
    
  }, error = function(e) {
    cat("Error accessing RNA assay:", e$message, "\n")
  })
}

# ==================== Now time for DE tests ====================
library(Seurat)
library(dplyr)

# Define exotic microglia color palette for consistent plotting
microglia_colors_exotic <- c(
  "Antigen-Presenting Microglia"      = "#E63946",  # Flamingo red
  "Disease-Associated Microglia"      = "#6A4C93",  # Mystic violet
  "Homeostatic Microglia"             = "#2A9D8F",  # Deep jade green
  "Interferon-Responsive Microglia"   = "#4CC9F0",  # Electric cyan
  "Proliferating Microglia"           = "#F4A261"   # Spicy mango
)

# ==================== PRIORITY 1: DISEASE BIOMARKERS ====================

# 1A. GBA1 vs KOLF across all microglia (broad disease signature)
Idents(PD) <- "Genotype"
DefaultAssay(PD) <- "RNA"
de_genotype_all <- FindMarkers(PD, 
                               ident.1 = "GBA1", 
                               ident.2 = "KOLF",
                               test.use = "wilcox",
                               logfc.threshold = 0.25,
                               min.pct = 0.1)

# Add gene names and filter significant
de_genotype_all$gene <- rownames(de_genotype_all)
de_genotype_sig <- de_genotype_all[de_genotype_all$p_val_adj < 0.05, ]

# =================DE in Homeostatic Microglia=================

celltype <- "Homeostatic Microglia"
# Set identity to Genotype
Idents(PD) <- "Genotype"
# Subset cells of interest
subset_cells <- WhichCells(PD, expression = microglia_type == celltype)
cells_gba1 <- WhichCells(PD, expression = Genotype == "GBA1" & microglia_type == celltype)
cells_kolf <- WhichCells(PD, expression = Genotype == "KOLF" & microglia_type == celltype)
if (length(cells_gba1) > 0 && length(cells_kolf) > 0) {
  
  PD_sub <- subset(PD, cells = subset_cells)
  Idents(PD_sub) <- "Genotype"
  
  de_homeo <- FindMarkers(PD_sub,
                          ident.1 = "GBA1",
                          ident.2 = "KOLF",
                          test.use = "wilcox",
                          logfc.threshold = 0.25,
                          min.pct = 0.1)
    # Save DE results
  write.csv(de_homeo, file = paste0("de_", gsub(" ", "_", celltype), ".csv"), row.names = FALSE)
  de_genotype_sig <- de_homeo[de_homeo$p_val_adj < 0.05, ]
  if (nrow(de_genotype_sig) > 0) {
    de_genotype_sig$gene <- rownames(de_genotype_sig)
    
    dbs <- "GO_Biological_Process_2023"
    top_genes <- head(de_genotype_sig[order(de_genotype_sig$p_val_adj), "gene"], 100)
    
    enrichr_results <- enrichr(genes = top_genes, databases = dbs)
    go_bp <- enrichr_results[[dbs]]
    
    write.csv(go_bp, file = paste0("go_bp_", gsub(" ", "_", celltype), ".csv"), row.names = FALSE)
    
    cat("\n=== TOP GO BP TERMS (GBA1 vs KOLF - Homeostatic Microglia) ===\n")
    print(head(go_bp[, c("Term", "Adjusted.P.value", "Combined.Score")], 10))
    
    top_terms <- head(go_bp[order(go_bp$Adjusted.P.value), ], 10)
    p <- ggplot(top_terms, aes(x = reorder(Term, Combined.Score), y = Combined.Score)) +
      geom_col(fill = "#6A4C93") +
      coord_flip() +
      labs(title = paste("GO BP Enrichment (Top 10): GBA1 vs KOLF in", celltype),
           x = "", y = "Combined Score") +
      theme_minimal(base_size = 14)
    
    ggsave(filename = paste0("go_bp_", gsub(" ", "_", celltype), ".png"), plot = p, width = 8, height = 5)
    
  } else {
    message("Not enough significant DE genes for enrichment in ", celltype)
  }
  
} else {
  message("Not enough GBA1 or KOLF cells in ", celltype)
}


# =================DE in Disease Associated Microglia=================

celltype <- "Disease-Associated Microglia"
# Set identity to Genotype
Idents(PD) <- "Genotype"
# Subset cells of interest
subset_cells <- WhichCells(PD, expression = microglia_type == celltype)
cells_gba1 <- WhichCells(PD, expression = Genotype == "GBA1" & microglia_type == celltype)
cells_kolf <- WhichCells(PD, expression = Genotype == "KOLF" & microglia_type == celltype)
if (length(cells_gba1) > 0 && length(cells_kolf) > 0) {
  
  PD_sub <- subset(PD, cells = subset_cells)
  Idents(PD_sub) <- "Genotype"
  
  de_homeo <- FindMarkers(PD_sub,
                          ident.1 = "GBA1",
                          ident.2 = "KOLF",
                          test.use = "wilcox",
                          logfc.threshold = 0.25,
                          min.pct = 0.1)
  # Save DE results
  write.csv(de_homeo, file = paste0("de_", gsub(" ", "_", celltype), ".csv"), row.names = FALSE)
  de_genotype_sig <- de_homeo[de_homeo$p_val_adj < 0.05, ]
  if (nrow(de_genotype_sig) > 0) {
    de_genotype_sig$gene <- rownames(de_genotype_sig)
    
    dbs <- "GO_Biological_Process_2023"
    top_genes <- head(de_genotype_sig[order(de_genotype_sig$p_val_adj), "gene"], 100)
    
    enrichr_results <- enrichr(genes = top_genes, databases = dbs)
    go_bp <- enrichr_results[[dbs]]
    
    write.csv(go_bp, file = paste0("go_bp_", gsub(" ", "_", celltype), ".csv"), row.names = FALSE)
    
    cat("\n=== TOP GO BP TERMS (GBA1 vs KOLF - Disease-Associated Microglia) ===\n")
    print(head(go_bp[, c("Term", "Adjusted.P.value", "Combined.Score")], 10))
    
    top_terms <- head(go_bp[order(go_bp$Adjusted.P.value), ], 10)
    p <- ggplot(top_terms, aes(x = reorder(Term, Combined.Score), y = Combined.Score)) +
      geom_col(fill = "#6A4C93") +
      coord_flip() +
      labs(title = paste("GO BP Enrichment (Top 10): GBA1 vs KOLF in", celltype),
           x = "", y = "Combined Score") +
      theme_minimal(base_size = 14)
    
    ggsave(filename = paste0("go_bp_", gsub(" ", "_", celltype), ".png"), plot = p, width = 8, height = 5)
    
  } else {
    message("Not enough significant DE genes for enrichment in ", celltype)
  }
  
} else {
  message("Not enough GBA1 or KOLF cells in ", celltype)
}



# =================DE in Antigen-Presenting Microglia=================

celltype <- "Antigen-Presenting Microglia"
# Set identity to Genotype
Idents(PD) <- "Genotype"
# Subset cells of interest
subset_cells <- WhichCells(PD, expression = microglia_type == celltype)
cells_gba1 <- WhichCells(PD, expression = Genotype == "GBA1" & microglia_type == celltype)
cells_kolf <- WhichCells(PD, expression = Genotype == "KOLF" & microglia_type == celltype)
if (length(cells_gba1) > 0 && length(cells_kolf) > 0) {
  
  PD_sub <- subset(PD, cells = subset_cells)
  Idents(PD_sub) <- "Genotype"
  
  de_homeo <- FindMarkers(PD_sub,
                          ident.1 = "GBA1",
                          ident.2 = "KOLF",
                          test.use = "wilcox",
                          logfc.threshold = 0.25,
                          min.pct = 0.1)
  # Save DE results
  write.csv(de_homeo, file = paste0("de_", gsub(" ", "_", celltype), ".csv"), row.names = FALSE)
  de_genotype_sig <- de_homeo[de_homeo$p_val_adj < 0.05, ]
  if (nrow(de_genotype_sig) > 0) {
    de_genotype_sig$gene <- rownames(de_genotype_sig)
    
    dbs <- "GO_Biological_Process_2023"
    top_genes <- head(de_genotype_sig[order(de_genotype_sig$p_val_adj), "gene"], 100)
    
    enrichr_results <- enrichr(genes = top_genes, databases = dbs)
    go_bp <- enrichr_results[[dbs]]
    
    write.csv(go_bp, file = paste0("go_bp_", gsub(" ", "_", celltype), ".csv"), row.names = FALSE)
    
    cat("\n=== TOP GO BP TERMS (GBA1 vs KOLF - Antigen-Presenting Microglia) ===\n")
    print(head(go_bp[, c("Term", "Adjusted.P.value", "Combined.Score")], 10))
    
    top_terms <- head(go_bp[order(go_bp$Adjusted.P.value), ], 10)
    p <- ggplot(top_terms, aes(x = reorder(Term, Combined.Score), y = Combined.Score)) +
      geom_col(fill = "#6A4C93") +
      coord_flip() +
      labs(title = paste("GO BP Enrichment (Top 10): GBA1 vs KOLF in", celltype),
           x = "", y = "Combined Score") +
      theme_minimal(base_size = 14)
    
    ggsave(filename = paste0("go_bp_", gsub(" ", "_", celltype), ".png"), plot = p, width = 8, height = 5)
    
  } else {
    message("Not enough significant DE genes for enrichment in ", celltype)
  }
  
} else {
  message("Not enough GBA1 or KOLF cells in ", celltype)
}



# =================DE in Interferon-Responsive Microglia=================

celltype <- "Interferon-Responsive Microglia"
# Set identity to Genotype
Idents(PD) <- "Genotype"
# Subset cells of interest
subset_cells <- WhichCells(PD, expression = microglia_type == celltype)
cells_gba1 <- WhichCells(PD, expression = Genotype == "GBA1" & microglia_type == celltype)
cells_kolf <- WhichCells(PD, expression = Genotype == "KOLF" & microglia_type == celltype)
if (length(cells_gba1) > 0 && length(cells_kolf) > 0) {
  
  PD_sub <- subset(PD, cells = subset_cells)
  Idents(PD_sub) <- "Genotype"
  
  de_homeo <- FindMarkers(PD_sub,
                          ident.1 = "GBA1",
                          ident.2 = "KOLF",
                          test.use = "wilcox",
                          logfc.threshold = 0.25,
                          min.pct = 0.1)
  # Save DE results
  write.csv(de_homeo, file = paste0("de_", gsub(" ", "_", celltype), ".csv"), row.names = FALSE)
  de_genotype_sig <- de_homeo[de_homeo$p_val_adj < 0.05, ]
  if (nrow(de_genotype_sig) > 0) {
    de_genotype_sig$gene <- rownames(de_genotype_sig)
    
    dbs <- "GO_Biological_Process_2023"
    top_genes <- head(de_genotype_sig[order(de_genotype_sig$p_val_adj), "gene"], 100)
    
    enrichr_results <- enrichr(genes = top_genes, databases = dbs)
    go_bp <- enrichr_results[[dbs]]
    
    write.csv(go_bp, file = paste0("go_bp_", gsub(" ", "_", celltype), ".csv"), row.names = FALSE)
    
    cat("\n=== TOP GO BP TERMS (GBA1 vs KOLF - Interferon-Responsive Microglia) ===\n")
    print(head(go_bp[, c("Term", "Adjusted.P.value", "Combined.Score")], 10))
    
    top_terms <- head(go_bp[order(go_bp$Adjusted.P.value), ], 10)
    p <- ggplot(top_terms, aes(x = reorder(Term, Combined.Score), y = Combined.Score)) +
      geom_col(fill = "#6A4C93") +
      coord_flip() +
      labs(title = paste("GO BP Enrichment (Top 10): GBA1 vs KOLF in", celltype),
           x = "", y = "Combined Score") +
      theme_minimal(base_size = 14)
    
    ggsave(filename = paste0("go_bp_", gsub(" ", "_", celltype), ".png"), plot = p, width = 8, height = 5)
    
  } else {
    message("Not enough significant DE genes for enrichment in ", celltype)
  }
  
} else {
  message("Not enough GBA1 or KOLF cells in ", celltype)
}



# =================DE in Proliferating Microglia=================

celltype <- "Proliferating Microglia"
# Set identity to Genotype
Idents(PD) <- "Genotype"
# Subset cells of interest
subset_cells <- WhichCells(PD, expression = microglia_type == celltype)
cells_gba1 <- WhichCells(PD, expression = Genotype == "GBA1" & microglia_type == celltype)
cells_kolf <- WhichCells(PD, expression = Genotype == "KOLF" & microglia_type == celltype)
if (length(cells_gba1) > 0 && length(cells_kolf) > 0) {
  
  PD_sub <- subset(PD, cells = subset_cells)
  Idents(PD_sub) <- "Genotype"
  
  de_homeo <- FindMarkers(PD_sub,
                          ident.1 = "GBA1",
                          ident.2 = "KOLF",
                          test.use = "wilcox",
                          logfc.threshold = 0.25,
                          min.pct = 0.1)
  # Save DE results
  write.csv(de_homeo, file = paste0("de_", gsub(" ", "_", celltype), ".csv"), row.names = FALSE)
  de_genotype_sig <- de_homeo[de_homeo$p_val_adj < 0.05, ]
  if (nrow(de_genotype_sig) > 0) {
    de_genotype_sig$gene <- rownames(de_genotype_sig)
    
    dbs <- "GO_Biological_Process_2023"
    top_genes <- head(de_genotype_sig[order(de_genotype_sig$p_val_adj), "gene"], 100)
    
    enrichr_results <- enrichr(genes = top_genes, databases = dbs)
    go_bp <- enrichr_results[[dbs]]
    
    write.csv(go_bp, file = paste0("go_bp_", gsub(" ", "_", celltype), ".csv"), row.names = FALSE)
    
    cat("\n=== TOP GO BP TERMS (GBA1 vs KOLF - Proliferating Microglia) ===\n")
    print(head(go_bp[, c("Term", "Adjusted.P.value", "Combined.Score")], 10))
    
    top_terms <- head(go_bp[order(go_bp$Adjusted.P.value), ], 10)
    p <- ggplot(top_terms, aes(x = reorder(Term, Combined.Score), y = Combined.Score)) +
      geom_col(fill = "#6A4C93") +
      coord_flip() +
      labs(title = paste("GO BP Enrichment (Top 10): GBA1 vs KOLF in", celltype),
           x = "", y = "Combined Score") +
      theme_minimal(base_size = 14)
    
    ggsave(filename = paste0("go_bp_", gsub(" ", "_", celltype), ".png"), plot = p, width = 8, height = 5)
    
  } else {
    message("Not enough significant DE genes for enrichment in ", celltype)
  }
  
} else {
  message("Not enough GBA1 or KOLF cells in ", celltype)
}


# ==================== TREATMENT RESPONSES IN GBA1 CELLS ====================
# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(enrichR)

cat("=== ANALYZING TREATMENT EFFECTS IN GBA1 DAMs ===\n")

# Check if microglia_type column exists
if(!"microglia_type" %in% colnames(PD@meta.data)) {
  stop("Column 'microglia_type' not found in metadata. Please check column name.")
}

# Subset to GBA1 DAMs only
gba1_dam_cells <- WhichCells(PD, expression = Genotype == "GBA1" & microglia_type == "Disease-Associated Microglia")
PD_gba1_dam <- subset(PD, cells = gba1_dam_cells)

cat("Total GBA1 DAM cells found:", length(gba1_dam_cells), "\n")

if(length(gba1_dam_cells) == 0) {
  stop("No GBA1 DAM cells found. Check your microglia_type annotations and Genotype labels.")
}

# Check treatment distribution in GBA1 DAMs
treatment_counts <- table(PD_gba1_dam@meta.data$Treatment)
cat("Treatment distribution in GBA1 DAMs:\n")
print(treatment_counts)

# Test each treatment vs control in GBA1
treatments <- c("LPS", "aSyn", "debris")
de_treatments_gba1 <- list()

# Set identity to Treatment for the GBA1 DAM subset
Idents(PD_gba1_dam) <- "Treatment"

for(treatment in treatments) {
  cat("\n=== Processing", treatment, "vs Control in GBA1 DAMs ===\n")
  
  tryCatch({
    # Check if we have enough cells for each condition
    treatment_cells <- WhichCells(PD_gba1_dam, expression = Treatment == treatment)
    ctrl_cells <- WhichCells(PD_gba1_dam, expression = Treatment == "ctrl")
    
    cat("  -", treatment, "DAM cells:", length(treatment_cells), "\n")
    cat("  - Control DAM cells:", length(ctrl_cells), "\n")
    
    if(length(treatment_cells) >= 10 && length(ctrl_cells) >= 10) {
      
      cat("  - Running DE analysis...\n")
      de_result <- FindMarkers(PD_gba1_dam,
                               ident.1 = treatment,
                               ident.2 = "ctrl",
                               test.use = "wilcox",
                               logfc.threshold = 0.25,
                               min.pct = 0.1,
                               verbose = FALSE)
      
      # Add metadata
      de_result$gene <- rownames(de_result)
      de_result$treatment <- treatment
      
      # Filter for significant results
      de_sig <- de_result[de_result$p_val_adj < 0.05, ]
      de_treatments_gba1[[treatment]] <- de_sig
      
      # Save all results
      write.csv(de_result, file = paste0("de_gba1_dam_", treatment, "_vs_ctrl_all.csv"), row.names = TRUE)
      
      # Save only significant results
      if(nrow(de_sig) > 0) {
        write.csv(de_sig, file = paste0("de_gba1_dam_", treatment, "_vs_ctrl_sig.csv"), row.names = TRUE)
      }
      
      cat("  - Total DE genes:", nrow(de_result), "\n")
      cat("  - Significant DE genes (padj < 0.05):", nrow(de_sig), "\n")
      cat("  - Upregulated in", treatment, ":", sum(de_sig$avg_log2FC > 0), "\n")
      cat("  - Downregulated in", treatment, ":", sum(de_sig$avg_log2FC < 0), "\n")
      
      # Show top upregulated genes
      if(nrow(de_sig) > 0) {
        top_up <- head(de_sig[de_sig$avg_log2FC > 0, ][order(de_sig[de_sig$avg_log2FC > 0, ]$p_val_adj), ], 5)
        top_down <- head(de_sig[de_sig$avg_log2FC < 0, ][order(de_sig[de_sig$avg_log2FC < 0, ]$p_val_adj), ], 5)
        
        if(nrow(top_up) > 0) {
          cat("  - Top upregulated genes:", paste(rownames(top_up), collapse = ", "), "\n")
        }
        if(nrow(top_down) > 0) {
          cat("  - Top downregulated genes:", paste(rownames(top_down), collapse = ", "), "\n")
        }
        
        # Pathway enrichment analysis
        cat("  - Running pathway enrichment...\n")
        
        # Set enrichR databases
        dbs <- "GO_Biological_Process_2023"
        
        # Get top 100 most significant genes
        top_genes <- head(de_sig[order(de_sig$p_val_adj), "gene"], 100)
        
        tryCatch({
          enrichr_results <- enrichr(genes = top_genes, databases = dbs)
          go_bp <- enrichr_results[[dbs]]
          
          # Save enrichment results
          write.csv(go_bp, file = paste0("go_bp_gba1_dam_", treatment, "_vs_ctrl.csv"), row.names = FALSE)
          
          if(nrow(go_bp) > 0 && any(go_bp$Adjusted.P.value < 0.05)) {
            # Filter for significant terms
            go_bp_sig <- go_bp[go_bp$Adjusted.P.value < 0.05, ]
            
            cat("  - Significant GO BP terms found:", nrow(go_bp_sig), "\n")
            cat("  - Top GO BP terms for", treatment, "in DAMs:\n")
            
            top_5_terms <- head(go_bp_sig[, c("Term", "Adjusted.P.value", "Combined.Score")], 5)
            for(i in 1:nrow(top_5_terms)) {
              cat("    ", i, ". ", top_5_terms$Term[i], " (padj=", 
                  format(top_5_terms$Adjusted.P.value[i], digits=3), ")\n", sep="")
            }
            
            # Create visualization
            top_terms <- head(go_bp_sig[order(go_bp_sig$Adjusted.P.value), ], 10)
            
            p <- ggplot(top_terms, aes(x = reorder(Term, Combined.Score), y = Combined.Score)) +
              geom_col(fill = "#6A4C93") +
              coord_flip() +
              labs(title = paste("GO BP Enrichment (Top 10):", treatment, "vs Control in GBA1 DAMs"),
                   x = "GO Terms", y = "Combined Score") +
              theme_minimal(base_size = 12) +
              theme(axis.text.y = element_text(size = 10))
            
            ggsave(filename = paste0("go_bp_gba1_dam_", treatment, "_vs_ctrl.png"), 
                   plot = p, width = 12, height = 8, dpi = 300)
            
            cat("  - Enrichment plot saved for", treatment, "in DAMs\n")
            
          } else {
            cat("  - No significant GO BP terms found for", treatment, "in DAMs\n")
          }
          
        }, error = function(e) {
          cat("  - Enrichment analysis failed for", treatment, ":", e$message, "\n")
        })
        
      }
      
    } else {
      cat("  - Insufficient DAM cells for", treatment, "analysis (need ≥10 per group)\n")
    }
    
  }, error = function(e) {
    cat("  - ERROR with", treatment, ":", e$message, "\n")
  })
}

# Summary of all results
cat("\n=== SUMMARY OF TREATMENT EFFECTS IN GBA1 DAMs ===\n")
for(treatment in names(de_treatments_gba1)) {
  n_sig <- nrow(de_treatments_gba1[[treatment]])
  cat(treatment, ": ", n_sig, " significant DE genes in DAMs\n")
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Output files generated:\n")
cat("- de_gba1_dam_[treatment]_vs_ctrl_all.csv : All DE results in DAMs\n")
cat("- de_gba1_dam_[treatment]_vs_ctrl_sig.csv : Only significant results (padj < 0.05)\n")
cat("- go_bp_gba1_dam_[treatment]_vs_ctrl.csv : GO enrichment results for DAMs\n")
cat("- go_bp_gba1_dam_[treatment]_vs_ctrl.png : GO enrichment plots for DAMs\n")




# ==================== TIMEPOINT RESPONSES IN GBA1 DAMs ====================

cat("=== ANALYZING TIMEPOINT EFFECTS IN GBA1 DAMs ===\n")

# Check necessary metadata columns
if(!"microglia_type" %in% colnames(PD@meta.data)) {
  stop("Column 'microglia_type' not found in metadata.")
}
if(!"Timepoint" %in% colnames(PD@meta.data)) {
  stop("Column 'Timepoint' not found in metadata.")
}

# Subset to GBA1 DAMs only
gba1_dam_cells <- WhichCells(PD, expression = Genotype == "GBA1" & microglia_type == "Disease-Associated Microglia")
PD_gba1_dam <- subset(PD, cells = gba1_dam_cells)

cat("Total GBA1 DAM cells found:", length(gba1_dam_cells), "\n")

if(length(gba1_dam_cells) == 0) {
  stop("No GBA1 DAM cells found.")
}

# Check timepoint distribution
timepoint_counts <- table(PD_gba1_dam@meta.data$Timepoint)
cat("Timepoint distribution in GBA1 DAMs:\n")
print(timepoint_counts)

# Compare 7d vs 24h in GBA1 DAMs
Idents(PD_gba1_dam) <- "Timepoint"
de_timepoint_gba1 <- list()

# Define timepoints
tp1 <- "7d"
tp2 <- "24h"

cat("\n=== Processing", tp1, "vs", tp2, "in GBA1 DAMs ===\n")

tryCatch({
  tp1_cells <- WhichCells(PD_gba1_dam, expression = Timepoint == tp1)
  tp2_cells <- WhichCells(PD_gba1_dam, expression = Timepoint == tp2)
  
  cat("  -", tp1, "DAM cells:", length(tp1_cells), "\n")
  cat("  -", tp2, "DAM cells:", length(tp2_cells), "\n")
  
  if(length(tp1_cells) >= 10 && length(tp2_cells) >= 10) {
    
    cat("  - Running DE analysis...\n")
    de_result <- FindMarkers(PD_gba1_dam,
                             ident.1 = tp1,
                             ident.2 = tp2,
                             test.use = "wilcox",
                             logfc.threshold = 0.25,
                             min.pct = 0.1,
                             verbose = FALSE)
    
    de_result$gene <- rownames(de_result)
    de_result$comparison <- paste0(tp1, "_vs_", tp2)
    
    # Filter significant
    de_sig <- de_result[de_result$p_val_adj < 0.05, ]
    de_timepoint_gba1[[paste0(tp1, "_vs_", tp2)]] <- de_sig
    
    # Save results
    write.csv(de_result, paste0("de_gba1_dam_", tp1, "_vs_", tp2, "_all.csv"), row.names = TRUE)
    if(nrow(de_sig) > 0) {
      write.csv(de_sig, paste0("de_gba1_dam_", tp1, "_vs_", tp2, "_sig.csv"), row.names = TRUE)
    }
    
    cat("  - Total DE genes:", nrow(de_result), "\n")
    cat("  - Significant DE genes (padj < 0.05):", nrow(de_sig), "\n")
    cat("  - Upregulated in", tp1, ":", sum(de_sig$avg_log2FC > 0), "\n")
    cat("  - Downregulated in", tp1, ":", sum(de_sig$avg_log2FC < 0), "\n")
    
    if(nrow(de_sig) > 0) {
      top_up <- head(de_sig[de_sig$avg_log2FC > 0, ][order(de_sig[de_sig$avg_log2FC > 0, ]$p_val_adj), ], 5)
      top_down <- head(de_sig[de_sig$avg_log2FC < 0, ][order(de_sig[de_sig$avg_log2FC < 0, ]$p_val_adj), ], 5)
      
      if(nrow(top_up) > 0) {
        cat("  - Top upregulated genes:", paste(rownames(top_up), collapse = ", "), "\n")
      }
      if(nrow(top_down) > 0) {
        cat("  - Top downregulated genes:", paste(rownames(top_down), collapse = ", "), "\n")
      }
      
      # Enrichment analysis
      cat("  - Running pathway enrichment...\n")
      dbs <- "GO_Biological_Process_2023"
      top_genes <- head(de_sig[order(de_sig$p_val_adj), "gene"], 100)
      
      tryCatch({
        enrichr_results <- enrichr(genes = top_genes, databases = dbs)
        go_bp <- enrichr_results[[dbs]]
        
        write.csv(go_bp, file = paste0("go_bp_gba1_dam_", tp1, "_vs_", tp2, ".csv"), row.names = FALSE)
        
        if(nrow(go_bp) > 0 && any(go_bp$Adjusted.P.value < 0.05)) {
          go_bp_sig <- go_bp[go_bp$Adjusted.P.value < 0.05, ]
          
          cat("  - Significant GO BP terms found:", nrow(go_bp_sig), "\n")
          top_5_terms <- head(go_bp_sig[, c("Term", "Adjusted.P.value", "Combined.Score")], 5)
          for(i in 1:nrow(top_5_terms)) {
            cat("    ", i, ". ", top_5_terms$Term[i], " (padj=", 
                format(top_5_terms$Adjusted.P.value[i], digits=3), ")\n", sep="")
          }
          
          # Plot top terms
          top_terms <- head(go_bp_sig[order(go_bp_sig$Adjusted.P.value), ], 10)
          p <- ggplot(top_terms, aes(x = reorder(Term, Combined.Score), y = Combined.Score)) +
            geom_col(fill = "#FF6B6B") +
            coord_flip() +
            labs(title = paste("GO BP Enrichment:", tp1, "vs", tp2, "in GBA1 DAMs"),
                 x = "GO Terms", y = "Combined Score") +
            theme_minimal(base_size = 12) +
            theme(axis.text.y = element_text(size = 10))
          
          ggsave(filename = paste0("go_bp_gba1_dam_", tp1, "_vs_", tp2, ".png"), 
                 plot = p, width = 12, height = 8, dpi = 300)
          
          cat("  - Enrichment plot saved for", tp1, "vs", tp2, "\n")
          
        } else {
          cat("  - No significant GO BP terms found for", tp1, "vs", tp2, "\n")
        }
        
      }, error = function(e) {
        cat("  - Enrichment analysis failed:", e$message, "\n")
      })
    }
    
  } else {
    cat("  - Not enough cells per timepoint (need ≥10 per group)\n")
  }
  
}, error = function(e) {
  cat("  - ERROR:", e$message, "\n")
})

# Final summary
cat("\n=== SUMMARY OF TIMEPOINT EFFECTS IN GBA1 DAMs ===\n")
n_sig <- nrow(de_timepoint_gba1[[paste0(tp1, "_vs_", tp2)]])
cat(tp1, "vs", tp2, ":", n_sig, "significant DE genes\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Output files generated:\n")
cat("- de_gba1_dam_", tp1, "_vs_", tp2, "_all.csv : All DE results\n")
cat("- de_gba1_dam_", tp1, "_vs_", tp2, "_sig.csv : Significant DE results (padj < 0.05)\n")
cat("- go_bp_gba1_dam_", tp1, "_vs_", tp2, ".csv : GO enrichment results\n")
cat("- go_bp_gba1_dam_", tp1, "_vs_", tp2, ".png : GO enrichment plot\n")



# ==================== SEPARATE GENOTYPES DE FOR TREATMENTS ====================


celltypes <- unique(PD$microglia_type)
treatments <- c("LPS", "aSyn", "debris")
dbs <- "GO_Biological_Process_2023"
DefaultAssay(PD) <- "RNA"

for (geno in c("GBA1", "KOLF")) {
  cat("\n=== Processing genotype:", geno, "===\n")
  
  PD_geno <- subset(PD, subset = Genotype == geno)
  
  for (celltype in celltypes) {
    cat("\n--- Cell type:", celltype, "---\n")
    
    subset_cells <- WhichCells(PD_geno, expression = microglia_type == celltype)
    if (length(subset_cells) == 0) next
    
    PD_sub <- subset(PD_geno, cells = subset_cells)
    Idents(PD_sub) <- "Treatment"
    
    for (treatment in treatments) {
      treat_cells <- WhichCells(PD_sub, expression = Treatment == treatment)
      ctrl_cells <- WhichCells(PD_sub, expression = Treatment == "ctrl")
      
      if (length(treat_cells) < 10 || length(ctrl_cells) < 10) next
      
      de_res <- FindMarkers(PD_sub,
                            ident.1 = treatment,
                            ident.2 = "ctrl",
                            test.use = "wilcox",
                            logfc.threshold = 0.25,
                            min.pct = 0.1,
                            verbose = FALSE)
      de_res$gene <- rownames(de_res)
      
      # Save all DE results
      write.csv(de_res, paste0("de_", geno, "_", gsub(" ", "_", celltype), "_", treatment, "_vs_ctrl_all.csv"),
                row.names = TRUE)
      
      # Significant genes
      de_sig <- de_res[de_res$p_val_adj < 0.05, ]
      if (nrow(de_sig) == 0) next
      write.csv(de_sig, paste0("de_", geno, "_", gsub(" ", "_", celltype), "_", treatment, "_vs_ctrl_sig.csv"),
                row.names = TRUE)
      
      # Enrichment
      top_genes <- head(de_sig[order(de_sig$p_val_adj), "gene"], 100)
      enrichr_results <- enrichr(genes = top_genes, databases = dbs)
      go_bp <- enrichr_results[[dbs]]
      write.csv(go_bp, paste0("go_bp_", geno, "_", gsub(" ", "_", celltype), "_", treatment, "_vs_ctrl.csv"),
                row.names = FALSE)
      
      if (nrow(go_bp) > 0 && any(go_bp$Adjusted.P.value < 0.05)) {
        top_terms <- head(go_bp[order(go_bp$Adjusted.P.value), ], 10)
        p <- ggplot(top_terms, aes(x = reorder(Term, Combined.Score), y = Combined.Score)) +
          geom_col(fill = "#6A4C93") +
          coord_flip() +
          labs(title = paste("GO BP Enrichment:", treatment, "vs ctrl in", geno, celltype),
               x = "", y = "Combined Score") +
          theme_minimal(base_size = 12)
        ggsave(paste0("go_bp_", geno, "_", gsub(" ", "_", celltype), "_", treatment, "_vs_ctrl.png"),
               plot = p, width = 8, height = 5)
      }
    }
  }
}


library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(pheatmap)
library(tibble)



# Set your folder containing CSV DE results
folder_path <- "/data/outputs/cellranger_outputs/Analysis_by_Mo_Multi/DEGs_and_EnrichmentPlots/Separated_Genotype_by_Treatment"

# Set working directory to the folder path
setwd(folder_path)

# List GO BP csv files (exclude PNGs)
go_bp_files <- list.files(pattern = "^go_bp_.*\\.csv$", full.names = FALSE)

# Check if files were found
if(length(go_bp_files) == 0) {
  stop("No GO BP CSV files found. Check your file pattern and directory.")
}

print(paste("Found", length(go_bp_files), "GO BP files:"))
print(go_bp_files)

# Read and combine GO BP csv files
go_bp_combined <- go_bp_files %>%
  set_names() %>%
  map_dfr(~ {
    cat("Reading file:", .x, "\n")
    read_csv(.x)
  }, .id = "filename")

# Parse metadata from GO BP filenames (handle hyphens in cell types)
go_bp_combined <- go_bp_combined %>%
  mutate(
    filename_trim = filename %>% str_remove("^go_bp_") %>% str_remove("\\.csv$")
  ) %>%
  # Split filename: genotype_celltype_treatment_vs_ctrl
  separate(filename_trim, 
           into = c("genotype", "remaining"), 
           sep = "_", 
           extra = "merge") %>%
  # Extract treatment (last part before _vs_ctrl)
  mutate(
    treatment = str_extract(remaining, "(aSyn|debris|LPS)(?=_vs_ctrl)"),
    # Remove treatment and _vs_ctrl to get cell type
    microglia_type = str_remove(remaining, "_(aSyn|debris|LPS)_vs_ctrl$"),
    comparison = "vs_ctrl"
  ) %>%
  dplyr::select(-remaining)

# Save combined GO BP table
write_csv(go_bp_combined, "combined_GO_BP_results.csv")

# Load data
go_data <- go_bp_combined

# Inspect the data structure
cat("Data dimensions:", dim(go_data), "\n")
cat("Column names:\n")
print(colnames(go_data))

# Create comprehensive condition identifier (cell_type + treatment)
go_data <- go_bp_combined

# Clean and verify the parsed data
cat("Data dimensions:", dim(go_data), "\n")
cat("Column names:\n")
print(colnames(go_data))

# Check parsing results
cat("Parsing verification:\n")
parsing_check <- go_data %>% 
  dplyr::select(filename, genotype, microglia_type, treatment, comparison) %>% 
  distinct() %>%
  arrange(genotype, microglia_type, treatment)
print(parsing_check)

go_data <- go_data %>%
  mutate(
    # Clean up individual components
    microglia_type_clean = ifelse(is.na(microglia_type) | microglia_type == "", "Unknown", microglia_type),
    treatment_clean = ifelse(is.na(treatment) | treatment == "", "Control", treatment),
    
    # Create cell_treatment combination for x-axis
    cell_treatment = paste(microglia_type_clean, treatment_clean, sep = "_"),
    
    # Remove trailing underscores and clean up
    cell_treatment = str_remove(cell_treatment, "_$")
  )

# Print unique combinations to see what we have
cat("Unique cell_treatment combinations:\n")
unique_combinations <- go_data %>% 
  dplyr::select(microglia_type_clean, treatment_clean, cell_treatment) %>% 
  distinct() %>%
  arrange(microglia_type_clean, treatment_clean)
print(unique_combinations)

# Check available score and term columns
score_column <- NULL
possible_score_cols <- c("Combined.Score", "Combined_Score", "NES", "enrichmentScore", 
                         "pvalue", "p_value", "padj", "p.adjust", "qvalue")

for(col in possible_score_cols) {
  if(col %in% colnames(go_data)) {
    score_column <- col
    break
  }
}

if(is.null(score_column)) {
  cat("Available columns:\n")
  print(colnames(go_data))
  stop("No recognized score column found.")
}

term_column <- NULL
possible_term_cols <- c("Term", "term", "GO_term", "Description", "pathway", "ID")

for(col in possible_term_cols) {
  if(col %in% colnames(go_data)) {
    term_column <- col
    break
  }
}

if(is.null(term_column)) {
  cat("Available columns:\n")
  print(colnames(go_data))
  stop("No recognized term column found.")
}

cat("Using score column:", score_column, "\n")
cat("Using term column:", term_column, "\n")

# Convert score to numeric and remove NAs
go_data[[score_column]] <- as.numeric(go_data[[score_column]])
go_data <- go_data %>% filter(!is.na(.data[[score_column]]))

# Get unique genotypes
genotypes <- unique(go_data$genotype)
genotypes <- genotypes[!is.na(genotypes)]
cat("Found genotypes:", paste(genotypes, collapse = ", "), "\n")

# Function to create heatmap for each genotype
create_genotype_heatmap <- function(data, geno, top_n = 20) {
  
  cat("\n=== Processing genotype:", geno, "===\n")
  
  # Filter data for this genotype
  geno_data <- data %>%
    filter(genotype == geno) %>%
    dplyr::select(all_of(c(term_column, "cell_treatment", score_column)))
  
  if(nrow(geno_data) == 0) {
    cat("No data found for genotype:", geno, "\n")
    return(NULL)
  }
  
  # Get top enriched terms (by absolute score)
  top_terms <- geno_data %>%
    group_by(.data[[term_column]]) %>%
    summarise(max_abs_score = max(abs(.data[[score_column]]), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_abs_score)) %>%
    head(top_n) %>%
    pull(.data[[term_column]])
  
  cat("Selected top", length(top_terms), "terms for", geno, "\n")
  
  # Filter for top terms and pivot to wide format
  heatmap_data <- geno_data %>%
    filter(.data[[term_column]] %in% top_terms) %>%
    # Handle duplicates by taking mean
    group_by(.data[[term_column]], cell_treatment) %>%
    summarise(score = mean(.data[[score_column]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      names_from = cell_treatment,
      values_from = score,
      values_fill = 0
    )
  
  # Convert to matrix
  heatmap_matrix <- heatmap_data %>%
    column_to_rownames(var = term_column) %>%
    as.matrix()
  
  # Remove all-zero rows and columns
  heatmap_matrix <- heatmap_matrix[rowSums(abs(heatmap_matrix)) > 0, , drop = FALSE]
  heatmap_matrix <- heatmap_matrix[, colSums(abs(heatmap_matrix)) > 0, drop = FALSE]
  
  cat("Final matrix dimensions for", geno, ":", dim(heatmap_matrix), "\n")
  
  if(nrow(heatmap_matrix) == 0 || ncol(heatmap_matrix) == 0) {
    cat("No data to plot for genotype:", geno, "\n")
    return(NULL)
  }
  
  # Create heatmap (display in plot, don't save)
  tryCatch({
    # Create column annotations for treatments
    col_data <- data.frame(
      cell_treatment = colnames(heatmap_matrix),
      stringsAsFactors = FALSE
    ) %>%
      separate(cell_treatment, into = c("cell_type", "treatment"), sep = "_", remove = FALSE, extra = "merge") %>%
      column_to_rownames("cell_treatment")
    
    # Create annotation colors
    unique_treatments <- unique(col_data$treatment)
    unique_celltypes <- unique(col_data$cell_type)
    
    treatment_colors <- rainbow(length(unique_treatments))
    names(treatment_colors) <- unique_treatments
    
    celltype_colors <- rainbow(length(unique_celltypes))
    names(celltype_colors) <- unique_celltypes
    
    ann_colors <- list(
      treatment = treatment_colors,
      cell_type = celltype_colors
    )
    
    p <- pheatmap(
      heatmap_matrix,
      scale = "row",  # Scale by GO terms (rows)
      clustering_method = "complete",
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 8,
      fontsize_col = 8,
      angle_col = 45,
      main = paste("GO BP Enrichment -", geno, "Genotype"),
      annotation_col = col_data,
      annotation_colors = ann_colors
    )
    
    cat("Heatmap displayed for genotype:", geno, "\n")
    return(p)
    
  }, error = function(e) {
    cat("Error creating annotated heatmap for", geno, ":", e$message, "\n")
    
    # Try simplified version without annotations
    cat("Trying simplified heatmap for", geno, "...\n")
    p_simple <- pheatmap(
      heatmap_matrix,
      scale = "row",
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 8,
      fontsize_col = 8,
      angle_col = 45,
      main = paste("GO BP Enrichment -", geno, "Genotype (Simplified)")
    )
    
    cat("Simplified heatmap displayed for genotype:", geno, "\n")
    return(p_simple)
  })
}

# Create heatmaps for each genotype
# Create heatmaps for each genotype
heatmap_plots <- list()
for(geno in genotypes) {
  cat("\n" , paste(rep("=", 50), collapse=""), "\n")
  cat("Creating heatmap for genotype:", geno, "\n")
  cat(paste(rep("=", 50), collapse=""), "\n")
  
  heatmap_plots[[geno]] <- create_genotype_heatmap(go_data, geno, top_n = 20)
}

# Display summary
cat("\n=== SUMMARY ===\n")
cat("Created heatmaps for genotypes:", paste(names(heatmap_plots), collapse = ", "), "\n")

# Function to redisplay any specific genotype heatmap
display_genotype_heatmap <- function(genotype_name) {
  if(genotype_name %in% names(heatmap_plots)) {
    create_genotype_heatmap(go_data, genotype_name, top_n = 20)
    cat("Displayed heatmap for:", genotype_name, "\n")
  } else {
    cat("Genotype", genotype_name, "not found. Available:", paste(names(heatmap_plots), collapse = ", "), "\n")
  }
}

# Display both heatmaps
cat("\n=== DISPLAYING ALL HEATMAPS ===\n")
for(geno in names(heatmap_plots)) {
  cat("Displaying heatmap for:", geno, "\n")
  display_genotype_heatmap(geno)
  cat("------------------------\n")
}

cat("\nTo redisplay a specific heatmap, use:\n")
cat("display_genotype_heatmap('GBA1')\n")
cat("display_genotype_heatmap('KOLF')\n")


# To see GBA1 heatmap
display_genotype_heatmap('GBA1')

# To see KOLF heatmap  
display_genotype_heatmap('KOLF')




# ====================  TREATMENT-SPECIFIC GENOTYPE ANALYSIS (SEURAT v5) ====================

library(Seurat)
library(dplyr)
library(ggplot2)

cat("=== PREPARING SEURAT OBJECT ===\n")

# 0. Align metadata with true cell names
cat("Fixing metadata rownames...\n")
rownames(PD@meta.data) <- Cells(PD)

# 1. Clear graphs (optional, prevents conflicts)
cat("Removing graphs...\n")
PD@graphs <- list()

# 2. Use RNA assay for DE
cat("Setting default assay to RNA...\n")
DefaultAssay(PD) <- "RNA"
cat("Default assay:", DefaultAssay(PD), "\n")

# 3. Verify RNA assay layers
rna_counts <- PD[["RNA"]]@layers$counts
rna_data   <- PD[["RNA"]]@layers$data
cat("RNA counts:", dim(rna_counts), "\n")
cat("RNA data:", dim(rna_data), "\n")

# ==================== ANALYSIS ====================

cat("\n=== STARTING TREATMENT-SPECIFIC DE ANALYSIS ===\n")

treatments <- c("ctrl", "LPS", "aSyn", "debris")
celltypes  <- unique(PD$microglia_type)

de_results_by_treatment <- list()

for (treatment in treatments) {
  cat("\n", paste(rep("=", 50), collapse=""), "\n")
  cat("ANALYZING TREATMENT:", treatment, "\n")
  cat(paste(rep("=", 50), collapse=""), "\n")
  
  # Subset by treatment
  PD_treatment <- subset(PD, subset = Treatment == treatment)
  cat("Cells in", treatment, "condition:", ncol(PD_treatment), "\n")
  
  # Genotype counts
  geno_counts <- table(PD_treatment$Genotype)
  print(geno_counts)
  
  if (!all(c("GBA1","KOLF") %in% names(geno_counts))) {
    cat("Missing genotypes - skipping\n")
    next
  }
  if (any(geno_counts < 20)) {
    cat("Insufficient cells - skipping\n")
    next
  }
  
  Idents(PD_treatment) <- "Genotype"
  de_results_by_treatment[[treatment]] <- list()
  
  for (celltype in celltypes) {
    cat("\n--- Analyzing", celltype, "in", treatment, "---\n")
    
    # Subset by celltype
    PD_sub <- subset(PD_treatment, subset = microglia_type == celltype)
    if (ncol(PD_sub) < 20) {
      cat("Too few", celltype, "cells (", ncol(PD_sub), ") - skipping\n")
      next
    }
    
    geno_ct <- table(PD_sub$Genotype)
    print(geno_ct)
    
    if (!all(c("GBA1","KOLF") %in% names(geno_ct))) {
      cat("Missing genotypes - skipping\n")
      next
    }
    if (any(geno_ct < 10)) {
      cat("Too few cells per genotype - skipping\n")
      next
    }
    
    Idents(PD_sub) <- "Genotype"
    
    cat("Running DE...\n")
    de_result <- FindMarkers(PD_sub,
                             ident.1 = "GBA1",
                             ident.2 = "KOLF",
                             assay = "RNA",
                             test.use = "wilcox",
                             logfc.threshold = 0.25,
                             min.pct = 0.1,
                             verbose = FALSE)
    
    de_result$gene <- rownames(de_result)
    de_result$treatment <- treatment
    de_result$celltype  <- celltype
    
    # Save results
    filename_prefix <- paste0("de_", treatment, "_", gsub(" ", "_", celltype), "_GBA1_vs_KOLF")
    write.csv(de_result, paste0(filename_prefix, "_all.csv"), row.names = TRUE)
    
    # Filter significant
    de_sig <- de_result[de_result$p_val_adj < 0.05, ]
    de_results_by_treatment[[treatment]][[celltype]] <- de_sig
    
    cat("Results:", nrow(de_result), "total,", nrow(de_sig), "significant\n")
    if (nrow(de_sig) > 0) {
      write.csv(de_sig, paste0(filename_prefix, "_sig.csv"), row.names = TRUE)
      cat("- Upregulated in GBA1:", sum(de_sig$avg_log2FC > 0), "\n")
      cat("- Downregulated in GBA1:", sum(de_sig$avg_log2FC < 0), "\n")
    }
  }
}

# ==================== SUMMARY ====================

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("ANALYSIS SUMMARY\n")
cat(paste(rep("=", 60), collapse=""), "\n")

summary_data <- data.frame()
total_comparisons <- 0
successful_comparisons <- 0

for (treatment in names(de_results_by_treatment)) {
  for (celltype in names(de_results_by_treatment[[treatment]])) {
    n_sig <- nrow(de_results_by_treatment[[treatment]][[celltype]])
    if (n_sig > 0) {
      n_up   <- sum(de_results_by_treatment[[treatment]][[celltype]]$avg_log2FC > 0)
      n_down <- sum(de_results_by_treatment[[treatment]][[celltype]]$avg_log2FC < 0)
      
      summary_data <- rbind(summary_data, data.frame(
        Treatment = treatment,
        CellType  = celltype,
        Significant_DE_genes = n_sig,
        Upregulated_GBA1 = n_up,
        Downregulated_GBA1 = n_down
      ))
      successful_comparisons <- successful_comparisons + 1
    }
    total_comparisons <- total_comparisons + 1
  }
}

cat("Successful:", successful_comparisons, "of", total_comparisons, "\n")

if (nrow(summary_data) > 0) {
  print(summary_data)
  write.csv(summary_data, "summary_treatment_genotype_effects.csv", row.names = FALSE)
} else {
  cat("No significant DE found.\n")
}

cat("\nOutput files:\n")
cat("- de_[treatment]_[celltype]_GBA1_vs_KOLF_all.csv\n")
cat("- de_[treatment]_[celltype]_GBA1_vs_KOLF_sig.csv\n")
cat("- summary_treatment_genotype_effects.csv\n")



# ==================== ENRICHR ANALYSIS (FILTERED) ====================

library(enrichR)

# List available databases
dbs <- listEnrichrDbs()
if(!"GO_Biological_Process_2023" %in% dbs$libraryName) {
  stop("GO_Biological_Process_2023 not available in your EnrichR installation")
}

cat("\n=== RUNNING ENRICHR (GO BP 2023, FDR < 0.05) ===\n")

all_enrich_results <- data.frame()

for (treatment in names(de_results_by_treatment)) {
  for (celltype in names(de_results_by_treatment[[treatment]])) {
    
    de_sig <- de_results_by_treatment[[treatment]][[celltype]]
    if (is.null(de_sig) || nrow(de_sig) == 0) next
    
    # Take significant DE genes
    sig_genes <- de_sig$gene
    if (length(sig_genes) < 5) {
      cat("Too few genes for", treatment, "-", celltype, "skipping enrichR\n")
      next
    }
    
    # Run enrichR
    enr <- enrichr(sig_genes, databases = c("GO_Biological_Process_2023"))
    go_bp <- enr[["GO_Biological_Process_2023"]]
    
    if (!is.null(go_bp) && nrow(go_bp) > 0) {
      # Filter by adjusted p-value
      go_bp <- go_bp[go_bp$Adjusted.P.value < 0.05, ]
      
      if (nrow(go_bp) == 0) next  # skip if nothing passes filter
      
      # Optional: keep top 10 terms by Combined Score
      go_bp <- go_bp[order(-go_bp$Combined.Score), ]
      go_bp <- head(go_bp, 10)
      
      # Add metadata
      go_bp$treatment <- treatment
      go_bp$celltype  <- celltype
      all_enrich_results <- rbind(all_enrich_results, go_bp)
      
      cat("EnrichR done for", treatment, "-", celltype, "\n")
    }
  }
}

# Save filtered enrichment results
if (nrow(all_enrich_results) > 0) {
  write.csv(all_enrich_results, "enrichr_GO_BP2023_filtered_results.csv", row.names = FALSE)
  cat("\nFiltered enrichment results saved to enrichr_GO_BP2023_filtered_results.csv\n")
} else {
  cat("\nNo enrichment results passed the FDR < 0.05 filter.\n")
}






# ====================  TREATMENT-SPECIFIC TIMEPOINT ANALYSIS (SEURAT v5) ====================

library(Seurat)
library(dplyr)
library(ggplot2)

cat("=== PREPARING SEURAT OBJECT ===\n")

# 0. Align metadata with true cell names
cat("Fixing metadata rownames...\n")
rownames(PD@meta.data) <- Cells(PD)

# 1. Clear graphs (optional, prevents conflicts)
cat("Removing graphs...\n")
PD@graphs <- list()

# 2. Use RNA assay for DE
cat("Setting default assay to RNA...\n")
DefaultAssay(PD) <- "RNA"
cat("Default assay:", DefaultAssay(PD), "\n")

# 3. Verify RNA assay layers
rna_counts <- PD[["RNA"]]@layers$counts
rna_data   <- PD[["RNA"]]@layers$data
cat("RNA counts:", dim(rna_counts), "\n")
cat("RNA data:", dim(rna_data), "\n")

# ==================== ANALYSIS ====================

cat("\n=== STARTING TREATMENT-SPECIFIC DE ANALYSIS (24h vs 7d) ===\n")

treatments <- c("ctrl", "LPS", "aSyn", "debris")
celltypes  <- unique(PD$microglia_type)

de_results_by_treatment <- list()

for (treatment in treatments) {
  cat("\n", paste(rep("=", 50), collapse=""), "\n")
  cat("ANALYZING TREATMENT:", treatment, "\n")
  cat(paste(rep("=", 50), collapse=""), "\n")
  
  # Subset by treatment
  PD_treatment <- subset(PD, subset = Treatment == treatment)
  cat("Cells in", treatment, "condition:", ncol(PD_treatment), "\n")
  
  # Timepoint counts
  time_counts <- table(PD_treatment$Timepoint)
  print(time_counts)
  
  if (!all(c("24h","7d") %in% names(time_counts))) {
    cat("Missing timepoints - skipping\n")
    next
  }
  if (any(time_counts < 20)) {
    cat("Insufficient cells - skipping\n")
    next
  }
  
  Idents(PD_treatment) <- "Timepoint"
  de_results_by_treatment[[treatment]] <- list()
  
  for (celltype in celltypes) {
    cat("\n--- Analyzing", celltype, "in", treatment, "---\n")
    
    # Subset by celltype
    PD_sub <- subset(PD_treatment, subset = microglia_type == celltype)
    if (ncol(PD_sub) < 20) {
      cat("Too few", celltype, "cells (", ncol(PD_sub), ") - skipping\n")
      next
    }
    
    time_ct <- table(PD_sub$Timepoint)
    print(time_ct)
    
    if (!all(c("24h","7d") %in% names(time_ct))) {
      cat("Missing timepoints - skipping\n")
      next
    }
    if (any(time_ct < 10)) {
      cat("Too few cells per timepoint - skipping\n")
      next
    }
    
    Idents(PD_sub) <- "Timepoint"
    
    cat("Running DE...\n")
    de_result <- FindMarkers(PD_sub,
                             ident.1 = "24h",
                             ident.2 = "7d",
                             assay = "RNA",
                             test.use = "wilcox",
                             logfc.threshold = 0.25,
                             min.pct = 0.1,
                             verbose = FALSE)
    
    de_result$gene <- rownames(de_result)
    de_result$treatment <- treatment
    de_result$celltype  <- celltype
    
    # Save results
    filename_prefix <- paste0("de_", treatment, "_", gsub(" ", "_", celltype), "_24h_vs_7d")
    write.csv(de_result, paste0(filename_prefix, "_all.csv"), row.names = TRUE)
    
    # Filter significant
    de_sig <- de_result[de_result$p_val_adj < 0.05, ]
    de_results_by_treatment[[treatment]][[celltype]] <- de_sig
    
    cat("Results:", nrow(de_result), "total,", nrow(de_sig), "significant\n")
    if (nrow(de_sig) > 0) {
      write.csv(de_sig, paste0(filename_prefix, "_sig.csv"), row.names = TRUE)
      cat("- Upregulated in 24h:", sum(de_sig$avg_log2FC > 0), "\n")
      cat("- Downregulated in 24h:", sum(de_sig$avg_log2FC < 0), "\n")
    }
  }
}

# ==================== SUMMARY ====================

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("ANALYSIS SUMMARY (24h vs 7d)\n")
cat(paste(rep("=", 60), collapse=""), "\n")

summary_data <- data.frame()
total_comparisons <- 0
successful_comparisons <- 0

for (treatment in names(de_results_by_treatment)) {
  for (celltype in names(de_results_by_treatment[[treatment]])) {
    n_sig <- nrow(de_results_by_treatment[[treatment]][[celltype]])
    if (n_sig > 0) {
      n_up   <- sum(de_results_by_treatment[[treatment]][[celltype]]$avg_log2FC > 0)
      n_down <- sum(de_results_by_treatment[[treatment]][[celltype]]$avg_log2FC < 0)
      
      summary_data <- rbind(summary_data, data.frame(
        Treatment = treatment,
        CellType  = celltype,
        Significant_DE_genes = n_sig,
        Upregulated_24h = n_up,
        Downregulated_24h = n_down
      ))
      successful_comparisons <- successful_comparisons + 1
    }
    total_comparisons <- total_comparisons + 1
  }
}

cat("Successful:", successful_comparisons, "of", total_comparisons, "\n")

if (nrow(summary_data) > 0) {
  print(summary_data)
  write.csv(summary_data, "summary_treatment_timepoint_effects.csv", row.names = FALSE)
} else {
  cat("No significant DE found.\n")
}

cat("\nOutput files:\n")
cat("- de_[treatment]_[celltype]_24h_vs_7d_all.csv\n")
cat("- de_[treatment]_[celltype]_24h_vs_7d_sig.csv\n")
cat("- summary_treatment_timepoint_effects.csv\n")









