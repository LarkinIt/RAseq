# Load required libraries
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

cat("Set output and loading seurat...\n")
seurat_f <- "/ix/djishnu/Aaron_F/RA/Results/EMATB8322/integrated_samples/integrated_seurat_DimReduced_and_Clustered.rds"
# Set your data directory path
out_dir <- "/ix/djishnu/Aaron_F/RA/Results/EMATB8322/integrated_samples/out"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)

seurat_obj <- LoadSeuratRds(seurat_f)
cat("Done!\n")

# cat("Clustering and Dim Reduction...\n")
# # Step 1: Run UMAP based on PCA embeddings using first 12 PCs (as specified in paper)
# seurat_obj <- RunUMAP(seurat_obj, 
#                       reduction = "pca", 
#                       dims = 1:12,
#                       verbose = TRUE)

# # Step 2: Perform clustering
# # Find neighbors using the first 12 PCs
# seurat_obj <- FindNeighbors(seurat_obj, 
#                             reduction = "pca", 
#                             dims = 1:12)

# # Find clusters
# seurat_obj <- FindClusters(seurat_obj, 
#                           resolution = 0.5)  # Adjust resolution as needed
# cat("Done!\n")


cat("Get cluster markers...\n")
# Find markers for all clusters using MAST test
# Using non-batch corrected counts as recommended
DefaultAssay(seurat_obj) <- "RNA"  # Use original counts, not integrated
# First, join the layers in your RNA assay
seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])

cluster_markers <- FindAllMarkers(
  seurat_obj,
  test.use = "MAST",
  min.pct = 0.4,
  logfc.threshold = 0.1,
  only.pos = FALSE,
  verbose = TRUE
)

# Filter for significantly DE genes (adjusted p < 0.05)
significant_markers <- cluster_markers %>%
  filter(p_val_adj < 0.05) %>%
  arrange(cluster, desc(avg_log2FC))

cat("Done!\n")

cat("Saving...\n")
# save
write.csv(significant_markers, "all_clusterMarkers_DEres.csv")
write.csv(significant_markers, "sig_clusterMarkers_DEres.csv")
cat("Done!\n")