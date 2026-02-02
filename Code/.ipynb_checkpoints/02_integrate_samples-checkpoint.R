# Load required libraries
library(Seurat)
library(Matrix)
library(dplyr)

#### OUTPUT ####
out_dir <- "/ix/djishnu/Aaron_F/RA/Results/EMATB8322/integrated_samples"
dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

#### READ IN ALL SAMPLES ####
patient_info <- read.csv("/ix/djishnu/Aaron_F/RA/Data/EMATB8322/patient_info_das28.csv")
# Set your data directory path
data_dir <- "/ix/djishnu/Aaron_F/RA/Data/EMATB8322/processed_data"

seurat_objects_list <- list()
for (i in seq(dim(patient_info)[1])){
  sample_name <- patient_info[i,"Sample_name"]
  pass_qc <- patient_info[i, "pass_qc"]
  
  # Only process samples that passed QCa
  if (pass_qc == TRUE) {
    sample_f <- paste(data_dir, sample_name, paste(sample_name, "processed_seurat.rds", sep="_"),sep="/")
    
    # Add error checking for file existence
    if (file.exists(sample_f)) {
      seurat_obj <- LoadSeuratRds(sample_f)
      
      # Add sample name as metadata for tracking
      seurat_obj$sample_id <- sample_name
      seurat_obj$batch <- sample_name
      
      seurat_objects_list[[sample_name]] <- seurat_obj
      cat("Loaded sample:", sample_name, "\n")
    } else {
      cat("Warning: File not found for sample:", sample_name, "\n")
    }
  }
}

# Check how many samples were loaded
cat("Total samples loaded:", length(seurat_objects_list), "\n")

#### INTEGRATE SAMPLES ####
# Find integration anchors
cat("Running FindIntegrationAnchors...\n")
anchors <- FindIntegrationAnchors(object.list = seurat_objects_list)
cat("done!")

# Integrate using all common genes
cat("Running IntegrateData...\n")
# common_genes <- Reduce(intersect, lapply(seurat_objects_list, rownames))
# integrated_seurat <- IntegrateData(anchorset = anchors, 
#                                   features.to.integrate = common_genes)
integrated_seurat <- IntegrateData(anchorset = anchors)
cat("done!")

# Set integrated assay as the default assay
DefaultAssay(integrated_seurat) <- "integrated"

# Scale data and run PCA (as mentioned in your paper)
integrated_seurat <- ScaleData(integrated_seurat, verbose = FALSE)
integrated_seurat <- RunPCA(integrated_seurat, npcs = 30, verbose = FALSE)

#### SAVE ####
saveRDS(integrated_seurat, file = "integrated_seurat.rds")

# Save some summary info
cat("Final integrated object:\n")
print(integrated_seurat)
cat("Sample distribution:\n")
print(table(integrated_seurat$sample_id))