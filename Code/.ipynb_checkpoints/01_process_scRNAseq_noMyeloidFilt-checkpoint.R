# Load required libraries
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(dplyr)

# load patient info for samples to process and qc threshold to use
patient_info <- read.csv("/ihome/rgottschalk/cil8/RAseq/data/Caroline/Data/EMATB8322/patient_info_das28.csv")

# Initialize summary dataframe
processing_summary <- data.frame(
  sample_name = character(),
  preQC_cells = numeric(),
  postQC_cells = numeric(),
  percent_retained_afterQC = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq(dim(patient_info)[1])){
  # qc metrics to use for sample processing:
  sample_name <- patient_info[i,"Sample_name"]
  min_cells <- patient_info[i,"min_cells"]
  min_features <- patient_info[i, "min_features"]
  max_features <- patient_info[i, "max_features"]
  max_mt_exp <- patient_info[i, "max_mt_exp"]
  pass_qc <- patient_info[i, "pass_qc"]
  print(paste("Processing sample", i, "of", nrow(patient_info), ":", sample_name))

  # Set your data directory path
  data_dir <- paste("/ihome/rgottschalk/cil8/RAseq/data/Caroline/Data/EMATB8322", sample_name, sep="/")
  setwd(data_dir)

  # Read the 10X format data
  # The function expects barcode.tsv, features.tsv (or genes.tsv), and matrix.mtx
  data_10x <- Read10X(data.dir = data_dir)

  # Create Seurat object with min.cells = 5 (genes expressed in at least 5 cells)
  seurat_obj <- CreateSeuratObject(
    counts = data_10x,
    min.cells = min_cells
  )

  # Display initial object statistics
  print(paste("Initial dimensions:", dim(seurat_obj)[1], "genes x", dim(seurat_obj)[2], "cells"))
  preQC_count <- ncol(seurat_obj)

  # Calculate QC metrics
  # Calculate mitochondrial gene percentages
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

  # Calculate ribosomal gene percentages (optional but useful)
  seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")

  # Display QC metrics summary
  print("Pre-QC Metrics Summary:")

  
  # Visualize QC metrics
  # Create violin plots for QC metrics
  #pdf("preQC_plots.pdf")
  qc_plots <- VlnPlot(seurat_obj, 
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                  combine = FALSE
  )
  #print(qc_plots)

  # Create scatter plots to visualize relationships between metrics
  scatter1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  scatter2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  #print(scatter1)
  #print(scatter2)
  #dev.off()

  # Filter cells with min/max number features expressed and below mitochondrial expr threshold
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > min_features)
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA < max_features)
  seurat_obj <- subset(seurat_obj, subset = percent.mt < max_mt_exp)

  print(paste("After QC filtering:", dim(seurat_obj)[1], "genes x", dim(seurat_obj)[2], "cells"))
  postQC_count <- ncol(seurat_obj)
  percent_retained <- round(100 * postQC_count / preQC_count, 1)

  # Visualize QC metrics
  # Create violin plots for QC metrics
  pdf("postQC_plots.pdf")
  qc_plots <- VlnPlot(seurat_obj, 
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                  combine = FALSE
  )
  print(qc_plots)

  # Create scatter plots to visualize relationships between metrics
  scatter1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  scatter2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(scatter1)
  print(scatter2)
  dev.off()

  # Normalize data
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

  # Find variable genes
  # Identify the top 2,000 most variable genes
  seurat_obj <- FindVariableFeatures(seurat_obj, 
                                      selection.method = "vst", 
                                      nfeatures = 2000)

  # Display the top 10 most variable genes
  top10 <- head(VariableFeatures(seurat_obj), 10)

  pdf("variable_feature_plot.pdf")
  # Plot variable features
  variable_plot <- VariableFeaturePlot(seurat_obj)
  print(variable_plot)

  # Label the top 10 most variable genes on the plot
  labeled_plot <- LabelPoints(plot = variable_plot, points = top10, repel = TRUE)
  print(labeled_plot)
  dev.off()

  # Save the Seurat object
  # saveRDS(seurat_obj, file = paste(sample_name, "processed_seurat.rds", sep = "_"))

  # save metadata table:
  seurat_obj@meta.data$barcode <- rownames(seurat_obj@meta.data)
  write.csv(seurat_obj@meta.data, file='processed_metadata.csv', quote=F, row.names=F)
  
  # write expression counts matrix
  counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='counts')
  writeMM(counts_matrix, file='processed_counts.mtx')
  
  # write gene names
  write.table(data.frame('gene'=rownames(counts_matrix)),file='processed_gene_names.csv', quote=F,row.names=F,col.names=F)
  processing_summary <- rbind(
    processing_summary,
    data.frame(
      sample_name = sample_name,
      preQC_cells = preQC_count,
      postQC_cells = postQC_count,
      percent_retained_afterQC = percent_retained
    )
  )

  print("done!")

}

# Save the processing summary
write.csv(processing_summary, file = "~/RAseq/data/Caroline/Data/EMATB8322/processed_data/processing_summary.csv", 
          row.names = FALSE)