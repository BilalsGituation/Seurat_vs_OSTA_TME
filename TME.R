# Setup: Import Packages and Query Dataset ####

# My usual options for seurat analysis
options(timeout=10000000) # no timeouts
options(future.globals.maxSize = 1e12) # allow large variables

# pkg list contains some leftovers from OSTA
pacman::p_load(OSTA.data, # ebook data
               tidyverse, 
               Seurat, 
               patchwork,
               visiumStitched, # Interoperation
               CARDspa, # CARD deconvolution for future version
               pheatmap,
               SpatialExperiment, # read in data
               VisiumIO, # for reading in a spatial dataset
               hdf5r, # for reading in a spatial dataset
               devtools, # install remote pkgs
               RNiftyReg, # Quantify img similarity
               png, # visual comparison
               pals, # colour scheme
               caret#, # confusion matrix
               
)
#devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
library(spacexr) # RCTD deconvolution for future version


# retrieve dataset from OSF repository
id <- "Visium_HumanBreast_Janesick"
pa <- OSTA.data_load(id)
dir.create(td <- tempfile())
unzip(pa, exdir=td)

# Make Seurat obj out of Query Dataset ####

visium_breast <- 
  Load10X_Spatial(
    file.path(td, "outs"),
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "slice1",
    bin.size = NULL,
    filter.matrix = TRUE,
    to.upper = FALSE,
    image = NULL,
    image.name = "tissue_lowres_image.png",
    segmentation.type = NULL,
    compact = TRUE
  )

visium_breast <- UpdateSeuratObject(visium_breast)

# Ground truth (spot-level type annotation) for later ####
df <- read.csv(file.path(td, "annotation.csv"))
cs <- match(colnames(visium_breast), df$Barcode)
visium_breast$anno <- factor(df$Annotation[cs])

# Get QC metrics ####

# I used the QC values in OSTA so I could benchmark the downstream outputs against theirs

# Num RNA molecules per spot
plot1 <- VlnPlot(visium_breast, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(visium_breast, features = "nCount_Spatial", pt.size.factor = 3) + theme(legend.position = "right")
wrap_plots(plot1, plot2) 

# num unique genes per spot
plot3 <- VlnPlot(visium_breast, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
wrap_plots(plot3, plot2)

# calculate mitochondrial genes per spot and plot them
visium_breast[["percent.mt"]] <- PercentageFeatureSet(visium_breast, pattern = "^MT-")   
plot4 <- VlnPlot(visium_breast, features = "percent.mt", pt.size = 0.1) + NoLegend()
wrap_plots(plot4, plot2)

# Remove Spots that fail QC ####

dim(visium_breast)
# 1000 RNA molecules per spot is considered low for this dataset
# mito counts look ok
# ignore the warnings. 30 spots should be removed
visium_breast <- subset(visium_breast, 
                        subset = nCount_Spatial > 1000 # & percent.mt < 5 
)
dim(visium_breast)

# SCTransform, cluster Query dataset ####

visium_breast <- SCTransform(visium_breast, assay = "Spatial")
visium_breast <- visium_breast %>% RunPCA(assay = "SCT", verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(reduction = "pca", dims = 1:30)

# Visualise clusters ####

# p1 = visualise clustering in UMAP space, p2 is visualised as image overlay
p1 <- UMAPPlot(visium_breast, label = TRUE) + ggtitle("UMAP plot")
p2 <- SpatialDimPlot(visium_breast, label = TRUE, label.size = 5, pt.size.factor = 3) + ggtitle("UMAP clusters over Image")
p3 <- SpatialDimPlot(visium_breast, label = TRUE, label.size = 5, pt.size.factor = 3, image.alpha = 0) + ggtitle("UMAP clusters - No Image")
p1 + p2 + p3

# Check clusters against ground truth ####

Idents(visium_breast) <- "anno"
p1 <- SpatialDimPlot(visium_breast, label = TRUE, label.size = 0, pt.size.factor = 3) + ggtitle("Ground Truth with img")
p2 <- SpatialDimPlot(visium_breast, label = TRUE, image.alpha = 0, label.size = 4, pt.size.factor = 3) + ggtitle("Ground Truth without img")
p1+p2+p3

# Import, prepare deconvolution reference ####

# We are provided with single cell reference data from 10X for deconvolution 
# using the algorithms FindTransferAnchors, RCTD and CARD
# This is the CHROMIUM set, the experimental data is the VISIUM set
# retrieve dataset from OSF repo
id <- "Chromium_HumanBreast_Janesick"
pa <- OSTA.data_load(id)
dir.create(td <- tempfile())
unzip(pa, exdir=td)

chromium_breast <- Read10X_h5(file.path(td, "filtered_feature_bc_matrix.h5"), 
                              use.names = TRUE, unique.features = TRUE)
chromium_metadata <- read.csv(file.path(td, "cell_metadata.csv"), row.names = 1)

# Streamline demo by removing mixture spots,
# consolidating and abbreviating metadata columns
chromium_ann <- chromium_metadata$Annotation
chromium_ann[grepl("Hyb", chromium_ann)] <- NA
metadata_short <- c(
  "B Cell"="B", "T Cell"="T", "Mac"="macro", "Mast"="mast", 
  "DCs"="dendritic", "Peri"="perivas", "End"="endo", 
  "Str"="stromal", "Inv"="tumor", "Myo"="myoepi")
for (pattern_name in names(metadata_short)) {
  chromium_ann[grep(pattern_name, chromium_ann)] <- metadata_short[pattern_name]
}
chromium_ann <- gsub("\\s", "", chromium_ann)
chromium_metadata$ShortName <- chromium_ann

# Are all rownames and colnames equal?
all(rownames(chromium_metadata) == colnames(chromium_breast))

# Reorder metadata to exactly match cell barcodes
#chromium_metadata <- chromium_metadata[colnames(chromium_breast),]

chromium_breast <- CreateSeuratObject(counts = chromium_breast, meta.data = chromium_metadata)
chromium_breast <- UpdateSeuratObject(chromium_breast)

Idents(chromium_breast) <- "ShortName"

chromium_breast <- SCTransform(chromium_breast, assay = "RNA")
chromium_breast <- chromium_breast %>% RunPCA(assay = "SCT", verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(reduction = "pca", dims = 1:30)

# Deconvolute, ID cells with FindTransferAnchors() ####

anchors <- FindTransferAnchors(reference = chromium_breast, 
                               query = visium_breast, 
                               normalization.method = "SCT",
                               npcs = 50)

predictions.assay <- TransferData(anchorset = anchors, 
                                  refdata = chromium_breast$ShortName, 
                                  prediction.assay = TRUE,
                                  weight.reduction = visium_breast[["pca"]], 
                                  dims = 1:50)

visium_breast[["predictions"]] <- predictions.assay

DefaultAssay(visium_breast) <- "predictions"

# Visualise Cell type predictions ####

cell.types.of.interest <- c("B", "DCIS1", "DCIS2", "T", "dendritic", "endo", "macro", "mast", "myoepi", "perivas", "stromal", "tumor")

SpatialFeaturePlot(visium_breast, features = cell.types.of.interest, alpha = c(0.1, 1), pt.size.factor = 3)

SpatialFeaturePlot(visium_breast, 
                   features = cell.types.of.interest, 
                   image.alpha = 0, 
                   pt.size.factor = 3.5) & 
  scale_fill_gradientn(colors = pals::jet(), 
                       limits = c(0, 1), 
                       oob = scales::squish) 

# Calculate match rate for TME elements ####

# Extract predicted cell type (highest score)
visium_breast$predicted_id <- rownames(visium_breast[["predictions"]])[
  apply(visium_breast[["predictions"]]@data, 2, which.max)
]

# Define clear, unambiguous mappings only
clear_mapping <- c(
  "DCIS #1" = "DCIS1",
  "DCIS #2" = "DCIS2",
  "invasive" = "tumor",
  "stromal" = "stromal"
)

# Subset to spots with clear ground truth
visium_clear <- visium_breast[, visium_breast$anno %in% names(clear_mapping)]

# Map ground truth to prediction names
# (!!! captures the name mapping scheme)
visium_clear$anno_mapped <- recode(visium_clear$anno, !!!clear_mapping)

# Extract predictions for these spots
visium_clear$predicted_id <- rownames(visium_clear[["predictions"]])[
  apply(visium_clear[["predictions"]]@data, 2, which.max)
]

# Calculate match rate per cell type
match_rate <- sapply(unique(visium_clear$anno_mapped), function(celltype) {
  spots <- which(visium_clear$anno_mapped == celltype)
  correct <- sum(visium_clear$predicted_id[spots] == celltype, na.rm = TRUE)
  return(correct / length(spots) * 100)
})

# Print match rates
for(celltype in names(match_rate)) {
  total <- sum(visium_clear$anno_mapped == celltype)
  correct <- sum(visium_clear$predicted_id[visium_clear$anno_mapped == celltype] == celltype)
  cat(sprintf("%-10s: %5.1f%% (%d/%d)\n", 
              celltype, match_rate[celltype], correct, total))
}

# Print success rate and spot total
overall_accuracy <- sum(visium_clear$predicted_id == visium_clear$anno_mapped) / ncol(visium_clear) * 100
cat(sprintf("\nOverall Accuracy: %.1f%%\n", overall_accuracy))
cat(sprintf("Total spots: %d\n", ncol(visium_clear)))

# Confusion matrix of clear, unambiguous mappings
cm_clear <- table(Predicted = visium_clear$predicted_id, 
                  Ground_Truth = visium_clear$anno_mapped)
# Only keep non-empty columns
cm_clear <- cm_clear[, colSums(cm_clear) > 0, drop = FALSE]
cm_clear <- cm_clear[rowSums(cm_clear) > 0, , drop = FALSE]
print(cm_clear)

# Column percentages for heatmap
col_pct <- prop.table(cm_clear, 2) * 100

pheatmap(col_pct,
         main = paste("Deconvolution Performance - Clear Cell Types\n(Overall Accuracy:", 
                      round(overall_accuracy, 1), "%)"),
         color = colorRampPalette(c("white", "#2E86C1"))(50),
         display_numbers = TRUE,
         number_format = "%.1f",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         angle_col = 45)

# Quantify similarity between Seurat and Bioconductor deconvolutions ####

# Using NMI and NiftyReg

# Read images in.
# You can replace data/seurat_FindTransferAnchors_img.png with your own image
seurat_FindTransferAnchors_img <- readPNG("data/seurat_FindTransferAnchors_img.png")
OSTA_RCTD_img <- readPNG("data/Screenshot_2026-06-10_18-47-50.png")

# NMI before NiftyReg affine registrations
print(paste("NMI before:", similarity(seurat_FindTransferAnchors_img, OSTA_RCTD_img, interpolation = 3L,
                                      threads = getOption("RNiftyReg.threads")))) 

# Array with affine registrations (see readme)
nmi_registered <- niftyreg(seurat_FindTransferAnchors_img, 
                           OSTA_RCTD_img, scope = "affine",
                           interpolation = 3L)

# NMI after NiftyReg affine registrations
score_after <- similarity(nmi_registered$image, OSTA_RCTD_img, 
                          interpolation = 3L,
                          threads = getOption("RNiftyReg.threads"))
print(paste("NMI after:", score_after)) 

# Overlayed image of both sets of plots
writePNG(nmi_registered$image, "data/aligned_seurat_to_osta.png")



# References ####
sink(paste0("references.txt"))
print("HL Crowell*°, Y Dong*, I Billato, P Cai, M Emons, S Gunz, B Guo, M Li, A Mahmoud, A Manukyan, H Pagès, P Panwar, S Rao, CJ Sargeant, L Shepherd Kern, M Ramos, J Sun, M Totty, VJ Carey, Y Chen, L Collado-Torres, S Ghazanfar, KD Hansen, K Martinowich, KR Maynard, E Patrick, D Righelli, D Risso, S Tiberi, L Waldron, R Gottardo†°, MD Robinson†°, SC Hicks†°, LM Weber†°. Orchestrating spatial transcriptomics analysis with Bioconductor. bioRxiv (2025). DOI: 10.1101/2025.11.20.688607")
print("Janesick, A., Shelansky, R., Gottscho, A.D. et al. High resolution mapping of the tumor microenvironment using integrated single-cell, spatial and in situ analysis. Nat Commun 14, 8353 (2023). https://doi.org/10.1038/s41467-023-43458-x")
citation("OSTA.data")
print("")
citation("tidyverse") 
print("")
citation("Seurat")
print("")
citation("patchwork")
print("")
citation("visiumStitched")
print("")
citation("CARDspa")
print("")
citation("pheatmap")
print("")
citation("SpatialExperiment")
print("")
citation("VisiumIO")
print("")
citation("hdf5r")
print("")
citation("devtools")
print("")
citation("RNiftyReg")
print("")
citation("png")
print("")
citation("pals")
print("")
citation("spacexr")
sink()

# Session Info ####

sink(paste0("session_info_", Sys.Date(), ".txt"))
sessionInfo()
sink()

# Use this link to compare your predictions with the predictions the
# ebook authors got with RCTD and CARD algorithms
# (their images are flipped vertically)
# https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#visualization