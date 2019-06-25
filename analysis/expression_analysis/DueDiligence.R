#!/usr/bin/env Rscript
#
# Since Alex's Treg data did not have any published findings, I would like
# to do a due-diligence analysis of the expression data, to better understand
# confounding factors affecting the variant calling results.
#
# Usage:
#   run R script from the home directory for appropriate file location definitions
#
# NOTE:
#   Using Seurat v2, some function calls and features may not be portable to older
#     versions of Seurat.
#
# Integrated analysis performed by following Seurat vignette: https://satijalab.org/seurat/v3.0/immune_alignment.html
#

library(dplyr)
library(RPostgreSQL)
library(Seurat)
library(Matrix)
library(stringr)
library(hashmap)
library(ggplot2)
library(cowplot)

file_names <- list.files('./raw_data/matrices/');
sample_details <- str_match(file_names, "^(\\d+)([BF])");

# load donor information from database

drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, dbname = "ms", host = "localhost", port = 5432)
donor_df <- dbGetQuery(con, 'SELECT * FROM donor');

donorNum2scRNAseqDate <- hashmap(donor_df$Sample, donor_df$scRNAseqDate)
donorNum2Group <- hashmap(donor_df$Sample, donor_df$Group)

filenamePrefix2scRNAseqDate <- hashmap(
  sample_details[,1],
  donorNum2scRNAseqDate[[as.integer(sample_details[,2])]]
)

filenamePrefix2Group <- hashmap(
  sample_details[,1],
  donorNum2Group[[as.integer(sample_details[,2])]]
)

# load matrix data
srt_objs <- c();
for (i in seq(length(file_names))) {
  print(paste("PROCESSING: ", file_names[i], sep=''));
  nd <- read.table(paste('./raw_data/matrices/', file_names[i], sep=''), 
          sep = "\t", row.names = 1, header = TRUE);
  srt_obj <- UpdateSeuratObject(CreateSeuratObject(nd));
  srt_obj <- NormalizeData(srt_obj, verbose=FALSE);
  srt_obj <- FindVariableFeatures(srt_obj, selection.method = "vst", nfeatures = 2000);
  
  srt_objs <- c(srt_objs, srt_obj);
}

for (i in seq(length(file_names))) {
  srt_objs[[i]] <- AddMetaData(object = srt_objs[[i]], col.name='group',
              metadata = filenamePrefix2Group[[sample_details[,1]]][[i]]);
  srt_objs[[i]] <- AddMetaData(object = srt_objs[[i]], col.name='batch',
              metadata = as.character(filenamePrefix2scRNAseqDate[[sample_details[,1]]][[i]]))
  srt_objs[[i]] <- AddMetaData(object = srt_objs[[i]], col.name='source.tissue',
              metadata = sample_details[,3][[i]])
  srt_objs[[i]] <- AddMetaData(object = srt_objs[[i]], col.name='sample.id',
              metadata = sample_details[,1][[i]])
}

# 4Fat (file 22) has very low yield, only 16 cells, will throw out to allow better data integration
# 
# sapply(srt_objs, function(x){ ncol(x) } )

treg.anchors <- FindIntegrationAnchors(object.list = srt_objs[-c(22)], dims=1:20, k.filter = 50)
treg.combined <- IntegrateData(anchorset = treg.anchors, dims = 1:20)
treg.combined.raw <- treg.combined;

# Seurat v3 analysis
DefaultAssay(treg.combined) <- "RNA"

# Basic QC
treg.combined[["percent.mt"]] <- PercentageFeatureSet(treg.combined, pattern = "^ENSG\\d+-MT-")
treg.combined[["percent.ercc"]] <- PercentageFeatureSet(treg.combined, pattern = "^ERCC-")
treg.combined[["cell.type.treg"]] <- grepl("Treg", colnames(treg.combined))

VlnPlot(treg.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ercc"), ncol = 4)

treg.combined <- subset(treg.combined,
                        nFeature_RNA > 200 & nFeature_RNA < 2500
                        & percent.mt < 5
                        & percent.ercc < 1
                        & cell.type.treg)

# DefaultAssay(treg.combined) <- "integrated"
# NOTE: not sure this is correct 
# using the 'integrated' assay causes analysis to fail
DefaultAssay(treg.combined) <- "RNA"

treg.combined <- FindVariableFeatures(object = treg.combined)

# Run the standard workflow for visualization and clustering
treg.combined <- ScaleData(treg.combined, verbose = FALSE, 
                           vars.to.regress = c('batch', 
                                               'nCount_RNA', 
                                               'percent.mt', 
                                               'percent.ercc'))
treg.combined <- RunPCA(treg.combined, npcs = 30, verbose = FALSE)

# UMAP and Clustering
treg.combined <- RunUMAP(treg.combined, reduction = "pca", dims = 1:20)
treg.combined <- FindNeighbors(treg.combined, reduction = "pca", dims = 1:20)
treg.combined <- FindClusters(treg.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(treg.combined, reduction = "umap", group.by = "group")
p2 <- DimPlot(treg.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

p1 <- DimPlot(treg.combined, reduction = "umap", group.by = "batch")
p2 <- DimPlot(treg.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

p1 <- DimPlot(treg.combined, reduction = "umap", group.by = "source.tissue")
p2 <- DimPlot(treg.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

p1 <- DimPlot(treg.combined, reduction = "umap", group.by = "sample.id")
p2 <- DimPlot(treg.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

############################################################
## Separate Fat
############################################################

fat.treg.combined <- subset(treg.combined, source.tissue == 'F');
fat.treg.combined <- RunPCA(fat.treg.combined, npcs = 30, verbose = FALSE);

# UMAP and Clustering
fat.treg.combined <- RunUMAP(fat.treg.combined, reduction = "pca", dims = 1:20);
fat.treg.combined <- FindNeighbors(fat.treg.combined, reduction = "pca", dims = 1:20);
fat.treg.combined <- FindClusters(fat.treg.combined, resolution = 0.5);

# Visualization
p1 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "group");
p2 <- DimPlot(fat.treg.combined, reduction = "umap", label = TRUE);
plot_grid(p1, p2);

markers.0 <- FindMarkers(fat.treg.combined, ident.1 = 0, ident.2 = 1);

p1 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "batch")
p2 <- DimPlot(fat.treg.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

p1 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "source.tissue")
p2 <- DimPlot(fat.treg.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

p1 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "sample.id")
p2 <- DimPlot(fat.treg.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

# Genes

rownames(treg.combined)[grepl("CCL21", rownames(treg.combined))]

p1 <- FeaturePlot(fat.treg.combined, c("ENSG00000111537-IFNG"))
p2 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(fat.treg.combined, c("ENSG00000170345-FOS"))
p2 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(fat.treg.combined, c("ENSG00000169429-IL8"))
p2 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(fat.treg.combined, c("ENSG00000112486-CCR6"))
p2 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(fat.treg.combined, c("ENSG00000161570-CCL5"))
p2 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(fat.treg.combined, c("ENSG00000168329-CX3CR1"))
p2 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(fat.treg.combined, c("ENSG00000126353-CCR7"))
p2 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(fat.treg.combined, c("ENSG00000153563-CD8A"))
p2 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(fat.treg.combined, c("ENSG00000010610-CD4"))
p2 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(fat.treg.combined, c("ENSG00000167286-CD3D"))
p2 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(fat.treg.combined, c("ENSG00000148773-MKI67"))
p2 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(fat.treg.combined, c("ENSG00000132646-PCNA"))
p2 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(fat.treg.combined, c("ENSG00000128016-ZFP36"))
p2 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(fat.treg.combined, c("ENSG00000196154-S100A4"))
p2 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

# CD45
p1 <- FeaturePlot(fat.treg.combined, c("ENSG00000179820-MYADM"))
p2 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(fat.treg.combined, c("ENSG00000153234-NR4A2"))
p2 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(fat.treg.combined, c("ENSG00000271204-RP11-138A9.1"))
p2 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(fat.treg.combined, c("ENSG00000049768-FOXP3"))
p2 <- FeaturePlot(fat.treg.combined, c("ENSG00000134460-IL2RA"))
p3 <- FeaturePlot(fat.treg.combined, c("ENSG00000168685-IL7R"))
p4 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2, p3, p4, ncol = 4)

p1 <- FeaturePlot(fat.treg.combined, c("ENSG00000049768-FOXP3"))
p2 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(fat.treg.combined, c("ENSG00000134460-IL2RA"))
p2 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(fat.treg.combined, c("ENSG00000168685-IL7R"))
p2 <- DimPlot(fat.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

fat.treg.combined[['cell.names']] <- Idents(fat.treg.combined);
Idents(fat.treg.combined) <- "group"
fat.diff.group.markers <- FindMarkers(fat.treg.combined, min.diff.pct = 0.20,
                                      ident.1 = 'HC', ident.2 = 'MS', verbose = F)
head(fat.diff.group.markers, n=15)

############################################################
## Separate Blood
############################################################

blood.treg.combined <- subset(treg.combined, source.tissue == 'B');
blood.treg.combined <- RunPCA(blood.treg.combined, npcs = 30, verbose = FALSE);

blood.treg.combined <- RunUMAP(blood.treg.combined, reduction = "pca", dims = 1:20);
blood.treg.combined <- FindNeighbors(blood.treg.combined, reduction = "pca", dims = 1:20);
blood.treg.combined <- FindClusters(blood.treg.combined, resolution = 0.5);

# Visualization
p1 <- DimPlot(blood.treg.combined, reduction = "umap", group.by = "group");
p2 <- DimPlot(blood.treg.combined, reduction = "umap", label = TRUE);
plot_grid(p1, p2);

p1 <- DimPlot(blood.treg.combined, reduction = "umap", group.by = "batch")
p2 <- DimPlot(blood.treg.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

p1 <- DimPlot(blood.treg.combined, reduction = "umap", group.by = "source.tissue")
p2 <- DimPlot(blood.treg.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

p1 <- DimPlot(blood.treg.combined, reduction = "umap", group.by = "sample.id")
p2 <- DimPlot(blood.treg.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

p1 <- FeaturePlot(blood.treg.combined, c("ENSG00000049768-FOXP3"))
p2 <- FeaturePlot(blood.treg.combined, c("ENSG00000134460-IL2RA"))
p3 <- FeaturePlot(blood.treg.combined, c("ENSG00000168685-IL7R"))
p4 <- DimPlot(blood.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2, p3, p4, ncol = 4)

p1 <- FeaturePlot(blood.treg.combined, c("ENSG00000148773-MKI67"))
p2 <- DimPlot(blood.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(blood.treg.combined, c("ENSG00000128016-ZFP36"))
p2 <- DimPlot(blood.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(blood.treg.combined, c("ENSG00000196154-S100A4"))
p2 <- DimPlot(blood.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

# CD45
p1 <- FeaturePlot(blood.treg.combined, c("ENSG00000179820-MYADM"))
p2 <- DimPlot(blood.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(blood.treg.combined, c("ENSG00000153234-NR4A2"))
p2 <- DimPlot(blood.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(blood.treg.combined, c("ENSG00000271204-RP11-138A9.1"))
p2 <- DimPlot(blood.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(blood.treg.combined, c("ENSG00000049768-FOXP3"))
p2 <- DimPlot(blood.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(blood.treg.combined, c("ENSG00000134460-IL2RA"))
p2 <- DimPlot(blood.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

p1 <- FeaturePlot(blood.treg.combined, c("ENSG00000168685-IL7R"))
p2 <- DimPlot(blood.treg.combined, reduction = "umap", group.by = "group");
plot_grid(p1, p2);

blood.treg.combined[['cell.names']] <- Idents(blood.treg.combined);
Idents(blood.treg.combined) <- "group";
blood.diff.group.markers <- FindMarkers(blood.treg.combined,
                                      ident.1 = 'HC', ident.2 = 'MS', verbose = F)
head(blood.diff.group.markers, n=15)

#########################################################
## EXAMINE CHANGES TO GENE EXPRESSION BETWEEN GROUPS
#########################################################

blood.diff.group.markers$Feature <- rownames(blood.diff.group.markers)
blood.diff.group.markers$NegLogP <- -log10(blood.diff.group.markers$p_val)
  
fat.diff.group.markers$Feature <- rownames(fat.diff.group.markers)
fat.diff.group.markers$NegLogP <- -log10(fat.diff.group.markers$p_val)

# Write summary files
write.csv(blood.diff.group.markers, file = "./Blood_Diff_Group_Markers.csv", row.names = F)
write.csv(fat.diff.group.markers, file = "./Fat_Diff_Group_Markers.csv", row.names = F)

# ad hoc code
p_threshold <- 1e-4;
bd <- subset(blood.diff.group.markers, p_val_adj < p_threshold)
ft <- subset(fat.diff.group.markers, p_val_adj < p_threshold)

intersect(rownames(bd), rownames(ft))
setdiff(rownames(ft), rownames(bd))
setdiff(rownames(bd), rownames(ft))

joined.diff.markers <- fat.diff.group.markers %>% 
                          merge(blood.diff.group.markers, 
                                    by='row.names',
                                    suffixes=c('.fat', '.blood'))

# volcano plots
plot(blood.diff.group.markers$avg_logFC, -log10(blood.diff.group.markers$p_val))
plot(fat.diff.group.markers$avg_logFC, -log10(fat.diff.group.markers$p_val))


