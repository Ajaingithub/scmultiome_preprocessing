### Neeed to have a object that comprises of both RNA and ATAC
# The commented part should be already being done
# old <- Read10X("/research/labs/immunology/goronzy_weyand/GoronzyLab_SCG/Huimin/Single_cell_Multiome/analysis/Old/outs/filtered_feature_bc_matrix/")
# young <- Read10X("/research/labs/immunology/goronzy_weyand/GoronzyLab_SCG/Huimin/Single_cell_Multiome/analysis/Young/outs/filtered_feature_bc_matrix/")
#
# GEX_old = old$`Gene Expression`
# GEX_young = young$`Gene Expression`
#
# colnames(GEX_old) <- paste("Old#",colnames(GEX_old),sep = "")
# colnames(GEX_young) <- paste("Young#",colnames(GEX_young),sep = "")
# GEX=cbind(GEX_old,GEX_young)
#
# ATAC_old <- old$Peaks
# ATAC_young <- young$Peaks
#
# ATAC_old.peaks <- StringToGRanges(rownames(ATAC_old), sep = c(":", "-"))
# ATAC_young.peaks <- StringToGRanges(rownames(ATAC_young), sep = c(":", "-"))
# all.peaks <- c(ATAC_old.peaks, ATAC_young.peaks)
# #####################
# reduced.peaks <- reduce(all.peaks, ignore.strand=TRUE)
#
# young_fragment <- CreateFragmentObject("/research/labs/immunology/goronzy_weyand/GoronzyLab_SCG/Huimin/Single_cell_Multiome/analysis/Young/outs/atac_fragments.tsv.gz")
# Young <- FeatureMatrix(young_fragment,features = reduced.peaks, cells = colnames(ATAC_young))
#
# Old_fragment <- CreateFragmentObject("/research/labs/immunology/goronzy_weyand/GoronzyLab_SCG/Huimin/Single_cell_Multiome/analysis/Old/outs/atac_fragments.tsv.gz")
# Old <- FeatureMatrix(Old_fragment,features = reduced.peaks, cells = colnames(ATAC_old))
#
# colnames(Young) <- paste("Young#",colnames(Young),sep = "")
# colnames(Old) <- paste("Old#",colnames(Old),sep = "")
#
# ### Checking the GEX column are the same as the column of the ATAC
# ATAC=cbind(Old,Young)
#
# ## Creating the Seurat Object
# resting <- CreateSeuratObject(counts = GEX, project = "CD4_resting")
# resting[["percent.mt"]] <- PercentageFeatureSet(resting, pattern = "^MT-")
#
# ### Now we will add the ATAC seq assay to the object
# ## We will use only the standard chromosomes
# grange.counts <- StringToGRanges(rownames(ATAC), sep = c(":", "-"))
# grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
# atac_counts <- ATAC[as.vector(grange.use),]
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) # EnsDb.Hsapiens.v75 is for hg19 while EnsDb.Hsapiens.v86 for hg38
# seqlevelsStyle(annotations) <- 'UCSC'
# genome(annotations) <- "hg38"
#
# write.table(colnames(atac_counts), file="/research/labs/immunology/goronzy_weyand/GoronzyLab_SCG/Huimin/Single_cell_Multiome/analysis/cellname.txt",
#             sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
#
# ## In bash
# ## [m256617@mforgehn1 outs]$ zcat atac_fragments.tsv.gz > /research/labs/immunology/goronzy_weyand/GoronzyLab_SCG/Huimin/Single_cell_Multiome/analysis/young_atac_fragment.ts
# ### awk '{print $1"\t"$2"\t"$3"\tOld#"$4"\t"$5}' old_atac_fragment.tsv > old_atac_fragment_2.tsv
# ## grep -wFf cellname.txt old_and_young_atac_fragment_2.tsv > old_and_young_atac_fragment_3.tsv
# # (seurat_new) bgzip old_and_young_atac_fragment_cells_req.tsv
# # (seurat_new) tabix -p bed old_and_young_atac_fragment_cells_req.tsv.gz
#
#
# frag.file <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_SCG/Huimin/Single_cell_Multiome/analysis/old_and_young_atac_fragment_cells_req_sorted.tsv.gz"
# chrom_assay <- CreateChromatinAssay(
#   counts = atac_counts,
#   sep = c(":", "-"),
#   genome = 'hg38',
#   fragments = frag.file,
#   min.cells = 10,
#   annotation = annotations
# )
# resting[["ATAC"]] <- chrom_assay

scmultiome_pre <- function(obj, savedir, objname, subset_cellnames){
  library(BiocParallel)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(Seurat)
  library(SeuratDisk)
  library(Signac)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(EnsDb.Hsapiens.v86) # for hg38
  # library(EnsDb.Hsapiens.v75) # for hg19
  library(BSgenome.Hsapiens.UCSC.hg38)
  # library(BSgenome.Hsapiens.UCSC.hg19)
  library(chromVAR)
  library(JASPAR2020)
  library(TFBSTools)
  library(motifmatchr)
  library(chromVAR)
  library(JASPAR2020)
  library(TFBSTools)
  library(motifmatchr)
  library(BSgenome.Hsapiens.UCSC.hg38)

  message(paste("Subsetting the obj for",objname))
  obj_2 <- subset(obj, cells=subset_cellnames)

  message(paste("Filtering the object"))
  ### Filtering
  resting_2 <- subset(x = obj_2, subset = nCount_ATAC < 7e4 & nCount_ATAC > 5e3 & nCount_RNA < 25000 & nCount_RNA > 1000 & percent.mt < 20)


  # We next perform pre-processing and dimensional reduction on both assays independently, using standard approaches for RNA and ATAC-seq data.
  message("Performing the RNA analysis")
  DefaultAssay(resting_2) <- "RNA"
  resting_2 <- SCTransform(resting_2, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

  message("Performing the ATAC analysis")

  # ATAC analysis
  # We exclude the first dimension as this is typically correlated with sequencing depth
  # We then reduced the dimensionality of these large binary matrices using a term frequency-inverse document frequency (‘‘TF-IDF’’)
  # transformation. To do this, we first weighted all the sites for individual cells by the total number of sites accessible in that
  # cell (‘‘term frequency’’). We then multiplied these weighted values by log(1 + the inverse frequency of each site across all cells),
  # the ‘‘inverse document frequency.’’ We then used singular value decomposition on the TF-IDF matrix to generate a lower dimensional representation of the data by only retaining the 2nd through 10th dimensions
  DefaultAssay(resting_2) <- "ATAC"

  resting_2 <- RunTFIDF(resting_2)
  resting_2 <- FindTopFeatures(resting_2, min.cutoff = 'q0')
  resting_2 <- RunSVD(resting_2)
  resting_2 <- RunUMAP(resting_2, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

  message("Performing Multimodal")
  resting_2 <- FindMultiModalNeighbors(resting_2, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
  resting_2 <- RunUMAP(resting_2, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  resting_2 <- FindClusters(resting_2, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.8)

  p1 <- DimPlot(resting_2, reduction = "umap.rna", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("RNA")
  p2 <- DimPlot(resting_2, reduction = "umap.atac", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("ATAC")
  p3 <- DimPlot(resting_2, reduction = "wnn.umap", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("WNN")


  pdf(paste(savedir,"RNA_and_ATAC.pdf",sep = ""), width = 14, height = 6)
  print(p1 + p2 + p3 & theme(plot.title = element_text(hjust = 0.5)))
  dev.off()

  p1 <- DimPlot(resting_2, reduction = "umap.rna", group.by = "Age", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("RNA")
  p2 <- DimPlot(resting_2, reduction = "umap.atac",group.by = "Age", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("ATAC")
  p3 <- DimPlot(resting_2, reduction = "wnn.umap", group.by = "Age", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("WNN")

  pdf(paste(savedir,"Age.pdf",sep = ""), width = 14, height = 6)
  print(p1 + p2 + p3 & theme(plot.title = element_text(hjust = 0.5)))
  dev.off()

  p1 <- FeaturePlot(resting_2, reduction = "umap.rna", features = c("rna_CD69")) + ggtitle("RNA")
  p3 <- FeaturePlot(resting_2, reduction = "wnn.umap", features = c("rna_CD69")) + ggtitle("WNN")
  pdf(paste(savedir,"Activation.pdf",sep = ""), width = 9, height = 6)
  print(p1 + p3 & theme(plot.title = element_text(hjust = 0.5)))
  dev.off()

  p <- CoveragePlot(resting_2, region = 'CD69', features = 'CD69', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE)

  pdf(paste(savedir,"CD69_ATAC_RNA.pdf",sep = ""), width = 9, height = 6)
  print(p)
  dev.off()

  message(paste("Please find the UMAP for RNA and ATAC Coverage at this location ",savedir))

  saveRDS(resting_2,file=paste(savedir,objname,".RDS",sep = ""))

  message("performing ChromVar")
  ### Find the Enriched Motif. we will use chromVar package from Greenleaf
  # This calculates a per-cell accessibility score for known motifs, and adds these scores as a third assay (chromvar) in the Seurat object.
  # chromVAR for the analysis of motif accessibility in scATAC-seq
  # presto for fast differential expression analyses.
  # TFBSTools for TFBS analysis
  # JASPAR2020 for JASPAR motif models
  # motifmatchr for motif matching
  # BSgenome.Hsapiens.UCSC.hg38 for chromVAR

  # Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
  DefaultAssay(resting_2) <- "ATAC"
  pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
  motif.matrix <- CreateMotifMatrix(features = granges(resting_2), pwm = pwm_set, genome = 'hg38', use.counts = FALSE)
  motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
  resting_2 <- SetAssayData(resting_2, assay = 'ATAC', slot = 'motifs', new.data = motif.object)

  # Note that this step can take 30-60 minutes
  resting_2 <- RunChromVAR(
    object = resting_2,
    genome = BSgenome.Hsapiens.UCSC.hg38
  )

  # We’d like to quantify this relationship, and search across all cell types to find similar examples. To do so, we will use the
  # presto package to perform fast differential expression. We run two tests: one using gene expression data, and the other using
  # chromVAR motif accessibilities. presto calculates a p-value based on the Wilcox rank sum test, which is also the default test in
  # Seurat, and we restrict our search to TFs that return significant results in both tests.
  # presto also calculates an “AUC” statistic, which reflects the power of each gene (or motif) to serve as a marker of cell type.

  # resting <- resting_2
  # markers_rna <- presto:::wilcoxauc.Seurat(X = resting, group_by = 'seurat_clusters', assay = 'data', seurat_assay = 'SCT')
  # markers_motifs <- presto:::wilcoxauc.Seurat(X = resting, group_by = 'seurat_clusters', assay = 'data', seurat_assay = 'chromvar')
  # motif.names <- markers_motifs$feature
  # colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))
  # colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))
  # markers_rna$gene <- markers_rna$RNA.feature
  # markers_motifs$gene <- ConvertMotifID(resting_2, id = motif.names)
  #
  # # a simple function to implement the procedure above
  # topTFs <- function(celltype, padj.cutoff = 1e-2) {
  #   ctmarkers_rna <- dplyr::filter(
  #     markers_rna, RNA.group == celltype, RNA.padj < padj.cutoff, RNA.logFC > 0) %>%
  #     arrange(-RNA.auc)
  #   ctmarkers_motif <- dplyr::filter(
  #     markers_motifs, motif.group == celltype, motif.padj < padj.cutoff, motif.logFC > 0) %>%
  #     arrange(-motif.auc)
  #   top_tfs <- inner_join(
  #     x = ctmarkers_rna[, c(2, 11, 6, 7)],
  #     y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
  #   )
  #   top_tfs$avg_auc <- (top_tfs$RNA.auc + top_tfs$motif.auc) / 2
  #   top_tfs <- arrange(top_tfs, -avg_auc)
  #   return(top_tfs)
  # }
  #
  # head(topTFs("1"), 3)
  #

  saveRDS(resting_2,file=paste(savedir,objname,".RDS",sep = ""))

  message("Please find the object and UMAP at this location")
}





