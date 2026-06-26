suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(harmony)
  library(ggplot2)
  library(ggsci)
  library(hdf5r)
  library(ggsignif)
  library(RColorBrewer)
  library(ggrepel)
})

## 1. Environment & Paths
data_dir <- "./data"
out_dir <- "./results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

aging_samples <- c("Aging_1", "Aging_2", "Aging_3", "Aging_4", "Aging_5", "Aging_6", "Aging_7")
young_samples <- c("Young_1", "Young_2", "Young_3", "Young_4", "Young_5")

plot_theme <- theme(
  plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
  panel.background = element_blank(),
  panel.border = element_rect(color = "black", fill = NA)
)

## 2. Data Loading
h5_files <- list.files(data_dir, pattern = "\\.h5$", full.names = TRUE)
sample_names <- tools::file_path_sans_ext(basename(h5_files))
keep_idx <- sample_names %in% c(aging_samples, young_samples)
h5_files <- h5_files[keep_idx]
sample_names <- sample_names[keep_idx]

seurat_list <- lapply(seq_along(h5_files), function(i) {
  obj <- CreateSeuratObject(counts = Read10X_h5(h5_files[i]), project = sample_names[i], min.cells = 3, min.features = 200)
  obj$original_sample <- sample_names[i]
  obj$group <- ifelse(sample_names[i] %in% aging_samples, "Aging", "Young")
  return(obj)
})
names(seurat_list) <- sample_names

merged_seurat <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = sample_names)
merged_seurat$original_sample <- factor(merged_seurat$original_sample, levels = c(aging_samples, young_samples))
rm(seurat_list)

## 3. QC & Filtering
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")
HB_genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
HB_matched <- HB_genes[HB_genes %in% rownames(merged_seurat)]
merged_seurat[["percent.HB"]] <- PercentageFeatureSet(merged_seurat, features = HB_matched)

qc1 <- VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB"), group.by = "original_sample", pt.size = 0, ncol = 4) + geom_boxplot(width = 0.1, outlier.shape = NA)
ggsave(file.path(out_dir, "QC_1_before_filtering.pdf"), qc1, width = 14, height = 7)

filtered_seurat <- subset(merged_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10 & percent.HB < 5)
qc2 <- VlnPlot(filtered_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB"), group.by = "original_sample", pt.size = 0, ncol = 4) + geom_boxplot(width = 0.1, outlier.shape = NA)
ggsave(file.path(out_dir, "QC_2_after_filtering.pdf"), qc2, width = 14, height = 7)
rm(merged_seurat)

## 4. Normalization & Cell Cycle
filtered_seurat <- NormalizeData(filtered_seurat) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% ScaleData(features = rownames(filtered_seurat))
data(cc.genes.updated.2019)
filtered_seurat <- CellCycleScoring(filtered_seurat, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
filtered_seurat <- ScaleData(filtered_seurat, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(filtered_seurat))

## 5. Dimensionality Reduction & Clustering
filtered_seurat <- RunPCA(filtered_seurat, features = VariableFeatures(filtered_seurat), npcs = 50, verbose = FALSE)
filtered_seurat <- RunHarmony(object = filtered_seurat, group.by.vars = "original_sample", reduction.use = "pca", dims.use = 1:30, plot_convergence = FALSE)

filtered_seurat <- FindNeighbors(filtered_seurat, reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.03) %>%
  RunUMAP(reduction = "harmony", dims = 1:20)

## 6. Cluster Markers & Annotation
cluster_markers <- FindAllMarkers(filtered_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
ggsave(file.path(out_dir, "Cluster_Marker_Heatmap.pdf"), DoHeatmap(subset(filtered_seurat, downsample = 200), features = top10$gene), width = 12, height = 10)

cluster_ids <- c("0"="Stromal", "1"="Immune", "2"="Epithelial", "3"="Perivascular", "4"="Endothelial", "5"="Immune", "6"="Epithelial")
filtered_seurat$cell_type <- factor(cluster_ids[as.character(filtered_seurat$seurat_clusters)], levels = c("Stromal", "Immune", "Epithelial", "Perivascular", "Endothelial"))
Idents(filtered_seurat) <- "cell_type"

## 7. Differential Expression (Volcano)
all_markers <- FindAllMarkers(filtered_seurat, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all_markers, file.path(out_dir, "All_5_CellTypes_Markers_Bidirectional.csv"), row.names = FALSE)

sig_data <- all_markers %>% filter(p_val_adj < 0.05) %>% mutate(label = ifelse(avg_log2FC > 0, "sigUp", "sigDown"))
set.seed(42)
sig_data$x_jitter <- as.numeric(factor(sig_data$cluster, levels = levels(filtered_seurat$cell_type))) + runif(nrow(sig_data), -0.35, 0.35)
top_genes <- sig_data %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5) %>% bind_rows(sig_data %>% group_by(cluster) %>% slice_min(avg_log2FC, n = 5)) %>% mutate(gene = toupper(gene))

p_volcano <- ggplot() +
  geom_point(data = sig_data, aes(x = x_jitter, y = avg_log2FC, color = label), size = 0.8) +
  scale_color_manual(values = c("sigDown" = "#0077c0", "sigUp" = "#c72d2e")) +
  geom_text_repel(data = top_genes, aes(x = x_jitter, y = avg_log2FC, label = gene), size = 3.5, fontface = "italic", segment.color = NA) +
  labs(x = "Cell type", y = expression("Average log"[2]*" fold change"), title = "Differentially expressed genes") +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "Volcano_Plot_5_CellTypes_Refined.pdf"), p_volcano, width = 11, height = 6)

## 8. UMAP Visualization
p_clust <- DimPlot(filtered_seurat, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("UMAP by clusters") + plot_theme
p_type <- DimPlot(filtered_seurat, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) + ggtitle("UMAP by cell types") + plot_theme
p_grp <- DimPlot(filtered_seurat, reduction = "umap", group.by = "group", cols = c("#FF7F00", "#1F78B4")) + ggtitle("UMAP by groups") + plot_theme

ggsave(file.path(out_dir, "UMAP_1_Clusters.pdf"), p_clust, width = 7, height = 6)
ggsave(file.path(out_dir, "UMAP_2_CellTypes.pdf"), p_type, width = 7, height = 6)
ggsave(file.path(out_dir, "UMAP_3_Groups.pdf"), p_grp, width = 7, height = 6)
ggsave(file.path(out_dir, "UMAP_4_Combined.pdf"), p_type + p_grp, width = 14, height = 6)

p_split <- DimPlot(filtered_seurat, reduction = "umap", group.by = "cell_type", split.by = "group", label = TRUE, repel = TRUE) + plot_theme + ggtitle("UMAP by cell types (split by group)")
ggsave(file.path(out_dir, "UMAP_5_Split_By_Group.pdf"), p_split, width = 14, height = 6)

## 9. DotPlot
marker_list <- list("Stromal"=c("IGF1", "DCN", "PCOLCE", "LUM", "SFRP4"), "Immune"=c("PTPRC", "NKG7", "CD2", "TYROBP", "CCL5"), "Epithelial"=c("WFDC2", "KRT18", "KRT8", "CAPS", "EPCAM"), "Perivascular"=c("NOTCH3", "RGS5", "MYH11", "OLFML2B", "GUCY1A2"), "Endothelial"=c("EGFL7", "PECAM1", "VWF", "CD34", "CLDN5"))
p_dot <- DotPlot(filtered_seurat, features = marker_list, group.by = "cell_type", cols = c("lightgrey", "red"), dot.scale = 5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"), plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "Marker gene expression by cell type", x = "Features", y = "Identity")
ggsave(file.path(out_dir, "DotPlot_All_Markers.pdf"), p_dot, width = 12, height = 6)