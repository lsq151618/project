## 0. Initialization
suppressPackageStartupMessages({
  library(pheatmap)
  library(openxlsx)
  library(grid)
  library(gtable)
})

data_dir <- "./"
diff_file_path <- file.path(data_dir, "DEGs_significant.csv")

## 1. Data Loading & Preprocessing
diff_data <- read.csv(diff_file_path, stringsAsFactors = FALSE, check.names = FALSE)
diff_data <- as.data.frame(diff_data)

colnames(diff_data)[colnames(diff_data) == "gene name"] <- "gene name"

required_cols <- c("gene name", "logFC", "pvalue", "padj", "BZBS1", "BZBS2", "BZBS3", "D-gal1", "D-gal2", "D-gal3")
if (!all(required_cols %in% colnames(diff_data))) {
  stop("Input CSV missing required columns or expression data.")
}

diff_data$logFC <- as.numeric(as.character(diff_data$logFC))
diff_data$padj <- as.numeric(as.character(diff_data$padj))

if (any(duplicated(diff_data$`gene name`))) {
  diff_data <- diff_data[!duplicated(diff_data$`gene name`), ]
}
rownames(diff_data) <- diff_data$`gene name`

## 2. Filtering & AMPK Selection
diff_data$classification <- "None"
diff_data$classification[diff_data$padj < 0.05 & diff_data$logFC > 0.5] <- "Up"
diff_data$classification[diff_data$padj < 0.05 & diff_data$logFC < -0.5] <- "Down"

sig_genes <- diff_data[diff_data$classification != "None", ]

ampk_targets <- c("CREB3L3", "CREB3L4", "FOXO1", "INSR", "IRS2", "LIPE", "MLYCD", 
                  "PCK2", "PFKFB2", "ACACA", "AKT3", "CCND1", "CD36", "EEF2K", 
                  "FASN", "HMGCR", "LEPR", "PFKFB3", "PFKL", "PIK3R3", "PPARG", 
                  "PPP2R3A", "SCD", "SREBF1", "TSC2")

found_genes <- rownames(sig_genes)[toupper(rownames(sig_genes)) %in% toupper(ampk_targets)]
ampk_plot_data <- sig_genes[found_genes, ]

## 3. Visualization
expr_cols <- c("BZBS1", "BZBS2", "BZBS3", "D-gal1", "D-gal2", "D-gal3")

if (nrow(ampk_plot_data) > 1) {
  mat <- as.matrix(ampk_plot_data[, expr_cols])
  mat <- apply(mat, 2, as.numeric)
  rownames(mat) <- rownames(ampk_plot_data)
  
  mat_transposed <- t(scale(t(mat)))
  
  gene_status_info <- data.frame(Status = factor(ampk_plot_data$classification, levels = c("Up", "Down")))
  rownames(gene_status_info) <- rownames(ampk_plot_data)
  
  ann_colors <- list(Status = c(Up = "#D6604D", Down = "#4393C3"))
  
  p <- pheatmap(t(mat_transposed),
                annotation_col = gene_status_info,
                annotation_colors = ann_colors,
                color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                cluster_cols = TRUE,
                cluster_rows = FALSE,
                show_colnames = TRUE,
                show_rownames = TRUE,
                fontsize_col = 10,
                fontsize_row = 11,
                fontface_col = "bold",
                fontface_row = "bold",
                angle_col = 45,
                main = "AMPK pathway differential genes (Z-score)",
                border_color = "white",
                silent = TRUE)
  
  ## Legend adjustment
  g <- p$gtable
  ann_leg_idx <- which(g$layout$name == "annotation_legend")
  if(length(ann_leg_idx) > 0) {
    ann_legend_grob <- g$grobs[[ann_leg_idx]]
    ann_legend_grob$vp <- viewport(x = unit(0.95, "npc"), y = unit(0.3, "npc"), just = c("right", "top"))
    g$grobs[[ann_leg_idx]] <- ann_legend_grob
  }
  
  pdf_filename <- file.path(data_dir, "AMPK_Pathway_Heatmap_Final_Fixed.pdf")
  pdf(pdf_filename, width = 12, height = 8)
  grid.newpage()
  grid.draw(g)
  dev.off()
  
  cat(">>> Heatmap saved to:", normalizePath(pdf_filename, mustWork = FALSE), "\n")
} else {
  warning("Insufficient genes matched for heatmap.")
}