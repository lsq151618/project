## 0. Initialization
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(readxl)
  library(pheatmap)
  library(openxlsx)
})

## 1. Data Loading & Preprocessing
data_dir <- "./"
diff_file_path <- file.path(data_dir, "differential_expression.xlsx")
diff_data <- as.data.frame(read_excel(diff_file_path))

required_cols <- c("gene name", "logFC", "pvalue", "padj")
if (!all(required_cols %in% colnames(diff_data))) {
  stop("Excel file must contain 'gene name', 'logFC', 'pvalue', and 'padj' columns.")
}

if (!is.numeric(diff_data$logFC)) {
  diff_data$logFC <- as.numeric(as.character(diff_data$logFC))
}

if (!is.numeric(diff_data$padj)) {
  diff_data$padj <- as.numeric(as.character(diff_data$padj))
}

if (any(duplicated(diff_data$`gene name`))) {
  diff_data <- diff_data[!duplicated(diff_data$`gene name`), ]
}

rownames(diff_data) <- diff_data$`gene name`

## 2. Statistical Thresholding
diff_data$classification <- ifelse(
  diff_data$pvalue < 0.05 & abs(diff_data$logFC) > 1.0,
  ifelse(diff_data$logFC > 1.0, "Up", "Down"),
  "None"
)

sig_genes <- diff_data[diff_data$classification != "None", ]

## 3. Data Export & Summary
write.csv(diff_data, file = file.path(data_dir, "all_genes_analysis_results.csv"), row.names = TRUE)
write.csv(sig_genes, file = file.path(data_dir, "DEGs_significant.csv"), row.names = TRUE)

cat("\n--- Differential Expression Summary (D-gal vs D-gal+BZBS) ---\n")
cat("Up-regulated genes  :", sum(diff_data$classification == "Up", na.rm = TRUE), "\n")
cat("Down-regulated genes:", sum(diff_data$classification == "Down", na.rm = TRUE), "\n")
cat("Total DEGs          :", sum(diff_data$classification != "None", na.rm = TRUE), "\n")
cat("Total input genes   :", nrow(diff_data), "\n")
cat("-------------------------------------------------------------\n")

## 4. Visualization
volcano_data <- diff_data[-log10(diff_data$pvalue) <= 300, ]

y_upper_limit <- 300
y_break_interval <- 50

p <- ggplot(volcano_data, aes(x = logFC, y = -log10(pvalue), color = classification)) +
  geom_point(alpha = 0.6, size = 2) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("Up" = "#FF6767", "Down" = "#3C8DAD", "None" = "gray60")) +
  xlab(expression(Log[2]*" fold change")) + 
  ylab(expression(-log[10](italic("p")))) +
  ggtitle("Volcano plot: D-gal vs D-gal+BZBS") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_vline(xintercept = c(-1.0, 1.0), linetype = "dashed", color = "black", linewidth = 0.7) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank()) +
  scale_y_continuous(limits = c(0, y_upper_limit),
                     breaks = seq(0, y_upper_limit, by = y_break_interval),
                     expand = expansion(mult = c(0, 0.05)))

output_pdf <- file.path(data_dir, "volcano_plot_filtered.pdf")
pdf(file = output_pdf, width = 10, height = 8)
print(p)
invisible(dev.off())

cat("Volcano plot saved to:", normalizePath(output_pdf, mustWork = FALSE), "\n")