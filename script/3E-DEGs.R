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
expr_data <- as.data.frame(read_excel(diff_file_path))

required_cols <- c("gene name", "C1", "C2", "C3", "G1", "G2", "G3")
if (!all(required_cols %in% colnames(expr_data))) {
  stop("Missing required columns: 'gene name', 'C1','C2','C3','G1','G2','G3'")
}

if (any(duplicated(expr_data$`gene name`))) {
  expr_data <- expr_data[!duplicated(expr_data$`gene name`), ]
}

rownames(expr_data) <- expr_data$`gene name`

expr_cols <- c("C1", "C2", "C3", "G1", "G2", "G3")
for (col in expr_cols) {
  if (!is.numeric(expr_data[[col]])) {
    expr_data[[col]] <- as.numeric(as.character(expr_data[[col]]))
  }
}

## 2. Statistical Analysis
expr_data$Con_mean <- rowMeans(expr_data[, c("C1", "C2", "C3")], na.rm = TRUE)
expr_data$Dgal_mean <- rowMeans(expr_data[, c("G1", "G2", "G3")], na.rm = TRUE)

pseudo_count <- 0.001
expr_data$logFC <- log2((expr_data$Con_mean + pseudo_count) / (expr_data$Dgal_mean + pseudo_count))

calculate_pvalue <- function(con_values, dgal_values) {
  if (length(na.omit(con_values)) < 2 || length(na.omit(dgal_values)) < 2) return(NA)
  tryCatch({
    t.test(con_values, dgal_values, var.equal = TRUE)$p.value
  }, error = function(e) NA)
}

expr_data$pvalue <- sapply(1:nrow(expr_data), function(i) {
  calculate_pvalue(as.numeric(expr_data[i, c("C1", "C2", "C3")]),
                   as.numeric(expr_data[i, c("G1", "G2", "G3")]))
})

expr_data$padj <- p.adjust(expr_data$pvalue, method = "fdr")

expr_data$classification <- ifelse(
  expr_data$pvalue < 0.05 & abs(expr_data$logFC) > 1.0,
  ifelse(expr_data$logFC > 1.0, "Up", "Down"),
  "None"
)

results_data <- data.frame(
  `gene name` = expr_data$`gene name`,
  Con_mean = expr_data$Con_mean,
  Dgal_mean = expr_data$Dgal_mean,
  logFC = expr_data$logFC,
  pvalue = expr_data$pvalue,
  padj = expr_data$padj,
  classification = expr_data$classification,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
rownames(results_data) <- expr_data$`gene name`

## 3. Data Export & Summary
write.csv(results_data, file = file.path(data_dir, "DEGs_Con_vs_Dgal_calculated_results.csv"), row.names = TRUE)

sig_genes <- results_data[results_data$classification != "None", ]
write.csv(sig_genes, file = file.path(data_dir, "DEGs_significant_Con_vs_Dgal.csv"), row.names = TRUE)

cat("\n--- Differential Expression Summary (Ctrl vs D-gal) ---\n")
cat("Up-regulated genes  :", sum(results_data$classification == "Up", na.rm = TRUE), "\n")
cat("Down-regulated genes:", sum(results_data$classification == "Down", na.rm = TRUE), "\n")
cat("Total DEGs          :", sum(results_data$classification != "None", na.rm = TRUE), "\n")
cat("Total input genes   :", nrow(results_data), "\n")
cat("-------------------------------------------------------\n")

## 4. Visualization
volcano_data <- results_data[-log10(results_data$pvalue) <= 10, ]

p <- ggplot(volcano_data, aes(x = logFC, y = -log10(pvalue), color = classification)) +
  geom_point(alpha = 0.6, size = 2) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("Up" = "#FF6767", "Down" = "#3C8DAD", "None" = "gray60")) +
  xlab(expression(Log[2]*" fold change")) + 
  ylab(expression(-log[10](italic("p")))) +
  ggtitle("Volcano plot: Ctrl vs D-gal") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_vline(xintercept = c(-1.0, 1.0), linetype = "dashed", color = "black", linewidth = 0.7) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank()) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))

output_pdf <- file.path(data_dir, "volcano_plot_Con_vs_Dgal.pdf")
pdf(file = output_pdf, width = 10, height = 8)
print(p)
invisible(dev.off())

cat("Volcano plot saved to:", normalizePath(output_pdf, mustWork = FALSE), "\n")