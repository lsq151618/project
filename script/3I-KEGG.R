## 0. Initialization
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(tidyr)
})

## 1. Data Loading & Preprocessing
data_dir <- "./"
kegg_data <- read.csv(file.path(data_dir, "KEGG.csv"), header = TRUE, stringsAsFactors = FALSE)

if (!"classify" %in% colnames(kegg_data)) {
  stop("Missing 'classify' column in KEGG.csv")
}

## 2. Data Processing & Factor Reordering
data <- kegg_data %>%
  mutate(
    Adjusted_FoldEnrichment = ifelse(classify == "Activated", FoldEnrichment, -FoldEnrichment)
  ) %>%
  arrange(Adjusted_FoldEnrichment) %>%
  mutate(Description = factor(Description, levels = Description))

## 3. Visualization Configuration
custom_colors <- colorRampPalette(c("#3C8DAD", "white", "#FF6767"))(n = 100)
levelcolor <- c("Activated" = "#f1c7e1", "Suppressed" = "#99b3c6")

## 4. Plotting
p1 <- ggplot(data, aes(x = Adjusted_FoldEnrichment, y = Description)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_col(aes(fill = classify), width = 0.2) +
  geom_point(aes(size = Count, color = pvalue), shape = 16) +
  scale_size_continuous(range = c(4.5, 11.25), name = "Gene count") +
  scale_fill_manual(values = levelcolor) +
  scale_color_gradientn(
    colors = custom_colors,
    name = expression(italic("p") * " value")
  ) +
  scale_x_continuous(
    limits = c(-max(abs(data$Adjusted_FoldEnrichment)) * 1.1, 
               max(abs(data$Adjusted_FoldEnrichment)) * 1.1)
  ) +
  labs(
    x = "Fold enrichment", 
    y = NULL, 
    size = "Gene count", 
    fill = "Classification" 
  ) +
  theme(
    text = element_text(family = "sans"),
    panel.background = element_rect(fill = NA, color = NA), 
    panel.grid.major.x = element_line(color = "grey90", size = 0.2), 
    panel.grid.minor.x = element_blank(), 
    panel.grid.major.y = element_blank(), 
    panel.border = element_rect(colour = "Black", fill = NA, linewidth = 1), 
    axis.text = element_text(color = "Black", size = 10, face = "bold"), 
    axis.title = element_text(size = 12, face = "bold"), 
    legend.text = element_text(size = 12), 
    legend.title = element_text(size = 12) 
  )

## 5. Export Output
output_pdf <- file.path(data_dir, "KEGG_classified_lollipop_plot.pdf")
ggsave(output_pdf, p1, width = 10, height = 8)

cat(">>> KEGG analysis and visualization successfully completed!\n")