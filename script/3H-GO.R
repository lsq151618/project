## 0. Initialization
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(tidyr)
})

## 1. Data Loading & Preprocessing
data_dir <- "./"
go_data <- read.csv(file.path(data_dir, "GO.csv"), header = TRUE, stringsAsFactors = FALSE)

if (!"classify" %in% colnames(go_data)) {
  stop("Missing 'classify' column in GO.csv")
}
if (!"ONTOLOGY" %in% colnames(go_data)) {
  stop("Missing 'ONTOLOGY' column in GO.csv")
}

go_data$classify <- ifelse(go_data$classify == "up", "Activated", 
                           ifelse(go_data$classify == "down", "Suppressed", go_data$classify))

## 2. Data Processing & Factor Reordering
data <- go_data %>%
  mutate(
    Adjusted_FoldEnrichment = ifelse(classify == "Activated", FoldEnrichment, -FoldEnrichment)
  ) %>%
  arrange(Adjusted_FoldEnrichment) %>%
  mutate(Description = factor(Description, levels = Description)) %>%
  mutate(Shape = case_when(
    ONTOLOGY == "BP" ~ 16,
    ONTOLOGY == "CC" ~ 17,
    ONTOLOGY == "MF" ~ 15,
    TRUE ~ 16
  ))

## 3. Visualization Configuration
custom_colors <- colorRampPalette(c("#FF6767", "white", "#3C8DAD"))(n = 100)
levelcolor <- c("Activated" = "#f1c7e1", "Suppressed" = "#99b3c6")

## 4. Plotting
p1 <- ggplot(data, aes(x = Adjusted_FoldEnrichment, y = Description)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_col(aes(fill = classify), width = 0.2) +
  geom_point(aes(size = Count, color = pvalue, shape = ONTOLOGY)) +
  scale_size_continuous(range = c(3, 7.5), name = "Gene count") +
  scale_fill_manual(values = levelcolor) +
  scale_color_gradientn(
    colors = custom_colors,
    name = expression(italic("p") * " value")   
  ) +
  scale_shape_manual(
    values = c("BP" = 16, "CC" = 17, "MF" = 15),
    name = "GO category"
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
    legend.text = element_text(size = 10),  
    legend.title = element_text(size = 12)  
  )

## 5. Export Output
output_pdf <- file.path(data_dir, "GO_classified_lollipop_plot.pdf")
ggsave(output_pdf, p1, width = 10, height = 8)

cat(">>> GO analysis and visualization successfully completed!\n")