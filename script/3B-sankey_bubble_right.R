## 0. Configuration
data_dir <- "./" 

plot_width <- 10 
plot_height <- 8 

sankey_width <- 3 
dot_width <- 2 

sankey_text_size_gene <- 5 
sankey_text_size_pathway <- 5 
bubble_x_label <- "Gene ratio" 
sankey_x_label <- "Gene-pathway relationship" 

id_colors <- c(
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
  "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"
)

desc_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4",
  "#91D1C2", "#C8DC96", "#7E6148", "#B09C85", "#6A5ACD", "#A0522D"
)

bubble_size_range <- c(3, 8) 
spacer_freq <- 1 

bubble_low_color <- "#0000FF" 
bubble_high_color <- "#FF0000" 

windowsFonts(`Times New Roman` = windowsFont("Times New Roman"))
plot_font <- "Times New Roman" 

axis_order <- c("Gene", "Pathway") 
legend_justification <- c(0, 0.15) 

## 1. Package Initialization
source("install_dependencies.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggalluvial)
  library(patchwork)
  library(RColorBrewer)
  library(Cairo)
})

cols <- list(
  id = "geneID", 
  desc = "Description", 
  count = "Count", 
  sig = "pvalue", 
  ratio = "GeneRatio" 
)

## 2. Data Processing
enrichData <- read_tsv(file.path(data_dir, "enrich.tsv")) |>
  select(!!!cols) |>
  mutate(
    count = as.numeric(count),
    ratio = sapply(ratio, function(x) {
      if (grepl("/", x)) eval(parse(text = x)) else as.numeric(x)
    })
  )

## 3. Sankey Data Preparation
dfForLodes <- enrichData |>
  separate_rows(id, convert = TRUE, sep = "/") |>
  count(Gene = id, Pathway = desc, order = count, name = "Freq") |>
  arrange(desc(order))

sankeyData <- to_lodes_form(dfForLodes,
                            key = "axis",
                            axes = axis_order
) |>
  mutate(
    alluvium = as.character(alluvium),
    stratum = as.character(stratum)
  )

insert_spacer_nodes <- function(d, spacer_freq) {
  result <- d |>
    group_by(axis) |>
    summarise(
      nodes = list(unique(stratum)),
      spacers = list(paste0("spacer_", cur_group_id(), "_", seq_len(n_distinct(stratum) - 1))),
      .groups = "drop"
    ) |>
    mutate(
      levels = map2(nodes, spacers, ~ c(.x, .y)[order(c(seq_along(.x), seq_along(.y) + 0.5))])
    )
  
  spacer_rows <- result |>
    select(axis, spacers) |>
    unnest(cols = spacers) |>
    rename(stratum = spacers) |>
    mutate(Freq = spacer_freq, alluvium = stratum)
  
  bind_rows(d, spacer_rows) |>
    mutate(stratum = factor(stratum, levels = unlist(result$levels)))
}

sankeyData <- insert_spacer_nodes(sankeyData, spacer_freq)

## 4. Color Mapping
gene_colors <- setNames(
  colorRampPalette(id_colors)(length(unique(dfForLodes$Gene))),
  unique(dfForLodes$Gene)
)

pathway_colors <- setNames(
  colorRampPalette(desc_colors)(length(unique(dfForLodes$Pathway))),
  unique(dfForLodes$Pathway)
)

spacer_strata <- grep("spacer_", unique(sankeyData$stratum), value = TRUE)
spacer_colors <- setNames(rep("transparent", length(spacer_strata)), spacer_strata)

nodeColors <- c(gene_colors, pathway_colors, spacer_colors)

sankeyData <- sankeyData |>
  mutate(axis = factor(axis, levels = axis_order)) |>
  mutate(node_color = nodeColors[as.character(stratum)]) |>
  group_by(alluvium) |>
  mutate(
    to_node_name = lead(as.character(stratum), order_by = axis),
    to_node_name = ifelse(is.na(to_node_name), as.character(stratum), to_node_name)
  ) |>
  ungroup() |>
  mutate(flow_color = nodeColors[to_node_name])

## 5. Plot Generation

sankeyPlot <- ggplot(
  sankeyData,
  aes(
    x = axis, stratum = stratum, alluvium = alluvium,
    y = Freq, label = stratum
  )
) +
  geom_stratum(aes(fill = node_color), color = NA, width = 0.05) +
  geom_flow(
    aes(fill = flow_color),
    alpha = 0.3, width = 0.05,
    knot.pos = 0.3, color = "transparent"
  ) +
  geom_text(
    stat = "stratum",
    data = function(x) filter(x, axis == "Gene"),
    aes(label = ifelse(
      grepl("spacer_", as.character(after_stat(stratum))),
      "",
      as.character(after_stat(stratum))
    )),
    hjust = 1, nudge_x = -0.03,
    size = sankey_text_size_gene, family = plot_font, fontface = "bold"
  ) +
  geom_text(
    stat = "stratum",
    data = function(x) filter(x, axis == "Pathway"),
    aes(label = str_wrap(
      ifelse(
        grepl("spacer_", as.character(after_stat(stratum))),
        "",
        as.character(after_stat(stratum))
      ),
      50
    )),
    hjust = 1, nudge_x = -0.03,
    size = sankey_text_size_pathway, family = plot_font, fontface = "bold"
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_identity() +
  guides(fill = "none") +
  theme_void() +
  labs(x = sankey_x_label) +
  theme(
    text = element_text(family = plot_font, face = "bold", color = "black"),
    axis.title.x = element_text(margin = margin(t = 2), size = 16)
  )

sankeyPlotData <- ggplot_build(sankeyPlot)

rightNodes <- sankeyPlotData$data[[1]] |>
  filter(x == max(x)) |>
  mutate(
    node_name = as.character(stratum),
    node_ymin = ymin,
    node_ymax = ymax,
    node_center_y = (ymin + ymax) / 2
  ) |>
  filter(!grepl("spacer_", node_name)) |>
  select(node_name, node_center_y, ymin, ymax)

dotData <- enrichData |>
  distinct(desc, .keep_all = TRUE) |>
  mutate(GeneRatio = ratio) |>
  left_join(rightNodes, by = c("desc" = "node_name"))

dotPlot <- ggplot(dotData, aes(x = GeneRatio, y = node_center_y, color = -log10(sig))) +
  geom_point(aes(size = count), stroke = 0.5) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0.005)) +
  scale_color_gradient(low = bubble_low_color, high = bubble_high_color) +
  scale_size_continuous(
    range = bubble_size_range, 
    name = "Count", 
    breaks = function(x) unique(round(seq(min(x), max(x), length.out = 3)))
  ) +
  guides(
    color = guide_colorbar(order = 1, barwidth = 0.8, barheight = 4),
    size = guide_legend(order = 2, keywidth = 0.8, keyheight = 0.8)
  ) +
  labs(
    size = "Count",
    color = expression(bold("-log")[bold("10")] * bold("(") * bolditalic("p") * bold(")")), 
    x = bubble_x_label
  ) +
  theme_void() +
  theme(
    text = element_text(family = plot_font, face = "bold", color = "black"),
    axis.text.x = element_text(margin = margin(t = 4), size = 12),
    axis.title.x = element_text(margin = margin(t = 6), size = 16),
    axis.ticks.x = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length.x = unit(0.1, "cm"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.position = "right",
    legend.justification = legend_justification,
    legend.margin = margin(0, 0, 0, 0)
  )

## 6. Plot Assembly and Output
yRange <- sankeyPlotData$layout$panel_params[[1]]$y.range

sankeyPlot <- sankeyPlot +
  coord_cartesian(clip = "off", ylim = yRange) +
  theme(plot.margin = margin(0, 0, 0, 2.5, "cm"))

dotPlot <- dotPlot +
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = max(rightNodes$ymax),
           fill = NA, color = "black", linewidth = 0.5
  ) +
  coord_cartesian(ylim = yRange)

combinedPlot <- (sankeyPlot + dotPlot) +
  plot_layout(widths = c(sankey_width, dot_width))

ggsave(file.path(data_dir, "sankey_bubble_right.pdf"), combinedPlot,
       width = plot_width, height = plot_height, device = cairo_pdf
)

cat("Output saved to directory: ", normalizePath(data_dir, mustWork = FALSE), "\n", sep = "")
cat("  - sankey_bubble_right.pdf\n")