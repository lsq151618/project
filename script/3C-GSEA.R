## 0. Initialization
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(openxlsx)
})

## 1. Configuration & Data Loading
data_dir <- "./"

data <- read.table(file.path(data_dir, "gene.txt"), 
                   header = TRUE, 
                   sep = "\t", 
                   stringsAsFactors = FALSE)

colnames(data) <- c("Gene", "Pi_score")

## 2. ID Conversion & Ranking
gene_list <- bitr(data$Gene, 
                  fromType = "SYMBOL",
                  toType = "ENTREZID", 
                  OrgDb = org.Hs.eg.db)

data <- merge(data, gene_list, by.x = "Gene", by.y = "SYMBOL")

gene_ranks <- data$Pi_score
names(gene_ranks) <- data$ENTREZID
gene_ranks <- sort(gene_ranks, decreasing = TRUE) 

## 3. GSEA Analysis
gsea_result <- gseKEGG(
  geneList = gene_ranks,   
  organism = "hsa",        
  keyType = "ncbi-geneid", 
  exponent = 0.5,          
  minGSSize = 10,          
  maxGSSize = 500,         
  pvalueCutoff = 0.5,      
  verbose = FALSE          
)

gsea_result <- setReadable(gsea_result,
                           OrgDb = 'org.Hs.eg.db',
                           keyType = 'ENTREZID')

## 4. Statistics Extraction
p_val <- gsea_result@result[gsea_result@result$ID == "hsa04152", "pvalue"]
adj_p_val <- gsea_result@result[gsea_result@result$ID == "hsa04152", "p.adjust"]

cat("\n=== AMPK Signaling Pathway Statistics ===\n")
cat(sprintf("Nominal P-value : %g\n", p_val))
cat(sprintf("Adjusted P-value: %g\n", adj_p_val))
cat("=========================================\n\n")

## 5. Visualization
g <- gseaplot2(gsea_result,
               geneSetID = "hsa04152", 
               title = "AMPK signaling pathway", 
               subplots = 1:3,
               ES_geom = 'line',
               pvalue_table = FALSE)

g[[1]] <- g[[1]] + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Running enrichment score")

g[[3]] <- g[[3]] + 
  ylab("Ranked list metric") +
  xlab("Rank in ordered dataset")

## 6. Output Generation
output_file <- file.path(data_dir, "AMPK_signaling_pathway.pdf")

ggsave(filename = output_file, 
       plot = g, width = 6, height = 5, dpi = 300)

cat(">>> GSEA plot successfully saved to PDF format!\n")