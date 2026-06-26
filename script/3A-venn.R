## 0. Initialization
suppressPackageStartupMessages({
  library(VennDiagram)
  library(openxlsx)
  library(grid)
})

## 1. Configuration
data_dir <- "./" 

aged_file <- file.path(data_dir, "Aged.xlsx")
bzbs_file <- file.path(data_dir, "BZBS.xlsx")

output_pdf <- file.path(data_dir, "Aged_BZBS_VennDiagram.pdf")
output_common_genes <- file.path(data_dir, "Aged_BZBS_Intersect.xlsx")

clean_genes <- function(genes) {
  as.character(unique(na.omit(genes[genes != ""])))
}

## 2. Data Processing
aged_genes <- clean_genes(read.xlsx(aged_file, colNames = FALSE)[[1]])
bzbs_genes <- clean_genes(read.xlsx(bzbs_file, colNames = FALSE)[[1]])

common_genes <- intersect(aged_genes, bzbs_genes)

write.xlsx(data.frame(Intersect_Genes = common_genes), 
           file = output_common_genes, 
           rowNames = FALSE)

if (length(aged_genes) == 0 && length(bzbs_genes) == 0) {
  stop("Error: All gene lists are empty. Cannot plot Venn diagram.")
}

## 3. Visualization
pdf(output_pdf, width = 8, height = 8)

if (length(aged_genes) > 0 && length(bzbs_genes) > 0) {
  venn.plot <- venn.diagram(
    x = list(
      "Aged" = aged_genes,
      "BZBS" = bzbs_genes
    ),
    filename = NULL,
    col = "#696969",        
    lty = "dashed",         
    lwd = 0.75,             
    fill = c("#FFB6C1", "#B0E2FF"), 
    alpha = 0.6,
    cex = 1.4,              
    cat.cex = 1.4,          
    cat.fontface = "bold",
    cat.col = c("#FFB6C1", "#B0E2FF"), 
    cat.pos = c(-20, 20),  
    cat.dist = c(0.05, 0.05),
    margin = 0.1,
    force.unique = TRUE,    
    euler.d = FALSE,        
    scaled = FALSE,
    print.mode = c("raw", "percent"), 
    sigdigs = 3
  )
  grid.draw(venn.plot)
} else {
  grid.text("Cannot plot Venn diagram:\nOne or more gene lists are empty.", 
            gp = gpar(fontsize = 16, col = "red"))
}

invisible(dev.off())

## 4. Execution Summary
cat("--- Execution Complete ---\n")
cat("Venn diagram saved to:", normalizePath(output_pdf, mustWork = FALSE), "\n")
cat("Intersect genes saved to:", normalizePath(output_common_genes, mustWork = FALSE), "\n\n")

cat("Aged gene count:", length(aged_genes), "\n")
cat("BZBS gene count:", length(bzbs_genes), "\n")
cat("Total intersecting genes:", length(common_genes), "\n")

if (length(common_genes) > 0) {
  cat("\nTop 10 intersecting genes:\n")
  show_n <- min(10, length(