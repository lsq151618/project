## 0. Initialization
suppressPackageStartupMessages({
  library(VennDiagram)
  library(openxlsx)
  library(grid)
})

## 1. Configuration
data_dir <- "./"

con_file <- file.path(data_dir, "Con vs Dgal.csv")
bzbs_file <- file.path(data_dir, "BZBS vs Dgal.csv")

output_pdf <- file.path(data_dir, "Con_BZBS_VennDiagram.pdf")
output_common_genes <- file.path(data_dir, "Con_BZBS_Intersect.xlsx")

clean_genes <- function(genes) {
  as.character(unique(na.omit(genes[genes != ""])))
}

## 2. Data Processing
con_data <- read.csv(con_file, check.names = FALSE, stringsAsFactors = FALSE)
bzbs_data <- read.csv(bzbs_file, check.names = FALSE, stringsAsFactors = FALSE)

con_genes <- clean_genes(con_data$`gene name`)
bzbs_genes <- clean_genes(bzbs_data$`gene name`)

common_genes <- intersect(con_genes, bzbs_genes)

write.xlsx(data.frame(Intersect_Genes = common_genes), 
           file = output_common_genes, 
           rowNames = FALSE)

if (length(con_genes) == 0 && length(bzbs_genes) == 0) {
  stop("Gene lists are empty. Cannot plot Venn diagram.")
}

## 3. Visualization
pdf(output_pdf, width = 8, height = 8)

if (length(con_genes) > 0 && length(bzbs_genes) > 0) {
  custom_colors <- c("#FF8396", "#7DCFFF")
  
  venn.plot <- venn.diagram(
    x = list(
      "Ctrl vs D-gal" = con_genes,
      "D-gal+BZBS vs D-gal" = bzbs_genes
    ),
    filename = NULL,
    col = "#696969",          
    lty = "dashed",           
    lwd = 0.75,               
    fill = custom_colors,     
    alpha = 0.6,
    cex = 1.4,                
    cat.cex = 1.2,            
    cat.fontface = "bold",
    cat.col = custom_colors,  
    cat.pos = c(-15, 15),   
    cat.dist = c(0.06, 0.06),
    margin = 0.12,            
    force.unique = TRUE,    
    euler.d = FALSE,          
    scaled = FALSE          
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

cat("Ctrl vs D-gal gene count        :", length(con_genes), "\n")
cat("D-gal+BZBS vs D-gal gene count  :", length(bzbs_genes), "\n")
cat("Total intersecting genes        :", length(common_genes), "\n")

if (length(common_genes) > 0) {
  cat("\nTop 10 intersecting genes:\n")
  show_n <- min(10, length(common_genes))
  cat(paste(common_genes[1:show_n], collapse = ", "))
  if (length(common_genes) > 10) cat("\n... (Total", length(common_genes), "genes)\n")
  cat("\n")
} else {
  cat("\nNo intersecting genes found between the two groups.\n")
}