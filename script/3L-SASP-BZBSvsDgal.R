## 0. Initialization
pkgs <- c("readxl", "circlize")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(readxl)
library(circlize)

## 1. Data Loading & Preprocessing
data_dir <- "./"
df <- as.data.frame(read_excel(file.path(data_dir, "BZBSvsDgal.xlsx"), sheet = 1))
rownames(df) <- df[[1]]
mat_raw <- as.matrix(df[, -1, drop = FALSE])
storage.mode(mat_raw) <- "numeric"

mat <- t(scale(t(mat_raw)))

## 2. Sample & Gene Ordering
colnames(mat) <- trimws(colnames(mat))
dgal_cols <- grep("D-gal", colnames(mat), value = TRUE, ignore.case = TRUE)
bzbs_cols <- grep("BZBS", colnames(mat), value = TRUE, ignore.case = TRUE)
mat <- mat[, c(dgal_cols, bzbs_cols), drop = FALSE]

hc <- hclust(dist(mat), method = "complete")
ordered_genes <- rownames(mat)[hc$order]
target_genes <- c("SERPINE1", "IL6", "CXCL3")
idx <- which(ordered_genes %in% target_genes)
if(length(idx) > 0) ordered_genes[idx] <- rev(ordered_genes[idx])
mat <- mat[ordered_genes, , drop = FALSE]

## 3. Visualization Constants
min_val <- floor(min(mat) * 10) / 10
max_val <- ceiling(max(mat) * 10) / 10
col_fun <- colorRamp2(c(min_val, 0, max_val), c("#8CA1E0", "#FFFFFF", "#F48B8B"))

GAP_DEG <- 30; GAP_POS <- 2; START_DEG <- 90
TRACK_H <- 0.40; CEX_GENE <- 1.0; CEX_SAMPLE <- 1.4
DEND_RLEAF <- 0.36; DEND_RIN <- 0.17

## 4. Helper Functions
draw_radial_dend <- function(hc, leaf_angle_deg, r_leaf, r_inner, lwd = 0.8) {
  rad <- function(d) d * pi / 180; th <- rad(leaf_angle_deg); hmax <- max(hc$height)
  node_th <- node_r <- numeric(nrow(hc$merge))
  getpos <- function(idx) if (idx < 0) c(th[-idx], r_leaf) else c(node_th[idx], node_r[idx])
  draw_arc <- function(t1, t2, r) {
    d <- t2 - t1; while(d > pi) d <- d - 2*pi; while(d < -pi) d <- d + 2*pi
    tt <- seq(t1, t1 + d, length.out = 80); lines(r * cos(tt), r * sin(tt), lwd = lwd)
    t1 + d / 2
  }
  for (i in seq_len(nrow(hc$merge))) {
    a <- getpos(hc$merge[i, 1]); b <- getpos(hc$merge[i, 2])
    ri <- r_leaf - (hc$height[i] / hmax) * (r_leaf - r_inner)
    segments(a[2]*cos(a[1]), a[2]*sin(a[1]), ri*cos(a[1]), ri*sin(a[1]), lwd=lwd)
    segments(b[2]*cos(b[1]), b[2]*sin(b[1]), ri*cos(b[1]), ri*sin(b[1]), lwd=lwd)
    mid <- draw_arc(a[1], b[1], ri); node_th[i] <- mid; node_r[i] <- ri
  }
}

draw_circos <- function() {
  split <- factor(rownames(mat), levels = rownames(mat))
  gaps <- rep(2, nrow(mat)); gaps[GAP_POS] <- GAP_DEG
  circos.clear(); circos.par(start.degree = START_DEG, gap.after = gaps, track.height = TRACK_H)
  circos.heatmap(mat, col = col_fun, split = split, cluster = FALSE, dend.side = "none", bg.border = NA)
  
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    if (CELL_META$sector.index == levels(split)[GAP_POS]) {
      cn <- colnames(mat); n <- length(cn)
      y_pos <- n - seq_len(n) + 0.5
      circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1.5, "mm"),
                  y_pos, cn, cex = CEX_SAMPLE, adj = c(0, 0.5), facing = "inside", niceFacing = TRUE)
    }
  }, bg.border = NA)
  
  for(sn in get.all.sector.index()) {
    set.current.cell(sector.index = sn, track.index = 1)
    circos.text(CELL_META$xcenter, CELL_META$ycenter, labels = sn, facing = "bending.inside", 
                niceFacing = TRUE, cex = CEX_GENE * 1.5, col = "black")
  }
  
  sec_center <- sapply(get.all.sector.index(), function(s) (get.cell.meta.data("cell.start.degree", s) + get.cell.meta.data("cell.end.degree", s)) / 2)
  leaf_angle <- as.numeric(sec_center[rownames(mat)])
  circos.clear()
  draw_radial_dend(hc, leaf_angle, DEND_RLEAF, DEND_RIN)
  
  xL <- -0.040; xR <- -0.005; yB <- -0.13; yT <- 0.09
  cols <- col_fun(seq(min_val, max_val, length.out = 200))
  for (i in 1:200) rect(xL, yB + (i-1)*(yT-yB)/200, xR, yB + i*(yT-yB)/200, col = cols[i], border = NA)
  rect(xL, yB, xR, yT, border = "black", lwd = 0.6)
  ticks <- ceiling(min_val):floor(max_val)
  ty <- yB + (ticks - min_val) / (max_val - min_val) * (yT - yB)
  segments(xR, ty, xR + 0.012, ty, lwd = 0.6)
  text(xR + 0.017, ty, labels = ticks, adj = c(0, 0.5), cex = CEX_SAMPLE * 0.8)
}

## 5. Export Output
cairo_pdf(file.path(data_dir, "SASP-BZBSvsDgal.pdf"), width = 7, height = 7, pointsize = 10.5)
draw_circos(); dev.off()

png(file.path(data_dir, "SASP-BZBSvsDgal.png"), width = 2100, height = 2100, res = 300, pointsize = 10.5)
draw_circos(); dev.off()

cat(">>> Execution completed successfully.\n")