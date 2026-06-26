#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://mirrors.westlake.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor/")

log_msg <- function(level, fmt, ...) {
  cat(sprintf("[%s] %s ", level, format(Sys.time(), "%H:%M:%S")),
      sprintf(fmt, ...), "\n",
      sep = ""
  )
}

cran_packages <- c(
  "Cairo",
  "ggalluvial",
  "patchwork",
  "RColorBrewer",
  "tidyverse"
)

bioc_packages <- c()
github_packages <- c()
version_packages <- c()

ipk <- unique(rownames(installed.packages(fields = "Package")))

need <- function(p, src = "CRAN") {
  miss <- setdiff(p, ipk)
  skip <- intersect(p, ipk)
  if (length(skip)) {
    log_msg("INFO", "Skipping installed %s package(s): %s", src, toString(skip))
  }
  miss
}

install_cran <- function(pkgs) {
  if (!length(pkgs)) return(invisible())
  
  tryCatch({
    install.packages(pkgs, dependencies = TRUE)
    log_msg("SUCCESS", "CRAN batch installation complete: %s", toString(pkgs))
  }, error = function(e) {
    log_msg("WARN", "Batch installation failed, retrying individually: %s", e$message)
    for (p in pkgs) {
      tryCatch(install.packages(p, dependencies = TRUE),
        error = function(e) log_msg("ERROR", "%s installation failed: %s", p, e$message)
      )
    }
  })
}

install_bioc <- function(pkgs) {
  if (!length(pkgs)) return(invisible())
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install(pkgs, update = FALSE, ask = FALSE)
  log_msg("SUCCESS", "Bioconductor installation complete: %s", toString(pkgs))
}

install_github <- function(repos) {
  if (!length(repos)) return(invisible())
  
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }

  for (repo in repos) {
    pkg <- basename(repo)
    if (pkg %in% ipk) {
      log_msg("INFO", "Skipping installed GitHub package: %s", pkg)
      next
    }
    tryCatch({
      devtools::install_github(repo, upgrade = "never")
      log_msg("SUCCESS", "GitHub package installation complete: %s", pkg)
      ipk <<- c(ipk, pkg)
    }, error = function(e) {
      log_msg("ERROR", "%s installation failed: %s", pkg, e$message)
    })
  }
}

install_version <- function(named_vec) {
  if (!length(named_vec)) return(invisible())
  
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }

  for (pkg in names(named_vec)) {
    tgt_ver <- named_vec[[pkg]]
    local_ver <- tryCatch(as.character(utils::packageVersion(pkg)),
      error = function(e) NA
    )
    
    if (!is.na(local_ver) && local_ver == tgt_ver) {
      log_msg("INFO", "Skipping installed versioned package: %s %s", pkg, tgt_ver)
      next
    }
    tryCatch({
      remotes::install_version(pkg, version = tgt_ver, upgrade = "never")
      log_msg("SUCCESS", "%s version %s installation complete", pkg, tgt_ver)
    }, error = function(e) {
      log_msg("ERROR", "%s version %s installation failed: %s", pkg, tgt_ver, e$message)
    })
  }
}

log_msg("INFO", "========== Starting Dependency Installation ==========")
install_cran(need(cran_packages, "CRAN"))
install_bioc(need(bioc_packages, "Bioconductor"))
install_github(need(github_packages, "GitHub"))
install_version(version_packages)
log_msg("INFO", "========== Dependency Installation Finished ==========")