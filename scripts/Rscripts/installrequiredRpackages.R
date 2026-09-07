# Install and load required R packages (pinned versions).
#
# Singularity mounts the user's HOME by default, which makes R pick up packages
# from ~/. This line forces use of the container's library instead.
  .libPaths("/usr/local/lib/R/site-library")

# ---------------------------------------------------------------------------
# CRAN packages  (dated snapshot -> deterministic dependency resolution)
# ---------------------------------------------------------------------------
  CRAN_SNAPSHOT <- "https://packagemanager.posit.co/cran/2025-08-07"
  options(repos = c(CRAN = CRAN_SNAPSHOT))

  cran_pkgs <- c(
    openxlsx     = "4.2.8",
    dplyr        = "1.1.4",
    stringr      = "1.5.1",
    ggplot2      = "3.5.2",
    plotly       = "4.11.0",
    crayon       = "1.5.3",
    htmlwidgets  = "1.6.4",
    scales       = "1.4.0",
    withr        = "3.0.2",
    data.table   = "1.17.8",
    BiocManager  = "1.30.26"
  )

  if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

  for (pkg in names(cran_pkgs)) {
    want <- cran_pkgs[[pkg]]
    have <- tryCatch(as.character(packageVersion(pkg)), error = function(e) NA)
    if (is.na(have) || have != want)
      remotes::install_version(pkg, version = want, upgrade = "never")
  }

# ---------------------------------------------------------------------------
# Bioconductor packages  (pinned release 3.21)
# ---------------------------------------------------------------------------
  bioc_pkgs <- c(
    Rsubread       = "2.22.1",
    ComplexHeatmap = "2.24.1"
  )

  BiocManager::install(version = "3.21", update = FALSE, ask = FALSE)
  for (pkg in names(bioc_pkgs)) {
    if (!requireNamespace(pkg, quietly = TRUE))
      BiocManager::install(pkg, version = "3.21", update = FALSE, ask = FALSE)
  }

# ---------------------------------------------------------------------------
# Record the exact environment for the record
# ---------------------------------------------------------------------------
  writeLines(capture.output(sessionInfo()), "/scripts/Rscripts/sessionInfo.txt")
  cat("R package installation complete (pinned versions).\n")
