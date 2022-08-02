# Install and load required packages
# Singularity mounts the HOME directory of the user by default.
# This leads to R using libraries installed in the HOME folder of the user. 
# This line assures that libraries from the Singularity container are used instead.
  .libPaths("/usr/local/lib/R/site-library")
      
list.of.packages <- c("openxlsx","dplyr","stringr",
                      "ggplot2","grid","plotly","crayon",
                      "htmlwidgets","scales","withr",
                      "data.table","BiocManager")
bioconductor.packages <- c("Rsubread","ComplexHeatmap")
# General packages
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
for(n in new.packages){
  install.packages(n, repos = "http://cran.us.r-project.org")
}
# Bioconductor packages
for(b in bioconductor.packages){
  if (!requireNamespace(b, quietly = TRUE))
	BiocManager::install(b)
}
