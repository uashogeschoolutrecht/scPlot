## helper script for Shiny App
#library(BiocManager)
#options(repos = BiocManager::repositories())
library(purrr)

functions <- list.files(pattern = "Ge", full.names = TRUE)

purrr::map(functions, source)

