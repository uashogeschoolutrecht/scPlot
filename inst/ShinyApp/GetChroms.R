#' Get a list of the chromosomes of a genome
#'
#' @param g_url         Url to the ".gtf" file of the researched genome on the ncbi website.
#' @param gen_name      Name to store the downloaded files, and title name of the plot.
#' @param dl_folder     Path to the folder to download the files in.
#'
#' @return List of chromosomes
#'
#' @examples
#' GetChroms(g_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/277/835/GCF_002277835.1_ASM227783v1/GCF_002277835.1_ASM227783v1_genomic.gtf.gz",
#'            dl_folder = "C:/Users/Downloads",
#'            gen_name = Flavobacterium_ir1)
#'
#' @export

GetChroms <-
  function(g_url,
           dl_folder,
           gen_name = NULL) {

    library(R.utils)
    library(tidyverse)

    file_url <- g_url
    pathD <- dl_folder

    pathD2 <- paste0(pathD, "/", gen_name)

    pathD5 <- paste0(pathD2, "/", gen_name, ".gtf.gz")
    pathD6 <- paste0(pathD2, "/", gen_name, ".gtf")

    pathD7 <- paste0(pathD2, "/", gen_name, ".gff.gz")
    pathD8 <- paste0(pathD2, "/", gen_name, ".gff")

    is_gtf <- stringr::str_detect(string = g_url, pattern = ".gtf")
    if (is_gtf == TRUE) {
      pathD9 <- pathD5
      pathD10 <- pathD6
    } else if (is_gtf == FALSE) {
      pathD9 <- pathD7
      pathD10 <- pathD8
    }


    FE <- file.exists(pathD10)

    if (FE == FALSE) {

      dir.create(pathD2)

      download.file(url = file_url,
                    destfile = pathD9,
                    cacheOK = TRUE,
                    method = "curl",
                    quiet = TRUE,
                    mode = "a")

      # Unzip the files

      R.utils::gunzip(filename = pathD9, destname = pathD10, remove = TRUE)

    } else if (FE == TRUE) {}

    gendat <-
      readr::read_delim(pathD10,
                        delim = "\t",
                        comment = "#",
                        col_names = FALSE)
    gendat_group <- gendat %>% dplyr::select(X1, X4, X5)
    gendat_group <- dplyr::group_by(gendat_group, X1)
    gendat_group <- dplyr::summarise(gendat_group)
    as.character(gendat_group$X1)

    print_chr <- paste(gendat_group$X1, sep = ", ")

    return(chrom = print_chr)

  }
