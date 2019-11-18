#' Make Genomic Plots Function
#'
#' @param g_url         Url to the ".gtf" file of the researched genome on the ncbi website.
#'                      when the url is the same as the .fna file (exept for the ".fna"), only the fnaUrl is enough.
#' @param f_url         Url to the ".fna" file of the researched genome on the ncbi website.
#' @param chr_id        Chromosome ID name of the researched chromosome, as named in
#'                      the gff and gtf files. As character. Default is the first chromosome.
#' @param gen_name      Name to store the downloaded files, and title name of the plot.
#'                      The name of the organism without spaces is recommended. As character.
#' @param start         Base start of plot.
#' @param end           Base end of plot.
#' @param gene_only     Show only the Gene or also start and stop codons? Logical.
#' @param show_gene_id  Show the gene ID instead of the product in the plot. Logical.
#'                      The gene ID is needed to get the sequence of a gene in function GetSequenceNow.
#' @param dl_folder     Path to the folder to download the files in.
#'
#' @return A plot of the inserted genome files.
#'
#' @examples
#' MakeGenPlot(g_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/277/835/GCF_002277835.1_ASM227783v1/GCF_002277835.1_ASM227783v1_genomic.gtf.gz",
#'            dl_folder = "C:/Users/boazk/Downloads",
#'            chr_id = "NZ_NQOT01000381.1",
#'            gen_name = "Flavobacterium_ir1",
#'            start = 80000,
#'            end = 90000,
#'            show_gene_id = FALSE,
#'            f_url = NULL)
#'
#' @export

MakeGenPlot <-
  function(g_url = NULL,
           f_url = NULL,
           dl_folder = NULL,
           chr_id = NULL,
           gen_name = NULL,
           start = 1,
           end = 10000,
           show_gene_id = FALSE) {

    library(tidyverse)
    library(Biostrings)
    library(IRanges)
    library(GenomicRanges)
    library(S4Vectors)
    library(Gviz)
    library(rtracklayer)
    library(biomaRt)
    library(ggbio)
    library(GenomicFeatures)
    library(seqRFLP)
    library(R.utils)
    library(stringr)

    # With the help of this tutorial:
    # http://www.sthda.com/english/wiki/gviz-visualize-genomic-data


    # Get the gff and gtf files and make them GRanges objects.
    # Choose between only CDS and CDS, start and stopcodon display
    # Download files

    file_url <- g_url
    pathD <- dl_folder

    pathD2 <- paste0(pathD, "/", gen_name)

    pathD3 <- paste0(pathD2, "/", gen_name, ".fna.gz")
    pathD4 <- paste0(pathD2, "/", gen_name, ".fna")

    pathD5 <- paste0(pathD2, "/", gen_name, ".gtf.gz")
    pathD6 <- paste0(pathD2, "/", gen_name, ".gtf")

    pathD7 <- paste0(pathD2, "/", gen_name, ".gff.gz")
    pathD8 <- paste0(pathD2, "/", gen_name, ".gff")

    urlSame <- str_sub(file_url, end = -8)
    urlSame2 <- paste0(urlSame, ".fna.gz")

    is_gtf <- stringr::str_detect(string = g_url, pattern = ".gtf")

    # Is it a .gtf or .gff file?

    if (is_gtf == TRUE) {
      pathD9 <- pathD5
      pathD10 <- pathD6
    } else if (is_gtf == FALSE) {
      pathD9 <- pathD7
      pathD10 <- pathD8
    }


    if (is.null(f_url)) {
      fnaLink <- urlSame2
    } else if (!is.null(f_url)) {
      fnaLink <- f_url
    }

    file_urlG <- fnaLink

    # If the files already exist. Don't download them again and inform user.
    # If they don't exist. Download them and inform user.

    FE <- file.exists(pathD4)

    if (FE == FALSE) {

      dir.create(pathD2)

      download.file(url = file_url,
                    destfile = pathD9,
                    cacheOK = TRUE,
                    method = "curl",
                    quiet = TRUE,
                    mode = "a")

      download.file(url = file_urlG,
                    destfile = pathD3,
                    cacheOK = TRUE,
                    method = "curl",
                    quiet = TRUE,
                    mode = "a")

      # Unzip the files

      R.utils::gunzip(filename = pathD9, destname = pathD10, remove = TRUE)
      R.utils::gunzip(filename = pathD3, destname = pathD4, remove = TRUE)

      print_info <- print(paste0("Files were downloaded to folder: ", pathD2))

    } else if (FE == TRUE) {
      print_info <- print(paste0("Files already existed in: ", pathD2, "  ---- No files were downloaded."))
    }

    # Make dataframe

    gendat <-
      readr::read_delim(pathD10,
                        delim = "\t",
                        comment = "#",
                        col_names = FALSE)

    # chromosome ID default

    if (is.null(chr_id)) {
      chrID2 <- gendat$X1[[1]]
    } else {
      chrID2 <- chr_id
    }

    GD = gendat$X3[[1]]

    if (GD == "gene") {
      SG <- "gene"
    } else if (GD == "CDS") {
      SG <- "CDS"
    } else if (GD == "cds") {
      SG <- "cds"
    } else {
      SG <- "CDS"
    }

    if (is_gtf == TRUE) {
      SG <- "CDS"
    } else {}

      gendat2 <- dplyr::filter(gendat, gendat$X3 == SG)


      # Ranges

      gendat_range <- dplyr::filter(gendat2, gendat2$X1 == chrID2)
      range_start <- 1
      range_end <- dplyr::last(gendat_range$X5)

    new_gtf <- str_sub(pathD10, end = -5)

    if (is_gtf == TRUE) {
      gff_gtf <- "_new.gtf"
    } else if (is_gtf == FALSE) {
      gff_gtf <- "_new.gff"
    }



    new_gtf2 <- paste0(new_gtf, gff_gtf, sep = "")

    utils::write.table(
      gendat2,
      file = new_gtf2,
      col.names = FALSE,
      row.names = FALSE,
      sep = "\t",
      quote = FALSE
    )

    gtf <- rtracklayer::readGFFAsGRanges(new_gtf2)

    # Plotting annotation track, genomic coordinates and gene region track

    chr <- as.character(unique(seqnames(gtf)))
    gen <- GenomeInfoDb::genome(gtf)

    options(ucscChromosomeNames = FALSE)
    gtrack <- Gviz::GenomeAxisTrack(name = gen_name, range = gtf)
    atrack <-
      Gviz::AnnotationTrack(
        gtf,
        strand = gtf@strand@elementMetadata,
        name = gen_name,
        genome = gen,
        chromosome = chr,
        featureAnnotation = "id"
      )


    # What to show (gene_id or products)

    if (is_gtf == TRUE) {
      gffI <- gtf$gene_id
    } else if (is_gtf == FALSE) {
      gffI <- gtf$ID
    }


    if (show_gene_id == TRUE) {
      sym <- gffI
    } else if (show_gene_id == FALSE) {
      sym <- gtf$product
    }

    # In .gff are no products, so always show gene ID:

    if (is_gtf == TRUE) {
    } else if (is_gtf == FALSE) {
      sym <- gtf$ID
    }

    grtrack <- Gviz::GeneRegionTrack(
      gtf,
      chromosome = chr,
      symbol = sym,
      name = gen_name,
      transcriptAnnotation = "symbol",
      background.panel = "#FFFEDB",
      background.title = "lightblue",
      shape = "arrow"
    )

    # Adding a sequence track.
    # First tidy fasta file to right chromosome names.

    FA <- pathD4
    fastaFile <- Biostrings::readDNAStringSet(FA)
    seq_name = names(fastaFile)
    sequence = paste(fastaFile)
    datair <- data.frame(seq_name, sequence)

    # remove first characters of chr_id
    NewID <- str_split(chrID2,
                       n = 2,
                       pattern = "_",
                       simplify = TRUE)

    # find matching row in fasta with chr_id
    NewID2 <- as.data.frame(NewID)
    whichRow <- which(grepl(NewID2$V2, datair$seq_name))

    # change seqname in fasta to right chr_id
    datair <- datair %>% dplyr::slice(whichRow)
    datair <- dplyr::transmute(datair, sequence, seqnames = chrID2)

    names <- datair$seqnames
    seq <- datair$sequence
    df <- data.frame(names, seq)

    FA2 <- str_sub(FA, end = -5)
    pathF <- paste0(FA2, "_newSeqnames.fa")

    df.fasta = seqRFLP::dataframe2fas(df, file = pathF)

    strack <-
      Gviz::SequenceTrack(sequence = pathF,
                          chromosome = chrID2,
                          name = gen_name)

    names(strack) <- chrID2

    # Plot tracks.

   my_gen_plot <- Gviz::plotTracks(
      list(gtrack, grtrack, strack),
      chromosome = chrID2,
      from = start,
      to = end,
      col = "grey",
      add53 = TRUE,
      add35 = TRUE,
      littleTicks = TRUE,
      cex = 0.8
    )

    #Show chromosomes?
    gendat_group <- gendat %>% dplyr::select(X1)
    gendat_group <- dplyr::group_by(gendat_group, X1)
    gendat_group <- dplyr::summarise(gendat_group)
    new_row <- data.frame(X1 = "NULL")
    gendat_group <- rbind(new_row, gendat_group)

    print_chr <- list(gendat_group$X1)


    return(list(chrom = print_chr,
                plot = my_gen_plot,
                info = print_info,
                RS = range_start,
                RE = range_end))

  }


# install.packages("BiocManager")
# library("BiocManager")
# BiocManager::install(c("IRanges", "GenomicRanges", "Biostrings", "Gviz", "rtracklayer", "ggbio", "seqRFLP"))

 

