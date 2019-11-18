#' Function to get the sequence of a gene in fasta format
#'
#' @param gen_name    Name of the genome/folder. Has to be the same as in MakeGenPlot
#' @param dl_folder   Path to the the download folder. Has to be the same as in MakeGenPlot
#' @param gen_id      Name of the wanted gene (The gene name of a gene in a plot of MakeGenPlot when
#'                    show_gen_id = TRUE)
#' @param fasta_name  The first line in the fasta file "> gen_id, strand = ... , fasta_name"
#' @param protein     Logical. To get the protein sequence instead of the nucleotide sequence: protein = TRUE
#'
#' @return A fasta file with the sequence of the given gene.
#'
#' @examples
#' GetSequenceNow(gen_name = "Flavobacterium_ir1",
#'                dl_folder = "C:/Users/boazk/Downloads",
#'                gen_id = "B4N84_RS04735",
#'                fasta_name = "Bladibladibla",
#'                protein = TRUE)
#'
#' @export

GetSequenceNow <-
  function(gen_name = NULL,
           dl_folder,
           gen_id,
           fasta_name,
           protein = FALSE) {

    library(tidyverse)
    library(Biostrings)
    library(rtracklayer)
    library(stringr)
    library(plyr)
    library(seqinr)

    pathG <- dl_folder
    pathG2 <- paste0(pathG, "/", gen_name)
    pathG3 <- paste0(pathG2, "/", gen_name, ".fna")
    pathG4 <- paste0(pathG2, "/", gen_name, ".gtf")
    pathG5 <- paste0(pathG2, "/", gen_name, ".gff")
    pathG7 <- paste0(pathG2, "/", gen_name, "_chr.fa")

    FE <- file.exists(pathG4)

    if (FE == TRUE) {
      pathG6 <- pathG4
    } else if (FE == FALSE) {
      pathG6 <- pathG5
    }

    gtfGen <- pathG6

    options(warn = -1)
    GF <- gtfGen
    gendata <-
      readr::read_delim(GF,
                        delim = "\t",
                        comment = "#",
                        col_names = FALSE)

    GID <- gen_id

    # Filter chromosome names, ranges, strand and gene ID.

    G = gendata$X3[[1]]

    if (G == "gene") {
      GOC <- "gene"
    } else if (G == "CDS") {
      GOC <- "CDS"
    } else if (G == "cds") {
      GOC <- "cds"
    } else {
      GOC <- "gene"
    }

    gendata <- filter(gendata, X3 == GOC)

    # Make gene_id names the same as the MakeGenPlot gen_id names.
    # .gff and .gtf need different sparates when transforming

    if (FE == TRUE) {
      sepp <- 8
    } else if (FE == FALSE) {
      sepp <- 3
    }

    # Transform

    gendata <- gendata %>% separate(X9, into = "gene_id", sep = ";")
    gendata <-
      gendata %>% separate(gene_id,
                           into = c("delete", "gene_id"),
                           sep = sepp)
    gendata <- dplyr::transmute(gendata, X1, X4, X5, X7, gene_id)

    # Remove quotes
    del <- colwise(function(gendata)
      str_replace_all(gendata, '\"', ""))
    gendata2 <- del(gendata)

    # Set ranges and strand for fasta.
    row <- dplyr::filter(gendata2, gendata2$gene_id == GID)
    Sname <- row$gene_id[[1]]
    Sstart <- row$X4[[1]]
    Sstart2 <- as.integer(Sstart)
    Send <- row$X5[[1]]
    Send2 <- as.integer(Send)
    Sstrand <- row$X7[[1]]

    # Filter chromosome sequence by specific row

    fnaGen <- pathG3

    FF <- fnaGen
    fastaFile <- Biostrings::readDNAStringSet(FF)

    seq_name = names(fastaFile)
    sequence = paste(fastaFile)
    datair <- data.frame(seq_name, sequence)
    sequN <- row %>% separate(X1, into = c("extra", "seqname"), sep = 3)
    whichRow <- which(grepl(sequN$seqname[[1]], datair$seq_name))
    datair <- datair %>% dplyr::slice(whichRow)

    # Make fasta from data frame
    df.fasta = seqRFLP::dataframe2fas(datair, file = pathG7)

    fastaFile3 <- Biostrings::readDNAStringSet(pathG7)
    fastaFile2 <- as.character(fastaFile3)

    FastaAA <- as.SeqFastaAA(fastaFile2, name = GID, Annot = Sname)

    # Filter gene sequence with getFrag with the aquired ranges
    seqfrag <-
      seqinr::getFrag(FastaAA[[1]],
                      begin = Sstart2,
                      end = Send2,
                      name = Sname)
    seqfrag2 <- dplyr::as_data_frame(seqfrag)

    # Protein or nucleotide
    z = protein

    if (z == TRUE) {
      P <- "_P"
    } else if (z == FALSE) {
      P <- "_N"
    }

    # Write sequence in fasta file and make it nice and clean
    FF2 <- str_sub(FF, end = -5)
    pathS <- paste0(FF2, "_", gen_id, P, ".fna")
    utils::write.table(
      seqfrag2,
      file = pathS,
      sep = "",
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE
    )
    seqfrag4 <- utils::read.table(pathS, sep = "\t", quote = "")
    seqfrag5 <- as.alignment(seq = seqfrag4$V1)
    Fname <- paste0(GID, ", strand = ", Sstrand, " , ", fasta_name)
    seqinr::write.fasta(seqfrag5, names = Fname, file.out = pathS)

    seqfrag6 <- utils::read.table(pathS, sep = "\t", quote = "")
    seqfrag6 <- seqfrag6[!(seqfrag6$V1 == ">NA"), ]

    utils::write.table(
      seqfrag6,
      file = pathS,
      sep = "",
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE
    )

    # Strand reverse/foward:

    S = Sstrand

    if (S == "-") {
      Nstrand <- "R"
    } else if (S == "+") {
      Nstrand <- "F"
    }

    # Protein or nucleotide?
    z = protein

    if (z == TRUE) {
      Fas <- "Protein"
      seqfrag7 <- seqinr::translate(seqfrag, sens = Nstrand)
      seqfrag8 <- dplyr::as_data_frame(seqfrag7)
      utils::write.table(
        seqfrag8,
        file = pathS,
        sep = "",
        col.names = FALSE,
        row.names = FALSE,
        quote = FALSE
      )
      seqfrag9 <- utils::read.table(pathS, sep = "\t", quote = "")
      seqfrag10 <- as.alignment(seq = seqfrag9$V1)
      Fname <- paste0(GID, ", strand = ", Sstrand, " , ", fasta_name)
      seqinr::write.fasta(seqfrag10, names = Fname, file.out = pathS)
      seqfrag11 <- utils::read.table(pathS, sep = "\t", quote = "")
      seqfrag11 <- seqfrag11[!(seqfrag11$V1 == ">NA"), ]
      utils::write.table(
        seqfrag11,
        file = pathS,
        sep = "",
        col.names = FALSE,
        row.names = FALSE,
        quote = FALSE
      )

    } else if (z == FALSE) {
      Fas <- "Nucleotide"
      utils::write.table(
        seqfrag6,
        file = pathS,
        sep = "",
        col.names = FALSE,
        row.names = FALSE,
        quote = FALSE
      )
    }

    seqfrag99 <- utils::read.table(pathS, sep = "\t", quote = "")

    # Print information

      outp <- paste0(
        "Ranges of gene ",
        gen_id,
        ": ",
        Sstart2,
        " - ",
        Send2,
        ". ",
        "Strand: ",
        Sstrand,
        ". ",
        Fas,
        " fasta sequence saved to path: ",
        pathS)


    outs <- as.character(seqfrag99$V1)

    text <- print(outp)

    return(text)


  }
