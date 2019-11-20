# scPlot

# Overview
scPlot is a tool for plotting genes and researching them in an interactive way. You can also easily get protein or nucleotide sequences of any gene on the plot. Just type in the gene identification and the tool will automaticly provide an fasta file with the sequence. 

# Functions

- MakeGenPlot() Plots the genome and shows products or gene identifications of the genes.
- GetSequenceNow() Gets the sequence of a gene and makes an fasta file for it.


- GetChroms() Gets the different chromosomes of an genome. (meant only for the shiny app)

# Input

As input you need 3 things at least:

- The download folder (dl_folder) to store the files in. The app/function will automaticly  make a new folder in the download folder with the "gen_name" as folder name. The generated fasta files will also go to this folder.

- The genome name (gen_name) is the name of the folder and files. With this name the functions MakeGenPlot and GetSequenceNow are linked. Always first make a plot with MakeGenPlot and then you can use GetSequenceNow with the same gen_name as MakeGenPlot to get a sequence. 

- Preferably a "_genomic.gtf.gz" url from an RefSeq or GenBank assembly from the NCBI website. "_genomic.gff.gz" will also work, but won't show any products in the plot. With the ".gtf" of ".gff" url it will automaticly also download and read the "_genomic.fna.gz" file for the sequences. 


