library(shiny)

shinyUI(fluidPage(
  titlePanel(title = "Make Genomic Plots & Get Gene Sequences"),
  sidebarLayout(
    sidebarPanel(("Input: (Hit 'Plot my genome' to create/update your plot.)"),

                 textInput(inputId = "dl_folder", label = "Please enter the path to the download folder."),

                 textInput(inputId = "gen_name", label = "Please enter genome name."),

                 textInput(inputId = "g_url", label = "Please enter the .gtf or .gff assembly link from the NCBI website."),

                 actionButton("button", "Plot my genome"),
                 
                 textInput(inputId = "chr_id", label = "Enter chromosome. After filling in the three boxes above, you can see the chromosomes to choose from in the 'Chromosome List' tab.", ""),
                 
                 radioButtons(inputId = "endR", label = "Select range length",
                              choices = list("0 - 10.000", "0 - 50.000", "0 - 200.000", "100.000 - 1 million", "800.000 - 10 million")),

                 sliderInput(
                   inputId = "ranges", label = "Please select the plot ranges",
                   min = 0, max = 10000, value = c(1, 5000), step = 5),

                 radioButtons(inputId = "show_gene_id", label = "Show gene identification (TRUE) or products (FALSE)?",
                              choices = list("FALSE", "TRUE")),
                 
                 radioButtons(inputId = "prot", label = "Do you want the protein (TRUE) or the nucleotide (FALSE) sequence?",
                              choices = list("FALSE", "TRUE")),
                 
                 textInput(inputId = "getseq", label = "Enter the gene identification of a gene in the plot to get the fasta output.", ""),
                 
                 actionButton("button2", "Get sequence")),



    mainPanel(("Output"),
              tabsetPanel(type = "tab",
                          tabPanel("My plot", plotOutput("myplot")),
                          tabPanel("Chromosome List", chromz <- textOutput("chrom")),
                          tabPanel("Get sequence", textOutput("fasta"))))

)))

