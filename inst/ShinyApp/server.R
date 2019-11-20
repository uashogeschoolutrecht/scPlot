
source("helpers.R")
library(shiny)
library(grDevices)
library(BiocGenerics)
library(utils)

shinyServer(function(input, output, session){

    observe({
    range <- input$endR

    if (range == "0 - 10.000"){
      updateSliderInput(session, "ranges", label = "Please select the plot ranges",
                        min = 0, max = 10000, value = c(1, 1000), step = 1)
    }
    if (range == "0 - 50.000"){
      updateSliderInput(session, "ranges", label = "Please select the plot ranges",
                        min = 0, max = 50000, value = c(1, 5000), step = 1)
    }
    if (range == "0 - 200.000"){
      updateSliderInput(session, "ranges", label = "Please select the plot ranges",
                        min = 0, max = 200000, value = c(1, 20000), step = 1)
    }
    if (range == "100.000 - 1 million"){
      updateSliderInput(session, "ranges", label = "Please select the plot ranges",
                        min = 100000, max = 1000000, value = c(100001, 190000), step = 1)
    }
    if (range == "800.000 - 10 million"){
      updateSliderInput(session, "ranges", label = "Please select the plot ranges",
                        min = 800000, max = 10000000, value = c(800001, 1720000), step = 1)
    }
    
  })

  plotnow <- eventReactive(input$button, {MakeGenPlot(g_url = input$g_url,
                                                     dl_folder = ".",
                                                     chr_id = input$chr_id,
                                                     gen_name = input$gen_name,
                                                     start = input$ranges[1],
                                                     end = input$ranges[2],
                                                     show_gene_id = input$show_gene_id)})
  

  plotnu <- reactive({ MakeGenPlot(g_url = input$g_url,
                                   dl_folder = ".",
                                   chr_id = input$chr_id,
                                   gen_name = input$gen_name,
                                   start = input$ranges[1],
                                   end = input$ranges[2],
                                   show_gene_id = input$show_gene_id)
   
  })

  
  chromnow <- eventReactive(input$button, {GetChroms(g_url = input$g_url,
                                                     dl_folder = ".",
                                                     gen_name = input$gen_name)})
  
  seqnow <- eventReactive(input$button2, {GetSequenceNow(gen_name = input$gen_name,
                                                        dl_folder = ".",
                                                        gen_id = input$getseq,
                                                        fasta_name = input$gen_name,
                                                        protein = input$prot)})
   
  
   output$myplot <- renderPlot({plotnow()})

   output$chrom <- renderText({chromnow()})

   output$fasta <- renderText({seqnow()})
   
   
   #----------------Download button for plot-------------------------------
   
   output$Download <- downloadHandler(
     filename = function(){
         base::paste0(input$gen_name, "_", input$chr_id, "_", input$ranges[1], "-", input$ranges[2], ".png")},
                                    content = function(file){
                                        png(file)
                                        plotnu()
                                        dev.off()
                                    })
   
   #----------------Download button for fasta-------------------------------
   
   output$Download2 <- downloadHandler(
     filename = function(){
     base::paste(input$gen_name, "_", input$getseq, ".fna", sep = "")
       },
     content = function(con){
       file.copy(from = file.path(input$gen_name, paste0( 
         input$gen_name, "_", 
         input$getseq, ".fna"))
         , to = con)
     })
   
})

