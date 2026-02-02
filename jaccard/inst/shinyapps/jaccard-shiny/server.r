library(shiny)
library(shinyjs)
library(jaccard)
library(qvalue)
library(tidyverse)
load("table.Rdata")

shinyServer(function(input, output) {
  output$subTabs <- renderUI({
    tabsetPanel(id = "subTabPanel1",
                tabPanel("Summary", verbatimTextOutput("summary")),
                tabPanel("Scatter Plot",plotOutput("scatterPlot", height="720px", width="720px"), downloadButton('downloadScatter', 'Download Scatter Plot')),
                tabPanel("P-value Histogram", plotOutput("pvalueHist", height="720px", width="720px"), downloadButton('downloadPvalueHist', 'Download P-value Histogram')),
                tabPanel("P-value Heatmap", plotOutput("pvalueHeat", height="720px", width="720px"), downloadButton('downloadPvalueHeat', 'Download P-value Heatmap')),
                tabPanel("Q-value Plots", plotOutput("qvalueHist", height="720px", width="720px"), downloadButton('downloadQvalueHist', 'Download Q-value Plots'))
    )
  })
  # Return the requested dataset

  # input$file1 will be NULL initially. After the user selects and uploads a
  # file, it will be a data frame with 'name', 'size', 'type', and 'datapath'
  # columns. The 'datapath' column will contain the local filenames where the
  # data can be found.
  methodInput <- reactive({
    switch(input$method,
           "exact" = "exact",
           "asymptotic" = "asymptotic",
           "mca" = "mca",
           "bootstrap" = "bootstrap")
  })
  datasetInput <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) return()
    dataset <- read.csv(inFile$datapath, header=TRUE, sep=input$sep, row.names = 1)
    # generate an rnorm distribution and plot it
    #dataset <-unlist(dataset,use.names=FALSE)
    dataset <- as.matrix(dataset)

    if(methodInput() == "exact") {
      output <- jaccard.test.pairwise(dataset, method=methodInput(), verbose=TRUE, compute.qvalue=input$compute.qvalue)
    } else if(methodInput() == "asymptotic") {
      output <- jaccard.test.pairwise(dataset, method=methodInput(), verbose=TRUE, compute.qvalue=input$compute.qvalue)
    } else if(methodInput() == "mca") {
      output <- jaccard.test.pairwise(dataset,
                                       method=methodInput(),
                                       verbose=FALSE,
                                       compute.qvalue=input$compute.qvalue,
                                       error.type=input$error.type,
                                       accuracy=input$accuracy)
    } else if(methodInput() == "bootstrap") {
      output <- jaccard.test.pairwise(dataset,
                                       method=methodInput(),
                                       verbose=FALSE,
                                       compute.qvalue=input$compute.qvalue,
                                       B=input$B,
                                       seed=input$seed)
    }

  })
  # Generate a summary of the dataset
  output$summary <- renderPrint({
    out <- datasetInput()
    if (is.null(out)) return()
    cat("Computation Completed! \n \n")
    cat(paste("The uploaded dataset has", nrow(out$statistics) ,"variables (rows). \n"))
    cat(paste("All pair-wise comparisons of", nrow(out$statistics) ,"variables are evaluated. \n \n"))
    cat("The resulting output is consisted of \n")
    print(summary(out))
    cat("\n Summary of pair-wise p-values: \n")
    print(summary(as.vector(out$pvalues),na.rm=TRUE))
    cat("\n Summary of pair-wise statistics: \n")
    print(summary(as.vector(out$statistics),na.rm=TRUE))
  })
  # function/render/download a scatterplot
  scatterPlot <- function(){
    out <- datasetInput()
    if (is.null(out)) return()
    plot(out$pvalues, out$statistics, xlab="P-values", ylab="Statistics", pch=20)
    abline(v=input$threshold,col="red",lty="dashed")
  }
  output$scatterPlot <- renderPlot({ scatterPlot() })
  output$downloadScatter <- downloadHandler(
    filename<- paste0('scatter', Sys.Date(),'.pdf'),
    content <- function(file) {
      pdf(file, width=10, height=10)
      scatterPlot()
      dev.off()
    },
    contentType = 'image/pdf'
  )

  # function/render/download a pvalueHist
  pvalueHist <- function(){
    out <- datasetInput()
    if (is.null(out)) return()
    hist(out$pvalues,col="black", xlim=c(0,1), xlab="P-values")
    abline(v=input$threshold,col="red",lty="dashed")
  }
  output$pvalueHist <- renderPlot({ pvalueHist() })
  output$downloadPvalueHist <- downloadHandler(
    filename<- paste0('pvalueHist', Sys.Date(),'.pdf'),
    content <- function(file) {
      pdf(file, width=10, height=10)
      pvalueHist()
      dev.off()
    },
    contentType = 'image/pdf'
  )

  # function/render/download a qvalueHist
  qvalueHist <- function(){
    out <- datasetInput()
    if (is.null(out)) return()
    plot(out$qvalues)
  }
  output$qvalueHist <- renderPlot({ qvalueHist() })
  output$downloadQvalueHist <- downloadHandler(
    filename<- paste0('qvalueHist', Sys.Date(),'.pdf'),
    content <- function(file) {
      pdf(file, width=10, height=10)
      qvalueHist()
      dev.off()
    },
    contentType = 'image/pdf'
  )
  
  # function/render/download a pvalueHist
  pvalueHeat <- function(){
    out <- datasetInput()
    if (is.null(out)) return()
    m <- nrow(out$pvalues)
    pvalues.tidy <- out$pvalues %>%
      as_tibble() %>%
      rownames_to_column('Var1') %>%
      gather(Var2, value, -Var1) %>%
      drop_na() %>%
      mutate(
        Var1 = factor(Var1, levels=1:m),
        Var2 = factor(gsub("V", "", Var2), levels=1:m)
      )
    
    ggplot(pvalues.tidy, aes(Var1, Var2))+
      geom_tile(aes(fill = value))+
      scale_fill_gradient2(low = "blue", high = "darkorange", mid = "darkgreen", space = "Lab", 
                           name="P-values") +
      theme_minimal() +  coord_fixed() + xlab("") + ylab("") +
      theme(strip.text.x = element_blank(),
            strip.background = element_rect(colour="white", fill="white"))
    
  }
  output$pvalueHeat <- renderPlot({ pvalueHeat() })
  output$downloadPvalueHeat <- downloadHandler(
    filename<- paste0('pvalueHeat', Sys.Date(),'.pdf'),
    content <- function(file) {
      ggsave(file, plot=pvalueHeat(), device="pdf", width=10, height=10)
    },
    contentType = 'image/pdf'
  )
  
  output$downloadPvalues <- downloadHandler(
    filename = paste0("pvalues_", methodInput() , "_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(datasetInput()$pvalues, file)
    },
    contentType = "text/csv")
  output$downloadStatistics <- downloadHandler(
    filename = paste0("statistics_", methodInput() , "_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(datasetInput()$statistics, file)
    },
    contentType = "text/csv")  
  output$downloadExpectation <- downloadHandler(
    filename = paste0("expectation_", methodInput() , "_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(datasetInput()$expectation, file)
    },
    contentType = "text/csv")  
  output$help <- renderTable({
    if (is.null(table)) {return()}
    print(table)
  }, 'include.rownames' = FALSE
  , 'include.colnames' = TRUE
  , 'sanitize.text.function' = function(x){x}
  )
  output$about <- renderUI({
    wellPanel(
      p("This is an interactive web interface for the ", a("jaccard R package", href="https://github.com/ncchung/jaccard"), "accompanying a manuscript titled `Jaccard/Tanimoto binary similarity test, estimation methods, and implementations to probabilistically evaluate species co-occurrences'. Generally, competition between different operational taxonomic units (OTUs) can be evaluated by a Jaccard/Tanimoto coefficient applied to their absence/presence vectors across multiple bioregions."),
      p(a("Jaccard similarity coefficient or index", href="https://en.wikipedia.org/wiki/Jaccard_index"), "is one of the most popular similarity measures for binary biological data; namely, the ratio of their intersection to their union. It's alternatively called Tanimoto coefficient and Intersection over Union. This app performs pair-wise similarity tests among binary variables (i.e., rows), outputting p-values and q-values."),
      #p("Competition between different operational taxonomic units (OTUs) can be evaluated by a Jaccard/Tanimoto coefficient applied to their absence/presence vectors across multiple bioregions."),
      p("Although this app and example target the absense/presence dataset in biogeography, this test is readily applicable to biochemical fingerprints, genomic intervals, and other binary datasets."),
      p(""),
      p("To obtain p-values for similarity:"),
      p("1. Prepare a dataset in a plain text file with one binary variable per row. This app will evaluate pair-wise similarity among variables. The first column and row are assumed to be the row names and column names (headers), respectively."),
      p("2. Upload this dataset using the 'Choose File' button located on the left side panel."),
      p("3. Choose appropriate options used to estimate p-values. See the 'Help' tab for explanations of these options."),
      p("4. View and download selected visualizations of the results in the 'Output' tab."),
      p("5. Save the p-values and other outputs under the 'Download Outputs' panel on the left side panel."),
      p(" "),
      p("For finer controls and further analyses, please use the R package available in the ", a("Github Repository.", href="https://github.com/ncchung/jaccard")),
      p("Created by Neo Christopher Chung -- ", a("Contact Details", href="http://ncc.name"))
    )
  })

  # download sample data file
  output$sampleDataFile <- downloadHandler(
    filename = function() { "ExampleDataset_BirdSpecies.csv" },
    content = function(file) {
      file.copy(from = "BirdSpecies.csv", to = file, overwrite = TRUE)
    }
  )
})

