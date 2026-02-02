library(shiny)
library(shinyjs)
library(jaccard)
library(qvalue)
library(tidyverse)

source(file.path("ui", "helpPopup.R"))

shinyUI(pageWithSidebar(
  headerPanel("Jaccard/Tanimoto Test for Similarity between Binary Variables"),
   sidebarPanel(
     
    fileInput('file1',
              div("Upload a Dataset",
                  helpPopup(title=NULL, content="Browse and select a plain text file with variables as rows.", trigger="hover"),
                  downloadLink("sampleDataFile", "Example")
              ),
              accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),

    selectInput("sep", "Dataset Separator:",
                choices = c(Comma=',',
                            Semicolon=';',
                            Tab='\t',
                            Space=' ')),
    #checkboxInput('header', 'Header', TRUE),
    
    ##### HTML("&pi;<sub>0</sub>")
    wellPanel(
      selectInput("method", p("Computational Method"), 
                   choices = c("bootstrap","mca","asymptotic", "exact"),
                   width = "70%"),
      checkboxInput('compute.qvalue', 'Compute q-values?', TRUE)
    ),
  conditionalPanel(
      condition = "input.method == 'bootstrap'",
      wellPanel(
        p(strong("Arguments for the Bootstrap")),
        numericInput("B", "Iterations:", 1000, width = "70%"),
        numericInput("seed", "Seed:", 1234, width = "70%")
      )
    ),
  conditionalPanel(
      condition = "input.method == 'mca'",
      wellPanel(
        p(strong("Arguments for the MCA")),
        selectInput("error.type", "Error type:", 
                   choices = c("average", "upper", "lower"), width = "70%"),
        numericInput("accuracy", "Error bound:", 1e-05, width = "70%")
      )
  ),
    wellPanel(
      p(strong("Download Outputs")),
      #sliderInput("threshold", "P-value threshold:", step = 0.01, value = 0.05, min = 0, max = 1),
      downloadButton('downloadPvalues', 'P-values'),
      downloadButton('downloadStatistics', 'Stats'),
      downloadButton('downloadExpectation', 'Expectation')
    )
  ),

  mainPanel(
    tabsetPanel(id="tabSelected",
      tabPanel("About", h4("Quick Introduction"), uiOutput("about")),
      tabPanel("Output", uiOutput("subTabs")),
      tabPanel("Help", h4("Descriptions of Input Arguments"), uiOutput("help")))

  )
))

