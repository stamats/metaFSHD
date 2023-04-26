load("resData_metaFSHD.RData")
source("forestplotFSHD.R")
library(grid)
library(forestploter)
library(metafor)
library(shiny)
library(shinythemes)
library(ggplot2)

ui <- fluidPage(
  theme = shinytheme("darkly"),

  # Application title
  titlePanel("Random-Effects Meta-Analysis of FSHD Omics-Data"),
  
  sidebarLayout(
    sidebarPanel(
      h3('Provide ID', align = "left"),
      radioButtons("ID", 
                   label = ("Select type of ID"), 
                   choices = list("Gene Name" = 1, 
                                  "ENSEMBL-ID" = 2,
                                  "ENTREZ-ID" = 3), 
                   selected = NULL),
      conditionalPanel(
        condition = "input.ID == 1",
        selectizeInput(inputId = "GeneName",
                       label = "Select gene name",
                       choices = NULL
        )
      ),
      conditionalPanel(
        condition = "input.ID == 2",
        selectizeInput(inputId = "ENSEMBL",
                       label = "Select ENSEMBL ID",
                       choices = NULL
        )
      ),
      conditionalPanel(
        condition = "input.ID == 3",
        selectizeInput(inputId = "ENTREZ",
                       label = "Select ENTREZ ID",
                       choices = NULL
        )
      ),
      downloadButton('savePlot', 'Save Forestplot (pdf)')
    ),
    mainPanel(
      plotOutput("plot", width = "100%", height = "720px")
    )
  )
)

server <- function(input, output, session) {
  ensembl_list <- sort(resData$ENSEMBL)
  gene_name_list <- resData$gene_name[resData$gene_name != ""]
  gene_name_list <- sort(unique(gene_name_list))
  entrez_list <- resData$entrezid[resData$entrezid != ""]
  entrez_list <- unique(entrez_list)
  ## first entry for sorting
  entrez_list1 <- sapply(strsplit(entrez_list, " \\| "), "[", 1)
  entrez_list <- entrez_list[order(as.numeric(entrez_list1))]
  updateSelectizeInput(session, 'GeneName', choices = gene_name_list, selected = "7SK", server = TRUE)
  updateSelectizeInput(session, 'ENSEMBL', choices = ensembl_list, selected = "ENSG00000000003", server = TRUE)
  updateSelectizeInput(session, 'ENTREZ', choices = entrez_list, selected = "1", server = TRUE)
  
  
  output$plot <- renderPlot({
    if(input$ID == 1)
      gg <- forestplotFSHD(resData, gene.name = input$GeneName, called.in.app = TRUE)
    if(input$ID == 2)
      gg <- forestplotFSHD(resData, ensemblid = input$ENSEMBL, called.in.app = TRUE)
    if(input$ID == 3)
      gg <- forestplotFSHD(resData, entrezid = input$ENTREZ, called.in.app = TRUE)
    p.wh <- get_wh(plot = gg, unit = "in")
    gg.plot <<- gg
    gg
  })
  output$savePlot <- downloadHandler(
    filename = function() { 
      if(input$ID == 1)
        filename <- paste("metaFSHD_", input$GeneName, ".pdf", sep = "")
      if(input$ID == 2)
        filename <- paste("metaFSHD_", input$ENSEMBL, ".pdf", sep = "")
      if(input$ID == 3)
        filename <- paste("metaFSHD_", input$ENTREZ, ".pdf", sep = "")
      filename
    },
    content = function(file) {
      ggsave(file, plot = gg.plot, width = 12, height = 8)
    }
  )
}

shinyApp(ui, server)