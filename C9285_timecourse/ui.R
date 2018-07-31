
library(shiny)

shinyUI(fluidPage(
  
  # Application title
  titlePanel("Visualize C9285 and T65 RNAseq data"),
  
  sidebarLayout(
    sidebarPanel(
      tags$strong("Which genes should be plotted ?"),
      tags$br(),
      tags$em("Do not use more than 50 genes for plotting !"),
      tags$br(),
      tags$textarea(id = "genesToPlotDirect",
                    paste("DWR_OS06G001006",
                          "DWR_OS03G001620",
                          "DWR_OS08G002109",
                          "DWR_OS12G001081",
                          "DWR_OS03G001140",
                          "DWR_OS02G002918",
                          "DWR_OS01G004400",
                          "DWR_OS05G001135", sep = "\n")),
      tags$br(),
      selectInput("facetswitch",
                  label = "What should be in a single panel ?",
                  choices = list("A single transcript" = "transcript",
                                 "One of the four genotype/leaf stages datasets" = "geno_LS")),
      checkboxGroupInput("whichGenoLS", label = "Which genotype/leafstages should be plotted ?",
                         choices = list("C9285 4-leaf stage" = "C9285_4LS",
                                        "C9285 6-leaf stage" = "C9285_6LS",
                                        "T65 6-leaf stage" = "T65_6LS"),
                         selected = c("C9285_6LS","T65_6LS")),
      numericInput("numFacetCol", label = "How many panel next to each other ?", value = "4"),
      selectInput("plotSort",
                  label = "How should the expression plots be sorted ?",
                  choices = list("by input order" = "inOrder",
                                 "alphanumeric order" = "alphaOrder")),
      checkboxGroupInput("whichQTLregion" , label = "Which QTL region should be displayed ? (select at least one)",
                         choices = list("none","qGTIL3","qGLEI3","qGNEI3","qGLEI8",
                                        "qGTIL9","qGLEI9","qGNEI9","qGNEI10","qGTIL12",
                                        "qTIL1","qRIE1","qLEI3","qTIL12","qNEI12",
                                        "qLEI12","qTIL2","qTIL4"),
                         selected = "none",inline = TRUE
      ),
      numericInput("exprPlotHeight", label = "Height of the expression plot in pixel ?", value = "800")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("expression Plot", plotOutput("expressionPlot", height = 800),
                 dataTableOutput("expr_transcript_info")),
        tabPanel("chromsome overview",
                 plotOutput("chromPlot",
                            brush = "transcript_brush"),
                 tags$p("Click and drag on the chromosome plot to select genes."),
                 dataTableOutput("transcript_info")),
        tabPanel("data table", dataTableOutput("expression_table")),
        tabPanel("MAPMAN level 1", dataTableOutput("MM_lev1_out")),
        tabPanel("MAPMAN level 2", dataTableOutput("MM_lev2_out")),
        tabPanel("MAPMAN level 3", dataTableOutput("MM_lev3_out"))
      )
    )
  )
)
)
