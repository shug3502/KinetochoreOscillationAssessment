library(shiny)
library(ggplot2)

#  ui <- fluidPage(
#    textInput("caption", "Caption", ""),
#    verbatimTextOutput("value")
#  )
#  server <- function(input, output) {
#    output$value <- renderText("FRED") #{ input$caption })
#  }
#  shinyApp(ui, server)
# 
# # creating sample files to upload
# 
# ui <- fluidPage(
#   sidebarLayout(
#     sidebarPanel(
#       fileInput(
#         inputId = "files", 
#         label = "Choose CSV File", 
#         multiple = TRUE,
#         accept = c("text/csv",
#                    "text/comma-separated-values,text/plain",
#                    ".csv")
#       )
#     ),
#     mainPanel(
#       tableOutput("contents")
#     )
#   )
# )
# 
# server <- function(input, output) {
#   output$contents <- renderTable({
#     req(input$files)
#     upload = list()
#     
#     for(nr in 1:length(input$files[, 1])){
#       upload[[nr]] <- read.csv(
#         file = input$files[[nr, 'datapath']]
#       )
#     }
#     
#     return(upload)
#   })
# }
# 
# shinyApp(ui, server)
# 
# 
# ########

source('assess.R')

ui <- fluidPage(
  
  titlePanel("Tabsets"),
  
  sidebarLayout(
    
    sidebarPanel(
      textInput("lineage", "Name of cell line", "eg. MC191"),
      verbatimTextOutput("value"),
      fileInput(
        inputId = "files", 
        label = "Choose CSV File containing tracking from KiT", 
        multiple = TRUE,
        accept = c("text/csv",
                   "text/comma-separated-values,text/plain",
                   ".csv")
      )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", plotOutput("plot")), 
        tabPanel("Summary", verbatimTextOutput("summary"))
      )
    )
  )
)

server <- function(input, output) {
  output$plot <- renderPlot({
    req(input$files)
    req(input$lineage)
    upload_paths = c()
    for(nr in 1:length(input$files[, 1])){
      upload_paths[nr] <- input$files[[nr, 'datapath']]
    }
    print(upload_paths)
    assess_and_compare_files(input$lineage,upload_paths)
  })
  output$summary <- renderText("THANK YOU FOR USING THIS APP")
}

shinyApp(ui, server)

