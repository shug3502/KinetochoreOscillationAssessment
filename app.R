library(shiny)
library(ggplot2)


source('assess.R')

ui <- fluidPage(
  
  titlePanel("Knetochore Oscsillation Assessment App"),
  
  sidebarLayout(
    
    sidebarPanel(
      textInput("lineage", "Input name of cell line", "eg. MC191"),
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
  output$summary <- renderText("THANK YOU FOR USING THIS APP! (Created by Jonathan U Harrison 2022-04-11)")
}

shinyApp(ui, server)

