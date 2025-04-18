library(shiny)
library(bslib)

ui <- page_sidebar(
  title = "Simple Calculator",
  
  sidebar = sidebar(
    numericInput("numb1", "Enter first value:", value = 0),
    numericInput("numb2", "Enter second value:", value = 0),
    
    selectInput("operatorxxx", "Select an operation to be Performed:",
                choices = c("Add", "Subtract", "Multiply", "Divide")),
    actionButton("calculate123", "Calculate")
  ),

  textOutput("result")
)

server <- function(input, output) {
  observeEvent(input$calculate123, {                        # observes when "Calculate" button is clicked
    theresult <- switch(input$operatorxxx,                  #switch() chooses operation based on the selected input
                     "Add" = input$numb1 + input$numb2,
                     "Subtract" = input$numb1 - input$numb2,
                     "Multiply" = input$numb1 * input$numb2,
                     "Divide" = {
                       if (input$numb2 == 0) {
                         "Error: Cannot divide by zero"
                       } else {
                         input$numb1 / input$numb2
                       }
                     })
    
    output$result <- renderText({                        #rendertext generates output dynamically when button is clicked
      paste("Result:", theresult)                        #paste() combines label with the computed results
    })
  })
}


shinyApp(ui = ui, server = server)
