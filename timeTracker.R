library(shiny)
library(bslib)

ui<- page_sidebar(
  title = "Time Tracker",
  
  sidebar = sidebar (       
    selectInput ( "chooseMode", "Choose the mode:", 
                  choices = c("work", "sleep")),
    
  sliderInput("hours", "How many hours spent:", min = 0, max = 12, value = 5),
    
    actionButton("calculateButton", "Calculate time")
    
), 

#selectInput("chooseMode", ) creates a dropdown. The ID is "chooseMode", can be accessed as input$chooseMode in server.
#sliderInput("hours", )	creates a slider. ID is "hours". 
#actionButton("calculateButton",)	creates button to trigger calculation. ID is "calculateButton", Used with eventReactive() in server.

card(
  card_header("Daily Activity Time Tracker"),
  
  value_box(
    title = "Hours Spent", 
    value = textOutput ("hoursWorked")
     ), 
  
 
  
  value_box(
    title = "Productivity level",
    value = textOutput("productivity"),
    theme = "blue" 
    
),

  card_footer("Balance is the key.", "Make sure to sleep for atleast 7 to 8", "Without compromising on your work :)")
)
  )


server <- function (input, output) {
  
  calculated_hours <- eventReactive(input$calculateButton, {
    input$hours
  })                            #calculate only when button is clicked 
  
 # output$hoursWorked <- renderText ({
 #   paste( input$hours, "hour(s)")
 #   }) this syntax caused error 
  
  output$hoursWorked <- renderText({
    req(calculated_hours())                  # Wait until there's a value
    paste(calculated_hours(), "hour(s)")
  })

 output$productivity <- renderText ({
   req(calculated_hours())
   h <- calculated_hours()
   level <- if (h == 0){ 
     "Disappointing Behaviour"
   } else if (h >= 0 & h <= 2) { 
       "Bare Minimum"
   } else if (h >= 2 & h <= 4) {
       "Only satisfactory"
   } else if (h >= 4 & h <= 6) {
       " Doing Better "
   } else if (h >= 6 & h <= 8) {
       "Excellent Job"
   } else {
      "ERROR: It is Enough " 
   } 
   
   level 
   
   })
}


shinyApp(ui, server)