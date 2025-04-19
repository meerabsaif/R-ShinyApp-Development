library(shiny)
library(bslib)

ui <- page_sidebar(
  title = "Workout Tracker",
  sidebar = sidebar (
    selectInput("exertype", "Choose Type of Activity:",
                choices = c("cardio", "strength training", "pilates", "yoga"),
                ),
    sliderInput("timespent1", "Time active in hours:" , min = 0, max = 5, value = 2),
    
    actionButton("calcbttn", "Calculate")
),

card(
  card_header("Workout tracker"),
  
  value_box( 
    title = "Time Active:" , 
    value = textOutput ("timespent2")
    
    ),
  
  value_box(
    title = "Total calories burnt:" , 
    value = textOutput ("calories") , 
    theme = "yellow"
  ),
  
  card_footer(" A healthy body= A healthy brain")
)
 )

server <- function (input, output) {

   calculated_time <- eventReactive (input$calcbttn, {
        input$timespent1
})
  
 output$timespent2 <- renderText({
   req(calculated_time())
   paste (calculated_time(), "hours")
   })
  
 output$calories <- renderText ({
   req(calculated_time())
   t <- calculated_time()
   level <- if (t == 0) {
     "zero calories"
   } else if ( t >= 0 & t<= 1 ){
      "100 calories burnt"
   } else if ( t >= 1 & t<= 2 ) {
      "200 calories"
   } else if ( t >= 2 & t<= 3 ) {
     "300 calories"
   } else if ( t >= 3 & t<= 4 ) {
     "400 calories"
   } else {
     "limit exceeded, not recommended"
   } 
   
   level
     
   })
 
 }
  
shinyApp(ui, server)  