#user interface 
#card and value box
#Building first shiny app from scratch in detail 
#laying out the user interface 

library(shiny)           #loading libraries 
library(bslib)           #package in R that helps us style and design our Shiny apps


# Defining UI
ui <- page_sidebar(      #page_sidebar() to create a page with a sidebar
title = "Creating Shiny app in Detail",
sidebar = sidebar("Input contents"),  #by default sidebar appears at left, we can change allingment by giving the argument: position="right"

value_box(                      #Value boxes, useful UI component. Used to highlight important values in our app.
  title = "Value box",
  value = 100,
  theme = "teal"
),

card (                                #using card() function like: card(Contents) displays it in a separate rectangular box
  card_header( "Contents"),           #Card elements 
  card_body("Main contents"),
  card_footer("Extra contents")
)

  
)

# Defining server logic 
server <- function(input, output) {
  
}

# Running the app
shinyApp(ui = ui, server = server)
