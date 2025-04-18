library(shiny)                              #for creating app
library(bslib)                              #for layout

# UI ----
ui <- page_sidebar(
  title = "Meerub's Coffee Tracker"  ,
  
  sidebar = sidebar(                                   #creates dropdown menu
    selectInput("coffeeType", "Choose your brew:",     #"coffeeType" is input ID, refers to the user's selection in the server.
                choices = c("Espresso", "Latte", "Cappuccino", "decaf", "iced latte")),
    
    sliderInput("cups", "How many cups today?", min = 0, max = 10, value = 4),
    "DONT EXCEED 7!!!" , 
    
    actionButton("brewButton", "Brew Coffee pleaseee", class = "btn-success") , 
    "Storyline : It's 3am and instead of sleeping, I built an app to track how much coffee I drank to stay awake building it. "
  ),
  
  card(
    card_header("Daily Coffee Dashboard"),
    
    value_box(
      title = "Cups Consumed",
      value = textOutput("cupCount"), #this is where server writes the number of cups
      theme = "brown",
      "I may look normal, but trust me I am high on caffeine"
    ),
    
    value_box(
      title = "Productivity Level",
      value = textOutput("productivity"),
      theme = "orange"
    ),
    
    card_footer("Stay caffeinated = Stay productive")
  )
)


server <- function(input, output, session) {
  

  
  output$cupCount <- renderText({
    paste(input$cups, "cups
          ")
  })
  
  output$productivity <- renderText({
    level <- if (input$cups == 0) {
      "Sleeping "
    } else if (input$cups <= 2) {
      "Warming up "
    } else if (input$cups <= 5) {
      "Focused "
    } else if (input$cups <= 8) {
      "High on Caffiene"
    } else {
      "ERROR: Too Much Coffee! 🫣"
    }
    level                          #result was stored in level variable and returned
  })
}

shinyApp(ui, server)
