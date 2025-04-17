#Widgets 
#its a web element that users interact with. It provides a way for users to send messages to the Shiny app.
#Shiny widgets collect a value from the user.When we change the widget, the value also changes. 
#shiny comes with multiple pre-built widgets
#we can add widgets functions into our web page
#they require several arguments, first two for each widget are:
#1. Name for the widget: User cant see this name, but we can use it to access the widget’s value.The name should be a character string.
#2.A label:label appears with the widget in our app. It should be a character string, but it can be an empty string "".
#The other arguments depend widget to widget depending on there job. they may include(initial values, ranges and increments.)

library(shiny)

# Defining UI
ui <- page_fluid(
  titlePanel("Trying out widgets"),
  layout_columns(
    col_width = 4,
    card(
      card_header("Buttons"),
      "Select your gender before starting the survey",
      actionButton("action", "male"),
      actionButton("action", "female"),
      submitButton("Done")
    ),
    card(
      card_header("Single checkbox"),
      "Are you hungry?",
      checkboxInput("checkbox", "yes", value = TRUE) #value can also be set to false
    ),
    card(
      card_header("Multiple Checkbox group"),
      checkboxGroupInput(
        "checkGroup",
        "Choose from the options given below",
        choices = list("Apple" = 1, "Grapes" = 2, "Pasta" = 3),
        selected = 3
      )
    ),
    card(
      card_header("Date Input"),
      dateInput("date", "Date of Birth", value = "2004-03-18")
    ),
    card(
      card_header("Date range input"),
      dateRangeInput("dates", "Check dates for availability")
    ),
    card(
      card_header("Insert File"),
      fileInput("file", label = NULL)
    ),
    card(
      card_header("Help text"),
      helpText(
        "Note: Not a true widget,",
        "provides an easy way to add text to",
        "accompany other widgets."
      )
    ),
    card(
      card_header("Numeric input"),
      numericInput("num", "My age:", value = 21)
    ),
    card(
      card_header("Radio buttons"),
      radioButtons(
        "radio",
        "Is programming fun",
        choices = list("yes" = 1, "no" = 2, "obv" = 3),
        selected = 3
      )
    ),
    card(
      card_header("Drop down option"),
      selectInput(
        "select",
        "Select country",
        choices = list("italy" = 1, "turkey" = 2, "egypt" = 3),
        selected = 2
      )
    ),
    card(
      card_header("Sliders"),
      sliderInput(
        "slider1",
        "Set value",
        min = 0,
        max = 50,
        value = 25
      ),
      sliderInput(
        "slider2",
        "Set value range",
        min = 0,
        max = 50,
        value = c(15, 35)
      )
    ),
    card(
      card_header("Text input"),
      textInput("text", label = NULL, value = "Enter address...")
    ),
    
    card(
      card_header("last question"),
      checkboxGroupInput( 
        "checkGroup" , 
        "Did you answer all the questions?",
        choices = list ( "yes" = 1, "no" = 2 ),
        selected = 2
        
      )
    ),
    
  )
)

server <- function(input, output) {
}

shinyApp(ui = ui, server = server)
