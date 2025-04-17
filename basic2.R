#test 1 

library(shiny)
library(bslib)

ui <- page_sidebar(
  title = "My Shiny App",
  sidebar = sidebar ("Shiny is available on CRAN, so you can install it in the usual way from your R console:",
  code('install.packages("shiny")'),
),

card(
  card_header ("Introducing Shiny"),
  "Shiny is a package from Posit that makes it incredibly easy to build interactive web applications with R.",
  card_footer("Shiny is a product of postit")
  
)
) 
  
server <- function (input, output)  {

}
shinyApp( ui = ui , server = server)  
