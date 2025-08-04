library(shiny)


ui <- fluidPage(
  #App title
  h2("Plotting thermal performace curves of mosquito traits"),
  plotOutput(outputId = "moz_tpc"),
  fluidRow(column(width = 4,
                  checkboxGroupInput(inputId = "trait_options",
                                label = "trait to plot" )))
)
server <- function(input, output, session) {
}


shinyApp(ui, server)