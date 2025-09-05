library(shiny)
library(bslib)
library(bsicons)
library(shinyWidgets)
library(shinythemes)

# Define UI ----
ui <- page_sidebar(
  #theme = shinythemes::shinytheme("flatly"),
  theme = bs_theme(preset = "shiny"),
  title = "These models are run using the rTPC package in R.",
  sidebar = sidebar(
    helpText("Select what models to fit and what data to include"),
    title = "Model Controls",
    
    ##choose trait -----
     layout_columns(
      radioButtons(inputId = "trait",
                   "Trait:",
                   levels(dat$trait.name)
      ),
     # tooltip(
     #    bsicons::bs_icon("info-circle", title = "About tooltips"),
     #    "Text shown in the tooltip.",
     #    placement = "bottom"
     #  ),
      tooltip( bsicons::bs_icon("info-circle", title = "About tooltips"),
        "Text shown in the tooltip2.",
        placement = "left"
      ),
      col_widths = c(6,6,6)
    ),
   
    
    ##choose species -----
    pickerInput(
      inputId = "species", 
      label = "Species to include:", 
      choices = levels(dat$sp_name), 
      selected = levels(dat$sp_name)[1:10], #this selects 1st 10 initially
      options = pickerOptions(
        actionsBox = TRUE,  #this adds the select all / deselect all
        size = 5,
        selectedTextFormat = "count > 3",
        liveSearch = T
      ), 
      multiple = TRUE
    ),
    
    ##choose model -----
    pickerInput(
      inputId = "model", 
      label = "Models to fit:", 
      choices = model_options, 
      selected = model_options[1],
      options = pickerOptions(
        actionsBox = T, 
        size = 5,
        selectedTextFormat = "count > 3", 
        liveSearch = T
      ), 
      multiple = TRUE
    ),
    
    ##best model tickbox -----
    checkboxInput(inputId = "highlight_best", "Highlight TPC with lowest AIC?", FALSE)
  ),

  ##output cards -----
  layout_columns(
  navset_card_tab(
    nav_panel("Plot", plotOutput('plot1')), 
    nav_panel("Table", tableOutput("table1"),
              # Button
              downloadButton("downloadData", "Download")
              )
    ),
  card(textOutput("selected_species"),
  textOutput("selected_model")
  ),
  col_widths = c(9,3))
)




# Define server logic ----
server <- function(input, output) { 
  

  output$plot1 <- renderPlot({      
    validate(
        need(input$species != "", "Please select a species"),
        need(input$model != "", "Please select a model")
      )
    ifelse(input$highlight_best == FALSE,

      shiny_TPC(input$trait, input$species, input$model)[1],
      shiny_TPC(input$trait, input$species, input$model)[3]
    )
    })


output$table1 <- renderTable({
  validate(
    need(input$species != "", "Please select a species"),
    need(input$model != "", "Please select a model")
  )
  shiny_TPC(input$trait, input$species, input$model)[2]
})


#this was the only way i got it to work, but feels wrong- calling shinytpc again rather than table1
output$downloadData <- downloadHandler(
  filename = function() {
    paste("data-", Sys.Date(), ".csv", sep="")
  },
  content = function(file) {
    write.csv(shiny_TPC(input$trait, input$species, input$model)[2], file)
  }
)

output$selected_species <- renderText({
  validate(
    need(input$species != "", "Please select a species")
  )
  paste("The species you have selected are: ", paste(input$species, collapse = ", "))
})

output$selected_model <- renderText({
  validate(
    need(input$model != "", "Please select a model")
  )
  paste("The models you have selected are: ", paste(input$model, collapse = ", "))
  
})

}

# Run the app ----
shinyApp(ui = ui, server = server)