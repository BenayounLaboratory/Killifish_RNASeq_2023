library(shiny)

# Loading Data
brain.age <- read.csv("data/shiny_brain_age.csv", row.names = NULL)
brain.sex <- read.csv("data/shiny_brain_sex.csv", row.names = NULL)
heart.age <- read.csv("data/shiny_heart_age.csv", row.names = NULL)
heart.sex <- read.csv("data/shiny_heart_sex.csv", row.names = NULL)
muscle.age <- read.csv("data/shiny_muscle_age.csv", row.names = NULL)
muscle.sex <- read.csv("data/shiny_muscle_sex.csv", row.names = NULL)
spleen.age <- read.csv("data/shiny_spleen_age.csv", row.names = NULL)
spleen.sex <- read.csv("data/shiny_spleen_sex.csv",row.names = NULL)

# Define UI for application
ui <- fluidPage(
  
  # Application title
  titlePanel("Interactive N. furzeri Transcriptomic Database"),
  sidebarLayout(
    sidebarPanel(
      selectInput("var", 
                  label = "Which Tissue is Of Interest?",
                  choices = list("Brain with Age as Variable" = 1,
                                 "Brain with Sex as Variable" = 2,
                                 "Heart with Age as Variable" = 3,
                                 "Heart with Sex as Variable" = 4,
                                 "Muscle with Age as Variable" = 5, 
                                 "Muscle with Sex as Variable" = 6,
                                 "Spleen with Age as Variable" = 7,
                                 "Spleen with Sex as Variable" = 8),
                  selected = 1),
      
      textInput("p",
                   h3("Export All Genes With Adjusted p-value Cutoff Of")),
      
      textInput("log",
                   h3("Export All Genes With Adjusted Log Cutoff Of")),
      
      textInput("human.s",
                   h3("Find Human Gene Symbol Of")),
      
      textInput("org.s",
                h3("Find Original Gene Nomenclature Of")),
      
      downloadButton('download', 'Download the Displayed Table')
    ),
    
    mainPanel(
      dataTableOutput("db")
    )
  )
  
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  data <- reactive({
    if (input$var == 1){
      return(brain.age)
    } else if (input$var == 2){
      return(brain.sex)
    } else if (input$var == 3) {
      return(heart.age)
    } else if (input$var == 4){
      return(heart.sex)
    } else if (input$var == 5){
      return(muscle.age)
    } else if (input$var == 6){
      return(muscle.sex)
    } else if (input$var == 7){
      return(spleen.age)
    } else {
      return(spleen.sex)
    }
  })
  
  data.p <- reactive({
    if(input$p == ""){
      return(data())
    } else{
      validate(need(!is.na(as.numeric(input$p)), "Enter a number for the p value!"))
      temp <- subset(data(), padj <= as.numeric(input$p))
      return(temp)
    }
  })
  
  data.p.log <- reactive({
    if(input$log == ""){
      return(data.p())
    } else {
      validate(need(!is.na(as.numeric(input$log)), "Enter a number for the log value!"))
      v.log <- as.numeric(input$log)
      temp <- data.p()
      if (v.log >= 0){
        temp <- subset(temp, log2FoldChange >= v.log)
      } else {
        temp <- subset(temp, log2FoldChange <= v.log)
      }
      return(temp)
    }
  })
  
  data.p.log.human <- reactive({
    if(input$human.s == ""){
      return(data.p.log())
    } else{
      temp <- subset(data.p.log(), Human_Homolog == input$human.s)
      return(temp)
    }
  })
  
  data.p.log.human.org <- reactive({
    if (input$org.s == ""){
      return(data.p.log.human())
    } else {
      temp <- subset(data.p.log.human(), Gene_ID == input$org.s)
      return(temp)
    }
  })
  
  output$db <- renderDataTable(data.p.log.human.org())
  
  output$download <- downloadHandler(
    filename = paste0(Sys.Date(), "_exported_nf_table",".csv"),
    content = function(file) {
      write.csv(data.p.log.human.org(), file, row.names = FALSE)
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
