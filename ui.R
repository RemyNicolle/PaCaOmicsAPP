library(shiny)


shinyUI(bootstrapPage(

textInput(inputId="gene",
value="NPC1L1",
label="Enter gene symbol"
),

selectInput(inputId="data",label="Select dataset:",
choices=c("Tumor","Stroma","TumorVsStroma"),selected="Tumor"),

plotOutput(outputId = "main_plot")


))
