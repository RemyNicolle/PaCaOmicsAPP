#library(shiny)


navbarPage("PaCaOmics",




tabPanel("RNAseq",
sidebarLayout(


sidebarPanel(
textInput(inputId="gene",
value="NPC1L1",
label="Enter gene symbol"
),
selectInput(inputId="data",label="Select dataset:",
choices=c("Tumor","Stroma","TumorVsStroma"),selected="Tumor")
),



mainPanel(
plotOutput("rnaseq_plot")
)

)),






tabPanel("Methylation",
sidebarLayout(


sidebarPanel(
textInput(inputId="methgene",
value="NPC1L1",
label="Enter gene symbol or CpG ID:"
),
selectInput(inputId="island",label="Select island status:",
choices=c("any","Island","Shore","Shelf","non-CGI"),selected="any"),

selectInput(inputId="group",label="Select position in gene:",
choices=c("any","TSS1500",  "TSS200","1stExon","5'UTR",  "Body" ,"3'UTR"),selected="any"),

tags$p("Other possible CpG:"),
tableOutput("cpglist")
),



mainPanel(
plotOutput("methylation_plot")
)

)))




