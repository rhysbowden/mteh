#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(plotly)
library(ggplot2)
library(reshape2)
library(scales)
library(dplyr)

source("msteh_helper.R", local=TRUE)

# Define UI for application that draws a histogram
ui <- pageWithSidebar(
    h3("Estimating variance components when treatment effect heterogeneity is incorrectly excluded from the model."),
    sidebarPanel(
        actionButton("updatebutton", "Draw plot"),
        selectInput("design", label = ("Which cluster randomised design would you like to consider?"), 
                    choices = list("Stepped wedge" = 1, "CRXO" = 2, "SW/CXO ratio"=3), selected = 2),
        selectInput("panelvar", label = ("Vary clusters per period (c) or number of individuals per cluster period (N)?"), 
                    choices = list("c" = 1, "N" = 2), selected = 2),
        #selectInput("model", label = ("Compare Model 3 (correlation decay model) to:"), 
        #            choices = list("Model 1 (one-way model)" = 1, "Model 2 (two-way model)" = 2), selected = 1),
        
        # fileInput('file1', 'Upload a design matrix instead:',
        #           accept=c('text/plain', '.txt', '.csv')),
        # helpText("The file must be a comma separated .csv or .txt file consisting of 0s and 1s, with a column for each time period. Do not include row or column names."),
        # actionButton('reset', 'Clear file'),
        
        conditionalPanel(condition = "input.panelvar==1",
                         numericInput("m",
                                      "Number of subjects/cluster-period (N):",
                                      min = 1,
                                      max=1000,
                                      step = 1,
                                      value = 10)
        ),
        conditionalPanel(condition= "input.panelvar==2",
                         numericInput("c",
                                      "Number of clusters per sequence (c):",
                                      min = 1,
                                      max=1000,
                                      step = 1,
                                      value = 1)
        ),
        # sliderInput("rho0",
        #             label = HTML(paste("Intra-cluster correlation, &rho;",tags$sub(0),":", sep="")),
        #             #"rho 0:"
        #             min = 0.001,
        #             max = 0.2,
        #             step = 0.001,
        #             value = 0.05),
        sliderInput("rhoTTrange",
                    label = HTML(paste("Range of values of treatment-treatment intra-cluster correlation, &rho;",tags$sub("TT"),":", sep="")),
                    #"rho_TT:"
                    min = 0,
                    max = 1,
                    value = c(0.01,0.05),
                    step=0.001 ),
        sliderInput("rhoCC",
                    label = HTML(paste("Control-control intra-cluster correlation, &rho;",tags$sub("CC"),":", sep="")),
                    #"rho_CC:"
                    min = 0.001,
                    max = 0.2,
                    step = 0.001,
                    value = 0.01),
        sliderInput("rhoCT",
                    label = HTML(paste("Control-treatment intra-cluster correlation, &rho;",tags$sub("CT"),":", sep="")),
                    #"rho 0:"
                    min = -1,
                    max = 0.2,
                    step = 0.001,
                    value = 0.01),
        sliderInput("rhoCAC",
                    label = HTML(paste("Cluster auto-correlation, &rho;",tags$sub("CAC"),":", sep="")),
                    #"rho 0:"
                    min = 0.01,
                    max = 1,
                    step = 0.01,
                    value = 0.8),
        
        # helpText(HTML(paste("Note: calculations assume &sigma;",tags$sup(2),tags$sub(HTML(paste("1&epsilon;"))), "+", "&sigma;",
        #                     tags$sup(2),tags$sub(HTML(paste("1&alpha;"))), " = 1",  " so  &rho; = &sigma;",tags$sup(2),tags$sub(HTML(paste("1&alpha;"))), "/(&sigma;",tags$sup(2),tags$sub(HTML(paste("1&alpha;"))), "+", "&sigma;",
        #                     tags$sup(2),tags$sub(HTML(paste("1&epsilon;"))), ")  = &sigma;",tags$sup(2),tags$sub(HTML(paste("1&alpha;"))), sep="")))
        
        
    ),
    mainPanel(
        tabsetPanel(
            tabPanel("Ratio of variances of treatment effect estimator", value=1, plotOutput("ratiovarplot"),
            #htmlOutput("Plotexplan")
            ),
            tabPanel("Design Schematics", value=2, 
                     textOutput("text1"),
                     tableOutput("SWxmat"),
                     textOutput("text4"),
                     tableOutput("crxoxmat"),
                     # textOutput("text2"),
                     # tableOutput("pllelxmat"),
                     # textOutput("text3"),
                     # tableOutput("pllelbxmat")
                     ),
            #tabPanel("Ratios for user-input design", value=3,  plotlyOutput("mymat_ratiovarplot")),
            tabPanel("Ref. and contact details", value=4, htmlOutput("Contactdetails")) )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    observeEvent(input$updatebutton,{
    output$ratiovarplot <- renderPlot({
        #vartreat_versus3_plot(input$design, input$m, input$rho0, input$model)
        # design, panel_var, cc, ni, rhoCC, rhoCT, rhoCAC
        my_outer_graph(design=isolate(input$design), panel_var=isolate(input$panelvar), cc=isolate(input$c), ni=isolate(input$m), rhoCC=isolate(input$rhoCC), rhoCT=isolate(input$rhoCT), rhoCAC=isolate(input$rhoCAC), rhoTTrange=isolate(input$rhoTTrange))
    })
    })
    # observeEvent(input$rhoTTmin, {
    #     updateSliderInput(session, "rhoTTmax", min = input$rhoTTmin)
    # })
    
    observeEvent(input$rhoCC|input$rhoTTrange, {
        updateSliderInput(session, "rhoCT", max = rhoct_max(input$rhoCC,input$rhoTTrange[1],input$rhoTTrange[2])$objective)
    })
    
    # mymatrix <- reactiveValues(data=NULL)
    # 
    # # observe({
    # #     req(input$file1)
    # #     mymatrix$data <- read.csv(input$file1$datapath)
    # # })
    # 
    # observeEvent(input$reset, {
    #     mymatrix$data <- NULL
    #     reset('file1')
    # })
    # 
    # 
    # output$mymat_ratiovarplot <- renderPlotly({
    #     if(!is.null(mymatrix$data))   mymat <- mymatrix$data
    #     # else mymat <- SWdesmat(7)
    #     myplot <- vartreat_versus3_plotDESMAT(mymat, input$m, input$rho0, input$model)
    #     myplot
    # })
    
    output$Plotexplan <- renderUI({
        HTML(paste("The plot displays V&#770",
                   tags$sub(1), "/V",tags$sub(3),  " or  V&#770",tags$sub(2), 
                   "/V",tags$sub(3), ". V&#770", tags$sub(1), " is the variance of the 
        treatment effect estimator obtained using Model 1 if the expected values of 
        the variance components obtained using the Model 1 ANOVA formulas
        (with expectation under Model 3, the autoregressive model) 
        are used to estimate variance  components. V",tags$sub(3), " is the variance of 
        the treatment effect estimator under Model 3 using the correct Model 3 
        within-cluster correlation structure and the true value of the 
        decay parameters. Results are displayed for a range of
        study lengths and  a range of decay values.", sep=""))
    })
    
    output$text1 <- renderText({ 
        "An example of a stepped wedge design matrix:"
    })
    # output$text2 <- renderText({ 
    #     "An example of a parallel design matrix:"
    # })
    # output$text3 <- renderText({ 
    #     "An example of a parallel w/ baseline design matrix:"
    # })
    output$text4 <- renderText({ 
        "An example of a CRXO design matrix:"
    })
    
    
    output$SWxmat <- renderTable({
        head(SWdesmat(7), n=6)
    },digits=0)
    # output$pllelxmat <- renderTable({
    #     head(plleldesmat(6), n=6)
    # },digits=0)
    # output$pllelbxmat <- renderTable({
    #     head(pllelBLdesmat(7), n=6)
    # },digits=0)
    output$crxoxmat <- renderTable({
        head(CRXOdesmat(7), n=7)
    },digits=0)
    
    
    
    output$Contactdetails <- renderUI({
        HTML(paste("This Shiny app accompanies the paper &quot;Inference for the treatment effect in longitudinal cluster-randomized trials when treatment effect heterogeneity is ignored&quot; by Rhys Bowden, Andrew Forbes and Jessica Kasza.
             For questions or comments, please contact Rhys Bowden: 
            rhys.bowden &quot; at &quot; monash.edu", sep=""))
    })
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
