if (!require("shiny")) install.packages("shiny")
if (!require("shinydashboard")) install.packages("shinydashboard")
library(shiny)
library(shinydashboard)

shinyUI(
  dashboardPage(
    dashboardHeader(title ="Time Series Analysis"),
    dashboardSidebar( 
      sidebarMenu(
        
        fileInput("file1", "Choose CSV File",
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        menuItem("Home", tabName = "Home", icon = icon("arrow-right")),
        
        
        menuItem("Normalization", tabName = "Normalization", icon = icon("arrow-right")),
        
        
        menuItem("Reapeted Measures ANOVA", tabName = "Anova", icon = icon("arrow-right")),
        
        
        menuItem("Post-Hoc ANOVA", tabName = "panova", icon = icon("arrow-right")),
        
        menuItem("Reapeted Measures ANOVA Plots", tabName = "plots", icon = icon("arrow-right")),
        
        menuItem("Friedman Test", tabName = "Friedman", icon = icon("arrow-right")),
        
        menuItem("Friedman Test Plots", tabName = "fplot", icon = icon("arrow-right")),       
        
        menuItem("Post-Hoc Friedman Test", tabName = "pfriedman",icon = icon("arrow-right")),
        
        
        
        checkboxInput("checkbox1", label = "Normalization", value = FALSE),
        checkboxInput("checkbox2", label = "Reapeted Measures ANOVA", value = FALSE),
        checkboxInput("checkbox5", label = "ANOVA Plots", value = FALSE),
        checkboxInput("checkbox3", label = "Friedman Test", value = FALSE),
        checkboxInput("checkbox4", label = "Post hoc", value = FALSE),
        checkboxInput("checkbox6", label = "Friedman Test Plots", value = FALSE)
        
        
        
        
        
        
      )),
    
    dashboardBody(
      tabItems(
        
        tabItem(tabName = "Home", 
                h2("Mass Spectrometry Intensity Scores"),
                fluidRow(
                  box(title = "Time series analysis for repeated measures",  status = "success", solidHeader = T, width = "50px", 
                      img(src="uni.png"),
                      p(strong("Discription")),
                      p("This is a web-based Shiny application, built during Module 6 of MSc Bioinformatics at University of Birmingham."),
                      p("The main purpose of that application is to perform time series analysis for mass spectrometry-based metabolomics data."),
                      p("The application gives the opportunity to perform:"),
                      p("-Normalization of dataset"),
                      p("-Parametric statistical test (repeated measures ANOVA)"),
                      p("-Non-parametric statistical test (Friedman test)"),
                      p("-Post hoc test for ANOVA models)"),
                      p("(-Post hoc test for Freidman tests is under construction...)"),
                      p("-Plots and graphs for significantly changed metabolites collected from ANOVA"),
                      p("-Plots and graphs for significantly changed metabolites collected from Friedman test"),
                      p("For more details about how options work and their outputs, please visit GitHub link below"),
                      p("- - - - -"),
                      p(strong("Demo dataset")),
                      p("We bulit the application making the stastistical analysis for a given metabolomics dataset."),
                      p("The demo dataset was provided by Professor Warwick Dunn at University of Birmingham and was an untargeted metabolomics dataset."),
                      p("The biological study investigated 38 adults in three different age ranges (<30, 30-50 and >50) who undertook an exercise-related"),
                      p("intervention with muscle biopsies collected at six different time points, 3 time points before a 20-week exercise intervention and"),
                      p("3 time points after the 20-week exercise intervention. The 3 time points were fasted state, immediately after exercise and 1.5 hours"), 
                      p("after feeding. Muscle biopsies were analysed applying ultra-performance liquid chromatography-mass spectrometry (LC/MS) and raw data"),
                      p("was processed applying the open source software XCMS to generate the data matrix applied for the demonstration."),
                      p("For our purposes, the factors that taken into consideration were TIMEPOINTS and different SUBJECTS."),
                      p("- - - - -"),
                      p("The application was built using 'shiny' R package. (R version 3.5.3)"),
                      p("The code of the statistical analysis and the application are available in GitHub:"),
                      p("https://github.com/module6group-msc/Siny-web-application-for-mass-spectrometry-based-metabolomics-data-analysis"),
                      p("- - - - -"),
                      p("No programming skills required to use the application.")
                      )
                  
                  
                  
                ),
                fluidRow(
                  box(title = "Your Data Set", status = "success", solidHeader = T, width = "50px", tableOutput("contents" ))
                )
                
        ),
        tabItem(tabName = "Normalization", 
                h2("Mass Spectrometry intensity Scores Normalized"),
                
                fluidRow(
                  box(title = "Normalization", status = "success", solidHeader = T, width = "50px", tableOutput("contents1" ))
                ),
                fluidRow(
                  box(title = "Download Normalized Values of Your Data Set", status = "success", solidHeader = T,  downloadButton("download Norm","Download Data"))
                )
                
        ),
        tabItem(tabName = "Anova", 
                h2("Reapeted Measures ANOVA"),
                
                fluidRow(
                  box(title = "Significant Metabolites", status = "success", solidHeader = T, tableOutput("contents2" )),
                  
                  box(title = "Download Significant Data Set", status = "success", solidHeader = T,  downloadButton("download Anova","Download Data"))
                  
                )
                
        ),
        
        tabItem(tabName = "panova", 
                h1("Reapeted Measures ANOVA Post-Hoc"),
                
                fluidRow(
                  box(title = "Files of One Way ANOVA Analysis", status = "success", solidHeader = T, tableOutput("contents3")),
                  
                  box(title = "Download Post-Hoc for ANOVA Analysis zip ", status = "success", solidHeader = T,  downloadButton("download panova","Download Zip file"))
                  
                )
        ),
        
        tabItem(tabName = "plots", 
                h1("Reapeted Measures ANOVA Plots"),
                
                fluidRow(
                  box(title = "Number of Significant Metabolites Plots", status = "success", solidHeader = T, tableOutput("contents4")),
                  
                  box(title = "Download Plots for ANOVA Analysis ", status = "success", solidHeader = T,  downloadButton("download plots","Download Plots"))
                  
                )
        ),  
        
        tabItem(tabName = "Friedman",
                h1 ("Friedman Test"),
                fluidRow(
                  box(title = "Friedman Test", status = "success", solidHeader = T, tableOutput("contents5" )),
                  
                  box(title = "Download Friedman Test ", status = "success", solidHeader = T,  downloadButton("download friedman","Download Data"))
                  
                )),
        
        tabItem(tabName = "fplot",
                h1 ("Friedman Test Plots"),
                fluidRow(
                  box(title = "Friedman Test Plots", status = "success", solidHeader = T, tableOutput("contents6" )),
                  
                  
                  box(title = "Download Friedman Plots zip ", status = "success", solidHeader = T,  downloadButton("download fplots","Download Data"))
                )),
        
        tabItem(tabName = "pfriedman", 
                h1("Friedman Test Post-Hoc"),
                
                
                fluidRow(
                  box(title = "Under Construction",  status = "warning", solidHeader = T, 
                      "Coming Soon!!!")
                  #box(title = "Examble2", status = "success", solidHeader = T, tableOutput("contents5" )),
                  
                  # box(title = "Download Post-Hoc for Friedman's Test zip ", status = "success", solidHeader = T,  downloadButton("download pfriedman","Download Data"))
                )
        )
        
      ))
    
    
    
    
    
  )
)
  

