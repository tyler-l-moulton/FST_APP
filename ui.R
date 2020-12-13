#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


library(shiny)
library(Matrix)
library(tidyverse)
library(magrittr)
library(plotly)
library(adegenet)
library(poppr)
library(hierfstat)
library(DT)
library(shinythemes)
library(shinybusy)
library(shinyalert)
library(stringr)
library(leaflet)
library(maptools)
library(sf)
date <- parse_date("2020/10/21")
vers <- 1.0


walleye_gen <- read_rds("./Data/Bowles_walleye_2020_10_20.rds")
# walleye_gen <- read_rds("./Data/Bowles_walleye_2020_10_20_no_outlier.rds")

# Define UI for application 
shinyUI(fluidPage(

    useShinyalert(),
    
    
    add_busy_gif(
        src = "https://media2.giphy.com/media/1xneO4pA1mgnx1GKnx/giphy.gif?cid=ecf05e479yy84pgylkullpsmlqoa3po31tcfaz39p2wf1u9p&rid=giphy.gif"),   
    
    navbarPage("Fish Population Differentiation", theme = shinytheme("cosmo"),
               
               
               
               # browser(),        
               tabPanel("Data Import & Filtering",fluid = TRUE,
                        verticalLayout(
                            
                            titlePanel("How do you quantify genetic relatedness of populations?"),
                            uiOutput("FST_exp"),
                            
                            br(), br(), br(), hr(),
                            
                            
                            sidebarLayout(
                                sidebarPanel(
                                    radioButtons("species",
                                                 "Select Fish Species!",
                                                 choices = "Walleye" ),
                                    actionButton("GetData", "Get Data")
                                ),
                                
                                
                                mainPanel(
                                    
                                    titlePanel("Selecting your data!"),
                                    
                                    br(),br(),
                                    
                                    uiOutput("sp_sel_exp"),
                                    
                                    br(), br(),
                                    
                                    div(dataTableOutput(outputId = "tbl_head"), style = "font-size:80%"),
                                    
                                )
                            ),   #sidebarLayout
                            
                            hr(), ## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
                            
                            sidebarLayout(
                                sidebarPanel(
                                    actionButton("locmiss_eval",
                                                 "Find Fish Missing from Loci"),
                                    br(),
                                    sliderInput(inputId = "locmiss_threshold", step = .1,
                                                label = "Tolerance % Fish Missing from Loci",
                                                min = 0, max = 10,
                                                value = 2),
                                    br(),
                                    
                                    actionButton("make_loc_filter_plots", "Filter Out Loci")
                                    
                                ),
                                
                                mainPanel(
                                    titlePanel("Percent of Fish Genotyped by Each Locus"),
                                    br(),
                                    uiOutput("filter_loci_exp"),
                                    br(),
                                    splitLayout(
                                        
                                        dataTableOutput("locmiss_tbl"),
                                        
                                        plotOutput("locmiss_plot")
                                    )
                                )
                            ), #sidebarPanel
                            
                            hr(), ## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
                            
                            sidebarLayout(
                                sidebarPanel(
                                    actionButton("indmiss_eval",
                                                 "Find Loci Missing from Fish"),
                                    br(),
                                    
                                    sliderInput(inputId = "indmiss_threshold",
                                                label = "Tolerance % Loci Missing from Fish",
                                                min = 0, max = 50,
                                                value = 25),
                                    br(),
                                    
                                    actionButton("make_ind_filter_plots", "Filter Out Fish")
                                    
                                ),
                                mainPanel(
                                    titlePanel("Percent of Loci Genotyped for Each Fish"),
                                    
                                    br(),
                                    uiOutput("filter_fish_exp"),
                                    br(),
                                    
                                    
                                    splitLayout(
                                        
                                        dataTableOutput("indmiss_tbl"),
                                        
                                        plotOutput("indmiss_plot")
                                        
                                        
                                        
                                    )
                                )
                            ), #sidebarLayout
                            
                            hr(),## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
                            
                            sidebarLayout(
                                
                                sidebarPanel(
                                    actionButton("filter_gen_obj",
                                                 "Finalize Filters"),
                                    # actionButton("stop_app", "Interrupt Finalization!")
                                ),
                                
                                mainPanel(
                                    titlePanel("Summary of Filtering Steps"),
                                    
                                    
                                    tableOutput("summary_tbl"),
                                    
                                    uiOutput("filter_exp"),
                                    br(),
                                    
                                    p("The computation to determine genetic differentiation is very expensive. While you wait, you can watch the helpful video below!"),
                                    tags$a(href = "https://www.youtube.com/watch?v=I8RCOI7n4XI&ab_channel=iBiologyTechniques",
                                           "Genetic Differentiaion Video", target = "_blank"),
                                    p("when the computation is finished, scroll to the plot and click the 'Plots' tab."),
                                    br(),
                                    br(),
                                    br(),
                                    br()
                                    
                                    
                                )
                            )   #sidebarLayout
                            
                            
                        )#verticalLayout
               ), #tabPanel
               
               tabPanel("Plots", fluid = TRUE,
                        verticalLayout(
                            sidebarLayout(
                                
                                sidebarPanel(
                                    checkboxGroupInput(inputId = "location", label = "Select Locations", choices= character(), selected = character()),
                                    
                                    br(),
                                    
                                    checkboxGroupInput(inputId = "year", label = "Select Years", choices= character(), selected = character()),
                                    
                                    actionButton("make_plot",
                                                 "Make FST Plot")
                                ),
                                
                                mainPanel(
                                    titlePanel("Pairwise Genetic Difference Among Populations"),
                                    br(),
                                    uiOutput("FST_plot_exp"),
                                    br(),
                                    plotOutput("fst_heat")
                                    
                                ) #end mainPanel
                                
                                
                            ), #end sidebarLayout
                            
                            hr(),
                            
                            sidebarLayout(
                                sidebarPanel(
                                    actionButton("plot_yrs", 
                                                 HTML("Plot Genetic Distance <br/> Between Years")),
                                    br(),
                                    br(),
                                    actionButton("plot_sites", 
                                                 HTML("Plot Genetic Distance <br/> Among Sites"))
                                ),
                                mainPanel( 
                                    titlePanel("Genetic Distances Map"),
                                    br(),
                                    uiOutput("FST_map_exp"),
                                    br(),
                                    leafletOutput("map", width = 800, height = 800),
                                    br(), br(), br()
                                    
                                )
                            )
                            
                            
                        )
                        
               ), # end tab panel
               tabPanel("Documentation", fluid = TRUE,
                        verticalLayout(
                            titlePanel("Documentation"),
                            hr(),
                            br(),
                            h1("Data Availability and Associated Publications"),
                            br(),
                            p("Walleye Publications:"),
                            tags$a(href="https://onlinelibrary-wiley-com.lib-ezproxy.concordia.ca/doi/full/10.1111/eva.12987",
                                   "Bowles et al. 2020: Size reductions and genomic changes within two generations in wild walleye populations: associated with harvest? Evolutionary Applications.",
                                   target="_blank"),
                            br(), br(),
                            p("Walleye Data:"),
                            tags$a(href = "https://datadryad.org/stash/dataset/doi:10.5061/dryad.5tb2rbp1z",
                                   "Bowles et al. 2020 Evol. Appl. Data Dryad", target = "_blank"),
                            p("select: Bowles_et_al-EvolApp-2020-r12dWithOutOddBallSamps_includ-outlier_loci.genepop"),
                            p("This file was accessed October 20, 2020."),
                            br(),
                            hr(),
                            h1("Source Code:"),
                            tags$a(href = "https://github.com/tyler-l-moulton/FST_APP",
                                   "GitHub", target = "_blank"),
                            
                            br(), br(),
                            hr(),
                            br(),
                            p("App written by Tyler Moulton, December 12, 2020")
                            
                            
                            
                            
                            
                            
                            
                            
                        )
               ) 
    )
))
