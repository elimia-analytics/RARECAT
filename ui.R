#' ---
#' title: NatureServe Rapid Rank Review Tool
#' ---
#'
#' # Load libraries
library(shiny)
library(leaflet)
library(leaflet.extras)
library(purrr)
library(shinyjs)
library(sf)
library(terra)
library(dplyr)
library(plotly)
library(htmltools)
library(htmlwidgets)
library(shinyWidgets)
library(sortable)
library(shinydashboard)
library(shinycssloaders)
library(shinyBS)
library(DT)
library(natserv)
library(rgbif)
library(flexdashboard)
library(shinybusy)
library(leafpm)
library(red)
library(dygraphs)
library(RWmisc)
library(units)

scr <- tags$script(HTML(
  "
Shiny.addCustomMessageHandler(
  'removeleaflet',
  function(x){
    console.log('deleting',x)
    // get leaflet map
    var map = HTMLWidgets.find('#' + x.elid).getMap();
    // remove
    map.removeLayer(map._layers[x.layerid])
  })
"
))
#'
#' # User Interface
navbarPage(title = HTML("<span style='float: left; display: inline-block; padding-left: 20px;'><img src = 'ns_logo.png', height = '45'></span><span style='display: inline-block; padding: 13px 5px 15px 15px;'><p style = 'font-size: 28px; font-family: Archivo !important;'><strong>RARECAT</strong></p></span>"), 
           windowTitle = "RARECAT", 
           id="nav", theme = "style.css", collapsible = TRUE,
           
           scr,
           
           tags$head(
             HTML("<link href='https://fonts.googleapis.com/css2?family=Roboto&display=swap' rel='stylesheet'>"),
             HTML("<link href='https://fonts.googleapis.com/css2?family=Archivo&display=swap' rel='stylesheet'>"),
             HTML("<meta name='viewport' content='width=device-width, initial-scale=1'>"),
           ),
           
           tabPanel("SINGLE SPECIES", height = "100%", 
                    
                    div(class="outer",
                        
                        useShinyjs(),     ## Call to use shinyJS
                        
                        scr,
                        
                        absolutePanel(id = "cond_inputs_panel",
                                      class = "panel panel-default",
                                      top = 57, left = 50, right = "auto", bottom = "auto",
                                      width = "25em",
                                      height = "3em",
                                      style = "margin-top: 0; padding: 0em 1.8em 1em 3em; border: none; border-bottom: none; border-color: transparent; background-color: rgba(169, 169, 169, 0); z-index: 1 !important; overflow-y: hidden !important; overflow-x: hidden;",
                                      textInput(inputId = "search_taxon", label = "", placeholder = "Select target taxon", width = "100%")
                        ),
                        
                        hidden(
                          absolutePanel(id = "taxon_search_panel", 
                                        class = "panel panel-default",
                                        top = 100, left = 90, right = "auto", bottom = "auto",
                                        width = "33vw",
                                        height = "25vh",
                                        style = "padding: 0.5em; border-bottom: none; border-color: transparent; background-color: rgba(255, 255, 255, 0.8); z-index: 1 !important; overflow-y: scroll; overflow-x: hidden;",
                                        shinycssloaders::withSpinner(DT::dataTableOutput("taxon_NS_table", width = "100%"), type = 7, proxy.height = "50px")
                          )
                        ),
                        
                        hidden(
                          absolutePanel(id = "taxon_options_panel", 
                                        class = "panel panel-default",
                                        top = 100, left = 90, right = "auto", bottom = "auto",
                                        width = "33vw",
                                        height = "25vh",
                                        style = "padding: 0.5em; border-bottom: none; border-color: transparent; background-color: rgba(255, 255, 255, 0.8); z-index: 1 !important; overflow-y: scroll; overflow-x: hidden;",
                                        div(style = "padding-left: 10px;", p(em("Click on all taxonomic concepts to include in assessment: "))),
                                        shinycssloaders::withSpinner(DT::dataTableOutput("taxon_options_table", width = "100%"), type = 7, proxy.height = "50px"),
                                        div(style = "padding-top: 20px; position: relative; float: right;", actionButton(inputId = "begin_assessment", label = "Start assessment", icon = icon("close"), block = TRUE, class = "btn-primary btn-sm", width = "100%"))
                          )
                        ),

                        absolutePanel(id = "cond_inputs_panel", 
                                      class = "panel panel-default", 
                                      top = 63, left = 380, right = "auto", bottom = "auto",
                                      width = "13em",
                                      height = "2.5em",
                                      style = "margin-top: 0; padding: 0em 1.8em 1em 3em; border-bottom: none; border-color: transparent; background-color: rgba(169, 169, 169, 0); z-index: 10 !important; overflow-y: hidden !important; overflow-x: hidden;", 
                                      div(style = "float: left;",
                                      actionButton(inputId = "clear_map", label = "Start fresh", icon = icon("close"), block = TRUE, class = "btn-primary btn-sm", width = "100%")
                                      )
                                      ),
                        
                        shinycssloaders::withSpinner(leafletOutput("main_map", height="60vh", width = "106vw"), type = 7),
                        
                          div(id = "analysis_panel",
                              absolutePanel(id = "cond_inputs_panel2", 
                                            class = "panel panel-default", 
                                            top = 58, left = "auto", right = 14, bottom = "auto",
                                            width = "33vw",
                                            height = "58vh",
                                            style = "margin-top: 0; padding: 0.3em 1.2em 1em 1.2em; background-color: white; box-shadow: -5px 5px 5px rgba(169, 169, 169, .8); z-index: 1000 !important; overflow-y: scroll; scrollbar-color: white;", 
                                            fluidRow(style = "padding-left: 20px",
                                                            htmlOutput("species_name")
                                            ),

                                            fluidRow(
                                              column(width = 6, style = "margin-bottom: -15px; padding-bottom: 0px;",
                                                     fluidRow(
                                                       column(width = 4, style = "padding-top: 25px; padding-right: 0; margin-right: 0;",
                                                              materialSwitch(inputId = "map_uploads", 
                                                                             label = "", 
                                                                             value = FALSE, 
                                                              )
                                                       ),
                                                       column(width = 8, style = "text-align: left;",
                                                              h4("Add records from CSV")                                                                     
                                                       )
                                                     ),
                                                     fileInput(inputId = "filedata", 
                                                               label = "",
                                                               accept = c("text/csv",
                                                                          "text/comma-separated-values,text/plain",
                                                                          ".csv"),
                                                               placeholder = "",
                                                               multiple = TRUE
                                                     )
                                              ),
                                              hidden(
                                                div(id = "load_data_panel", style = "padding-left: 10px;",
                                                    column(width = 6, style = "padding-left: 10px;",
                                                           fluidRow(
                                                             column(width = 4, style = "padding-top: 25px; padding-right: 0; margin-right: 0;",
                                                                    materialSwitch(inputId = "load_gbif_data", 
                                                                                   label = "", 
                                                                                   value = FALSE
                                                                    )
                                                             ),
                                                             column(width = 8, style = "text-align: left;",
                                                                    h4("Add GBIF records")
                                                             )
                                                           ),
                                                           fluidRow(style = "padding-left: 20px;",
                                                                    column(width = 4, style = "position: relative; float: left; padding-right: 0; padding-top: 0; padding-left: 0.3em;", 
                                                                           textInput(inputId = "number_gbif_occurrences", label = "", value = 1000)
                                                                    ),
                                                                    column(width = 5, style = "padding-top: 1em; padding-left: 0.5em;",
                                                                           p("records max.")
                                                                    )
                                                           )
                                                    )
                                                )
                                              )
                                                     ),

                                            hidden(
                                              div(id = "data_panel", style = "padding-left: 10px;",
                                                  fluidRow(style = "padding-left: 5px;",
                                                           column(width = 12, style = "padding-top: 0.4em;",
                                                                  shinycssloaders::withSpinner(htmlOutput("number_occurrences"), type = 7, proxy.height = "0px")
                                                           )
                                                  ),
                                                  br(),
                                                  fluidRow(style = "padding-top: 0.6em;",
                                                           column(width = 2, style = "position: relative; float: left; padding-top: 1em;",
                                                                  materialSwitch(inputId = "range_extent",
                                                                                 label = "",
                                                                                 value = FALSE
                                                                  )
                                                           ),
                                                           column(width = 5, style = "position: relative; float: left; padding-left: 0.3em;",
                                                                  fluidRow(
                                                                    column(width = 12, 
                                                                           p("Range Extent")
                                                                    )
                                                                  )
                                                           ),
                                                           hidden(
                                                             column(id = "EOO_panel", width = 5, style = "margin-left: 0px; padding-left: 0px;",
                                                                    fluidRow(style = "padding-top: 0;",
                                                                             column(width = 12, style = "padding-top: 0;",
                                                                                    htmlOutput("species_range_value")
                                                                             )
                                                                    )
                                                             )
                                                           )
                                                  ),
                                                  fluidRow(style = "padding-top: 0.6em;",
                                                           column(width = 2, style = "position: relative; float: left; padding-top: 1em;",
                                                                  materialSwitch(inputId = "area_of_occupancy",
                                                                                 label = "",
                                                                                 value = FALSE
                                                                  )
                                                           ),
                                                           column(width = 5, style = "position: relative; float: left; padding-left: 0.3em;",
                                                                  fluidRow(
                                                                    column(width = 12, 
                                                                           p("Area of Occupancy")
                                                                    )
                                                                  )
                                                           ),
                                                           hidden(
                                                             column(id = "AOO_panel", width = 5, style = "margin-left: 0px; padding-left: 0px;",
                                                                    fluidRow(style = "padding-top: 0;",
                                                                             column(width = 12, style = "padding-top: 0;",
                                                                                    shinycssloaders::withSpinner(htmlOutput("AOO_value"), type = 7, proxy.height = "0px")
                                                                             )
                                                                    )
                                                             )
                                                           )
                                                  ),
                                                  fluidRow(
                                                    column(width = 3, style = "position: relative; float: left; padding-left: 2em;",
                                                           selectInput(inputId = "grid_cell_size", label = "", choices = list("2 x 2 km" = 2, "1 x 1 km" = 1))
                                                    ),
                                                    column(width = 9, style = "position: relative; float: left; padding: 1em 0 1.5em 0em;",
                                                           HTML(paste0(p("Grid cell size")))
                                                    )
                                                  ),
                                                  fluidRow(style = "padding-top: 0.6em;",
                                                           column(width = 2, style = "position: relative; float: left; padding-top: 1em;",
                                                                  materialSwitch(inputId = "number_EOs",
                                                                                 label = "",
                                                                                 value = FALSE
                                                                  )
                                                           ),
                                                           column(width = 5, style = "position: relative; float: left; padding-left: 0.3em;",
                                                                  fluidRow(
                                                                    column(width = 12, 
                                                                           p("Number of Occurrences")
                                                                    )
                                                                  )
                                                           ),
                                                           hidden(
                                                             column(id = "EOcount_panel", width = 5, style = "margin-left: 0px; padding-left: 0px;",
                                                                    fluidRow(style = "padding-top: 0;",
                                                                             column(width = 12, style = "padding-top: 0;",
                                                                                    shinycssloaders::withSpinner(htmlOutput("EOcount_value"), type = 7, proxy.height = "0px")                                                                                   )
                                                                    )
                                                             )
                                                           )
                                                  ),
                                                  fluidRow(
                                                    column(width = 3, style = "position: relative; float: left; padding-left: 2em;",
                                                           textInput(inputId = "separation_distance", label = "", value = 1000)
                                                    ),
                                                    column(width = 9, style = "position: relative; float: left; padding: 1em 0 1.5em 0em;",
                                                           p("Separation distance (m)")
                                                    )
                                                  ),
                                                  fluidRow(style = "padding-top: 2px;",
                                                           column(width = 6,
                                                                  downloadButton(outputId = "download_occurrence_data", label = "Download analysis records", class = "btn-primary btn-sm", style = "width: 100%;")
                                                           ),
                                                           column(width = 6, 
                                                                  downloadButton(outputId = "download_rank_data", label = "Download rank data", class = "btn-primary btn-sm", style = "width: 100%;")
                                                           )
                                                  )
                                              )
                                            )
                              )
                          )
                        ),
                        
                        hidden(
                          div(id = "species_occurrences_table", style = "padding-right: 20px; padding-left: 20px",
                              
                              fluidRow(
                                column(width = 2, style = "padding-bottom: 20px; margin-left: 0; background-color: rgba(249, 249, 249, 1);",
                                       h3("Filters", style = "padding-bottom: 12px;"),
                                       br(),
                                       materialSwitch(inputId = "clean_occ", 
                                                      label = "Clean up GBIF records", 
                                                      value = TRUE, 
                                                      right = TRUE
                                                      
                                       ),
                                       br(),
                                       materialSwitch(inputId = "centroid_filter", 
                                                      label = "Remove centroids identified", 
                                                      value = FALSE, 
                                                      right = TRUE
                                                      
                                       ),
                                       br(),
                                       br(),
                                       dateRangeInput("year_filter", "Set time frame of records", format = "yyyy", start = "1900-01-01", end = "2023-01-01"),
                                       br(),
                                       br(),
                                       materialSwitch(inputId = "no_year", 
                                                      label = "Select only records with no year", 
                                                      value = FALSE, 
                                                      right = TRUE
                                       ),
                                       br(),
                                       br(),
                                       textInput( 
                                         inputId = "uncertainty_filter", 
                                         label = "Set spatial uncertainty of records (in meters):",
                                         value = "", 
                                         width = "100%"
                                       ),
                                       br(),
                                       fluidRow(
                                         column(width = 6, 
                                                selectizeInput(inputId = "nation_filter",
                                                               label = "Filter records by nation(s)",
                                                               choices = NULL,
                                                               multiple = TRUE, 
                                                               width = "100%"
                                                )
                                                ),
                                         column(width = 6, 
                                                selectizeInput(inputId = "states_filter",
                                                               label = "Filter records by subnation(s)",
                                                               choices = NULL,
                                                               multiple = TRUE, 
                                                               width = "100%"
                                                )      
                                         )
                                       ),
                                       br(),
                                       fluidRow(style = "padding-left: 15px; padding-right: 15px;",
                                                selectizeInput(inputId = "synonyms_filter",
                                                               label = "Select taxon concepts included",
                                                               choices = NULL,
                                                               multiple = TRUE, 
                                                               width = "100%"
                                                )  
                                       ),
                                       br(),
                                       fluidRow(style = "padding-left: 15px; padding-right: 15px;",
                                                selectizeInput(inputId = "sources_filter",
                                                               label = "Select data sources included",
                                                               choices = NULL,
                                                               multiple = TRUE, 
                                                               width = "100%"
                                                )  
                                       ),
                                       br(),
                                       ),
                                
                                
                                
                                column(width = 6, style = "background-color: transparent;", 
                                       fluidRow(style = "padding: 10px 10px 10px 20px; margin-right: 0px;",
                                                column(width = 8,
                                                       materialSwitch(inputId = "remove_selections",
                                                                      label = "",
                                                                      value = FALSE,
                                                                      inline = TRUE,
                                                                      right = TRUE
                                                       ),
                                                       h4("Remove selected records from dataset", style = "display: inline-block;")
                                                ),
                                                column(width = 4, style = "padding-right: 0; padding-left: 0;",
                                                       div(style = "float: right !important;", actionButton(inputId = "clear_selected_records", label = "Unselect records", icon = icon("close"), block = TRUE, class = "btn-primary btn-sm", width = "12em"))
                                                )
                                       ),
                                       fluidRow(style = "padding-left: 15px; padding-right: 15px; overflow-x: scroll; overflow-y: hidden; scrollbar-color: white;",
                                         DT::dataTableOutput("occurrences_table", height="40vh")
                                       )
                                ),
                                column(width = 4, 
                                       fluidRow(style = "padding-left: 20px;",
                                         h3("Temporal trend"),
                                       ),
                                       fluidRow(style = "padding: 25px 20px 2px 20px;",
                                                dygraphs::dygraphOutput("occurrences_barchart_full", height = "16vh", width = "100%")
                                       ),
                                       br(),
                                       fluidRow(style = "padding-left: 20px; padding-top: 5px; padding-bottom: 5px; background-color: rgba(249, 249, 249, 1);",
                                                column(width = 4,
                                                       fluidRow(style = "padding-bottom: 20px;",
                                                         h3("Summarize by time period", style = "padding-top: 0;"),
                                                       ),
                                                       fluidRow(style = "padding-bottom: 20px;",
                                                                dateRangeInput("period1", "Time period 1", format = "yyyy", start = "1980-01-01", end = "1994-01-01")
                                                       ),
                                                       fluidRow(style = "padding-bottom: 20px;",
                                                                dateRangeInput("period2", "Time period 2", format = "yyyy", start = "1995-01-01", end = "2009-01-01")
                                                       ),
                                                       fluidRow(style = "padding-bottom: 20px;",
                                                                dateRangeInput("period3", "Time period 3", format = "yyyy", start = "2010-01-01", end = Sys.Date())
                                                       )
                                                ),
                                                column(width = 8,
                                                       plotly::plotlyOutput("occurrences_barchart_period", height = "30vh")
                                                       )
                                       ) 
                                )
                              )
                          )
                        )
                        

),

tabPanel("DOCUMENTATION", height = "100%",
         fluidRow(style = "padding: 20px 50px 20px 50px; background-color: rgba(249, 249, 249, 1);",
           h2(paste0("You are using RARECAT version 1.1.1 (", Sys.Date(), "). For more information: "),
           strong(a(
             "NatureServe RARECAT Documentation.",
             target = "_blank",
             href = "https://natureserve01.sharepoint.com/:w:/g/teamsites/element_ranking/ETf6mA-3TvtMpo6iYO45JQ4BEhhuQFRKYfbUEbnlm_u7xA"
           ))),
           br(),
           h2("Suggested citation:"),
           h2(paste0("NatureServe (", substr(Sys.Date(), 1, 4), "). ", "RARECAT version 1.1.1. ",
                     "Available from https://natureserve.shinyapps.io/RARECAT. Accessed [Date].")
           )
         )
)
)
