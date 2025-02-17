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

#' ## Load NatureServe Network subnation polygons for subnation overlay
network_polys <- readRDS("data/subnation_polys.rds")

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
           
           tabPanel("SINGLE SPECIES MODE", height = "100%", 

                    tags$head(
                      HTML("<link href='https://fonts.googleapis.com/css2?family=Roboto&display=swap' rel='stylesheet'>"),
                      HTML("<link href='https://fonts.googleapis.com/css2?family=Archivo&display=swap' rel='stylesheet'>"),
                      HTML("<meta name='viewport' content='width=device-width, initial-scale=1'>"),
                    ),
                    
                    div(class="outer",
                        
                        useShinyjs(),     ## Call to use shinyJS
                        
                        scr,
                        
                        absolutePanel(id = "cond_inputs_panel",
                                      class = "panel panel-default",
                                      top = 57, left = 50, right = "auto", bottom = "auto",
                                      width = "25em",
                                      height = "3em",
                                      style = "margin-top: 0; padding: 0em 1.8em 1em 3em; border: none; border-bottom: none; border-color: transparent; background-color: rgba(169, 169, 169, 0); z-index: 1 !important; overflow-y: hidden !important; overflow-x: hidden; box-shadow: none !important;",
                                      textInput(inputId = "search_taxon", label = "", placeholder = "Select target taxon", width = "100%")
                        ),
                        
                        hidden(
                          absolutePanel(id = "taxon_search_panel", 
                                        class = "panel panel-default",
                                        top = 100, left = 90, right = "auto", bottom = "auto",
                                        width = "33vw",
                                        height = "25vh",
                                        style = "padding: 0.5em; border-bottom: none; border-color: transparent; background-color: rgba(255, 255, 255, 0.8); z-index: 1 !important; overflow-y: scroll; overflow-x: hidden; scrollbar-color: #fff !important;",
                                        shinycssloaders::withSpinner(DT::dataTableOutput("taxon_NS_table", width = "100%"), type = 7, proxy.height = "150px")
                          )
                        ),
                        
                        hidden(
                          absolutePanel(id = "taxon_options_panel", 
                                        class = "panel panel-default",
                                        top = 100, left = 90, right = "auto", bottom = "auto",
                                        width = "33vw",
                                        height = "25vh",
                                        style = "padding: 0.5em; border-bottom: none; border-color: transparent; background-color: rgba(255, 255, 255, 0.8); z-index: 1 !important; overflow-y: scroll; overflow-x: hidden; scrollbar-color:  rgba(255, 255, 255, 0.8) !important;",
                                        shinycssloaders::withSpinner(DT::dataTableOutput("taxon_options_table", width = "100%"), type = 7, proxy.height = "50px"),
                                        div(style = "padding-top: 20px; position: relative; text-align: center;", actionButton(inputId = "begin_assessment", label = "Start assessment", block = TRUE, class = "btn-primary btn-lg", width = "60%"))
                          )
                        ),
                        
                        shinycssloaders::withSpinner(leafletOutput("main_map", height="60vh", width = "106vw"), type = 7),
                        
                          div(id = "analysis_panel",
                              absolutePanel(id = "cond_inputs_panel2", 
                                            class = "panel panel-default", 
                                            top = 58, left = "auto", right = 14, bottom = "auto",
                                            width = "33vw",
                                            height = "58vh",
                                            style = "margin: 0; padding: 0.3em 0.2em 0.3em 0.3em; background-color: white; box-shadow: -5px 5px 5px rgba(169, 169, 169, .8); z-index: 1000 !important; overflow-y: scroll; scrollbar-color: #fff !important;", 
                                            fluidRow(style = "padding: 0px 10px 10px 10px",
                                                     column(width = 9,
                                                            htmlOutput("species_name")
                                                            ),
                                                     column(width = 3, style = "padding-top: 15px;",
                                                            actionButton(inputId = "clear_map", label = "Clear data", icon = icon("close"), block = TRUE, class = "btn-primary btn-sm", width = "100%")
                                                            )
                                            ),
                                            shinyBS::bsCollapse(id = "inputs_single", open = "Add assessment data",
                                                                shinyBS::bsCollapsePanel(title = "Add assessment data", style = "primary", 
                                            fluidRow(style = "padding-top: 10px;",
                                              column(width = 6, style = "margin-bottom: -15px; padding-bottom: 0px;",
                                                     fluidRow(
                                                       column(width = 4, style = "width: 21%; padding: 10px 10px 10px 20px; margin-right: 0;",
                                                              materialSwitch(inputId = "map_uploads", 
                                                                             label = "", 
                                                                             value = FALSE, 
                                                              )
                                                       ),
                                                       column(width = 8, style = "text-align: left; padding-top: 5px;",
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
                                                div(id = "load_data_panel",
                                                    column(width = 6, style = "padding-left: 10px;",
                                                           fluidRow(
                                                             column(width = 4, style = "width: 21%; padding: 10px 10px 10px 20px; margin-right: 0;",
                                                                    materialSwitch(inputId = "load_gbif_data", 
                                                                                   label = "", 
                                                                                   value = FALSE
                                                                    )
                                                             ),
                                                             column(width = 8, style = "text-align: left; padding-top: 5px;",
                                                                    h4("Add records from")
                                                             )
                                                             # column(width = 3, style = "position: relative; float: left; padding-right: 0; padding-top: 0; padding-left: 0.3em;", 
                                                             #        textInput(inputId = "number_gbif_occurrences", label = "", value = 1000)
                                                             # ),
                                                           ),
                                                           fluidRow(style = "padding: 0px 5px 0px 5px;",
                                                                    column(width = 12,
                                                                           selectizeInput("input_sources", label = "", choices = c("gbif", "inat", "ebird"), selected = c("gbif", "inat", "ebird"), multiple = TRUE)
                                                                           )
                                                           )
                                                    )
                                                )
                                              )
                                            )
                                                                )
                                                     ),

                                            hidden(
                                              div(id = "data_panel", style = "padding-top: 0; margin-top: 0;",
                                                  shinyBS::bsCollapse(id = "outputs_single", open = "Rank factor calculations",
                                                                      shinyBS::bsCollapsePanel(title = "Rank factor calculations", style = "primary", 
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
                                                  )
                                                                      )
                                                  ),
                                                  fluidRow(style = "padding-top: 2px;",
                                                           column(width = 4, style = "padding-right: 4px;",
                                                                  downloadButton(outputId = "download_occurrence_data", label = "Download records", class = "btn-primary btn-sm", style = "width: 100%;")
                                                           ),
                                                           column(width = 4, style = "padding-right: 2px; padding-left: 2px;", 
                                                                  downloadButton(outputId = "download_rank_data", label = "Download rank data", class = "btn-primary btn-sm", style = "width: 100%;")
                                                           ),
                                                           column(width = 4, style = "padding-left: 4px;",
                                                                  actionButton(inputId = "send_to_batch_mode", label = "Send to multispecies mode", block = TRUE, class = "btn-primary btn-sm", style = "width: 100%;")
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
                                column(width = 2, style = "margin-top: 5px; padding-bottom: 20px; margin-left: 0; background-color: white; border: 1px solid #347AB7;",
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
                                       dateRangeInput("year_filter", "Set time frame of records", format = "yyyy", start = "1900-01-01", end = Sys.Date(), startview = "decade", ),
                                       br(),
                                       br(),
                                       materialSwitch(inputId = "no_year", 
                                                      label = "Select only records with no year", 
                                                      value = FALSE, 
                                                      right = TRUE
                                       ),
                                       br(),
                                       h5("Select time of year", style = "color: #333333"),
                                       div(style = "padding: 0 10px 10px 10px;",
                                         sliderTextInput(
                                           inputId = "seasonality",
                                           label = "",
                                           choices = substr(month.name, 1, 3),
                                           selected = substr(month.name, 1, 3)[c(1, 12)], 
                                           grid = TRUE, 
                                           width = "100%"
                                         # column(width = 4, h5("Start date", style = "float: right; padding-top: 5px;")),
                                         # column(width = 5, style = "padding-right: 0;", selectInput(inputId = "seasonality_month1", label = "", choices = (purrr::map(1:12, function(x) x) %>% purrr::set_names(month.name)))),
                                         # column(width = 3, style = "padding-left: 0;", selectInput(inputId = "seasonality_day1", label = "", choices = 1:31, selected = 1))
                                       )
                                       ),
                                       # fluidRow(
                                       #   column(width = 4, h5("End date", style = "float: right; padding-top: 5px;")),
                                       #   column(width = 5, style = "padding-right: 0;", selectInput(inputId = "seasonality_month2", label = "", choices = (purrr::map(1:12, function(x) x) %>% purrr::set_names(month.name)), selected = 12)),
                                       #   column(width = 3, style = "padding-left: 0;", selectInput(inputId = "seasonality_day2", label = "", choices = 1:31, selected = 31))
                                       # ),
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
                                
                                
                                
                                column(width = 6, style = "background-color: transparent; width: 48vw; padding-left: 30px;", 
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
                                                       div(style = "float: right !important;", actionButton(inputId = "clear_selected_records", label = "Unselect records", block = TRUE, class = "btn-primary btn-sm", width = "12em"))
                                                )
                                       ),
                                       fluidRow(style = "padding-right: 15px; overflow-x: scroll; overflow-y: hidden; scrollbar-color: #fff !important;",
                                         DT::dataTableOutput("occurrences_table", height="40vh")
                                       )
                                ),
                                column(width = 4, style = "padding-top: 5px; width: 34vw;",
                                       tabsetPanel(type = "tabs",
                                         tabPanel("Change over time", style = "background-color: white; border-bottom: 1px solid #347AB7; border-left: 1px solid #347AB7; border-right: 1px solid #347AB7;",
                                                  fluidRow(style = "padding: 5px 35px 2px 35px;",
                                                           h3("Records per year"),
                                                           dygraphs::dygraphOutput("occurrences_barchart_full", height = "16vh", width = "100%")
                                                  ),
                                                  br(),
                                                  fluidRow(style = "padding: 0px 35px 20px 35px;",
                                                           h3("Rarity change by time period"),
                                                           br(),
                                                           column(width = 4,
                                                                  fluidRow(style = "padding-bottom: 20px;",
                                                                           dateRangeInput("period1", "Set time period 1", format = "yyyy", start = "1980-01-01", end = "1994-01-01")
                                                                  ),
                                                                  fluidRow(style = "padding-bottom: 20px;",
                                                                           dateRangeInput("period2", "Set time period 2", format = "yyyy", start = "1995-01-01", end = "2009-01-01")
                                                                  ),
                                                                  fluidRow(style = "padding-bottom: 20px;",
                                                                           dateRangeInput("period3", "Set time period 3", format = "yyyy", start = "2010-01-01", end = Sys.Date())
                                                                  )
                                                           ),
                                                           column(width = 8,
                                                                  selectInput("barchart_metric", label = "Select metric:", choices = c(list("Number of Records" = "rec_count", "Range Extent" = "eoo", "Area of Occupancy" = "aoo", "Number of Occurrences" = "eo_count")), multiple = FALSE, width = "70%"),
                                                                  plotly::plotlyOutput("metric_barchart_period", height = "30vh", width = "105%")
                                                           )
                                                  )
                                                  ),
                                         tabPanel("Temporal trend analysis", style = "background-color: white; border-bottom: 1px solid #347AB7; border-left: 1px solid #347AB7; border-right: 1px solid #347AB7;",
                                                  fluidRow(style = "padding-left: 5px; padding-bottom: 10px;",
                                                           column(width = 4,
                                                                  p("Select reference taxon ", style = "padding-top: 19px; float: right;")
                                                                  ),
                                                           column(width = 4, style = "padding-top: 5px; padding-left: 0;",
                                                                  selectizeInput(inputId = "select_reference_taxon",
                                                                                 label = "",
                                                                                 choices = c("genus", "family", "order", "class", "phylum", "kingdom"),
                                                                                 selected = "genus",
                                                                                 multiple = FALSE, 
                                                                                 width = "90%"
                                                                  )
                                                           ),
                                                           column(width = 4, style = "padding: 12px 30px 5px 0px; float: right;",
                                                                  actionButton(inputId = "temporal_trend", label = "Calculate temporal trend", block = TRUE, class = "btn-primary btn-sm", width = "100%")
                                                           ),
                                                           hidden(
                                                           fluidRow(id = "temporal_trend_plots", style = "padding: 25px 25px 10px 25px;",
                                                             plotlyOutput("temporal_trends_output", height = "800px")
                                                           )
                                                           )
                                                  )
                                                  )
                                       ),
                                       
                                )
                              )
                          )
                        )
                        

),

tabPanel("MULTISPECIES MODE", height = "100%", 
         
         tags$head(
           HTML("<link href='https://fonts.googleapis.com/css2?family=Roboto&display=swap' rel='stylesheet'>"),
           HTML("<link href='https://fonts.googleapis.com/css2?family=Archivo&display=swap' rel='stylesheet'>"),
           HTML("<meta name='viewport' content='width=device-width, initial-scale=1'>"),
         ),
         
         div(class="outer",
             
             useShinyjs(),     ## Call to use shinyJS
             
             scr,
             shinyBS::bsCollapse(id = "batch_parameters", open = "Rank assessment parameters", 
               shinyBS::bsCollapsePanel(title = "Rank assessment parameters", style = "primary", 
             fluidRow(style = "padding: 20px 20px 0px 20px; margin-bottom: 0;",
                      column(width = 2,
                             h1("1. Set assessment geography")
                      ),
                      column(width = 2,
                             h3("Select assessment type"),
                             selectizeInput(inputId = "batch_assessment_type",
                                            label = "",
                                            choices = c(list("Global (G)" = "global", "National (N)" = "national", "Subnational (S)" = "subnational")),
                                            multiple = FALSE, 
                                            width = "100%"
                             )
                      ),
                      column(width = 2, 
                             h3("Select assessment nation(s)"),
                             selectizeInput(inputId = "batch_nation_filter",
                                            label = "",
                                            choices = c(list("Canada" = "CA", "United States" = "US")),
                                            multiple = TRUE, 
                                            width = "100%"
                             )
                      ),
                      column(width = 2, 
                             h3("Select assessment subnation(s)"),
                             selectizeInput(inputId = "batch_states_filter",
                                            label = "",
                                            choices = (network_polys$Admin_abbr %>% na.omit() %>% as.character()) %>% set_names(network_polys$ADMIN_NAME%>% na.omit() %>% as.character()) %>% sort(),
                                            multiple = TRUE, 
                                            width = "100%"
                             )      
                      )
             ),
             fluidRow(style = "padding: 20px 20px 0px 20px;",
                      column(width = 2,
                             h1("2. Set assessment taxa")
                      ),
                      column(width = 3,
                             h3("Add taxon names from Rank Calculator file"),
                             fileInput(inputId = "batch_filedata_rank", 
                                       label = "",
                                       accept = c("text/csv",
                                                  "text/comma-separated-values,text/plain",
                                                  ".csv",
                                                  ".xslx"),
                                       placeholder = "",
                                       multiple = TRUE, 
                                       width = "100%"
                             )
                      ),
                      column(width = 3, 
                             h3("Add taxon names from observations CSV file"),
                             fileInput(inputId = "batch_filedata_obs", 
                                       label = "",
                                       accept = c("text/csv",
                                                  "text/comma-separated-values,text/plain",
                                                  ".csv"),
                                       placeholder = "",
                                       multiple = TRUE,
                                       width = "100%"
                             )
                      ),
                      column(width = 2, style = "padding: 0 10px 10px 10px;",
                             h3("Type or paste species list"),
                             textAreaInput(inputId = "typed_list", 
                                           label = "",
                                           height = "80px",
                                           width = "100%",
                                           resize = "vertical"
                                           
                             )
                      )
             ),
             fluidRow(style = "padding: 20px 20px 0px 20px;",
                      column(width = 2,
                             h1("3. Set additional filters")
                      ),
                      column(width = 2, style = "padding: 50px 10px 10px 10px;",
                             materialSwitch(inputId = "batch_clean_occ", 
                                            label = "Clean up GBIF records", 
                                            value = TRUE, 
                                            right = TRUE
                                            
                             ),
                             br(),
                             materialSwitch(inputId = "batch_centroid_filter", 
                                            label = "Remove centroids identified", 
                                            value = FALSE, 
                                            right = TRUE
                                            
                             )
                      ),
                      column(width = 2, style = "padding: 0 10px 10px 10px;",
                             h3("Set max. spatial uncertainty (m)"),
                             textInput( 
                               inputId = "batch_uncertainty_filter", 
                               label = "",
                               value = "", 
                               width = "90%"
                             )
                      ),
                      column(width = 2, style = "padding: 0 0 10px 10px;",
                             h3("Select data sources to include"),
                             selectizeInput(inputId = "batch_sources_filter",
                                            label = "",
                                            choices = c("gbif", "inat", "ebird", "uploaded"),
                                            multiple = TRUE, 
                                            width = "90%"
                             )
                      ),
                      column(width = 2,
                             h3("Select time frame"),
                             fluidRow(style = "padding-left: 10px;",
                             dateRangeInput("batch_year_filter", "", format = "yyyy", start = "1900-01-01", end = Sys.Date(), width = "90%"),
                             ),
                             fluidRow(style = "padding: 35px 0 0 20px;",
                             materialSwitch(inputId = "batch_no_year", 
                                     label = "Remove records with no year", 
                                     value = FALSE, 
                                     right = TRUE
                                     )
                             )
                      ),
                      column(width = 2,
                             h3("Select time of year"),
                             div(style = "padding-left: 10px;",
                               sliderTextInput(
                                 inputId = "seasonality",
                                 label = "",
                                 choices = substr(month.name, 1, 3),
                                 selected = substr(month.name, 1, 3)[c(1, 12)], 
                                 grid = TRUE, 
                                 width = "100%"
                             )
                      )
                      )
             ),
             fluidRow(style = "padding: 10px 20px 0px 20px;",
                      column(width = 2,
                             h1("4. Set calculation parameters")
                      ),
                      column(width = 2,
                             h3("Select AOO grid cell size"),
                             selectInput(inputId = "batch_grid_cell_size", label = "", choices = list("2 x 2 km" = 2, "1 x 1 km" = 1))
                      ),
                      column(width = 2,
                             h3("Select occurrence separation distance"),
                             textInput(inputId = "batch_separation_distance", label = "", value = 1000)
                      )
             )
               )
             ),
             fluidRow(style = "padding: 0px 20px 0px 20px; margin-top: 0;",
                      column(width = 2,
                             div(style = "padding: 0px 0px 10px 0px; position: relative; text-align:center; display:block;", actionButton(inputId = "batch_assessment", label = "Start assessment", block = TRUE, class = "btn-primary btn-lg", width = "100%")),
                             ),
                      column(width = 1,
                             div(style = "padding: 0px 10px 10px 0px; position: relative; text-align:center; display:block;", actionButton(inputId = "batch_clear", label = "Clear data", block = TRUE, class = "btn-primary btn-lg", width = "100%")),
                             )
             ),
             fluidRow(style = "padding: 10px 20px 0px 20px;",
                      column(width = 8,
                             hidden(
                             div(id = "batch_output",
                             DT::dataTableOutput("batch_run_results_table"),
                             downloadButton(outputId = "download_rank_data_batch", label = "Download rank data", class = "btn-primary btn-lg", style = "width: 30%; float: left;")
                             )
                             )
                      )
             )
         )
),

tabPanel("DOCUMENTATION", height = "100%",
         fluidRow(style = "padding: 20px 50px 20px 50px; background-color: rgba(249, 249, 249, 1);",
           h2(paste0("You are using RARECAT version 1.1.1 (2024-10-01). For more information: "),
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
