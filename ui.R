#' ---
#' title: NatureServe Rapid Rank Review Tool
#' ---
#'
#' # Load libraries
library(shiny)
library(leaflet)
library(leaflet.extras)
library(tidyverse)
library(shinyjs)
library(sf)
library(terra)
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

#' ## Add custom Javascript
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
tooltip_js <- "
$(function () {
  $('[data-toggle=tooltip]').tooltip()
})
"
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
                      tags$script(HTML(tooltip_js))
                    ),
                    
                    div(class="outer",
                        
                        useShinyjs(),     ## Call to use shinyJS
                        
                        scr,
                        
                        fluidRow(style = "padding: 0 0 10px 30px; margin-right: 0;",
                                 column(width = 2, style = "width: 15%; padding-left: 0; padding-right: 2px; z-index: 1005;",
                                        selectizeInput(inputId = "single_assessment_type",
                                                       label = "",
                                                       choices = c(list("Global Assessment (G)" = "global", "National Assessment (N)" = "national", "Subnational Assessment (S)" = "subnational")),
                                                       multiple = FALSE,
                                                       width = "100%",
                                                       options = list(placeholder = "Select assessment geography")
                                        )
                                 ),
                                 span(style = "float: left; padding: 0; margin-right: 5px;",
                                      `data-toggle` = "tooltip", `data-placement` = "right",  `data-animation` = "true",
                                      title = "Begin by selecting the geography of the assessment; if selecting a national or subnational assessment, also make a selection for the specific nation and subnation(s) from the corresponding dropdowns.",
                                      icon("info-circle", style = "color: #1F417D;")
                                 ),
                                 column(id = "single_nation", width = 2, style = "width: 14%; padding-left: 5px; padding-right: 2px;",
                                        selectizeInput(inputId = "single_assessment_nation",
                                                       label = "",
                                                       choices = "", 
                                                       multiple = FALSE,
                                                       width = "100%",
                                                       options = list(placeholder = "Select assessment nation")
                                        )
                                 ),
                                 column(id = "single_subnation", width = 2, style = "width: 14%; padding-left: 5px; padding-right: 2px;",
                                        selectizeInput(inputId = "single_assessment_subnation",
                                                       label = "",
                                                       choices = "",
                                                       multiple = TRUE,
                                                       width = "100%",
                                                       options = list(placeholder = "Select assessment subnation(s)")
                                        )
                                 ),
                                 column(width = 2, style = "width: 19%; padding-left: 5px; margin-left: 0; padding-right: 2px; display: inline !important; float: left;",
                                        div(textInput(inputId = "search_taxon", label = "", placeholder = "Select assessment taxon", width = "100%"))
                                 ),
                                 span(style = "float: left; padding: 0; margin-right: 5px;",
                                      `data-toggle` = "tooltip", `data-placement` = "right",  `data-animation` = "true",
                                      title = "Select the taxon to be assessed from the NatureServe or Global Biodiversity Information Facility (GBIF) taxonomic backbones, as well as the GBIF data to be included in the assessment.",
                                      icon("info-circle", style = "color: #1F417D;")
                                 ),
                                 column(id = "single_clear", width = 2, style = "padding-top: 5px; padding-left: 10px; font-size: 13px !important;",
                                        actionButton(inputId = "single_assessment_clear", label = "Start over", icon = icon("close"), block = TRUE, class = "btn-primary btn-lg", width = "90%", style = "font-size: 14px !important;")
                                 )
                        ),
                        
                        hidden(
                          absolutePanel(id = "taxon_search_panel", 
                                        class = "panel panel-default",
                                        top = 94, left = 15, right = "auto", bottom = "auto",
                                        width = "63%",
                                        height = "35vh",
                                        style = "padding: 0.5em; border-bottom: none; border-color: transparent; background-color: rgba(255, 255, 255, 1); z-index: 1001 !important; overflow-y: scroll; overflow-x: hidden; scrollbar-color: #C7C7C7 rgba(255, 255, 255, 1) !important; border: 1px solid rgba(0, 0, 0, 0.15); border-radius: 4px; box-shadow: 0 6px 12px rgba(0, 0, 0, 0.175);",
                                        DT::dataTableOutput("taxon_NS_table", width = "100%")
                          )
                        ),
                        
                        hidden(
                          absolutePanel(id = "taxon_options_panel", 
                                        class = "panel panel-default",
                                        top = 94, left = 15, right = "auto", bottom = "auto",
                                        width = "63%",
                                        height = "35vh",
                                        style = "padding: 0.5em; border-bottom: none; border-color: transparent; background-color: rgba(255, 255, 255, 1); z-index: 1001 !important; overflow-y: scroll; overflow-x: hidden; scrollbar-color: #C7C7C7 rgba(255, 255, 255, 1) !important; border: 1px solid rgba(0, 0, 0, 0.15); border-radius: 4px; box-shadow: 0 6px 12px rgba(0, 0, 0, 0.175);",
                                        span(style = "float: right; padding: 0; margin-right: 5px;",
                                             `data-toggle` = "tooltip", `data-html`="true", `data-placement` = "bottom",  `data-animation` = "true", 
                                             title = "- Select the GBIF taxon concepts and associated occurrence data relevant to assessment taxon <br/> - Click on corresponding row to select (in blue) or unselect (white) <br/> - Add to assessment either all GBIF records available for the assessment taxon and geography (up to 5000) or proceed to select specific datasets",
                                             icon("info-circle", style = "color: #1F417D;")
                                        ),
                                        DT::dataTableOutput("taxon_options_table", width = "100%"),
                                        br(),
                                        fluidRow(style = "padding-top: 20px; ", 
                                                 column(width = 6, style = "position: absolute; bottom: 20px; left: 10px;", actionButton(inputId = "begin_assessment_coarse", label = paste0("Start assessment with all records (up to 5000)"), block = TRUE, class = "btn-primary btn-lg", width = "95%", style = "font-size: 14px !important;")),
                                                 column(width = 6, style = "position: absolute; bottom: 20px; right: 5px;", actionButton(inputId = "select_datasets", label = "Select specific datasets", block = TRUE, class = "btn-primary btn-lg", width = "95%", style = "font-size: 14px !important;"))
                                        )
                          )
                        ),
                        
                        hidden(
                          absolutePanel(id = "taxon_datasets_panel",
                                        class = "panel panel-default",
                                        top = 94, left = 15, right = "auto", bottom = "auto",
                                        width = "63%",
                                        height = "55vh",
                                        style = "padding: 0.5em; border-bottom: none; border-color: transparent; background-color: rgba(255, 255, 255, 1); z-index: 1001 !important; overflow-y: scroll; overflow-x: hidden; scrollbar-color: #C7C7C7 rgba(255, 255, 255, 1) !important; border: 1px solid rgba(0, 0, 0, 0.15); border-radius: 4px; box-shadow: 0 6px 12px rgba(0, 0, 0, 0.175);",
                                        tabsetPanel(type = "tabs",
                                                    tabPanel("Occurrences", style = "padding: 10px; background-color: rgba(230, 239, 240, 0.5); border-color: transparent;",
                                                             fluidRow(
                                                               column(width = 2, style = "padding-top: 1em;",
                                                                      checkboxInput("select_all_occ", label = "Include all datasets", value = FALSE)
                                                               ),
                                                               column(width = 3, style = "padding-top: 1em;",
                                                                      checkboxInput("deselect_all_occ", label = "Exclude all datasets", value = FALSE)
                                                               ),
                                                               column(width = 5, offset = 2,
                                                                      p("Maximum records to include per dataset: ", style = "padding-top: 1em; padding-left: 30px; color: black; display: inline-block !important;"),
                                                                      span(style = "float: right; padding: 0; margin-right: 5px;",
                                                                           `data-toggle` = "tooltip", `data-html`="true", `data-placement` = "bottom",  `data-animation` = "true", 
                                                                           title = "- Select GBIF records to include in assessment from specific datasets <br/> - The 'Occurrences' tab lists biocollection datasets while the 'Human Observations' lists citizen and community science datasets <br/> - Click on corresponding row to select (in blue) or unselect (white) <br/> - The maximum number of records to include per datasets is a shared limit of records which will be used across all datasets of the corresponding type",
                                                                           icon("info-circle", style = "color: #1F417D;")
                                                                      ),
                                                                      span(textInput(inputId = "number_occ", label = "", value = 1000, width = "100%"), style = "display: inline-block !important; float:right; margin-left: 0; padding-right: 5px;")
                                                               )
                                                             ),
                                                             DT::dataTableOutput("taxon_datasets_occ_table", width = "100%"),
                                                    ),
                                                    tabPanel("Human Observations", style = "padding: 10px; background-color: rgba(230, 239, 240, 0.5); border-color: transparent;",
                                                             fluidRow(
                                                               column(width = 2, style = "padding-top: 1em;",
                                                                      checkboxInput("select_all_humobs", label = "Include all datasets", value = FALSE)
                                                               ),
                                                               column(width = 3, style = "padding-top: 1em;",
                                                                      checkboxInput("deselect_all_humobs", label = "Exclude all datasets", value = FALSE)
                                                               ),
                                                               column(width = 5, offset = 2,
                                                                      p("Maximum records to include per dataset: ", style = "padding-top: 1em; padding-left: 30px; color: black; display: inline-block !important;"),
                                                                      span(textInput(inputId = "number_humobs", label = "", value = 1000, width = "100%"), style = "display: inline-block !important; float:right; margin-left: 0;")
                                                               )
                                                             ),
                                                             DT::dataTableOutput("taxon_datasets_humobs_table", width = "100%"),                                                    ),
                                                    div(style = "padding-top: 20px; text-align: center; position: absolute; bottom: 20px; left: 37%; right: auto;", actionButton(inputId = "begin_assessment_detailed", label = "Start assessment with selected data", block = TRUE, class = "btn-primary btn-lg", width = "130%", style = "font-size: 14px !important;"))
                                        )
                          )
                        ),
                        
                        leafletOutput("main_map", height="60vh", width = "100vw"),
                        
                        add_busy_spinner(spin = "circle", color = "#1F417D", margins = c("40vh", "50vw"), height = "75px", width = "75px"),
                        
                        div(id = "analysis_panel",
                            absolutePanel(id = "cond_inputs_panel2", 
                                          class = "panel panel-default", 
                                          top = 115, left = "auto", right = 0, bottom = "auto",
                                          width = "35vw",
                                          height = "58vh",
                                          style = "margin: 0; padding: 0.3em 0 0 0.3em; background-color: white; z-index: 1000 !important; overflow-y: scroll; overflow-x: hidden; scrollbar-color: #C7C7C7 rgba(255, 255, 255, 1) !important; border: 1px solid rgba(0, 0, 0, 0.15); border-radius: 4px; box-shadow: 0 6px 12px rgba(0, 0, 0, 0.175);", 
                                          fluidRow(style = "padding: 0px 0px 10px 10px",
                                                   column(width = 9, style = "padding-top: 0; margin-top: 0;",
                                                          htmlOutput("species_name")
                                                   )
                                          ),
                                          shinyBS::bsCollapse(id = "inputs_single", open = "Add assessment data",
                                                              shinyBS::bsCollapsePanel(title = "Add assessment data", style = "primary", 
                                                                                       fluidRow(style = "padding-top: 10px;",
                                                                                                column(width = 6, style = "margin-bottom: -15px; padding-bottom: 0px;",
                                                                                                       fluidRow(
                                                                                                         column(width = 4, style = "width: 23%; padding: 10px 10px 10px 20px; margin-right: 0;",
                                                                                                                materialSwitch(inputId = "map_uploads", 
                                                                                                                               label = "", 
                                                                                                                               value = FALSE, 
                                                                                                                )
                                                                                                         ),
                                                                                                         column(width = 8, style = "text-align: left; padding-top: 5px;",
                                                                                                                h4("Add records from CSV", style = "display: inline !important; float: left; padding-right: 2px;"),
                                                                                                                span(style = "float: right; padding: 0; margin-right: 1px;",
                                                                                                                     `data-toggle` = "tooltip", `data-html`="true", `data-placement` = "bottom",  `data-animation` = "true", 
                                                                                                                     title = "- Select records from a formatted CSV file by browsing your local folders <br/> - Slide switch to the right to add records to the map <br/> - See Documentation for addictional details on how CSV files should be formatted",
                                                                                                                     icon("info-circle", style = "color: #1F417D;")
                                                                                                                ),
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
                                                                                                               column(width = 4, style = "width: 23%; padding: 10px 10px 10px 20px; margin-right: 0;",
                                                                                                                      materialSwitch(inputId = "load_gbif_data", 
                                                                                                                                     label = "", 
                                                                                                                                     value = FALSE
                                                                                                                      )
                                                                                                               ),
                                                                                                               column(width = 8, style = "text-align: left; padding-top: 5px;",
                                                                                                                      h4("Add records from GBIF", style = "display: inline !important; float: left; padding-right: 2px;"),
                                                                                                                      span(style = "float: right; padding: 0; margin-right: 1px;",
                                                                                                                           `data-toggle` = "tooltip", `data-html`="true", `data-placement` = "bottom",  `data-animation` = "true", 
                                                                                                                           title = "Add selected GBIF records to map by sliding switch to the right",
                                                                                                                           icon("info-circle", style = "color: #1F417D;")
                                                                                                                      ),
                                                                                                               )
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
                                            div(id = "data_panel", style = "padding-top: 0; margin-top: 0; padding-bottom: 0; margin-bottom: 2px !important;",
                                                shinyBS::bsCollapse(id = "outputs_single", open = "Rank factor calculations",
                                                                    shinyBS::bsCollapsePanel(title = "Rank factor calculations", style = "primary", 
                                                                                             fluidRow(style = "padding-left: 5px;",
                                                                                                      column(width = 12, style = "padding-top: 0;",
                                                                                                             column(width = 10, style = "padding-left: 0;", shinycssloaders::withSpinner(htmlOutput("number_occurrences"), type = 7, proxy.height = "0px"), style = "display: inline !important; float: left;"),
                                                                                                             span(style = "float: right; padding: 0; margin-right: 5px;",
                                                                                                                  `data-toggle` = "tooltip", `data-html`="true", `data-placement` = "bottom",  `data-animation` = "true", 
                                                                                                                  title = "- Slide each toggle to calculate the corresponding rarity metric using the records added to the map <br/> - Where applicable, change parameters for the calculations using the corresponding dropdown menu/text box <br/> - Use buttons at the bottom of this panel to download outputs and feed data back to multispecies mode",
                                                                                                                  icon("info-circle", style = "color: #1F417D;")
                                                                                                             ),
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
                                                                                                                               shinycssloaders::withSpinner(htmlOutput("species_range_value"), type = 7, proxy.height = "0px")
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
                                                                                               column(width = 4, style = "position: relative; float: left; padding-left: 2em;",
                                                                                                      selectInput(inputId = "grid_cell_size", label = "", choices = list("2 x 2 km" = 2, "1 x 1 km" = 1))
                                                                                               ),
                                                                                               column(width = 8, style = "position: relative; float: left; padding: 1em 0 0.5em 0em;",
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
                                                                                             fluidRow(style = "padding-bottom: 0; margin-bottom: 0;",
                                                                                                      column(width = 4, style = "position: relative; float: left; padding-left: 2em;",
                                                                                                             textInput(inputId = "separation_distance", label = "", value = 1000)
                                                                                                      ),
                                                                                                      column(width = 8, style = "position: relative; float: left; padding: 1em 0 0 0;",
                                                                                                             p("Separation distance (m)")
                                                                                                      )
                                                                                             )
                                                                    )
                                                ),
                                                fluidRow(style = "padding-top: 0;",
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
                      div(id = "species_occurrences_table", style = "padding-right: 20px; padding-left: 20px; height: 80vh;",
                          
                          fluidRow(
                            column(width = 2, style = "height: 80vh; margin-top: 5px; padding-top: 10px; padding-bottom: 20px; margin-left: 0; background-color: rgba(230, 239, 240, 0.5);",
                                   fluidRow(
                                     h3("Filters", style = "padding: 10px 2px 12px 10px; display: inline !important;"),
                                     span(style = "float: right; padding: 5px; margin-right: 5px;",
                                          `data-toggle` = "tooltip", `data-html`="true", `data-placement` = "bottom",  `data-animation` = "true", 
                                          title = "Use dropdown menus and text boxes below to apply filters to the data mapped and used to calculate rarity metrics and estimate temporal change",
                                          icon("info-circle", style = "color: #1F417D;")
                                     )
                                   ),
                                   div(style = "height: 70vh; padding: 0.4em; overflow-y: scroll; overflow-x: hidden; scrollbar-color: #C7C7C7 rgba(230, 239, 240, 0) !important;",
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
                                       selectizeInput(
                                         inputId = "seasonality",
                                         label = "",
                                         choices = substr(month.name, 1, 3),
                                         selected = substr(month.name, 1, 3), 
                                         multiple = TRUE,
                                         width = "100%"
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
                                       fluidRow(style = "padding-left: 15px; padding-right: 15px;",
                                                selectizeInput(inputId = "type_filter",
                                                               label = "Select record types included",
                                                               choices = NULL,
                                                               multiple = TRUE, 
                                                               width = "100%"
                                                )  
                                       ),
                                       br(),
                                       fluidRow(style = "padding-left: 15px; padding-right: 15px;",
                                                selectizeInput(inputId = "rank_filter",
                                                               label = "Select element occurrences by rank",
                                                               choices = NULL,
                                                               multiple = TRUE, 
                                                               width = "100%"
                                                )  
                                       )
                                   )
                            ),
                            
                            
                            
                            column(width = 6, style = "background-color: transparent; width: 47vw; padding-left: 30px; height: 80vh; ", 
                                   fluidRow(style = "padding: 10px 10px 10px 10px; margin-right: 0px;",
                                            column(width = 8, style = "padding-left: 5px;",
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
                                   fluidRow(style = "padding-right: 15px; overflow-x: scroll; overflow-y: hidden; scrollbar-color: #C7C7C7 rgba(255, 255, 255, 1) !important;",
                                            shinycssloaders::withSpinner(DT::dataTableOutput("occurrences_table", height="75vh"), type = 7, proxy.height = "0px")
                                   )
                            ),
                            column(width = 4, style = "padding-top: 5px; padding-right: 0; margin-right: 0; width: 34vw; height: 80vh !important;",
                                   tabsetPanel(type = "tabs",
                                               tabPanel("Change over time", style = "height: 76vh !important; width: 34vw; background-color: rgba(230, 239, 240, 0.5);",
                                                        fluidRow(style = "padding: 5px 35px 2px 35px;",
                                                                 h3("Records per year"),
                                                                 shinycssloaders::withSpinner(dygraphs::dygraphOutput("occurrences_barchart_full", height = "16vh", width = "100%"), type = 7, proxy.height = "0px")
                                                        ),
                                                        br(),
                                                        fluidRow(style = "padding: 0px 35px 20px 35px;",
                                                                 fluidRow(style = "padding-left: 15px;",
                                                                   h3("Rarity change by time period", style = "display: inline !important; float: left; padding-right: 2px;"),
                                                                   span(style = "float: right; padding: 10px; margin-right: 5px;",
                                                                        `data-toggle` = "tooltip", `data-html`="true", `data-placement` = "bottom",  `data-animation` = "true", 
                                                                        title = "Click on year text boxes to update the range of each time period and calculate the corresponding change in number of occurrences and rarity metrics across time periods.",
                                                                        icon("info-circle", style = "color: #1F417D;")
                                                                   )
                                                                 ),
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
                                                                        shinycssloaders::withSpinner(plotly::plotlyOutput("metric_barchart_period", height = "30vh", width = "105%"), type = 7)
                                                                 )
                                                        )
                                               ),
                                               tabPanel("Temporal analysis", style = "height: 76vh !important; width: 34vw; background-color: rgba(230, 239, 240, 0.5);",
                                                        fluidRow(style = "padding: 10px 10px 0px 10px;",
                                                          span(style = "float: right; display: inline !important; padding-right: 30px; margin-right: 5px;",
                                                               `data-toggle` = "tooltip", `data-html`="true", `data-placement` = "bottom",  `data-animation` = "true", 
                                                               title = "- Select the taxonomic scale against which assessment taxon records will be compared over time <br/> - Click the 'Calculate temporal change' button to estimate bias-corrected changes in the number of records and the probability of observing the assessment taxon across years <br/> - See Documentation for additional details on the methods used to run this analysis",
                                                               icon("info-circle", style = "color: #1F417D;")
                                                          )
                                                        ),
                                                        fluidRow(style = "padding: 0 20px 20px 20px;",
                                                                 # column(width = 3, style = "width: 20%;",
                                                                 #        p("Reference taxon", style = "padding-top: 19px; float: right;")
                                                                 # ),
                                                                 column(width = 4, style = "padding-top: 10px; padding-left: 10px;",
                                                                        selectizeInput(inputId = "select_reference_taxon",
                                                                                       label = "Reference taxon: ",
                                                                                       choices = c("genus", "family", "order", "class", "phylum", "kingdom"),
                                                                                       selected = "genus",
                                                                                       multiple = FALSE, 
                                                                                       width = "120%"
                                                                        )
                                                                 ),
                                                                 # column(width = 2, style = "width: 14%;",
                                                                 #        p("Start year ", style = "padding-top: 19px; padding-left: 0; padding-right: 2px; float:right;")
                                                                 # ),
                                                                 column(width = 3, style = "padding-top: 10px; padding-left: 10px; padding-right: 5px;",
                                                                        selectizeInput(inputId = "select_start_year",
                                                                                       label = "Start year:",
                                                                                       choices = 1900:2000,
                                                                                       selected = 1980,
                                                                                       multiple = FALSE, 
                                                                                       width = "100%"
                                                                        )
                                                                 ),
                                                                 column(width = 4, style = "width: 35%; padding: 12px 0 10px 5px;",
                                                                        div(actionButton(inputId = "temporal_trend", label = "Calculate temporal change", block = TRUE, class = "btn-primary btn-sm", width = "100%"), style = "padding: 10px 0px 5px 5px; display: inline !important; float:right;")
                                                                 ),
                                                                 hidden(
                                                                   fluidRow(id = "temporal_trend_plots", style = "padding: 50px 35px 10px 35px;",
                                                                            plotlyOutput("temporal_trends_output", height = "63vh")
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
                        
                        add_busy_spinner(spin = "circle", color = "#1F417D", margins = c("40vh", "50vw"), height = "75px", width = "75px"),
                        
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
                                                                              hidden(
                                                                                column(id = "batch_nation", width = 2,  
                                                                                       h3("Select assessment nation(s)"),
                                                                                       selectizeInput(inputId = "batch_nation_filter",
                                                                                                      label = "",
                                                                                                      choices = c(list("Canada" = "CA", "United States" = "US")),
                                                                                                      multiple = TRUE, 
                                                                                                      width = "100%"
                                                                                       )
                                                                                )
                                                                              ),
                                                                              hidden(
                                                                                column(id = "batch_subnation", width = 2, 
                                                                                       h3("Select assessment subnation(s)"),
                                                                                       selectizeInput(inputId = "batch_states_filter",
                                                                                                      label = "",
                                                                                                      choices = (network_polys$Admin_abbr %>% na.omit() %>% as.character()) %>% set_names(network_polys$ADMIN_NAME%>% na.omit() %>% as.character()) %>% sort(),
                                                                                                      multiple = TRUE, 
                                                                                                      width = "100%"
                                                                                       )
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
                                                                              column(width = 2, style = "padding: 0 0 10px 10px;",
                                                                                     h3("Select data sources to include"),
                                                                                     selectizeInput(inputId = "batch_sources_filter",
                                                                                                    label = "", 
                                                                                                    choices = c("OCCURRENCE", "HUMAN_OBSERVATION"),  
                                                                                                    selected = c("OCCURRENCE", "HUMAN_OBSERVATION"),
                                                                                                    multiple = TRUE, 
                                                                                                    width = "90%"
                                                                                     )
                                                                              ),
                                                                              column(width = 2, style = "padding: 50px 0px 10px 10px;",
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
                                                                              column(width = 2,
                                                                                     h3("Select time frame"),
                                                                                     fluidRow(style = "padding-left: 10px;",
                                                                                              dateRangeInput("batch_year_filter", "", format = "yyyy", start = "1900-01-01", end = Sys.Date(), startview = "decade", width = "90%"),
                                                                                     ),
                                                                                     fluidRow(style = "padding: 35px 0 0 20px;",
                                                                                              materialSwitch(inputId = "batch_no_year", 
                                                                                                             label = "Remove records with no year", 
                                                                                                             value = FALSE, 
                                                                                                             right = TRUE
                                                                                              )
                                                                                     )
                                                                              ),
                                                                              column(width = 2, style = "width: 16%",
                                                                                     h3("Select time of year"),
                                                                                     selectizeInput(
                                                                                       inputId = "batch_seasonality",
                                                                                       label = "",
                                                                                       choices = substr(month.name, 1, 3),
                                                                                       selected = substr(month.name, 1, 3), 
                                                                                       multiple = TRUE,
                                                                                       width = "100%"
                                                                                     )
                                                                              )
                                                                     ),
                                                                     fluidRow(style = "padding: 10px 20px 10px 20px;",
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
                                                                              ),
                                                                              column(width = 4,
                                                                                     h3("Calculate temporal change in rarity metrics"),
                                                                                     br(),
                                                                                     column(width = 6, style = "padding-left: 0;",
                                                                                            dateRangeInput("batch_period1", "Set time period 1", format = "yyyy", start = "1900-01-01", end = "2025-12-31")
                                                                                     ),
                                                                                     column(width = 6, style = "padding-left: 0;",
                                                                                            dateRangeInput("batch_period2", "Set time period 2", format = "yyyy", start = "1985-01-01", end = "2025-12-31")
                                                                                     )
                                                                              )
                                                                     )
                                            )
                        ),
                        fluidRow(style = "padding: 0px 20px 0px 20px; margin-top: 0;",
                                 column(width = 2,
                                        div(style = "padding: 0px 0px 10px 0px; position: relative; text-align:center; display:block;", actionButton(inputId = "batch_assessment", label = "Start assessment", block = TRUE, class = "btn-primary btn-lg", width = "100%")),
                                 ),
                                 column(width = 2,
                                        div(style = "padding: 0px 10px 10px 0px; position: relative; text-align:center; display:block;", actionButton(inputId = "batch_clear", label = "Clear data", block = TRUE, class = "btn-primary btn-lg", width = "100%")),
                                 )
                        ),
                        fluidRow(style = "padding: 10px 20px 0px 20px;",
                                 column(width = 12,
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
           tabPanel("TUTORIALS", height = "100%",
                    fluidRow(style = "padding: 20px 10px 20px 10px;",
                      column(width = 6, style = "padding-right: 0;",
                             tags$iframe(style='height:500px; width:95%;', src="https://www.youtube.com/embed/NTRjjfIb_wM?si=YLGXGqmwa61rcEgT#frameborder=0")
                      ),
                      column(width = 6, style = "padding-left: 0;",
                             tags$iframe(style='height:500px; width:95%;', src="https://www.youtube.com/embed/9eG5T_tWeME?si=-rZePIyEr5JrvQMZ#frameborder=0")
                      )
                    )
                    ),
           tabPanel("DOCUMENTATION", height = "100%",
                             # h2(paste0("You are using RARECAT version 2.1.1 (2025-03-31). For more information: view/download "),
                             #    strong(a(
                             #      "NatureServe RARECAT Documentation",
                             #      target = "_blank",
                             #      href = "https://natureserve01.sharepoint.com/:w:/g/teamsites/element_ranking/ETf6mA-3TvtMpo6iYO45JQ4BEhhuQFRKYfbUEbnlm_u7xA"
                             #    )), " or navigate pdf below."),
                             # br(),
                             # h2(paste0("To report issues or provide feedback, please submit tickets by either by logging into the"), 
                             #           strong(a("RARECAT Help Desk", 
                             #                    target = "_blank",
                             #                    href = "https://rarecatsupport.natureserve.org/"
                             #           )), paste0("or by emailing "), strong(a("rarecatsupport@natureserve.org", target = "_blank", href = "rarecatsupport@natureserve.org")), paste0("which will automatically create a ticket.")),
                             # br(),
                             # h2("Suggested citation:"),
                             # h2(paste0("NatureServe (", substr(Sys.Date(), 1, 4), "). ", "RARECAT version 1.1.1. ",
                             #           "Available from https://natureserve.shinyapps.io/RARECAT. Accessed [Date].")
                             # ),
                             tags$iframe(style='height:1000px; width:100%; scrolling=yes;', src="NatureServe RARECAT v2.1.1 - Documentation.pdf#zoom=125")
           )
)
