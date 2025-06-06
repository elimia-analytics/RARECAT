#' ---
#' title: NatureServe RARECAT: Rapid Assessment of Rarity and Endangerment Conservation Assessment Tool
#' ---
#'
#' # Server setup
#' ## Load libraries
library(shiny)
library(leaflet)
library(leaflet.extras)
library(tidyverse)
library(shinyjs)
library(sf)
library(terra)
library(plotly)
library(ggpubr)
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
library(spocc)
library(flexdashboard)
library(shinybusy)
library(leafpm)
library(red)
library(dygraphs)
library(RWmisc)
library(units)
library(lme4)
library(writexl)
#'
#' ## Load NatureServe Network subnation polygons for subnation overlay
network_polys <- readRDS("data/network_polys.rds")
gadm_df <- readRDS("data/gadm_df.rds")
#'
#' ## Load Rank Calculator Domain Table Lookup
rank_factor_definitions <- read.csv("data/rank_factor_definitions.csv", header = TRUE)
#'
#' ## Load custom functions to enable app functionality 
source("RARECAT_functions.R")
#' 
#' ## Begin server source code
function(input, output, session) {
  
  ## Create web map and add basic elements and functionality
  output$main_map <- renderLeaflet({
    
    leaflet::leaflet(options = leafletOptions(zoomDelta = 0.5, zoomSnap = 0, attributionControl = FALSE, worldCopyJump = FALSE)) %>% # Open new leaflet web map
      leaflet::setView(lng = mean(c(-104.4474, -67.27911)), lat = 50, zoom = 3.5) %>%  # Zoom in on North America
      leaflet::addMapPane("basemap1", zIndex = -100) %>% # Add basemap 1
      leaflet::addProviderTiles(providers$Esri.WorldTerrain, group = "Esri World Terrain", options = list(pathOptions(pane = "basemap1")),
                                providerTileOptions(
                                  updateWhenZooming = FALSE,      # map won't update tiles until zoom is done
                                  updateWhenIdle = TRUE           # map won't load new tiles when panning
                                )) %>%
      leaflet::addMapPane("basemap2", zIndex = -100) %>% # Add basemap 2
      leaflet::addProviderTiles(providers$Esri.WorldTopoMap, group = "Esri Topographic Map", options = list(pathOptions(pane = "basemap2")),
                                providerTileOptions(
                                  updateWhenZooming = FALSE,      # map won't update tiles until zoom is done
                                  updateWhenIdle = TRUE           # map won't load new tiles when panning
                                )) %>%
      leaflet::addMapPane("basemap3", zIndex = -100) %>% # Add basemap 3
      leaflet::addProviderTiles(providers$Esri.WorldImagery, group = "Esri World Imagery", options = list(pathOptions(pane = "basemap3"))) %>%
      leaflet::addMapPane("basemap4", zIndex = -100) %>% # Add basemap 4
      leaflet::addProviderTiles(providers$OpenStreetMap, group = "Open Street Map", options = list(pathOptions(pane = "basemap4"))) %>%
      leaflet::addMapPane("basemap6", zIndex = -100) %>% # Add basemap 5
      leaflet::addProviderTiles(providers$Esri.WorldStreetMap, group = "Esri World Street Map", options = list(pathOptions(pane = "basemap5"))) %>%
      leaflet::addScaleBar(position = "bottomleft") %>% # Add scale bar
      leaflet.extras::addResetMapButton() %>% # Add button to reset map bounds
      leaflet.extras::addSearchOSM() %>% # Add functionality to search for specific location using Open Street Map
      leaflet::addLayersControl(baseGroups = c("Esri World Street Map", "Esri Topographic Map", "Open Street Map", "Esri World Terrain", "Esri World Imagery"), # Add layers control widget
                                options = layersControlOptions(collapsed = TRUE), position = "topleft") %>% 
      leafpm::addPmToolbar(toolbarOptions = leafpm::pmToolbarOptions(drawCircle = FALSE, drawPolyline = FALSE, editMode = FALSE, cutPolygon = FALSE, removalMode = FALSE), # Add point/polygon drawing tools
                           drawOptions = leafpm::pmDrawOptions(snappable = FALSE, markerStyle = list(draggable = FALSE))
      )
  })
  
  #' ## Ensure markers are not drawn by default
  observeEvent(input$main_map_draw_new_feature, {
    session$sendCustomMessage("removeleaflet", list(elid="main_map", layerid=input$main_map_draw_new_feature$properties$`_leaflet_id`))
  })
  #'
  #' ## Create static objects
  #' ### Specify minimumcommon fields that  all (up)loaded data should share
  minimum_fields <- c("key", "scientificName", "prov", "longitude", "latitude", "coordinateUncertaintyInMeters", "stateProvince", "countryCode", "year", "month", "datasetName", "institutionCode", "basisOfRecord", "EORANK", "references")
  #'
  #' ## Create reactive objects to store output
  #' ### Object to store NatureServe element selected
  selected_taxon <- reactiveValues(name = NULL)
  #' ### Object to store all data and outputs for the focal taxon
  taxon_data <- reactiveValues(
    info = data.frame(scientificName = "New taxon"),
    info_extended = NULL,
    synonyms = NULL,
    synonyms_selected = NULL,
    datasets = NULL,
    datasets_selected = NULL,
    assessment_polygon = NULL,
    gbif_occurrences_raw = NULL,
    gbif_occurrences = NULL,
    uploaded_occurrences = NULL,
    drawn_occurrences = NULL,
    all_occurrences = NULL,
    shifted = FALSE,
    sf = NULL,
    sf_filtered = NULL,
    filtered_occurrences = NULL,
    selected_points = data.frame("Key" = character(), "Scientific name" = character(), "Source" = character(), "Institution code" = character(), "Year" = numeric(), "Coordinate Uncertainty" = numeric(), "Place" = character(), "URL" = character())[NULL, ],
    removed_points = NULL,
    nations = NULL,
    states = NULL,
    records_over_time = NULL,
    species_range_value = NULL,
    species_range_map = NULL,
    species_range_factor = NULL,
    AOO_value = NULL,
    AOO_map = NULL,
    AOO_factor = NULL,
    EOcount_map = NULL,
    EOcount_value = NULL,
    EOcount_factor = NULL,
    rank_factor_comparison = NULL,
    temporal_change = NULL
  )
  #' ### Object to store all clicked point IDs
  clicks <- reactiveValues(IDs = vector(mode = "character"))
  #' ### Object to store "Begin assessment" button presses
  assessment_start_coarse_button_presses <- reactiveValues(values = 0)
  assessment_start_detailed_button_presses <- reactiveValues(values = 0)
  select_datasets_button_presses <- reactiveValues(values = 0)
  #'
  #' ## Set up reactive expressions
  #' ### Load NatureServe API information for taxon entered in search bar
  taxon_NS_options <- reactive({
    
    ns_table <- natserv::ns_search_spp(text_adv = list(searchToken = input$search_taxon, matchAgainst = "allScientificNames", operator="contains"))$results
    gbif_table <- rgbif::name_suggest(q = input$search_taxon, rank = c("species", "subspecies", "variety", "infraspecific_name"), limit = 10)$data
    
    out <- NULL
    
    if (nrow(ns_table) > 0){
      ns_table <- ns_table %>% 
        dplyr::mutate(Source = "NatureServe", synonyms = ns_table$speciesGlobal$synonyms, phylum = ns_table$speciesGlobal$phylum, kingdom = ns_table$speciesGlobal$kingdom) %>% 
        dplyr::select(scientificName, Source, elementGlobalId, primaryCommonName, roundedGRank, elcode, uniqueId, synonyms, phylum, kingdom)
      out <- rbind(out, ns_table)
    } 

      if (nrow(gbif_table) > 0){
        gbif_table <- gbif_table %>% 
          dplyr::rename(scientificName = canonicalName, elementGlobalId = key) %>% 
          dplyr::mutate(Source = "GBIF", primaryCommonName = NA, roundedGRank = NA, elcode = NA, synonyms = NA, uniqueId = NA, phylum = NA, kingdom = NA) %>% 
          dplyr::select(scientificName, Source, elementGlobalId, primaryCommonName, roundedGRank, elcode, uniqueId, synonyms, phylum, kingdom)
      out <- rbind(out, gbif_table)
      }
    
    out
    
  })
  
  download_gbif_data <- reactive({
    
    if (!is.null(selected_taxon$name)) {
      
      gbif_download <- purrr::safely(get_gbif_data)(
        taxa_metadata = taxon_data$synonyms_selected, 
        datasets_metadata = taxon_data$datasets_selected,
        query_polygon = taxon_data$assessment_polygon,
        all_occ_data = input$select_all_occ,
        all_humobs_data = input$select_all_humobs,
        shift_occurrences = TRUE
        )
      
      if(!is.null(gbif_download$result)){
      
      gbif_download <- gbif_download$result
      
      if (nrow(gbif_download$sp_occurrences) > 0){
        
        gbif_download
      
      } else {
        
        # shinybusy::remove_modal_spinner() # show the modal window
        sendSweetAlert(session, type = "warning", title = "Oops!", text = "There are no valid occurrences on GBIF for this taxon or any of its synonyms recognized by NatureServe", closeOnClickOutside = TRUE)
      }
      } else {
        if (!is.null(gbif_download$error)){
          if (grepl("gbif", gbif_download$error)){
            sendSweetAlert(session, type = "warning", title = "Oops!", text = "RARECAT failed to get a response from the GBIF API; please wait a few minutes and try again", closeOnClickOutside = TRUE)
          }
        } else {
          sendSweetAlert(session, type = "warning", title = "Oops!", text = "RARECAT could not load GBIF data at this time; please wait a few minutes and try again", closeOnClickOutside = TRUE)
        }
      } 
    } else {
      # shinybusy::remove_modal_spinner() # show the modal window
      sendSweetAlert(session, type = "warning", title = "Oops!", text = "You need to select a taxon before loading GBIF data!", closeOnClickOutside = TRUE)
    }
    
  })
  
  #' ### Load user-uploaded data
  uploaded_data <- reactive({
    
    
    out <- purrr::map(input$filedata$datapath, process_user_data, minimum_fields = c(minimum_fields, "scientificName_Assessment")) %>% 
      dplyr::bind_rows() 

    out
    
  })
  #'
  #' ## Specify what happens when user enters text in search bar
  observeEvent(input$search_taxon, {
    
    if (input$search_taxon != ""){
      
      shinyjs::show(id = "taxon_search_panel")
      shinyjs::hide(id = "taxon_options_panel")
      shinyjs::hide(id = "taxon_datasets_panel")
      updateMaterialSwitch(session = session, inputId = "range_extent", value = FALSE)
      updateMaterialSwitch(session = session, inputId = "area_of_occupancy", value = FALSE)
      updateMaterialSwitch(session = session, inputId = "number_EOs", value = FALSE)
      
      selected_taxon_info <- taxon_NS_options() 
      
      if (!is.null(selected_taxon_info)){
        
        output$taxon_NS_table <- DT::renderDataTable({
          
          selected_taxon_info <- selected_taxon_info %>% 
            dplyr::mutate(elementGlobalId = purrr::map(1:length(elementGlobalId), function(i){
              if (Source[i] == "NatureServe"){
                paste0("<a href='", paste0("https://explorer.natureserve.org/Taxon/", uniqueId[i]), "' target='_blank'>", elementGlobalId[i], "</a>")
              } else if (Source[i] == "GBIF"){
                paste0("<a href='", paste0("https://www.gbif.org/species/", elementGlobalId[i]), "' target='_blank'>", elementGlobalId[i], "</a>")
              } else {
                elementGlobalId
              }
            }) %>% unlist()) %>% 
            dplyr::select(-synonyms, -elcode, -uniqueId) %>% 
            dplyr::rename("Scientific name" = scientificName, 
                          "ID" = elementGlobalId, 
                          "Common name" = primaryCommonName, 
                          "G rank" = roundedGRank,
                          "Phylum" = phylum,
                          "Kingdom" = kingdom
                          ) %>% 

            DT::datatable(options = list(dom = 't', pageLength = 10, autoWidth = TRUE), selection = list(mode = 'single', target = 'row'), escape = FALSE, rownames = FALSE)
          
        })
        
      }
      
    } else {
      shinyjs::hide(id = "taxon_search_panel")
      # shinyjs::hide(id = "taxon_options_panel")
      # shinyjs::hide(id = "taxon_datasets_panel")
    }
    
  })
  
  
  observeEvent({
    input$taxon_NS_table_rows_selected
  }, {
    
    # shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
    
    updateTextInput(inputId = "search_taxon", value = "")
    
    taxon_data$info <- taxon_NS_options()[input$taxon_NS_table_rows_selected, ]
    
    if (taxon_data$info$Source == "NatureServe"){
      taxon_data$info_extended <- natserv::ns_id(uid = taxon_data$info$uniqueId)
      selected_taxon$name <- c(taxon_data$info$scientificName, unlist(taxon_data$info$synonyms))
      selected_taxon$name <- gsub("ssp. |var. ", "", selected_taxon$name)
    } else if (taxon_data$info$Source == "GBIF"){
      selected_taxon$name <- c(taxon_data$info$scientificName)
      taxon_data$info_extended <- NULL
    }
    
    if (length(selected_taxon$name) == 1){
      taxon_data$synonyms <- rgbif::name_usage(name = selected_taxon$name)$data
    } else {
      taxon_data$synonyms <- purrr::map(selected_taxon$name, function(sp) rgbif::name_usage(name = sp)) %>%
        purrr::map("data") %>%
        bind_rows()
    }
    
    taxon_data$synonyms <- taxon_data$synonyms %>% 
      dplyr::distinct(., .keep_all = TRUE) 

      gbif_counts <- purrr::map_dbl(taxon_data$synonyms$key, 
                                    function(x) rgbif::occ_count(
                                      taxonKey = x,
                                      gadmGid = taxon_data$assessment_polygon,
                                      hasCoordinate = TRUE)
                                    )

    taxon_data$synonyms <- taxon_data$synonyms %>% 
      dplyr::mutate(occurrence_count = gbif_counts) %>% 
      dplyr::filter(occurrence_count > 0)
    
    shinyInput <- function(FUN, len, id, ...) {
      inputs <- character(len)
      for (i in seq_len(len)) {
        inputs[i] <- as.character(FUN(paste0(id, i), ...))# label = NULL, ...))
      }
      inputs
    }
    
    DT::dataTableProxy("taxon_NS_table") %>% 
      selectRows(selected = NULL)
    
    shinyjs::show(id = "taxon_options_panel")
    shinyjs::show(id = "analysis_panel")
    updateCollapse(session = session, id = "inputs_single", open = "Add assessment data")
    shinyjs::hide(id = "taxon_search_panel")
    shinyjs::hide(id = "taxon_datasets_panel")
    
    m <- leafletProxy("main_map") %>%
      clearShapes() %>% 
      clearMarkers() %>% 
      clearMarkerClusters() 
    taxon_data$datasets = NULL
    taxon_data$datasets_selected = NULL
    taxon_data$gbif_occurrences_raw = NULL
    taxon_data$gbif_occurrences = NULL
    taxon_data$uploaded_occurrences = NULL
    taxon_data$drawn_occurrences = NULL
    taxon_data$all_occurrences = NULL
    taxon_data$shifted = FALSE
    taxon_data$sf = NULL
    taxon_data$sf_filtered = NULL
    taxon_data$filtered_occurrences = NULL
    taxon_data$selected_points = data.frame("Key" = character(), "Scientific name" = character(), "Source" = character(), "Institution code" = character(), "Year" = numeric(), "Coordinate Uncertainty" = numeric(), "Place" = character(), "URL" = character())[NULL, ]
    taxon_data$removed_points = NULL
    taxon_data$nations = NULL
    taxon_data$states = NULL
    taxon_data$records_over_time = NULL
    taxon_data$species_range_value = NULL
    taxon_data$species_range_map = NULL
    taxon_data$species_range_factor = NULL
    taxon_data$AOO_value = NULL
    taxon_data$AOO_map = NULL
    taxon_data$AOO_factor = NULL
    taxon_data$EOcount_map = NULL
    taxon_data$EOcount_value = NULL
    taxon_data$EOcount_factor = NULL
    taxon_data$rank_factor_comparison = NULL
    taxon_data$temporal_change = NULL
    
    updateMaterialSwitch(session = session, inputId = "load_gbif_data", value = FALSE)
    updateMaterialSwitch(session = session, inputId = "map_uploads", value = FALSE)    
    updateMaterialSwitch(session = session, inputId = "range_extent", value = FALSE)
    updateMaterialSwitch(session = session, inputId = "area_of_occupancy", value = FALSE)
    updateMaterialSwitch(session = session, inputId = "number_EOs", value = FALSE)
    
    shinyjs::show(id = "analysis_panel")
    shinyjs::hide(id = "data_panel")
    
    # shinybusy::remove_modal_spinner()

  })
  
  output$taxon_options_table <- DT::renderDataTable({
    
    taxon_data$synonyms %>% 
      dplyr::mutate(key = paste0("<a href='", paste0("https://www.gbif.org/species/", key), "' target='_blank'>", key, "</a>")) %>%         
      dplyr::select(scientificName, key, occurrence_count) %>% 
      # dplyr::mutate(add = ifelse(add == FALSE, as.character(icon("xmark", "fa-2x", style = "color: #ef8a62;")), as.character(icon("check", "fa-2x", style = "color: #67a9cf;")))) %>% 
      dplyr::rename("Scientific name" = scientificName, 
                    "Key" = key,
                    "Number of GBIF records" = occurrence_count
      ) %>% 
      DT::datatable(options = list(dom = 't', pageLength = 10, autoWidth = TRUE,
                                   columnDefs = list(list(className = "text-left", width = '400px', targets = c(0)),
                                                     list(className = "text-left", targets = c(1, 2))
                                   )
      ), 
      selection = list(mode = 'multiple', target = 'row', selected = 1:nrow(taxon_data$synonyms)), 
      escape = FALSE, 
      rownames = FALSE
      )
    
  })

  observeEvent({
    input$begin_assessment_coarse
  }, {
    
    assessment_start_coarse_button_presses$values <- c(assessment_start_coarse_button_presses$values, input$begin_assessment_coarse)
    
    if (assessment_start_coarse_button_presses$values[length(assessment_start_coarse_button_presses$values)] != assessment_start_coarse_button_presses$values[length(assessment_start_coarse_button_presses$values)-1]){
      taxon_data$synonyms_selected <- taxon_data$synonyms[input$taxon_options_table_rows_selected, ]
      taxon_data$datasets_selected <- data.frame(datasetKey = "gbif", count = sum(taxon_data$synonyms_selected$occurrence_count))
        
      shinyjs::hide(id = "taxon_options_panel")
      shinyjs::show(id = "load_data_panel")
      
      checkbox_choices <- paste0("gbif (", ifelse(taxon_data$datasets_selected$count < 5000, taxon_data$datasets_selected$count, 5000), ")")
      
      updateSelectizeInput(session = session,
                           "input_sources", 
                           choices = checkbox_choices, 
                           selected = checkbox_choices
      )
    }
    
  })
  
  observeEvent({
    input$select_datasets
    input$taxon_options_table_rows_selected
    }, {
      
    select_datasets_button_presses$values <- c(select_datasets_button_presses$values, input$select_datasets)
    
    if (select_datasets_button_presses$values[length(select_datasets_button_presses$values)] != select_datasets_button_presses$values[length(select_datasets_button_presses$values)-1]){
    
    # shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
      
    taxon_data$synonyms_selected <- taxon_data$synonyms[input$taxon_options_table_rows_selected, ]

    gbif_counts_occ <- purrr::map(taxon_data$synonyms_selected$key, function(k){
      out <- rgbif::occ_count(
        taxonKey = k, 
        basisOfRecord = "OCCURRENCE;PRESERVED_SPECIMEN;OBSERVATION;MACHINE_OBSERVATION", 
        gadmGid = taxon_data$assessment_polygon, 
        hasCoordinate = TRUE, facet = "datasetKey", facetLimit = 100)
      out
    }) %>% bind_rows()
      
      # rgbif::occ_count(taxonKey = taxon_keys, basisOfRecord = c("OCCURRENCE", "PRESERVED_SPECIMEN", "OBSERVATION", "MACHINE_OBSERVATION"), gadmGid = taxon_data$assessment_polygon, hasCoordinate = TRUE, facet = "datasetKey", facetLimit = 100)
    
    if (length(taxon_data$synonyms_selected$key) == 1){
      taxon_keys <- taxon_data$synonyms_selected$key
    } else {
      taxon_keys <- paste0(taxon_data$synonyms_selected$key, collapse = ";")
    }
    
    gbif_counts_humobs <- rgbif::occ_count(
        taxonKey = taxon_keys, 
        basisOfRecord = "HUMAN_OBSERVATION", 
        gadmGid = taxon_data$assessment_polygon, 
        hasCoordinate = TRUE, facet = "datasetKey", facetLimit = 100
        )
    
    gbif_counts_humobs <- rgbif::occ_count(taxonKey = taxon_keys, basisOfRecord = c("HUMAN_OBSERVATION"), gadmGid = taxon_data$assessment_polygon, hasCoordinate = TRUE, facet = "datasetKey", facetLimit = 100)
    
    gbif_counts <- rbind(
      gbif_counts_occ %>% dplyr::mutate(basisOfRecord = "OCCURRENCE"),
      gbif_counts_humobs %>% dplyr::mutate(basisOfRecord = "HUMAN_OBSERVATION")
    )
    
    gbif_counts$datasetName <- purrr::map_chr(gbif_counts$datasetKey, function(k) rgbif::dataset_get(k)$title)
    taxon_data$datasets <- gbif_counts %>% 
      # dplyr::mutate(recordsMax = count) %>% 
      dplyr::select(datasetName, datasetKey, count, basisOfRecord) #, recordsMax)

    shinyjs::hide(id = "taxon_options_panel")
    shinyjs::hide(id = "taxon_search_panel")
    shinyjs::show(id = "taxon_datasets_panel")
    
    # shinybusy::remove_modal_spinner()
    
    }
    
  })
  
  output$taxon_datasets_occ_table <- DT::renderDataTable({

    dat <- taxon_data$datasets %>% dplyr::filter(basisOfRecord == "OCCURRENCE")
    
    DT::dataTableProxy("taxon_datasets_occ_table") %>% selectRows(selected = 1:nrow((taxon_data$datasets %>% dplyr::filter(basisOfRecord == "OCCURRENCE"))))
    
    dat %>% 
      dplyr::mutate(datasetName = paste0("<a href='", paste0("https://www.gbif.org/dataset/", datasetKey), "' target='_blank'>", datasetName, "</a>")) %>%         
      dplyr::select(-datasetKey, -basisOfRecord) %>% 
      dplyr::rename("Dataset" = datasetName,
                    "Number of records available" = count# ,
                    # "Number of records to include" = recordsMax
      ) %>%
      DT::datatable(options = list(dom = 'tp', pageLength = 5, autoWidth = TRUE, language = list(emptyTable = 'There are no records of this type for this taxon')),
                    selection = list(mode = 'multiple', target = 'row', selected = 1:nrow(dat)),
                    escape = FALSE,
                    rownames = FALSE
      )
    
  })
  
  output$taxon_datasets_humobs_table <- DT::renderDataTable({
    
    dat <- taxon_data$datasets %>% dplyr::filter(basisOfRecord == "HUMAN_OBSERVATION")
    
    DT::dataTableProxy("taxon_datasets_humobs_table") %>% selectRows(selected = 1:nrow((taxon_data$datasets %>% dplyr::filter(basisOfRecord == "HUMAN_OBSERVATION"))))
    
    dat %>% 
      dplyr::mutate(datasetName = paste0("<a href='", paste0("https://www.gbif.org/dataset/", datasetKey), "' target='_blank'>", datasetName, "</a>")) %>%         
      dplyr::select(-datasetKey, -basisOfRecord) %>% 
      dplyr::rename("Dataset" = datasetName,
                    "Number of records available" = count# ,
                    # "Number of records to include" = recordsMax
      ) %>%
      DT::datatable(options = list(dom = 'tp', pageLength = 5, autoWidth = TRUE, language = list(emptyTable = 'There are no records of this type for this taxon')),
                    selection = list(mode = 'multiple', target = 'row', selected = 1:nrow(dat)),
                    escape = FALSE,
                    rownames = FALSE
      )
    
  })

  observeEvent(input$taxon_datasets_occ_table_rows_selected, {
    
    if (length(input$taxon_datasets_occ_table_rows_selected) < (taxon_data$datasets %>% dplyr::filter(basisOfRecord == "OCCURRENCE") %>% nrow())){
      
      updateCheckboxInput(session = session, "select_all_occ", value = FALSE)
      
    } else {
      
      updateCheckboxInput(session = session, "select_all_occ", value = TRUE)
      
    }
    
  }, ignoreInit = TRUE, ignoreNULL = FALSE)
  
  observeEvent(input$select_all_occ, {
    
    if (input$select_all_occ){
      DT::dataTableProxy("taxon_datasets_occ_table") %>% selectRows(selected = 1:nrow((taxon_data$datasets %>% dplyr::filter(basisOfRecord == "OCCURRENCE"))))
      updateCheckboxInput(session = session, "deselect_all_occ", value = FALSE)
    }
    
  })
  
  observeEvent(input$deselect_all_occ, {
    
    if (input$deselect_all_occ){
      DT::dataTableProxy("taxon_datasets_occ_table") %>% selectRows(selected = NULL)
      updateCheckboxInput(session = session, "select_all_occ", value = FALSE)
    }
    
  })

  observeEvent(input$taxon_datasets_humobs_table_rows_selected, {
    
    if (length(input$taxon_datasets_humobs_table_rows_selected) < (taxon_data$datasets %>% dplyr::filter(basisOfRecord == "HUMAN_OBSERVATION") %>% nrow())){
      
      updateCheckboxInput(session = session, "select_all_humobs", value = FALSE)
      
    } else {
      
      updateCheckboxInput(session = session, "select_all_humobs", value = TRUE)
      
    }
    
  }, ignoreInit = TRUE, ignoreNULL = FALSE)
  
  observeEvent(input$select_all_humobs, {
    
    if (input$select_all_humobs){
      DT::dataTableProxy("taxon_datasets_humobs_table") %>% selectRows(selected = 1:nrow((taxon_data$datasets %>% dplyr::filter(basisOfRecord == "HUMAN_OBSERVATION"))))
      updateCheckboxInput(session = session, "deselect_all_humobs", value = FALSE)
    }
    
  })
  
  observeEvent(input$deselect_all_humobs, {
    
    if (input$deselect_all_humobs){
      DT::dataTableProxy("taxon_datasets_humobs_table") %>% selectRows(selected = NULL)
      updateCheckboxInput(session = session, "select_all_humobs", value = FALSE)
    }
    
  })
  
  observeEvent(input$single_assessment_type, {
    
    if (input$single_assessment_type == "global"){
      updateSelectizeInput(session = session,
                           inputId = "single_assessment_nation",
                           choices = "",
                           selected = NULL
      )
      
      updateSelectizeInput(session = session,
                           inputId = "single_assessment_subnation",
                           choices = "",
                           selected = NULL
      )
    }
    
    if (input$single_assessment_type == "national"){
      updateSelectizeInput(session = session,
                           inputId = "single_assessment_nation",
                           choices = c(list("Canada" = "CA", "United States" = "US"))
      )
      
      updateSelectizeInput(session = session,
                           inputId = "single_assessment_subnation",
                           choices = "",
                           selected = NULL
      )
    
    }
    
    if (input$single_assessment_type == "subnational"){
      updateSelectizeInput(session = session,
                           inputId = "single_assessment_nation",
                           choices = c(list("Canada" = "CA", "United States" = "US"))
      )
      
    }

  })
  
  observeEvent(input$single_assessment_nation, {
    
    if (input$single_assessment_type == "subnational"){
      
    nation_subset <- network_polys %>% dplyr::filter(FIPS_CNTRY %in% input$single_assessment_nation)
    
    subnations_already_selected <- input$single_assessment_subnation
    
    updateSelectizeInput(session = session,
                         inputId = "single_assessment_subnation",
                         # choices = (nation_subset$Admin_abbr %>% na.omit() %>% as.character()) %>% set_names(nation_subset$ADMIN_NAME%>% na.omit() %>% as.character()) %>% sort(),
                         choices = nation_subset$ADMIN_NAME %>% na.omit() %>% as.character() %>% sort(),
                         selected = subnations_already_selected
    )
    
    }
    
    if (input$single_assessment_type == "national"){
    
    updateSelectizeInput(session = session,
                         inputId = "nation_filter",
                         selected = input$single_assessment_nation
    )
        
      # query_poly <- network_polys %>% dplyr::filter(
      #   FIPS_CNTRY %in% input$single_assessment_nation
      # )
      
      # query_poly_bbox <- sf::st_bbox(query_poly) %>% sf::st_as_sfc()
      # p <- terra::vect(query_poly_bbox)
      # pcc <- terra::forceCCW(p)
      # query_poly_bbox_wkt <- query_poly_bbox %>%
      #   terra::vect() %>%
      #   terra::forceCCW() %>%
      #   geom(wkt = TRUE)
      
      taxon_data$assessment_polygon <- ifelse(input$single_assessment_nation == "US", "USA", ifelse(input$single_assessment_nation == "CA", "CAN", NULL))
      
      } 
    
    
  })
  
  observeEvent(input$single_assessment_subnation, {
    
    relevant_nation <- network_polys %>% dplyr::filter(ADMIN_NAME %in% input$single_assessment_subnation) %>% dplyr::pull(FIPS_CNTRY)
    
    updateSelectizeInput(session = session,
                         inputId = "single_assessment_nation", 
                         selected = relevant_nation %>% set_names(relevant_nation)
    )
    
    nation_subset <- network_polys %>% dplyr::filter(FIPS_CNTRY %in% input$single_assessment_nation)
    
    subnations_already_selected <- input$single_assessment_subnation
    
    # updateSelectizeInput(session = session,
    #                      inputId = "single_assessment_subnation",
    #                      choices = (nation_subset$Admin_abbr %>% na.omit() %>% as.character()) %>% set_names(nation_subset$ADMIN_NAME%>% na.omit() %>% as.character()) %>% sort(),
    #                      selected = subnations_already_selected
    # )
    
    updateSelectizeInput(session = session,
                         inputId = "single_assessment_subnation",
                         choices = nation_subset$ADMIN_NAME %>% na.omit() %>% as.character() %>% sort(),
                         selected = subnations_already_selected
    )
    
    updateSelectizeInput(session = session,
                         inputId = "states_filter", 
                         selected = input$single_assessment_subnation
    )
    
    if (!is.null(input$single_assessment_subnation)){
      
    # query_poly <- network_polys %>% dplyr::filter(
    #   FIPS_CNTRY %in% input$single_assessment_nation,
    #   Admin_abbr %in% input$single_assessment_subnation
    # )
    # 
    # query_poly_bbox <- sf::st_bbox(query_poly) %>% sf::st_as_sfc()
    # p <- terra::vect(query_poly_bbox)
    # pcc <- terra::forceCCW(p)
    # query_poly_bbox_wkt <- query_poly_bbox %>%
    #   terra::vect() %>%
    #   terra::forceCCW() %>%
    #   geom(wkt = TRUE)
    
    taxon_data$assessment_polygon <- gadm_df %>% 
      dplyr::filter(NAME_1 %in% gsub(" ", "", input$single_assessment_subnation)) %>% 
      dplyr::pull(GID_1) %>% 
      unique()
    
    } else {
      taxon_data$assessment_polygon <- NULL
    }
    
  })

  observeEvent({
    input$begin_assessment_detailed
  }, {
    
    assessment_start_detailed_button_presses$values <- c(assessment_start_detailed_button_presses$values, input$begin_assessment_detailed)
    
    if (assessment_start_detailed_button_presses$values[length(assessment_start_detailed_button_presses$values)] != assessment_start_detailed_button_presses$values[length(assessment_start_detailed_button_presses$values)-1]){
      
      if (!is.null(input$taxon_datasets_occ_table_rows_selected)){
        taxon_data$datasets_selected <- rbind(
          taxon_data$datasets_selected,
          taxon_data$datasets %>% 
            dplyr::filter(basisOfRecord == "OCCURRENCE") %>% 
            dplyr::slice(input$taxon_datasets_occ_table_rows_selected) %>% 
            dplyr::mutate(recordsMax = ifelse(count < as.numeric(input$number_occ), count, as.numeric(input$number_occ)) %>% as.numeric())
        ) 
      } 
      
      if (!is.null(input$taxon_datasets_humobs_table_rows_selected)){
        taxon_data$datasets_selected <- rbind(
          taxon_data$datasets_selected,
          taxon_data$datasets %>% 
            dplyr::filter(basisOfRecord == "HUMAN_OBSERVATION") %>% 
            dplyr::slice(input$taxon_datasets_humobs_table_rows_selected) %>% 
            dplyr::mutate(recordsMax = ifelse(count < as.numeric(input$number_humobs), count, as.numeric(input$number_humobs)) %>% as.numeric())
        ) 
      } 
      
      shinyjs::hide(id = "taxon_options_panel")
      shinyjs::hide(id = "taxon_search_panel")
      shinyjs::hide(id = "taxon_datasets_panel")
      shinyjs::show(id = "load_data_panel")
      
      datasets_grouped <- taxon_data$datasets_selected %>% 
        dplyr::group_by(datasetName) %>% 
        dplyr::summarise(recordsMax = sum(recordsMax))
      
      checkbox_choices <- paste0(datasets_grouped$datasetName, " (", 
                                 datasets_grouped$recordsMax, ")"
                                 )
      
      updateSelectizeInput(session = session,
                           "input_sources", 
                           choices = checkbox_choices, 
                           selected = checkbox_choices
      )
    }
    
  })
  
  observeEvent({
    input$load_gbif_data
    input$clean_occ
    input$centroid_filter
    }, {
    
    if (isTRUE(input$load_gbif_data)){
      
      if (is.null(taxon_data$gbif_occurrences_raw)){
        
        # shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
        
        gbif_download <- download_gbif_data()
        
        taxon_data$gbif_occurrences_raw <- gbif_download$sp_occurrences
        taxon_data$shifted <- gbif_download$shifted
        
        uploaded_long_max <- ifelse(!is.null(taxon_data$uploaded_occurrences), max(taxon_data$uploaded_occurrences$longitude, na.rm = TRUE), NA)
        uploaded_long_min <- ifelse(!is.null(taxon_data$uploaded_occurrences), min(taxon_data$uploaded_occurrences$longitude, na.rm = TRUE), NA)
        
        long_range <- abs(max(c(max(taxon_data$gbif_occurrences_raw$longitude, na.rm = TRUE), uploaded_long_max), na.rm = TRUE)) +
          abs(min(c(min(taxon_data$gbif_occurrences_raw$longitude, na.rm = TRUE), uploaded_long_min), na.rm = TRUE))
        
        if (long_range > 360){
          taxon_data$gbif_occurrences_raw$longitude[taxon_data$gbif_occurrences_raw$longitude < 0] <- taxon_data$gbif_occurrences_raw$longitude[taxon_data$gbif_occurrences_raw$longitude < 0] + 360
        }
      }

      taxon_data$gbif_occurrences <- taxon_data$gbif_occurrences_raw %>% 
        clean_gbif_data(clean = input$clean_occ, remove_centroids = input$centroid_filter, minimum_fields = minimum_fields) %>% 
        dplyr::mutate(scientificName_Assessment = taxon_data$info$scientificName)

      print(names(taxon_data$gbif_occurrences))
      
      shinyjs::show(id = "data_panel")
      
      # shinybusy::remove_modal_spinner()
      
    } 
    
  })
  
  observeEvent({
    input$filedata
    input$map_uploads
  }, {
    
    if (isTRUE(input$map_uploads)){
      
      if (!is.null(input$filedata)){
        
        # shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
        
        taxon_data$uploaded_occurrences <- uploaded_data()
        # taxon_data$uploaded_occurrences <- taxon_data$uploaded_occurrences %>% 
        #   dplyr::mutate(key = paste(prov, 1:nrow(taxon_data$uploaded_occurrences), sep = "_"))
        
        if (length(is.na(taxon_data$uploaded_occurrences$key)) > 0){
          taxon_data$uploaded_occurrences$key[is.na(taxon_data$uploaded_occurrences$key)] <- paste(taxon_data$uploaded_occurrences$prov, 1:length(taxon_data$uploaded_occurrences$key[is.na(taxon_data$uploaded_occurrences$key)]), sep = "_")
        }
        
        # taxon_data$uploaded_occurrences$longitude[taxon_data$uploaded_occurrences$longitude > 180] <- taxon_data$uploaded_occurrences$longitude[taxon_data$uploaded_occurrences$longitude > 180] - 360
        taxon_data$uploaded_occurrences$longitude[taxon_data$uploaded_occurrences$longitude > 180] <- taxon_data$uploaded_occurrences$longitude[taxon_data$uploaded_occurrences$longitude > 180] - 360
        #
        max_long <- max(taxon_data$uploaded_occurrences$longitude, na.rm = TRUE)/2
        shifted_long <- taxon_data$uploaded_occurrences$longitude

        if (length(taxon_data$uploaded_occurrences$longitude[taxon_data$uploaded_occurrences$longitude > max_long]) > 0){
          shifted_long[shifted_long > max_long] <- shifted_long[shifted_long > max_long] - 360
          shifted_long <- shifted_long + 360
        }

        if ((max(shifted_long)-min(shifted_long)) < (max(taxon_data$uploaded_occurrences$longitude) - min(taxon_data$uploaded_occurrences$longitude))){
          taxon_data$uploaded_occurrences$longitude <- shifted_long
          taxon_data$shifted <- TRUE
        }
        
        gbif_long_max <- ifelse(!is.null(taxon_data$gbif_occurrences), max(taxon_data$gbif_occurrences$longitude, na.rm = TRUE), NA)
        gbif_long_min <- ifelse(!is.null(taxon_data$gbif_occurrences), min(taxon_data$gbif_occurrences$longitude, na.rm = TRUE), NA)
        
        long_range <- abs(max(c(max(taxon_data$uploaded_occurrences$longitude, na.rm = TRUE), gbif_long_max), na.rm = TRUE)) +
                      abs(min(c(min(taxon_data$uploaded_occurrences$longitude, na.rm = TRUE), gbif_long_min), na.rm = TRUE))
        
        if (long_range > 360){
          taxon_data$uploaded_occurrences$longitude[taxon_data$uploaded_occurrences$longitude < 0] <- taxon_data$uploaded_occurrences$longitude[taxon_data$uploaded_occurrences$longitude < 0] + 360
        }

        selected_taxon$NS <- c(selected_taxon$NS, unique(taxon_data$uploaded_occurrences$scientificName)) %>% unique()
        
        if (taxon_data$info$scientificName == "New taxon"){
          taxon_data$info <- data.frame(
            scientificName = ifelse(("scientificName_Assessment" %in% names(taxon_data$uploaded_occurrences)), taxon_data$uploaded_occurrences$scientificName_Assessment[1], taxon_data$uploaded_occurrences$scientificName[1]),
            elementGlobalId = NA, Source = "Uploaded", primaryCommonName = NA, roundedGRank = NA, elcode = NA, synonyms = NA, uniqueId = NA
          )
        }
        
        shinyjs::show(id = "data_panel")
        
        # shinybusy::remove_modal_spinner()
        
      }
    }
    
  })
  
  observeEvent(input$main_map_draw_new_feature, {
    
      if (input$main_map_draw_new_feature$geometry$type == "Polygon" & !is.null(taxon_data$sf_filtered)){
        
        # shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
        
        drawn_shape_coordinates <- input$main_map_draw_new_feature$geometry$coordinates[[1]]
        
        pol <- st_polygon(
          list(
            cbind(
              purrr::map(drawn_shape_coordinates, 1) %>% unlist() %>% as.numeric(),
              purrr::map(drawn_shape_coordinates, 2) %>% unlist() %>% as.numeric()
            )
          )
        )
        
        pol_bbox <- st_bbox(pol)
        
        taxon_data$selected_points <- taxon_data$sf_filtered[which(purrr::map_int(st_intersects(taxon_data$sf_filtered, pol), length) > 0), ]
        
        leafletProxy("main_map") %>%
          flyToBounds(pol_bbox[[1]], pol_bbox[[2]], pol_bbox[[3]], pol_bbox[[4]], options = list(animate = TRUE, duration = 1, easeLinearity = 0.1, noMoveStart = TRUE)) %>% 
          addCircleMarkers( 
            data = taxon_data$selected_points,
            lng = ~longitude,
            lat = ~latitude,
            layerId = ~key,
            fillColor = "#4169E1",
            fillOpacity = 0.75,
            color = "#4169E1"
          )
        
        clicks$IDs <- c(clicks$IDs, taxon_data$selected_points$key)
        
        # shinybusy::remove_modal_spinner()
        
      }
      
      if (input$main_map_draw_new_feature$geometry$type == "Point"){
        
        # shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
        
        drawn_point <- data.frame(longitude = input$main_map_draw_new_feature$geometry$coordinates[[1]],
                                  latitude = input$main_map_draw_new_feature$geometry$coordinates[[2]],
                                  year = substr(Sys.Date(), 1, 4),
                                  scientificName = taxon_data$info$scientificName
        )
        
        taxon_data$drawn_occurrences <- list(taxon_data$drawn_occurrences, drawn_point) %>% bind_rows()
        
        taxon_data$drawn_occurrences <- taxon_data$drawn_occurrences %>%
          cbind(matrix(NA, nrow = nrow(taxon_data$drawn_occurrences), ncol = length(setdiff(minimum_fields, names(taxon_data$drawn_occurrences)))) %>%
                  as.data.frame() %>%
                  set_names(setdiff(minimum_fields, names(taxon_data$drawn_occurrences)))
          ) %>%
          dplyr::mutate(prov = "drawn",
                        key = paste("user_created", 1:nrow(taxon_data$drawn_occurrences), sep = "_")
          ) %>%
          dplyr::select(all_of(minimum_fields))

        shinyjs::show(id = "data_panel")
        
        # shinybusy::remove_modal_spinner()
        
      }
    
  })
  
  observe({
    
    ### Combine all occurrences
    taxon_data$all_occurrences <- rbind(
      taxon_data$gbif_occurrences,
      taxon_data$uploaded_occurrences,
      taxon_data$drawn_occurrences
    )
    
    if (!is.null(taxon_data$all_occurrences)){

      ### Create simple features object for geospatial calculations
      taxon_data$sf <- taxon_data$all_occurrences %>% 
        dplyr::filter(complete.cases(longitude, latitude)) %>% 
        dplyr::mutate(lon = longitude,
                      lat = latitude) %>% 
        sf::st_as_sf(
          coords = c("lon", "lat"),
          crs = 4326
        )
      
      selected_taxon$NS <- c(selected_taxon$NS, unique(taxon_data$sf$scientificName)) %>% unique()
      
      taxon_data$nations <- network_polys[which(purrr::map_int(st_intersects(network_polys, taxon_data$sf), length) > 0), ]$FIPS_CNTRY %>% na.omit() %>% as.character()
      taxon_data$states <- network_polys[which(purrr::map_int(st_intersects(network_polys, taxon_data$sf), length) > 0), ]$Admin_abbr %>% na.omit() %>% as.character()

      shinyjs::show("species_occurrences_table")
      updateCollapse(session = session, id = "inputs_single", close = "Add assessment data")
      
      m <- leafletProxy("main_map") %>%
        clearShapes() %>% 
        clearMarkers() %>% 
        clearMarkerClusters() 
      
      if (nrow(taxon_data$sf) >= 2){

        if (is.null(taxon_data$drawn_occurrences)){
          
        species_occurrences_bbox <- sf::st_bbox(taxon_data$sf)
        
        m <- m %>%
          flyToBounds(species_occurrences_bbox[[1]],
                      species_occurrences_bbox[[2]]-0.25*(species_occurrences_bbox[[4]] - species_occurrences_bbox[[2]]),
                      species_occurrences_bbox[[3]]+0.75*(species_occurrences_bbox[[3]] - species_occurrences_bbox[[1]]),
                      species_occurrences_bbox[[4]]+0.25*(species_occurrences_bbox[[4]] - species_occurrences_bbox[[2]]),
                      options = list(animate = TRUE, duration = 1, easeLinearity = 0.1, noMoveStart = TRUE)
          )
        
        m
        
        }
        
      }
      
      updateDateRangeInput(session = session, inputId = "year_filter",
                           start = paste0(min(unique(taxon_data$sf$year), na.rm = TRUE), "-01-01"),
                           end = paste0(max(unique(taxon_data$sf$year), na.rm = TRUE), "-01-01")
      )
      
      updateMaterialSwitch(session = session, inputId = "range_extent", value = FALSE)
      updateMaterialSwitch(session = session, inputId = "area_of_occupancy", value = FALSE)
      updateMaterialSwitch(session = session, inputId = "number_EOs", value = FALSE)
      # updateSelectizeInput(session = session, inputId = "nation_filter", choices = ifelse(is.null(input$single_assessment_nation), taxon_data$nations, input$single_assessment_nation), selected = ifelse(is.null(input$single_assessment_nation), NULL, input$single_assessment_nation))
      # updateSelectizeInput(session = session, inputId = "states_filter", choices = ifelse(is.null(input$single_assessment_subnation), taxon_data$states, input$single_assessment_subnation), selected = ifelse(is.null(input$single_assessment_subnation), NULL, input$single_assessment_subnation)) # , taxon_data$states)
      updateSelectizeInput(session = session, inputId = "nation_filter", choices = taxon_data$nations, selected = NULL)
      updateSelectizeInput(session = session, inputId = "states_filter", choices = taxon_data$states, selected = NULL) # , taxon_data$states)
      
      updateSelectizeInput(session = session, inputId = "synonyms_filter", choices = unique(taxon_data$sf$scientificName), selected = unique(taxon_data$sf$scientificName))
      
      if (!identical(NA, unique(taxon_data$sf$basisOfRecord))){
        updateSelectizeInput(session = session, inputId = "type_filter", choices = unique(taxon_data$sf$basisOfRecord), selected = unique(taxon_data$sf$basisOfRecord))
      }
      if (!identical(NA, unique(taxon_data$sf$EORANK))){
        updateSelectizeInput(session = session, inputId = "rank_filter", choices = setdiff(unique(taxon_data$sf$EORANK), NA), selected = setdiff(unique(taxon_data$sf$EORANK), NA))
      }
      if (!identical(NA, unique(taxon_data$sf$datasetName))){
      updateSelectizeInput(session = session, inputId = "sources_filter", choices = unique(taxon_data$sf$datasetName) %>% sort(), selected = unique(taxon_data$sf$datasetName) %>% sort())
      }
      updateDateRangeInput(session = session, inputId = "year_filter",
                           start = paste0(min(unique(taxon_data$sf$year), na.rm = TRUE), "-01-01"),
                           end = paste0(max(unique(taxon_data$sf$year), na.rm = TRUE), "-01-01")
      )
      
      if (!is.null(batch_taxon_focus$taxon)){
        batch_taxon_filters <- batch_run_output$results[[batch_taxon_focus$taxon]]$filters_selected

        updateMaterialSwitch(session = session, inputId = "clean_occ", value = batch_taxon_filters$clean_occ)
        updateMaterialSwitch(session = session, inputId = "centroid_filter", value = batch_taxon_filters$centroid_filter)

        if (!is.null(batch_taxon_filters$nations_filter)){
          updateSelectizeInput(session = session, inputId = "nation_filter", selected = batch_taxon_filters$nations_filter)
        }
        if (!is.null(batch_taxon_filters$states_filter)){
          updateSelectizeInput(session = session, inputId = "states_filter", selected = batch_taxon_filters$states_filter) # , taxon_data$states)
        }
        if (!is.null(batch_taxon_filters$sources_filter)){
          if (!identical(NA, unique(batch_run_output$results[[batch_taxon_focus$taxon]]$sf_filtered$datasetName))){
          updateSelectizeInput(session = session, inputId = "sources_filter", choices = unique(batch_run_output$results[[batch_taxon_focus$taxon]]$sf_filtered$datasetName), selected = unique(batch_run_output$results[[batch_taxon_focus$taxon]]$sf_filtered$datasetName))
          }
        }
        if (!is.null(batch_taxon_filters$type_filter)){
          if (!identical(NA, unique(batch_run_output$results[[batch_taxon_focus$taxon]]$sf_filtered$basisOfRecord))){
          updateSelectizeInput(session = session, inputId = "type_filter", choices = unique(batch_run_output$results[[batch_taxon_focus$taxon]]$sf_filtered$basisOfRecord), selected = unique(batch_run_output$results[[batch_taxon_focus$taxon]]$sf_filtered$basisOfRecord))
          }
        }
        if (!is.null(batch_taxon_filters$rank_filter)){
          if (!identical(NA, unique(batch_run_output$results[[batch_taxon_focus$taxon]]$sf_filtered$EORANK))){
            updateSelectizeInput(session = session, inputId = "rank_filter", choices = unique(batch_run_output$results[[batch_taxon_focus$taxon]]$sf_filtered$EORANK), selected = unique(batch_run_output$results[[batch_taxon_focus$taxon]]$sf_filtered$EORANK))
          }
        }
        if (batch_taxon_filters$uncertainty_filter != ""){
          updateTextInput(session = session, inputId = "uncertainty_filter", selected = batch_taxon_filters$uncertainty_filter)
        }
        updateDateRangeInput(session = session, inputId = "year_filter",
                             start = batch_taxon_filters$date_start,
                             end = batch_taxon_filters$date_end
        )
        updateSelectizeInput(session = session,
                             inputId = "seasonality",
                             selected = batch_taxon_filters$months
        )
        updateSelectInput(session = session, inputId = "grid_cell_size", selected = batch_taxon_filters$grid_cell_size)
        updateTextInput(session = session, inputId = "separation_distance", value = batch_taxon_filters$separation_distance)
        updateMaterialSwitch(session = session, inputId = "range_extent", value = TRUE)
        updateMaterialSwitch(session = session, inputId = "area_of_occupancy", value = TRUE)
        updateMaterialSwitch(session = session, inputId = "number_EOs", value = TRUE)
        updateDateRangeInput(session = session, inputId = "period1", start = input$batch_period1[1], end = input$batch_period1[2])
        updateDateRangeInput(session = session, inputId = "period2", start = input$batch_period2[1], end = input$batch_period2[2])
        updateDateRangeInput(session = session, inputId = "period3", start = NA, end = NA)

      }
      
      observeEvent({
        input$year_filter
        input$no_year
        input$seasonality
        input$uncertainty_filter
        input$nation_filter
        input$states_filter
        input$synonyms_filter
        input$sources_filter
        input$type_filter
        input$rank_filter
        input$remove_selections
      }, {
        
        # shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
        
        if (!is.null(taxon_data$sf)){
          
          taxon_data$sf_filtered <- taxon_data$sf
          
          long_range <- abs(max(taxon_data$sf_filtered$longitude, na.rm = TRUE))+abs(min(taxon_data$sf_filtered$longitude, na.rm = TRUE))
          
          if (long_range > 360){
            taxon_data$sf_filtered$longitude[taxon_data$sf_filtered$longitude < 0] <- taxon_data$sf_filtered$longitude[taxon_data$sf_filtered$longitude < 0] + 360
          }
          
          # if (length(taxon_data$sf_filtered$longitude > 180) > 0){
            
          # taxon_data$sf_filtered$longitude[taxon_data$sf_filtered$longitude > 180] <- taxon_data$sf_filtered$longitude[taxon_data$sf_filtered$longitude > 180] - 360
          # 
          # max_long <- max(taxon_data$sf_filtered$longitude, na.rm = TRUE)/2
          # shifted_long <- taxon_data$sf_filtered$longitude
          # 
          # if (length(taxon_data$sf_filtered$longitude[taxon_data$sf_filtered$longitude > max_long]) > 0){
          #   shifted_long[shifted_long > max_long] <- shifted_long[shifted_long > max_long] - 360
          #   shifted_long <- shifted_long + 360
          # }
          # 
          # if ((max(shifted_long)-min(shifted_long)) < (max(taxon_data$sf_filtered$longitude) - min(taxon_data$sf_filtered$longitude))){
          #   taxon_data$sf_filtered$longitude <- shifted_long
          #   taxon_data$shifted <- TRUE
          # }
          # 
          
          # }
          
          months <- purrr::map_dbl(1:12, function(x) x) %>% purrr::set_names(substr(month.name, 1, 3))
          season <- which(names(months) %in% input$seasonality)
          # season <- ifelse(nchar(season) == 1, paste0("0", season), season)
            
          taxon_data$sf_filtered <- taxon_data$sf_filtered %>%
            dplyr::filter(year >= substr(input$year_filter[1], 1, 4) & year <= substr(input$year_filter[2], 1, 4) | is.na(year),
                          month %in% season | is.na(month),
                          key %in% setdiff(key, taxon_data$removed_points$key)
            )
          
          if (input$uncertainty_filter != ""){
            
            taxon_data$sf_filtered <- taxon_data$sf_filtered %>%
              dplyr::filter(coordinateUncertaintyInMeters <= as.numeric(input$uncertainty_filter) | is.na(coordinateUncertaintyInMeters))
          }
          
          if (!is.null(input$nation_filter)){
            taxon_data$sf_filtered <- taxon_data$sf_filtered %>%
              dplyr::filter(purrr::map_int(st_intersects(taxon_data$sf_filtered, network_polys %>% dplyr::filter(FIPS_CNTRY %in% input$nation_filter)), length) > 0)
          }
          
          if (!is.null(input$states_filter)){
            taxon_data$sf_filtered <- taxon_data$sf_filtered %>%
              dplyr::filter(purrr::map_int(st_intersects(taxon_data$sf_filtered, network_polys %>% dplyr::filter(Admin_abbr %in% input$states_filter)), length) > 0)
          }
          
          if (!is.null(input$rank_filter)){
            taxon_data$sf_filtered <- taxon_data$sf_filtered %>%
              dplyr::filter(EORANK %in% input$rank_filter | is.na(EORANK))
          }
          
          if (isTRUE(input$remove_selections)){
            
            taxon_data$removed_points <- rbind(taxon_data$removed_points, taxon_data$selected_points)
            
          }
          
          name_exclusions <- setdiff(taxon_data$sf_filtered$scientificName, input$synonyms_filter)
          if (length(name_exclusions) > 0){
            taxon_data$sf_filtered <- taxon_data$sf_filtered %>%
              dplyr::filter(!(scientificName %in% name_exclusions))
          }
          
          type_exclusions <- setdiff(unique(taxon_data$sf_filtered$basisOfRecord), input$type_filter)
          
          if (length(type_exclusions) > 0){
            if (!identical(NA, type_exclusions)){
            taxon_data$sf_filtered <- taxon_data$sf_filtered %>%
              dplyr::filter(!(basisOfRecord %in% type_exclusions) | is.na(basisOfRecord))
             updateSelectizeInput(session = session, inputId = "sources_filter", selected = unique(taxon_data$sf_filtered$datasetName) %>% sort())
            }
          }
          #  } else {
          #    updateSelectizeInput(session = session, inputId = "sources_filter", selected = unique(taxon_data$sf$datasetName) %>% sort())
          # }
          
          # 
          # if (!is.null(taxon_data$datasets_selected)){
          # if ("gbif" %in% taxon_data$datasets$datasetKey){
          #   source_exclusions <- setdiff(taxon_data$sf_filtered$prov, input$sources_filter)
          #   if (length(source_exclusions) > 0){
          #     taxon_data$sf_filtered <- taxon_data$sf_filtered %>%
          #       dplyr::filter(!(prov %in% source_exclusions) | is.na(prov))
          #   }
          # } else {
          
          source_exclusions <- setdiff(unique(taxon_data$sf_filtered$datasetName), input$sources_filter)
          
          if (length(source_exclusions) > 0){
            if (!identical(NA, source_exclusions)){
              taxon_data$sf_filtered <- taxon_data$sf_filtered %>%
                dplyr::filter(!(datasetName %in% source_exclusions) | is.na(datasetName))
              # updateSelectizeInput(session = session, inputId = "sources_filter", selected = unique(taxon_data$sf_filtered$datasetName) %>% sort())
            }
          }
          
          #   source_exclusions <- setdiff(taxon_data$datasets_selected$datasetName, input$sources_filter)
          #   if (length(source_exclusions) > 0){
          #     source_exclusions <- taxon_data$datasets_selected %>%
          #       dplyr::filter(datasetName %in% source_exclusions) %>%
          #       dplyr::pull(datasetKey)
          #     records_to_exclude <- taxon_data$gbif_occurrences_raw %>%
          #       dplyr::filter(datasetKey %in% source_exclusions) %>%
          #       dplyr::pull(key)
          #     taxon_data$sf_filtered <- taxon_data$sf_filtered %>%
          #       dplyr::filter(!(key %in% records_to_exclude))
          #   }
          # }
          # }

          m <- leafletProxy("main_map") %>%
            clearMarkers() %>%
            clearMarkerClusters() %>%
            leaflet::addMapPane("species_observations", zIndex = 300) %>%
            addCircleMarkers(
              data = taxon_data$sf_filtered,
              lng = ~longitude,
              lat = ~latitude,
              clusterOptions = markerClusterOptions(),
              layerId = ~key,
              fillColor = "#4169E1",
              fillOpacity = 0.25,
              color = grey(0.25),
              options = pathOptions(pane = "species_observations")
            )
          
          if (!is.null(taxon_data$selected_points)){
            
            taxon_data$selected_points <- taxon_data$sf_filtered %>%
              dplyr::filter(key %in% taxon_data$selected_points$key)
            
            m <- m %>%
              clearMarkerClusters() %>%
              clearMarkers() %>%
              leaflet::addMapPane("selected_observations", zIndex = 400) %>%
              addCircleMarkers(
                data = taxon_data$selected_points,
                lng = ~longitude,
                lat = ~latitude,
                layerId = ~key,
                fillColor = "#4169E1",
                fillOpacity = 0.75,
                color = "#4169E1",
                options = pathOptions(pane = "selected_observations")
              ) %>%
              addCircleMarkers(
                data = taxon_data$sf_filtered %>% dplyr::filter(key %in% setdiff(key, taxon_data$selected_points$key)),
                lng = ~longitude,
                lat = ~latitude,
                # clusterOptions = markerClusterOptions(),
                layerId = ~key,
                fillColor = "#4169E1",
                fillOpacity = 0.25,
                color = grey(0.25),
                options = pathOptions(pane = "species_observations"),
                group = "Point observations"
              )
            
          }
          
          if (isTRUE(input$range_extent)){
            updateMaterialSwitch(session = session, inputId = "range_extent", value = FALSE)
            updateMaterialSwitch(session = session, inputId = "range_extent", value = TRUE)
          }
          
          if (isTRUE(input$area_of_occupancy)){
            updateMaterialSwitch(session = session, inputId = "area_of_occupancy", value = FALSE)
            updateMaterialSwitch(session = session, inputId = "area_of_occupancy", value = TRUE)
          }
          
          if (isTRUE(input$number_EOs)){
            updateMaterialSwitch(session = session, inputId = "number_EOs", value = FALSE)
            updateMaterialSwitch(session = session, inputId = "number_EOs", value = TRUE)
          }
          
          if (isTRUE(input$remove_selections)){
            shinyjs::delay(
              1000,
              updateMaterialSwitch(session = session, inputId = "remove_selections", value = FALSE)
            )
          }
          
          m
          
          # shinybusy::remove_modal_spinner() # remove the modal window
          
        }
        
      }) 
      
    }
  })
 
  observeEvent(input$single_assessment_subnation, {
    
    if (!is.null(input$single_assessment_subnation)){
      
      subnation_bbox <- network_polys %>% 
        dplyr::filter(Admin_abbr %in% input$single_assessment_subnation) %>% 
        sf::st_bbox()
      
      m <- leafletProxy("main_map") %>%
        flyToBounds(subnation_bbox[[1]],
                    subnation_bbox[[2]]-0.25*(subnation_bbox[[4]] - subnation_bbox[[2]]),
                    subnation_bbox[[3]]+0.75*(subnation_bbox[[3]] - subnation_bbox[[1]]),
                    subnation_bbox[[4]]+0.25*(subnation_bbox[[4]] - subnation_bbox[[2]]),
                    options = list(animate = TRUE, duration = 1, easeLinearity = 0.1, noMoveStart = TRUE)
        )
      
      m
      
    }
    
  })
  
  observeEvent(input$main_map_marker_click, {
    
    # shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
    
    clicks$IDs <- c(clicks$IDs, input$main_map_marker_click$id)
    
    if (sum(clicks$IDs %in% input$main_map_marker_click$id) > 1) clicks$IDs <- setdiff(clicks$IDs, input$main_map_marker_click$id)
    
    taxon_data$selected_points <- taxon_data$sf_filtered %>% 
      dplyr::filter(key %in% clicks$IDs)
    
    leafletProxy("main_map") %>%
      clearMarkerClusters() %>% 
      clearMarkers() %>% 
      leaflet::addMapPane("selected_observations", zIndex = 400) %>% 
      addCircleMarkers( 
        data = taxon_data$selected_points,
        lng = ~longitude,
        lat = ~latitude,
        layerId = ~key,
        fillColor = "#4169E1",
        fillOpacity = 0.75,
        color = "#4169E1",
        options = pathOptions(pane = "selected_observations")
      ) %>% 
      addCircleMarkers( 
        data = taxon_data$sf_filtered %>% dplyr::filter(key %in% setdiff(key, taxon_data$selected_points$key)),
        lng = ~longitude,
        lat = ~latitude,
        layerId = ~key,
        fillColor = "#4169E1",
        fillOpacity = 0.25,            
        color = grey(0.25),
        options = pathOptions(pane = "species_observations")
      )
    
    # shinybusy::remove_modal_spinner() # remove the modal window
    
  })
  
  output$occurrences_barchart_full <- dygraphs::renderDygraph({
    
    dat <- taxon_data$sf_filtered %>%
      dplyr::filter(complete.cases(year)) %>% 
      dplyr::group_by(year) %>%
      dplyr::summarise(number_records = n()) %>%
      dplyr::mutate(year = paste0(year, "-01-01") %>% as.Date())
    
    taxon_data$records_over_time <- xts::xts(x = dat$number_records, order.by = dat$year)
    
    start_window <- input$year_filter[1]
    end_window <- input$year_filter[2]
    
    dygraph(taxon_data$records_over_time, ylab = "") %>%
      dyBarChart() %>%
      dySeries("V1", label = "Number of records", color = "#1f417d") %>%
      dyAxis(
        "y",
        axisLabelWidth = 0
      ) %>% 
      dyRangeSelector(dateWindow = c(start_window %>% as.Date(), end_window %>% as.Date()))
    
  })
  
  
  output$metric_barchart_period <- plotly::renderPlotly({
    
    dat <- calculate_rarity_change(taxon_data = taxon_data, 
                                   period1 = input$period1,
                                   period2 = input$period2,
                                   period3 = input$period3,
                                   aoo_grid_cell_size = input$grid_cell_size, 
                                   occ_sep_distance = input$separation_distance
                                   ) %>% 
      dplyr::filter(complete.cases(.))
    
    taxon_data$temporal_change <- dat 
    
    p <- ggplot(data = dat) + 
      geom_bar(mapping = aes(x = .data[["period"]], y = .data[[input$barchart_metric]]), stat = "identity", fill = "#B2BBCB") +
      theme_bw() +
      theme(panel.background = element_rect(fill = "#F9F9F9")) +
      xlab("Time period") +
      ylab(
        which(c("Number of Records" = "rec_count", "Range Extent" = "eoo", "Area of Occupancy" = "aoo", "Number of Occurrences" = "eo_count") == input$barchart_metric) %>% names()
      ) +
      xlim(1-0.5, nrow(dat)+0.5)

      p <- p + annotate("text", x = dat$period[2], y = max(dat[[input$barchart_metric]], na.rm = TRUE)/2, label = paste0(dat[[paste0(input$barchart_metric, "_change")]][2], "%"), parse = TRUE, colour = grey(.1)) 
      p <- p + annotate("text", x = dat$period[3], y = max(dat[[input$barchart_metric]], na.rm = TRUE)/2, label = paste0(dat[[paste0(input$barchart_metric, "_change")]][3], "%"), parse = TRUE, colour = grey(.1))
    
    gg <- plotly_build(p) %>%
      config(displayModeBar = FALSE) %>%
      layout(plot_bgcolor  = "rgba(0, 0, 0, 0)",
             paper_bgcolor = "rgba(0, 0, 0, 0)",
             font = list(family = "Helvetica", size = 14), 
             xaxis = list(titlefont = list(size = 13),
                          tickfont = list(size = 13)),
             yaxis = list(titlefont = list(size = 13),
                          tickfont = list(size = 13))
      )
    
    gg 
    
  })
 
  observeEvent(input$temporal_trend, {
    
    shinyjs::show("temporal_trend_plots")
    
    # shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
    
    temporal_trends_output <- purrr::safely(get_temporal_trends)(
      taxon_data = taxon_data,
      referenceTaxon = input$select_reference_taxon,
      start_year = input$select_start_year
    )

    if (!is.null(temporal_trends_output$result)){
      output$temporal_trends_output <- plotly::renderPlotly(
        temporal_trends_output$result
      )
    } else {
      sendSweetAlert(session, type = "warning", title = "Oops!", text = "Temporal trend could not be calculated; try a different reference taxon or start year", closeOnClickOutside = TRUE)
    }

    # shinybusy::remove_modal_spinner()
  })
  
  output$occurrences_table <- DT::renderDataTable({
    
    dat <- taxon_data$selected_points

    if ("sf" %in% class(dat)){
      dat <- dat %>%
        sf::st_set_geometry(NULL)
    }

    dat <- dat %>%
      dplyr::select(all_of(minimum_fields)) %>%
      dplyr::mutate(
        key = paste0("<a href='", references, "' target='_blank'>", key, "</a>"),
        references = paste0("<a href='", references, "' target='_blank'>", references, "</a>"),
        latitude = round(latitude, 4),
        longitude = round(longitude, 4)
        ) %>% 
      dplyr::select(key, prov, scientificName, EORANK, stateProvince,	countryCode, year, month, basisOfRecord, datasetName,	institutionCode, coordinateUncertaintyInMeters, longitude, latitude, references) %>%
      dplyr::rename("Key" = key,
                    "Scientific name" = scientificName,
                    "Source" = prov,
                    "Latitude" = latitude,
                    "Longitude" = longitude,
                    "Dataset Name" = datasetName,
                    "Institution" = institutionCode,
                    "Year" = year,
                    "Month" = month,
                    "Coordinate Uncertainty" = coordinateUncertaintyInMeters,
                    "Subnation" = stateProvince,
                    "Nation" = countryCode,
                    "Record type" = basisOfRecord,
                    "EO Rank" = EORANK,
                    "References" = references
      ) 
    
    dat %>%
      DT::datatable(options = list(dom = 'tp',
                                   pageLength = 15,
                                   columnDefs = list(list(width = "10%", className = 'dt-left', targets = c(1,2))),
                                   language = list(emptyTable = 'You have not selected any records from the map'),
                                   width = "100%"
      ),
      # filter = list(position = 'top'),
      selection = list(mode = 'multiple', target = 'row', selected = 0:nrow(taxon_data$selected_points)),
      escape = FALSE,
      rownames = FALSE
      ) %>% 
      formatStyle(1:ncol(dat),"white-space"="nowrap")
    
  })
  
  observeEvent(input$occurrences_table_rows_selected, {
    
    taxon_data$selected_points <- taxon_data$selected_points %>% 
      dplyr::filter(taxon_data$selected_points$key %in% taxon_data$selected_points$key[input$occurrences_table_rows_selected])
    
    clicks$IDs <- intersect(clicks$IDs, taxon_data$selected_points$key)
    
    m <- leafletProxy("main_map") %>%
      clearMarkerClusters() %>%
      clearMarkers() %>%
      leaflet::addMapPane("selected_observations", zIndex = 400) %>% 
      addCircleMarkers(
        data = taxon_data$selected_points,
        lng = ~longitude,
        lat = ~latitude,
        layerId = ~key,
        fillColor = "#4169E1",
        fillOpacity = 0.75,
        color = "#4169E1",
        options = pathOptions(pane = "selected_observations"),
        group = "Point observations"
      ) %>%
      addCircleMarkers(
        data = taxon_data$sf_filtered %>% dplyr::filter(key %in% setdiff(key, taxon_data$selected_points$key)),
        lng = ~longitude,
        lat = ~latitude,
        # clusterOptions = markerClusterOptions(),
        layerId = ~key,
        fillColor = "#4169E1",
        fillOpacity = 0.25,            
        color = grey(0.25),
        options = pathOptions(pane = "species_observations")
      )
    
    m
    
  })
  
  observeEvent(input$clear_selected_records, { 
    
    # shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
    
    taxon_data$selected_points <- taxon_data$selected_points %>% 
      dplyr::filter(key %in% taxon_data$selected_points$key[NULL])  
    
    # shinybusy::remove_modal_spinner()
    
  })
  
  observeEvent(input$no_year, { 
      
      if (isTRUE(input$no_year)){
        
        no_year_keys <- taxon_data$sf_filtered %>% dplyr::filter(is.na(year)) %>% dplyr::pull(key)
        
        if (length(no_year_keys) > 0){
          taxon_data$selected_points <- taxon_data$sf_filtered %>% 
            dplyr::filter(key %in% no_year_keys)
          
        } else {
          sendSweetAlert(session, type = "warning", title = "Oops!", text = "No records seem to be missing an assigned year value", closeOnClickOutside = TRUE)
        }
      } else {
        
        taxon_data$selected_points <- data.frame("Key" = character(), "Scientific name" = character(), "Source" = character(), "Institution code" = character(), "Year" = numeric(), "Coordinate Uncertainty" = numeric(), "Place" = character(), "URL" = character())[NULL, ]
        
      }
    
  })
  
  observeEvent(input$range_extent, {

    if (isTRUE(input$range_extent)){

      # shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
      
      shinyjs::show(id = "EOO_panel")

      # taxon_data$filtered_occurrences <- taxon_data$sf_filtered
      # 
      # taxon_data$filtered_occurrences <- taxon_data$filtered_occurrences %>%
      #   st_set_geometry(NULL) %>%
      #   dplyr::select(longitude, latitude) %>%
      #   as.data.frame()
      
      if (nrow(taxon_data$sf_filtered) >= 3){

        if (taxon_data$shifted){
          eoo_output <- taxon_data$sf_filtered %>% safe_eoo(shifted = TRUE)
        } else {
          eoo_output <- taxon_data$sf_filtered %>% safe_eoo(shifted = FALSE)
        }
        
        if (!is.null(eoo_output$result)){

        eoo_output <- eoo_output$result
        taxon_data$species_range_value <- eoo_output$EOO
        taxon_data$species_range_map <- eoo_output$hull
        taxon_data$species_range_factor <- eoo_output$factor %>% as.character()
          
        m <- leafletProxy("main_map") %>%
          clearGroup("Range Extent") %>%
          leaflet::addMapPane("species_range", zIndex = 200) %>%
          addPolygons(data = taxon_data$species_range_map,
                      color = grey(.2),
                      fillOpacity = 0.1,
                      fill = TRUE,
                      options = pathOptions(pane = "species_range"),
                      group = "Range Extent"
          )

        m
        } else {
          
          taxon_data$species_range_value <- NULL
          
          leafletProxy("main_map") %>%
            leaflet::clearGroup("Range Extent")
          
          sendSweetAlert(session, type = "warning", title = "Oops!", text = "Range Extent could not be calculated correctly for this taxon", closeOnClickOutside = TRUE)

        }
        } else {
          
          taxon_data$species_range_value <- NULL
          
          leafletProxy("main_map") %>%
            leaflet::clearGroup("Range Extent")
          
          sendSweetAlert(session, type = "warning", title = "Oops!", text = "Range Extent could not be calculated correctly for this taxon", closeOnClickOutside = TRUE)
          
        }

    } else {
      
      taxon_data$species_range_value <- NULL
      
      shinyjs::hide(id = "EOO_panel")
      leafletProxy("main_map") %>%
        leaflet::clearGroup("Range Extent")
    }

    # shinybusy::remove_modal_spinner()
    
  })
  
  observeEvent({
    input$area_of_occupancy
    input$grid_cell_size
  }, {
    
    if (isTRUE(input$area_of_occupancy)){
      
      # shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
      
      shinyjs::show(id = "AOO_panel")
      
      taxon_data$filtered_occurrences <- taxon_data$sf_filtered
      taxon_data$filtered_occurrences <- taxon_data$filtered_occurrences %>%
        st_set_geometry(NULL) %>%
        dplyr::select(longitude, latitude) %>%
        as.data.frame()
      taxon_data$filtered_occurrences$longitude[taxon_data$filtered_occurrences$longitude > 180] <- taxon_data$filtered_occurrences$longitude[taxon_data$filtered_occurrences$longitude > 180] - 360
      
      if (nrow(taxon_data$filtered_occurrences) >= 1){
        
      taxon_data$AOO_value <- purrr::safely(aoo2)(taxon_data$filtered_occurrences, as.numeric(input$grid_cell_size)*1000)
      if (!is.null(taxon_data$AOO_value$result)){
        taxon_data$AOO_value <- taxon_data$AOO_value$result/4
        taxon_data$AOO_factor <- purrr::safely(get_aoo_factor)(taxon_data$AOO_value, grid_cell_size = as.numeric(input$grid_cell_size))
        taxon_data$AOO_factor <- ifelse(!is.null(taxon_data$AOO_factor$result), as.character(taxon_data$AOO_factor$result), NULL)
        taxon_data$AOO_map <- purrr::safely(get_aoo_polys)(taxon_data$sf_filtered, as.numeric(input$grid_cell_size))      
        if (!is.null(taxon_data$AOO_map$result)){
          taxon_data$AOO_map <- taxon_data$AOO_map$result
        } 
      } 
      }
      
      if (!is.null(taxon_data$AOO_map)){
        if (nrow(taxon_data$filtered_occurrences) >= 2){
          
      m <- leafletProxy("main_map") %>%
        clearGroup("Occupancy") %>% 
        leaflet::addMapPane("aoo", zIndex = 200) %>% 
        addPolygons(data = taxon_data$AOO_map,
                    color = "#2c7bb6",
                    fill = "#2c7bb6",
                    fillOpacity = .2,
                    options = pathOptions(pane = "aoo"),
                    group = "Occupancy"
        )
        } else {
          m <- leafletProxy("main_map") %>%
            clearShapes()
          sendSweetAlert(session, type = "warning", title = "Oops!", text = "AOO could not be mapped for this taxon", closeOnClickOutside = TRUE)
        }
      } else {
        m <- leafletProxy("main_map") %>%
          clearShapes()
        sendSweetAlert(session, type = "warning", title = "Oops!", text = "AOO could not be mapped for this taxon", closeOnClickOutside = TRUE)
      }
      
      m
      
      # shinybusy::remove_modal_spinner()
      
    } else {
      
      shinyjs::hide(id = "AOO_panel")
      leafletProxy("main_map") %>%
        leaflet::clearGroup("Occupancy")
      
    }
  })
  
  observeEvent({
    input$number_EOs
    input$separation_distance
  }, {
    
    if (isTRUE(input$number_EOs)){
      
      # shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
      
      shinyjs::show(id = "EOcount_panel")
      
      if (nrow(taxon_data$sf_filtered) >= 1){
        
      ##### Calculate numbers of EOs
      number_EOs <- purrr::safely(calculate_number_occurrences)(taxon_data$sf_filtered, separation_distance = input$separation_distance %>% as.numeric(), added_distance = 0)
      if (!is.null(number_EOs$result)){
        number_EOs <- number_EOs$result
        taxon_data$EOcount_value <- number_EOs$eo_count
        taxon_data$EOcount_map <- number_EOs$buffered_occurrences
        taxon_data$EOcount_factor <- number_EOs$factor %>% as.character()
      }
      }

      if (!is.null(taxon_data$EOcount_map)){
        
      m <- leafletProxy("main_map") %>%
        clearGroup("Occurrences") %>%
        leaflet::addMapPane("eos", zIndex = 200) %>% 
        addPolygons(data = taxon_data$EOcount_map, 
                    fill = TRUE,
                    fillColor = "#8b0000",
                    fillOpacity = .5,
                    opacity = 0,
                    options = pathOptions(pane = "eos"),
                    group = "Occurrences"
        )
      
      } else {
        m <- leafletProxy("main_map") %>%
          clearShapes()
        sendSweetAlert(session, type = "warning", title = "Oops!", text = "Occurrences could not be mapped for this taxon", closeOnClickOutside = TRUE)
      }
      
      m
      
      # shinybusy::remove_modal_spinner()
      
    } else {
      
      shinyjs::hide(id = "EOcount_panel")
      leafletProxy("main_map") %>%
        leaflet::clearGroup("Occurrences")
    }
  })
  
  output$species_name <- renderUI({
    
    if (!is.null(taxon_data$info)){
      
      if (taxon_data$info$scientificName == "New taxon"){
        span(style = "; display-inline: block;",
             h2(strong(em("New taxon", style = "padding-top: 0; padding-bottom: 7px; margin-top: 0; display-inline: block; float:left; padding-right:10px;")))
        )
      } else {
        if (taxon_data$info$Source == "NatureServe"){
          display_ID <- paste0("(NatureServe ID: ", taxon_data$info$elementGlobalId, ")")
          display_href <- paste0("https://explorer.natureserve.org/Taxon/", taxon_data$info$uniqueId)
        } else if (taxon_data$info$Source == "GBIF"){
          display_ID <- paste0("(GBIF ID: ", taxon_data$info$elementGlobalId, ")")
          display_href <- paste0("https://www.gbif.org/species/", taxon_data$info$elementGlobalId)
        } else if (taxon_data$info$Source == "Uploaded"){
          display_ID <- "(Detected from uploaded file)"
          display_href <- paste0("http://resolver.globalnames.org/name_resolvers?names=", gsub(" ", "+", taxon_data$info$scientificName))
        }
        span(style = "; display-inline: block;",
             h2(strong(em(req(taxon_data$info$scientificName)), style = "padding-top: 0; margin-top: 0; display-inline: block; float:left; padding-right:10px;")),
             h2(strong(a(display_ID, href = display_href, target = "_blank")))
        )
      }
      
    } 
  })
  
  output$species_range_value <- renderUI({
    
    if (!is.null(taxon_data$species_range_value)){
      
      HTML(
        paste0(
          h2(format(as.numeric(req(taxon_data$species_range_value)), big.mark=",", scientific = FALSE), style = "padding-top: 0; margin-top: 0; display: inline;"), 
          h2(" km", tags$sup("2"), style = "padding-top: 0; margin-top: 0; display: inline;"),
          h1(strong(taxon_data$species_range_factor), style = "padding-top: 0; margin-top: 0; padding-left: 1em; display: inline;")
        )
      )
    } else {
      HTML(paste0(
        h4("Could not be calculated", style = "padding-top: 0; margin-top: 0; display: inline;")
      ))
      }
    
  })
  
  output$AOO_value <- renderUI({
    
    if (!is.null(taxon_data$AOO_value)){
      
    if (input$grid_cell_size == 1) out <- base::cut(req(as.numeric(taxon_data$AOO_value)), breaks = c(0, 0.999, 4.999, 10.999, 20.999, 100.999, 500.999, 2000.999, 10000.999, 50000.999, 1000000000), labels = c("Z", LETTERS[1:9]))
    
    if (input$grid_cell_size > 1) out <- base::cut(req(as.numeric(taxon_data$AOO_value)), breaks = c(0, 0.999, 1.999, 2.999, 5.999, 25.999, 125.999, 500.999, 2500.999, 12500.999, 1000000000), labels = c("Z", LETTERS[1:9]))
    
    HTML(
      paste0(
        h2(format(as.numeric(req(taxon_data$AOO_value)), big.mark=",", scientific = FALSE), style = "padding-top: 0; margin-top: 0; display: inline;"), 
        h2(paste0(" cells"), style = "padding-top: 0; margin-top: 0; display: inline;"),
        h1(strong(taxon_data$AOO_factor), style = "padding-top: 0; margin-top: 0; padding-left: 1em; display: inline;")
      )
    ) 
    } else {
      HTML(paste0(
        h4("Could not be calculated", style = "padding-top: 0; margin-top: 0; display: inline;")
        ))
    }
    
  })
  
  output$EOcount_value <- renderUI({

    if (!is.null(taxon_data$EOcount_value)){    
      HTML(
        paste0(
          h2(format(as.numeric(req(taxon_data$EOcount_value)), big.mark=",", scientific = FALSE), style = "padding-top: 0; margin-top: 0; display: inline;"), 
          h2("  EOs", style = "padding-top: 0; margin-top: 0; display: inline;"),
          h1(strong(taxon_data$EOcount_factor), style = "padding-top: 0; margin-top: 0; padding-left: 1em; display: inline;")
        )
      )
    } else {
      HTML(paste0(
        h4("Could not be calculated", style = "padding-top: 0; margin-top: 0; display: inline;")
      ))
    }
    
  })
  
  output$number_occurrences <- renderUI({
    
    if (!is.null(taxon_data$sf_filtered)){
      sources_count <- taxon_data$sf_filtered %>% dplyr::count(prov)
      sources_filter_labels <- paste0(paste0(sources_count$prov, ": ", sources_count$n), collapse = "; ")
      records_count_text <- ifelse(sum(sources_count$n) > 0, paste0("  records included after filtering (", sources_filter_labels, ")"), "  records included after filtering")
      
      HTML(
        paste0(
          strong(h3(format(as.numeric(req(nrow(taxon_data$sf_filtered))), big.mark=",", scientific = FALSE), style = "padding-top: 0; margin-top: 0; display: inline;")), 
          strong(h4(records_count_text, style = "padding-top: 0; margin-top: 0; display: inline;"))
        )
      )
    }

  })
  
  output$download_rank_data <- downloadHandler(
    filename = function() {
      paste(gsub(" ", "_", ifelse(!is.null(taxon_data$info), taxon_data$info$scientificName, "New taxon")), "-rank_factor_values-", Sys.Date(), ".xlsx", sep="")
    },
    content = function(file) {
      out <- matrix(nrow = 1, ncol = 42, data = "") %>%
        as.data.frame()
      names(out)[c(1:3, 6:9, 11, 13, 15, 16:17, 19:20, 22:28, 30, 32:42)] <- c("Calc Rank", "Assigned Rank", "Species or Community Scientific Name*", "Element ID",
                                                                               "Elcode*", "Common Name*", "Classification*", "Range Extent", "Area of Occup 4-km2 grid cells",
                                                                               "# Occur", "Pop Size", "# Occur Good Viab", "Environm Specif (opt.)", "Overall Threat Impact",
                                                                               "Intrinsic Vulner (opt.)", "Short-term Trend", "Long-term Trend", "Rank Adjustment Reasons", "Assigned Rank Reasons", "Rank Factors Author", 
                                                                               "Rank Factors Date", "Rank Review Date", "Range Extent Comments", "Area of Occupancy Comments",
                                                                               "# of Occurrences Comments", "Population Size Comments", "Good Viability/Integrity Comments", "Environmental Specificity Comments",
                                                                               "Threat Impact Comments", "Threat Impact Adjustment Reasons", "Intrinsic Vulnerability Comments", "Short-term Trend Comments", "Long-term Trend Comments"
      )
      
      out[, 3] <-  ifelse(!is.null(taxon_data$info), taxon_data$info$scientificName, "")
      if (!is.null(taxon_data$species_range_value)){
        out[, 11] <- taxon_data$species_range_factor
      }
      if (!is.null(taxon_data$AOO_value)){
        out[, 13] <- taxon_data$AOO_factor
      }
      if (!is.null(taxon_data$EOcount_value)){
        out[, 15] <- taxon_data$EOcount_factor
      }
      
      if (input$single_assessment_type == "global"){
        
      out[, 2] <-  ifelse(!is.null(taxon_data$info), taxon_data$info$roundedGRank, "")
      out[, 6] <-  ifelse(!is.null(taxon_data$info), taxon_data$info$elementGlobalId, "")
      out[, 7] <-  ifelse(!is.null(taxon_data$info), taxon_data$info$elcode, "")
      out[, 8] <-  ifelse(!is.null(taxon_data$info), taxon_data$info$primaryCommonName, "")
      out[, 9] <-  ifelse(!is.null(taxon_data$info_extended), taxon_data$info_extended$classificationStatus$classificationStatusDescEn, "")

      if (!is.null(taxon_data$info_extended$rankInfo$popSize)){
        rank_def <- rank_factor_definitions$population_size_description[grep(taxon_data$info_extended$rankInfo$popSize$popSizeDescEn, rank_factor_definitions$population_size_description, fixed = TRUE)]
        out[, 16] <- rank_factor_definitions$population_size_value[rank_factor_definitions$population_size_description == rank_def]
      }
      if (!is.null(taxon_data$info_extended$rankInfo$numberGoodEos$numberGoodEosDescEn)){
        rank_def <- rank_factor_definitions$good_viability_description[grep(taxon_data$info_extended$rankInfo$numberGoodEos$numberGoodEosDescEn, rank_factor_definitions$good_viability_description, fixed = TRUE)]
        out[, 17] <- rank_factor_definitions$good_viability_value[rank_factor_definitions$good_viability_description == rank_def]
      }
      if (!is.null(taxon_data$info_extended$rankInfo$enviromentalSpecificity)){
        rank_def <- rank_factor_definitions$environmental_specificity_description[grep(gsub(" |\\.", "", taxon_data$info_extended$rankInfo$enviromentalSpecificity$enviromentalSpecificityDescEn), gsub(" |\\.", "", rank_factor_definitions$environmental_specificity_description), fixed = TRUE)]
        out[, 19] <- rank_factor_definitions$environmental_specificity_value[rank_factor_definitions$environmental_specificity_description == rank_def]
      }
      if (!is.null(taxon_data$info_extended$rankInfo$threatImpactAssigned)){
        rank_def <- rank_factor_definitions$threat_impact_description[grep(taxon_data$info_extended$rankInfo$threatImpactAssigned$threatImpactAssignedDescEn, rank_factor_definitions$threat_impact_description, fixed = TRUE)]
        out[, 20] <- rank_factor_definitions$threat_impact_value[rank_factor_definitions$threat_impact_description == rank_def]
      }
      if (!is.null(taxon_data$info_extended$rankInfo$intrinsicVulnerability)){
        rank_def <- rank_factor_definitions$intrinsic_vulnerability_description[grep(taxon_data$info_extended$rankInfo$intrinsicVulnerability$intrinsicVulnerabilityDescEn, rank_factor_definitions$intrinsic_vulnerability_description, fixed = TRUE)]
        out[, 22] <- rank_factor_definitions$intrinsic_vulnerability_value[rank_factor_definitions$intrinsic_vulnerability_description == rank_def]
        
      }
      if (!is.null(taxon_data$info_extended$rankInfo$shortTermTrend)){
        rank_def<- rank_factor_definitions$temporal_trend_description[grep(gsub(" ", "", taxon_data$info_extended$rankInfo$shortTermTrend$shortTermTrendDescEn), gsub(" ", "", rank_factor_definitions$temporal_trend_description), fixed = TRUE)]
        out[, 23] <- rank_factor_definitions$temporal_trend_value[rank_factor_definitions$temporal_trend_description == rank_def]
        
      }
      if (!is.null(taxon_data$info_extended$rankInfo$longTermTrend)){
        rank_def <- rank_factor_definitions$temporal_trend_description[grep(gsub(" ", "", taxon_data$info_extended$rankInfo$longTermTrend$longTermTrendDescEn), gsub(" ", "", rank_factor_definitions$temporal_trend_description), fixed = TRUE)]
        out[, 24] <- rank_factor_definitions$temporal_trend_value[rank_factor_definitions$temporal_trend_description == rank_def]
        
      }
      out[, 25] <- ifelse(!is.null(taxon_data$info_extended$grankAdjustmentReasons), taxon_data$info_extended$grankAdjustmentReasons, "")
      out[, 26] <- ifelse(!is.null(taxon_data$info_extended$grankReasons), taxon_data$info_extended$grankReasons, "")
      out[, 27] <- ifelse(!is.null(taxon_data$info_extended$conservationStatusFactorsEditionAuthors), taxon_data$info_extended$conservationStatusFactorsEditionAuthors, "")
      out[, 28] <- ifelse(!is.null(taxon_data$info_extended$conservationStatusFactorsEditionDate), taxon_data$info_extended$conservationStatusFactorsEditionDate, "")
      out[, 30] <- ifelse(!is.null(taxon_data$info_extended$grankReviewDate), taxon_data$info_extended$grankReviewDate, "")
      out[, 32] <- ifelse(!is.null(taxon_data$info_extended$rankInfo$rangeExtentComments), taxon_data$info_extended$rankInfo$rangeExtentComments, "")
      out[, 33] <- ifelse(!is.null(taxon_data$info_extended$rankInfo$areaOfOccupancyComments), taxon_data$info_extended$rankInfo$areaOfOccupancyComments, "")
      out[, 34] <- ifelse(!is.null(taxon_data$info_extended$rankInfo$numberEosComments), taxon_data$info_extended$rankInfo$numberEosComments, "")
      out[, 35] <- ifelse(!is.null(taxon_data$info_extended$rankInfo$popSizeComments), taxon_data$info_extended$rankInfo$popSizeComments, "")
      out[, 36] <- ifelse(!is.null(taxon_data$info_extended$rankInfo$viabilityComments), taxon_data$info_extended$rankInfo$viabilityComments, "")
      out[, 37] <- ifelse(!is.null(taxon_data$info_extended$rankInfo$enviromentalSpecificityComments), taxon_data$info_extended$rankInfo$enviromentalSpecificityComments, "")
      out[, 38] <- ifelse(!is.null(taxon_data$info_extended$rankInfo$threatImpactComments), taxon_data$info_extended$rankInfo$threatImpactComments, "")
      out[, 39] <- ifelse(!is.null(taxon_data$info_extended$rankInfo$grankAdjustmentReasons), taxon_data$info_extended$rankInfo$grankAdjustmentReasons, "")
      out[, 40] <- ifelse(!is.null(taxon_data$info_extended$rankInfo$intrinsicVulnerabilityComments), taxon_data$info_extended$rankInfo$intrinsicVulnerabilityComments, "")
      out[, 41] <- ifelse(!is.null(taxon_data$info_extended$rankInfo$shortTermTrendComments), taxon_data$info_extended$rankInfo$shortTermTrendComments, "")
      out[, 42] <- ifelse(!is.null(taxon_data$info_extended$rankInfo$longTermTrendComments), taxon_data$info_extended$rankInfo$longTermTrendComments, "")
      
      }
      
      out2_names <- c("NatureServe accepted name", "NatureServe synonyms", "GBIF taxonomic concepts with GBIF IDs", "EGT ID", "EGT UID", "ELCODE", 
                      "Assessment Type", "Nations included", "Subnations included", "New Range Extent value (sq km)", "New Range Extent letter",
                      "Previous Range Extent letter", "Compare Range Extent letter (new vs. previous)", "New Area of Occupancy grid cell size",
                      "New Area of Occupancy value", "New Area of Occupancy letter", "Previous Area of Occupancy letter", "Compare Area of Occupancy letter (new vs. previous)",
                      "RARECAT occurrence separation distance (m)", "New Number of Occurrences value", "New Number of Ocurrences letter", "Previous Number of Occurrences letter",
                      "Compare Number of Occurrences letter (new vs. previous)", "Time Period 1", "Time Period 2", "Range Extent temporal change analysis", 
                      "Area of Occupancy temporal change analysis", "Number of Occurrences temporal change analysis", "Number of records included",
                      "Data sources included", "Date range of records included", "Months of records included", "Locational uncertainty cutoff", "Other filters", "New assessment date"
      )
      
      out2 <- matrix(nrow = 1, ncol = length(out2_names), data = "") %>%
        as.data.frame()
      names(out2) <- out2_names
      
      if (input$single_assessment_type == "global"){
        out2[, 4] <- taxon_data$info_extended$elementGlobalId
        out2[, 5] <- taxon_data$info_extended$uniqueId
        out2[, 6] <- taxon_data$info_extended$elcode
      }
      out2[, 1] <- ifelse(!is.null(taxon_data$info), taxon_data$info$scientificName, "")
      out2[, 2] <- ifelse(length(taxon_data$info$synonyms %>% unlist() %>% na.omit() %>% as.character()) > 0, paste0(taxon_data$info$synonyms %>% unlist() %>% na.omit() %>% as.character(), collapse = "; "), "")
      out2[, 3] <- ifelse(!is.null(taxon_data$synonyms$scientificName), 
                      paste0(taxon_data$synonyms$scientificName %>% na.omit() %>% as.character(), collapse = "; "), 
                      ""
                      )
      out2[, 7] <- input$single_assessment_type
      out2[, 8] <- ifelse(!is.null(input$nation_filter), paste0(input$nation_filter, collapse = "; "), "")
      out2[, 9] <- ifelse(!is.null(input$states_filter), paste0(input$states_filter, collapse = "; "), "")
      out2[, 10] <- ifelse(!is.null(taxon_data$species_range_value), taxon_data$species_range_value %>% as.numeric(), "")
      out2[, 11] <- ifelse(!is.null(taxon_data$species_range_factor), taxon_data$species_range_factor %>% as.character(), "")
      out2[, 14] <- input$grid_cell_size %>% as.numeric()
      out2[, 15] <- ifelse(!is.null(taxon_data$AOO_value), taxon_data$AOO_value %>% as.numeric(), "")
      out2[, 16] <- ifelse(!is.null(taxon_data$AOO_factor), taxon_data$AOO_factor %>% as.character(), "")
      out2[, 19] <- input$separation_distance %>% as.numeric()
      out2[, 20] <- ifelse(!is.null(taxon_data$EOcount_value), taxon_data$EOcount_value %>% as.numeric(), "")
      out2[, 21] <- ifelse(!is.null(taxon_data$EOcount_factor), taxon_data$EOcount_factor %>% as.character(), "")
      out2[, 24] <- paste0(input$period1, collapse = " - ")
      out2[, 25] <- paste0(input$period2, collapse = " - ")
      out2[, 26] <- taxon_data$temporal_change$eoo_change[2]
      out2[, 27] <- taxon_data$temporal_change$aoo_change[2]
      out2[, 28] <- taxon_data$temporal_change$eo_count_change[2]
      out2[, 29] <- nrow(taxon_data$sf_filtered)
      out2[, 30] <- paste0(input$sources_filter, collapse = "; ")
      out2[, 31] <- paste0(input$year_filter, collapse = "; ")
      out2[, 32] <- paste0(input$seasonality, collapse = "; ")  
      out2[, 33] <- input$uncertainty_filter
      out2[, 34] <- paste0(
        c(ifelse(input$clean_occ, "Basic GBIF data cleanup implemented", ""),
          ifelse(input$clean_occ, "Putative centroids removed", ""),
          paste0("Data types included: ", input$type_filter),
          ifelse(length(taxon_data$sf_filtered$EORANK %>% complete.cases(.)) > 0, paste0("Element occurrence ranks included: ", paste0(unique(taxon_data$sf_filtered$EORANK), collapse = "|")), "")
        ), collapse = "; ")
      out2[, 35] <- Sys.Date()

      taxon_data$rank_factor_comparison <- compare_rank_factors(taxon_data)
      out2[, 12] <- taxon_data$rank_factor_comparison$previous_species_range_letter
      out2[, 13] <- taxon_data$rank_factor_comparison$new_previous_species_range_comparison
      out2[, 17] <- taxon_data$rank_factor_comparison$previous_aoo_letter
      out2[, 18] <- taxon_data$rank_factor_comparison$new_previous_aoo_comparison
      out2[, 22] <- taxon_data$rank_factor_comparison$previous_eocount_letter
      out2[, 23] <- taxon_data$rank_factor_comparison$new_previous_eocount_comparison
      
      writexl::write_xlsx(
        list(
          "Rank Calculator" = out,
          "Analysis Summary" = out2
        ), 
        path = file
      )
      
    }
  )
  
  output$download_occurrence_data <- downloadHandler(
    filename = function() {
      paste(gsub(" ", "_", ifelse(!is.null(taxon_data$info), taxon_data$info$scientificName, "New taxon")), "-RARECAT_assessment_records-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      out <- taxon_data$sf_filtered %>% 
        st_set_geometry(NULL) %>% 
        dplyr::mutate(egt_uid = "",
                      el_code = "",
                      scientificName_Assessment = "New taxon",
                      scientificName_Source = "New taxon"
        ) 
      
      if (!is.null(taxon_data$info$Source)){
        out <- out %>% 
          dplyr::mutate(egt_uid = ifelse(taxon_data$info$Source == "NatureServe", taxon_data$info$uniqueId, ""),
                        el_code = ifelse(taxon_data$info$Source == "NatureServe", taxon_data$info$elcode, ""),
                        scientificName_Assessment = taxon_data$info$scientificName,
                        scientificName_Source = scientificName
          ) 
      }
      
      extra_gbif_fields <- c("key", "collectionCode", "recordedBy", "recordNumber", "samplingProtocol", "accessRights", "taxonKey", "samplingProtocol", "occurrenceID", "bibliographicCitation", "datasetID", "license", "rightsHolder", "eventDate", "occurrenceStatus", "georeferenceProtocol")
      if (!is.null(taxon_data$gbif_occurrences_raw)){
        out <- out %>% 
          dplyr::left_join(taxon_data$gbif_occurrences_raw %>% dplyr::select(any_of(extra_gbif_fields)), by = "key")
      }

      out <- out %>% 
        cbind(matrix("", nrow = nrow(out), ncol = length(setdiff(extra_gbif_fields, names(out)))) %>%
                as.data.frame() %>%
                set_names(setdiff(extra_gbif_fields, names(out)))
        ) 
      
      out <- out %>%
        dplyr::select(all_of(c("key", "scientificName_Assessment", "scientificName_Source", out %>% dplyr::select(-key, -scientificName_Assessment, -scientificName_Source, -scientificName) %>% names())))
      
      write.csv(out, file, row.names = FALSE, na = "")
    }
  )
  
  observeEvent(input$single_assessment_clear, {
    
    leafletProxy("main_map") %>%
      clearShapes() %>% 
      clearMarkers() %>% 
      clearMarkerClusters() 
    taxon_data$info <- data.frame(scientificName = "New taxon")
    taxon_data$gbif_occurrences_raw <- NULL
    taxon_data$gbif_occurrences <- NULL
    synonyms <- NULL
    synonyms_selected <- NULL
    taxon_data$uploaded_occurrences <- NULL
    taxon_data$drawn_occurrences <- NULL
    taxon_data$all_occurrences <- NULL
    taxon_data$circumboreal <- FALSE
    taxon_data$sf <- NULL
    taxon_data$sf_filtered <- NULL
    taxon_data$filtered_occurrences <- NULL
    taxon_data$selected_points <- data.frame("Key" = character(), "Scientific name" = character(), "Source" = character(), "Institution code" = character(), "Year" = numeric(), "Coordinate Uncertainty" = numeric(), "Place" = character(), "URL" = character())[NULL, ]
    taxon_data$removed_points <- NULL
    taxon_data$nations <- NULL
    taxon_data$states <- NULL
    taxon_data$records_over_time <- NULL
    taxon_data$species_range_value <- NULL
    taxon_data$species_range_map <- NULL
    taxon_data$AOO_value <- NULL
    taxon_data$AOO_map <- NULL
    taxon_data$EOcount_map <- NULL
    taxon_data$EOcount_value <- NULL
    
    updateSelectizeInput(session = session, inputId = "single_assessment_type", selected = "global")
    updateMaterialSwitch(session = session, inputId = "load_gbif_data", value = FALSE)
    updateMaterialSwitch(session = session, inputId = "map_uploads", value = FALSE)    
    updateMaterialSwitch(session = session, inputId = "range_extent", value = FALSE)
    updateMaterialSwitch(session = session, inputId = "area_of_occupancy", value = FALSE)
    updateMaterialSwitch(session = session, inputId = "number_EOs", value = FALSE)
    
    shinyjs::hide(id = "data_panel")
    shinyjs::hide(id = "taxon_search_panel")
    shinyjs::hide(id = "taxon_options_panel")
    shinyjs::hide(id = "data_panel")
    shinyjs::hide(id = "load_data_panel")
    shinyjs::hide(id = "EOO_panel")
    shinyjs::hide(id = "AOO_panel")
    shinyjs::hide(id = "EOcount_panel")
    updateCollapse(session = session, id = "inputs_single", open = "Add assessment data")
    
    # if (!is.null(input$filedata)){
    unlink(input$filedata$datapath)
    reset(id = "filedata")
    # }
  })
  
  batch_assessment_polygon <- reactiveValues(gadmGid = NULL)
  batch_run_taxon_list <- reactiveValues(names = NULL)
  batch_run_output <- reactiveValues(
    results = NULL,
    table = NULL
  )
  
  
  #' ### Load user-uploaded data
  uploaded_obs_data_batch <- reactive({

    out <- purrr::map(input$batch_filedata_obs$datapath, process_user_data, minimum_fields = c(minimum_fields, "scientificName_Source", "scientificName_Assessment")) %>% 
      dplyr::bind_rows() 
    
    out 
  })
  
  observeEvent(input$batch_assessment_type, {
    
    if (input$batch_assessment_type == "global"){
      shinyjs::hide(id = "batch_nation")
      shinyjs::hide(id = "batch_subnation")
    }
    
    if (input$batch_assessment_type == "national"){
      shinyjs::show(id = "batch_nation")
      
    }
    
    if (input$batch_assessment_type == "subnational"){
      shinyjs::show(id = "batch_nation")
      shinyjs::show(id = "batch_subnation")
    }
    
  })
  
  observeEvent(input$batch_nation_filter, {

      nation_subset <- network_polys %>% dplyr::filter(FIPS_CNTRY %in% input$batch_nation_filter)
      
      subnations_already_selected <- input$batch_states_filter
      
      updateSelectizeInput(session = session,
                           inputId = "batch_states_filter",
                           choices = nation_subset$ADMIN_NAME %>% na.omit() %>% as.character() %>% sort(),
                           selected = subnations_already_selected
      )
      
      if (input$batch_assessment_type == "national"){
        batch_assessment_polygon$gadmGid <- ifelse(input$batch_nation_filter == "US", "USA", ifelse(input$batch_nation_filter == "CA", "CAN", NULL))
      }
      
  })
  
  observeEvent(input$batch_states_filter, {

    relevant_nation <- network_polys %>% dplyr::filter(ADMIN_NAME %in% input$batch_states_filter) %>% dplyr::pull(FIPS_CNTRY)
    
    updateSelectizeInput(session = session,
                         inputId = "batch_nation_filter", 
                         selected = relevant_nation %>% set_names(relevant_nation)
    )
    
    if (input$batch_assessment_type == "subnational"){
    batch_assessment_polygon$gadmGid <- gadm_df %>% 
      dplyr::filter(NAME_1 %in% gsub(" ", "", input$batch_states_filter)) %>% 
      dplyr::pull(GID_1) %>% 
      unique()
    }
  })
  
  observeEvent(input$batch_assessment, {

    # shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
    
    if (!is.null(input$batch_filedata_obs$datapath)){
      batch_uploaded_occurrences <- uploaded_obs_data_batch()
      # batch_uploaded_occurrences <- batch_uploaded_occurrences %>% 
      #   dplyr::mutate(key = paste(prov, 1:nrow(batch_uploaded_occurrences), sep = "_"))
      # uploaded_names <- purrr::map(batch_uploaded_occurrences$scientificName %>% unique(), function(sp){
      #   out <- rgbif::name_suggest(q = sp, rank = c("species", "subspecies", "variety", "infraspecific_name"), limit = 10)$data
      #   out$uploaded_name <- sp
      #   out
      if ("scientificName_Assessment" %in% names(batch_uploaded_occurrences)){
        if (length(complete.cases(batch_uploaded_occurrences$scientificName_Assessment)) > 0){
          uploaded_names <- c(batch_uploaded_occurrences$scientificName_Assessment) %>% unique()
        }
      } else {
        uploaded_names <- c(batch_uploaded_occurrences$scientificName_Assessment) %>% unique()
      }

    } else {
      batch_uploaded_occurrences <- NULL
      uploaded_names <- NULL
    }
    
    if (!is.null(input$batch_filedata_rank$datapath)){
      batch_rank_factor_file <- read.csv(input$batch_filedata_rank$datapath, header = TRUE) 
      rank_names <- batch_rank_factor_file[, 3] %>% as.character()
    } else {
      rank_names <- NULL
    }
    
    typed_names <- strsplit(input$typed_list, "\n|,|;|, |; ")[[1]] 
    
    batch_run_taxon_list$names <- c(rank_names, uploaded_names, typed_names) %>% 
      unique() %>% 
      as.data.frame() %>% 
      set_names("user_supplied_name")
    
    if (nrow(batch_run_taxon_list$names) > 0){
      
    batch_run_output$results <- vector("list", length(batch_run_taxon_list$names$user_supplied_name))
    
    updateCollapse(session = session, id = "batch_parameters", close = "Rank assessment parameters")

    time1 <- Sys.time()
      batch_run_output$results <- safe_batch_run(
        taxon_names = batch_run_taxon_list$names$user_supplied_name,
        minimum_fields = c("key", "scientificName", "prov", "longitude", "latitude", "coordinateUncertaintyInMeters", "stateProvince", "countryCode", "year", "month", "datasetName", "institutionCode", "basisOfRecord", "EORANK", "references"),
        max_number_observations = 10000,
        uploaded_data = batch_uploaded_occurrences,
        rank_factor_upload = input$batch_filedata_rank$datapath,
        clean_occ = input$batch_clean_occ,
        centroid_filter = input$batch_centroid_filter,
        date_start = input$batch_year_filter[1],
        date_end = input$batch_year_filter[2],
        months = input$batch_seasonality,
        uncertainty_filter = input$batch_uncertainty_filter,
        query_polygon = batch_assessment_polygon$gadmGid,
        nations_filter = input$batch_nation_filter,
        states_filter = input$batch_states_filter,
        sources_filter = input$batch_sources_filter,
        grid_cell_size = input$batch_grid_cell_size,
        sep_distance = input$batch_separation_distance,
        trends_period1 = input$batch_period1,
        trends_period2 = input$batch_period2
        )
    time2 <- Sys.time()
    
    print(paste0("batch run took ", difftime(time2, time1, units = "mins")))
    
    batch_run_output$results <- batch_run_output$results$result
    
    # batch_run_output$results <- batch_run_output$results %>%
      # set_names(batch_run_taxon_list$names$user_supplied_name)
    
    batch_run_output$table <- data.frame(
      taxon = names(batch_run_output$results), # batch_run_taxon_list$names$user_supplied_name,
      total_observations_used = purrr::map(batch_run_output$results, function(out) ifelse(!is.null(out$sf_filtered), nrow(out$sf_filtered), nrow(out$all_occurrences))) %>% unlist(),
      range_value = purrr::map(batch_run_output$results, function(out){
        ifelse(!is.null(out$species_range_value)&!is.na(out$species_range_value), paste0(out$species_range_value, " (", out$species_range_factor, ")"), NA)
        }) %>% unlist(),
      current_range_value = purrr::map(batch_run_output$results, function(out){
        if (!is.na(out$rank_factor_comparison$new_previous_species_range_comparison)){
          x <- out$rank_factor_comparison$new_previous_species_range_comparison
          if (x == "equal") res <- as.character(icon("equals", "fa-2x", style = "color: #BEBEBE;"))
          if (x == "lower") res <- as.character(icon("arrow-down", "fa-2x", style = "color: #ef8a62;"))
          if (x == "higher") res <- as.character(icon("arrow-up", "fa-2x", style = "color: #67a9cf;"))
        } else {
          res <- NA
        }
        res
      }) %>% unlist(),
      range_value_trend = purrr::map(batch_run_output$results, function(out){
        out <- ifelse(!is.na(out$temporal_change$eoo_change[2]), paste0(out$temporal_change$eoo_change[2], "%"), NA)
        out
      }) %>% unlist(),
      AOO_value = purrr::map(batch_run_output$results, function(out){
        ifelse(!is.null(out$AOO_value)&!is.na(out$AOO_value), paste0(out$AOO_value, " (", out$AOO_factor, ")"), NA)
        }) %>% unlist(),
      current_AOO_value = purrr::map(batch_run_output$results, function(out){
        if (!is.na(out$rank_factor_comparison$new_previous_aoo_comparison)){
          x <- out$rank_factor_comparison$new_previous_aoo_comparison
          if (x == "equal") res <- as.character(icon("equals", "fa-2x", style = "color: #BEBEBE;"))
          if (x == "lower") res <- as.character(icon("arrow-down", "fa-2x", style = "color: #ef8a62;"))
          if (x == "higher") res <- as.character(icon("arrow-up", "fa-2x", style = "color: #67a9cf;"))
        } else {
          res <- NA
        }
        res
      }) %>% unlist(),
      AOO_value_trend = purrr::map(batch_run_output$results, function(out){
        out <- ifelse(!is.na(out$temporal_change$aoo_change[2]), paste0(out$temporal_change$aoo_change[2], "%"), NA)
        out
      }) %>% unlist(),
      EOcount_value = purrr::map(batch_run_output$results, function(out){
        ifelse(!is.null(out$EOcount_value)&!is.na(out$EOcount_value), paste0(ifelse(out$EOcount_value > 300, ">300", out$EOcount_value), " (", out$EOcount_factor, ")"), NA)
        # ifelse(!is.null(out$EOcount_value), paste0(out$EOcount_value, " (", out$EOcount_factor, ")"), NA)
        }) %>% unlist(),
      current_EOcount_value = purrr::map(batch_run_output$results, function(out){
        if (!is.na(out$rank_factor_comparison$new_previous_eocount_comparison)){
          x <- out$rank_factor_comparison$new_previous_eocount_comparison
        if (x == "equal") res <- as.character(icon("equals", "fa-2x", style = "color: #BEBEBE;"))
        if (x == "lower") res <- as.character(icon("arrow-down", "fa-2x", style = "color: #ef8a62;"))
        if (x == "higher") res <- as.character(icon("arrow-up", "fa-2x", style = "color: #67a9cf;"))
        } else {
          res <- NA
        }
        res
      }) %>% unlist(),
      EOcount_value_trend = purrr::map(batch_run_output$results, function(out){
        out <- ifelse(!is.na(out$temporal_change$eo_count_change[2]), paste0(out$temporal_change$eo_count_change[2], "%"), NA)
        out
      }) %>% unlist(),
      Reviewed = FALSE
    )
    
    if (length(setdiff(batch_run_taxon_list$names$user_supplied_name, batch_run_output$table$taxon)) > 0){
      
    batch_run_output$table <- batch_run_output$table %>% 
      rbind(
        data.frame(
          taxon = setdiff(batch_run_taxon_list$names$user_supplied_name, batch_run_output$table$taxon),
          total_observations_used = "Not available",
          range_value = NA,
          current_range_value = NA,
          range_value_trend = NA,
          AOO_value = NA,
          current_AOO_value = NA,
          AOO_value_trend = NA,
          EOcount_value = NA,
          current_EOcount_value = NA,
          EOcount_value_trend = NA,
          Reviewed = FALSE
        )
      )
    }
    
    shinyjs::show("batch_output")
    
    # shinybusy::remove_modal_spinner()
    
    } else {
      
      sendSweetAlert(session, type = "warning", title = "Oops!", text = "You have not selected any taxa to assess!", closeOnClickOutside = TRUE)
      
    }
    
  })
  
  batch_taxon_focus <- reactiveValues(taxon = NULL)
  
  output$batch_run_results_table <- DT::renderDataTable({
    
    batch_run_table <- batch_run_output$table %>%
      dplyr::rename("Scientific Name (Assessment)" = taxon,
                    "Number Of Records Included" = total_observations_used,
                    "Range Extent (New)" = range_value,
                    "Range Extent (Previous)" = current_range_value,
                    "AOO (New)" = AOO_value,
                    "AOO (Previous)" = current_AOO_value,
                    "Number of Occurrences (New)" = EOcount_value,
                    "Number of Occurrences (Previous)" = current_EOcount_value,
                    "Range Extent Change" = range_value_trend,
                    "AOO Change" = AOO_value_trend,
                    "Number of Occurrences Change" = EOcount_value_trend,
      ) %>% 
      dplyr::mutate(" " = purrr::map(1:nrow(batch_run_output$table), function(i){
        shinyInput(actionButton, 1, paste0(`Scientific Name (Assessment)`[i], "_"), label = "Review assessment", onclick = 'Shiny.onInputChange(\"select_button\",  this.id + "_" + Date.now());')
      }),
      Reviewed = ifelse(Reviewed == FALSE, as.character(icon("xmark", "fa-2x", style = "color: #ef8a62;")), as.character(icon("check", "fa-2x", style = "color: #67a9cf;")))
      ) %>% 
      dplyr::select("Scientific Name (Assessment)", "Number Of Records Included", "Range Extent (New)", "Range Extent (Previous)", "AOO (New)", "AOO (Previous)", "Number of Occurrences (New)", "Number of Occurrences (Previous)", "Range Extent Change", "AOO Change", "Number of Occurrences Change", "Reviewed", " ")
    
    if (input$batch_assessment_type != "global"){
      batch_run_table$`Range Extent (Previous)` <- NA
      batch_run_table$`AOO (Previous)` <- NA
      batch_run_table$`Number of Occurrences (Previous)` <- NA
    }  
    
      out_tab <- batch_run_table %>% 
      DT::datatable(options = list(dom = 'tp',
                                   pageLength = 10,
                                   # columnDefs = list(list(width = "10%", className = 'dt-left', targets = c(1,2))),
                                   width = "100%"
      ),
      # filter = list(position = 'top'),
      selection = "none", 
      escape = FALSE, 
      rownames = FALSE
      ) 

      range_extent_trend_value <- purrr::map(batch_run_table$`Range Extent Change`, function(x){
        ifelse(!is.na(x), gsub("%", "", x) %>% as.numeric(), NA)
      }) %>% unlist()
      if (sum(range_extent_trend_value <= -30 & range_extent_trend_value >= -50, na.rm = TRUE) > 0) out_tab <- out_tab %>% formatStyle(columns = "Range Extent Change", backgroundColor = styleRow(rows = which(range_extent_trend_value <= -30 & range_extent_trend_value >= -50), values = "#fee08b"))
      if (sum(range_extent_trend_value <= -51 & range_extent_trend_value >= -70, na.rm = TRUE) > 0) out_tab <- out_tab %>%  formatStyle(columns = "Range Extent Change", backgroundColor = styleRow(rows = which(range_extent_trend_value <= -51 & range_extent_trend_value >= -70), values = "#fdae61"))
      if (sum(range_extent_trend_value <= -71, na.rm = TRUE) > 0) out_tab <- out_tab %>% formatStyle(columns = "Range Extent Change", backgroundColor = styleRow(rows = which(range_extent_trend_value <= -71), values = "#f46d43"))
      
      aoo_trend_value <- purrr::map(batch_run_table$`AOO Change`, function(x){
        ifelse(!is.na(x), gsub("%", "", x) %>% as.numeric(), NA)
      }) %>% unlist()
      if (sum(aoo_trend_value <= -30 & aoo_trend_value >= -50, na.rm = TRUE) > 0) out_tab <- out_tab %>% formatStyle(columns = "AOO Change", backgroundColor = styleRow(rows = which(aoo_trend_value <= -30 & aoo_trend_value >= -50), values = "#fee08b"))
      if (sum(aoo_trend_value <= -51 & aoo_trend_value >= -70, na.rm = TRUE) > 0) out_tab <- out_tab %>% formatStyle(columns = "AOO Change", backgroundColor = styleRow(rows = which(aoo_trend_value <= -51 & aoo_trend_value >= -70), values = "#fdae61"))
      if (sum(aoo_trend_value <= -71, na.rm = TRUE) > 0) out_tab <- out_tab %>% formatStyle(columns = "AOO Change", backgroundColor = styleRow(rows = which(aoo_trend_value <= -71), values = "#f46d43"))
      
      eo_count_trend_value <- purrr::map(batch_run_table$`Number of Occurrences Change`, function(x){
        ifelse(!is.na(x), gsub("%", "", x) %>% as.numeric(), NA)
      }) %>% unlist()
      if (sum(eo_count_trend_value <= -30 & eo_count_trend_value >= -50, na.rm = TRUE) > 0) out_tab <- out_tab %>% formatStyle(columns = "Number of Occurrences Change", backgroundColor = styleRow(rows = which(eo_count_trend_value <= -30 & eo_count_trend_value >= -50), values = "#fee08b"))
      if (sum(eo_count_trend_value <= -51 & eo_count_trend_value >= -70, na.rm = TRUE) > 0) out_tab <- out_tab %>% formatStyle(columns = "Number of Occurrences Change", backgroundColor = styleRow(rows = which(eo_count_trend_value <= -51 & eo_count_trend_value >= -70), values = "#fdae61"))
      if (sum(eo_count_trend_value <= -71, na.rm = TRUE) > 0) out_tab <- out_tab %>% formatStyle(columns = "Number of Occurrences Change", backgroundColor = styleRow(rows = which(eo_count_trend_value <= -71), values = "#f46d43"))

      out_tab
    
  })
  
  observeEvent(input$send_to_batch_mode, {
    
    # shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
    
    if (!is.null(batch_run_output$results)){
      
    updateTabsetPanel(inputId = "nav", selected = "MULTISPECIES MODE")

    batch_run_output$results[[batch_taxon_focus$taxon]]$info[1, ] <- taxon_data$info
    batch_run_output$results[[batch_taxon_focus$taxon]]$info_extended <- taxon_data$info_extended
    batch_run_output$results[[batch_taxon_focus$taxon]]$family <- taxon_data$family
    batch_run_output$results[[batch_taxon_focus$taxon]]$synonyms <- taxon_data$synonyms
    batch_run_output$results[[batch_taxon_focus$taxon]]$synonyms_selected <- taxon_data$synonyms_selected 
    batch_run_output$results[[batch_taxon_focus$taxon]]$gbif_occurrences_raw <- taxon_data$gbif_occurrences_raw
    batch_run_output$results[[batch_taxon_focus$taxon]]$gbif_occurrences <- taxon_data$gbif_occurrences
    batch_run_output$results[[batch_taxon_focus$taxon]]$uploaded_occurrences <- taxon_data$uploaded_occurrences
    batch_run_output$results[[batch_taxon_focus$taxon]]$all_occurrences <- taxon_data$all_occurrences
    batch_run_output$results[[batch_taxon_focus$taxon]]$shifted <- taxon_data$shifted
    batch_run_output$results[[batch_taxon_focus$taxon]]$sf <- taxon_data$sf
    batch_run_output$results[[batch_taxon_focus$taxon]]$sf_filtered <- taxon_data$sf_filtered
    batch_run_output$results[[batch_taxon_focus$taxon]]$filters_selected <- 
    list(
      clean_occ = input$clean_occ,
      centroid_filter = input$centroid_filter,
      date_start = input$year_filter[1],
      date_end = input$year_filter[2],
      months = input$seasonality,
      uncertainty_filter = input$uncertainty_filter,
      nations_filter = input$nation_filter,
      states_filter = input$states_filter,
      sources_filter = input$sources_filter,
      grid_cell_size = input$grid_cell_size,
      sep_distance = input$separation_distance
    )

    if (!is.null(input$batch_filedata_rank$datapath)){
      batch_rank_factor_file <- read.csv(input$batch_filedata_rank$datapath, header = TRUE) 
      if (taxon_data$info$scientificName %in% (batch_rank_factor_file[, 3] %>% as.character())){
        batch_rank_factor_sp <- batch_rank_factor_file[which((batch_rank_factor_file[, 3] %>% as.character()) %in% taxon_data$info$scientificName), ]
        batch_run_output$results[[batch_taxon_focus$taxon]]$rank_factor_comparison <- compare_rank_factors(taxon_data, rank_factor_upload = batch_rank_factor_sp)
      } else {
        batch_run_output$results[[batch_taxon_focus$taxon]]$rank_factor_comparison <- compare_rank_factors(taxon_data)
      }
    } else {
      batch_run_output$results[[batch_taxon_focus$taxon]]$rank_factor_comparison <- compare_rank_factors(taxon_data)
    }
    
    # batch_run_output$results[[batch_taxon_focus$taxon]]$rank_factor_comparison <- compare_rank_factors(taxon_data)
    
    batch_run_output$results[[batch_taxon_focus$taxon]]$temporal_change <- calculate_rarity_change(
      taxon_data = taxon_data,
      period1 = input$period1,
      period2 = input$period2,
      period3 = input$period3, 
      aoo_grid_cell_size = input$grid_cell_size, 
      occ_sep_distance = input$separation_distance
    )
    
    updated_batch_table <- data.frame(
      taxon = taxon_data$info$scientificName,
      total_observations_used = nrow(taxon_data$sf_filtered),
      range_value = paste0(taxon_data$species_range_value, " (", taxon_data$species_range_factor, ")"),
      current_range_value = purrr::map(1, function(out){
        if (!is.na(batch_run_output$results[[batch_taxon_focus$taxon]]$rank_factor_comparison$new_previous_species_range_comparison)){
          x <- batch_run_output$results[[batch_taxon_focus$taxon]]$rank_factor_comparison$new_previous_species_range_comparison
          if (x == "equal") res <- as.character(icon("equals", "fa-2x", style = "color: #BEBEBE;"))
          if (x == "lower") res <- as.character(icon("arrow-down", "fa-2x", style = "color: #ef8a62;"))
          if (x == "higher") res <- as.character(icon("arrow-up", "fa-2x", style = "color: #67a9cf;"))
        } else {
          res <- NA
        }
        res
      }) %>% unlist(),
      range_value_trend = purrr::map_chr(1, function(i){
        out <- ifelse(!is.na(batch_run_output$results[[batch_taxon_focus$taxon]]$temporal_change$eoo_change[2]), batch_run_output$results[[batch_taxon_focus$taxon]]$temporal_change$eoo_change[2], NA)
        paste0(out, "%")
      }) %>% unlist(),
      AOO_value = paste0(taxon_data$AOO_value, " (", taxon_data$AOO_factor, ")"),
      current_AOO_value = purrr::map(1, function(out){
        if (!is.na(batch_run_output$results[[batch_taxon_focus$taxon]]$rank_factor_comparison$new_previous_aoo_comparison)){
          x <- batch_run_output$results[[batch_taxon_focus$taxon]]$rank_factor_comparison$new_previous_aoo_comparison
          if (x == "equal") res <- as.character(icon("equals", "fa-2x", style = "color: #BEBEBE;"))
          if (x == "lower") res <- as.character(icon("arrow-down", "fa-2x", style = "color: #ef8a62;"))
          if (x == "higher") res <- as.character(icon("arrow-up", "fa-2x", style = "color: #67a9cf;"))
        } else {
          res <- NA
        }
        res
      }) %>% unlist(),
      AOO_value_trend = purrr::map(1, function(i){
        out <- ifelse(!is.na(batch_run_output$results[[batch_taxon_focus$taxon]]$temporal_change$aoo_change[2]), batch_run_output$results[[batch_taxon_focus$taxon]]$temporal_change$aoo_change[2], NA)
        paste0(out, "%")
      }) %>% unlist(),
      EOcount_value = paste0(taxon_data$EOcount_value, " (", taxon_data$EOcount_factor, ")"),
      current_EOcount_value = purrr::map(1, function(out){
        if (!is.na(batch_run_output$results[[batch_taxon_focus$taxon]]$rank_factor_comparison$new_previous_eocount_comparison)){
          x <- batch_run_output$results[[batch_taxon_focus$taxon]]$rank_factor_comparison$new_previous_eocount_comparison
          if (x == "equal") res <- as.character(icon("equals", "fa-2x", style = "color: #BEBEBE;"))
          if (x == "lower") res <- as.character(icon("arrow-down", "fa-2x", style = "color: #ef8a62;"))
          if (x == "higher") res <- as.character(icon("arrow-up", "fa-2x", style = "color: #67a9cf;"))
        } else {
          res <- NA
        }
        res
      }) %>% unlist(),
      EOcount_value_trend = purrr::map(1, function(i){
        out <- ifelse(!is.na(batch_run_output$results[[batch_taxon_focus$taxon]]$temporal_change$eo_count_change[2]), batch_run_output$results[[batch_taxon_focus$taxon]]$temporal_change$eo_count_change[2], NA)
        paste0(out, "%")
      }) %>% unlist(),
      Reviewed = TRUE
    )

    if (updated_batch_table$taxon %in% batch_run_output$table$taxon){
      batch_run_output$table <- batch_run_output$table %>%
        dplyr::filter(taxon != updated_batch_table$taxon) %>%
        rbind(
          updated_batch_table
        )
    } else {
      batch_run_output$table <- batch_run_output$table %>% 
        rbind(
          updated_batch_table
        )
    }
    } else {
      sendSweetAlert(session, type = "warning", title = "Oops!", text = "You need to run a multispecies analysis before you are able to send data back to it!", closeOnClickOutside = TRUE)
    }
    
    # shinybusy::remove_modal_spinner()
 
  })
  
  observeEvent(input$select_button, {
    
    # shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
    
    updateTabsetPanel(session = session, inputId = "nav", selected = "SINGLE SPECIES MODE")
    
    Sys.sleep(1)
    # selectedRow <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
    # batch_taxon_focus$taxon <- batch_run_output$table[selectedRow, ]$taxon
    
    batch_taxon_focus$taxon <-  str_split_1(input$select_button, "_")[1]
    
    # updateTextInput(inputId = "search_taxon", value = batch_taxon_focus$taxon)
    
    taxon_data$info <- batch_run_output$results[[batch_taxon_focus$taxon]]$info[1, ]
    taxon_data$info_extended <- batch_run_output$results[[batch_taxon_focus$taxon]]$info_extended
    taxon_data$family <- batch_run_output$results[[batch_taxon_focus$taxon]]$family
    taxon_data$synonyms <- batch_run_output$results[[batch_taxon_focus$taxon]]$synonyms
    taxon_data$synonyms_selected <- batch_run_output$results[[batch_taxon_focus$taxon]]$synonyms_selected
    taxon_data$gbif_occurrences_raw <- batch_run_output$results[[batch_taxon_focus$taxon]]$gbif_occurrences_raw
    taxon_data$gbif_occurrences <- batch_run_output$results[[batch_taxon_focus$taxon]]$gbif_occurrences
    taxon_data$uploaded_occurrences <- batch_run_output$results[[batch_taxon_focus$taxon]]$uploaded_occurrences
    taxon_data$all_occurrences <- batch_run_output$results[[batch_taxon_focus$taxon]]$all_occurrences
    taxon_data$shifted <- batch_run_output$results[[batch_taxon_focus$taxon]]$shifted
    taxon_data$sf <- batch_run_output$results[[batch_taxon_focus$taxon]]$sf
    taxon_data$sf_filtered <- batch_run_output$results[[batch_taxon_focus$taxon]]$sf_filtered
    
    # taxon_data$species_range_value <- batch_run_output$results[[batch_taxon_focus$taxon]]$species_range_value
    # taxon_data$species_range_map <- batch_run_output$results[[batch_taxon_focus$taxon]]$species_range_map
    # taxon_data$AOO_value <- batch_run_output$results[[batch_taxon_focus$taxon]]$AOO_value
    # taxon_data$AOO_map <- batch_run_output$results[[batch_taxon_focus$taxon]]$AOO_map
    # taxon_data$EOcount_map <- batch_run_output$results[[batch_taxon_focus$taxon]]$EOcount_map
    # taxon_data$EOcount_value <- batch_run_output$results[[batch_taxon_focus$taxon]]$EOcount_value
    taxon_data$info$scientificName <- batch_taxon_focus$taxon
    selected_taxon$name <- batch_taxon_focus$taxon[1]
    
    # updateTextInput(session = session, inputId = "number_gbif_occurrences", label = "", value = nrow(taxon_data$all_occurrences))

    shinyjs::show(id = "data_panel")
    
    # shinybusy::remove_modal_spinner()
    
  })
  
  output$download_rank_data_batch <- downloadHandler(
    filename = function() {
      paste("RARECAT-multispecies_run-rank_factor_values-", Sys.Date(), ".xlsx", sep="")
    },
    content = function(file) {
      if (!is.null(input$batch_filedata_rank$datapath)){
        batch_rank_factor_file <- read.csv(input$batch_filedata_rank$datapath, header = TRUE) 
        # names(batch_rank_factor_file) <- gsub("\\.", " ", names(batch_rank_factor_file))
      } else {
        batch_rank_factor_file <- NULL
      }
      batch_out_tab1 <- purrr::map(1:length(batch_run_output$results), function(i){
        taxon_data <- batch_run_output$results[[i]]
        out <- matrix(nrow = 1, ncol = 42, data = "") %>%
          as.data.frame()
        names(out)[c(1:3, 6:9, 11, 13, 15, 16:17, 19:20, 22:28, 30, 32:42)] <- c("Calc Rank", "Assigned Rank", "Species or Community Scientific Name*", "Element ID",
                                                                                 "Elcode*", "Common Name*", "Classification*", "Range Extent", "Area of Occup 4-km2 grid cells",
                                                                                 "# Occur", "Pop Size", "# Occur Good Viab", "Environm Specif (opt.)", "Overall Threat Impact",
                                                                                 "Intrinsic Vulner (opt.)", "Short-term Trend", "Long-term Trend", "Rank Adjustment Reasons", "Assigned Rank Reasons", "Rank Factors Author", 
                                                                                 "Rank Factors Date", "Rank Review Date", "Range Extent Comments", "Area of Occupancy Comments",
                                                                                 "# of Occurrences Comments", "Population Size Comments", "Good Viability/Integrity Comments", "Environmental Specificity Comments",
                                                                                 "Threat Impact Comments", "Threat Impact Adjustment Reasons", "Intrinsic Vulnerability Comments", "Short-term Trend Comments", "Long-term Trend Comments"
        )
    if (!is.null(batch_rank_factor_file)){
      batch_rank_factor_sp <- batch_rank_factor_file %>% 
        dplyr::filter(Species.or.Community.Scientific.Name. == names(batch_run_output$results)[i])
      if (nrow(batch_rank_factor_sp) > 0){
        out <- batch_rank_factor_sp
        names(out)[c(1:3, 6:9, 11, 13, 15, 16:17, 19:20, 22:28, 30, 32:42)] <- c("Calc Rank", "Assigned Rank", "Species or Community Scientific Name*", "Element ID",
                                                                                 "Elcode*", "Common Name*", "Classification*", "Range Extent", "Area of Occup 4-km2 grid cells",
                                                                                 "# Occur", "Pop Size", "# Occur Good Viab", "Environm Specif (opt.)", "Overall Threat Impact",
                                                                                 "Intrinsic Vulner (opt.)", "Short-term Trend", "Long-term Trend", "Rank Adjustment Reasons", "Assigned Rank Reasons", "Rank Factors Author", 
                                                                                 "Rank Factors Date", "Rank Review Date", "Range Extent Comments", "Area of Occupancy Comments",
                                                                                 "# of Occurrences Comments", "Population Size Comments", "Good Viability/Integrity Comments", "Environmental Specificity Comments",
                                                                                 "Threat Impact Comments", "Threat Impact Adjustment Reasons", "Intrinsic Vulnerability Comments", "Short-term Trend Comments", "Long-term Trend Comments"
        )
        if (!is.null(taxon_data$species_range_value)){
          out[, 11] <- taxon_data$species_range_factor
        }
        if (!is.null(taxon_data$AOO_value)){
            out[, 13] <- taxon_data$AOO_factor
        }
        if (!is.null(taxon_data$EOcount_value)){
          out[, 15] <- taxon_data$EOcount_factor
        }
      } else {
        batch_rank_factor_sp <- NULL
      }
    } else {
      batch_rank_factor_sp <- NULL
      out[, 3] <-  ifelse(!is.null(taxon_data$info), taxon_data$info$scientificName, "")
      if (!is.null(taxon_data$species_range_value)){
        # out[, 11] <- cut(as.numeric(taxon_data$species_range_value), breaks = c(0, 0.999, 99.999, 249.999, 999.999, 4999.999, 19999.999, 199999.999, 2499999.999, 1000000000), labels = c("Z", LETTERS[1:8]))
        out[, 11] <- taxon_data$species_range_factor
      }
      if (!is.null(taxon_data$AOO_value)){
        out[, 13] <- taxon_data$AOO_factor
      }
      if (!is.null(taxon_data$EOcount_value)){
        # out[, 15] <- cut(as.numeric(taxon_data$EOcount_value), breaks = c(0, 0.999, 5.999, 19.999, 79.999, 299.999, 1000000000), labels = c("Z", LETTERS[1:5]))
        out[, 15] <- taxon_data$EOcount_factor
      }
      if (input$batch_assessment_type == "global"){
      out[, 2] <-  ifelse(!is.null(taxon_data$info), taxon_data$info$roundedGRank, "")
      out[, 6] <-  ifelse(!is.null(taxon_data$info), taxon_data$info$elementGlobalId, "")
      out[, 7] <-  ifelse(!is.null(taxon_data$info), taxon_data$info$elcode, "")
      out[, 8] <-  ifelse(!is.null(taxon_data$info), taxon_data$info$primaryCommonName, "")
      out[, 9] <-  ifelse(!is.null(taxon_data$info_extended), taxon_data$info_extended$classificationStatus$classificationStatusDescEn, "")

      if (!is.null(taxon_data$info_extended$rankInfo$popSize)){
        rank_def <- rank_factor_definitions$population_size_description[grep(taxon_data$info_extended$rankInfo$popSize$popSizeDescEn, rank_factor_definitions$population_size_description, fixed = TRUE)]
        out[, 16] <- rank_factor_definitions$population_size_value[rank_factor_definitions$population_size_description == rank_def]
      }
      if (!is.null(taxon_data$info_extended$rankInfo$numberGoodEos$numberGoodEosDescEn)){
        rank_def <- rank_factor_definitions$good_viability_description[grep(taxon_data$info_extended$rankInfo$numberGoodEos$numberGoodEosDescEn, rank_factor_definitions$good_viability_description, fixed = TRUE)]
        out[, 17] <- rank_factor_definitions$good_viability_value[rank_factor_definitions$good_viability_description == rank_def]
      }
      if (!is.null(taxon_data$info_extended$rankInfo$enviromentalSpecificity)){
        rank_def <- rank_factor_definitions$environmental_specificity_description[grep(gsub(" |\\.", "", taxon_data$info_extended$rankInfo$enviromentalSpecificity$enviromentalSpecificityDescEn), gsub(" |\\.", "", rank_factor_definitions$environmental_specificity_description), fixed = TRUE)]
        out[, 19] <- rank_factor_definitions$environmental_specificity_value[rank_factor_definitions$environmental_specificity_description == rank_def]
      }
      if (!is.null(taxon_data$info_extended$rankInfo$threatImpactAssigned)){
        rank_def <- rank_factor_definitions$threat_impact_description[grep(taxon_data$info_extended$rankInfo$threatImpactAssigned$threatImpactAssignedDescEn, rank_factor_definitions$threat_impact_description, fixed = TRUE)]
        out[, 20] <- rank_factor_definitions$threat_impact_value[rank_factor_definitions$threat_impact_description == rank_def]
      }
      if (!is.null(taxon_data$info_extended$rankInfo$intrinsicVulnerability)){
        rank_def <- rank_factor_definitions$intrinsic_vulnerability_description[grep(taxon_data$info_extended$rankInfo$intrinsicVulnerability$intrinsicVulnerabilityDescEn, rank_factor_definitions$intrinsic_vulnerability_description, fixed = TRUE)]
        out[, 22] <- rank_factor_definitions$intrinsic_vulnerability_value[rank_factor_definitions$intrinsic_vulnerability_description == rank_def]
        
      }
      if (!is.null(taxon_data$info_extended$rankInfo$shortTermTrend)){
        rank_def<- rank_factor_definitions$temporal_trend_description[grep(gsub(" ", "", taxon_data$info_extended$rankInfo$shortTermTrend$shortTermTrendDescEn), gsub(" ", "", rank_factor_definitions$temporal_trend_description), fixed = TRUE)]
        out[, 23] <- rank_factor_definitions$temporal_trend_value[rank_factor_definitions$temporal_trend_description == rank_def]
        
      }
      if (!is.null(taxon_data$info_extended$rankInfo$longTermTrend)){
        rank_def <- rank_factor_definitions$temporal_trend_description[grep(gsub(" ", "", taxon_data$info_extended$rankInfo$longTermTrend$longTermTrendDescEn), gsub(" ", "", rank_factor_definitions$temporal_trend_description), fixed = TRUE)]
        out[, 24] <- rank_factor_definitions$temporal_trend_value[rank_factor_definitions$temporal_trend_description == rank_def]
        
      }
      out[, 25] <- ifelse(!is.null(taxon_data$info_extended$grankAdjustmentReasons), taxon_data$info_extended$grankAdjustmentReasons, "")
      out[, 26] <- ifelse(!is.null(taxon_data$info_extended$grankReasons), taxon_data$info_extended$grankReasons, "")
      out[, 27] <- ifelse(!is.null(taxon_data$info_extended$conservationStatusFactorsEditionAuthors), taxon_data$info_extended$conservationStatusFactorsEditionAuthors, "")
      out[, 28] <- ifelse(!is.null(taxon_data$info_extended$conservationStatusFactorsEditionDate), taxon_data$info_extended$conservationStatusFactorsEditionDate, "")
      out[, 30] <- ifelse(!is.null(taxon_data$info_extended$grankReviewDate), taxon_data$info_extended$grankReviewDate, "")
      out[, 32] <- ifelse(!is.null(taxon_data$info_extended$rankInfo$rangeExtentComments), taxon_data$info_extended$rankInfo$rangeExtentComments, "")
      out[, 33] <- ifelse(!is.null(taxon_data$info_extended$rankInfo$areaOfOccupancyComments), taxon_data$info_extended$rankInfo$areaOfOccupancyComments, "")
      out[, 34] <- ifelse(!is.null(taxon_data$info_extended$rankInfo$numberEosComments), taxon_data$info_extended$rankInfo$numberEosComments, "")
      out[, 35] <- ifelse(!is.null(taxon_data$info_extended$rankInfo$popSizeComments), taxon_data$info_extended$rankInfo$popSizeComments, "")
      out[, 36] <- ifelse(!is.null(taxon_data$info_extended$rankInfo$viabilityComments), taxon_data$info_extended$rankInfo$viabilityComments, "")
      out[, 37] <- ifelse(!is.null(taxon_data$info_extended$rankInfo$enviromentalSpecificityComments), taxon_data$info_extended$rankInfo$enviromentalSpecificityComments, "")
      out[, 38] <- ifelse(!is.null(taxon_data$info_extended$rankInfo$threatImpactComments), taxon_data$info_extended$rankInfo$threatImpactComments, "")
      out[, 39] <- ifelse(!is.null(taxon_data$info_extended$rankInfo$grankAdjustmentReasons), taxon_data$info_extended$rankInfo$grankAdjustmentReasons, "")
      out[, 40] <- ifelse(!is.null(taxon_data$info_extended$rankInfo$intrinsicVulnerabilityComments), taxon_data$info_extended$rankInfo$intrinsicVulnerabilityComments, "")
      out[, 41] <- ifelse(!is.null(taxon_data$info_extended$rankInfo$shortTermTrendComments), taxon_data$info_extended$rankInfo$shortTermTrendComments, "")
      out[, 42] <- ifelse(!is.null(taxon_data$info_extended$rankInfo$longTermTrendComments), taxon_data$info_extended$rankInfo$longTermTrendComments, "")
      }
    }
      out
    }) %>% bind_rows()
      
    batch_out_tab2 <- purrr::map(1:length(batch_run_output$results), function(i){
        taxon_data <- batch_run_output$results[[i]]
        out2_names <- c("NatureServe accepted name", "NatureServe synonyms", "GBIF taxonomic concepts with GBIF IDs", "EGT ID", "EGT UID", "ELCODE", 
                        "Assessment Type", "Nations included", "Subnations included", "New Range Extent value (sq km)", "New Range Extent letter",
                        "Previous Range Extent letter", "Compare Range Extent letter (new vs. previous)", "New Area of Occupancy grid cell size",
                        "New Area of Occupancy value", "New Area of Occupancy letter", "Previous Area of Occupancy letter", "Compare Area of Occupancy letter (new vs. previous)",
                        "RARECAT occurrence separation distance (m)", "New Number of Occurrences value", "New Number of Ocurrences letter", "Previous Number of Occurrences letter",
                        "Compare Number of Occurrences letter (new vs. previous)", "Time Period 1", "Time Period 2", "Range Extent temporal change analysis", 
                        "Area of Occupancy temporal change analysis", "Number of Occurrences temporal change analysis", "Number of records included",
                        "Data sources included", "Date range of records included", "Months of records included", "Locational uncertainty cutoff", "Other filters", "New assessment date"
        )
        
        out2 <- matrix(nrow = 1, ncol = length(out2_names), data = "") %>%
          as.data.frame()
        names(out2) <- out2_names
        
        if (input$batch_assessment_type == "global"){
          out2[, 4] <- taxon_data$info_extended$elementGlobalId
          out2[, 5] <- taxon_data$info_extended$uniqueId
          out2[, 6] <- taxon_data$info_extended$elcode
        }
        out2[, 1] <- ifelse(!is.null(taxon_data$info), taxon_data$info$scientificName, "")
        out2[, 2] <- ifelse(length(taxon_data$info$synonyms %>% unlist() %>% na.omit() %>% as.character()) > 0, paste0(taxon_data$info$synonyms %>% unlist() %>% na.omit() %>% as.character(), collapse = "; "), "")
        out2[, 3] <- ifelse(!is.null(taxon_data$synonyms$scientificName), 
                            paste0(taxon_data$synonyms$scientificName %>% na.omit() %>% as.character(), collapse = "; "), 
                            ""
        )
        out2[, 7] <- input$batch_assessment_type
        out2[, 8] <- ifelse(!is.null(input$nation_filter), paste0(input$nation_filter, collapse = "; "), "")
        out2[, 9] <- ifelse(!is.null(input$states_filter), paste0(input$states_filter, collapse = "; "), "")
        out2[, 10] <- ifelse(!is.null(taxon_data$species_range_value), taxon_data$species_range_value %>% as.character(), "")
        out2[, 11] <- ifelse(!is.null(taxon_data$species_range_factor), taxon_data$species_range_factor %>% as.character(), "")
        out2[, 14] <- input$grid_cell_size %>% as.numeric()
        out2[, 15] <- ifelse(!is.null(taxon_data$AOO_value), taxon_data$AOO_value %>% as.character(), "")
        out2[, 16] <- ifelse(!is.null(taxon_data$AOO_factor), taxon_data$AOO_factor %>% as.character(), "")
        out2[, 19] <- input$separation_distance %>% as.numeric()
        out2[, 20] <- ifelse(!is.null(taxon_data$EOcount_value), taxon_data$EOcount_value %>% as.character(), "")
        out2[, 21] <- ifelse(!is.null(taxon_data$EOcount_factor), taxon_data$EOcount_factor %>% as.character(), "")
        out2[, 24] <- paste0(input$period1, collapse = " - ")
        out2[, 25] <- paste0(input$period2, collapse = " - ")
        out2[, 26] <- taxon_data$temporal_change$eoo_change[2]
        out2[, 27] <- taxon_data$temporal_change$aoo_change[2]
        out2[, 28] <- taxon_data$temporal_change$eo_count_change[2]
        out2[, 29] <- nrow(taxon_data$sf_filtered)
        out2[, 30] <- paste0(input$sources_filter, collapse = "; ")
        out2[, 31] <- paste0(input$year_filter, collapse = "; ")
        out2[, 32] <- paste0(input$seasonality, collapse = "; ")  
        out2[, 33] <- input$uncertainty_filter
        out2[, 34] <- paste0(
          c(ifelse(input$clean_occ, "Basic GBIF data cleanup implemented", ""),
            ifelse(input$clean_occ, "Putative centroids removed", ""),
            paste0("Data types included: ", input$type_filter),
            ifelse(length(taxon_data$sf_filtered$EORANK %>% complete.cases(.)) > 0, paste0("Element occurrence ranks included: ", paste0(unique(taxon_data$sf_filtered$EORANK), collapse = "|")), "")
          ), collapse = "; ")
        out2[, 35] <- Sys.Date()
        taxon_data$rank_factor_comparison <- data.frame(
          previous_species_range_letter = NA,
          new_previous_species_range_comparison = NA,
          previous_aoo_letter = NA,
          new_previous_aoo_comparison = NA,
          previous_eocount_letter = NA,
          new_previous_eocount_comparison = NA
        )
        if (!is.null(batch_rank_factor_file)){
          batch_rank_factor_sp <- batch_rank_factor_file %>%
            dplyr::filter(Species.or.Community.Scientific.Name. == names(batch_run_output$results)[i])
          if (nrow(batch_rank_factor_sp) > 0){
            taxon_data$rank_factor_comparison <- compare_rank_factors(taxon_data, rank_factor_upload = batch_rank_factor_sp)
          }
        }
        if (input$batch_assessment_type == "global" & is.null(batch_rank_factor_file)){
          taxon_data$rank_factor_comparison <- compare_rank_factors(taxon_data)
        }
        out2[, 12] <- taxon_data$rank_factor_comparison$previous_species_range_letter %>% as.character()
        out2[, 13] <- taxon_data$rank_factor_comparison$new_previous_species_range_comparison %>% as.character()
        out2[, 17] <- taxon_data$rank_factor_comparison$previous_aoo_letter %>% as.character()
        out2[, 18] <- taxon_data$rank_factor_comparison$new_previous_aoo_comparison %>% as.character()
        out2[, 22] <- taxon_data$rank_factor_comparison$previous_eocount_letter %>% as.character()
        out2[, 23] <- taxon_data$rank_factor_comparison$new_previous_eocount_comparison %>% as.character()
      out2
      
    }) %>% bind_rows()
      
  writexl::write_xlsx(
    list(
      "Rank Calculator" = batch_out_tab1,
      "Assessment Details" = batch_out_tab2
    ), 
    path = file
  )

})
  
  observeEvent(input$batch_clear, {
    
    # updateTabsetPanel(inputId = "nav", selected = "MULTISPECIES MODE")
    reset('batch_filedata_rank')
    reset('batch_filedata_obs')
    updateSelectizeInput(session = session, inputId = 'batch_assessment_type', selected = "global")
    updateTextInput(session = session, inputId = "typed_list", value = "")
    shinyjs::hide("batch_output")
    batch_run_output <- reactiveValues(
      results = NULL,
      table = NULL
    )
    updateCollapse(session = session, id = "batch_parameters", open = "Rank assessment parameters")
    
  })
  
}
