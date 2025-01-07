#' ---
#' title: NatureServe RARECAT: Rapid Assessment of Rarity and Endangerment Conservation Assessment Tool
#' ---
#'
#' # Server setup
#' ## Load libraries
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
#'
#' ## Load NatureServe Network subnation polygons for subnation overlay
network_polys <- readRDS("data/subnation_polys.rds")
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
      leaflet::addProviderTiles(providers$Esri.WorldTerrain, group = "Esri World Terrain", options = list(pathOptions(pane = "basemap1"))) %>%
      leaflet::addMapPane("basemap2", zIndex = -100) %>% # Add basemap 2
      leaflet::addProviderTiles(providers$Esri.WorldImagery, group = "Esri World Imagery", options = list(pathOptions(pane = "basemap2"))) %>%
      leaflet::addMapPane("basemap3", zIndex = -100) %>% # Add basemap 3
      leaflet::addProviderTiles(providers$OpenStreetMap, group = "Open Street Map", options = list(pathOptions(pane = "basemap3"))) %>%
      leaflet::addMapPane("basemap4", zIndex = -100) %>% # Add basemap 4
      leaflet::addProviderTiles(providers$Esri.WorldStreetMap, group = "Esri World Street Map", options = list(pathOptions(pane = "basemap4"))) %>%
      leaflet::addScaleBar(position = "bottomleft") %>% # Add scale bar
      leaflet.extras::addResetMapButton() %>% # Add button to reset map bounds
      leaflet.extras::addSearchOSM() %>% # Add functionality to search for specific location using Open Street Map
      leaflet::addLayersControl(baseGroups = c("Esri World Street Map", "Open Street Map", "Esri World Terrain", "Esri World Imagery"), # Add layers control widget
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
  minimum_fields <- c("key", "scientificName", "prov", "longitude", "latitude", "coordinateUncertaintyInMeters", "stateProvince", "countryCode", "year", "institutionCode", "references")
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
    AOO_value = NULL,
    AOO_map = NULL,
    EOcount = NULL,
    EOcount_value = NULL
  )
  #' ### Object to store all clicked point IDs
  clicks <- reactiveValues(IDs = vector(mode = "character"))
  #' ### Object to store "Begin assessment" button presses
  assessment_start_button_presses <- reactiveValues(values = 0)
  #'
  #' ## Set up reactive expressions
  #' ### Load NatureServe API information for taxon entered in search bar
  taxon_NS_options <- reactive({
    
    ns_table <- natserv::ns_search_spp(text_adv = list(searchToken = input$search_taxon, matchAgainst = "allScientificNames", operator="contains"))$results
    gbif_table <- rgbif::name_suggest(q = input$search_taxon, rank = c("species", "subspecies"), limit = 10)$data
    
    out <- NULL
    
    if (nrow(ns_table) > 0){
      ns_table <- ns_table %>% 
        dplyr::mutate(Source = "NatureServe", synonyms = ns_table$speciesGlobal$synonyms) %>% 
        dplyr::select(scientificName, Source, elementGlobalId, primaryCommonName, roundedGRank, elcode, uniqueId, synonyms)
      out <- rbind(out, ns_table)
    } 

      if (nrow(gbif_table) > 0){
        gbif_table <- gbif_table %>% 
          dplyr::rename(scientificName = canonicalName, elementGlobalId = key) %>% 
          dplyr::mutate(Source = "GBIF", primaryCommonName = NA, roundedGRank = NA, elcode = NA, synonyms = NA, uniqueId = NA) %>% 
          dplyr::select(scientificName, Source, elementGlobalId, primaryCommonName, roundedGRank, elcode, uniqueId, synonyms)
      out <- rbind(out, gbif_table)
      }
    
    out
    
  })
  
  download_gbif_data <- reactive({
    
    if (!is.null(selected_taxon$name)) {
      
      gbif_download <- get_gbif_data(sp_data = taxon_data$synonyms_selected, number_observations = input$number_gbif_occurrences)
      
      if (nrow(gbif_download$sp_occurrences) > 0){
        
        gbif_download
        
      } else {
        
        shinybusy::remove_modal_spinner() # show the modal window
        sendSweetAlert(session, type = "warning", title = "Oops!", text = "There are no valid occurrences on GBIF for this taxon or any of its synonyms recognized by NatureServe", closeOnClickOutside = TRUE)
      }
      
    } else {
      shinybusy::remove_modal_spinner() # show the modal window
      sendSweetAlert(session, type = "warning", title = "Oops!", text = "You need to select a taxon before loading GBIF data!", closeOnClickOutside = TRUE)
    }
    
  })
  
  #' ### Load user-uploaded data
  uploaded_data <- reactive({
    
    
    out <- purrr::map(input$filedata$datapath, process_user_data, minimum_fields = minimum_fields) %>% 
      dplyr::bind_rows() 
    
    out 
  })
  #'
  #' ## Specify what happens when user enters text in search bar
  observeEvent(input$search_taxon, {
    
    if (input$search_taxon != ""){
      
      shinyjs::show(id = "taxon_search_panel")
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
                paste0("<a href='", paste0("https://www.gbif.org/species", elementGlobalId[i]), "' target='_blank'>", elementGlobalId[i], "</a>")
              } else {
                elementGlobalId
              }
            }) %>% unlist()) %>% 
            dplyr::select(-synonyms, -elcode, -uniqueId) %>% 
            dplyr::rename("Scientific name" = scientificName, 
                          "ID" = elementGlobalId, 
                          "Common name" = primaryCommonName, 
                          "G rank" = roundedGRank
                          ) %>% 

            DT::datatable(options = list(dom = 't', pageLength = 100, autoWidth = TRUE), selection = list(mode = 'single', target = 'row'), escape = FALSE, rownames = FALSE)
          
        })
        
      }
      
    } else {
      shinyjs::hide(id = "taxon_search_panel")
    }
    
  })
  
  observeEvent({
    input$taxon_NS_table_rows_selected
  }, {
    
    shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
    
    updateTextInput(inputId = "search_taxon", value = "")
    
    taxon_data$info <- taxon_NS_options()[input$taxon_NS_table_rows_selected, ]
    
    if (taxon_data$info$Source == "NatureServe"){
      taxon_data$info_extended <- natserv::ns_id(uid = taxon_data$info$uniqueId)
      selected_taxon$name <- c(taxon_data$info$scientificName, unlist(taxon_data$info$synonyms))
      selected_taxon$name <- gsub("ssp. |var. ", "", selected_taxon$name)
    } else if (taxon_data$info$Source == "GBIF"){
      selected_taxon$name <- c(taxon_data$info$scientificName)
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
    
    gbif_counts <- purrr::map_dbl(taxon_data$synonyms$key, function(x) rgbif::occ_count(taxonKey = x, hasCoordinate=TRUE))
    
    taxon_data$synonyms <- taxon_data$synonyms %>% 
      dplyr::mutate(occurrence_count = gbif_counts) %>% 
      dplyr::filter(occurrence_count > 0)
    
    shinyInput <- function(FUN, len, id, ...) {
      inputs <- character(len)
      for (i in seq_len(len)) {
        inputs[i] <- as.character(FUN(paste0(id, i), label = NULL, ...))
      }
      inputs
    }

    output$taxon_options_table <- DT::renderDataTable({
      
      taxon_data$synonyms %>% 
        dplyr::mutate(key = paste0("<a href='", paste0("https://www.gbif.org/species/", key), "' target='_blank'>", key, "</a>")) %>%         
        dplyr::select(scientificName, key, occurrence_count) %>% 
        dplyr::rename("Scientific name" = scientificName, 
                      "Key" = key,
                      "Number of occurrences" = occurrence_count
        ) %>% 
        DT::datatable(options = list(dom = 't', pageLength = 100, autoWidth = TRUE,
                                     columnDefs = list(list(width = '400px', targets = c(0)))
        ), 
        selection = list(mode = 'multiple', target = 'row'), 
        escape = FALSE, 
        rownames = FALSE
        )
      
    })
    
    DT::dataTableProxy("taxon_NS_table") %>% 
      selectRows(selected = NULL)
    
    shinyjs::show(id = "taxon_options_panel")
  
    
    m <- leafletProxy("main_map") %>%
      clearShapes() %>% 
      clearMarkers() %>% 
      clearMarkerClusters() 
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
    taxon_data$EOcount <- NULL
    taxon_data$EOcount_value <- NULL

    updateMaterialSwitch(session = session, inputId = "load_gbif_data", value = FALSE)
    updateMaterialSwitch(session = session, inputId = "map_uploads", value = FALSE)    
    updateMaterialSwitch(session = session, inputId = "range_extent", value = FALSE)
    updateMaterialSwitch(session = session, inputId = "area_of_occupancy", value = FALSE)
    updateMaterialSwitch(session = session, inputId = "number_EOs", value = FALSE)
    
    shinyjs::show(id = "analysis_panel")
    shinyjs::hide(id = "data_panel")
    
    shinybusy::remove_modal_spinner()
    
  })
  
  observeEvent({
    input$begin_assessment
    input$taxon_options_table_rows_selected
  }, {
    
    assessment_start_button_presses$values <- c(assessment_start_button_presses$values, input$begin_assessment)
    
    if (assessment_start_button_presses$values[length(assessment_start_button_presses$values)] != assessment_start_button_presses$values[length(assessment_start_button_presses$values)-1]){
      taxon_data$synonyms_selected <- taxon_data$synonyms[input$taxon_options_table_rows_selected, ]
      shinyjs::hide(id = "taxon_options_panel")
      
      updateTextInput(session = session, inputId = "number_gbif_occurrences", label = "", value = sum(taxon_data$synonyms_selected$occurrence_count))
      shinyjs::show(id = "load_data_panel")
      
    }
    
    
  })
  
  observeEvent({
    input$load_gbif_data
    input$clean_occ
    input$centroid_filter
    }, {
    
    if (isTRUE(input$load_gbif_data)){
      
      if (is.null(taxon_data$gbif_occurrences_raw)){
        
        shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
        
        gbif_download <- download_gbif_data()
        taxon_data$gbif_occurrences_raw <- gbif_download$sp_occurrences
        taxon_data$shifted <- gbif_download$shifted
      }

      taxon_data$gbif_occurrences <- taxon_data$gbif_occurrences_raw %>% 
        clean_gbif_data(clean = input$clean_occ, remove_centroids = input$centroid_filter, minimum_fields = minimum_fields)

      shinyjs::show(id = "data_panel")
      
      shinybusy::remove_modal_spinner()
      
    } 
    
  })
  
  observeEvent({
    input$filedata
    input$map_uploads
  }, {
    
    if (isTRUE(input$map_uploads)){
      
      if (!is.null(input$filedata)){
        
        shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
        
        taxon_data$uploaded_occurrences <- uploaded_data()
        taxon_data$uploaded_occurrences <- taxon_data$uploaded_occurrences %>% 
          dplyr::mutate(key = paste(prov, 1:nrow(taxon_data$uploaded_occurrences), sep = "_"))

        selected_taxon$NS <- c(selected_taxon$NS, unique(taxon_data$uploaded_occurrences$scientificName)) %>% unique()
        
        if (taxon_data$info$scientificName == "New taxon"){
          taxon_data$info <- data.frame(
            scientificName = taxon_data$uploaded_occurrences$scientificName[1],
            elementGlobalId = NA, Source = "Uploaded", primaryCommonName = NA, roundedGRank = NA, elcode = NA, synonyms = NA, uniqueId = NA
          )
        }
        
        shinyjs::show(id = "data_panel")
        
        shinybusy::remove_modal_spinner()
        
      }
    }
    
  })
  
  observeEvent(input$main_map_draw_new_feature, {
    
      if (input$main_map_draw_new_feature$geometry$type == "Polygon" & !is.null(taxon_data$sf_filtered)){
        
        shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
        
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
        
        shinybusy::remove_modal_spinner()
        
      }
      
      if (input$main_map_draw_new_feature$geometry$type == "Point"){
        
        shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
        
        drawn_point <- data.frame(longitude = input$main_map_draw_new_feature$geometry$coordinates[[1]],
                                  latitude = input$main_map_draw_new_feature$geometry$coordinates[[2]],
                                  year = 2024,
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
        
        shinybusy::remove_modal_spinner()
        
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
      
      m <- leafletProxy("main_map") %>%
        clearShapes() %>% 
        clearMarkers() %>% 
        clearMarkerClusters() 
      
      if (nrow(taxon_data$sf) >= 2){

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
      
      updateDateRangeInput(session = session, inputId = "year_filter",
                           start = paste0(min(unique(taxon_data$sf$year), na.rm = TRUE), "-01-01"),
                           end = paste0(max(unique(taxon_data$sf$year), na.rm = TRUE), "-01-01")
      )
      updateMaterialSwitch(session = session, inputId = "range_extent", value = FALSE)
      updateMaterialSwitch(session = session, inputId = "area_of_occupancy", value = FALSE)
      updateMaterialSwitch(session = session, inputId = "number_EOs", value = FALSE)
      updateSelectizeInput(session = session, inputId = "nation_filter", choices = taxon_data$nations, selected = NULL)
      updateSelectizeInput(session = session, inputId = "states_filter", choices = taxon_data$states, selected = NULL) # , taxon_data$states)
      updateSelectizeInput(session = session, inputId = "synonyms_filter", choices = unique(taxon_data$sf$scientificName), selected = unique(taxon_data$sf$scientificName))
      updateSelectizeInput(session = session, inputId = "sources_filter", choices = unique(taxon_data$sf$prov), selected = unique(taxon_data$sf$prov))
      
      shinybusy::remove_modal_spinner()
      
      observeEvent({
        input$year_filter
        input$no_year
        input$uncertainty_filter
        input$nation_filter
        input$states_filter
        input$synonyms_filter
        input$sources_filter
        input$remove_selections
      }, {
        
        shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
        
        if (!is.null(taxon_data$sf)){
          
          taxon_data$sf_filtered <- taxon_data$sf
          

          taxon_data$sf_filtered <- taxon_data$sf_filtered %>%
            dplyr::filter(year >= substr(input$year_filter[1], 1, 4) & year <= substr(input$year_filter[2], 1, 4) | is.na(year),
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
          
          if (isTRUE(input$remove_selections)){
            
            taxon_data$removed_points <- rbind(taxon_data$removed_points, taxon_data$selected_points)
            
          }
          
          name_exclusions <- setdiff(taxon_data$sf_filtered$scientificName, input$synonyms_filter)
          if (length(name_exclusions) > 0){
            taxon_data$sf_filtered <- taxon_data$sf_filtered %>%
              dplyr::filter(!(scientificName %in% name_exclusions))
          }
          
          source_exclusions <- setdiff(taxon_data$sf_filtered$prov, input$sources_filter)
          
          if (length(source_exclusions) > 0){
            taxon_data$sf_filtered <- taxon_data$sf_filtered %>%
              dplyr::filter(!(prov %in% source_exclusions))
          }
          
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
          
          shinybusy::remove_modal_spinner() # remove the modal window
          
        }
        
      }) 
      
    }
  })
 
  observeEvent(input$main_map_marker_click, {
    
    shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
    
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
    
    shinybusy::remove_modal_spinner() # remove the modal window
    
  })
  
  output$occurrences_barchart_full <- dygraphs::renderDygraph({
    
    dat <- taxon_data$sf_filtered %>%
      dplyr::filter(complete.cases(year)) %>% 
      dplyr::group_by(year) %>%
      dplyr::summarise(number_records = n()) %>%
      dplyr::mutate(year = paste0(year, "-01-01") %>% as.Date())
    
    taxon_data$records_over_time <- xts::xts(x = dat$number_records, order.by = dat$year)
    
    start_window <- "1980-01-01"
    end_window <- Sys.Date()
    
    dygraph(taxon_data$records_over_time, ylab = "") %>%
      dyBarChart() %>%
      dySeries("V1", label = "Number of records", color = "#1f417d") %>%
      dyAxis(
        "y",
        axisLabelWidth = 0
      ) %>% 
      dyRangeSelector(dateWindow = c(start_window %>% as.Date(), end_window %>% as.Date()))
    
  })
  
  
  output$occurrences_barchart_period <- plotly::renderPlotly({
    
    dat <- taxon_data$sf_filtered %>%
      dplyr::filter(complete.cases(year)) %>% 
      dplyr::mutate(period = case_when(
        year >= substr(input$period1[1], 1, 4) & year < substr(input$period1[2], 1, 4) ~ "1",
        year >= substr(input$period2[1], 1, 4) & year < substr(input$period2[2], 1, 4) ~ "2",
        year >= substr(input$period3[1], 1, 4) & year < substr(input$period3[2], 1, 4) ~ "3"
        # year >= substr(input$period4[1], 1, 4) & year < substr(input$period4[2], 1, 4) ~ 4
      )
      ) %>%
      dplyr::filter(complete.cases(period)) %>% 
      dplyr::group_by(period) %>%
      dplyr::summarise(number_records = n()) %>%
      dplyr::mutate(year = paste0(period, "-01-01") %>% as.Date())
    
    p <- ggplot(data = dat) + 
      geom_bar(mapping = aes(x = period, y = number_records), stat = "identity", fill = "#1f417d") +
      theme_bw() +
      theme(panel.background = element_rect(fill = "#F9F9F9")) +
      xlab("Time period") +
      ylab("Number of records")
    
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
  
  output$occurrences_table <- DT::renderDataTable({
    
    dat <- taxon_data$selected_points

    if ("sf" %in% class(dat)){
      dat <- dat %>%
        sf::st_set_geometry(NULL)
    }

    dat %>%
      dplyr::select(all_of(minimum_fields)) %>%
      dplyr::mutate(references = paste0("<a href='", references, "' target='_blank'>", references, "</a>"),
        latitude = round(latitude, 4),
        longitude = round(longitude, 4)
        ) %>%
      dplyr::rename("Key" = key,
                    "Scientific name" = scientificName,
                    "Source" = prov,
                    "Latitude" = latitude,
                    "Longitude" = longitude,
                    "Institution" = institutionCode,
                    "Year" = year,
                    "Coordinate Uncertainty" = coordinateUncertaintyInMeters,
                    "Subnation" = stateProvince,
                    "Nation" = countryCode,
                    "References" = references
      ) %>%
      DT::datatable(options = list(dom = 'tp',
                                   pageLength = 10,
                                   columnDefs = list(list(width = "10%", className = 'dt-left', targets = c(1,2))),
                                   language = list(emptyTable = 'You have not selected any records from the map'),
                                   width = "100%"
      ),
      # filter = list(position = 'top'),
      selection = list(mode = 'multiple', target = 'row', selected = 0:nrow(taxon_data$selected_points)),
      escape = FALSE,
      rownames = FALSE
      )
    
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
    
    shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
    
    taxon_data$selected_points <- taxon_data$selected_points %>% 
      dplyr::filter(key %in% taxon_data$selected_points$key[NULL])  
    
    shinybusy::remove_modal_spinner()
    
  })
  
  observeEvent(input$no_year, { 
      
      if (isTRUE(input$no_year)){
        
        no_year_keys <- taxon_data$sf_filtered %>% dplyr::filter(is.na(year)) %>% dplyr::pull(key)
        
        print(no_year_keys)
        
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

      shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
      
      shinyjs::show(id = "EOO_panel")

      taxon_data$filtered_occurrences <- taxon_data$sf_filtered

      taxon_data$filtered_occurrences <- taxon_data$filtered_occurrences %>%
        st_set_geometry(NULL) %>%
        dplyr::select(longitude, latitude) %>%
        as.data.frame()

      if (nrow(taxon_data$filtered_occurrences) >= 3){

        
        if (taxon_data$shifted){
          eoo_output <- taxon_data$sf_filtered %>% calculate_eoo(shifted = TRUE)
        } else {
          eoo_output <- taxon_data$sf_filtered %>% calculate_eoo(shifted = FALSE)
        }
        
        print("Range calculation went through")
        
        taxon_data$species_range_value <- eoo_output$EOO
        taxon_data$species_range_map <- eoo_output$hull

        m <- leafletProxy("main_map") %>%
          clearShapes() %>%
          leaflet::addMapPane("species_range", zIndex = 200) %>%
          addPolygons(data = taxon_data$species_range_map,
                      color = grey(.2),
                      fillOpacity = 0.1,
                      fill = TRUE,
                      options = pathOptions(pane = "species_range"),
                      group = "Range Extent"
          )

        if (input$area_of_occupancy & input$number_EOs){
          m <- leafletProxy("main_map") %>%
            leaflet::addMapPane("aoo", zIndex = 200) %>%
            addPolygons(data = taxon_data$AOO_map,
                        color = "#2c7bb6",
                        fill = "#2c7bb6",
                        fillOpacity = .2,
                        options = pathOptions(pane = "aoo"),
                        group = "Occupancy"
            ) %>%
            leaflet::addMapPane("eos", zIndex = 200) %>%
            addPolygons(data = taxon_data$EOcount_map,
                        fill = TRUE,
                        fillColor = "#8b0000",
                        fillOpacity = .5,
                        opacity = 0,
                        options = pathOptions(pane = "eos"),
                        group = "Occurrences"
            )
        }

        if (input$area_of_occupancy & isFALSE(input$number_EOs)){
          m <- leafletProxy("main_map") %>%
            leaflet::addMapPane("aoo", zIndex = 200) %>%
            addPolygons(data = taxon_data$AOO_map,
                        color = "#2c7bb6",
                        fill = "#2c7bb6",
                        fillOpacity = .2,
                        options = pathOptions(pane = "aoo"),
                        group = "Occupancy"
            )
        }

        if (input$number_EOs & isFALSE(input$area_of_occupancy)){
          m <- leafletProxy("main_map") %>%
            leaflet::addMapPane("eos", zIndex = 200) %>%
            addPolygons(data = taxon_data$EOcount_map,
                        fill = TRUE,
                        fillColor = "#8b0000",
                        fillOpacity = .5,
                        opacity = 0,
                        options = pathOptions(pane = "eos"),
                        group = "Occurrences"
            )
        }

        m

        shinybusy::remove_modal_spinner()
        
      }

    } else {

      taxon_data$species_range_value <- NA

      shinyjs::hide(id = "EOO_panel")
      leafletProxy("main_map") %>%
        leaflet::clearGroup("Range Extent")
    }

  })
  
  observeEvent({
    input$area_of_occupancy
    input$grid_cell_size
  }, {
    
    if (isTRUE(input$area_of_occupancy)){
      
      shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
      
      shinyjs::show(id = "AOO_panel")
      
      taxon_data$filtered_occurrences <- taxon_data$sf_filtered
      taxon_data$filtered_occurrences <- taxon_data$filtered_occurrences %>%
        st_set_geometry(NULL) %>%
        dplyr::select(longitude, latitude) %>%
        as.data.frame()
      taxon_data$filtered_occurrences$longitude[taxon_data$filtered_occurrences$longitude > 180] <- taxon_data$filtered_occurrences$longitude[taxon_data$filtered_occurrences$longitude > 180] - 360
      
      taxon_data$AOO_value <- (aoo2(taxon_data$filtered_occurrences, as.numeric(input$grid_cell_size)*1000))/4
      
      taxon_data$AOO_map <- get_aoo_polys(taxon_data$sf_filtered, as.numeric(input$grid_cell_size))
      # taxon_data$AOO_value <- taxon_data$AOO_map %>% nrow()
      
      if (!is.null(taxon_data$AOO_map)){
        
      m <- leafletProxy("main_map") %>%
        clearShapes() %>% 
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
        sendSweetAlert(session, type = "warning", title = "Oops!", text = "There are too many AOO cells to be mapped efficiently", closeOnClickOutside = TRUE)
      }
      
      if (input$range_extent & input$number_EOs){
        m <- leafletProxy("main_map") %>%
          leaflet::addMapPane("species_range", zIndex = 200) %>% 
          addPolygons(data = taxon_data$species_range_map,
                      color = grey(.2),
                      fillOpacity = 0.1,
                      fill = TRUE,
                      options = pathOptions(pane = "species_range")
                      # group = "Range Extent"
          ) %>% 
          leaflet::addMapPane("eos", zIndex = 200) %>% 
          addPolygons(data = taxon_data$EOcount_map, 
                      fill = TRUE,
                      fillColor = "#8b0000",
                      fillOpacity = .5,
                      opacity = 0,
                      options = pathOptions(pane = "eos"),
                      group = "Occurrences"
          )
      } 
      
      if (input$range_extent & isFALSE(input$number_EOs)){
        m <- leafletProxy("main_map") %>%
          leaflet::addMapPane("species_range", zIndex = 200) %>% 
          addPolygons(data = taxon_data$species_range_map,
                      color = grey(.2),
                      fillOpacity = 0.1,
                      fill = TRUE,
                      options = pathOptions(pane = "species_range"),
                      group = "Range Extent"
          )
      } 
      
      if (input$number_EOs & isFALSE(input$range_extent)){
        m <- leafletProxy("main_map") %>%
          leaflet::addMapPane("eos", zIndex = 200) %>% 
          addPolygons(data = taxon_data$EOcount_map, 
                      fill = TRUE,
                      fillColor = "#8b0000",
                      fillOpacity = .5,
                      opacity = 0,
                      options = pathOptions(pane = "eos"),
                      group = "Occurrences"
          )
      } 
      
      m
      
      shinybusy::remove_modal_spinner()
      
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
      
      shinybusy::show_modal_spinner("circle", color = "#024b6c") # show the modal window
      
      shinyjs::show(id = "EOcount_panel")
      
      ##### Calculate numbers of EOs
      number_EOs <- calculate_number_occurrences(taxon_data$sf_filtered, separation_distance = input$separation_distance %>% as.numeric(), added_distance = 0)
      taxon_data$EOcount_value <- number_EOs$eo_count
      taxon_data$EOcount_map <- number_EOs$buffered_occurrences
      
      m <- leafletProxy("main_map") %>%
        clearShapes() %>%
        leaflet::addMapPane("eos", zIndex = 200) %>% 
        addPolygons(data = taxon_data$EOcount_map, 
                    fill = TRUE,
                    fillColor = "#8b0000",
                    fillOpacity = .5,
                    opacity = 0,
                    options = pathOptions(pane = "eos"),
                    group = "Occurrences"
        )
      
      if (input$area_of_occupancy & input$range_extent){
        m <- leafletProxy("main_map") %>%
          leaflet::addMapPane("aoo", zIndex = 200) %>% 
          addPolygons(data = taxon_data$AOO_map,
                      color = "#2c7bb6",
                      fill = "#2c7bb6",
                      fillOpacity = .2,
                      options = pathOptions(pane = "aoo"),
                      group = "Occupancy"
          ) %>% 
          leaflet::addMapPane("species_range", zIndex = 200) %>% 
          addPolygons(data = taxon_data$species_range_map,
                      color = grey(.2),
                      fillOpacity = 0.1,
                      fill = TRUE,
                      options = pathOptions(pane = "species_range"),
                      group = "Range Extent"
          )
      } 
      
      if (input$area_of_occupancy & isFALSE(input$range_extent)){
        m <- leafletProxy("main_map") %>%
          leaflet::addMapPane("aoo", zIndex = 200) %>% 
          addPolygons(data = taxon_data$AOO_map,
                      color = "#2c7bb6",
                      fill = "#2c7bb6",
                      fillOpacity = .2,
                      options = pathOptions(pane = "aoo"),
                      group = "Occupancy"
          )
      } 
      
      if (input$range_extent & isFALSE(input$area_of_occupancy)){
        m <- leafletProxy("main_map") %>%
          leaflet::addMapPane("species_range", zIndex = 200) %>% 
          addPolygons(data = taxon_data$species_range_map,
                      color = grey(.2),
                      fillOpacity = 0.1,
                      fill = TRUE,
                      options = pathOptions(pane = "species_range"),
                      group = "Range Extent"
          )
      } 
      
      m
      
      shinybusy::remove_modal_spinner()
      
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
    
    if (!is.na(taxon_data$species_range_value)){
      
      HTML(
        paste0(
          h2(format(as.numeric(req(taxon_data$species_range_value)), big.mark=",", scientific = FALSE), style = "padding-top: 0; margin-top: 0; display: inline;"), 
          h2(" km", tags$sup("2"), style = "padding-top: 0; margin-top: 0; display: inline;"),
          h1(strong(base::cut(req(as.numeric(taxon_data$species_range_value)), breaks = c(0, 0.999, 99.999, 249.999, 999.999, 4999.999, 19999.999, 199999.999, 2499999.999, 1000000000), labels = c("Z", LETTERS[1:8]))), style = "padding-top: 0; margin-top: 0; padding-left: 1em; display: inline;")
        )
      )
    } else {
      HTML(paste0(h3("NA")))
    }
    
  })
  
  output$AOO_value <- renderUI({
    
    if (input$grid_cell_size == 1) out <- base::cut(req(as.numeric(taxon_data$AOO_value)), breaks = c(0, 0.999, 4.999, 10.999, 20.999, 100.999, 500.999, 2000.999, 10000.999, 50000.999, 1000000000), labels = c("Z", LETTERS[1:9]))
    
    if (input$grid_cell_size > 1) out <- base::cut(req(as.numeric(taxon_data$AOO_value)), breaks = c(0, 0.999, 1.999, 2.999, 5.999, 25.999, 125.999, 500.999, 2500.999, 12500.999, 1000000000), labels = c("Z", LETTERS[1:9]))
    
    HTML(
      paste0(
        h2(format(as.numeric(req(taxon_data$AOO_value)), big.mark=",", scientific = FALSE), style = "padding-top: 0; margin-top: 0; display: inline;"), 
        h2(paste0(" cells"), style = "padding-top: 0; margin-top: 0; display: inline;"),
        h1(strong(out), style = "padding-top: 0; margin-top: 0; padding-left: 1em; display: inline;")
      )
    )
    
  })
  
  output$EOcount_value <- renderUI({
    
    HTML(
      paste0(
        h2(format(as.numeric(req(taxon_data$EOcount_value)), big.mark=",", scientific = FALSE), style = "padding-top: 0; margin-top: 0; display: inline;"), 
        h2("  EOs", style = "padding-top: 0; margin-top: 0; display: inline;"),
        h1(strong(cut(req(as.numeric(taxon_data$EOcount_value)), breaks = c(0, 0.999, 5.999, 19.999, 79.999, 299.999, 1000000000), labels = c("Z", LETTERS[1:5]))), style = "padding-top: 0; margin-top: 0; padding-left: 1em; display: inline;")
      )
    )
    
  })
  
  output$number_occurrences <- renderUI({
    
    if (!is.null(taxon_data$sf_filtered)){
      sources_count <- taxon_data$sf_filtered %>% dplyr::count(prov)
      sources_filter_labels <- paste0(paste0(sources_count$prov, ": ", sources_count$n), collapse = "; ")
      records_count_text <- paste0("  records added (", sources_filter_labels, ")")
      
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
      paste(gsub(" ", "_", ifelse(!is.null(taxon_data$info), taxon_data$info$scientificName, "New taxon")), "-rank_factor_values-", Sys.Date(), ".csv", sep="")
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
      out[, 2] <-  ifelse(!is.null(taxon_data$info), taxon_data$info$roundedGRank, "")
      out[, 3] <-  ifelse(!is.null(taxon_data$info), taxon_data$info$scientificName, "")
      out[, 6] <-  ifelse(!is.null(taxon_data$info), taxon_data$info$elementGlobalId, "")
      out[, 7] <-  ifelse(!is.null(taxon_data$info), taxon_data$info$elcode, "")
      out[, 8] <-  ifelse(!is.null(taxon_data$info), taxon_data$info$primaryCommonName, "")
      out[, 9] <-  ifelse(!is.null(taxon_data$info_extended), taxon_data$info_extended$classificationStatus$classificationStatusDescEn, "")
      if (!is.null(taxon_data$species_range_value)){
        out[, 11] <- cut(as.numeric(taxon_data$species_range_value), breaks = c(0, 0.999, 99.999, 249.999, 999.999, 4999.999, 19999.999, 199999.999, 2499999.999, 1000000000), labels = c("Z", LETTERS[1:8]))
      }
      if (!is.null(taxon_data$AOO_value)){
        if (input$grid_cell_size == 2){
          out[, 13] <- base::cut(as.numeric(taxon_data$AOO_value), breaks = c(0, 0.999, 1.999, 2.999, 5.999, 24.999, 124.999, 499.999, 2499.999, 12499.999, 1000000000), labels = c("Z", LETTERS[1:9]))
        } else if (input$grid_cell_size == 1){
          out[, 13] <- base::cut(as.numeric(taxon_data$AOO_value), breaks = c(0, 0.999, 4.999, 10.999, 20.999, 100.999, 500.999, 2000.999, 10000.999, 50000.999, 1000000000), labels = c("Z", LETTERS[1:9]))
        } 
      }
      if (!is.null(taxon_data$EOcount_value)){
      out[, 15] <- cut(as.numeric(taxon_data$EOcount_value), breaks = c(0, 0.999, 5.999, 19.999, 79.999, 299.999, 1000000000), labels = c("Z", LETTERS[1:5]))
      }
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

      write.csv(out, file, row.names = FALSE, na = "")
      
    }
  )
  
  output$download_occurrence_data <- downloadHandler(
    filename = function() {
      paste(gsub(" ", "_", ifelse(!is.null(taxon_data$info), taxon_data$info$scientificName, "New taxon")), "-RARECAT_assessment_records-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      out <- taxon_data$sf_filtered %>% 
        st_set_geometry(NULL)
      write.csv(out, file, row.names = FALSE)
    }
  )
  
  observeEvent(input$clear_map, {
    
    session$reload()
    
  })
  
  batch_run_taxon_list <- reactiveValues(names = NULL)
  
  batch_run_output <- reactiveValues(
    results = NULL,
    table = NULL
  )
  
  observeEvent(input$batch_assessment, {
    
    batch_run_taxon_list$names <- strsplit(input$typed_list, "\n|,|;|, |; ")[[1]] %>% as.data.frame() %>% set_names("user_supplied_name")
    
    batch_run_output$results <- vector("list", length(batch_run_taxon_list$names$user_supplied_name))
    
    print(batch_run_taxon_list$names)
    
    # Create a Progress object
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    
    progress$set(message = "Running Batch Assessment", value = 0)
    
    n <- length(batch_run_taxon_list$names$user_supplied_name)
    
    for (i in 1:n) {
      
      taxon_name <- batch_run_taxon_list$names$user_supplied_name[i]
      
      batch_run_output$results[[i]] <- run_rank_assessment(taxon_name = taxon_name)
      
      # Increment the progress bar, and update the detail text.
      progress$inc(1/n, detail = paste("Taxon", i))
      
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
    }
    
    batch_run_output$results <- batch_run_output$results %>% 
      set_names(batch_run_taxon_list$names$user_supplied_name)
    
    batch_run_output$table <- data.frame(
      taxon = batch_run_taxon_list$names$user_supplied_name,
      total_observations_used = purrr::map(batch_run_output$results, function(out) nrow(out$sf_filtered)) %>% unlist(),
      range_value = purrr::map(batch_run_output$results, function(out) out$species_range_value) %>% unlist(),
      AOO_value = purrr::map(batch_run_output$results, function(out) out$AOO_value) %>% unlist(),
      EOcount_value = purrr::map(batch_run_output$results, function(out) out$EOcount_value) %>% unlist()
    )
    
  })
  
  output$batch_run_results_table <- DT::renderDataTable({
    
    batch_run_output$table %>%
      # dplyr::rename("Scientific name" = taxon,
      #               "Total Observations Used" = total_observations_used,
      #               "Range Extent" = range_value,
      #               "AOO" = AOO_value,
      #               "Occurrence Count" = EOcount_value
      # ) %>%
      DT::datatable(options = list(dom = 'tp',
                                   pageLength = 10,
                                   columnDefs = list(list(width = "10%", className = 'dt-left', targets = c(1,2))),
                                   width = "100%"
      ),
      # filter = list(position = 'top'),
      selection = list(mode = 'multiple', target = 'row'),
      escape = FALSE,
      rownames = FALSE
      )
    
  })
  
}
