run_rank_assessment <- function(taxon_name){
  
  taxon_data <- list(
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
  
  ns_table <- natserv::ns_search_spp(text_adv = list(searchToken = taxon_name, matchAgainst = "allScientificNames", operator="contains"))$results
  gbif_table <- rgbif::name_suggest(q = taxon_name, rank = c("species", "subspecies"), limit = 10)$data
  
  taxon_info <- NULL
  if (nrow(ns_table) > 0){
    ns_table <- ns_table %>% 
      dplyr::mutate(Source = "NatureServe", synonyms = ns_table$speciesGlobal$synonyms) %>% 
      dplyr::select(scientificName, Source, elementGlobalId, primaryCommonName, roundedGRank, elcode, uniqueId, synonyms)
    taxon_info <- rbind(taxon_info, ns_table)
  } 
  if (nrow(gbif_table) > 0){
    gbif_table <- gbif_table %>% 
      dplyr::rename(scientificName = canonicalName, elementGlobalId = key) %>% 
      dplyr::mutate(Source = "GBIF", primaryCommonName = NA, roundedGRank = NA, elcode = NA, synonyms = NA, uniqueId = NA) %>% 
      dplyr::select(scientificName, Source, elementGlobalId, primaryCommonName, roundedGRank, elcode, uniqueId, synonyms)
    taxon_info <- rbind(taxon_info, gbif_table)
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
  
  gbif_download <- get_gbif_data(sp_data = taxon_data$synonyms_selected, number_observations = input$number_gbif_occurrences)
  
  taxon_data$gbif_occurrences <- taxon_data$gbif_occurrences_raw %>% 
    clean_gbif_data(clean = input$clean_occ, remove_centroids = input$centroid_filter, minimum_fields = minimum_fields)
  
  out <- purrr::map(input$filedata$datapath, process_user_data, minimum_fields = minimum_fields) %>% 
    dplyr::bind_rows() 
  
  taxon_data$uploaded_occurrences <- uploaded_data()
  
  taxon_data$all_occurrences <- rbind(
    taxon_data$gbif_occurrences,
    taxon_data$uploaded_occurrences
  )
  
  ### Create simple features object for geospatial calculations
  taxon_data$sf <- taxon_data$all_occurrences %>% 
    dplyr::filter(complete.cases(longitude, latitude)) %>% 
    dplyr::mutate(lon = longitude,
                  lat = latitude) %>% 
    sf::st_as_sf(
      coords = c("lon", "lat"),
      crs = 4326
    )
  
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
  
  eoo_output <- taxon_data$sf_filtered %>% calculate_eoo(shifted = TRUE)
  
  taxon_data$AOO_value <- (aoo2(taxon_data$filtered_occurrences, as.numeric(input$grid_cell_size)*1000))/4
  
  number_EOs <- calculate_number_occurrences(taxon_data$sf_filtered, separation_distance = input$separation_distance %>% as.numeric(), added_distance = 0)
  
  
}
