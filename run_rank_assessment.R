run_rank_assessment <- function(taxon_name, 
                                minimum_fields =  c("key", "scientificName", "prov", "longitude", "latitude", "coordinateUncertaintyInMeters", "stateProvince", "countryCode", "year", "institutionCode", "references"),
                                max_number_observations = 10000,
                                clean_occ = TRUE,
                                centroid_filter = FALSE,
                                date_start = "1900-01-01",
                                date_end = "2025-01-01",
                                uncertainty_filter = "",
                                nations_filter = NULL,
                                states_filter = NULL,
                                sources_filter = NULL,
                                grid_cell_size = 2,
                                sep_distance = 1000
                                
){
  
  taxon_data <- list(
    info = data.frame(scientificName = "New taxon"),
    family = NULL,
    synonyms = NULL,
    synonyms_selected = NULL,
    gbif_occurrences_raw = NULL,
    gbif_occurrences = NULL,
    uploaded_occurrences = NULL,
    all_occurrences = NULL,
    shifted = FALSE,
    sf = NULL,
    sf_filtered = NULL,
    species_range_value = NULL,
    species_range_map = NULL,
    AOO_value = NULL,
    AOO_map = NULL,
    EOcount = NULL,
    EOcount_value = NULL
  )
  
  ns_table_full <- natserv::ns_search_spp(text_adv = list(searchToken = taxon_name, matchAgainst = "allScientificNames", operator="contains"))$results
  gbif_table <- rgbif::name_suggest(q = taxon_name, rank = c("species", "subspecies"), limit = 10)$data
  
  taxon_info <- NULL
  if (nrow(ns_table_full) > 0){
    ns_table <- ns_table_full %>% 
      dplyr::mutate(Source = "NatureServe", synonyms = ns_table_full$speciesGlobal$synonyms) %>% 
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
  taxon_data$info <- taxon_info
  taxon_data$family <- ns_table_full$speciesGlobal$family
  selected_taxon <- c(taxon_info$scientificName, unlist(taxon_info$synonyms)) %>% na.omit() %>% unique()
  selected_taxon <- gsub("ssp. |var. ", "", selected_taxon)
  
  if (length(selected_taxon) == 1){
    taxon_data$synonyms <- rgbif::name_usage(name = selected_taxon)$data
  } else {
    taxon_data$synonyms <- purrr::map(selected_taxon, function(sp) rgbif::name_usage(name = sp)) %>%
      purrr::map("data") %>%
      bind_rows()
  }
  
  taxon_data$synonyms <- taxon_data$synonyms %>% 
    dplyr::distinct(., .keep_all = TRUE) 
  
  taxon_data$gbif_occurrences_raw <- get_gbif_data(sp_data = taxon_data$synonyms, number_observations = max_number_observations)$sp_occurrences
  
  taxon_data$gbif_occurrences <- taxon_data$gbif_occurrences_raw %>% 
    clean_gbif_data(clean = clean_occ, remove_centroids = centroid_filter, minimum_fields = minimum_fields)
  
  # out <- purrr::map(input$filedata$datapath, process_user_data, minimum_fields = minimum_fields) %>% 
  #   dplyr::bind_rows() 
  # 
  # taxon_data$uploaded_occurrences <- uploaded_data()
  # 
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
    dplyr::filter(year >= substr(date_start, 1, 4) & year <= substr(date_end, 1, 4) | is.na(year))
  
  if (uncertainty_filter != ""){
    
    taxon_data$sf_filtered <- taxon_data$sf_filtered %>%
      dplyr::filter(coordinateUncertaintyInMeters <= as.numeric(uncertainty_filter) | is.na(coordinateUncertaintyInMeters))
  }
  
  if (!is.null(nations_filter)){
    taxon_data$sf_filtered <- taxon_data$sf_filtered %>%
      dplyr::filter(purrr::map_int(st_intersects(taxon_data$sf_filtered, network_polys %>% dplyr::filter(FIPS_CNTRY %in% nations_filter)), length) > 0)
  }
  
  if (!is.null(states_filter)){
    taxon_data$sf_filtered <- taxon_data$sf_filtered %>%
      dplyr::filter(purrr::map_int(st_intersects(taxon_data$sf_filtered, network_polys %>% dplyr::filter(Admin_abbr %in% states_filter)), length) > 0)
  }
  
  if (!is.null(sources_filter)){
    source_exclusions <- setdiff(taxon_data$sf_filtered$prov, sources_filter)
    
    if (length(source_exclusions) > 0){
      taxon_data$sf_filtered <- taxon_data$sf_filtered %>%
        dplyr::filter(!(prov %in% source_exclusions))
    }
  }

  eoo_output <- taxon_data$sf_filtered %>% calculate_eoo(shifted = TRUE)
  taxon_data$species_range_value <- eoo_output$EOO
  taxon_data$species_range_map <- eoo_output$hull
  taxon_data$AOO_value <- (aoo2(taxon_data$sf_filtered, as.numeric(grid_cell_size)*1000))/4
  taxon_data$AOO_map <- get_aoo_polys(taxon_data$sf_filtered, as.numeric(grid_cell_size))
  number_EOs <- calculate_number_occurrences(taxon_data$sf_filtered, separation_distance = sep_distance %>% as.numeric(), added_distance = 0)
  taxon_data$EOcount_value <- number_EOs$eo_count
  taxon_data$EOcount_map <- number_EOs$buffered_occurrences
  
  return(taxon_data)
  
}

taxon_names <- c("Acacia millefolia", "Juncus abortivus", "Artemisia pattersonii")
Acacia millefolia
Juncus abortivus
Artemisia pattersonii
x <- purrr::map(taxon_names, run_rank_assessment) %>% 
  set_names(taxon_names)
