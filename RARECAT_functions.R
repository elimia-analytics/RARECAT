# Function to extract observations from GBIF for target taxon
get_gbif_data <- function(sp_data, number_observations, gbif = TRUE, inat = TRUE, ebird = TRUE){
  
  gbif_occurrences <- NULL
  
  if (gbif){
  # Get all non-human observations
  ## Extract species occurrences across all relevant scientific names from GBIF using the SPOCC package
  gbif_occurrences_nh <- purrr::map(sp_data$key, function(sp){
    gbif_data <- spocc::occ(from = "gbif", 
               gbifopts = list(
                 taxonKey = sp, 
                 basisOfRecord = c("OCCURRENCE", "PRESERVED_SPECIMEN", "OBSERVATION", "MACHINE_OBSERVATION")
               ),
               limit = 10000, # as.integer(number_observations), 
               has_coords = TRUE
    )
    gbif_data$gbif$data %>% dplyr::bind_rows() 
  }) %>% bind_rows()

  if (nrow(gbif_occurrences_nh) > 0){
  ### Bind results across all relevant scientific names
    gbif_occurrences_nh <- gbif_occurrences_nh %>%
    dplyr::filter(complete.cases(latitude, longitude, basisOfRecord)) %>% # Only keep records that have at a minimum a value for longitude, latitude, and basisOfRecord
    dplyr::filter(latitude != 0 | longitude != 0) # Exclude records that have latitude and longitude values of 0

  if (!is.null(gbif_occurrences)){
    shared_names <- intersect(names(gbif_occurrences), names(gbif_occurrences_nh))
    gbif_occurrences <- gbif_occurrences %>% 
      dplyr::select(all_of(shared_names)) %>% 
      rbind(gbif_occurrences_nh %>% dplyr::select(all_of(shared_names)))
  } else {
    gbif_occurrences <- gbif_occurrences %>% 
      rbind(gbif_occurrences_nh)
  }
  
  }
  }
  
  if (inat){
    ## Extract species occurrences across all relevant scientific names from GBIF using the SPOCC package
    inat_data <- spocc::occ(from = "gbif", 
                             gbifopts = list(
                               taxonKey = sp_data$key, 
                               institutionCode = "iNaturalist"
                             ),
                             limit = 10000, # as.integer(number_observations), 
                             has_coords = TRUE
    )
    
    inat_occurrences <- inat_data$gbif$data %>%
      dplyr::bind_rows()
    
    if (nrow(inat_occurrences) > 0){
    ### Bind results across all relevant scientific names
    inat_occurrences <- inat_occurrences %>% 
      dplyr::filter(complete.cases(latitude, longitude, basisOfRecord)) %>% # Only keep records that have at a minimum a value for longitude, latitude, and basisOfRecord
      dplyr::filter(latitude != 0 | longitude != 0) %>% # Exclude records that have latitude and longitude values of 0
      dplyr::mutate(prov = "inat")
    
    if (!is.null(gbif_occurrences)){
      shared_names <- intersect(names(gbif_occurrences), names(inat_occurrences))
      print(shared_names)
      gbif_occurrences <- gbif_occurrences %>% 
        dplyr::select(all_of(shared_names)) %>% 
        rbind(inat_occurrences %>% dplyr::select(all_of(shared_names)))
    } else {
      gbif_occurrences <- gbif_occurrences %>% 
        rbind(inat_occurrences)
    }
    }
  }
  
  if (ebird){
    ## Extract species occurrences across all relevant scientific names from GBIF using the SPOCC package
    ebird_data <- spocc::occ(from = "gbif", 
                            gbifopts = list(
                              taxonKey = sp_data$key, 
                              datasetKey = "4fa7b334-ce0d-4e88-aaae-2e0c138d049e"
                            ),
                            limit = 10000, # as.integer(number_observations), 
                            has_coords = TRUE
    )
    
    ebird_occurrences <- ebird_data$gbif$data %>%
      dplyr::bind_rows()
    
    if (nrow(ebird_occurrences) > 0){
    ### Bind results across all relevant scientific names
      ebird_occurrences <- ebird_occurrences %>% 
      dplyr::filter(complete.cases(latitude, longitude, basisOfRecord)) %>% # Only keep records that have at a minimum a value for longitude, latitude, and basisOfRecord
      dplyr::filter(latitude != 0 | longitude != 0) %>% # Exclude records that have latitude and longitude values of 0
      dplyr::mutate(prov = "ebird")

      if (!is.null(gbif_occurrences)){
        shared_names <- intersect(names(gbif_occurrences), names(ebird_occurrences))
        gbif_occurrences <- gbif_occurrences %>% 
          dplyr::select(all_of(shared_names)) %>% 
          rbind(ebird_occurrences %>% dplyr::select(all_of(shared_names)))
      } else {
        gbif_occurrences <- gbif_occurrences %>% 
          rbind(ebird_occurrences)
      }
      
    }

  }
  
  # Clean up references field
  gbif_occurrences <- gbif_occurrences %>% 
    dplyr::mutate(references = paste0("https://www.gbif.org/occurrence/", key))

  sp_occurrences <- gbif_occurrences
  
  shifted <- FALSE
  
  sp_occurrences$longitude[sp_occurrences$longitude > 180] <- sp_occurrences$longitude[sp_occurrences$longitude > 180] - 360
  max_long <- max(sp_occurrences$longitude, na.rm = TRUE)/2
  shifted_long <- sp_occurrences$longitude
  
  if (length(sp_occurrences$longitude[sp_occurrences$longitude > max_long]) > 0){
    shifted_long[shifted_long > max_long] <- shifted_long[shifted_long > max_long] - 360
    shifted_long <- shifted_long + 360
    # shifted <- TRUE
  }
  
  if ((max(shifted_long)-min(shifted_long)) < (max(sp_occurrences$longitude) - min(sp_occurrences$longitude))){
    sp_occurrences$longitude <- shifted_long
    shifted <- TRUE
  }

  out <- list(sp_occurrences = sp_occurrences, shifted = shifted)
  
  return(out)
  
}

# Function to clean up GBIF records
clean_gbif_data <- function(gbif_occurrences, clean = TRUE, minimum_fields = minimum_fields, remove_centroids = FALSE){
  
  if (isTRUE(clean)){
    
    if (nrow(gbif_occurrences) > 0){
      
      gbif_occurrences <- gbif_occurrences %>% 
        dplyr::filter(!(basisOfRecord %in% c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN", "MATERIAL_SAMPLE"))) # Exclude these record types
      if ("occurrenceStatus" %in% names(gbif_occurrences)) gbif_occurrences <- gbif_occurrences %>% dplyr::filter(occurrenceStatus == "PRESENT")
      if ("coordinateUncertaintyInMeters" %in% names(gbif_occurrences)) gbif_occurrences <- gbif_occurrences %>% dplyr::filter(!coordinateUncertaintyInMeters %in% c(999, 9999))
      if ("samplingProtocol" %in% names(gbif_occurrences)) gbif_occurrences <- gbif_occurrences %>% dplyr::filter(!samplingProtocol %in% c("from a cultivated plant of known (indirect) wild origin", "grown"))

      if (nrow(gbif_occurrences) > 0){
        
        gbif_occurrences <- gbif_occurrences %>%
          distinct(longitude, latitude, speciesKey, datasetKey, .keep_all = TRUE)
        
      }
      
    } 
    
  }
  
  if (isTRUE(remove_centroids)){

    if ("coordinateUncertaintyInMeters" %in% names(gbif_occurrences)) gbif_occurrences <- gbif_occurrences %>% dplyr::filter(!coordinateUncertaintyInMeters %in% c(301, 3036))
    if ("georeferenceRemarks" %in% names(gbif_occurrences)) gbif_occurrences <- gbif_occurrences %>% dplyr::filter(!grepl("centroid|Centroid|CENTROID", georeferenceRemarks))

  }

  gbif_occurrences <- gbif_occurrences %>%
    dplyr::select(intersect(names(gbif_occurrences), minimum_fields)) %>%
    cbind(matrix(NA, nrow = nrow(gbif_occurrences), ncol = length(setdiff(minimum_fields, names(gbif_occurrences)))) %>%
            as.data.frame() %>%
            set_names(setdiff(minimum_fields, names(gbif_occurrences)))
    ) %>%
    dplyr::select(minimum_fields)

  return(gbif_occurrences)
}

# Function to load and process user data
process_user_data <- function(user_file, minimum_fields){
  
  user_data <- read.csv(user_file, header = TRUE) 
  user_data <- user_data %>% select_if(~!(all(is.na(.)) | all(. == "")))
  processed_data <- NULL
  longitude_names <- c("longitude", "LONGITUDE", "decimalLongitude", "private_longitude", "Decimal Longitude", "Longitude", "lon", "Lon", "X", "x")
  longitude_names <- longitude_names %>% set_names(rep("longitude", length(longitude_names)))
  latitude_names <- c("latitude", "LATITUDE", "decimalLatitude", "private_latitude", "Decimal Latitude", "Latitude", "lat", "Lat", "Y", "y")
  latitude_names <- latitude_names %>% set_names(rep("latitude", length(latitude_names)))
  
  
  if (sum(grepl(paste0(longitude_names, collapse = "|"), names(user_data))) > 0 & sum(grepl(paste0(latitude_names, collapse = "|"), names(user_data))) > 0){
    
    scientificName_names <- c("scientificName", "Scientific name", "GNAME", "scientific_name", "SciName", "scientific name", "Scientific Name")
    scientificName_names <- scientificName_names %>% set_names(rep("scientificName", length(scientificName_names)))
    stateProvince_names <- c("stateProvince", "STATE", "GNAME", "place_state_name", "State", "State or Province", "state", "Province", "province", "PROVINCE", "place_guess")
    stateProvince_names <- stateProvince_names %>% set_names(rep("stateProvince", length(stateProvince_names)))
    countryCode_names <- c("countryCode", "NATN", "country", "place_country_name", "Country", "COUNTRY", "Nation", "NATION")
    countryCode_names <- countryCode_names %>% set_names(rep("countryCode", length(countryCode_names)))
    year_names <- c("year", "Year", "YEAR", "Year Collected", "SiteDate")
    year_names <- year_names %>% set_names(rep("year", length(year_names)))
    coordinateUncertaintyInMeters_names <- c("coordinateUncertaintyInMeters", "public_positional_accuracy")
    lookup <- c(longitude_names, latitude_names, scientificName_names, stateProvince_names, countryCode_names, year_names)
    
    if ("latitude" %in% names(user_data)){
      user_data <- user_data %>% 
        dplyr::select(-any_of(setdiff(latitude_names, "latitude")))
    }
    
    if ("longitude" %in% names(user_data)){
      user_data <- user_data %>% 
        dplyr::select(-any_of(setdiff(longitude_names, "longitude")))
    }
    
    if ("observed_on" %in% names(user_data)){
      user_data <- user_data %>% 
        dplyr::mutate(year = substr(observed_on, 1, 4) %>% as.numeric())
    }

    target_columns <- user_data %>% 
      dplyr::select(matches(paste0(longitude_names, collapse = "|")), 
                    matches(paste0(latitude_names, collapse = "|")), 
                    matches(paste0(scientificName_names, collapse = "|")),
                    matches(paste0(stateProvince_names, collapse = "|")),
                    matches(paste0(countryCode_names, collapse = "|")),
                    matches(paste0(year_names, collapse = "|"))
      ) %>% 
      dplyr::rename(any_of(lookup))
    
      
    processed_data <- target_columns # user_data %>% cbind(target_columns) 
    
    if ("references" %in% names(user_data)){
      processed_data <- processed_data %>% 
        dplyr::mutate(references = user_data$references)
    } else if ("url" %in% names(user_data)){
      processed_data <- processed_data %>% 
        dplyr::mutate(references = user_data$url)      
    } else {
      processed_data <- processed_data %>% 
        dplyr::mutate(references = NA)  
    }
    
    processed_data <- processed_data %>% 
      dplyr::select(intersect(names(processed_data), minimum_fields)) %>%
      cbind(matrix(NA, nrow = nrow(processed_data), ncol = length(setdiff(minimum_fields, names(processed_data)))) %>%
              as.data.frame() %>%
              set_names(setdiff(minimum_fields, names(processed_data)))
      ) %>%
      dplyr::select(minimum_fields) %>%
      dplyr::distinct(., .keep_all = TRUE) %>%    
      dplyr::filter(complete.cases(latitude, longitude), latitude != 0, longitude != 0) 
    
    processed_data <- processed_data %>% 
      dplyr::mutate(prov = "uploaded",
                    key = paste(prov, 1:nrow(processed_data), sep = "_"),
                    year = ifelse(!is.na(year), year, NA), #format(Sys.time(), "%Y")),
                    scientificName = ifelse(!is.na(scientificName), scientificName, "user-uploaded")
      )
    
    processed_data$longitude[processed_data$longitude > 180] <- processed_data$longitude[processed_data$longitude > 180] - 360
    
  }
  return(processed_data)
}

longlat2utm <- function(longlat){
  longlat_df <- longlat %>% 
    sf::st_set_geometry(NULL) %>% 
    dplyr::select(longitude, latitude) %>% 
    as.data.frame()
  longlat_df <- as.matrix(longlat_df)
  minlong = min(longlat_df[,1])
  zone = floor((minlong + 180) / 6) + 1
  zone <- sf::st_crs(paste("+proj=utm +zone=", zone," ellps=WGS84",sep=''))
  res <- sf::st_transform(longlat, zone) %>% 
    sf::st_coordinates()
  return(res)
}

aoo2 <- function (spData, cell_size) 
{
  
  spData_sf <- spData %>% 
    dplyr::filter(complete.cases(longitude, latitude)) %>% 
    dplyr::mutate(lon = longitude,
                  lat = latitude) %>% 
    sf::st_as_sf(
      coords = c("lon", "lat"),
      crs = 4326
    )
  spData <- spData_sf %>%
    sf::st_set_geometry(NULL) %>%
    dplyr::select(longitude, latitude) %>%
    as.data.frame()
  
  if (ncol(spData) == 2) {
    if (max(spData) <= 180) {
      spData <- longlat2utm(spData_sf)
      spData <- floor(spData/cell_size)
      ncells <- nrow(unique(spData))
      aoo <- ncells * 4
    }
    else {
      spData$longitude[spData$longitude > 180] <- spData$longitude[spData$longitude > 180] - 360
      spData <- floor(spData/cell_size)
      ncells <- nrow(unique(spData))
      aoo <- ncells * 4
    }
  }
  return(round(aoo))
}

# Function to calculate extent of occurrence (i.e. range)
calculate_eoo <- function(occurrences_sf, shifted = FALSE){
  
  # Rescale hull for calculation
  if (isTRUE(shifted)){
    scaling_factor <- max(occurrences_sf$longitude)-min(occurrences_sf$longitude)
    # occurrences_sf$geometry <- (sf::st_geometry(occurrences_sf) + c(360,90)) %% c(360) - c(0,90)
    occurrences_sf$geometry <- (sf::st_geometry(occurrences_sf)-c(scaling_factor, 0))
    st_crs(occurrences_sf) <- 4326
    
  }
  centreofpoints <- trueCOGll(st_coordinates(occurrences_sf) %>% as.data.frame() %>% set_names("longitude", "latitude"))
  mypointsxy <- simProjWiz(st_coordinates(occurrences_sf) %>% as.data.frame() %>% set_names(c("long", "lat")), centreofpoints)
  # vertices <- chull(mypointsxy$xy)
  # mypointsxy_sf <- mypointsxy$xy %>% 
  #   sf::st_as_sf(
  #     coords = c("x", "y"),
  #     crs = mypointsxy$crs_new
  #   )
  hull <- occurrences_sf %>%  #[vertices, ] %>% 
    terra::vect() %>% 
    terra::convHull() %>% 
    sf::st_as_sf()
  st_crs(hull) <- 4326
  
  # safe_area <- purrr::safely(sf::st_area)
  EOO <- red::eoo(mypointsxy$xy)
  # if (!is.null(EOO$result)){
  #   EOO <- EOO$result %>% units::set_units(km^2)
  # } else {
  #   EOO <- NA
  # }
  
  # Rescale hull for mapping
  if (isTRUE(shifted)){
    # occurrences_sf$geometry <- (sf::st_geometry(occurrences_sf) + c(360,90)) %% c(360) - c(0,90)
    hull$geometry <- (sf::st_geometry(hull)+c(scaling_factor, 0))
    st_crs(hull) <- 4326
  }
  return(list(hull = hull, EOO = EOO))
}

get_aoo_polys <- function(occurrences_sf, cell_size){
  ## Extract bounding box containing all observations
  sp_bounding_box <- sf::st_bbox(occurrences_sf) 
  ## Create raster with extent equal to focal taxon bounding box and desired resolution (i.e. cell size divided by 111 - the distance equivalent to 1 degree of longitude at the equator)
  aoo_raster <- terra::rast(xmin = sp_bounding_box[[1]], ymin = sp_bounding_box[[2]], xmax = sp_bounding_box[[3]], ymax = sp_bounding_box[[4]], resolution = cell_size/111) # crs = sf::st_as_text(zone), ncol = box_cols, nrow = box_rows) # 
  ## Extend raster by 2 cells in each direction
  aoo_raster <- terra::extend(aoo_raster, 2)
  
  if (terra::ncell(aoo_raster) < 100000000){
  ## Identify AOO grid cells that are overlapped by points
  aoo_cells <- terra::cellFromXY(object = aoo_raster, xy = sf::st_coordinates(occurrences_sf)) %>% unique()
  ## Set values of overlapping AOO grid cells to 1
  aoo_raster[aoo_cells] <- 1
  ## Transform overlapped AOO grid cells to polygons
  aoo_polys <- aoo_raster %>% terra::as.polygons(aggregate = FALSE, na.rm = TRUE) %>% sf::st_as_sf()
  } else {
    aoo_polys <- NULL
  }
  
  ## Return sf object representing overlapped AOO grid cells
  return(aoo_polys)  
}

# Function to calculate number of occurrences
calculate_number_occurrences <- function(occ, separation_distance = 1000, added_distance = 0){
  ##### Buffer occurrences by desired separation distance
  buffered_occurrences <- sf::st_buffer(occ, (separation_distance+added_distance)/2)

  ##### Calculate intersections among occurrences
  buffered_occurrences_intersection <- sf::st_intersects(buffered_occurrences)
    ##### Calculate numbers of EOs by identifying individual clusters of points
    groups <- purrr::map(1:length(buffered_occurrences_intersection), function(i){
      focal_intersection <- buffered_occurrences_intersection[[i]] # for each point, assess its intersections
      intersections <- purrr::map(focal_intersection, function(x){
        buffered_occurrences_intersection[[x]] # Identify all intersecting points for all points intersecting with the focal point
      }) %>%
        unlist() %>%
        unique()
      intersections <- paste0(paste0(" ", intersections, collapse = ", "), ", ")
      intersections
    })
    connections <- groups
    group_count <- 0
    while (length(connections) != group_count){
      group_count <- length(connections)
      connections <- purrr::map(1:length(connections), function(i){
        connections <- grep(paste0(paste0(strsplit(connections[[i]], ",")[[1]], ", "), collapse = "|"), groups)
        paste0(paste0(" ", connections, collapse = ", "), ", ")
      }) %>% unique()
    }
  
  eo_count <- connections %>% unique() %>% length()

  out <- list(buffered_occurrences = buffered_occurrences, eo_count = eo_count)
  return(out)
}

# Run full rank assessment
run_rank_assessment <- function(taxon_name, 
                                minimum_fields =  c("key", "scientificName", "prov", "longitude", "latitude", "coordinateUncertaintyInMeters", "stateProvince", "countryCode", "year", "month", "day", "institutionCode", "references"),
                                uploaded_data = NULL,
                                max_number_observations = 10000,
                                clean_occ = TRUE,
                                centroid_filter = FALSE,
                                date_start = "1900-01-01",
                                date_end = "2025-01-01",
                                uncertainty_filter = "",
                                nations_filter = NULL,
                                states_filter = NULL,
                                network_polys = network_polys,
                                sources_filter = NULL,
                                grid_cell_size = 2,
                                sep_distance = 1000
                                
){
  
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
    EOcount_map = NULL,
    EOcount_value = NULL
  )
  
  ns_table_full <- natserv::ns_search_spp(text_adv = list(searchToken = taxon_name, matchAgainst = "allScientificNames", operator="contains"))$results
  gbif_table <- rgbif::name_suggest(q = taxon_name, rank = c("species", "subspecies"), limit = 10)$data
  
  if (nrow(ns_table_full) == 0){
    ns_table_full <- natserv::ns_search_spp(text_adv = list(searchToken = gbif_table$canonicalName[1], matchAgainst = "allScientificNames", operator="contains"))$results
  }
  
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
  taxon_data$info_extended <- natserv::ns_id(uid = taxon_data$info$uniqueId %>% unique() %>% na.omit() %>% as.character() %>% head(1)) 
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
  
  taxon_data$gbif_occurrences_raw <- get_gbif_data(
    sp_data = taxon_data$synonyms, 
    number_observations = max_number_observations, 
    gbif = TRUE,
    inat = TRUE,
    ebird = TRUE
    )$sp_occurrences
  
  taxon_data$gbif_occurrences <- taxon_data$gbif_occurrences_raw %>% 
    clean_gbif_data(clean = clean_occ, remove_centroids = centroid_filter, minimum_fields = minimum_fields)

  taxon_data$uploaded_occurrences <- uploaded_data
  
  taxon_data$all_occurrences <- rbind(
    taxon_data$gbif_occurrences,
    taxon_data$uploaded_occurrences
  )
  
  shifted <- FALSE
  
  taxon_data$all_occurrences$longitude[taxon_data$all_occurrences$longitude > 180] <- taxon_data$all_occurrences$longitude[taxon_data$gbif_occurrences$longitude > 180] - 360
  max_long <- max(taxon_data$all_occurrences$longitude, na.rm = TRUE)/2
  shifted_long <- taxon_data$all_occurrences$longitude
  
  if (length(taxon_data$all_occurrences$longitude[taxon_data$all_occurrences$longitude > max_long]) > 0){
    shifted_long[shifted_long > max_long] <- shifted_long[shifted_long > max_long] - 360
    shifted_long <- shifted_long + 360
    # shifted <- TRUE
  }
  
  if ((max(shifted_long)-min(shifted_long)) < (max(taxon_data$all_occurrences$longitude) - min(taxon_data$all_occurrences$longitude))){
    taxon_data$all_occurrences$longitude <- shifted_long
    shifted <- TRUE
  }
  
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
    dplyr::filter(
      year >= substr(date_start, 1, 4) & year <= substr(date_end, 1, 4) | is.na(year),
      (month == substr(date_start, 6, 7) & day >= substr(date_start, 9, 10)) | (month > substr(date_start, 6, 7) & month < substr(date_end, 6, 7)) | (month == substr(date_end, 6, 7) & day <= substr(date_end, 9, 10)) | is.na(month) | is.na(day),
      )
  
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
  
  eoo_output <- taxon_data$sf_filtered %>% calculate_eoo(shifted = shifted)
  taxon_data$species_range_value <- eoo_output$EOO
  taxon_data$species_range_map <- eoo_output$hull
  taxon_data$AOO_value <- (aoo2(taxon_data$sf_filtered, as.numeric(grid_cell_size)*1000))/4
  taxon_data$AOO_map <- get_aoo_polys(taxon_data$sf_filtered, as.numeric(grid_cell_size))
  number_EOs <- calculate_number_occurrences(taxon_data$sf_filtered, separation_distance = sep_distance %>% as.numeric(), added_distance = 0)
  taxon_data$EOcount_value <- number_EOs$eo_count
  taxon_data$EOcount_map <- number_EOs$buffered_occurrences
  taxon_data$shifted <- shifted 
  
  return(taxon_data)
  
}

# Function to extract id from table button click
shinyInput <- function(FUN, len, id, ...) {
  inputs <- character(len)
  for (i in seq_len(len)) {
    inputs[i] <- as.character(FUN(paste0(id, i), ...))
  }
  inputs
}

# Get temporal trends
get_temporal_trends <- function(taxon_data = taxon_data, referenceTaxon = "kingdom", start_year = 1980){
  
  min_fields <- c("key", "scientificName", "prov", "longitude", "latitude", "coordinateUncertaintyInMeters", "stateProvince", "countryCode", "year", "institutionCode", "references")

  key_value <- taxon_data$synonyms %>% 
    dplyr::filter(datasetKey == "d7dddbf4-2cf0-4f39-9b2a-bb099caae36c") %>% 
    dplyr::select(paste0(referenceTaxon, "Key")) %>% 
    as.character()
  
  if (is.null(taxon_data$AOO_map)){
    taxon_data$AOO_map <- get_aoo_polys(taxon_data$sf_filtered, 2)
  }
  
  gbif_data <- spocc::occ(from = "gbif", gbifopts = list(
    taxonKey = key_value,
    geometry = taxon_data$AOO_map %>% wk::as_wkt(),
    year = paste0(start_year, ",", substr(Sys.Date(), 1, 4))
  ),
  limit = 1000000, 
  has_coords = TRUE
  )$gbif$data$custom_query
  
  reference_taxon_name <- gbif_data[[referenceTaxon]][1]
  
  gbif_data <- gbif_data %>% 
    dplyr::filter(
      complete.cases(longitude, latitude),
      !(scientificName %in% unique(taxon_data$sf_filtered$scientificName))
    ) %>% 
    dplyr::select(intersect(names(gbif_data), min_fields)) %>% 
    dplyr::mutate(lon = longitude,
                  lat = latitude) %>% 
    sf::st_as_sf(
      coords = c("lon", "lat"),
      crs = 4326
    )
  
  gbif_data <- gbif_data %>% 
    rbind(taxon_data$sf_filtered %>% 
            dplyr::filter(year >= start_year) %>% 
            dplyr::mutate(scientificName = "focal_taxon")) %>% 
    dplyr::distinct(key, .keep_all = TRUE) 
  
  gbif_data <- gbif_data %>% 
    dplyr::mutate(cellID = sf::st_intersects(gbif_data, taxon_data$AOO_map, sparse = TRUE) %>% purrr::map(function(x) ifelse(!is.null(x), x, NA)) %>% unlist())
  
  gbif_data <- gbif_data %>% 
    dplyr::filter(complete.cases(cellID, year)) %>% 
    dplyr::mutate(visitID = paste0(cellID, "_", year))
  
  ### Get recording data from observations across all species
  count_data <- gbif_data %>%
    sf::st_set_geometry(NULL) %>% 
    dplyr::group_by(visitID, cellID, year) %>%
    dplyr::count(scientificName) %>%
    dplyr::ungroup() %>%
    tidyr::spread(key = scientificName, value = n, fill = 0) 
  
  detection_data <- count_data %>% dplyr::mutate(across(4:ncol(count_data), ~as.numeric(. > 0)))
  
  focal_species_detection_data <- count_data %>% 
    dplyr::mutate(focal_species_count = .data[["focal_taxon"]],
                  focal_species_detection = ifelse(focal_species_count > 0, 1, 0),
                  total_number_observations = count_data %>% dplyr::select(-visitID, -cellID, -year) %>% rowSums(),
                  species_list_length = detection_data %>% dplyr::select(-visitID, -cellID, -year) %>% rowSums()
    ) %>% 
    dplyr::select(visitID, cellID, year, focal_species_count, focal_species_detection, total_number_observations, species_list_length)
  
  focal_species_detection_data <- focal_species_detection_data %>% 
    dplyr::mutate(
      proportion_observations = focal_species_count/total_number_observations,
      proportion_species = focal_species_detection/species_list_length,
    )
  
  # Mean yearly proportions
  temporal_trend_data <- focal_species_detection_data %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(
      focal_species_count = sum(focal_species_count),
      total_number_observations = sum(total_number_observations),
      proportion_observations = focal_species_count/total_number_observations, # mean(proportion_observations),
      proportion_species = mean(proportion_species),
      focal_species_detection = sum(focal_species_detection),
      n_visits = n()
    ) %>%
    dplyr::ungroup() %>% 
    dplyr::mutate(
      reporting_rate = focal_species_detection/n_visits
    )
  
  # Number of focal observations
  p1 <- ggplot(data = temporal_trend_data, aes(x = year, y = focal_species_count)) + 
    geom_col(col = "#2c7bb680", alpha = 0.5) + 
    ylab("") +
    xlab("") +
    ylim(c(0, max(temporal_trend_data$focal_species_count)+(0.2 *max(temporal_trend_data$focal_species_count)))) +
    annotate("text", x = median(temporal_trend_data$year), y = max(temporal_trend_data$focal_species_count)+(0.08 *max(temporal_trend_data$focal_species_count)), label = paste0("Observations of ", taxon_data$info$scientificName)) +
    theme_linedraw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 9),
          axis.text.x = element_blank(),
          axis.text = element_text(size = 8)
    )
  # Total number of observations of reference taxon
  p2 <- ggplot(data = temporal_trend_data, aes(x = year, y = total_number_observations)) + 
    geom_col(col = "#2c7bb680", alpha = 0.5) + 
    ylim(c(0, max(temporal_trend_data$total_number_observations)+(0.2 *max(temporal_trend_data$total_number_observations)))) +
    ylab("") +
    xlab("") +
    annotate("text", x = median(temporal_trend_data$year), y = max(temporal_trend_data$total_number_observations)+(0.08 *max(temporal_trend_data$total_number_observations)), label = paste0("Observations of ", referenceTaxon, " ", reference_taxon_name)) +
    theme_linedraw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 9),
          axis.text.x = element_blank(),
          axis.text = element_text(size = 8)        
    )
  # Mean Observation Rate
  p3 <- ggplot(data = temporal_trend_data, aes(x = year, y = proportion_observations)) +
    geom_col(col = "#2c7bb680", alpha = 0.5) +
    geom_smooth(size = 1.5, se = FALSE, col = "black") +
    ylab("") +
    xlab("") +
    annotate("text", x = median(temporal_trend_data$year), y = max(temporal_trend_data$proportion_observations)+(0.08 *max(temporal_trend_data$proportion_observations)), label = "Proportion of Observations") +
    ylim(c(0, max(temporal_trend_data$proportion_observations)+(0.2 *max(temporal_trend_data$proportion_observations)))) +
    theme_linedraw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 9),
          axis.text.x = element_blank(),
          axis.text = element_text(size = 8)
          )
  # Mean proportion of species
  p4 <- ggplot(data = temporal_trend_data, aes(x = year, y = proportion_species)) +
    geom_col(col = "#2c7bb680", alpha = 0.5) +
    geom_smooth(size = 1.5, se = FALSE, col = "black") +
    ylab("") +
    xlab("") +
    annotate("text", x = median(temporal_trend_data$year), y = max(temporal_trend_data$proportion_species)+(0.08 *max(temporal_trend_data$proportion_species)), label = "Proportion of Species") +
    ylim(c(0, max(temporal_trend_data$proportion_species)+(0.2 *max(temporal_trend_data$proportion_species)))) +
    theme_linedraw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 9),
          axis.text.x = element_blank(),
          axis.text = element_text(size = 8)
    )
  
  # Yearly reporting rate
  p5 <- ggplot(data = temporal_trend_data, aes(x = year, y = reporting_rate)) + 
    geom_col(col = "#2c7bb680", alpha = 0.5) + 
    geom_smooth(size = 1, se = FALSE, col = "black") + 
    ylim(c(0, max(temporal_trend_data$reporting_rate)+(0.2 *max(temporal_trend_data$reporting_rate)))) +
    ylab("") +
    annotate("text", x = median(temporal_trend_data$year), y = max(temporal_trend_data$reporting_rate)+(0.08 *max(temporal_trend_data$reporting_rate)), label = "Observation Rate") +
    theme_linedraw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 9),
          axis.text.x = element_blank(),
          axis.text = element_text(size = 8)
    )
  ## GLMER (with spatial random effect)
  options(na.action = "na.fail")
  glmer_result <- lme4::glmer(focal_species_detection ~ splines::bs(year, df = 3) + total_number_observations + species_list_length + (1 | cellID), 
                              data = focal_species_detection_data, family = binomial, control = lme4::glmerControl(tol = 1e-5, optimizer = "bobyqa", optCtrl=list(maxfun=2e5))
  )
  # Modeled probability of detection
  p6 <- data.frame(fitted = fitted(glmer_result), year = focal_species_detection_data$year) %>% 
    ggplot(aes(x = year, y = fitted)) + 
    geom_point(col = "#2c7bb680", alpha = 0.5) + 
    geom_smooth(size = 1, se = FALSE, col = "black") + 
    # ylim(c(0, 1.2)) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1.2)) +
    xlab("Year") +
    ylab("") +
    annotate("text", x = median(focal_species_detection_data$year), y = 1.08, label = "Modeled probability of detection") +
    theme_linedraw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 8)
    )
  
  p <- subplot(p1, p2, p3, p5, p6, nrows = 5, shareX = TRUE, titleX = TRUE, margin = 0.005, which_layout = 1)
  
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
  
  return(gg)
  
}
###########################################################
#Simple area Projection Wizard
###########################################################
#' Simple equal area projection wizard
#' @title Simple Projection Wizard
#' @description 
#' Projects any set of lat long points to a "suitable" area projection, based on it's "true centre of gravity"
#' @author Justin Moat. J.Moat@kew.org
#' @note
#' Based around a simple continental projection, using two sets of projections
#' equal area cylindrical = Cylindrical equal-area = 8287
#' equal area azimuthal for polar (above 70) = Lambert azimuthal equal-area
#'  
#' note these are not cartographically pleasing projections, they are just so we can get the data into something simple for areal analysis
#' See below for a more cartographically pleasing projection engine
#' 
#' Šavric, B., Jenny, B., Jenny, H., 2016. Projection Wizard – An Online Map Projection Selection Tool. Cartogr. J. 53, 1–9. doi:10.1080/00087041.2015.1131938
#' @param thepoints set of points in latitude and longtitude ie c(lat,long)
#' @param thecentre one point ie c(lat,long)
#' 
#' @return set of points in metres (x,y)
#' @examples 
#'lat <- runif (200,-24,-12)
#'long <- runif (200,43,51)
#'ll <- data.frame(lat,long)
#'cp <- trueCOGll(ll)
#'pointsprojected <- simProjWiz(ll,cp)
#' @references 
#' Šavric, B., Jenny, B., Jenny, H., 2016. Projection Wizard – An Online Map Projection Selection Tool. Cartogr. J. 53, 1–9. doi:10.1080/00087041.2015.1131938
#' 
#' Snyder, J.P., 1987. Map projections: A working manual, Professional Paper. Washington, D.C.
#' @export
#' @import sp
#' @import rgdal



######################################################################
#simple projection wizard#
######################################################################
#determining projection around on center of points
#based around a simplied continental scheme. Also see:
#Šavric, B., Jenny, B., Jenny, H., 2016. Projection Wizard – An Online Map Projection Selection Tool. Cartogr. J. 53, 
#1–9. doi:10.1080/00087041.2015.1131938
#two sets of projections
#equal area cylindrical = Cylindrical equal-area = 8287
#equal area azimuthal for polar (above 70) = Lambert azimuthal equal-area = 
#note these are not cartographically pleasing projections, they are just so we can get the data into something simple for areal analysis
######################################################################
## NEEDS TO BE UPDATED TO USE SF PACKAGE
simProjWiz <- function(thepoints, thecentre){
  #setup and set projection to WGS84
  sp::coordinates(thepoints) <- c("long", "lat")
  sp::proj4string(thepoints) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # geographic and WGS 84
  #depending on centre point
  if((thecentre$lat < 70) & (thecentre$lat > -70)){
    CRSstring <- paste("+proj=cea +lon_0=", thecentre$long,   " +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",sep = "")
  } else {
    CRSstring <- paste("+proj=laea +lat_0=", thecentre$lat," +lon_0=", thecentre$long, " +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",sep = "")
  }
  CRS.new <- sp::CRS(CRSstring)
  #reproject
  xysp <- sp::spTransform(thepoints, CRS.new)
  xy <- as.data.frame(xysp)
  #rename to x and y as not longer lat long
  colnames (xy) <- c("x","y")
  return(list(xy = xy, crs_new = CRS.new))
}

######################################################################
#calculates 'true' centre of gravity from a set of lat long points in decimal degrees   #
#note z from mean of cartesian give some idea of weighted spread on the globe#
######################################################################
#' @title True centre of gravity from a set of Lat longs
#' @description 
#' Calculates the "true" centre of gravity (weighted) from a set of lat longs, using cartesian geometry
#' @author Justin Moat. J.Moat@kew.org
#' @param thepoints set of points c(lat,long)
#' @return a point (lat,long) from centre
#' @examples 
#'lat <- runif (200,-24,-12)
#'long <- runif (200,43,51)
#'ll <- data.frame(lat,long)
#'cp <- trueCOGll(ll)
#' @references Descartes, R., 1637. Discours de la methode. A Leyde, De l’imprimerie de I. Maire, Paris.
#' @export
trueCOGll <-function(thepoints){
  
  llrad <- deg2rad(thepoints) #to radians
  cartp <- ll2cart(llrad$lat,llrad$long) #to cartesian
  mp <- data.frame(x=mean(cartp$x),y=mean(cartp$y),z=mean(cartp$z)) #central point
  pmp <- pro2sph(mp$x,mp$y,mp$z) #projection to surface
  pmprll <- cart2ll(pmp$x,pmp$y,pmp$z) #to ll in radians
  pmpll <- rad2deg(pmprll) #to degrees
  return(data.frame(lat=pmpll$latr,long=pmpll$longr))
  
}

######################################################################
#calculates the Cartesian cordinates (x,y,z) from lat long in radians#
######################################################################
#' @title Geographic coordinates to cartesian (x,y,z)
#' @description 
#' Calculates the Cartesian cordinates (x,y,z) from lat long in radians
#' @author Justin Moat. J.Moat@kew.org
#' @param latr latitude point in radians
#' @param longr longtitude point in radians
#' @return dataframe of x,y,z
#' @examples 
#'lat <- runif (200,-24,-12)
#'long <- runif (200,43,51)
#'thepoints <- data.frame(lat,long)
#'llrad <- deg2rad(thepoints)
#'cartp <- ll2cart(llrad$lat,llrad$long)
#' @references Descartes, R., 1637. Discours de la methode. A Leyde, De l’imprimerie de I. Maire, Paris.
#' @export

ll2cart <- function(latr,longr){
  x <- cos(latr) * cos(longr)
  y <- cos(latr) * sin(longr)
  z <- sin(latr)
  return(data.frame(x,y,z))
}

######################################################################
#calculates the lat long cordinates in radians from Cartesian (x,y,z)#
######################################################################
#' @title Cartesian (x,y,z) to Geographic coordinates
#' @description 
#' calculates the latitude and longtitude cordinates in radians from Cartesian coordinates (x,y,z)
#' @author Justin Moat. J.Moat@kew.org
#' @param x East to West coordinate in metres
#' @param y South to North coordinate in metres
#' @param z height coordinate in metres
#' @return dataframe of latitude,longtitude
#' @export

cart2ll <-function (x,y,z){
  latr <- asin(z)
  longr <- atan2(y,x)
  return(data.frame(latr,longr))
}


######################################################################
#calculates Cartesian (x,y,z), projected from the centre of the sphere 
#to the earth surface, returns cartesian (x,y,z)
#used to calculate "true" centre of set of lat longs
# http://stackoverflow.com/questions/9604132/how-to-project-a-point-on-to-a-sphere
######################################################################
#' @title Cartesian coordinate projection
#' @description 
#' calculates Cartesian (x,y,z), projected from the centre of the sphere 
#' to the earth surface, returns cartesian (x,y,z)
#' used to calculate "true" centre of set of lat longs
#' @author Justin Moat. J.Moat@kew.org
#' @note
#' http://stackoverflow.com/questions/9604132/how-to-project-a-point-on-to-a-sphere
#' @param x East to West coordinate in metres
#' @param y South to North coordinate in metres
#' @param z height coordinate in metres
#' @return x,y,z
#' @references Descartes, R., 1637. Discours de la methode. A Leyde, De l’imprimerie de I. Maire, Paris.
#' @export

pro2sph <- function (x,y,z){
  sc <- 1/sqrt(x^2 + y^2 + z^2)
  x <- x * sc
  y <- y * sc
  z <- z * sc
  return(data.frame(x,y,z))
}

######################################################################
#radians to degrees and degrees to radians
######################################################################
#' @title Radians to Degrees
#' @description 
#' Calculates radians from degrees or degrees from radians
#' @author Justin Moat. J.Moat@kew.org
#' @param rad number in radians
#' @return number
#' @examples 
#' b <- 0.392699
#' rad2deg(b)
#' @export

rad2deg <- function(rad) {(rad * 180) / (pi)}

######################################################################
#radians to degrees and degrees to radians
######################################################################
#' @title 
#' Degrees to radians
#' @description 
#' Calculates radians from degrees or degrees from radians
#' @author Justin Moat. J.Moat@kew.org
#' @param deg number in degrees
#' @return number
#' @examples 
#' a <- 30
#' deg2rad(a)
#' @export

deg2rad <- function(deg) {(deg * pi) / (180)}