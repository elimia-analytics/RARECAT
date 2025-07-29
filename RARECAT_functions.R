# Function to extract observations from GBIF for target taxon
get_gbif_data <- function(taxa_metadata, 
                          datasets_metadata,
                          query_polygon = taxon_data$assessment_polygon,
                          all_occ_data = FALSE,
                          all_humobs_data = FALSE,
                          shift_occurrences = FALSE
                          ){

  gbif_occurrences <- gbif_occurrences_occ <- gbif_occurrences_humobs <- NULL

  if ("gbif" %in% datasets_metadata$datasetKey){
    
    occ_counts <- taxa_metadata$occurrence_count

    max_occ <- 5000
    if (sum(occ_counts >= 5000) == 1 & median(occ_counts) < 1000) max_occ <- 5000/sum(occ_counts > 5000)
    if (sum(occ_counts >= 5000) > 1 & median(occ_counts) < 1000) max_occ <- 5000/sum(occ_counts > 5000)
    if (sum(occ_counts >= 5000) > 1 & median(occ_counts) >= 1000) max_occ <- 5000/sum(occ_counts > median(x3))

    if (max(as.numeric(occ_counts)) < 5000){
    gbif_occurrences <- spocc::occ(from = "gbif",
                                   gbifopts = list(
                                     taxonKey = taxa_metadata$key,
                                     gadmGid = query_polygon
                                   ),
                                   limit = max(as.numeric(occ_counts)), 
                                   has_coords = TRUE
    )
    gbif_occurrences <- gbif_occurrences$gbif$data$custom_query
    } else {
      gbif_occurrences <- purrr::map2(
        list(
          c(paste0((substr(Sys.Date(), 1, 4) %>% as.numeric()) - 10, "-01-01"), Sys.Date() %>% as.character()),
          c(paste0((substr(Sys.Date(), 1, 4) %>% as.numeric()) - 20, "-01-01"), paste0((substr(Sys.Date(), 1, 4) %>% as.numeric()) - 11, "-12-31")),
          c(paste0((substr(Sys.Date(), 1, 4) %>% as.numeric()) - 30, "-01-01"), paste0((substr(Sys.Date(), 1, 4) %>% as.numeric()) - 21, "-12-31")),
          c(paste0((substr(Sys.Date(), 1, 4) %>% as.numeric()) - 40, "-01-01"), paste0((substr(Sys.Date(), 1, 4) %>% as.numeric()) - 31, "-12-31")),
          c(paste0((substr(Sys.Date(), 1, 4) %>% as.numeric()) - 50, "-01-01"), paste0((substr(Sys.Date(), 1, 4) %>% as.numeric()) - 41, "-12-31"))
        ),
        c(0.4*max_occ, 0.2*max_occ, 0.1*max_occ, 0.1*max_occ, 0.2*max_occ),# c(2000, 1000, 500, 500, 1000),
        function(x, y){
        out <- spocc::occ(from = "gbif",
                                     gbifopts = list(
                                       taxonKey = taxa_metadata$key,
                                       gadmGid = query_polygon
                                     ),
                                     limit = y, date = x,
                                     has_coords = TRUE
      )
        out$gbif$data$custom_query
        }
      ) %>% bind_rows()
    }
    
    dataset_keys <- gbif_occurrences$datasetKey %>% table() %>% sort() %>% head(50) %>% names() %>% unique()

    datasets_details <- data.frame(
      datasetKey = dataset_keys,
      datasetName = purrr::map_chr(dataset_keys, function(k) rgbif::dataset_get(k)$title)
    )

    gbif_occurrences <- gbif_occurrences %>%
      dplyr::select(-any_of("datasetName")) %>% 
      dplyr::left_join(datasets_details, by = "datasetKey")

  } else {
    
    if (all_occ_data){
      
      gbif_occurrences_occ <- purrr::map(taxa_metadata$key, function(k){
        out <- spocc::occ(from = "gbif",
                          gbifopts = list(
                            taxonKey = k,
                            basisOfRecord = c("OCCURRENCE", "PRESERVED_SPECIMEN", "OBSERVATION", "MACHINE_OBSERVATION"),
                            gadmGid = query_polygon
                          ),
                          limit = ifelse(all_humobs_data, 2500, 5000),
                          has_coords = TRUE
        )
        out$gbif$data$custom_query
      }) %>% bind_rows()
      
      if (nrow(gbif_occurrences_occ) > 0){
        if (!is.null(gbif_occurrences)){
          shared_names <- intersect(names(gbif_occurrences), names(gbif_occurrences_occ))
          gbif_occurrences <- gbif_occurrences %>% 
            dplyr::select(all_of(shared_names)) %>% 
            rbind(gbif_occurrences_occ %>% dplyr::select(all_of(shared_names)))
        } else {
          gbif_occurrences <- rbind(gbif_occurrences, gbif_occurrences_occ)
        }
      }
    } else {
      
      if (!("detailed" %in% datasets_metadata$datasetKey)){
        
      occ_data <- datasets_metadata %>% dplyr::filter(basisOfRecord == "OCCURRENCE")
      
      if (nrow(occ_data) > 0){
        gbif_occurrences_occ <- purrr::map(taxa_metadata$key, function(k){
          out2 <- purrr::map(c("OCCURRENCE", "PRESERVED_SPECIMEN", "OBSERVATION", "MACHINE_OBSERVATION"), function(d){
            out1 <- spocc::occ(from = "gbif",
                               gbifopts = list(
                                 taxonKey = k,
                                 datasetKey = occ_data$datasetKey,
                                 basisOfRecord = d,
                                 gadmGid = query_polygon
                               ),
                               limit = max(occ_data$recordsMax) %>% as.numeric(),
                               has_coords = TRUE
            )
            out1$gbif$data$custom_query
          }) %>% bind_rows()
          out2
        }) %>% bind_rows()

      if (nrow(gbif_occurrences_occ) > 0){
        if (!is.null(gbif_occurrences)){
          shared_names <- intersect(names(gbif_occurrences), names(gbif_occurrences_occ))
          gbif_occurrences <- gbif_occurrences %>% 
            dplyr::select(all_of(shared_names)) %>% 
            rbind(gbif_occurrences_occ %>% dplyr::select(all_of(shared_names)))
        } else {
          gbif_occurrences <- rbind(gbif_occurrences, gbif_occurrences_occ)
        }
      }
      }
      
      } else {
        gbif_occurrences_occ <- NULL
        gbif_occurrences <- rbind(gbif_occurrences, gbif_occurrences_occ)
      }
    }
    
    if (all_humobs_data){
      
      gbif_occurrences_humobs <- spocc::occ(from = "gbif",
                                            gbifopts = list(
                                              taxonKey = taxa_metadata$key,
                                              basisOfRecord = "HUMAN_OBSERVATION"
                                            ),
                                            limit = ifelse(all_occ_data, 2500, 5000),
                                            has_coords = TRUE
      )
      
      gbif_occurrences_humobs <- gbif_occurrences_humobs$gbif$data$custom_query
      
      if (nrow(gbif_occurrences_humobs) > 0){
        if (!is.null(gbif_occurrences)){
          shared_names <- intersect(names(gbif_occurrences), names(gbif_occurrences_humobs))
          gbif_occurrences <- gbif_occurrences %>% 
            dplyr::select(all_of(shared_names)) %>% 
            rbind(gbif_occurrences_humobs %>% dplyr::select(all_of(shared_names)))
        } else {
          gbif_occurrences <- rbind(gbif_occurrences, gbif_occurrences_humobs)
        }
      }
      
    } else {
      
      if (!("detailed" %in% datasets_metadata$datasetKey)){
        
      humobs_data <- datasets_metadata %>% dplyr::filter(basisOfRecord == "HUMAN_OBSERVATION")
      
      if (nrow(humobs_data) > 0){
        gbif_occurrences_humobs <- purrr::map(taxa_metadata$key, function(k){
          out <- spocc::occ(from = "gbif",
                            gbifopts = list(
                              taxonKey = k,
                              datasetKey = humobs_data$datasetKey,
                              basisOfRecord = "HUMAN_OBSERVATION",
                              gadmGid = query_polygon
                            ),
                            limit = max(humobs_data$recordsMax) %>% as.numeric(),
                            has_coords = TRUE
          )
          out$gbif$data$custom_query
        }) %>% bind_rows()

    if (nrow(gbif_occurrences_humobs) > 0){
      if (!is.null(gbif_occurrences)){
        shared_names <- intersect(names(gbif_occurrences), names(gbif_occurrences_humobs))
        gbif_occurrences <- gbif_occurrences %>% 
          dplyr::select(all_of(shared_names)) %>% 
          rbind(gbif_occurrences_humobs %>% dplyr::select(all_of(shared_names)))
      } else {
        gbif_occurrences <- rbind(gbif_occurrences, gbif_occurrences_humobs)
      }
    }
      }
      } else {
        gbif_occurrences_humobs <- NULL
        gbif_occurrences <- rbind(gbif_occurrences, gbif_occurrences_humobs)
      }
    }
    
    if (!is.null(gbif_occurrences) & !("detailed" %in% datasets_metadata$datasetKey)){
      gbif_occurrences <- gbif_occurrences %>%
        dplyr::select(-any_of("datasetName")) %>% 
        dplyr::left_join(datasets_metadata %>% dplyr::select(datasetKey, datasetName), by = "datasetKey")

    }
  }

    gbif_occurrences <- gbif_occurrences %>% 
    # dplyr::filter(basisOfRecord %in% c("OCCURRENCE", "PRESERVED_SPECIMEN", "OBSERVATION", "MACHINE_OBSERVATION")) %>%
    dplyr::filter(complete.cases(latitude, longitude, basisOfRecord)) %>% # Only keep records that have at a minimum a value for longitude, latitude, and basisOfRecord
    dplyr::filter(latitude != 0 | longitude != 0) %>% # Exclude records that have latitude and longitude values of 0
    dplyr::mutate(references = paste0("https://www.gbif.org/occurrence/", key)) # Clean up references field

  sp_occurrences <- gbif_occurrences

  shifted <- FALSE
  
  if (shift_occurrences){

    sp_occurrences$longitude[sp_occurrences$longitude > 180] <- sp_occurrences$longitude[sp_occurrences$longitude > 180] - 360
    # sp_occurrences$longitude[sp_occurrences$longitude < -180] <- sp_occurrences$longitude[sp_occurrences$longitude < -180] + 360
    #
    max_long <- max(sp_occurrences$longitude, na.rm = TRUE)/2
    shifted_long <- sp_occurrences$longitude

    if (length(sp_occurrences$longitude[sp_occurrences$longitude > max_long]) > 0){
      shifted_long[shifted_long > max_long] <- shifted_long[shifted_long > max_long] - 360
      shifted_long <- shifted_long + 360
    }

    if ((max(shifted_long)-min(shifted_long)) < (max(sp_occurrences$longitude) - min(sp_occurrences$longitude))){
      sp_occurrences$longitude <- shifted_long
      shifted <- TRUE
    }

}

  out <- list(sp_occurrences = sp_occurrences, shifted = shifted)
  
  return(out)
  
}

# Function to clean up GBIF records
clean_gbif_data <- function(gbif_occurrences, clean = TRUE, minimum_fields = minimum_fields, remove_centroids = FALSE){
  
  if (isTRUE(clean)){
    
    if (nrow(gbif_occurrences) >= 2){
      
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
    if ("georeferenceProtocol" %in% names(gbif_occurrences)) gbif_occurrences <- gbif_occurrences %>% dplyr::filter(!grepl("centroid|Centroid|CENTROID", georeferenceProtocol))
    
  }

  gbif_occurrences <- gbif_occurrences %>%
    dplyr::select(all_of(intersect(names(gbif_occurrences), minimum_fields))) %>%
    cbind(matrix(NA, nrow = nrow(gbif_occurrences), ncol = length(setdiff(minimum_fields, names(gbif_occurrences)))) %>%
            as.data.frame() %>%
            set_names(setdiff(minimum_fields, names(gbif_occurrences)))
    ) %>%
    dplyr::select(all_of(minimum_fields))

  return(gbif_occurrences)
}

# Function to load and process user data
process_user_data <- function(user_file, minimum_fields){
  
  user_data <- read.csv(user_file, header = TRUE) 
  user_data <- user_data %>% select_if(~!(all(is.na(.)) | all(. == "")))
  user_data <- user_data %>% mutate_if(is.character, na_if, "")
  processed_data <- NULL
  longitude_names <- c("longitude", "LONGITUDE", "decimalLongitude", "private_longitude", "Decimal Longitude", "Longitude", "lon", "Lon", "X", "x")
  longitude_names <- longitude_names %>% set_names(rep("longitude", length(longitude_names)))
  latitude_names <- c("latitude", "LATITUDE", "decimalLatitude", "private_latitude", "Decimal Latitude", "Latitude", "lat", "Lat", "Y", "y")
  latitude_names <- latitude_names %>% set_names(rep("latitude", length(latitude_names)))
  
  if (sum(grepl(paste0(longitude_names, collapse = "|"), names(user_data))) > 0 & sum(grepl(paste0(latitude_names, collapse = "|"), names(user_data))) > 0){
    
    scientificName_Source_names <- c("scientificName_Source", "scientificName", "Scientific name", "GNAME", "scientific_name", "SciName", "scientific name", "Scientific Name")
    scientificName_Source_names <- scientificName_Source_names %>% set_names(rep("scientificName", length(scientificName_Source_names)))
    scientificName_Assessment_names <- c("scientificName_Assessment")
    scientificName_Assessment_names <- scientificName_Assessment_names %>% set_names(rep("scientificName_Assessment", length(scientificName_Assessment_names)))
    stateProvince_names <- c("stateProvince", "STATE", "GNAME", "place_state_name", "State", "State or Province", "state", "Province", "province", "PROVINCE", "place_guess")
    stateProvince_names <- stateProvince_names %>% set_names(rep("stateProvince", length(stateProvince_names)))
    countryCode_names <- c("countryCode", "NATN", "country", "place_country_name", "Country", "COUNTRY", "Nation", "NATION")
    countryCode_names <- countryCode_names %>% set_names(rep("countryCode", length(countryCode_names)))
    year_names <- c("year", "Year", "YEAR", "Year Collected", "SiteDate")
    year_names <- year_names %>% set_names(rep("year", length(year_names)))
    coordinateUncertaintyInMeters_names <- c("coordinateUncertaintyInMeters", "public_positional_accuracy")
    coordinateUncertaintyInMeters_names <- coordinateUncertaintyInMeters_names %>% set_names(rep("coordinateUncertaintyInMeters", length(coordinateUncertaintyInMeters_names)))
    EORANK_names <- c("EO RANK", "EORANK", "EORANK_CD", "EO_RANK")
    EORANK_names <- EORANK_names %>% set_names(rep("EORANK", length(EORANK_names)))
    lookup <- c(longitude_names, latitude_names, scientificName_Source_names, scientificName_Assessment_names, stateProvince_names, countryCode_names, year_names, EORANK_names)
    
    if (("coordinates_obscured" %in% names(user_data)) & ("private_latitude" %in% names(user_data))){
      
      user_data <- user_data %>% 
        dplyr::mutate(
          latitude = ifelse(coordinates_obscured == "true" & !is.null(private_latitude), private_latitude, latitude),
          longitude = ifelse(coordinates_obscured == "true" & !is.null(private_longitude), private_longitude, longitude)
        )
      
    }
    
    if ("latitude" %in% names(user_data)){
      user_data <- user_data %>% 
        dplyr::select(-any_of(setdiff(latitude_names, "latitude")))
    }
    
    if ("longitude" %in% names(user_data)){
      user_data <- user_data %>% 
        dplyr::select(-any_of(setdiff(longitude_names, "longitude")))
    }
    
    if (length(intersect(countryCode_names, names(user_data))) > 1){
      user_data <- user_data %>% 
        dplyr::select(-intersect(countryCode_names, names(user_data))[2])
    }
    
    if ("observed_on" %in% names(user_data)){
      user_data <- user_data %>% 
        dplyr::mutate(year = substr(observed_on, 1, 4) %>% as.numeric())
    }

    target_columns <- user_data %>% 
      dplyr::select(matches(paste0(longitude_names, collapse = "|")), 
                    matches(paste0(latitude_names, collapse = "|")), 
                    matches(paste0(scientificName_Source_names, collapse = "|")),
                    matches(paste0(scientificName_Assessment_names, collapse = "|")),
                    matches(paste0(stateProvince_names, collapse = "|")),
                    matches(paste0(countryCode_names, collapse = "|")),
                    matches(paste0(year_names, collapse = "|")),
                    matches(paste0(EORANK_names, collapse = "|"))
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
      dplyr::select(intersect(names(processed_data), minimum_fields)) 
    
    user_columns_to_keep <- user_data %>% 
      dplyr::select(setdiff(intersect(names(user_data), minimum_fields), names(processed_data))) 
    
    if (ncol(user_columns_to_keep) > 0){
      processed_data <- processed_data %>% 
        cbind(user_columns_to_keep)
    }
    
    processed_data <- processed_data %>%
      cbind(matrix(NA, nrow = nrow(processed_data), ncol = length(setdiff(minimum_fields, names(processed_data)))) %>%
              as.data.frame() %>%
              set_names(setdiff(minimum_fields, names(processed_data)))
      ) %>%
      dplyr::select(all_of(minimum_fields)) %>%
      dplyr::distinct(., .keep_all = TRUE) %>%    
      dplyr::filter(complete.cases(latitude, longitude), latitude != 0, longitude != 0) 
    
    processed_data <- processed_data %>% 
      dplyr::mutate(prov = "uploaded",
                    year = ifelse(!is.na(year), year, NA), #format(Sys.time(), "%Y")),
                    scientificName = ifelse(!is.na(scientificName), scientificName, "user-uploaded")
      )
    
    if ("EO_ID" %in% names(user_data)){
      processed_data$basisOfRecord <- "ELEMENT_OCCURRENCE"
    }
  
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
  if (shifted){
    scaling_factor <- max(occurrences_sf$longitude)-min(occurrences_sf$longitude)
    # occurrences_sf$geometry <- (sf::st_geometry(occurrences_sf) + c(360,90)) %% c(360) - c(0,90)
    occurrences_sf$geometry <- (sf::st_geometry(occurrences_sf)-c(scaling_factor, 0))
    st_crs(occurrences_sf) <- 4326
    print(summary(st_coordinates(occurrences_sf)[, 1]))
  }

  hull <- occurrences_sf %>%  #[vertices, ] %>% 
    terra::vect() %>% 
    terra::convHull() %>% 
    sf::st_as_sf()
  st_crs(hull) <- 4326
  
  if (min(st_coordinates(occurrences_sf)[,1]) > 0){
  long_range <- abs(max(st_coordinates(occurrences_sf)[,1], na.rm = TRUE))-abs(min(st_coordinates(occurrences_sf)[,1], na.rm = TRUE))
  } else {
    long_range <- abs(max(st_coordinates(occurrences_sf)[,1], na.rm = TRUE))+abs(min(st_coordinates(occurrences_sf)[,1], na.rm = TRUE))
  }
  
  if (shifted & long_range >= 180){
    EOO <- purrr::safely(red::eoo)(st_coordinates(occurrences_sf))
  } else {
    centreofpoints <- trueCOGll(st_coordinates(occurrences_sf) %>% as.data.frame() %>% set_names("longitude", "latitude"))
    mypointsxy <- simProjWiz(st_coordinates(occurrences_sf) %>% as.data.frame() %>% set_names(c("long", "lat")), centreofpoints)
    EOO <- purrr::safely(red::eoo)(mypointsxy$xy)
  }
  
  if (shifted){
  # Rescale hull for mapping
  hull$geometry <- (sf::st_geometry(hull)+c(scaling_factor, 0))
  st_crs(hull) <- 4326
  }
  
  if (!is.null(EOO$result)){
    EOO <- EOO$result %>% units::set_units(km^2)
  } else {
    EOO <- NA
  }
  
  if (!is.na(EOO)){
    eoo_factor <- cut(as.numeric(EOO), breaks = c(0, 0.999, 99.999, 249.999, 999.999, 4999.999, 19999.999, 199999.999, 2499999.999, 1000000000), labels = c("Z", LETTERS[1:8]))
  } else {
    eoo_factor <- NA
  }
    
  out <- list(hull = hull, EOO = EOO, factor = eoo_factor)
  
  return(out)
}

# Safe version of calculate EOO
safe_eoo <- purrr::safely(calculate_eoo)

# Extract AOO polygons
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

# Function to derive AOO factor value from AOO value for a given grid cell size
get_aoo_factor <- function(AOO_value, grid_cell_size = 2){
  
  if (grid_cell_size == 1) out <- base::cut(req(as.numeric(AOO_value)), breaks = c(0, 0.999, 4.999, 10.999, 20.999, 100.999, 500.999, 2000.999, 10000.999, 50000.999, 1000000000), labels = c("Z", LETTERS[1:9]))
  
  if (grid_cell_size > 1) out <- base::cut(req(as.numeric(AOO_value)), breaks = c(0, 0.999, 1.999, 2.999, 5.999, 25.999, 125.999, 500.999, 2500.999, 12500.999, 1000000000), labels = c("Z", LETTERS[1:9]))
  
  return(out)

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
  
  eo_factor <- cut(as.numeric(eo_count), breaks = c(0, 0.999, 5.999, 20.001, 80.001, 300.001, 1000000000), labels = c("Z", LETTERS[1:5]))

  out <- list(buffered_occurrences = buffered_occurrences, eo_count = eo_count, factor = eo_factor)
  
  return(out)
}

compare_rank_factors <- function(taxon_data, rank_factor_upload = NULL){
  
  out <- data.frame(
    previous_species_range_letter = NA,
    new_previous_species_range_comparison = NA,
    previous_aoo_letter = NA,
    new_previous_aoo_comparison = NA,
    previous_eocount_letter = NA,
    new_previous_eocount_comparison = NA
  )
  
  if (is.null(rank_factor_upload)){
  
    if (!is.na(taxon_data$species_range_factor)){
      if (!is.null(taxon_data$info_extended$rankInfo$rangeExtent$rangeExtentDescEn)){
        
        previous_species_range_value <- taxon_data$info_extended$rankInfo$rangeExtent$rangeExtentDescEn %>% str_split_1(" square") %>% head(1) %>% str_split_1("-") %>% head(2) %>% parse_number()      
        
        out$previous_species_range_letter <- cut(previous_species_range_value, breaks = c(0, 0.999, 99.999, 249.999, 999.999, 4999.999, 19999.999, 199999.999, 2499999.999, 1000000000), labels = c("Z", LETTERS[1:8])) %>% as.character() %>% unique() %>% paste0(collapse = "")
        
        previous_letters <- strsplit(out$previous_species_range_letter, "")[[1]]
        
        if (length(previous_letters) == 1){
          
          out$new_previous_species_range_comparison <- ifelse(identical(taxon_data$species_range_factor, previous_letters), "equal", 
                                                              ifelse(
                                                                which(LETTERS %in% taxon_data$species_range_factor) < which(LETTERS %in% previous_letters),
                                                                "lower", "higher"
                                                              ))
        } else {
          
          out$new_previous_species_range_comparison <- 
            ifelse(taxon_data$species_range_factor %in% previous_letters, "equal",
                   ifelse(
                     (which(LETTERS %in% taxon_data$species_range_factor) < which(LETTERS %in% previous_letters[1])),
                     "lower",
                     "higher"
                   )
            )
        }
      } 
    } else {
      taxon_data$species_range_factor <- NA
    }

    if (!is.na(taxon_data$AOO_factor)){
    if (!is.null(taxon_data$info_extended$rankInfo$areaOfOccupancy$areaOfOccupancyDescEn)){
      previous_aoo_value <- taxon_data$info_extended$rankInfo$areaOfOccupancy$areaOfOccupancyDescEn %>% str_split_1("-") %>% head(2) %>% parse_number()
      out$previous_aoo_letter <- purrr::map(previous_aoo_value, get_aoo_factor) %>% unlist() %>% as.character() %>% unique() %>% paste0(collapse = "")
      previous_letters <- strsplit(out$previous_aoo_letter, "")[[1]]
      if (length(previous_letters) == 1){
        out$new_previous_aoo_comparison <- ifelse(identical(taxon_data$AOO_factor, previous_letters), "equal", 
                                                            ifelse(
                                                              which(LETTERS %in% taxon_data$AOO_factor) < which(LETTERS %in% previous_letters),
                                                              "lower", "higher"
                                                            ))
      } else {
        out$new_previous_aoo_comparison <- 
          ifelse(taxon_data$AOO_factor %in% previous_letters, "equal",
                 ifelse(
                   (which(LETTERS %in% taxon_data$AOO_factor) < which(LETTERS %in% previous_letters[1])),
                   "lower",
                   "higher"
                 )
          )
      }
    }
    } else {
      taxon_data$species_range_factor <- NA
    }

    if (!is.na(taxon_data$EOcount_factor)){
    if (!is.null(taxon_data$info_extended$rankInfo$numberEos$numberEosDescEn)){
      previous_EOcount_value <- taxon_data$info_extended$rankInfo$numberEos$numberEosDescEn %>% str_split_1("-") %>% head(2) %>% parse_number()
      out$previous_eocount_letter <- cut(as.numeric(previous_EOcount_value), breaks = c(0, 0.999, 5.999, 20.001, 80.001, 300.001, 1000000000), labels = c("Z", LETTERS[1:5])) %>% as.character() %>% unique() %>% paste0(collapse = "")
      previous_letters <- strsplit(out$previous_eocount_letter, "")[[1]]
      if (length(previous_letters) == 1){
        out$new_previous_eocount_comparison <- ifelse(identical(taxon_data$EOcount_factor, previous_letters), "equal", 
                                                  ifelse(
                                                    which(LETTERS %in% taxon_data$EOcount_factor) < which(LETTERS %in% previous_letters),
                                                    "lower", "higher"
                                                  ))
      } else {
        out$new_previous_eocount_comparison <- 
          ifelse(taxon_data$EOcount_factor %in% previous_letters, "equal",
                 ifelse(
                   (which(LETTERS %in% taxon_data$EOcount_factor) < which(LETTERS %in% previous_letters[1])),
                   "lower",
                   "higher"
                 )
          )
      }
    }
  } else {
    out$previous_species_range_letter <- rank_factor_upload$Range.Extent
    out$new_previous_species_range_comparison <- ifelse(identical(taxon_data$species_range_factor, out$previous_species_range_letter), "equal",
                                                        ifelse(
                                                          which(LETTERS %in% taxon_data$species_range_factor) < which(LETTERS %in% out$previous_species_range_letter),
                                                          "lower", "higher"
                                                        ))
    out$previous_aoo_letter <- rank_factor_upload[[names(rank_factor_upload)[grep("Area.of.Occup.", names(rank_factor_upload))][1]]]
    out$new_previous_aoo_comparison <- ifelse(identical(taxon_data$AOO_factor, out$previous_aoo_letter), "equal",
                                              ifelse(
                                                which(LETTERS %in% taxon_data$AOO_factor) < which(LETTERS %in% out$previous_aoo_letter),
                                                "lower", "higher"
                                              ))
    out$previous_eocount_letter <- rank_factor_upload$X..Occur
    out$new_previous_eocount_comparison <- ifelse(identical(taxon_data$EOcount_factor, out$previous_eocount_letter), "equal",
                                                  ifelse(
                                                    which(LETTERS %in% taxon_data$EOcount_factor) < which(LETTERS %in% out$previous_eocount_letter),
                                                    "lower", "higher"
                                                  ))
  }
  } else {
    taxon_data$EOcount_factor <- NA
  }

  return(out)
}
  
# Run full rank assessment
run_rank_assessment <- function(taxon_names, 
                                minimum_fields = c("key", "scientificName", "prov", "longitude", "latitude", "coordinateUncertaintyInMeters", "stateProvince", "countryCode", "year", "month", "datasetName", "institutionCode", "basisOfRecord", "EORANK", "references"),
                                uploaded_data = NULL,
                                rank_factor_upload = NULL,
                                max_number_observations = 10000,
                                clean_occ = TRUE,
                                centroid_filter = FALSE,
                                date_start = "1900-01-01",
                                date_end = "2025-01-01",
                                months = substr(month.name, 1, 3),
                                uncertainty_filter = "",
                                query_polygon = NULL,
                                nations_filter = NULL,
                                states_filter = NULL,
                                network_polys = network_polys,
                                sources_filter = c("OCCURRENCE", "HUMAN_OBSERVATION"),
                                rank_filter = NULL,
                                grid_cell_size = 2,
                                sep_distance = 1000,
                                trends_period1 = c("1900-01-01", "2025-12-31"),
                                trends_period2 = c("1985-01-01", "2025-12-31")
){
  
  # Create a Progress object
  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())

  progress$set(message = "Running Multispecies Assessment", value = 0)
  
  taxon_data_list <- purrr::map(1:length(taxon_names), function(i){
    list(
      info = data.frame(scientificName = "New taxon"),
      info_extended = NULL,
      synonyms = NULL,
      synonyms_selected = NULL,
      datasets = NULL,
      datasets_selected = NULL,
      gbif_occurrences_raw = NULL,
      gbif_occurrences = NULL,
      uploaded_occurrences = NULL,
      drawn_occurrences = NULL,
      all_occurrences = NULL,
      shifted = FALSE,
      sf = NULL,
      sf_filtered = NULL,
      filtered_occurrences = NULL,
      filters_selected = list(
        clean_occ = clean_occ,
        centroid_filter = centroid_filter,
        date_start = date_start,
        date_end = date_end,
        months = months,
        uncertainty_filter = uncertainty_filter,
        nations_filter = nations_filter,
        states_filter = states_filter,
        sources_filter = sources_filter,
        rank_filter = rank_filter,
        grid_cell_size = grid_cell_size,
        sep_distance = sep_distance
      ),
      selected_points = data.frame("Key" = character(), "Scientific name" = character(), "Source" = character(), "Institution code" = character(), "Year" = numeric(), "Coordinate Uncertainty" = numeric(), "Place" = character(), "URL" = character())[NULL, ],
      removed_points = NULL,
      nations = NULL,
      states = NULL,
      records_over_time = NULL,
      species_range_value = NA,
      species_range_map = NA,
      species_range_factor = NA,
      AOO_value = NA,
      AOO_map = NA,
      AOO_factor = NA,
      EOcount_map = NA,
      EOcount_value = NA,
      EOcount_factor = NA,
      rank_factor_comparison = data.frame(
        previous_species_range_letter = NA,
        new_previous_species_range_comparison = NA,
        previous_aoo_letter = NA,
        new_previous_aoo_comparison = NA,
        previous_eocount_letter = NA,
        new_previous_eocount_comparison = NA
      ),
      temporal_change = data.frame(period = c(1, 2), 
                   rec_count = NA, rec_count_change = NA,
                   eoo = NA, eoo_change = NA,
                   aoo = NA, aoo_change = NA,
                   eo_count = NA, eo_count_change = NA
        )
    )
  }) %>% 
    purrr::set_names(taxon_names)
  
  for (i in 1:length(taxon_names)){
    
    taxon_name <- taxon_names[i]

    # Increment the progress bar, and update the detail text.
    progress$inc(0.33/length(taxon_names), detail = paste0("Cross-referencing taxonomy for ", taxon_name))
    
    ns_table_full <- natserv::ns_search_spp(text_adv = list(searchToken = taxon_name, matchAgainst = "allScientificNames", operator="contains"))$results
    gbif_table <- rgbif::name_suggest(q = taxon_name, rank = c("species", "subspecies"), limit = 10)$data

    if (!is.null(ns_table_full) & !is.null(gbif_table)){
      if (nrow(ns_table_full) == 0 & nrow(gbif_table) > 0){
        ns_table_full <- natserv::ns_search_spp(text_adv = list(searchToken = gbif_table$canonicalName[1], matchAgainst = "allScientificNames", operator="contains"))$results
      }
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
    taxon_data_list[[taxon_name]]$info <- taxon_info
    
    if (!is.null(taxon_data_list[[taxon_name]]$info)){
      if (sum(!is.na(taxon_data_list[[taxon_name]]$info$uniqueId)) > 0){
      taxon_data_list[[taxon_name]]$info_extended <- natserv::ns_id(uid = taxon_data_list[[taxon_name]]$info$uniqueId %>% unique() %>% na.omit() %>% as.character() %>% head(1)) 
      taxon_data_list[[taxon_name]]$family <- ns_table_full$speciesGlobal$family
      }
    }
    
    selected_taxon <- c(taxon_info$scientificName, unlist(taxon_info$synonyms)) %>% na.omit() %>% unique()
    
    if (!is.null(uploaded_data)){
      if ("scientificName_Source" %in% names(uploaded_data)){
        selected_taxon <- c(selected_taxon, uploaded_data %>%
                              dplyr::filter(scientificName %in% selected_taxon) %>%
                              dplyr::pull(scientificName_Source) %>%
                              unique()
        ) %>% unique()
      }
    }
    
    selected_taxon <- gsub("ssp. |var. ", "", selected_taxon)
    
    if (length(selected_taxon) == 1){
      taxon_data_list[[taxon_name]]$synonyms <- rgbif::name_usage(name = selected_taxon)$data
    } else {
      taxon_data_list[[taxon_name]]$synonyms <- purrr::map(selected_taxon, function(sp) rgbif::name_usage(name = sp)) %>%
        purrr::map("data") %>%
        bind_rows()
    }

    taxon_data_list[[taxon_name]]$synonyms <- taxon_data_list[[taxon_name]]$synonyms %>% 
      dplyr::distinct(., .keep_all = TRUE)

    if (!is.null(taxon_data_list[[taxon_name]]$synonyms)){
      
    gbif_counts <- purrr::map_dbl(taxon_data_list[[taxon_name]]$synonyms$key,
                                    function(x) rgbif::occ_count(taxonKey = x,
                                                                 gadmGid = query_polygon,
                                                                 hasCoordinate = TRUE
                                                                 )
                                    )

    taxon_data_list[[taxon_name]]$synonyms <- taxon_data_list[[taxon_name]]$synonyms %>% 
      dplyr::mutate(occurrence_count = gbif_counts) %>% 
      dplyr::filter(occurrence_count > 0)

    total_count <- sum(taxon_data_list[[taxon_name]]$synonyms$occurrence_count)
    
    if (("OCCURRENCE" %in% sources_filter) & ("HUMAN_OBSERVATION" %in% sources_filter)){
      taxon_data_list[[taxon_name]]$datasets <- taxon_data_list[[taxon_name]]$datasets_selected <- data.frame(datasetKey = "gbif", count = total_count)
    } else {
      taxon_data_list[[taxon_name]]$datasets <- taxon_data_list[[taxon_name]]$datasets_selected <- data.frame(datasetKey = "detailed", count = total_count)
    }
    }
    
  }
  
  taxa_to_keep <- purrr::map_lgl(1:length(taxon_data_list), function(i){
    keep <- FALSE
    if (nrow(taxon_data_list[[i]]$synonyms) > 0){
      keep <- TRUE
    } else if (!is.null(uploaded_data)){
      keep <- (uploaded_data %>% dplyr::filter(scientificName %in% names(taxon_data_list)[i]) %>% nrow()) > 0
    } 
    keep
    })
  
  taxon_data_list <- taxon_data_list[taxa_to_keep]
  taxon_names <- names(taxon_data_list)
  
  # Increment the progress bar, and update the detail text.
  progress$inc(0.33, detail = paste0("Downloading and processing data..."))

  if ((length(sources_filter) > 0) & (length(taxon_data_list) >= 1)){
    gbif_occurrences_raw <- get_gbif_data(
      taxa_metadata = purrr::map(taxon_data_list, "synonyms") %>% bind_rows(), 
      datasets_metadata = taxon_data_list[[1]]$datasets_selected,
      query_polygon = query_polygon,
      all_occ_data = "OCCURRENCE" %in% sources_filter,
      all_humobs_data = "HUMAN_OBSERVATION" %in% sources_filter, 
      shift_occurrences = FALSE
    )$sp_occurrences
  } else {
    gbif_occurrences_raw <- NULL
  }
  
  taxon_assessment_output <- purrr::map(1:length(taxon_names), function(i){

    taxon_name <- taxon_names[i]
    
    taxon_data <- taxon_data_list[[taxon_name]]
    
    # Increment the progress bar, and update the detail text.
    progress$inc(0.33/length(taxon_names), detail = paste0("Calculating factors for ", taxon_name))
    
    if (!is.null(gbif_occurrences_raw)){
      taxon_data$gbif_occurrences_raw <- gbif_occurrences_raw %>% 
        dplyr::filter(taxonKey %in% taxon_data$synonyms$key)
      taxon_data$gbif_occurrences <- taxon_data$gbif_occurrences_raw %>% 
        clean_gbif_data(clean = clean_occ, remove_centroids = centroid_filter, minimum_fields = minimum_fields) %>% 
        dplyr::mutate(scientificName_Assessment = taxon_name)
      
      shifted <- FALSE
      
      taxon_data$gbif_occurrences$longitude[taxon_data$gbif_occurrences$longitude > 180] <- taxon_data$gbif_occurrences$longitude[taxon_data$gbif_occurrences$longitude > 180] - 360
      max_long <- max(taxon_data$gbif_occurrences$longitude, na.rm = TRUE)/2
      shifted_long <- taxon_data$gbif_occurrences$longitude
      if (length(taxon_data$gbif_occurrences$longitude[taxon_data$gbif_occurrences$longitude > max_long]) > 0){
        shifted_long[shifted_long > max_long] <- shifted_long[shifted_long > max_long] - 360
        shifted_long <- shifted_long + 360
        # shifted <- TRUE
      }

      if ((max(shifted_long)-min(shifted_long)) < (max(taxon_data$gbif_occurrences$longitude) - min(taxon_data$gbif_occurrences$longitude))){
        taxon_data$gbif_occurrences$longitude <- shifted_long
        shifted <- TRUE
      }
      
      taxon_data$shifted <- shifted

    }
    
   if (!is.null(uploaded_data)){

     # if ("scientificName_Assessment" %in% names(taxon_data$uploaded_occurrences)){
     #   if (length(complete.cases(taxon_data$uploaded_occurrences$scientificName_Assessment)) > 0){
     #     taxon_data$uploaded_occurrences <- uploaded_data %>%
     #       dplyr::filter(scientificName_Assessment == taxon_name) %>%
     #       dplyr::select(-scientificName_Source, -scientificName_Assessment)
     #   }
     # } else {
     #   taxon_data$uploaded_occurrences <- uploaded_data %>%
     #     dplyr::filter(scientificName == taxon_name) %>%
     #     dplyr::select(-scientificName_Source, -scientificName_Assessment)
     # }
     
     if ("scientificName_Assessment" %in% names(uploaded_data)){
       if (length(complete.cases(uploaded_data$scientificName_Assessment)) > 0){
       taxon_data$uploaded_occurrences <- uploaded_data %>%
         dplyr::filter(scientificName_Assessment == taxon_name) %>%
         dplyr::select(-scientificName_Source)
       }
     } else {
       taxon_data$uploaded_occurrences <- uploaded_data %>%
         dplyr::filter(scientificName == taxon_name) %>%
         dplyr::select(-scientificName_Source)
     }
     
    }

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
  
  months_numeric <- purrr::map_dbl(1:12, function(x) x) %>% purrr::set_names(substr(month.name, 1, 3))
  season <- which(names(months_numeric) %in% months)
  # season <- ifelse(nchar(season) == 1, paste0("0", season), season)
  
  taxon_data$sf_filtered <- taxon_data$sf_filtered %>%
    dplyr::filter(
      year >= substr(date_start, 1, 4) & year <= substr(date_end, 1, 4) | is.na(year),
      month %in% season | is.na(month)
    )
  
  if (uncertainty_filter != ""){
    
    taxon_data$sf_filtered <- taxon_data$sf_filtered %>%
      dplyr::filter(coordinateUncertaintyInMeters <= as.numeric(uncertainty_filter) | is.na(coordinateUncertaintyInMeters))
  }
  
  # if (!is.null(nations_filter)){
  #   taxon_data$sf_filtered <- taxon_data$sf_filtered %>%
  #     dplyr::filter(purrr::map_int(st_intersects(taxon_data$sf_filtered, network_polys %>% dplyr::filter(FIPS_CNTRY %in% nations_filter)), length) > 0)
  # }
  # 
  # if (!is.null(states_filter)){
  #   taxon_data$sf_filtered <- taxon_data$sf_filtered %>%
  #     dplyr::filter(purrr::map_int(st_intersects(taxon_data$sf_filtered, network_polys %>% dplyr::filter(Admin_abbr %in% states_filter)), length) > 0)
  # }
  
  if (nrow(taxon_data$sf_filtered) >= 1){
   
  # if (!is.null(sources_filter)){
  #   source_exclusions <- setdiff(taxon_data$sf_filtered$prov, sources_filter)
  #   
  #   if (length(source_exclusions) > 0){
  #     taxon_data$sf_filtered <- taxon_data$sf_filtered %>%
  #       dplyr::filter(!(prov %in% source_exclusions))
  #   }
  # }

  eoo_output <- taxon_data$sf_filtered %>% safe_eoo(shifted = shifted)

  if (!is.null(eoo_output$result)){
    eoo_output <- eoo_output$result
    taxon_data$species_range_value <- eoo_output$EOO
    taxon_data$species_range_map <- eoo_output$hull
    taxon_data$species_range_factor <- eoo_output$factor %>% as.character()
  }
  
  taxon_data$filtered_occurrences <- taxon_data$sf_filtered
  taxon_data$filtered_occurrences <- taxon_data$filtered_occurrences %>%
    st_set_geometry(NULL) %>%
    dplyr::select(longitude, latitude) %>%
    as.data.frame()
  taxon_data$filtered_occurrences$longitude[taxon_data$filtered_occurrences$longitude > 180] <- taxon_data$filtered_occurrences$longitude[taxon_data$filtered_occurrences$longitude > 180] - 360
  
  taxon_data$AOO_value <- purrr::safely(aoo2)(taxon_data$filtered_occurrences, as.numeric(grid_cell_size)*1000)
  
  if (!is.null(taxon_data$AOO_value$result)){
    taxon_data$AOO_value <- taxon_data$AOO_value$result/4
    taxon_data$AOO_factor <- purrr::safely(get_aoo_factor)(taxon_data$AOO_value, grid_cell_size = as.numeric(grid_cell_size))
    taxon_data$AOO_factor <- ifelse(!is.null(taxon_data$AOO_factor$result), as.character(taxon_data$AOO_factor$result), NULL)
    taxon_data$AOO_map <- purrr::safely(get_aoo_polys)(taxon_data$sf_filtered, as.numeric(grid_cell_size))      
    if (!is.null(taxon_data$AOO_map$result)){
      taxon_data$AOO_map <- taxon_data$AOO_map$result
    } 
  }

  if (taxon_data$AOO_factor != "H"){
  number_EOs <- purrr::safely(calculate_number_occurrences)(taxon_data$sf_filtered, separation_distance = sep_distance %>% as.numeric(), added_distance = 0)
  } else {
  # thinned_records <- spThin::thin(loc.data = taxon_data$sf_filtered %>% st_set_geometry(NULL) %>% dplyr::mutate(tax = "taxon"), lat.col = "latitude", long.col = "longitude", spec.col = "tax", thin.par = 1, reps = 1, locs.thinned.list.return = TRUE, write.files = FALSE)
  # new_dat <- taxon_data$sf_filtered[thinned_records[[1]] %>% row.names() %>% as.numeric(), ] %>% sample_n(1000)
  number_EOs <- purrr::safely(calculate_number_occurrences)(taxon_data$sf_filtered %>% sample_n(1000), separation_distance = sep_distance %>% as.numeric(), added_distance = 0)
  if (!is.null(number_EOs$result) & number_EOs$result$eo_count <= 300){
    number_EOs <- purrr::safely(calculate_number_occurrences)(taxon_data$sf_filtered %>% sample_n(2000), separation_distance = sep_distance %>% as.numeric(), added_distance = 0)
  }
  }
  if (!is.null(number_EOs$result)){
    number_EOs <- number_EOs$result
    taxon_data$EOcount_value <- number_EOs$eo_count
    taxon_data$EOcount_map <- number_EOs$buffered_occurrences
    taxon_data$EOcount_factor <- number_EOs$factor %>% as.character()
  } 

  if (!is.null(rank_factor_upload)){
    batch_rank_factor_file <- read.csv(rank_factor_upload, header = TRUE)
    if (taxon_name %in% (batch_rank_factor_file[, 3] %>% as.character())){
      batch_rank_factor_sp <- batch_rank_factor_file[which((batch_rank_factor_file[, 3] %>% as.character()) %in% taxon_name), ]
      taxon_data$rank_factor_comparison <- compare_rank_factors(taxon_data, rank_factor_upload = batch_rank_factor_sp)
    } else {
      taxon_data$rank_factor_comparison <- compare_rank_factors(taxon_data)
    }
    } else {
      taxon_data$rank_factor_comparison <- compare_rank_factors(taxon_data)
   }
  
  period1_dat <- taxon_data$sf_filtered %>%
    dplyr::filter(year >= substr(trends_period1[1], 1, 4) & year < substr(trends_period1[2], 1, 4))

  if (nrow(period1_dat) > 0){
    taxon_data$temporal_change$rec_count[1] <- nrow(period1_dat)
    taxon_data$temporal_change$rec_count_change[1] <- 0
    eoo_out <- period1_dat %>% safe_eoo(shifted = shifted)
    if (!is.null(eoo_out$result)){
      eoo_out <- eoo_out$result
      taxon_data$temporal_change$eoo[1] <- eoo_out$EOO
      taxon_data$temporal_change$eoo_change[1] <- 0
    }
    taxon_data$temporal_change$aoo[1] <- (aoo2(period1_dat, as.numeric(grid_cell_size)*1000))/4
    taxon_data$temporal_change$aoo_change[1] <- 0
    if (nrow(taxon_data$sf_filtered) <= 2000){
    taxon_data$temporal_change$eo_count[1] <- (calculate_number_occurrences(period1_dat, separation_distance = as.numeric(sep_distance), added_distance = 0))$eo_count
    taxon_data$temporal_change$eo_count_change[1] <- 0
    }
  }

  period2_dat <- taxon_data$sf_filtered %>%
    dplyr::filter(year >= substr(trends_period2[1], 1, 4) & year < substr(trends_period2[2], 1, 4))

  if (nrow(period2_dat) > 0){
    taxon_data$temporal_change$rec_count[2] <- nrow(period2_dat)
    taxon_data$temporal_change$rec_count_change[2] <- round(((taxon_data$temporal_change$rec_count[2]-taxon_data$temporal_change$rec_count[1])/taxon_data$temporal_change$rec_count[1])*100, 1)
    eoo_out2 <- period2_dat %>% safe_eoo(shifted = shifted)
    if (!is.null(eoo_out2$result)){
      eoo_out2 <- eoo_out2$result
      taxon_data$temporal_change$eoo[2] <- eoo_out2$EOO
      if (!is.na(taxon_data$temporal_change$eoo[1])){
      taxon_data$temporal_change$eoo_change[2] <- round(((taxon_data$temporal_change$eoo[2]-taxon_data$temporal_change$eoo[1])/taxon_data$temporal_change$eoo[1])*100, 1)
      }
    }
    taxon_data$temporal_change$aoo[2] <- (aoo2(period2_dat, as.numeric(grid_cell_size)*1000))/4
    taxon_data$temporal_change$aoo_change[2] <- round(((taxon_data$temporal_change$aoo[2]-taxon_data$temporal_change$aoo[1])/taxon_data$temporal_change$aoo[1])*100, 1)
    if (nrow(taxon_data$sf_filtered) <= 2000){
    taxon_data$temporal_change$eo_count[2] <- (calculate_number_occurrences(period2_dat, separation_distance = sep_distance %>% as.numeric(), added_distance = 0))$eo_count
    taxon_data$temporal_change$eo_count_change[2] <- round(((taxon_data$temporal_change$eo_count[2]-taxon_data$temporal_change$eo_count[1])/taxon_data$temporal_change$eo_count[1])*100, 1)
    }
  }
  }  

  taxon_data
  
  }) %>% 
    purrr::set_names(taxon_names)
  
  return(taxon_assessment_output)
  
}

# Run rank assessment safely without errors interrupting it
safe_batch_run <- purrr::safely(run_rank_assessment)

# Function to extract id from table button click
shinyInput <- function(FUN, len, id, ...) {
  inputs <- character(len)
  for (i in seq_len(len)) {
    inputs[i] <- as.character(FUN(paste0(id, i), ...))
  }
  inputs
}

# Calculate rarity change
calculate_rarity_change <- function(taxon_data = taxon_data, 
                                    period1 = c("1900-01-01", "2025-12-31"),
                                    period2 = c("1985-01-01", "2025-12-31"),
                                    period3 = c(NA, NA),
                                    aoo_grid_cell_size = 2,
                                    occ_sep_distance = 1000
                                    ){
  
  period_list <- list(period1, period2, period3)
  
  # Remove empty periods
  period_list <- period_list[purrr::map_lgl(period_list, function(x) !identical(x, c(NA, NA)))]
  
  number_periods <- length(period_list)
  
  dat <- data.frame(period = 1:number_periods, 
                    rec_count = NA, rec_count_change = NA,
                    eoo = NA, eoo_change = NA,
                    aoo = NA, aoo_change = NA,
                    eo_count = NA, eo_count_change = NA
  )
  
  for (period in 1:length(period_list)){
    
    period_dat <- taxon_data$sf_filtered %>% 
      dplyr::filter(year >= substr(period_list[[period]][1], 1, 4) & year < substr(period_list[[period]][2], 1, 4))
    
    if (nrow(period_dat) > 0){
      dat$rec_count[period] <- nrow(period_dat)
      dat$rec_count_change[period] <- ifelse(period > 1, round(((dat$rec_count[period]-dat$rec_count[period-1])/dat$rec_count[period-1])*100, 1), 0)
      dat$eoo[period] <- (period_dat %>% calculate_eoo(shifted = taxon_data$shifted))$EOO
      dat$eoo_change[period] <- ifelse(period > 1, round(((dat$eoo[period]-dat$eoo[period-1])/dat$eoo[period-1])*100, 1), 0)
      dat$aoo[period] <- (aoo2(period_dat, as.numeric(aoo_grid_cell_size)*1000))/4
      dat$aoo_change[period] <- ifelse(period > 1, round(((dat$aoo[period]-dat$aoo[period-1])/dat$aoo[period-1])*100, 1), 0)
      dat$eo_count[period] <- (calculate_number_occurrences(period_dat, separation_distance = occ_sep_distance %>% as.numeric(), added_distance = 0))$eo_count
      dat$eo_count_change[period] <- ifelse(period > 1, round(((dat$eo_count[period]-dat$eo_count[period-1])/dat$eo_count[period-1])*100, 1), 0)
    }

  }
  
  return(dat)

}

# Get temporal trends
get_temporal_trends <- function(taxon_data = taxon_data, referenceTaxon = "kingdom", start_year = 1980){
  
  # Create a Progress object
  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())
  
  progress$set(message = "Running Temporal Analysis", value = 0)
  
  min_fields <- c("key", "scientificName", "prov", "longitude", "latitude", "coordinateUncertaintyInMeters", "stateProvince", "countryCode", "year", "month", "datasetName", "institutionCode", "basisOfRecord", "EORANK", "references")

  key_value <- taxon_data$synonyms %>% 
    dplyr::filter(datasetKey == "d7dddbf4-2cf0-4f39-9b2a-bb099caae36c") %>% 
    dplyr::select(paste0(referenceTaxon, "Key")) %>% 
    as.character()
  
  if (is.null(taxon_data$AOO_map)){
    taxon_data$AOO_map <- get_aoo_polys(taxon_data$sf_filtered, 2)
  }
  
  if (nrow(taxon_data$AOO_map) >= 100){
    taxon_data$AOO_map <- taxon_data$AOO_map[which(st_intersects(taxon_data$sf_filtered, taxon_data$AOO_map) %>% unlist() %>% table() > 1), ]
  }
  
  query_poly <- taxon_data$AOO_map %>% 
    # sf::st_make_valid() %>%
    # sf::st_union() %>%
    # sf::st_make_valid() %>%
    terra::vect() %>%
    terra::forceCCW() %>%
    terra::geom(wkt = TRUE)
  
  # Increment the progress bar, and update the detail text.
  progress$inc(0.33/1, detail = paste0("Downloading GBIF data for reference taxon across target area"))
  
  gbif_data <- spocc::occ(from = "gbif", gbifopts = list(
    taxonKey = key_value,
    geometry = query_poly,
    year = paste0(start_year, ",", substr(Sys.Date(), 1, 4))
  ),
  limit = 100000, 
  has_coords = TRUE
  )
  
  if (nrow(gbif_data$gbif$data$custom_query) == 0){
    
    query_poly <- taxon_data$AOO_map %>% 
      # sf::st_make_valid() %>%
      # sf::st_union() %>%
      # sf::st_make_valid() %>%
      terra::vect() %>%
      terra::forceCCW() %>%
      terra::geom(wkt = TRUE)
    
    gbif_data <- spocc::occ(from = "gbif", gbifopts = list(
      taxonKey = key_value,
      geometry = query_poly,
      year = paste0(start_year, ",", substr(Sys.Date(), 1, 4))
    ),
    limit = 100000, 
    has_coords = TRUE
    )
    
  }
  
  # Increment the progress bar, and update the detail text.
  progress$inc(0.33/1, detail = paste0("Preparing GBIF data for temporal analysis"))
  
  gbif_data <- gbif_data$gbif$data$custom_query
  
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
            dplyr::select(names(gbif_data)) %>% 
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
  
  # Increment the progress bar, and update the detail text.
  progress$inc(0.33/1, detail = paste0("Creating yearly summaries"))
  
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
  
  annotation_x <- quantile(temporal_trend_data$year, .35)
  num_chars <- list(nchar(paste0("Observations of ", taxon_data$info$scientificName)),
                    nchar(paste0("Observations of ", referenceTaxon, " ", reference_taxon_name)),
                    nchar("Proportion of Observations"),
                    nchar("Modeled probability of detection")
                    )
  num_chars_missing <- purrr::map(num_chars, function(n) max(unlist(num_chars)) - n)

  # Number of focal observations
  p1 <- ggplot(data = temporal_trend_data, aes(x = year, y = focal_species_count)) + 
    geom_col(col = "#2c7bb680", alpha = 0.5) + 
    ylab(paste0("Records of \n", taxon_data$info$scientificName)) +
    xlab("") +
    # annotate("text", x = quantile(temporal_trend_data$year, .26), y = max(temporal_trend_data$focal_species_count)+(0.08 *max(temporal_trend_data$focal_species_count)), label = paste0("Observations of ", taxon_data$info$scientificName, rep(" ", num_chars_missing[[1]])), hjust = 1) +
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
    ylab(paste0("Records of \n ", referenceTaxon, " ", reference_taxon_name)) +
    xlab("") +
    # annotate("text", x = quantile(temporal_trend_data$year, .24), y = max(temporal_trend_data$total_number_observations)+(0.08 *max(temporal_trend_data$total_number_observations)), label = paste0("Observations of ", referenceTaxon, " ", reference_taxon_name, rep(" ", num_chars_missing[[2]])), hjust = 1) +
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
    geom_smooth(size = 1, se = FALSE, col = "black") +
    ylab(paste0("Proportion of \n records")) +
    xlab("Year") +
    # annotate("text", x = quantile(temporal_trend_data$year, .21), y = max(temporal_trend_data$proportion_observations)+(0.08 *max(temporal_trend_data$proportion_observations)), label = paste0("Proportion of Observations"), hjust = 1) +
    theme_linedraw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 8)
    )
  # # Mean proportion of species
  # p4 <- ggplot(data = temporal_trend_data, aes(x = year, y = proportion_species)) +
  #   # geom_col(col = "#2c7bb680", alpha = 0.5) +
  #   geom_smooth(size = 1.5, se = FALSE, col = "black") +
  #   ylab("") +
  #   xlab("") +
  #   # annotate("text", x = median(temporal_trend_data$year), y = max(temporal_trend_data$proportion_species)+(0.08 *max(temporal_trend_data$proportion_species)), label = "Proportion of Species") +
  #   # ylim(c(0, max(temporal_trend_data$proportion_species)+(0.2 *max(temporal_trend_data$proportion_species)))) +
  #   theme_linedraw() +
  #   theme(legend.position = "none",
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         axis.title = element_text(size = 9),
  #         axis.text.x = element_blank(),
  #         axis.text = element_text(size = 8)
  #   )
  # 
  # # Yearly reporting rate
  # p5 <- ggplot(data = temporal_trend_data, aes(x = year, y = reporting_rate)) + 
  #   geom_col(col = "#2c7bb680", alpha = 0.5) + 
  #   geom_smooth(size = 1, se = FALSE, col = "black") + 
  #   ylim(c(0, max(temporal_trend_data$reporting_rate)+(0.2 *max(temporal_trend_data$reporting_rate)))) +
  #   ylab("") +
  #   annotate("text", x = median(temporal_trend_data$year), y = max(temporal_trend_data$reporting_rate)+(0.08 *max(temporal_trend_data$reporting_rate)), label = "Observation Rate") +
  #   theme_linedraw() +
  #   theme(legend.position = "none",
  #         panel.grid.major = element_blank(), 
  #         panel.grid.minor = element_blank(),
  #         axis.title = element_text(size = 9),
  #         axis.text.x = element_blank(),
  #         axis.text = element_text(size = 8)
  #   )
  ## GLMER (with spatial random effect)
  # options(na.action = "na.fail")
  # glmer_result <- lme4::glmer(focal_species_detection ~ splines::bs(year, df = 3) + scale(total_number_observations) + scale(species_list_length) + (1 | cellID), 
  #                             data = focal_species_detection_data, family = binomial, control = lme4::glmerControl(tol = 1e-5, optimizer = "bobyqa", optCtrl=list(maxfun=2e5))
  # )
  # # Modeled probability of detection
  # p4 <- data.frame(fitted = fitted(glmer_result), year = focal_species_detection_data$year) %>% 
  #   ggplot(aes(x = year, y = fitted)) + 
  #   geom_point(col = "#2c7bb680", alpha = 0.5) + 
  #   geom_smooth(size = 1, se = FALSE, col = "black") + 
  #   # ylim(c(0, 1.2)) +
  #   scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1)) +
  #   xlab("Year") +
  #   ylab(paste0("Modeled probability \n of detection")) +
  #   # annotate("text", x = quantile(temporal_trend_data$year, .20), y = 1.08, label = paste0("Modeled probability of detection", rep(" ", num_chars_missing[[4]])), hjust = 1) +
  #   theme_linedraw() +
  #   theme(legend.position = "none",
  #         panel.grid.major = element_blank(), 
  #         panel.grid.minor = element_blank(),
  #         axis.title = element_text(size = 9),
  #         axis.text = element_text(size = 8)
  #   )
  # 
  # p <- subplot(p1, p2, p3, p4, nrows = 4, shareX = TRUE, titleX = TRUE, titleY = TRUE, margin = 0.01, which_layout = 1)

  p <- subplot(p1, p2, p3, nrows = 3, shareX = TRUE, titleX = TRUE, titleY = TRUE, margin = 0.01, which_layout = 1)
  
  gg <- plotly_build(p) %>%
    config(displayModeBar = FALSE) %>%
    layout(# plot_bgcolor  = "rgba(0, 0, 0, 0)",
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