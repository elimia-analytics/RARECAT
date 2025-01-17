get_temporal_trend <- function(referenceTaxon = "kingdom", start_year = 1980){
  
  key_value <- taxon_data$synonyms %>% 
    dplyr::filter(datasetKey == "d7dddbf4-2cf0-4f39-9b2a-bb099caae36c") %>% 
    dplyr::select(paste0(referenceTaxon, "Key")) %>% 
    as.character()
  
  gbif_data <- spocc::occ(from = "gbif", gbifopts = list(
    taxonKey = key_value,
    geometry = taxon_data$AOO_map %>% wk::as_wkt(),
    year = paste0(start_year, ",", substr(Sys.Date(), 1, 4))
  ),
  limit = 1000000, 
  has_coords = TRUE
  )$gbif$data$custom_query
  
  gbif_data <- gbif_data %>% 
    dplyr::filter(
      complete.cases(longitude, latitude),
      !(scientificName %in% unique(taxon_data$sf_filtered$scientificName))
    ) %>% 
    dplyr::select(intersect(names(gbif_data), minimum_fields)) %>% 
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
      proportion_observations = mean(proportion_observations),
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
    ylab("Observations of Target Taxon") +
    xlab("Year") +
    theme_linedraw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 9),
          axis.text.x = element_blank(),
          axis.text = element_text(size = 8),
    )
  # Total number of observations of reference taxon
  p2 <- ggplot(data = temporal_trend_data, aes(x = year, y = total_number_observations)) + 
    geom_col(col = "#2c7bb680", alpha = 0.5) + 
    ylab("Observations of Reference Taxon") +
    xlab("Year") +
    theme_linedraw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 9),
          axis.text.x = element_blank(),
          axis.text = element_text(size = 8),
    )
  # Mean Observation Rate
  p3 <- ggplot(data = temporal_trend_data, aes(x = year, y = proportion_observations)) +
    geom_col(col = "#2c7bb680", alpha = 0.5) +
    geom_smooth(size = 1.5, se = FALSE, col = "black") +
    ylab("Proportion of Observations") +
    xlab("Year") +
    ylim(c(0, max(mean_yearly_proportions$proportion_species)+0.05)) +
    theme_linedraw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 9),
          axis.text.x = element_blank(),
          axis.text = element_text(size = 8),
    )
  # Mean proportion of species
  p4 <- ggplot(data = temporal_trend_data, aes(x = year, y = proportion_species)) +
    geom_col(col = "#2c7bb680", alpha = 0.5) +
    geom_smooth(size = 1.5, se = FALSE, col = "black") +
    ylab("Proportion of Species") +
    xlab("Year") +
    ylim(c(0, max(mean_yearly_proportions$proportion_species)+0.05)) +
    theme_linedraw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 9),
          axis.text.x = element_blank(),
          axis.text = element_text(size = 8),
    )
  
  # Yearly reporting rate
  p5 <- ggplot(data = temporal_trend_data, aes(x = year, y = reporting_rate)) + 
    geom_col(col = "#2c7bb680", alpha = 0.5) + 
    geom_smooth(size = 1, se = FALSE, col = "black") + 
    ylim(c(0, max(reporting_rate_byYear$reporting_rate)+0.05)) +
    ylab("Observation Rate") +
    xlab("Year") +
    theme_linedraw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 9),
          axis.text.x = element_blank(),
          axis.text = element_text(size = 8),
    )
  ## GLMER (with spatial random effect)
  options(na.action = "na.fail")
  glmer_result <- lme4::glmer(focal_species_detection ~ splines::bs(year, df = 3) + total_number_observations + species_list_length + (1 | cellID), 
                              data = focal_species_detection_data, family = binomial, control = glmerControl(tol = 1e-5, optimizer = "bobyqa", optCtrl=list(maxfun=2e5))
  )
  # Modeled probability of detection
  p6 <- data.frame(fitted = fitted(glmer_result), year = focal_species_detection_data$year) %>% 
    ggplot(aes(x = year, y = fitted)) + 
    geom_point(col = "#2c7bb680", alpha = 0.5) + 
    geom_smooth(size = 1, se = FALSE, col = "black") + 
    ylim(c(0, max(fitted(glmer_result))+0.05)) +
    ylab("Modeled probability of detection") +
    xlab("Year") +
    theme_linedraw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 8),
    )
  
  p <- subplot(p1, p2, p3, p5, p6, nrows = 5, shareX = TRUE, titleX = TRUE, titleY = TRUE)
  
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