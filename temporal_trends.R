library(gbifdb)
library(tidyverse)

gbif_data <- gbifdb::gbif_remote(backend = "duckdb")
states <- unique(x$`Juncus abortivus`$sf_filtered$stateProvince)
gbif_data_subset <- gbif_data %>%
  filter(stateprovince %in% states, 
         family == x$`Juncus abortivus`$family
         ) 
system.time({
  gbif_data <- gbifdb::gbif_remote(backend = "duckdb")
  states <- unique(x$`Juncus abortivus`$sf_filtered$stateProvince)
  gbif_data_subset <- gbif_data %>%
    filter(stateprovince %in% states, 
           family == x$`Juncus abortivus`$family
    ) 
  control_count <- gbif_data_subset %>% 
  dplyr::count(year)
  control_count <- control_count %>% as.data.frame()
})

x <- gbif %>%
  filter(phylum == "Chordata", year > 1990) %>%
  count(class, year)


## Extract species occurrences across all relevant scientific names from GBIF using the SPOCC package
gbif_data <- spocc::occ(from = "gbif", gbifopts = list(kingdomKey = "6", geometry = x$`Juncus abortivus`$AOO_map %>% wk::as_wkt()),
                        limit = 1000000, 
                        has_coords = TRUE
)
