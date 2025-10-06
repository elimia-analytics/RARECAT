#' Please do not rename this file !
## Identify required package names
### CRAN packages
cran_packages <- c("shiny", "leaflet", "leaflet.extras", "tidyverse", "shinyjs", "sf", "terra", 
                   "plotly", "htmltools", "htmlwidgets", "shinyWidgets", "sortable", 
                   "shinydashboard", "shinycssloaders", "shinyBS", "DT", "natserv", "rgbif", "flexdashboard", 
                   "shinybusy", "leafpm", "red", "dygraphs", "RWmisc", "units", "shiny.exe", 
                   "writexl", "spocc")
installed_cran_packages <- cran_packages %in% rownames(installed.packages())
if (any(installed_cran_packages == FALSE)) {
  install.packages(cran_packages[!installed_cran_packages])
}
library(shiny.exe)
hostUnix(
    appDir = '/Users/giorap/Library/CloudStorage/Dropbox/elimia/code/apps/RARECAT',
    port = getOption('shiny.port'),
    launch.browser = TRUE,
    host = '127.0.0.1',
                  workerId = '',
                  quiet = FALSE,
                  display.mode = c('auto', 'normal', 'showcase'),
                  test.mode = getOption('shiny.testmode', FALSE))
