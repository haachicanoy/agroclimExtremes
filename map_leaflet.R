library(pacman)
pacman::p_load(terra,tidyverse,leaflet,shiny,htmlwidgets,leaflet.extras)

agex_sgn_poly <- terra::vect('D:/agex_features.gpkg')
my.palette    <- MetBrewer::met.brewer(name = 'Signac', n = nrow(agex_sgn_poly))
set.seed(1235)
my.palette    <- sample(x = my.palette, size = length(my.palette), replace = F)

int_map <- leaflet() |>
  addProviderTiles(providers$CartoDB.Positron) |>
  # addProviderTiles(providers$Esri.WorldImagery, group = "Aerial") |> # aerial map
  addResetMapButton() |> # reset map extent on click
  addSearchOSM() |> # search for towns/coords
  addPolygons(data    = agex_sgn_poly,
              color   = my.palette,
              opacity = 0.8,
              fill    = T,
              weight  = 2,
              label   = ~ extreme_signature,
              popup   = ~ paste("<div class='leaflet-popup-scrolled' style='max-width:250px;max-height:250px'", 
                                '<br>',
                                '<b>', 'Extreme signature: ', '</b>', extreme_signature, "<br>",
                                '<b>', 'Crops diversity:  ', '</b>', crop_types_diversity, "<br>",
                                '<b>', 'Value of production:  ', '</b>', value_of_production, "<br>",
                                '<b>', 'Population density:  ', '</b>', population_density, "<br>",
                                '<b>', 'Signature cohesion: ', '</b>', sgn_cohesion, "<br>",
                                '<b>', 'Signature contiguity:  ', '</b>', sgn_contiguity, "<br>"),
              highlightOptions = highlightOptions(fillColor = 'yellow', bringToFront = T))
int_map
