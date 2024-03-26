options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
pacman::p_load(shiny,terra,leaflet,htmlwidgets,leaflet.extras)

agex_sgn <- terra::vect('https://raw.githubusercontent.com/haachicanoy/agroclimExtremes/main/data/agex_features.gpkg')
agex_sgn <- sf::st_as_sf(agex_sgn)
col_pltt <- MetBrewer::met.brewer(name = 'Signac', n = nrow(agex_sgn))
set.seed(1235)
col_pltt <- sample(x = col_pltt, size = length(col_pltt), replace = F)

hl_opts <- leaflet::highlightOptions(fillColor = 'yellow', bringToFront = T)
int_map <- leaflet() |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addResetMapButton() |>
  addSearchOSM() |>
  addPolygons(data    = agex_sgn,
              color   = col_pltt,
              opacity = 0.8,
              fill    = T,
              weight  = 2,
              label   = ~ extreme_cluster,
              popup   = ~ paste("<div class='leaflet-popup-scrolled' style='max-width:250px;max-height:350px'",
                                '<br>',
                                '<b>','Extreme cluster: ','</b>',extreme_cluster,'<br>',
                                '<b>','Countries count:  ','</b>',countries_count,'<br>',
                                '<b>','Countries ISO3:  ','</b>',countries,'<br>',
                                '<b>','Continents:  ','</b>',continents,'<br>',
                                '<b>','Population density:  ','</b>',population_density,'<br>',
                                '<b>','Value of production:  ','</b>',value_of_production,'<br>',
                                '<b>','Crop types diversity:  ','</b>',crop_types_diversity,'<br>',
                                '<b>','Livestock units:  ','</b>',livestock_units,'<br>',
                                '<b>','SPEI-6 mean:  ','</b>',`SPEI-6_mean`,'<br>',
                                '<b>','SPEI-6 median:  ','</b>',`SPEI-6_median`,'<br>',
                                '<b>','SPEI-6 max:  ','</b>',`SPEI-6_max`,'<br>',
                                '<b>','SPEI-6 general trend:  ','</b>',`SPEI-6_slope`,'<br>',
                                '<b>','SPEI-6 95th pr trend:  ','</b>',`SPEI-6_slope_95th`,'<br>',
                                '<b>','Cluster cohesion: ','</b>',cls_cohesion,'<br>',
                                '<b>','Cluster contiguity:  ','</b>',cls_contiguity,'<br>',
                                '<b>','Cluster aggregation:  ','</b>',cls_aggregation,'<br>',
                                '<b>','Quality rank:  ','</b>',quality_rank,'<br>'),
              highlightOptions = hl_opts)
int_map
