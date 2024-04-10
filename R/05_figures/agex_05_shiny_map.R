## ------------------------------------------ ##
## Interactive map extreme clusters one season
## By: Harold Achicanoy
## WUR & ABC
## Mar. 2024
## ------------------------------------------ ##

options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
pacman::p_load(shiny,terra,leaflet,htmlwidgets,leaflet.extras)

agex_sgn <- terra::vect('https://raw.githubusercontent.com/haachicanoy/agroclimExtremes/main/data/agex_spei-6_one_s1_features.gpkg')
agex_sgn <- sf::st_as_sf(agex_sgn)
col_pltt <- MetBrewer::met.brewer(name = 'Signac', n = nrow(agex_sgn))
set.seed(1235)
col_pltt <- sample(x = col_pltt, size = length(col_pltt), replace = F)

ui <- shiny::bootstrapPage(
  tags$style(type='text/css','html, body {width:100%;height:100%}'),
  leaflet::leafletOutput('mymap', width='100%', height='100%')
)

server <- function(input, output, session) {
  
  rv <- shiny::reactiveValues()
  rv$selected <- NULL
  
  output$mymap <- leaflet::renderLeaflet({
    hl_opts <- leaflet::highlightOptions(fillColor = 'yellow', bringToFront = T)
    leaflet() |>
      addProviderTiles(providers$CartoDB.Positron) |>
      addResetMapButton() |>
      addSearchOSM() |>
      addPolygons(
        layerId   = ~ extreme_cluster,
        group     = 'cluster',
        data      = agex_sgn,
        color     = col_pltt,
        opacity   = 0.8,
        fill      = T,
        weight    = 2,
        label     = ~ extreme_cluster,
        popup     = ~ paste("<div class='leaflet-popup-scrolled' style='max-width:250px;max-height:350px'",
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
                            '<b>','SPEI-6 90th pr trend:  ','</b>',`SPEI-6_slope_90th`,'<br>',
                            '<b>','Cluster cohesion: ','</b>',cls_cohesion,'<br>',
                            '<b>','Cluster contiguity:  ','</b>',cls_contiguity,'<br>',
                            '<b>','Cluster aggregation:  ','</b>',cls_aggregation,'<br>',
                            '<b>','Quality rank:  ','</b>',quality_rank,'<br>'),
        highlightOptions = hl_opts)
  })
  
  output$click_on_shape <- shiny::renderPrint({ input$mymap_shape_click })
  
  shiny::observeEvent(c(input$mymap_shape_click, input$mymap_click), {
    new_selected <- shiny::req(input$mymap_shape_click)
    shiny::isolate(old_selected <- rv$selected)
    if(is.null(old_selected) || new_selected$.nonce != old_selected$.nonce){
      shiny::validate( shiny::need(new_selected$group!='selection', message=F) )
      rv$selected <- new_selected
      i <- which(agex_sgn$extreme_cluster==new_selected$id) 
      agex_sgn_filtered <- agex_sgn[i,]
      leaflet::leafletProxy('mymap') |>
        leaflet::clearGroup('selection') |>
        addPolygons(
          layerId     = ~ extreme_cluster,
          group       = 'selection',
          data        = agex_sgn_filtered,
          fillColor   = ~'cyan',
          weight      = 2,
          opacity     = 0.4,
          fillOpacity = 0.8)
    } else {
      rv$selected <- NULL
      leaflet::leafletProxy('mymap') |>
        leaflet::clearGroup('selection')
    }
  })
  
}

shiny::shinyApp(ui, server)