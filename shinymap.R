library(tidyverse)
library(shiny)
library(sf)
library(leaflet)
library(rnaturalearth)
library(rnaturalearthdata)

#data(africa)
world <- ne_countries(scale = "medium", returnclass = "sf")

africa <- world %>% filter(continent == "Africa") %>%
  # add rownumber
  dplyr::mutate(id = row_number()) %>%
  sf::st_transform(4326)


bbox <- st_bbox(africa$geometry) %>% 
  as.vector()

ui <- bootstrapPage(
  tags$style(type = "text/css", "html, body {width:100%;height:100%}"),
  leafletOutput("mymap", width = "100%", height = "100%")
)

server <- function(input, output, session) {
  
  rv <- reactiveValues()
  rv$selected <- NULL
  
  output$mymap <- renderLeaflet({
    hl_opts <- highlightOptions(
      color = "#CC0000", weight = 3, bringToFront = TRUE)
    leaflet() %>% addTiles() %>% 
      addPolygons(
        layerId = ~id,
        group = "countries",
        data = africa,
        label = ~name,
        fillColor = "#FCFFA4",
        weight = 1,
        color = "#666666",
        opacity = 0.4,
        fillOpacity = 0.8,
        highlightOptions = hl_opts)
  })
  
  output$click_on_shape <- renderPrint({
    input$mymap_shape_click
  })
  
  observeEvent(input$mymap_click, {
    new_selected <- req(input$mymap_shape_click)
    isolate(old_selected <- rv$selected)
    if (is.null(old_selected) || new_selected$.nonce != old_selected$.nonce) {
      validate(
        need(new_selected$group!="selection", message=FALSE)
      )
      rv$selected <- new_selected
      i <- which(africa$id==new_selected$id) 
      africa_filtered <- africa[i,]
      leafletProxy("mymap") %>%
        clearGroup("selection") %>%
        addPolygons(
          layerId = ~id,
          group = "selection",
          data = africa_filtered,
          fillColor = "cyan",
          weight = 1.2,
          color = "#666666",
          opacity = 0.4,
          fillOpacity = 0.8)
    } else {
      rv$selected <- NULL
      leafletProxy("mymap") %>%
        clearGroup("selection")
    }
  })
  
}

shinyApp(ui, server)