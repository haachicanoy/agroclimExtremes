## ------------------------------------------ ##
## Paper figures/maps
## By: Harold Achicanoy
## WUR & ABC
## Dec. 2023
## ------------------------------------------ ##

# R options and packages loading
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,lubridate,tidyverse))

## Key arguments
root   <- '//CATALOGUE/WFP_ClimateRiskPr1'
index  <- 'spei-6'
gs     <- 'one'
season <- 1

# ## Graph of one index-year
# wrl <- rnaturalearth::ne_countries(scale = 'large', returnclass = 'sv')
# wrl <- terra::aggregate(wrl)
# # my.palette <- RColorBrewer::brewer.pal(n = 20, name = 'Set1')
# my.palette <- MetBrewer::met.brewer(name = 'Tam', n = 20)
# tmp <- cdd[[terra::nlyr(cdd)]]
# tmp <- terra::trim(x = tmp)
# hist(terra::values(tmp))
# qnt <- quantile(x = terra::values(tmp), probs = 0.98, na.rm = T)
# terra::values(tmp)[terra::values(tmp) > qnt] <- qnt
# png(filename = "D:/test.png", width = 3132, height = 2359, units = 'px')
# plot(wrl, ext = terra::ext(tmp))
# plot(tmp, add = T, col = my.palette, plg = list(cex = 5), cex.main = 7)
# dev.off()

## Graph of global cluster One growing season
# World shapefile
wrl <- rnaturalearth::ne_countries(scale = 'large', returnclass = 'sv')
wrl <- terra::aggregate(wrl)
# Load extreme weather signatures
fls <- list.files(path = paste0(root,'/agroclimExtremes/agex_results/clusters'), pattern = paste0('agex_global_',index,'_',gs,'_s',season,'_fmadogram_*.*.tif$'), full.names = T)
agex_sgn <- terra::rast(fls)
names(agex_sgn) <- 'extreme_signature'
# Colour palette
n_col <- length(unique(terra::values(agex_sgn,na.rm = T)))
col_pltt <- MetBrewer::met.brewer(name = 'Signac', n = n_col)
set.seed(1); col_pltt <- sample(x = col_pltt, size = n_col, replace = F); rm(n_col)
agex_sgn <- terra::trim(x = agex_sgn)
# Figure
png(filename = paste0('D:/agex_global_',index,'_',gs,'_s',season,'_fmadogram.png'), width = 3132, height = 2359, units = 'px')
plot(wrl, ext = terra::ext(agex_sgn))
plot(agex_sgn, add = T, col = col_pltt, plg = list(cex = 5), cex.main = 7)
dev.off()

## Time series figure
# Load extreme weather signatures
fls <- list.files(path = paste0(root,'/agroclimExtremes/agex_results/clusters'), pattern = paste0('agex_global_',index,'_',gs,'_s',season,'_fmadogram_*.*.tif$'), full.names = T)
agex_sgn <- terra::rast(fls)
names(agex_sgn) <- 'extreme_signature'
# Load indices values
idx <- terra::rast(paste0(root,'/agroclimExtremes/agex_indices/agex_',index,'/agex_',index,'_25km/',gs,'_s',season,'_',index,'_25km.tif'))
idx_tsr <- terra::zonal(x = idx, z = agex_sgn, fun = 'median')

idx_tsr_lng <- idx_tsr |>
  tidyr::pivot_longer(cols = 2:ncol(idx_tsr), names_to = 'year', values_to = 'spei') |>
  base::as.data.frame()
idx_tsr_lng$year <- gsub('spei-6_','',idx_tsr_lng$year) |> as.numeric()
idx_tsr_lng$extreme_signature <- as.factor(idx_tsr_lng$extreme_signature)

tsr <- idx_tsr_lng |>
  ggplot2::ggplot(aes(x = year, y = spei, colour = extreme_signature)) +
  ggplot2::geom_line(alpha = 0.4) +
  ggplot2::scale_color_manual(values = col_pltt) +
  ggplot2::xlab('Year') +
  ggplot2::ylab('SPEI-6 median across signature') +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = 'blue') +
  ggplot2::theme_bw() +
  ggplot2::theme(text            = element_text(size = 17, colour = 'black'),
                 axis.text       = element_text(size = 16, colour = 'black'),
                 axis.title      = element_text(size = 20, colour = 'black'),
                 legend.text     = element_text(size = 13, colour = 'black'),
                 legend.title    = element_blank(),
                 plot.title      = element_text(size = 25, colour = 'black'),
                 plot.subtitle   = element_text(size = 17, colour = 'black'),
                 strip.text.x    = element_text(size = 17, colour = 'black'),
                 strip.text.y    = element_text(size = 17, colour = 'black'),
                 plot.caption    = element_text(size = 15, hjust = 0, colour = 'black'),
                 legend.position = 'none')
ggplot2::ggsave(filename = 'D:/time_series_per_signature.png', plot = tsr, device = 'png', units = 'in', width = 12, height = 8, dpi = 350)
