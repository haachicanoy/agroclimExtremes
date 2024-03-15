# ------------------------------------------ #
# Description of extreme weather signatures
# By: Harold Achicanoy
# WUR & ABC
# Dec. 2023
# ------------------------------------------ #

# R options and packages loading
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,geodata,NbClust,entropy,scales,dplyr,
                                ggplot2,RColorBrewer,landscapemetrics))

## Key arguments
root   <- '//CATALOGUE/WFP_ClimateRiskPr1'
index  <- 'spei-6'
gs     <- 'one'
season <- 1

# Load extreme weather signatures
fls <- list.files(path = paste0(root,'/agroclimExtremes/agex_results/clusters'), pattern = paste0('agex_global_',index,'_',gs,'_s',season,'_fmadogram_*.*.tif$'), full.names = T)
agex_sgn <- terra::rast(fls)
names(agex_sgn) <- 'extreme_signature'
crds_sgn <- terra::as.data.frame(x = agex_sgn, xy = T, cell = T, na.rm = T) # Coordinates

# Clean weather signatures from random miss-classification issues
vls <- terra::values(agex_sgn, na.rm = T) |> as.numeric() |> unique() |> sort()
agex_sgn_cln <- lapply(vls, FUN = function(i){
  agex_sgn_cln <- agex_sgn
  agex_sgn_cln[agex_sgn_cln != i] <- NA
  agex_sgn_cln[!is.na(agex_sgn_cln)] <- 1
  agex_sgn_cln <- terra::focal(x = agex_sgn_cln, w = 3, fun = 'mean')
  agex_sgn_cln[!is.na(agex_sgn_cln)] <- i
  return(agex_sgn_cln)
}) |> terra::sprc()
agex_sgn_cln <- terra::merge(agex_sgn_cln)
terra::writeRaster(x = agex_sgn_cln, filename = paste0(root,'/agroclimExtremes/agex_results/clusters/agex_global_',index,'_',gs,'_s',season,'_fmadogram_k192_clean.tif'))

crds_sgn_cln <- terra::as.data.frame(x = agex_sgn_cln, xy = T, cell = T, na.rm = T) # Coordinates
names(crds_sgn_cln)[ncol(crds_sgn_cln)] <- 'extreme_signature'

# Load global shapefile
wrl <- rnaturalearth::ne_countries(scale = 'large', returnclass = 'sv')
wrl <- wrl[,c('adm0_iso','name_en','continent','subregion')]

# Identify pixels within countries
crds_cty_cln <- terra::intersect(x = wrl, terra::vect(crds_sgn_cln, c('x','y'), crs = 'EPSG:4326')) |> base::as.data.frame()
crds_cty_cln <- unique(crds_cty_cln[,c('adm0_iso','continent','subregion','extreme_signature')])
crds_cty_cln <- dplyr::arrange(.data = crds_cty_cln, extreme_signature) |> base::as.data.frame()

# Identify collaboration opportunities among countries
vls <- sort(unique(crds_cty_cln$extreme_signature))
collaborations <- lapply(vls, function(i){
  smm <- data.frame(extreme_signature = i,
                    countries_count = length(crds_cty_cln[crds_cty_cln$extreme_signature == i,'adm0_iso']),
                    countries = paste0(crds_cty_cln[crds_cty_cln$extreme_signature == i,'adm0_iso'], collapse = ','),
                    continents = paste0(unique(crds_cty_cln[crds_cty_cln$extreme_signature == i,'continent']), collapse = ','),
                    regions = paste0(unique(crds_cty_cln[crds_cty_cln$extreme_signature == i,'subregion']), collapse = ','))
  return(smm)
}) |> dplyr::bind_rows()

collaborations <- collaborations[collaborations$extreme_signature != 3,]

# Missing signatures
base::setdiff(1:length(unique(as.numeric(terra::values(agex_sgn, na.rm = T)))), collaborations$extreme_signature)

collaborations <- rbind(collaborations,
                        data.frame(extreme_signature=155, countries_count=1, countries='IDN', continents='Asia', regions='South-Eastern Asia'))
collaborations <- dplyr::arrange(.data = collaborations, -countries_count) |> base::as.data.frame()
utils::write.csv(x = collaborations, file = paste0(root,'/agroclimExtremes/agex_results/collaborations_',index,'_',gs,'_s',season,'.csv'), row.names = F)

# Cluster quality measures
# They are computed to get rid of anomalous signatures
# Cohesion index
# 0 if patches of class become more isolated
# 100 if patches of class i become more aggregated
agex_chs <- landscapemetrics::lsm_c_cohesion(landscape = agex_sgn) |> base::as.data.frame()
agex_chs <- agex_chs[,c('class','value')]
names(agex_chs) <- c('extreme_signature','sgn_cohesion')
# Contiguity index mean: equals the mean of the contiguity index on class level for all patches
# 0 maximally non-contiguous
# 1 maximally contiguous
agex_ctg <- landscapemetrics::lsm_c_contig_mn(landscape = agex_sgn) |> base::as.data.frame()
agex_ctg <- agex_ctg[,c('class','value')]
names(agex_ctg) <- c('extreme_signature','sgn_contiguity')
# Aggregation index
# 0 maximally disaggregated
# 100 maximally aggregated
agex_agr <- landscapemetrics::lsm_c_ai(landscape = agex_sgn) |> base::as.data.frame()
agex_agr <- agex_agr[,c('class','value')]
names(agex_agr) <- c('extreme_signature','sgn_aggregation')

qlt_mtrcs <- dplyr::left_join(x = agex_chs, y = agex_ctg, by = 'extreme_signature')
qlt_mtrcs <- dplyr::left_join(x = qlt_mtrcs, y = agex_agr, by = 'extreme_signature')
rm(agex_chs, agex_ctg, agex_agr)

agex_qly_pca <- stats::prcomp(x = qlt_mtrcs[,-1], retx = T, center = T, scale. = T)
qlt_mtrcs$rank <- rank(agex_qly_pca$x[,1])

# Load global shapefile
wrl <- geodata::world(resolution = 1, level = 0, path = tempdir())

# Compute diversity of crop types
# High entropy means high variation of food/non-food classes
crp_fls <- list.files(path = paste0(root,'/agroclimExtremes/agex_raw_data'), pattern = '^agex_cropclass_', full.names = T)
crp_fls <- crp_fls[-grep2(pattern = c('fibres','stimulant'), x = crp_fls)]
crp_typ <- terra::rast(crp_fls)
crp_ntp <- terra::app(x = crp_typ, fun = function(i, ff) ff(i), cores = 20, ff = entropy::entropy)
names(crp_ntp) <- 'crop_types_diversity'
plot(crp_ntp)
crp_ntp_25km <- terra::resample(x = crp_ntp, y = agex_sgn, method = 'bilinear')

# Value of production
vop <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_spam2010V2r0_global_V_agg_VP_CROP_A.tif'))
names(vop) <- 'value_of_production'
vop_25km <- terra::resample(x = vop, y = agex_sgn, method = 'bilinear')

# tst <- terra::zonal(x = c(crp_ntp_25km, vop_25km), z = agex_sgn, fun = 'mean', na.rm = T)
# plot(tst$crop_types_diversity, tst$value_of_production, pch = 20)
# 
# tst |>
#   ggplot2::ggplot(aes(x = crop_types_diversity, y = value_of_production)) +
#   ggplot2::geom_point() +
#   ggplot2::geom_smooth(se = F)

# Index severity
idx <- terra::rast(paste0(root,'/agroclimExtremes/agex_indices/agex_',index,'/agex_',index,'_25km/',gs,'_s',season,'_',index,'_25km.tif'))
idx_avg <- mean(idx)
idx_max <- max(idx)
idx_mdn <- median(idx)
get_trend <- function(x){
  if(!all(is.na(x))){
    x <- na.omit(x) |> as.numeric()
    y <- trend::sens.slope(x)$estimates |> as.numeric()
  } else {
    y <- NA
  }
  return(y)
}
idx_slp <- terra::app(x = idx, fun = function(i, ff) ff(i), cores = 20, ff = get_trend)
plot(idx_slp)

get_p95_trend <- function(x){
  if(!all(is.na(x))){
    x <- as.numeric(na.omit(x))
    dfm <- data.frame(time = 1:length(x), ts = x)
    qrf <- quantreg::rq(formula = ts ~ time, tau = .95, data = dfm)
    y <- as.numeric(qrf$coefficients[2])
  } else {
    y <- NA
  }
  return(y)
}

tst <- terra::zonal(x = c(idx_slp, crp_ntp_25km), z = agex_sgn, fun = 'mean', na.rm = T)
names(tst)[-1] <- c('Response_speed','Food_diversity')

tst <- biscale::bi_class(.data = tst, x = Food_diversity, y = Response_speed, style = 'quantile', dim = 3)
colours <- data.frame(bi_class = sort(unique(tst$bi_class)))
colours$col <- pals::stevens.bluered(n = nrow(colours))
tst <- dplyr::left_join(x = tst, y = colours, by = 'bi_class') |> base::as.data.frame()

tst <- dplyr::arrange(.data = tst, -slope) |> base::as.data.frame()
n_col <- length(unique(terra::values(agex_sgn,na.rm = T)))
col_pltt <- MetBrewer::met.brewer(name = 'OKeeffe1', n = n_col)
tst$palette <- col_pltt
tst <- dplyr::arrange(.data = tst, extreme_signature) |> base::as.data.frame()

auxl_plt <- agex_sgn
auxl_plt <- terra::trim(x = auxl_plt)
# Figure
png(filename = paste0(root,'/agroclimExtremes/agex_global_',index,'_',gs,'_s',season,'_fmadogram.png'), width = 3132, height = 2359, units = 'px')
plot(wrl, ext = terra::ext(auxl_plt))
plot(auxl_plt, add = T, col = tst$palette, plg = list(cex = 5), cex.main = 7)
dev.off()

# Figure
png(filename = paste0(root,'/agroclimExtremes/agex_global_',index,'_',gs,'_s',season,'_fmadogram_bivariate.png'), width = 3132, height = 2359, units = 'px')
plot(wrl, ext = terra::ext(auxl_plt))
plot(auxl_plt, add = T, col = tst$col, plg = list(cex = 5), cex.main = 7)
dev.off()

# Legend box
tst$x <- as.numeric(substr(x = tst$bi_class, start=1, stop=1))
tst$y <- as.numeric(substr(x = tst$bi_class, start=3, stop=3))
tst %>%
  ggplot2::ggplot(aes(x = x, y = y)) +
  ggplot2::geom_tile(fill = tst$col, alpha = 1) +
  ggplot2::coord_equal() +
  ggplot2::theme_minimal() +
  ggplot2::xlab('Food classes diversity') +
  ggplot2::ylab('Response speed to drought') +
  ggplot2::theme(axis.text       = element_blank(),
                 axis.title      = element_text(size = 35),
                 legend.text     = element_text(size = 17),
                 legend.title    = element_blank(),
                 plot.title      = element_text(size = 25),
                 plot.subtitle   = element_text(size = 17),
                 plot.caption    = element_text(size = 15, hjust = 0),
                 legend.position = "bottom")

tst <- terra::zonal(x = c(crp_ntp_25km, vop_25km, cdd_avg, cdd_max, cdd_mdn), z = agex_sgn, fun = 'mean', na.rm = T)

tst |>
  ggplot2::ggplot(aes(x = median, y = crop_types_diversity)) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(se = F)

# Value of production per tertiles
vTert <- stats::quantile(tst$value_of_production, c(0:3/3))
tst$vop_tert = with(tst, cut(value_of_production, vTert, include.lowest = T, labels = c('Low', 'Medium', 'High')))
# CDD severity per tertiles
vTert <- stats::quantile(tst$median, c(0:3/3))
tst$cdd_tert = with(tst, cut(median, vTert, include.lowest = T, labels = c('Low', 'Medium', 'High')))

table(tst$cdd_tert, tst$vop_tert)

View(tst)

tst <- tst |> dplyr::arrange(value_of_production) |> base::as.data.frame()
tst$vop_colors <- NA
tst$vop_colors[tst$vop_tert == 'Low'] <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,'YlGn'))(sum(tst$vop_tert == 'Low'))
tst$vop_colors[tst$vop_tert == 'Medium'] <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,'Oranges'))(sum(tst$vop_tert == 'Medium'))
tst$vop_colors[tst$vop_tert == 'High'] <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,'Reds'))(sum(tst$vop_tert == 'High'))

tst <- tst |> dplyr::arrange(extreme_signature) |> base::as.data.frame()

# Population density
pop <- geodata::population(year = 2020, res = 10, path = tempdir())
pop_25km <- terra::resample(x = pop, y = agex_sgn, method = 'sum')

tst <- terra::zonal(x = c(crp_ntp_25km, vop_25km, idx_avg, idx_max, idx_mdn, pop_25km), z = agex_sgn, fun = 'mean', na.rm = T)
tst[,2:ncol(tst)] <- round(tst[,2:ncol(tst)], 2)

tst[,-1] |> cor(method = 'kendall') |> round(2)

## Graph of global cluster
wrl <- rnaturalearth::ne_countries(scale = 'large', returnclass = 'sv')
wrl <- terra::aggregate(wrl)
tmp <- agex_sgn
tmp <- terra::trim(x = tmp)
hist(terra::values(tmp))
png(filename = "D:/global_cluster_3palettes.png", width = 3132, height = 2359, units = 'px')
plot(wrl, ext = terra::ext(tmp))
plot(tmp, add = T, col = tst$vop_colors, plg = list(cex = 5), cex.main = 7)
dev.off()

biscale::bi_class(.data = tbl, x = Heat, y = Waterlogging, style = "quantile", dim = 3)

lapply(1:length(unique(crds_sgn$extreme_signature)), function(i){
  cat(paste0('>>> Extracting features from extreme signature: ',i,'\n'))
  cat('1. Verify coordinates within countries\n')
  crds_cty <- terra::intersect(x = wrl, terra::vect(crds_sgn[crds_sgn$extreme_signature == i,], c('x','y'), crs = 'EPSG:4326')) |> base::as.data.frame()
  cat('2. Get the full list of countries\n')
  countries_list <- unique(crds_cty$NAME_0)
  cat(paste0('   We got: ',length(countries_list),' countries\n'))
  cat('3. Count the number of pixels within each country\n')
  cnts_cty <- table(crds_cty$NAME_0) |> base::as.data.frame()
  names(cnts_cty)[1] <- 'Country'
  cat('4. Check adjacency among countries\n')
  if(length(countries_list) > 1){
    # Filtering global shapefile
    ctys_flt <- wrl[wrl$NAME_0 %in% countries_list,]
    # Find neighbor countries
    adjc_cty <- terra::relate(x = ctys_flt, relation = 'touches')
    adjc_cty <- 1 * adjc_cty
    colnames(adjc_cty) <- rownames(adjc_cty) <- countries_list
    adjc_cty <- adjc_cty[apply(adjc_cty, 2, sd) != 0, apply(adjc_cty, 2, sd) != 0]
    countries_adj_status <- ifelse(dim(adjc_cty)[1] == 0, 'Independent', 'Adjacent')
    # Best separation
    if(countries_adj_status == 'Adjacent' & nrow(adjc_cty) > 3){
      sprt_cty <- NbClust::NbClust(diss     = as.dist(1-adjc_cty),
                                   distance = NULL,
                                   min.nc   = 2,
                                   max.nc   = nrow(adjc_cty)-1,
                                   method   = 'ward.D2',
                                   index    = c('silhouette'))
      if(is.infinite(sprt_cty$Best.nc[1])){
        rownames(sprt_cty$Best.partition)
      } else {
        print(sort(sprt_cty$Best.partition))
        cnts_cty[cnts_cty$Country %in% split(names(sprt_cty$Best.partition), sprt_cty$Best.partition)[[1]],] |> dplyr::arrange(-Freq)
      }
    } else {
      # The country or two countries
    }
  } else {
    
  }
  
})

data.frame(Extreme_signature = i, Countries = paste0(countries_list, collapse = ','))



mtrcs <- dplyr::left_join(x = agex_chs, y = agex_ctg, by = 'extreme_signature')
mtrcs <- dplyr::left_join(x = mtrcs, y = agex_agr, by = 'extreme_signature')
rm(agex_agr, agex_chs, agex_ctg)

mtrcs[,-1] <- round(mtrcs[,-1], 2)

names(tst)[4:6] <- c('SPEI-6_mean','SPEI-6_max','SPEI-6_median')

tst <- dplyr::left_join(x = tst, y = mtrcs, by = 'extreme_signature')

agex_sgn_poly <- terra::as.polygons(x = agex_sgn)
agex_sgn_poly <- terra::merge(x = agex_sgn_poly, y = tst)

terra::writeVector(agex_sgn_poly,'D:/agex_features.gpkg')
