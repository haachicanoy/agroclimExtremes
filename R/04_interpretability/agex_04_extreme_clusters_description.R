# ------------------------------------------ #
# Description of extreme weather clusters
# By: Harold Achicanoy
# WUR & ABC
# May 2024
# ------------------------------------------ #

# R options and packages loading
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman')}else{library(pacman)})
suppressMessages(pacman::p_load(terra,geodata,NbClust,entropy,scales,dplyr,
                                ggplot2,RColorBrewer,landscapemetrics,
                                exactextractr,hrbrthemes,trend,quantreg,
                                OutliersO3))

list.files2 <- Vectorize(FUN = list.files, vectorize.args = 'path')
grep2 <- Vectorize(FUN = grep, vectorize.args = 'pattern')
get_trend <- function(x){
  if(!all(is.na(x))){
    x <- as.numeric(na.omit(x))
    y <- as.numeric(trend::sens.slope(x)$estimates)
  } else { y <- NA }
  return(y)
}
get_ext_trend <- function(x){
  if(!all(is.na(x))){
    x <- as.numeric(na.omit(x))
    dfm <- data.frame(time = 1:length(x), ts = x)
    qrf <- quantreg::rq(formula = ts ~ time, tau = .90, data = dfm)
    y <- as.numeric(qrf$coefficients[2])
  } else { y <- NA }
  return(y)
}
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

## Key arguments
root   <- '//CATALOGUE/WFP_ClimateRiskPr1'
index  <- 'spei-6'

# Load extreme weather clusters
agex_sgn_cln <- terra::rast(paste0(root,'/agroclimExtremes/agex_results/agex_results_clusters/agex_global_spei-6_combined_fmadogram_clean.tif'))
names(agex_sgn_cln) <- 'extreme_cluster'
crds_sgn_cln <- terra::as.data.frame(x = agex_sgn_cln, xy = T, cell = T, na.rm = T) # Coordinates

# Load global shapefile
wrl <- rnaturalearth::ne_countries(scale = 'large', returnclass = 'sv')
wrl <- wrl[,c('adm0_iso','name_en','continent','subregion')]

# Pixels within countries
crds_cty <- terra::intersect(x = wrl, y = terra::vect(crds_sgn_cln, c('x','y'), crs = 'EPSG:4326')) |> base::as.data.frame()

# Pixels number per country
crds_cnt <- crds_cty |>
  dplyr::group_by(adm0_iso, extreme_cluster) |>
  dplyr::count() |>
  dplyr::arrange(extreme_cluster) |>
  dplyr::ungroup() |>
  base::as.data.frame()

# Unique geographies
unq_geos <- unique(crds_cty[,c('adm0_iso','name_en','continent','subregion','extreme_cluster')])
unq_geos <- base::as.data.frame(dplyr::arrange(.data = unq_geos, extreme_cluster))

# Identify collaboration opportunities among countries
vls <- sort(unique(unq_geos$extreme_cluster))
collaborations <- lapply(vls, function(i){
  smm <- data.frame(extreme_cluster = i,
                    countries_count = length(unq_geos[unq_geos$extreme_cluster == i,'adm0_iso']),
                    isos = paste0(unq_geos[unq_geos$extreme_cluster == i,'adm0_iso'], collapse = ','),
                    countries = paste0(unq_geos[unq_geos$extreme_cluster == i,'name_en'], collapse = ','),
                    continents = paste0(unique(unq_geos[unq_geos$extreme_cluster == i,'continent']), collapse = ',') # ,
                    # regions = paste0(unique(unq_geos[unq_geos$extreme_cluster == i,'subregion']), collapse = ',')
  )
  return(smm)
}) |> dplyr::bind_rows() |>
  dplyr::arrange(-countries_count) |>
  base::as.data.frame()

utils::write.csv(x = collaborations, file = paste0(root,'/agroclimExtremes/agex_results/agex_country_collaborations.csv'), row.names = F)

# Cluster quality measures
# They are computed to provide a quality control of the produced clusters
# Cohesion index
# 0 if patches of class become more isolated
# 100 if patches of class i become more aggregated
agex_chs <- landscapemetrics::lsm_c_cohesion(landscape = agex_sgn_cln) |> base::as.data.frame()
agex_chs <- agex_chs[,c('class','value')]
names(agex_chs) <- c('extreme_cluster','cls_cohesion')
# Contiguity index mean: equals the mean of the contiguity index on class level for all patches
# 0 maximally non-contiguous
# 1 maximally contiguous
agex_ctg <- landscapemetrics::lsm_c_contig_mn(landscape = agex_sgn_cln) |> base::as.data.frame()
agex_ctg <- agex_ctg[,c('class','value')]
names(agex_ctg) <- c('extreme_cluster','cls_contiguity')
# Aggregation index
# 0 maximally disaggregated
# 100 maximally aggregated
agex_agr <- landscapemetrics::lsm_c_ai(landscape = agex_sgn_cln) |> base::as.data.frame()
agex_agr <- agex_agr[,c('class','value')]
names(agex_agr) <- c('extreme_cluster','cls_aggregation')

qlt_mtrcs <- dplyr::left_join(x = agex_chs, y = agex_ctg, by = 'extreme_cluster')
qlt_mtrcs <- dplyr::left_join(x = qlt_mtrcs, y = agex_agr, by = 'extreme_cluster')
qlt_mtrcs[,-1] <- round(qlt_mtrcs[,-1], 3)
rm(agex_chs, agex_ctg, agex_agr)

agex_qly_pca <- stats::prcomp(x = qlt_mtrcs[,-1], retx = T, center = T, scale. = T)
qlt_mtrcs$quality_rank <- rank(-1*agex_qly_pca$x[,1])

# Compute diversity of crop types
# High entropy means high variation of food/non-food classes
crp_fls <- list.files(path = paste0(root,'/agroclimExtremes/agex_raw_data'), pattern = '^agex_cropclass_', full.names = T)
crp_fls <- crp_fls[-grep2(pattern = c('fibres','stimulant'), x = crp_fls)]
crp_typ <- terra::rast(crp_fls)
crp_ntp <- terra::app(x = crp_typ, fun = function(i, ff) ff(i), cores = 20, ff = entropy::entropy)
names(crp_ntp) <- 'crop_types_diversity'
plot(crp_ntp)
crp_ntp_25km <- terra::resample(x = crp_ntp, y = agex_sgn_cln, method = 'cubicspline', threads = T)
crp_ntp_25km <- terra::mask(x = crp_ntp_25km, mask = agex_sgn_cln)

# Value of production
vop <- terra::rast(paste0(root,'/agroclimExtremes/agex_raw_data/agex_spam2010V2r0_global_V_agg_VP_CROP_A.tif'))
names(vop) <- 'value_of_production'
vop_25km <- terra::resample(x = vop, y = agex_sgn_cln, method = 'cubicspline', threads = T)
vop_25km <- terra::mask(x = vop_25km, mask = agex_sgn_cln)

# Compute diversity of livestock
lvs_dir <- 'D:/OneDrive - CGIAR/African_Crisis_Observatory/data/_global/livestock'
anmls   <- list.dirs(path = lvs_dir, full.names = F, recursive = F)
anmls   <- anmls[-grep('horses',anmls)]
lvstc_fls <- list.files2(path = paste0(lvs_dir,'/',anmls), pattern = '_Da.tif$', full.names = T); rm(lvs_dir, anmls)
names(lvstc_fls) <- NULL
lvstc_cnt <- terra::rast(lvstc_fls)

# Computing Livestock Units
# LU: Livestock Unit
# LU = Buffaloes * 1 + Cattle * 1 + Chickens * ((0.007 + 0.014)/2) +
#      Ducks * 0.01 + Goats * 0.1 +Pigs * ((0.5+0.027)/2) + Sheep * 0.1
# Source: https://ec.europa.eu/eurostat/statistics-explained/index.php?title=Glossary:Livestock_unit_(LSU)
lsu <- lvstc_cnt[[1]] + lvstc_cnt[[2]] + lvstc_cnt[[3]]*((0.007 + 0.014)/2) + lvstc_cnt[[4]]*0.01 + lvstc_cnt[[5]]*0.1 + lvstc_cnt[[6]]*((0.5+0.027)/2) + lvstc_cnt[[7]]*0.1
rm(lvstc_fls, lvstc_cnt)
names(lsu) <- 'livestock_units'
lsu_25km <- terra::resample(x = lsu, y = agex_sgn_cln, method = 'cubicspline', threads = T)
lsu_25km <- terra::mask(x = lsu_25km, mask = agex_sgn_cln)

# Population density
pop <- geodata::population(year = 2020, res = 10, path = tempdir())
pop_25km <- terra::resample(x = pop, y = agex_sgn_cln, method = 'cubicspline', threads = T)
pop_25km <- terra::mask(x = pop_25km, mask = agex_sgn_cln)

# Index severity. TO DO
idx <- terra::rast(paste0(root,'/agroclimExtremes/agex_indices/agex_',index,'/agex_',index,'_25km/',gs,'_s',season,'_',index,'_25km.tif'))
idx_avg <- mean(idx); names(idx_avg) <- 'SPEI-6_mean'
idx_mdn <- median(idx); names(idx_mdn) <- 'SPEI-6_median'
idx_max <- max(idx); names(idx_max) <- 'SPEI-6_max'
idx_slp <- terra::app(x = idx, fun = function(i, ff) ff(i), cores = 20, ff = get_trend)
names(idx_slp) <- 'SPEI-6_slope'
idx_sxt <- terra::app(x = idx, fun = function(i, ff) ff(i), cores = 20, ff = get_ext_trend)
names(idx_sxt) <- 'SPEI-6_slope_90th'

rst_mtrcs <- terra::zonal(x = c(pop_25km, vop_25km, crp_ntp_25km, lsu_25km,
                                idx_avg, idx_mdn, idx_max, idx_slp, idx_sxt),
                          z = agex_sgn_cln, fun = 'mean', na.rm = T)
rst_mtrcs[,-1] <- round(rst_mtrcs[,-1], 3)

dfm <- dplyr::left_join(x = collaborations, y = rst_mtrcs, by = 'extreme_cluster')
dfm <- dplyr::left_join(x = dfm, y = qlt_mtrcs, by = 'extreme_cluster')

agex_sgn_poly <- terra::as.polygons(x = agex_sgn_cln)
agex_sgn_poly <- terra::merge(x = agex_sgn_poly, y = dfm)

terra::writeVector(x = agex_sgn_poly, filename = paste0(root,'/agroclimExtremes/agex_results/clusters/vct_agex_global_',index,'_',gs,'_s',season,'_fmadogram_k',nrow(dfm),'.gpkg'), overwrite = T)








# Figure 1
tst <- biscale::bi_class(.data = dfm, x = crop_types_diversity, y = `SPEI-6_slope`, style = 'quantile', dim = 3)
colours <- data.frame(bi_class = sort(unique(tst$bi_class)))
colours$col <- pals::stevens.bluered(n = nrow(colours))
tst <- dplyr::left_join(x = tst, y = colours, by = 'bi_class') |> base::as.data.frame()

tst <- dplyr::arrange(.data = tst, -`SPEI-6_slope`) |> base::as.data.frame()
n_col <- length(unique(terra::values(agex_sgn,na.rm = T)))
col_pltt <- MetBrewer::met.brewer(name = 'OKeeffe1', n = n_col)
tst$palette <- col_pltt
tst <- dplyr::arrange(.data = tst, extreme_cluster) |> base::as.data.frame()

auxl_plt <- agex_sgn
auxl_plt <- terra::trim(x = auxl_plt)

png(filename = paste0(root,'/agroclimExtremes/agex_global_',index,'_',gs,'_s',season,'_fmadogram_bivariate_crops.png'), width = 3132, height = 2359, units = 'px')
plot(wrl, ext = terra::ext(auxl_plt))
plot(auxl_plt, add = T, col = tst$col, plg = list(cex = 5), cex.main = 7)
dev.off()

# Figure 2
tst <- biscale::bi_class(.data = dfm, x = livestock_units, y = `SPEI-6_slope`, style = 'quantile', dim = 3)
colours <- data.frame(bi_class = sort(unique(tst$bi_class)))
colours$col <- pals::stevens.bluered(n = nrow(colours))
tst <- dplyr::left_join(x = tst, y = colours, by = 'bi_class') |> base::as.data.frame()

tst <- dplyr::arrange(.data = tst, -`SPEI-6_slope`) |> base::as.data.frame()
n_col <- length(unique(terra::values(agex_sgn,na.rm = T)))
col_pltt <- MetBrewer::met.brewer(name = 'OKeeffe1', n = n_col)
tst$palette <- col_pltt
tst <- dplyr::arrange(.data = tst, extreme_cluster) |> base::as.data.frame()

auxl_plt <- agex_sgn
auxl_plt <- terra::trim(x = auxl_plt)

png(filename = paste0(root,'/agroclimExtremes/agex_global_',index,'_',gs,'_s',season,'_fmadogram_bivariate_livestock.png'), width = 3132, height = 2359, units = 'px')
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
  ggplot2::xlab('Livestock diversity') +
  ggplot2::ylab('Drought trend') +
  ggplot2::theme(axis.text       = element_blank(),
                 axis.title      = element_text(size = 35),
                 legend.text     = element_text(size = 17),
                 legend.title    = element_blank(),
                 plot.title      = element_text(size = 25),
                 plot.subtitle   = element_text(size = 17),
                 plot.caption    = element_text(size = 15, hjust = 0),
                 legend.position = "bottom")

# Figure 3
dfm |>
  dplyr::select(extreme_cluster,`SPEI-6_slope`,`SPEI-6_slope_90th`) |>
  tidyr::pivot_longer(cols = -1, names_to = 'variable', values_to = 'value') |>
  base::as.data.frame() |>
  ggplot2::ggplot(aes(x = value, fill = variable)) +
  ggplot2::geom_density(color="#e9ecef", alpha=0.6, position = 'identity') +
  ggplot2::scale_fill_manual(values=c("#69b3a2", "#404080")) +
  ggplot2::geom_vline(xintercept=0, linetype=3, colour='red') +
  ggplot2::theme_bw()


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
