# crops
crops_cls <- c(223, 420, 270, 33, 42, 444)
names(crops_cls) <- c('North America','South America','Africa','Europe','Asia','Australia')

# livestock
lvstc_cls <- c(227, 402, 289, 50, 171, 444)
names(lvstc_cls) <- c('North America','South America','Africa','Europe','Asia','Australia')

root <- '//CATALOGUE/AgroclimExtremes'
rpth <- paste0(root,'/agex_results')

agex_metrics <- read.csv(file.path(rpth,'agex_all_metrics.csv'))
agex_hvareas <- read.csv(file.path(rpth,'agex_harvested_areas_percentages.csv'))
agex_lvunits <- read.csv(file.path(rpth,'agex_livestock_units_percentages.csv'))

crops_cls_dfm <- data.frame(extreme_cluster = c(223, 420, 270, 33, 42, 444),
                            continent = c('North America','South America','Africa','Europe','Asia','Australia'))
lvstc_cls_dfm <- data.frame(extreme_cluster = c(227, 402, 289, 50, 171, 444),
                            continent = c('North America','South America','Africa','Europe','Asia','Australia'))

agex_crops_hotspots <- agex_hvareas[agex_hvareas$extreme_cluster %in% crops_cls,]
agex_crops_hotspots <- agex_crops_hotspots |>
  tidyr::pivot_longer(cols = -1, names_to = 'Crop', values_to = 'Percent') |>
  dplyr::filter(Percent > 5) |>
  dplyr::group_by(extreme_cluster) |>
  dplyr::arrange(extreme_cluster, -Percent) |>
  base::as.data.frame()
agex_crops_hotspots <- dplyr::left_join(x = agex_crops_hotspots, y = crops_cls_dfm, by = 'extreme_cluster')
(table(agex_crops_hotspots$continent)/34 * 100) |> round(1)

grp <- read.csv('https://raw.githubusercontent.com/wri/MAPSPAM/master/metadata_tables/4-Methodology-Crops-of-SPAM-2005-2015-02-26.csv')
grp <- grp[,c('SPAM.long.name','GROUP')]
grp$SPAM.long.name <- gsub(' ','.',stringr::str_to_title(grp$SPAM.long.name))
names(grp)[ncol(grp)] <- 'Class'
# Add classification column
agex_crops_hotspots <- dplyr::left_join(x = agex_crops_hotspots, y = grp,
                                        by = c('Crop' = 'SPAM.long.name'))
rm(grp)
utils::write.csv(x = agex_crops_hotspots, file = 'D:/crop_distribution_hotspots.csv', row.names = F)

agex_lvstc_hotspots <- agex_lvunits[agex_lvunits$extreme_cluster %in% lvstc_cls,]
agex_lvstc_hotspots <- dplyr::left_join(x = agex_lvstc_hotspots, y = lvstc_cls_dfm, by = 'extreme_cluster')
agex_lvstc_hotspots <- agex_lvstc_hotspots[,c('extreme_cluster','continent',base::setdiff(names(agex_lvstc_hotspots),c('extreme_cluster','continent')))]

utils::write.csv(x = agex_lvstc_hotspots, file = 'D:/livestock_distribution_hotspots.csv', row.names = F)
