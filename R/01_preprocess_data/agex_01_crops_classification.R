# ------------------------------------------ #
# Crops classification - croplands
# By: Harold Achicanoy
# WUR & ABC
# Feb. 2024
# ------------------------------------------ #

# R options and packages loading
options(warn = -1, scipen = 999)

# Load MapSPAM readme
mapspam_info <- readLines('D:/Data/Maps/spam2010/ReadMe_v2r0_Global.txt')

# Food crops
crps <- strsplit(x = mapspam_info[105:132], split = '\t', fixed = T)
nms <- crps[[1]][!(nchar(crps[[1]]) == 0)]
crps <- lapply(2:length(crps), function(i){
  dfm <- data.frame(t(crps[[i]][!(nchar(crps[[i]]) == 0)]))
  names(dfm) <- nms
  return(dfm)
}) |> dplyr::bind_rows()

# Non-food crops
ncrps <- strsplit(x = mapspam_info[136:151], split = '\t', fixed = T)
ncrps[[10]] <- c('31','other fibre crops', 'ofib')
ncrps <- lapply(2:length(ncrps), function(i){
  dfm <- data.frame(t(ncrps[[i]][!(nchar(ncrps[[i]]) == 0)]))
  names(dfm) <- nms
  return(dfm)
}) |> dplyr::bind_rows()

# Put them together
agr <- rbind(crps, ncrps)
agr$`crop #` <- as.numeric(agr$`crop #`)
agr <- agr |> dplyr::arrange(`crop #`) |> base::as.data.frame()

# Crops classification
grp <- read.csv('https://raw.githubusercontent.com/wri/MAPSPAM/master/metadata_tables/4-Methodology-Crops-of-SPAM-2005-2015-02-26.csv')
grp <- grp[,c('SPAM.short.name','GROUP')]

# Add classification column
agr <- dplyr::left_join(x = agr, y = grp, by = c('SPAM name' = 'SPAM.short.name'))
agr$GROUP[agr$`SPAM name` == 'rest'] <- 'rest of crops'
