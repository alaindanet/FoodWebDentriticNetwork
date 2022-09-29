library(tidyverse)
library(magrittr)
library(sizeTrophicInteractions)
library(targets)
library(sf)
lapply(list.files("./R", full.names = TRUE), source)

tar_load(c(network))
tar_load(c(station, op, op_desc, op_env))
tar_load(c(network_mat, com_trophic_species, metaweb))
tar_load(op_st_filtered)
names(op_st_filtered)
tar_load(c(riveratlas_station))

ti <- c(get_river_atlas_significant_var(), get_land_class_var())
tar_load(c(hydroriver_shp_files, snapped_site_river))


save(basinl12, file = "~/Téléchargements/basin_station_L12.rda")

basin_dce <- read_sf("~/Documents/post-these/mnhn/fishcom/data-raw/basin_dce/BassinDCE.shp")

load("~/Documents/post-these/mnhn/fishcom/data/the_8_hydrologic_basin.rda")
basin_dce <- the_8_hydrologic_basin
save(basin_dce, file = "~/Téléchargements/basin_dce.rda")

op_st_filtered$station

intersect_st_basin <- st_intersects(op_st_filtered$station %>%
  st_transform(crs = 2154),
  basin_dce)
station_basin_dce <- tibble(
  station = op_st_filtered$station$station, 
  basin = basin_dce$NomDistric[map_int(intersect_st_basin, ~.x[1])]
)


tar_load(riveratlas_station)
rv_at <- riveratlas_station %>%
  left_join(station_basin_dce, by = "station")

ti <- rv_at %>%
  group_by(basin, ord_stra) %>%
  summarise(n = n())

tar_load(c(com_species, op_st_filtered, station_basin_dce))





