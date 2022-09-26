library(tidyverse)
library(magrittr)
library(sizeTrophicInteractions)
library(targets)
library(sf)
lapply(list.files("./R", full.names = TRUE), source)

tar_load(c(network))

tar_load(c(station, op, op_desc, op_env))
station

tar_load(c(network_mat, com_trophic_species, metaweb))

station %>%
  filter(id %in% unique(network_mat$station))
library(mapview)
mapview(station)

tar_load(op_st_filtered)
op_st_filtered$station

tar_load(c(riveratlas_station))
hydroriver <- read_sf(hydroriver_shp_files)


hydroriver %>%
  filter(HYRIV_ID %in% snapped_site_river$riverid)

layer_name <- sf::st_layers(riveratlas_shp_files, do_count = TRUE)$name
to_select <- setNames(
  unique(c("HYRIV_ID", "LENGTH_KM", "DIST_UP_KM", "ORD_STRA", "ORD_FLOW",
    get_river_atlas_significant_var(), get_land_class_var())),
  NULL)


myquery <- paste0(
  "SELECT ",
  paste0(to_select[!to_select %in% c("length_km", "dist_up_km", "ord_stra", "ord_flow", "hft_ix_c9309_ratio", "hft_ix_c9309_log2_ratio")], collapse = ", "),
  " FROM ", layer_name, " WHERE HYRIV_ID IN ",
  "(",
  paste0(na.omit(snapped_site_river$riverid), collapse = ","),
  ")"
)

river <- sf::read_sf(dsn = riveratlas_shp_files, query = myquery)

