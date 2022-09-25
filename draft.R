library(tidyverse)
library(magrittr)
library(sizeTrophicInteractions)

tar_load(c(network))

library(sf)
tar_load(c(station, op, op_desc, op_env))
station

tar_load(c(network_mat, com_trophic_species, metaweb))

station %>%
  filter(id %in% unique(network_mat$station))
library(mapview)
mapview(station)

tar_load(op_st_filtered)
op_st_filtered$station
library(sf)

library("googledrive")
drive_find()
