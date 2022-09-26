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




