renv::init()
renv::activate()
renv::snapshot()

library(targets)
library(tidyverse)
library(magrittr)
##
#Make datasets for the workshop
dir.create("~/Téléchargements/NetworksInNetworks/data", recursive = TRUE)
dir.create("~/Téléchargements/NetworksInNetworks/data-raw/fishing_op_build", recursive = TRUE)

  
tar_load(c(op_sp_ind_file,
  fish_length_file,
  op_file,
  op_desc_file,
  op_env_file,
  station_file,
  pred_win_file,
  resource_diet_shift_file,
  fish_diet_shift_file))

file_cut <- map_chr(c(op_sp_ind_file,
  fish_length_file,
  op_file,
  op_desc_file,
  op_env_file,
  station_file,
  pred_win_file,
  resource_diet_shift_file,
  fish_diet_shift_file),
~str_remove(.x, "~/Documents.*com/")
)

fishcom_path <- "~/Documents/post-these/mnhn/fishcom/"
data_rep <- "~/Téléchargements/NetworksInNetworks/"

file.copy(
  from = paste0(fishcom_path, file_cut),
  to = paste0(data_rep, file_cut)
)

zip(
  zipfile = paste0(str_remove(data_rep, "/$"), ".zip"),
  files = dir(str_remove(data_rep, "/$"), full.names = TRUE))
