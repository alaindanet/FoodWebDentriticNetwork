## Load your packages, e.g. library(targets).
source("./packages.R")

## Load your R files
lapply(list.files("./R", full.names = TRUE), source)
# Building community metrics
source("https://raw.githubusercontent.com/alaindanet/fishcom/master/R/community_analysis.R")

# Where the fish data lives
fishcom_path <- "~/Documents/post-these/mnhn/fishcom/"

## tar_plan supports drake-style targets and also tar_target()
tar_plan(

  # Collect clean files from my postdoc project:
  #1. AFB database database extraction is in https://github.com/alaindanet/fishcom/blob/master/vignettes/aspe_database.Rmd
  #2. Data cleaning is in  https://github.com/alaindanet/fishcom/blob/master/data-raw/01_clean_fish_op_data.R
  tar_target(op_sp_ind_file,
    paste0(fishcom_path, "data/op_sp_ind.rda"),
    format = "file"),
  tar_target(fish_length_file,
    paste0(fishcom_path, "data/fish_length.rda"),
    format = "file"),
  tar_target(op_file,
    paste0(fishcom_path, "data-raw/fishing_op_build/op.rda"),
    format = "file"),
  tar_target(op_desc_file,
    paste0(fishcom_path, "data-raw/fishing_op_build/op_desc.rda"),
    format = "file"),
  tar_target(op_env_file,
    paste0(fishcom_path, "data-raw/fishing_op_build/op_env.rda"),
    format = "file"),
  tar_target(station_file,
    paste0(fishcom_path, "data-raw/fishing_op_build/station.rda"),
    format = "file"),

  # Collect ontogenic data diet and predation windows
  tar_target(pred_win_file,
    paste0(fishcom_path, "data/pred_win.rda"),
    format = "file"),
  tar_target(resource_diet_shift_file,
    paste0(fishcom_path, "data/resource_diet_shift.rda"),
    format = "file"),
  tar_target(fish_diet_shift_file,
    paste0(fishcom_path, "data/fish_diet_shift.rda"),
    format = "file"),

  # Load the data
  tar_target(op_sp_ind, return_data(op_sp_ind_file)),
  tar_target(op, return_data(op_file)),
  tar_target(op_desc, return_data(op_desc_file) %>%
    rename(opcod = ope_id)
    ),
  tar_target(op_env, return_data(op_env_file)),
  tar_target(station, return_data(station_file) %>%
  rename(station = id)),
  tar_target(pred_win, return_data(pred_win_file)),
  tar_target(resource_diet_shift, return_data(resource_diet_shift_file)),
  tar_target(fish_diet_shift, return_data(fish_diet_shift_file)),

  # Individual fish length
  tar_target(fish_length,
    return_data(fish_length_file) %>%
        filter(species %in% unique(fish_diet_shift$species)) %>%
        #remove because only 3 individuals (not enough for size class)
        filter(species != "APR") %>%
        left_join(select(op, station, opcod), by = "opcod") %>%
        # Get bodymass estimation
        mutate(weight = calc_fish_body_mass(length, unit = "gram"))
    ),

  # Network
  tar_target(metaweb,
    build_metaweb(
      data = fish_length,
      species = species,
      size = length,
      pred_win = pred_win,
      beta_min = beta_min,
      beta_max = beta_max,
      fish_diet_shift = fish_diet_shift,
      low_bound = size_min,
      upper_bound = size_max,
      fish = fish,
      resource_diet_shift = resource_diet_shift,
      class_method = "percentile",
      nb_class = 9,
      pred_win_method = "midpoint",
      fish_resource_method = "midpoint",
      na.rm = TRUE,
      replace_min_by_one = FALSE)
    ),
  tar_target(network,
    build_local_network(
      data = fish_length,
      species = species,
      var = length,
      group_var = opcod,
      metaweb = metaweb,
      classes = NULL,
      out_format = "igraph") %>%
    left_join(select(op, station, opcod), by = "opcod") %>%
    ungroup()
    ),

  # Community
  tar_target(com_species,
    fish_length %>%
      group_by(opcod, species) %>%
      summarise(
        nind = n(),
        biomass = sum(weight, na.rm =TRUE),
        length = mean(length, na.rm = TRUE)
      ) %>%
      left_join(select(op, opcod, surface), by = "opcod") %>%
      mutate(
        nind_std = nind / surface,
        bm_std = biomass / surface
      ) %>%
      ungroup()
  ),
  tar_target(com_metrics,
    com_species %>%
  group_by(opcod) %>%
  summarise(
    richness = n(),
    nind = sum(nind),
    biomass = sum(biomass)
  ) %>%
  left_join(select(op, opcod, surface), by = "opcod") %>%
  mutate(
    rich_std = richness / surface,
    nind_std = nind / surface,
    bm_std = biomass / surface
  )
    ),

  # Network metrics
  tar_target(com_trophic_species, {
    # Sanatize dataset according to the metaweb 
    tmp <- sanatize_metaweb(data = fish_length,
      species = species,
      fish_diet_shift = metaweb$size_class,
      nb_class = max(unique(metaweb$size_class$class_id))
      ) %>%
      filter(!is.na(length))

    # Get trophic species for each ind
    trophic_species_bm <- assign_size_class(
      tmp, species, var = length,
      classes = metaweb$size_class) %>%
      unite(sp_class, species, class_id, sep = "_") %>%
      dplyr::select(opcod, sp_class, weight)

    # Get summary metrics for each trophic species
    trophic_species_bm %>%
      group_by(opcod, sp_class) %>%
      summarise(
        nind = n(),
        biomass = sum(weight)
        ) %>%
      left_join(select(op, opcod, surface), by = "opcod") %>%
      mutate(
        nind_std = nind / surface,
        bm_std = biomass / surface
      ) %>%
      ungroup()
    }
    ),
  tar_target(network_mat,
    network %>%
      select(station, opcod, network) %>%
      mutate(
        network = future_map(network, igraph::graph_from_data_frame, directed = TRUE),
        network = future_map(network, igraph::as_adjacency_matrix, sparse = FALSE),
        troph = future_map(network, NetIndices::TrophInd),
        metrics = future_map(network, NetIndices::GenInd),
        obs_troph_level = map(troph,
          function (x) {
            out <- tibble(
              sp_class = row.names(x),
              obs_troph_level = x$TL
            )
            return(out)
          }
        )
      )
    ),
  tar_target(network_metrics, {
    tlvl <- network_mat %>%
      select(opcod, obs_troph_level) %>%
      unnest(obs_troph_level) %>%
      left_join(
        com_trophic_species %>%
          select(opcod, sp_class, bm_std),
        by = c("opcod", "sp_class")
        ) %>%
      # Excluding resources for which we do not have biomass
      filter(!sp_class %in% metaweb$resource) %>%
      group_by(opcod) %>%
      summarise(
        weighted_avg_obs_tlvl = sum(
          obs_troph_level * bm_std / sum(bm_std, na.rm = T),
          na.rm = T),
        non_weighted_avg_obs_tlvl = mean(obs_troph_level)
      )

      network_mat %>%
        mutate(
          mean_troph_level = map_dbl(obs_troph_level, ~mean(.x$obs_troph_level)),
          max_troph_lvl = map_dbl(obs_troph_level, ~max(.x$obs_troph_level)),
          connectance = map_dbl(metrics, "C"),
          nbnode = map_dbl(metrics, "N")
          ) %>%
      select(station, opcod, mean_troph_level, max_troph_lvl, connectance, nbnode) %>%
      left_join(tlvl, by = "opcod")
    }
  ),

  # Info station / sampling filtered according to network
  tar_target(op_st_filtered,
    list(
      op = op %>% filter(opcod %in% network$opcod),
      station = station %>%
        filter(station %in% unique(network$station)),
      op_desc = op_desc %>% filter(opcod %in% network$opcod),
      op_env = op_env %>% filter(opcod %in% network$opcod)
      )),
  # Rivers
  tar_target(hydroriver_shp_files,
    here::here("data-raw", "HydroRIVERS_v10_eu_shp", "HydroRIVERS_v10_eu.shp"),
    format = "file"),
  tar_target(snapped_site_river,
    target_snap_site_to_river(
      river_shp_filepath = hydroriver_shp_files,
      site_sf = op_st_filtered$station %>%
        st_transform(crs = 4326),
      proj_crs =  4326,
      length_chunk = 200,
      max_dist = 1000
        )
  )
)