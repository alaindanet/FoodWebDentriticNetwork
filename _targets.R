## Load your packages, e.g. library(targets).
source("./packages.R")

## Load your R files
lapply(list.files("./R", full.names = TRUE), source)
# Building community metrics
source("https://raw.githubusercontent.com/alaindanet/fishcom/master/R/community_analysis.R")
source("https://raw.githubusercontent.com/alaindanet/fishcom/master/R/codeWeb.R")

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
  tar_target(occ,{
    n_op_st <- op_st_filtered$op %>%
      group_by(station) %>%
      summarise(nop = n())

    com_species %>%
      left_join(op_st_filtered$op[, c("station", "opcod")], by = "opcod") %>%
      left_join(station_basin_dce, by = "station") %>%
      group_by(station, species) %>%
      summarise(occ = n(), .groups = "drop") %>%
      left_join(n_op_st, by = "station") %>%
      mutate(occ_prop = occ / nop)
    }),
  tar_target(occ_mat_basin, {
    occ %>%
      mutate(pres = occ_prop > .1) %>%
      filter(pres) %>%
      select(station, species, pres) %>%
      pivot_wider(names_from = "species", values_from = "pres") %>%
      mutate(across(where(is.logical), ~ ifelse(is.na(.), FALSE, .))) %>%
      left_join(station_basin_dce, by = "station") %>%
      select(basin, station, everything()) %>%
      filter(!is.na(basin)) %>%
      group_by(basin) %>%
      nest() %>%
      mutate(occ_mat = map(data,
          function(x) {
            mat_site_sp <- as.matrix(x[,!colnames(x) %in% "station"])
            mat_sp <- apply(mat_site_sp, 2, as.numeric)
            row.names(mat_sp) <- x$station
            mat_sp
          }
          ))
    }),
  tar_target(beta_div,
     occ_mat_basin %>%
       mutate(
         beta_tot = map(occ_mat, ~as_tibble(betapart::beta.multi(.x, index.family = "jaccard"))),
         beta_pair = map(occ_mat, ~betapart::beta.pair(.x, index.family = "jaccard")),
         beta_jac_st = map(beta_pair, ~enframe(colMeans(as.matrix(.x$beta.jac)), name = "station", value = "jaccard")),
         beta_tu_st = map(beta_pair, ~enframe(colMeans(as.matrix(.x$beta.jtu)), name = "station", value = "turnover")),
         beta_ne_st = map(beta_pair, ~enframe(colMeans(as.matrix(.x$beta.jne)), name = "station", value = "nestedness"))
       ) %>%
     select(-data, -occ_mat)
    ),
  tar_target(betadiv_site_by_basin,
    station_basin_dce %>%
      mutate(station = as.character(station)) %>%
      left_join(do.call(rbind, beta_div$beta_jac_st), by = "station") %>%
      left_join(do.call(rbind, beta_div$beta_tu_st), by = "station") %>%
      left_join(do.call(rbind, beta_div$beta_ne_st), by = "station")
    ),
  tar_target(beta_by_stra,
    beta_div %>%
      select(basin, beta_pair) %>%
      mutate(beta_by_stra = map(beta_pair, function(z) {
          map2_dfr(z, names(beta_div$beta_pair[[1]]), function(x, y) {
            x %>%
              as.matrix() %>%
              as.data.frame() %>%
              rownames_to_column(var = "station1") %>%
              pivot_longer(-station1, names_to = "station2", values_to = "beta") %>% # site pair lg df
              left_join(st_basin_width_ord_sta %>% #get ord_stra for each station of the pair
                  select(station1 = station, ord_stra1 = ord_stra),
                  by = "station1") %>%
              left_join(st_basin_width_ord_sta %>%
                select(station2 = station, ord_stra2 = ord_stra),
              by = "station2") %>%
              filter(ord_stra1 == ord_stra2, station1 != station2) %>% #garder les pairs de meme ord stra mais station differentes
              mutate(ord_stra = ord_stra1) %>%
                mutate(metric = y)
})
          })) %>%
      select(basin, beta_by_stra) %>%
      unnest(beta_by_stra) %>%
      mutate(metric = str_replace_all(metric, c("beta.jtu" = "turnover", "beta.jne" = "nestedness", "beta.jac" = "jaccard")))
    ),
  tar_target(beta_by_width,
    beta_div %>%
      select(basin, beta_pair) %>%
      mutate(beta_by_width = map(beta_pair, function(z) {
          map2_dfr(z, names(beta_div$beta_pair[[1]]), function(x, y) {
            x %>%
              as.matrix() %>%
              as.data.frame() %>%
              rownames_to_column(var = "station1") %>%
              pivot_longer(-station1, names_to = "station2", values_to = "beta") %>% # site pair lg df
              left_join(st_basin_width_ord_sta %>% #get ord_stra for each station of the pair
                select(station1 = station, mean.width_river_cat1 = mean.width_river_cat),
                by = "station1") %>%
              left_join(st_basin_width_ord_sta %>%
                select(station2 = station, mean.width_river_cat2 = mean.width_river_cat),
              by = "station2") %>%
              filter(mean.width_river_cat1 == mean.width_river_cat2, station1 != station2) %>% #garder les pairs de meme ord stra mais station differentes
              mutate(mean.width_river_cat = mean.width_river_cat1) %>%
                mutate(metric = y)
})
          })) %>%
      select(basin, beta_by_width) %>%
      unnest(beta_by_width) %>%
      mutate(metric = str_replace_all(metric, c("beta.jtu" = "turnover", "beta.jne" = "nestedness", "beta.jac" = "jaccard")))

    ),
  # Betadiv site versus all other site 
  tar_target(betadiv_site_by_basin_lg,
    betadiv_site_by_basin %>%
      left_join(st_basin_width_ord_sta, by = c("station", "basin")) %>%
      pivot_longer(c(jaccard, turnover, nestedness), names_to = "betadiv", values_to = "value")
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
    ),
  tar_target(riveratlas_shp_files,
    here::here("data-raw", "RiverATLAS_v10_shp", "RiverATLAS_v10_eu.shp"),
    format = "file"),
  tar_target(riveratlas_station, {

    to_select <- setNames(
      unique(c("HYRIV_ID", "LENGTH_KM", "DIST_UP_KM", "ORD_STRA", "ORD_FLOW",
          get_river_atlas_significant_var(), get_land_class_var())),
      NULL)

    layer_name <- sf::st_layers(riveratlas_shp_files, do_count = TRUE)$name
    myquery <- paste0(
      "SELECT ",
      paste0(to_select[!to_select %in% c("length_km", "dist_up_km", "ord_stra", "ord_flow", "hft_ix_c9309_ratio", "hft_ix_c9309_log2_ratio")], collapse = ", "),
      " FROM ", layer_name, " WHERE HYRIV_ID IN ",
      "(",
      paste0(na.omit(snapped_site_river$riverid), collapse = ","),
      ")"
    )
    sf::read_sf(dsn = riveratlas_shp_files, query = myquery) %>%
      mutate(riverid = HYRIV_ID) %>%
      janitor::clean_names() %>%
      left_join(select(st_drop_geometry(snapped_site_river), station, riverid), by = "riverid")
    }),
  tar_target(basin_dce,{
    load("~/Documents/post-these/mnhn/fishcom/data/the_8_hydrologic_basin.rda")
    the_8_hydrologic_basin
    }),
  tar_target(station_basin_dce, {
    intersect_st_basin <- st_intersects(op_st_filtered$station %>%
      st_transform(crs = 2154),
    basin_dce)
    tibble(
      station = op_st_filtered$station$station, 
      basin = basin_dce$NomDistric[map_int(intersect_st_basin, ~.x[1])]
    )
    }),
  tar_target(st_basin_width_ord_sta,
    op_st_filtered$op_desc %>%
      left_join(select(op_st_filtered$op, opcod, station), by = "opcod") %>%
      group_by(station) %>%
      summarise(
        mean.width_river = mean(width_river)
        ) %>%
      left_join(station_basin_dce, by = "station") %>%
      left_join(riveratlas_station %>%
        select(station, ord_stra) %>%
        st_drop_geometry(), by = "station") %>%
      mutate(station = as.character(station)) %>%
      mutate(
        mean.width_river_cat = cut(mean.width_river,
          c(0, 1, 2.5, 5, 10, 50, 100, 200, 400, 600)
          ),
        log10_mean.width_river = log10(mean.width_river),
        log10_mean.width_river_cat = cut(log10_mean.width_river,
          c(-1, 0, 1, 2, 3))
      )
    ),

  # report
  #tar_render(talk, "doc/slides.Rmd"),
  tar_render(betadiv_report, "doc/betadiv.Rmd"),


  #export
  tar_target(export, {
    save(fish_length, metaweb, network, com_species, com_metrics, com_trophic_species, network_mat, network_metrics, op_st_filtered, riveratlas_station,
      file = "~/Téléchargements/web_in_web.rda"
      )
    save(station_basin_dce, file = "~/Téléchargements/station_basin_dce.rda")
    save(betadiv_site_by_basin, file = "~/Téléchargements/betadiv_site_by_basin.rda")
    return("lololo") # trick for a happy targets
    }),

  # Plots

  ## Betadiv
  tar_target(p_beta_site_versus_all_ord_stra,
    betadiv_site_by_basin_lg %>%
      ggplot(aes(y = value, x = ord_stra, color = basin)) +
      geom_point() +
      geom_jitter() +
      geom_smooth(method = "lm", formula = as.formula(y ~ x + poly(x, 2))) +
      facet_grid(cols = vars(betadiv)) +
      theme(legend.position = "bottom")
    ),
  tar_target(p_beta_site_versus_all_width,
    betadiv_site_by_basin_lg %>%
      ggplot(aes(y = value, x = log(mean.width_river), color = basin)) +
      geom_point() +
      geom_jitter() +
      geom_smooth(method = "lm", formula = as.formula(y ~ x + poly(x, 2))) +
      facet_grid(cols = vars(betadiv)) +
      theme(legend.position = "bottom")
    ),
  tar_target(p_beta_site_by_stra,
    beta_by_stra %>%
      sample_n(2556) %>%
      ggplot(aes(y = beta, x = ord_stra, color = basin)) +
      geom_point() +
      geom_jitter() +
      geom_smooth(method = "lm", formula = as.formula(y ~ x + poly(x, 2))) +
      facet_grid(cols = vars(metric)) +
      theme(legend.position = "bottom")
    ),
  tar_target(p_beta_by_width,
    beta_by_width %>%
      sample_n(2556) %>%
      ggplot(aes(y = beta, x = as.numeric(mean.width_river_cat), color = basin)) +
      geom_point() +
      geom_jitter() +
      geom_smooth(method = "lm", formula = as.formula(y ~ x + poly(x, 2))) +
      facet_grid(cols = vars(metric)) +
      theme(legend.position = "bottom") +
      scale_x_continuous(
        name = "River width category (m)",
        breaks = seq_len(length(levels(beta_by_width$mean.width_river_cat))),
        labels = levels(beta_by_width$mean.width_river_cat)
      )
    )

)
