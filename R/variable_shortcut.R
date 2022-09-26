
get_river_atlas_significant_var <- function() {
  c(
    "Length of the reach" = "length_km",
    "Distance from source (km)" = "dist_up_km",
    "Strahler order" = "ord_stra",
    "Order flow" = "ord_flow",
    "Average elevation (m)" = "ele_mt_cav",
    "Annual average of temperature (°C)" = "tmp_dc_cyr",
    "Annual maximum of temperature (°C)" = "tmp_dc_cmx",
    "Annual minimum of temperature (°C)" = "tmp_dc_cmn",
    "Average slope (degree)" = "slp_dg_cav",
    "Annual average of discharge (m3/s)" = "dis_m3_pyr",
    "Annual maximum of discharge (m3/s)" = "dis_m3_pmx",
    "Annual minimum of discharge (m3/s)" = "dis_m3_pmn",
    "River area (reach segment, in ha)" = "ria_ha_csu",
    "River volume (reach segment, m3)" = "riv_tc_csu",
    "Protected area extent (%)" = "pac_pc_cse",
    "Urban extent (%)" = "urb_pc_cse",
    "Human footprint ratio 2009/1993" = "hft_ix_c9309_ratio",
    "Log2 Human footprint ratio (2009/1993)" = "hft_ix_c9309_log2_ratio",
    "Human footprint index 1993" = "hft_ix_c93",
    "Human footprint index 2009" = "hft_ix_c09"
  )
}

get_land_class_var <- function() {
  c(
    #  "Potential Natural Vegetation Extent (catchment)" = "pnv_pc_c01-c15",
    #  "Potential Natural Vegetation Extent (watershed)" = "pnv_pc_u01-u15",
    "Forest Cover Extent (%, catchment)" = "for_pc_cse",
    #"Forest Cover Extent (%, watershed)" = "for_pc_use",
    "Cropland Extent (%, catchment)" = "crp_pc_cse",
    #"Cropland Extent (%, watershed)" = "crp_pc_use",
    "Pasture Extent (%, catchment)" = "pst_pc_cse",
    #"Pasture Extent (%, watershed)" = "pst_pc_use",
    "Urban Extent (%, catchment)" = "urb_pc_cse",
    #"Urban Extent (%, watershed)" = "urb_pc_use",
    #"Population Density (watershed)" = "ppd_pk_uav",
    "Population Density (people per km2, catchment)" = "ppd_pk_cav"
  )
}
