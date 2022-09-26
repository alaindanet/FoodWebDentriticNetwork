
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FoodWebDentriticNetwork

<!-- badges: start -->

<!-- badges: end -->

The goal of FoodWebDentriticNetwork is to …

## Explain ouput web\_in\_web.rda

  - `fish_length`: individual fish bodysize
  - `metaweb`: list of `metaweb`, `species` and `resources`
  - `network`: networks for each sampling events
  - `com_species`, `com_trophic_species`: metric at scale of species and
    trophic species (e.g. biomass and number of individuals)
  - `com_metrics`, `network_metrics`: metrics at community and network
    scale
  - `network_mat`
  - `op_st_filtered`: list of “station”, “op”, “op\_desc”, “op\_env”
    (site information and sampling information)
  - `riveratlas_station`: collection of hydro-morpho, climatic variables
    for RiverAtlas

# load data

``` r
load("web_in_web.rda")
```
