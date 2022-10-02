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
tar_load(c(beta_by_stra, beta_by_width))
tar_load(c(p_beta_site_by_stra, p_beta_site_by_width))
beta_by_stra %>%
  group_by(basin) %>%
  sample_n(2556) %>%
  ggplot(aes(y = beta, x = ord_stra, color = basin)) +
  geom_point() +
  geom_jitter() +
  geom_smooth(method = "lm", formula = as.formula(y ~ x + poly(x, 2))) +
  facet_grid(cols = vars(metric)) +
  theme(legend.position = "bottom")

