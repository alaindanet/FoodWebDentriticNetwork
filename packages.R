## library() calls go here
library(conflicted)
library(dotenv)
library(targets)
library(tarchetypes)

library(tidyverse)
library(magrittr)
library(sf)

library(sizeTrophicInteractions)

library(furrr)

# Conflict
conflict_prefer("filter", "dplyr")
