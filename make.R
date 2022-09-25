library(targets)
library(tarchetypes)

tar_make()


tar_make_future(
  workers = min(future::availableCores() - 1, 24),
  names = !c(starts_with("beta_"), !!exclusion_vector)
  )

tar_make_future(
  workers = min(future::availableCores() - 1, 8),
  names = !c(starts_with("beta_"), !!exclusion_vector)
  )

tar_meta()
tar_visnetwork()

