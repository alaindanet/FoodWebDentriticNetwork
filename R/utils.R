return_data <- function(file_path = NULL) {
  load(file_path)
  return(get(ls()[! ls() %in% "file_path"]))
}
