clean <- function(x) {
  as.numeric(stringr::str_replace_all(stringr::str_extract(x, '[0-9,]*'),
                                      pattern = ',',
                                      replacement = ''))
}
