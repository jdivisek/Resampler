#' @noRd
dnearneigh_strat <- function(x, d1, d2, row.names, longlat){
  nb <- suppressWarnings(spdep::dnearneigh(x = x, d1 = d1, d2 = d2, row.names = row.names,
                                           bounds = c("GE", "LT"), longlat = longlat))
  lapply(nb, FUN = function(regions) {if(regions[1] > 0) {attr(nb, "region.id")[regions]} else {0}})
}
