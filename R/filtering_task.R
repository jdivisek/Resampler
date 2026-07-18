#' @noRd
filtering_task <- function(spec_sub, d, sim.threshold, sim.method){

  d_sub <- d[names(d) %chin% as.character(unique(spec_sub$id))]

  spec_sub <- Matrix::sparseMatrix(i = as.integer(factor(spec_sub$id), levels = names(d_sub)),
                                   j = as.integer(factor(spec_sub$Taxon_name)),
                                   x = spec_sub$cover,
                                   dimnames = list(names(d_sub), levels(factor(spec_sub$Taxon_name))))

  ds <- lapply(d_sub, FUN = calc_sim, spec_sub, sim.threshold, sim.method)

  pairs <- rbindlist(lapply(names(ds), function(p1) {
    p2 <- ds[[p1]][ds[[p1]] > sim.threshold & as.integer(names(ds[[p1]])) > as.integer(p1)]
    if (length(p2) == 0) return(NULL)
    data.table(p1 = p1, p2 = names(p2), sim = p2) }))

  if(nrow(pairs) > 0){
    pairs <- setorderv(pairs, c("sim"), c(-1))

    black <- NULL
    repeat
    {
      p <- pairs$p1[1]
      black <- c(black, as.integer(p))
      pairs <- pairs[!(pairs$p1 == p | pairs$p2 == p)]

      if(nrow(pairs) == 0){ break}
    }
    return(black)
  } else {
    return(NULL)
  }
}
