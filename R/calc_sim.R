#' @noRd
calc_sim <- function(plots, spec_sub, sim.threshold, sim.method){

  n <- length(plots)-1

  sim <- switch(sim.method,
                "simpson" = 1 - vegan::betadiver(as.matrix(spec_sub[as.character(plots), ]), method = "sim")[1:n],
                "sorensen" = ,
                "bray" = 1 - vegan::vegdist(spec_sub[as.character(plots), ], method = "bray")[1:n],
                "jaccard" = 1 - vegan::vegdist(spec_sub[as.character(plots), ], method = "jaccard")[1:n])
  names(sim) <- plots[-1]
  return(sim)
}
