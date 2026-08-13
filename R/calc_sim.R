#' @noRd
calc_sim <- function(plots, spec_sub, sim.method){

  n <- length(plots)-1

  sim <- switch(sim.method,
                "simpson" = 1 - vegan::betadiver(as.matrix(spec_sub[as.character(plots), ]), method = "sim")[1:n],
                "sorensen" = 1 - vegan::vegdist(spec_sub[as.character(plots), ], method = "bray", binary = TRUE)[1:n],
                "bray" = 1 - vegan::vegdist(spec_sub[as.character(plots), ], method = "bray", binary = FALSE)[1:n],
                "jaccard" = 1 - vegan::vegdist(spec_sub[as.character(plots), ], method = "jaccard", binary = TRUE)[1:n])
  names(sim) <- plots[-1]
  return(sim)
}
