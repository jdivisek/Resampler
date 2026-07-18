#' @title Distance and Similarity-Based Resampling of Vegetation Plots
#'
#' @description This function resamples vegetation plots (relevés) based on a
#' combination of geographic proximity and species composition similarity.
#' The goal is to reduce sampling density in areas where many plots are located
#' close to each other and record very similar or identical vegetation.
#'
#' @details
#' The function operates in the following steps:
#'
#' ### 1. Data Preparation and Validation
#' A series of checks is performed to verify that the input data meet all requirements regarding format, type, and the absence of NA values.
#'
#' ### 2. Data Preparation and Sorting
#'    * First, the \code{coord} table is sorted according to the rule defined in \code{remove}. For "random", the rows are randomly shuffled. For other methods, they are first randomly shuffled and then sorted by diversity or the value in \code{var.value} This serves as a universal tie-breaking rule in subsequent steps.
#'    * Based on this final order, a new, internal numeric identifier id (1, 2, 3...) is created and joined to the \code{spec} table. All further operations work with this id.
#'
#' ### 3. Neighbor Identification
#' A list of neighbors within the \code{dist.threshold} is found for each plot. If stratification is used, neighbors are only searched for within the same stratum.
#'
#' ### 4. Splitting into Groups
#' All plots are divided into independent, geographically contiguous groups
#'
#' ### 5. Resampling
#' The actual resampling process is performed for each group individually.
#'    * For the given group, similarity between neighboring plots is calculated using \code{sim.method} and a list of all unique pairs of plots that are both geographically close (i.e. neighboring) and exceed \code{sim.threshold} is created.
#'    * This list of "conflicting" pairs is sorted in descending order based on their compositional similarity.
#'    * The function then iterates through this sorted list, starting with the most similar pair. For the first pair in the list, it decides which plot to remove based on \code{remove} rule and adds it to a "blacklist". Subsequently, it removes ALL pairs from the list that contained this just-removed plot. The process is repeated on the reduced list until no conflicting pairs remain.
#'
#' ### 6. Result
#' The blacklists of removed plots from all groups are combined. Based on this final blacklist, the original \code{coord} table is filtered, and the thinned and cleaned table is returned.
#'
#' **Performance notes**
#'
#' The function can process very large datasets (lower hundreds of thousands of plots) but its actual performance critically depends on parameters set:
#'    * \code{dist.threshold} value is the most important. The larger the value, the larger the geographically contiguous groups of plots are processed. For example, when setting \code{dist.threshold = 1000} for grassland vegetation plots from the European Vegetation Archive, the largest group contained more than 29,000 plots! Depending on the actual size of the dataset, it is therefore recommended to set \code{dist.threshold} value no higher than 5,000 m. For smaller datasets, large distances can be set, but they are not ecologically very meaningful. For large datasets, high distance values increase processing time and can cause memory issues.
#'    * \code{sim.threshold} value is another important parameter that influences the performance of the function, but it is not as critical as \code{dist.threshold}. Setting a very low similarity value will result in fewer preserved plots, and vice versa.
#'    * \code{longlat} Although the function can handle geographical coordinates in degrees, it is highly recommended to provide coordinates in a projected coordinate system such as ETRS89 and set \code{longlat = FALSE}. In this case, the function uses Euclidean distance instead of Great Circle distance, which speeds up the identification of neighboring plots.
#'
#' The function was tested with a dataset containing 468,341 grassland vegetation
#' plots from the European Vegetation Archive using the following settings:
#' \code{longlat = FALSE}, \code{dist.threshold = 1000}, \code{sim.threshold = 0.8},
#' \code{sim.method = "simpson"}, \code{remove = "random"}, and \code{strata = NULL}.
#' On an older PC with 8 GB RAM and an Intel Core i5-9400F 2.8 GHz processor,
#' resampling took 20 hours and 53 minutes without any memory issues. Therefore,
#' you should be patient;-). The resampling removed 30.9% of plots (144,860 out of 468,341).
#'
#' The function was also tested with a smaller dataset of 114,854 plots from
#' the Czech Vegetation Database using the following settings: \code{longlat = FALSE},
#' \code{dist.threshold = 1000}, \code{sim.threshold = 0.5}, \code{sim.method = "simpson"},
#' \code{remove = "random"}, and \code{strata = NULL}. On a laptop with 32 GB RAM
#' and an Intel Core i7-11850H 2.5 GHz processor, the resampling procedure took 8 minutes.
#' When using environmental strata, the resampling took less than 3 minutes.
#'
#' @section Warning:
#' \code{coord} and \code{spec} tables must have the structure and column names described in the Arguments section!
#'
#' @param coord \code{data.table} with the following columns:
#'    * `PlotObservationID` A unique identifier for each plot (numeric, integer or character). Each PlotObservationID must be unique.
#'    * X-coordinate or Longitude. Must be of type numeric and must not contain any NA values.
#'    * Y-coordinate or Latitude. Must be of type numeric and must not contain any NA values.
#'    * Other optional column(s) containing, for example, environmental strata and/or a variable used for selecting vegetation plots (see \code{var.value} parameter)
#' @param spec \code{data.table} containing vegetation-composition data in "long format". Required columns are:
#'    * `PlotObservationID` The plot identifier, which corresponds to the IDs in \code{coord}.
#'    * `Taxon_name` The name of the recorded taxon (species).
#'    * `cover` The species cover value (usually in %). Must be of type numeric and must not contain any NA and zero values.
#' @param longlat A logical value (TRUE/FALSE). If TRUE, coordinates are treated
#' as latitude/longitude, and distances are calculated in kilometers. If FALSE
#' (default), a projected coordinate system (e.g., UTM) is assumed, and
#' distances are in meters.
#' @param dist.threshold A numeric value. The threshold for geographic distance.
#' From a pair of plots closer than this value, one will be removed (if they
#' also meet the \code{sim.threshold}). Units (meters/kilometers) depend on the
#' \code{longlat} parameter. Default is 1000.
#' @param sim.threshold A numeric value (0-1). The threshold for species
#' composition similarity. From a pair of plots more similar than this value,
#' one will be removed (if they also meet the \code{dist.threshold}). Default is 0.8.
#' @param sim.method The method for calculating similarity. Options: "simpson"
#' (Simpson), "sorensen" (Sørensen), "jaccard" (Jaccard) and "bray"
#' (Bray-Curtis). Default is "simpson". For all methods except "bray", cover data is automatically
#' converted to presence/absence. If only "presences" are provided (i.e. all
#' cover values are 1) and \code{sim.method = "bray"}, Sørensen index is calculated.
#' @param remove The rule that decides which plot from a conflicting pair will be removed. Default is "random".
#'    * `random` Randomly removes one of the two plots.
#'    * `less diverse` Removes the plot with the lower number of richness. Ties are broken by "random" order.
#'    * `more diverse` Removes the plot with the higher number of species. Ties are broken by "random" order.
#'    * `lower var.value` Removes the plot with the lower value in the column defined by the var.value parameter. NAs are allowed and plots with NA are removed first.
#'    * `higher var.value` Removes the plot with the higher value in the column defined by the var.value parameter. NAs are allowed and plots with NA are removed first.
#' @param var.value A character string. The name of a column in \code{coord} used for
#' decision-making with the "lower var.value" and "higher var.value" methods.
#' NAs are allowed in this variable and plots with NA are removed first.
#' @param strata A character string. The name of a column in \code{coord} that defines
#' plot stratification. If provided, resampling is performed separately within
#' each stratum (group).
#' @param seed A number. The seed value for the random number generator, which
#' ensures reproducibility of results. Default is 1234.
#'
#' @return Resampled \code{data.table} with PlotObservationIDs and geographical coordinates for selected plots.
#'
#' @author Jan Divíšek
#'
#' @seealso
#' To see a practical example of resampling vegetation plots from the Czech Vegetation
#' Database, visit the database website at \url{https://czechveg.github.io/DataProcessingTutorial/data_resampling.html}.
#'
#' @export
resample_plots <- function(coord, spec, longlat = FALSE, dist.threshold = 1000, sim.threshold = 0.8,
                           sim.method = c("simpson", "sorensen", "jaccard", "bray"),
                           remove = c("random", "less diverse", "more diverse", "lower var.value", "higher var.value"),
                           var.value = NULL, strata = NULL, seed = 1234) {

  start_time <- Sys.time()

  # --- 1. Check data ---
  sim.method <- match.arg(sim.method); remove <- match.arg(remove)
  if(!is.data.table(coord)) stop("Error: 'coord' must be a data.table")
  if(!is.data.table(spec)) stop("Error: 'spec' must be a data.table")

  if(colnames(coord)[1] != "PlotObservationID") stop("Error: First column in 'coord' must be 'PlotObservationID'")
  if(ncol(coord) < 3) stop("Error: 'coord' must contain at least three columns (PlotObservationID, X coordinate, Y coordinate)")
  if(!all(sapply(coord[, 2:3], is.numeric))) stop("Error: Coordinates (X, Y or Lon, Lat) must be numeric")
  if(any(is.na(coord[, 2:3]))) stop("Error: Coordinates contain NA values")
  if(remove %chin% c("lower var.value", "higher var.value") && (is.null(var.value) || !var.value %chin% colnames(coord))) stop("Error: Invalid or missing var.value")
  if(!is.null(strata) && !strata %chin% colnames(coord)) stop("Error: Invalid strata column name")
  if(!all(c("PlotObservationID", "Taxon_name", "cover") %chin% colnames(spec))) stop("Error: Check column names in 'spec'")
  if(any(spec$cover == 0)) stop("Error: Zero covers not allowed")
  if(any(is.na(spec$cover))) stop("Error: NAs in cover")
  if(!all(spec$PlotObservationID %in% coord$PlotObservationID)) stop("Error: Some PlotObservationID in 'spec' not found in 'coord'")
  if(!all(coord$PlotObservationID %in% spec$PlotObservationID)) stop("Error: Some PlotObservationID in 'coord' not found in 'spec'")

  # --- 2. Prepare data ---
  set.seed(seed)
  coord <- coord[sample(1:nrow(coord)), ]

  if(remove == "less diverse"){
    coord <- coord[spec[, .(.N), by = .(PlotObservationID)], on = c("PlotObservationID")][order(N)]}
  if(remove == "more diverse"){
    coord <- coord[spec[, .(.N), by = .(PlotObservationID)], on = c("PlotObservationID")][order(-N)]}
  if(remove == "lower var.value"){
    coord <- setorderv(coord, var.value, 1, na.last=FALSE)}
  if(remove == "higher var.value"){
    coord <- setorderv(coord, var.value, -1, na.last=FALSE)}

  if (sim.method != "bray") { spec[, cover := 1] }

  # --- 3. Indentification of neighbouring plots ---
  if (!is.null(strata)) {
    cat(paste("Searching for neighbouring plots within", uniqueN(coord[[strata]]), "strata:\n"))
    coord[, (strata) := as.character(.SD[[1]]), .SDcols = strata]

    coord <- setorderv(coord, strata, 1, na.last=FALSE)

    pb <- txtProgressBar(min = 0, max = uniqueN(coord[[strata]]), style = 3)
    d <- coord[, {setTxtProgressBar(pb, .GRP);
      data.table(NB = dnearneigh_strat(x = as.matrix(.SD), row.names = .I,
                                       d1 = 0, d2 = dist.threshold, longlat = longlat))},
      by = strata, .SDcols = 2:3]
    close(pb)

    d <- d$NB
  } else {
    cat("Searching for neighbouring plots. Please wait...\n")
    d <- spdep::dnearneigh(as.matrix(coord[, 2:3]), d1 = 0, d2 = dist.threshold, bounds = c("GE", "LT"), longlat = longlat)
  }

  ##set new ids based on ordered coord
  spec[coord[, .(.I, PlotObservationID)], on = "PlotObservationID", id := i.I]

  # --- 4. Split plots to groups ---
  g <- igraph::components(igraph::graph_from_adj_list(d, mode = "all"))

  spec[data.table(id = coord[,.I], group = g$membership), on = "id", grp := group]
  spec <- spec[grp %in% which(g$csize > 1), ]

  d <- mapply(append, seq_along(d), d, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  names(d) <- coord[,.I]
  d <- d[g$membership %in% which(g$csize > 1)]

  cat("Similarity-based resampling:\n")
  pb <- txtProgressBar(min = 0, max = uniqueN(spec$grp), style = 3)
  blacklist <-  spec[, {setTxtProgressBar(pb, .GRP);
    filtering_task(.SD, d, sim.threshold, sim.method)},
    by = .(grp),
    .SDcols = PlotObservationID:grp]

  close(pb)

  # --- 5. Return selected plots ---

  coord.filtered <- coord[-blacklist$V1, 1:3]

  cat(sprintf("Removed %.1f%% of plots (%d out of %d)\n", nrow(blacklist) / nrow(coord) * 100, nrow(blacklist), nrow(coord)))
  end_time <- Sys.time(); elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf("Elapsed time %02d:%02d:%02d\n", floor(elapsed / 3600), floor((elapsed %% 3600) / 60), round(elapsed %% 60)))

  return(coord.filtered)
}
