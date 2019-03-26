#' Data for performing a mini-benchmark
#'
#' A dataset containing ten triplets of trees, trait-values and models to
#' evaluate the likelihood calculation times for R and C++ implementations. 
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{X}{trait values}
#'   \item{tree}{phylogenetic tree (phylo) with set edge.regimes member}
#'   \item{model}{PCM model}
#' }
"miniBenchmarkData"
