#' Data for performing a mini-benchmark
#'
#' A dataset containing three triplets trees, trait-values and models to
#' evaluate the likelihood calculation times for R and C++ implementations. 
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{X}{trait values}
#'   \item{tree}{phylogenetic tree (phylo) with set edge.regimes member}
#'   \item{model}{MGPM model used to simulate the data in X}
#'   \item{ll}{log-likelihood value}
#'   \item{modelBM}{a random BM model}
#'   \item{llBM}{log-likelihood value form modelBM}
#'   \item{modelBM}{a random OU model}
#'   \item{llOU}{log-likelihood value for modelOU}
#' }
"benchmarkData"

#' Evaluate the likelihood calculation times for example trees and data
#' @param data a `data.frame` with at least the following columns: 
#' \itemize{
#' \item{tree: }{a list column of phylo objects with an edge.part member set.}
#' \item{X: }{a list column of k x N numerical matrices.}
#' \item{model: }{a list column of PCM objects.}
#' }
#' Defaults: to `benchmarkData`, which is small data.table included
#' with the PCMBaseCpp package.
#' @param nRepsCpp : number of repetitions for the cpp likelihood calculation 
#' calls: a bigger value increases the precision of time estimation at the 
#' expense of longer running time for the benchmark. Defaults to 10.
#' @param listOptions options to set before measuring the calculation times. 
#' Defaults to `list(PCMBase.Lmr.mode = 11, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9)`. 
#' `PCMBase.Lmr.mode` corresponds to the parallel traversal mode for the tree 
#' traversal algorithm (see 
#' \href{https://venelin.github.io/SPLITT/articles/SPLITTTraversalModes.html}{this page}
#' for possible values).
#' @return a data.frame.
#' @importFrom  PCMBase PCMInfo PCMLik PCMOptions MGPMDefaultModelTypes PCMTreeNumTips PCMTreeNumParts
#' @examples
#' \dontrun{
#' library(PCMBase)
#' library(PCMBaseCpp)
#' 
#' # original MGPM model
#' MiniBenchmarkRvsCpp()
#' 
#' # original MGPM model and parallel mode
#' MiniBenchmarkRvsCpp(
#' listOptions = list(PCMBase.Lmr.mode = 21, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9))
#' 
#' # single-trait data, original MGPM model and single mode and enabled option PCMBase.Use1DClasses
#' MiniBenchmarkRvsCpp(
#' data = PCMBaseCpp::benchmarkData[, list(
#'  tree, 
#'  X = lapply(X, function(x) x[1,, drop=FALSE]), 
#'  model = lapply(model, function(m) PCMExtractDimensions(m, dims = 1)))],
#' listOptions = list(PCMBase.Lmr.mode = 11, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9,
#' PCMBase.Use1DClasses = FALSE))
#' 
#' # random BM (non-MixedGaussian) model and parallel mode
#' MiniBenchmarkRvsCpp(
#' data = PCMBaseCpp::benchmarkData[, list(
#'  tree, 
#'  X,
#'  model = modelBM)],
#' listOptions = list(PCMBase.Lmr.mode = 21, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9))
#' 
#' # random OU (non-MixedGaussian) model and parallel mode
#' MiniBenchmarkRvsCpp(
#' data = PCMBaseCpp::benchmarkData[, list(
#'  tree, 
#'  X,, 
#'  model = modelOU)],
#' listOptions = list(PCMBase.Lmr.mode = 21, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9))
#'
#' # test on sinlge-trait data, BM (non-MGPM) model.
#' MiniBenchmarkRvsCpp(
#'  data = PCMBaseCpp::benchmarkData[, list(
#'  tree, 
#'  X = lapply(X, function(x) x[1,, drop=FALSE]), 
#'  model = lapply(modelBM, function(m) PCMExtractDimensions(m, dims = 1)))],
#'  listOptions = list(PCMBase.Lmr.mode = 11, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9, 
#'  PCMBase.Use1DClasses = FALSE))
#'  
#'  
#' # test on sinlge-trait data, OU (non-MGPM) model.
#' MiniBenchmarkRvsCpp(
#'  data = PCMBaseCpp::benchmarkData[, list(
#'  tree, 
#'  X = lapply(X, function(x) x[1,, drop=FALSE]), 
#'  model = lapply(modelOU, function(m) PCMExtractDimensions(m, dims = 1)))],
#'  listOptions = list(PCMBase.Lmr.mode = 11, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9))
#'  
#' # test on sinlge-trait data, OU (non-MGPM) model with enabled option PCMBase.Use1DClasses = TRUE
#' MiniBenchmarkRvsCpp(
#'  data = PCMBaseCpp::benchmarkData[, list(
#'  tree, 
#'  X = lapply(X, function(x) x[1,, drop=FALSE]), 
#'  model = lapply(model, function(m) PCMExtractDimensions(m, dims = 1)))],
#'  listOptions = list(PCMBase.Lmr.mode = 11, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9,
#'    PCMBase.Use1DClasses = FALSE))
#'
#' # test on sinlge-trait data, OU (non-MGPM) model with enabled option PCMBase.Use1DClasses = TRUE
#' MiniBenchmarkRvsCpp(
#'  data = PCMBaseCpp::benchmarkData[, list(
#'  tree, 
#'  X = lapply(X, function(x) x[1:2,, drop=FALSE]), 
#'  model = lapply(modelOU, function(m) PCMExtractDimensions(m, dims = 1:2)))],
#'  listOptions = list(PCMBase.Lmr.mode = 11, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9,
#'    PCMBase.Use1DClasses = TRUE))
#' 
#' MiniBenchmarkRvsCpp(
#'  data = PCMBaseCpp::benchmarkData[, list(
#'  tree, 
#'  X = lapply(X, function(x) x[1,, drop=FALSE]), 
#'  model = lapply(model, function(m) PCMExtractDimensions(m, dims = 1)))],
#'  listOptions = list(PCMBase.Lmr.mode = 11, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9))
#'
#' #' # test on sinlge-trait data using the PCMBase.Use1DClasses = TRUE option
#' MiniBenchmarkRvsCpp(
#'  data = PCMBaseCpp::benchmarkData[, list(
#'  tree, 
#'  X = lapply(X, function(x) x[1,, drop=FALSE]), 
#'  model = lapply(model, function(m) PCMExtractDimensions(m, dims = 1)))],
#'  listOptions = list(
#'    PCMBase.Lmr.mode = 11, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9,
#'    PCMBase.Use1DClasses = 1))
#'
#' # test on sinlge-trait data in parallel mode
#' MiniBenchmarkRvsCpp(
#'  data = PCMBaseCpp::benchmarkData[, list(
#'  tree, 
#'  X = lapply(X, function(x) x[1,, drop=FALSE]), 
#'  model = lapply(model, function(m) PCMExtractDimensions(m, dims = 1)))],
#'  listOptions = list(PCMBase.Lmr.mode = 21, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9))
#' 
#' MiniBenchmarkRvsCpp(listOptions = list(PCMBase.Lmr.mode = 21))
#' }
#' @export
MiniBenchmarkRvsCpp <- function(
  data = PCMBaseCpp::benchmarkData, nRepsCpp = 10L, 
  listOptions = list(PCMBase.Lmr.mode = 11, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9)) {
  
  listCurrentOptions <- options()
  
  do.call(options, listOptions)
  
  modelTypes <- MGPMDefaultModelTypes()
  
  
  res <- do.call(rbind, lapply(seq_len(nrow(data)), function(i) {
    tree <- data$tree[[i]]
    X <- data$X[[i]]
    model <- data$model[[i]]
    
    metaIR <- PCMInfo(X, tree, model)
    metaICpp <- PCMInfoCpp(X, tree, model, metaI = metaIR)
    
    valueR <- valueCpp <- as.double(NA)
    
    timeR <- system.time({
      valueR <- PCMLik(X, tree, model, metaI = metaIR)
    })[3]
    
    timeCpp <- system.time(
      valueCpp <-replicate(
        nRepsCpp, 
        PCMLik(X, tree, model, metaI = metaICpp))[1])[3] / nRepsCpp
    
    data.frame(
      N = PCMTreeNumTips(tree), 
      R = PCMTreeNumParts(tree),
      mapping = I(list(names(attr(model, "modelTypes"))[attr(model, "mapping")])), 
      type = I(list(class(model))),
      PCMBase.Lmr.mode = PCMOptions()$PCMBase.Lmr.mode,
      logLik = valueR, 
      logLikCpp = valueCpp, 
      timeR = unname(timeR), 
      timeCpp = unname(timeCpp))
  }))
  do.call(options, listCurrentOptions)
  res
}

#' A log-likelihood calculation time comparison for different numbers of traits and option-sets
#' @param ks a vector of positive integers, denoting different numbers of traits. 
#' Default: \code{c(1, 2, 4, 8)}.
#' @param optionSets a named list of lists of PCM-options. If NULL (the default) the 
#' option set is set to \code{
#' list(
#' `serial / 1D-multiv.` = list(PCMBase.Lmr.mode = 11, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9, PCMBase.Use1DClasses = FALSE),
#' `parallel / 1D-multiv.` = list(PCMBase.Lmr.mode = 21, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9, PCMBase.Use1DClasses = FALSE),
#' `serial / 1D-univar.` = list(PCMBase.Lmr.mode = 11, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9, PCMBase.Use1DClasses = TRUE),
#' `parallel / 1D-univar.` = list(PCMBase.Lmr.mode = 21, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9, PCMBase.Use1DClasses = TRUE)
#' )} for k = 1 and to \code{
#' list(
#' `serial / 1D-multiv.` = list(PCMBase.Lmr.mode = 11, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9, PCMBase.Use1DClasses = FALSE),
#' `parallel / 1D-multiv.` = list(PCMBase.Lmr.mode = 21, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9, PCMBase.Use1DClasses = FALSE))} for k > 1.
#' @param verbose logical indicating if log-messages should be printed to the console during the benchmark. Default FALSE.
#' @return a data.table for results similar to the data.table returned from \code{\link{MiniBenchmarkRvsCpp}} with 
#' additional columns for k, option-set and the type of model. 
#' @export
#' @importFrom data.table rbindlist data.table
#' @importFrom PCMBase PCMExtractDimensions
BenchmarkRvsCpp <- function(
  ks = c(1, 2, 4, 8),
  optionSets = NULL,
  verbose = FALSE) {

  findBiggestFactor <- function(k, fMax) {
    f <- fMax
    while(f >= 1) {
      if(k %% f == 0) {
        break
      }
      f <- f - 1
    }
    f
  }
  
  resultList <- lapply(ks, function(k) {
    if(is.null(optionSets)) {
      optionSets <- if(k == 1) {
        list(
          `serial / 1D-multiv.` = list(PCMBase.Lmr.mode = 11, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9, PCMBase.Use1DClasses = FALSE),
          `parallel / 1D-multiv.` = list(PCMBase.Lmr.mode = 21, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9, PCMBase.Use1DClasses = FALSE),
          `serial / 1D-univar.` = list(PCMBase.Lmr.mode = 11, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9, PCMBase.Use1DClasses = TRUE),
          `parallel / 1D-univar.` = list(PCMBase.Lmr.mode = 21, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9, PCMBase.Use1DClasses = TRUE))
      } else {
        list(
          `serial / 1D-multiv.` = list(PCMBase.Lmr.mode = 11, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9, PCMBase.Use1DClasses = FALSE),
          `parallel / 1D-multiv.` = list(PCMBase.Lmr.mode = 21, PCMBase.Threshold.EV = 1e-9, PCMBase.Threshold.SV = 1e-9, PCMBase.Use1DClasses = FALSE))
      }
    } 
    
    
    resultList <- lapply(names(optionSets), function(oset) {
      
      if(verbose) {
        cat("Performing benchmark for k: ", k, "; optionSet: ", oset, "...\n")
      }
      fk <- findBiggestFactor(k, 8) 
      ds <- seq_len(fk)
      nRB <- k / fk
      
      testData <- rbindlist(list(
        benchmarkData[, list(
          k = k,
          modelType = "MGPM (A-F)",
          options = oset,
          tree, 
          X = lapply(X, function(x) x[rep(ds, nRB),, drop=FALSE]), 
          model = lapply(model, function(m) PCMExtractDimensions(m, dims = ds, nRepBlocks = nRB)))],
        benchmarkData[, list(
          k = k,
          modelType = "BM (B)",
          options = oset,
          tree, 
          X = lapply(X, function(x) x[rep(ds, nRB),, drop=FALSE]), 
          model = lapply(modelBM, function(m) PCMExtractDimensions(m, dims = ds, nRepBlocks = nRB)))],
        benchmarkData[, list(
          k = k,
          modelType = "OU (E)",
          options = oset,
          tree, 
          X = lapply(X, function(x) x[rep(ds, nRB),, drop=FALSE]), 
          model = lapply(modelOU, function(m) PCMExtractDimensions(m, dims = ds, nRepBlocks = nRB)))]))
      
      resultData <- MiniBenchmarkRvsCpp(data = testData, listOptions = optionSets[[oset]])
      resultData <- cbind(testData[, list(k, modelType, options)], resultData)
      if(verbose) {
        print(resultData[, list(
          k, modelType, options, N, R, mapping, mode = PCMBase.Lmr.mode, logLik, 
          logLikCpp, timeR, timeCpp)])
      }
      resultData
    })
    rbindlist(resultList, use.names = TRUE)
  })
  rbindlist(resultList, use.names = TRUE)
}
