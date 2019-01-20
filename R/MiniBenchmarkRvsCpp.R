#' Evaluate the likelihood calculation times for example trees and data
#' @param data a `data.frame` with at least the following columns: 
#' \itemize{
#' \item{tree: }{a list column of phylo objects with an edge.regime member set.}
#' \item{X: }{a list column of k x N numerical matrices.}
#' \item{model: }{a list column of PCM objects.}
#' }
#' Defaults: to `miniBenchmarkData`, which is small data-table included
#' with the PCMBaseCpp package.
#' @param nRepsCpp : number of repetitions for the cpp likelihood calculation 
#' calls: a bigger value increases the precision of time estimation at the 
#' expense of longer running time for the benchmark. Defaults to 10.
#' @param listOptions options to set before measuring the calculation times. 
#' Defaults to `list(PCMBase.Lmr.mode = 11, PCMBase.Threshold.EV = 1e-7)`. 
#' `PCMBase.Lmr.mode` corresponds to the parallel traversal mode for the tree 
#' traversal algorithm (see 
#' \href{https://venelin.github.io/SPLITT/articles/SPLITTTraversalModes.html}{this page}
#' for possible values).
#' @return a data.frame.
#' @importFrom  PCMBase PCMInfo PCMLik PCMOptions MGPMDefaultModelTypes PCMTreeNumTips PCMTreeNumUniqueRegimes
#' @examples
#' \dontrun{
#' library(PCMBase)
#' library(PCMBaseCpp)
#' MiniBenchmarkRvsCpp()
#' 
#' MiniBenchmarkRvsCpp(listOptions = list(PCMBase.Lmr.mode = 21))
#' }
#' @export
MiniBenchmarkRvsCpp <- function(
  data = PCMBaseCpp::miniBenchmarkData, nRepsCpp = 10L, 
  listOptions = list(PCMBase.Lmr.mode = 11, PCMBase.Threshold.EV = 1e-7)) {
  
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
      R = PCMTreeNumUniqueRegimes(tree),
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