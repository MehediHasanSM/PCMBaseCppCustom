#'@export
PCMInfoCpp.BM <- function(X, tree, model, metaI = PCMInfo(X, tree, model, verbose), verbose = FALSE, ...) {
  
  res <- c(metaI, cppObject = QuadraticPolynomialBM$new(X, tree, metaI))
  res$TraverseTree = res$cppObject$TraverseTree
  
  class(res) <- c("PCMInfoCpp", class(metaI))
  res
}

loadModule( "QuadraticPolynomialBM", TRUE )
