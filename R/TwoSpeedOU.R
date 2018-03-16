#'@export
PCMInfoCpp.TwoSpeedOU <- function(X, tree, model, metaI = PCMInfo(X, tree, model, verbose), verbose = FALSE, ...) {
  
  res <- c(metaI, cppObject = QuadraticPolynomialTwoSpeedOU$new(X, tree, metaI))
  
  res$TraverseTree = res$cppObject$TraverseTree
  
  class(res) <- c("PCMInfoCpp", class(metaI))
  res
}

loadModule( "QuadraticPolynomialTwoSpeedOU", TRUE )
