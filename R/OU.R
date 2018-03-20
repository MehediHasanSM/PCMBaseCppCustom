#'@export
PCMInfoCpp.OU <- function(X, tree, model, metaI = PCMInfo(X, tree, model, verbose), verbose = FALSE, ...) {
  
  res <- c(metaI, cppObject = QuadraticPolynomialOU$new(X, tree, model, metaI))
  
  res$TraverseTree = res$cppObject$TraverseTree
  
  class(res) <- c("PCMInfoCpp", class(metaI))
  res
}

loadModule( "QuadraticPolynomialOU", TRUE )
