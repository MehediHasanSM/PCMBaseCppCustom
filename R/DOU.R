#'@export
PCMInfoCpp.DOU <- function(X, tree, model, metaI = PCMInfo(X, tree, model, verbose), verbose = FALSE, ...) {
  
  res <- c(metaI, cppObject = QuadraticPolynomialDOU$new(X, tree, model, metaI))
  
  res$TraverseTree = res$cppObject$TraverseTree
  
  class(res) <- c("PCMInfoCpp", class(metaI))
  res
}

loadModule( "QuadraticPolynomialDOU", TRUE )
