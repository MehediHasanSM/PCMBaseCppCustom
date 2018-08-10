#'@export
PCMInfoCpp.White <- function(X, tree, model, 
                             #metaI = PCMInfo(X, tree, model, verbose), 
                             metaI = PCMInfo(X, tree, model, verbose, preorder=PCMTreePreorderCpp(tree)), 
                             verbose = FALSE, ...) {
  
  if(is.Transformable(model)) {
    model <- PCMApplyTransformation(model)
  }
  
  metaI$pcListInt <- PCListInt(metaI$pc)
  
  res <- c(metaI, cppObject = QuadraticPolynomialWhite$new(X, tree, model, metaI))
  res$TraverseTree = res$cppObject$TraverseTree
  
  class(res) <- c("PCMInfoCpp", class(metaI))
  res
}

loadModule( "QuadraticPolynomialWhite", TRUE )
