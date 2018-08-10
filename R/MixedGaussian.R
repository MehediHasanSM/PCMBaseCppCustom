#'@export
PCMInfoCpp.MixedGaussian <- function(X, tree, model,
                           metaI = PCMInfo(X, tree, model, verbose), 
                           #metaI = PCMInfo(X, tree, model, verbose, preorder=PCMTreePreorderCpp(tree)), 
                           verbose = FALSE, ...) {
  
  if(is.Transformable(model)) {
    model <- PCMApplyTransformation(model)
  }
  
  spec <- attr(model, "spec", exact = TRUE)
  regimeModel <- unlist(sapply(spec, function(param) {
    if(is.PCM(param)) {
      if(startsWith(class(param)[1], "BM")) {
        "BM"
      } else if(startsWith(class(param)[1], "JOU")) {
        "JOU"
      } else if(startsWith(class(param)[1], "DOU")) {
        "DOU"
      } else if(startsWith(class(param)[1], "OU")) {
        "OU"
      } else 
        stop(paste0("ERR:03501:PCMBaseCpp:MixedGaussian.R:PCMInfoCpp.MixedGaussian:: Uknown model class ", class(param)[1]))
    } else {
      NULL
    }
  }))
  regimeModel <- regimeModel[!sapply(regimeModel, is.null)]
  
  metaI$pcListInt <- PCListInt(metaI$pc)
  
  res <- c(metaI, 
           cppObject = QuadraticPolynomialMixedGaussian$new(
             X, tree, model, metaI, 
             regimeModel))
  
  res$TraverseTree = res$cppObject$TraverseTree
  
  class(res) <- c("PCMInfoCpp", class(metaI))
  res
}

loadModule( "QuadraticPolynomialMixedGaussian", TRUE )
