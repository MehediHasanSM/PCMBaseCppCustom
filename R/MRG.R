#'@export
PCMInfoCpp.MRG <- function(X, tree, model, metaI = PCMInfo(X, tree, model, verbose), verbose = FALSE, ...) {
  specParams <- attr(model, "specParams")
  regimeModel <- unlist(sapply(specParams, function(param) {
    if(param$type[1] == "model") {
      if(startsWith(class(param$default)[1], "BM")) {
        "BM"
      } else if(startsWith(class(param$default)[1], "JOU")) {
        "JOU"
      } else if(startsWith(class(param$default)[1], "TwoSpeedOU")) {
        "TwoSpeedOU"
      } else if(startsWith(class(param$default)[1], "OU")) {
        "OU"
      } else 
        stop(paste0("ERR:03501:PCMBaseCpp:MRG.R:PCMInfoCpp.MRG:: Uknown model class ", class(param$default)[1]))
    } else {
      NULL
    }
  }))
  retimeModel <- regimeModel[!sapply(regimeModel, is.null)]
  
  res <- c(metaI, 
           cppObject = QuadraticPolynomialMRG$new(
             X, tree, model, metaI, 
             regimeModel))
  
  res$TraverseTree = res$cppObject$TraverseTree
  
  class(res) <- c("PCMInfoCpp", class(metaI))
  res
}

loadModule( "QuadraticPolynomialMRG", TRUE )
