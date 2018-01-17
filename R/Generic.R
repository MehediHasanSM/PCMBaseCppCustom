#' A generic method for creating C++ backend objects given a model, data and 
#' a tree.
#' @useDynLib PCMBaseCpp
#' @importFrom PCMBase validateModel
#' @export
newCppObject <- function(X, tree, model, metaInfo = validateModel(tree, model), ...) {
  UseMethod("newCppObject", model)
}

# loading the C++ modules
loadModule( "QuadraticPolynomialBM", TRUE )
loadModule( "QuadraticPolynomialOU", TRUE )
loadModule( "QuadraticPolynomialJOU", TRUE )
loadModule( "QuadraticPolynomialTwoSpeedOU", TRUE )
