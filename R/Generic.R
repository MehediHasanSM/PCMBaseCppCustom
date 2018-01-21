#' A generic method for creating C++ backend objects given a model, data and 
#' a tree.
#' @useDynLib PCMBaseCpp
#' @importFrom PCMBase validateModel
#' @export
newCppObject <- function(X, tree, model, metaInfo = validateModel(tree, model), ...) {
  UseMethod("newCppObject", model)
}

extractAbCdEfLmr <- function(pruneI) {
  tr <- pruneI$tree
  M <- tr$num_nodes
  # internal node ids from 1 to M
  nodeIds <- sapply(1:M, function(n) tr$FindIdOfNode(n) + 1)
  
  spec <- pruneI$spec
  specList <- list(
    A = spec$A, b = spec$b, C = spec$C, d = spec$d, E = spec$E, f = spec$f, 
    L = spec$L, m = spec$m, r = spec$r, V = spec$V, V_1 = spec$V_1
  )
  with(specList,
       list(A = abind(lapply(nodeIds, function(i) A[,, i, drop=TRUE]), along = -1),
            b = abind(lapply(nodeIds, function(i) b[, i, drop=TRUE]), along = -1),
            C = abind(lapply(nodeIds, function(i) C[,, i, drop=TRUE]), along = -1),
            d = abind(lapply(nodeIds, function(i) d[, i, drop=TRUE]), along = -1),
            E = abind(lapply(nodeIds, function(i) E[,, i, drop=TRUE]), along = -1),
            f = as.vector(abind(lapply(nodeIds, function(i) f[i, drop=TRUE]), along = -1)),
            V = abind(lapply(nodeIds, function(i) V[,, i, drop=TRUE]), along = -1),
            V_1 = abind(lapply(nodeIds, function(i) V_1[,, i, drop=TRUE]), along = -1),
            L = abind(lapply(nodeIds, function(i) L[,, i, drop=TRUE]), along = -1),
            m = abind(lapply(nodeIds, function(i) m[, i, drop=TRUE]), along = -1),
            r = as.vector(abind(lapply(nodeIds, function(i) r[i, drop=TRUE]), along = -1))
            ))
}

# loading the C++ modules
loadModule( "QuadraticPolynomialBM", TRUE )
loadModule( "QuadraticPolynomialOU", TRUE )
loadModule( "QuadraticPolynomialJOU", TRUE )
loadModule( "QuadraticPolynomialTwoSpeedOU", TRUE )
