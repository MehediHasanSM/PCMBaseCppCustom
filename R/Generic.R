#' A generic method for creating C++ backend objects given a model, data and 
#' a tree.
#' @useDynLib PCMBaseCpp
#' @importFrom PCMBase PCMInfo PCMTreeJumps
#' 
#' @export
PCMInfoCpp <- function(X, tree, model, metaI = PCMInfo(X, tree, model, verbose), verbose = FALSE, ...) {
  UseMethod("PCMInfoCpp", model)
}

PCMExtractAbCdEfLmr <- function(metaI) {
  tr <- metaI$cppObject$tree
  M <- tr$num_nodes
  # internal node ids from 1 to M
  nodeIds <- sapply(1:M, function(n) tr$FindIdOfNode(n) + 1)
  
  spec <- metaI$cppObject$spec
  specList <- list(
    A = spec$A, b = spec$b, C = spec$C, d = spec$d, E = spec$E, f = spec$f, 
    L = spec$L, m = spec$m, r = spec$r, V = spec$V, V_1 = spec$V_1
  )
  with(specList,
       list(A = abind(lapply(nodeIds, function(i) A[,, i, drop=TRUE]), along = 3),
            b = abind(lapply(nodeIds, function(i) b[, i, drop=TRUE]), along = 2),
            C = abind(lapply(nodeIds, function(i) C[,, i, drop=TRUE]), along = 3),
            d = abind(lapply(nodeIds, function(i) d[, i, drop=TRUE]), along = 2),
            E = abind(lapply(nodeIds, function(i) E[,, i, drop=TRUE]), along = 3),
            f = as.vector(abind(lapply(nodeIds, function(i) f[i, drop=TRUE]), along = 1)),
            V = abind(lapply(nodeIds, function(i) V[,, i, drop=TRUE]), along = 3),
            V_1 = abind(lapply(nodeIds, function(i) V_1[,, i, drop=TRUE]), along = 3),
            L = abind(lapply(nodeIds, function(i) L[,, i, drop=TRUE]), along = 3),
            m = abind(lapply(nodeIds, function(i) m[, i, drop=TRUE]), along = 2),
            r = as.vector(abind(lapply(nodeIds, function(i) r[i, drop=TRUE]), along = 1))
            ))
}

#' @importFrom PCMBase PCMLmr
#' @importFrom PCMBase PCMInfo
#'
#' @export
PCMLmr.PCMInfoCpp <- function(
  X, tree, model, 
  metaI = PCMInfo(X, tree, model), 
  root.only = FALSE, verbose = FALSE
) {
  
  par <- PCMGetVecParamsFull(model)
  
  PCMLmr_vec <- metaI$TraverseTree(par, mode=getOption("PCMBase.Lmr.mode", as.integer(11)))
  
  if(root.only) {
    # number of traits (variables)
    k <- metaI$k
    
    list(L = matrix(PCMLmr_vec[1:(k*k)], k, k),
         m = PCMLmr_vec[k*k+(1:k)],
         r = PCMLmr_vec[k*k+k+1]
    )  
  } else {
    PCMExtractAbCdEfLmr(metaI)
  }
}