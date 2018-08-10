loadModule( "PCMBaseCppTree", TRUE )
loadModule( "PCMBaseCppOrderedTree", TRUE )


#' Get a vector with all model parameters unrolled
#' @param model a PCM model object
#' @return a numerical vector
#' @export
PCMParamGetFullVector <- function(model, ...) {
  UseMethod("PCMParamGetFullVector", model)
}

#' @export
PCMParamGetFullVector.PCM <- function(model, ...) {
  R <- PCMNumRegimes(model)
  res <- do.call(c, lapply(names(model), function(name) {
    if(is.PCM(model[[name]])) {
      PCMParamGetFullVector(model[[name]], ...)
    } else if( is.Global(model[[name]]) ) {
      rep(as.vector(model[[name]]), R)
    } else {
      as.vector(model[[name]])
    }
  }))
  unname(res)
}

#' @export
PCMParamGetFullVector.MixedGaussian <- function(model, ...) {
  R <- PCMNumRegimes(model)
  k <- PCMNumTraits(model)
  spec <- attr(model, "spec", exact = TRUE)
  

  # X0 should be _Omitted in every sub-model
  vecX0 <- as.vector(model[["X0"]])
  
  # Sigmae_x can be _Omitted or not _Omitted in sub-models. If it is not _Omitted 
  # the Sigmae_x in the sub-model overwrites the Sigmae_x in the parent model 
  vecSigmae_x <- double(0)
  if(is.Omitted(spec[["Sigmae_x"]])) {
    # in this case Sigmae_x should be appended as fixed 0 (this is required by
    # the C++ classes)
    vecSigmae_x <- rep(0.0, k*k)
  } else {
    # will be appended to the Submodel's vectors only if they have Sigmae_x 
    # omitted.
    vecSigmae_x <- as.vector(model[["Sigmae_x"]])
  }
  
  # replicating the global parameters for each model
  res <- do.call(c, lapply(model, function(o) {
    if(is.PCM(o)) {
      vec <- c(vecX0, PCMParamGetFullVector(o, ...))
      if(! ("Sigmae_x" %in% names(o)) ) {
        vec <- c(vec, vecSigmae_x)
      }
      vec
    } else {
      double(0L)
    }
  }))
  
  unname(res)
}

#' A generic method for creating C++ backend objects given a model, data and 
#' a tree.
#' @useDynLib PCMBaseCpp
#' @importFrom PCMBase PCMInfo PCMTreeJumps PCMApplyTransformation is.Transformable
#' 
#' @export
PCMInfoCpp <- function(X, tree, model, 
                       metaI = PCMInfo(X, tree, model, verbose, preorder=PCMTreePreorderCpp(tree)), 
                       verbose = FALSE, ...) {
  UseMethod("PCMInfoCpp", model)
}

#' Fast preorder of the edges in a tree
#' @param tree a phylo object
#' @export
PCMTreePreorderCpp <- function(tree) {
  trCpp <- PCMBaseCppOrderedTree$new(tree)
  nodesInPreorder <- rev(trCpp$OrderNodes(1:PCMTreeNumNodes(tree)) + 1)[-1]
  match(nodesInPreorder, tree$edge[, 2])
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

#' @importFrom PCMBase PCMLmr PCMInfo PCMApplyTransformation is.Transformable
#'
#' @export
PCMLmr.PCMInfoCpp <- function(
  X, tree, model, 
  metaI = PCMInfo(X, tree, model), 
  root.only = FALSE, verbose = FALSE
) {
  
  if(is.Transformable(model)) {
    model <- PCMApplyTransformation(model)
  }
  
  par <- PCMParamGetFullVector(model)
  
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

#' Converts the logical matrix pc into a list of vectors denoting the (0-based) TRUE-indices in 
#' each column
PCListInt <- function(pc) {
  lapply(1:ncol(pc), function(i) (which(pc[, i]) - 1))
}
