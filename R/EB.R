#' @importFrom PCMBase PCMCond
#' @importFrom PCMBase PCMDescribe
#' @importFrom PCMBase PCMDescribeParameters
#' @importFrom PCMBase PCMInfo
#' @importFrom PCMBase PCMListDefaultParameterizations
#' @importFrom PCMBase PCMListParameterizations
#' @importFrom PCMBase PCMParentClasses
#' @importFrom PCMBase PCMSpecify
#'

#' @export
GetSigma_x <- function(
    o, name = "Sigma", r = 1,
    transpose = getOption("PCMBase.Transpose.Sigma_x", FALSE)) {
  
  name <- paste0(name, "_x")
  S <- if(is.Global(o[[name]])) as.matrix(o[[name]]) else as.matrix(o[[name]][,, r])
  
  if(transpose) {
    t(S)
  } else {
    S
  }
}


#' @export
PCMInfo.EB <- function(
    X, tree, model,
    SE = matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree)),
    verbose = FALSE, preorder = NULL, ...) {
  
  if(is.Transformable(model)) {
    model <- PCMApplyTransformation(model)
  }
  
  res <- NextMethod()
  
  res$nodeHeights <- nodeHeights(tree)
  
  #Removed if (res$RTree!=1) part
  
  res
}

#' @export
PCMParentClasses.EB <- function(model) {
  c("GaussianPCM", "PCM")
}


#' @export
PCMDescribe.EB <- function(model, ...) {
  "Early-Brust model"
}


#' @export
PCMCond.EB <- function(
    tree, model, r = 1, metaI = PCMInfo(NULL, tree, model, verbose = verbose),
    verbose=FALSE) {
  
  Sigma_x <- GetSigma_x(model, "Sigma", r)
  Sigma <- Sigma_x %*% t(Sigma_x)
  
  if(!is.null(model$Sigmae_x)) {
    Sigmae_x <- GetSigma_x(model, "Sigmae", r)
    Sigmae <- Sigmae_x %*% t(Sigmae_x)
  } else {
    Sigmae <- NULL
  }
  
  
  if(!is.null(model$rho)) {
    rho <- if(is.Global(model$rho)) model$rho else model$rho[r]
  }else{
    rho <- 1
  }
  
  
    # using the function V instead of PCMCondVOU_EB since H = 0
    V <- function(t, edgeIndex, metaI, ...) {
    
    # time from root to the start of the current branch
    t_s <- metaI$nodeHeights[edgeIndex, 1]
    
    # variance for a branch with EB exponential decay factor
    res <- t * Sigma * exp(-rho * (t + t_s) / 2)
    
    # Add non-heritable variance for tips of the tree
    if(!is.null(Sigmae) && metaI$edge[edgeIndex, 2] <= metaI$N) {
      res <- res + Sigmae
    }
    
    return(res)
  }
  omega <- function(t, edgeIndex, metaI) {
    rep(0, nrow(Sigma))
  }
  Phi <- function(t, edgeIndex, metaI, e_Ht = NULL) {
    diag(nrow(Sigma))
  }
  list(omega = omega, Phi = Phi, V = V)
}


#' @export
PCMDescribeParameters.EB <- function(model, ...) {
  list(
    X0 = "trait values at the root",
    rho = "time-dependent parameter affecting rates of evolution",
    Sigma_x = "Upper triangular factor of the unit-time variance rate",
    Sigmae_x = "Upper triangular factor of the non-heritable variance or the variance of the measurement error")
}


#' @export
PCMListParameterizations.EB <- function(model, ...) {
  list(
    X0 = list(
      c("VectorParameter", "_Global"),
      c("VectorParameter", "_Fixed", "_Global"),
      c("VectorParameter", "_AllEqual", "_Global"),
      c("VectorParameter", "_Omitted")),
    rho = list(
      c("ScalarParameter"),
      c("ScalarParameter", "_Fixed"),
      c("ScalarParameter", "_AllEqual"),
      c("ScalarParameter", "_Omitted")),
    Sigma_x = list(
      c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal")),
    
    Sigmae_x = list(
      c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
      c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
      c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
      c("MatrixParameter", "_Omitted"))
  )
}



#' @export
PCMListDefaultParameterizations.EB <- function(model, ...) {
  list(
    X0 = list(
      c("VectorParameter", "_Global"),
      c("VectorParameter", "_Omitted")
    ),
    rho = list(
      c("ScalarParameter")),
    
    Sigma_x = list(
      c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal")
    ),
    
    Sigmae_x = list(
      c("MatrixParameter", "_Omitted"))
  )
}


#' @export
PCMSpecify.EB <- function(model, ...) {
  spec <- list(
    X0 = structure(0.0, class = c('VectorParameter', '_Global'),
                   description = 'trait values at the root'),
    rho = structure(0.0, class = c('ScalarParameter'),
                    description = 'time-dependent parameter affecting rates of evolution'),
    Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
                        description = 'Upper triangular factor of the unit-time variance rate'),
    Sigmae_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', '_WithNonNegativeDiagonal'),
                         description = 'Upper triangular factor of the non-heritable variance or the variance of the measurement error'))
  attributes(spec) <- attributes(model)
  if(is.null(names(spec))) names(spec) <- c('X0', 'rho', 'Sigma_x', 'Sigmae_x')
  if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
  spec
}


#'@export
PCMInfoCpp.EB <- function(
    X, tree, model,
    SE = matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree)),
    metaI = PCMInfo(X, tree, model, SE, verbose, preorder=PCMTreePreorderCpp(tree)),
    verbose = FALSE, ...) {
  
  if(is.Transformable(model)) {
    model <- PCMApplyTransformation(model)
  }
  
  metaI$pcListInt <- PCListInt(metaI$pc)
  
  if(metaI$k == 1L && getOption("PCMBase.Use1DClasses", FALSE)) {
    res <- c(metaI, cppObject = PCMBaseCpp__QuadraticPolyEB1D$new(X, tree, model, metaI))
  } else {
    res <- c(metaI, cppObject = PCMBaseCpp__QuadraticPolyEB$new(X, tree, model, metaI))
  }
  
  res$TraverseTree <- res$cppObject$TraverseTree
  res$StateAtNode <- res$cppObject$StateAtNode
  
  class(res) <- c("PCMInfoCpp", class(metaI))
  res
}