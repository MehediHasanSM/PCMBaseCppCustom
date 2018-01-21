#'@export
newCppObject.JOU <- function(X, tree, model, metaInfo = validateModel(tree, model), ...) {
  QuadraticPolynomialJOU$new(
    X, tree, metaInfo,
    threshold_detV = getOption("PCMBase.Threshold.SV", 1e-6),
    thresholdLambda_ij = getOption("PCMBase.Threshold.Lambda_ij", 1e-8))
}

#' Calculate the coefficients L, m, r of the general
#' form (eq. 2) for the root node of a tree
#'
#' @param model parameters of the JOU process. This must be a
#' named list with the following elements:
#' Alpha: a k x k x R array, where R is the number of regimes of the
#' JOU process, k is the number of variables (traits), each Alpha[,,r]
#' containing the matrix Alpha for regime r;
#' Theta: a k x R matrix, row Theta[, r] containing the long-term
#' mean Theta for regime r;
#' Sigma: a k x k x R array, each Sigma[,,r] containing the
#' matrix Sigma for regime r;
#' Sigmae: a k x k x R array, each Sigmae[,,r] representing a diagonal matrix
#' with elements on the diagona corresponding to the environmental variances for
#' the k traits in regime r
#' @param pruneI an object of class QuadraticPolynomialJOU, which has been created using
#'   newCppObject.JOU.
#'
#' @importFrom abind abind
#' @importFrom PCMBase validateModel presentCoordinates
#'
#' @export
#' 
#' @return a named list containing the following elements:
#' L: a k x k matrix
#' m: a k vector
#' r: a number;
Lmr.Rcpp_QuadraticPolynomialJOU <- function(
  X, tree, model, 
  metaI = validateModel(tree, model), 
  pruneI = newCppObject(X, tree, model, metaI), 
  pc = presentCoordinates(X, tree), 
  root.only = FALSE, verbose = FALSE
  ) {
  # number of regimes
  R <- metaI$R
  
  # number of traits (variables)
  k <- metaI$k
  
  par <- c(model$Alpha, model$Theta, model$Sigma, model$Sigmae, model$mj, model$Sigmaj)
  
  Lmr_vec <- pruneI$TraverseTree(par, mode=getOption("PCMBase.Lmr.mode", 0))
  
  if(root.only) {
    list(L = matrix(Lmr_vec[1:(k*k)], k, k),
         m = Lmr_vec[k*k+(1:k)],
         r = Lmr_vec[k*k+k+1]
    )  
  } else {
    extractAbCdEfLmr(pruneI)
  }
}
