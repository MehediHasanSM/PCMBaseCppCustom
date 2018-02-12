#'@export
PCMCppPruningObject.OU <- function(X, tree, model, metaInfo = PCMValidate(tree, model), ...) {
  QuadraticPolynomialOU$new(
    X, tree, metaInfo,
    threshold_detV = getOption("PCMBase.Threshold.SV", 1e-6),
    thresholdLambda_ij = getOption("PCMBase.Threshold.Lambda_ij", 1e-8),
    internal_pc_full = getOption("PCMBase.Internal.PC.Full", TRUE))
}

#' Calculate the coefficients L, m, r of the general
#' form (eq. 2) for the root node of a tree
#'
#' @param model parameters of the OU process. This must be a
#' named list with the following elements:
#' H: a k x k x R array, where R is the number of regimes of the
#' OU process, k is the number of variables (traits), each H[,,r]
#' containing the matrix H for regime r;
#' Theta: a k x R matrix, row Theta[, r] containing the long-term
#' mean Theta for regime r;
#' Sigma: a k x k x R array, each Sigma[,,r] containing the
#' matrix Sigma for regime r;
#' Sigmae: a k x k x R array, each Sigmae[,,r] representing a diagonal matrix
#' with elements on the diagona corresponding to the environmental variances for
#' the k traits in regime r
#' @param pruneI an object of class QuadraticPolynomialOU, which has been created using
#'   PCMCppPruningObject.OU.
#'
#' @importFrom abind abind
#' @importFrom PCMBase PCMValidate PCMPresentCoordinates
#' 
#' @export
#' 
#' @return a named list containing the following elements:
#' L: a k x k matrix
#' m: a k vector
#' r: a number;
PCMLmr.Rcpp_QuadraticPolynomialOU <- function(
  X, tree, model, 
  metaI = PCMValidate(tree, model), 
  pruneI = PCMCppPruningObject(X, tree, model, metaI), 
  pc = PCMPresentCoordinates(X, tree), 
  root.only = FALSE, verbose = FALSE
  ) {
  
  # number of regimes
  R <- metaI$R
  
  # number of traits (variables)
  k <- metaI$k
  
  par <- c(model$H, model$Theta, model$Sigma, model$Sigmae)
  
  PCMLmr_vec <- pruneI$TraverseTree(par, mode=getOption("splittree.postorder.mode", as.integer(0)))
  
  if(root.only) {
    list(L = matrix(PCMLmr_vec[1:(k*k)], k, k),
         m = PCMLmr_vec[k*k+(1:k)],
         r = PCMLmr_vec[k*k+k+1]
    )  
  } else {
    PCMExtractAbCdEfLmr(pruneI)
  }
}
