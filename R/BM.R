#'@export
newCppObject.BM <- function(X, tree, model, metaInfo = validateModel(tree, model), ...) {
  QuadraticPolynomialBM$new(
    X, tree, metaInfo, threshold_detV = getOption("PCMBase.Threshold.SV", 1e-6),
    internal_pc_full = getOption("PCMBase.Internal.PC.Full", TRUE))
}

#' Calculate the coefficients L, m, r of the general
#' form (eq. 2) for the root node of a tree
#'
#' @param model parameters of the OU process. This must be a
#' named list with the following elements:
#' Sigma: a k x k x R array, each Sigma[,,r] containing the
#' matrix Sigma for regime r;
#' Sigmae: a k x k x R array, each Sigmae[,,r] representing a diagonal matrix
#' with elements on the diagona corresponding to the environmental variances for
#' the k traits in regime r
#' @param pruneI an object of class QuadraticPolynomialOU, which has been created using
#'   newCppObject.BM.
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
Lmr.Rcpp_QuadraticPolynomialBM <- function(
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
  
  par <- c(model$Sigma, model$Sigmae)
  
  Lmr_vec <- pruneI$TraverseTree(par, mode=getOption("splittree.postorder.mode", as.integer(0)))
  
  if(root.only) {
    list(L = matrix(Lmr_vec[1:(k*k)], k, k),
         m = Lmr_vec[k*k+(1:k)],
         r = Lmr_vec[k*k+k+1]
    )  
  } else {
    extractAbCdEfLmr(pruneI)
  }
}
