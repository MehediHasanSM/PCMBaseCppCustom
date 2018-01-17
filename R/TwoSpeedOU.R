#'@export
newCppObject.TwoSpeedOU <- function(X, tree, model, metaInfo = validateModel(tree, model), ...) {
  QuadraticPolynomialTwoSpeedOU$new(X, tree, metaInfo)
}

#' Calculate the coefficients L, m, r of the general
#' form (eq. 2) for the root node of a tree
#'
#' @param model parameters of the 2SpOU process. This must be a
#' named list with the following elements:
#' 
#' Alpha1: a R x k x k array, where R is the number of regimes of the
#' OU process, k is the number of variables (traits), each Alpha1[r,,]
#' containing the matrix Alpha1 for regime r;
#' 
#' Alpha2: a R x k x k array, where R is the number of regimes of the
#' OU process, k is the number of variables (traits), each Alpha2[r,,]
#' containing the matrix Alpha2 for regime r;
#' 
#' Theta: a R x k matrix, row Theta[r, ] containing the long-term
#' mean Theta for regime r;
#' 
#' Sigma: a R x k x k array, each Sigma[r,,] containing the
#' matrix Sigma for regime r;
#' 
#' Sigmae: a R x k x k array, each Sigmae[r,,] representing a diagonal matrix
#' with elements on the diagona corresponding to the environmental variances for
#' the k traits in regime r
#' @param pruneI an object of class QuadraticPolynomialOU, which has been created using
#'   newCppObject.OU.
#'
#' @export
#' 
#' @return a named list containing the following elements:
#' L: a k x k matrix
#' m: a k vector
#' r: a number;
Lmr.TwoSpeedOU <- function(model, metaI, pruneI) {
  # number of regimes
  R <- metaI$R
  
  # number of traits (variables)
  k <- metaI$k
  
  par <- c()
  for(r in 1:R) {
    par <- c(par, model$Alpha1[r,, , drop=FALSE])
  }
  for(r in 1:R) {
    par <- c(par, model$Alpha2[r,, , drop=FALSE])
  }
  for(r in 1:R) {
    par <- c(par, model$Theta[r, , drop=FALSE])
  }
  for(r in 1:R) {
    par <- c(par, model$Sigma[r,, , drop=FALSE])
  }
  for(r in 1:R) {
    par <- c(par, model$Sigmae[r,, , drop=FALSE])
  }
  
  Lmr_vec <- pruneI$TraverseTree(par, mode=getOption("PCMBase.Lmr.mode", 0))
  
  L <- matrix(Lmr_vec[1:(k*k)], k, k)
  m <- Lmr_vec[k*k+(1:k)]
  r <- Lmr_vec[k*k+k+1]
  
  list(L=L, m=m, r=r)
}
