# Copyright 2018 Venelin Mitov
#
# This file is part of PCMBaseCpp.
#
# PCMBase is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PCMBase is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PCMBase.  If not, see <http://www.gnu.org/licenses/>.

# Benchmark data for the likelihood calculation for the MGPM model

library(ape)
library(PCMBase)
library(data.table)
library(PCMBaseCpp)

if(FALSE) {
  set.seed(5, kind = "Mersenne-Twister", normal.kind = "Inversion")
  
  options(PCMBase.MaxLengthListCladePartitions = 1)
  
  Ns <- c(10, 100, 1000, 10000)
  
  trees <- lapply(Ns, function(N) {
    cat("Generating a tree of N=", N, "\n")
    tree <- PCMTree(rtree(N, TRUE, NULL, br = rexp, 1/4))
    tree
  })
  
  listsParts <- lapply(trees, function(tree) {
    N <- PCMTreeNumTips(tree)
    PCMTreeListCladePartitions(
      tree, 
      nNodes = if(N < 100) 1 else if(N < 1000) 3 else 10, 
      minCladeSize = if(N < 100) 2 else if(N < 1000) 15 else if(N < 10000) 25 else 250,
      verbose = TRUE)
  })
  
  for(i in seq_along(trees)) {
    PCMTreeSetPartition(trees[[i]], listsParts[[i]][[1]])
  }
}

if(TRUE) {
  options(PCMBase.Value.NA = -1e20)
  options(PCMBase.Lmr.mode = 11)
  
  simulatedModels <- MGPMDefaultModelTypes()
  argsMixedGaussian_SimulatedModels <- Args_MixedGaussian_MGPMDefaultModelTypes()
  
  
  # these options will affect the H matrix, so that we don't have to specify it
  # explicitly in the functions below
  options(PCMBase.ParamValue.LowerLimit = -1,
          PCMBase.ParamValue.UpperLimit = 1,
          PCMBase.ParamValue.LowerLimit.NonNegativeDiagonal = .2)
  
  PCMParamLowerLimit.BM <- function(o, k, R, ...) {
    o <- NextMethod()
    k <- attr(o, "k", exact = TRUE)
    R <- length(attr(o, "regimes", exact = TRUE))
    
    if(is.Global(o$Sigma_x)) {
      o$Sigma_x[1, 1] <- o$Sigma_x[1, 1] <- .05
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[1, 2] <- .0
      }
    } else {
      for(r in seq_len(R)) {
        o$Sigma_x[1, 1, r] <- o$Sigma_x[1, 1, r] <- .05
        if(!is.Diagonal(o$Sigma_x)) {
          o$Sigma_x[1, 2, r] <- .0
        }
      }
    }
    o
  }
  
  PCMParamUpperLimit.BM <- function(o, k, R, ...) {
    o <- NextMethod()
    k <- attr(o, "k", exact = TRUE)
    R <- length(attr(o, "regimes", exact = TRUE))
    
    if(is.Global(o$Sigma_x)) {
      o$Sigma_x[1, 1] <- o$Sigma_x[2, 2] <- .5
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[1, 2] <- .2
      }
    } else {
      for(r in seq_len(R)) {
        o$Sigma_x[1, 1, r] <- o$Sigma_x[2, 2, r] <- .5
        if(!is.Diagonal(o$Sigma_x)) {
          o$Sigma_x[1, 2, r] <- .2
        }
      }
    }
    o
  }
  
  PCMParamLowerLimit.OU <- function(o, k, R, ...) {
    o <- NextMethod()
    k <- attr(o, "k", exact = TRUE)
    R <- length(attr(o, "regimes", exact = TRUE))
    
    if(is.Global(o$Theta)) {
      o$Theta[1] <- 3.0
      o$Theta[2] <- 2.0
    } else {
      for(r in seq_len(R)) {
        o$Theta[1, r] <- 3.0
        o$Theta[2, r] <- 2.0
      }
    }
    if(is.Global(o$Sigma_x)) {
      o$Sigma_x[1, 1] <- o$Sigma_x[2, 2] <- .05
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[1, 2] <- -.0
      }
    } else {
      for(r in seq_len(R)) {
        o$Sigma_x[1, 1, r] <- o$Sigma_x[2, 2, r] <- .05
        if(!is.Diagonal(o$Sigma_x)) {
          o$Sigma_x[1, 2, r] <- .0
        }
      }
    }
    o
  }
  
  PCMParamUpperLimit.OU <- function(o, k, R, ...) {
    o <- NextMethod()
    k <- attr(o, "k", exact = TRUE)
    R <- length(attr(o, "regimes", exact = TRUE))
    
    if(is.Global(o$Theta)) {
      o$Theta[1] <- 3.0 + 3.0 * getOption("MGPMSimulations.Theta.factor", 1.0)
      o$Theta[2] <- 2.0 + 2.0 * getOption("MGPMSimulations.Theta.factor", 1.0)
    } else {
      for(r in seq_len(R)) {
        o$Theta[1, r] <- 3.0 + 3.0 * getOption("MGPMSimulations.Theta.factor", 1.0)
        o$Theta[2, r] <- 2.0 + 2.0 * getOption("MGPMSimulations.Theta.factor", 1.0)
      }
    }
    if(is.Global(o$Sigma_x)) {
      o$Sigma_x[1, 1] <- o$Sigma_x[2, 2] <- .5
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[1, 2] <- .2
      }
    } else {
      for(r in seq_len(R)) {
        o$Sigma_x[1, 1, r] <- o$Sigma_x[2, 2, r] <- .5
        if(!is.Diagonal(o$Sigma_x)) {
          o$Sigma_x[1, 2, r] <- .2
        }
      }
    }
    o
  }
}


if(FALSE) {
  set.seed(27, kind = "Mersenne-Twister", normal.kind = "Inversion")
  
  options(PCMBase.Threshold.EV = 1e-9)
  options(PCMBase.Threshold.SV = 1e-9)
  
  models <- lapply(trees, function(tree) {
    R <- PCMNumRegimes(tree)
    model <- MixedGaussian(
      k = 8, 
      modelTypes = simulatedModels, 
      #mapping = sample(seq_along(simulatedModels), R, TRUE),
      mapping = sample(c(2, 4, 5, 6), R, TRUE), 
      Sigmae_x = argsMixedGaussian_SimulatedModels$Sigmae_x)
    vec <- PCMParamRandomVecParams(model, k = PCMNumTraits(model), R = PCMNumRegimes(model))
    PCMParamLoadOrStore(model, vecParams = vec, k = PCMNumTraits(model), offset = 0, load = TRUE)
    model
  })
  
  benchmarkData <- data.table(tree = trees, model = models)
  
  benchmarkData[, X:=lapply(.I, function(i) {
    X <- PCMSim(tree[[i]], model = model[[i]], X0 = rep(0, PCMNumTraits(model[[i]])))
    X[, seq_len(PCMTreeNumTips(tree[[i]]))]
  })]
  
  benchmarkData[, ll:=lapply(.I, function(i) {
    PCMLik(X[[i]], tree[[i]], model[[i]])
  })]
  
  benchmarkData[, list(tree, model, ll)]
  
  usethis::use_data(benchmarkData, overwrite = TRUE)

}

if(FALSE) {
    
  set.seed(12, kind = "Mersenne-Twister", normal.kind = "Inversion")
  
  options(PCMBase.Threshold.EV = 1e-9)
  options(PCMBase.Threshold.SV = 1e-9)
  benchmarkData2 <- benchmarkData
  
  benchmarkData2[, modelOU:=lapply(.I, function(i) {
    k <- PCMNumTraits(model[[i]])
    R <- PCMNumRegimes(tree[[i]])
    model <- PCM(
      model = PCMDefaultModelTypes()[[5]], 
      modelTypes = PCMDefaultModelTypes()[[5]], 
      k = k, 
      regimes = PCMRegimes(tree[[i]]))
    vec <- PCMParamRandomVecParams(model, k = PCMNumTraits(model), R = PCMNumRegimes(model))
    PCMParamLoadOrStore(model, vecParams = vec, k = PCMNumTraits(model), offset = 0, load = TRUE)
    model
  })]
  
  benchmarkData2[, llOU:=lapply(.I, function(i) {
    PCMLik(X[[i]], tree[[i]], modelOU[[i]])
  })]
  benchmarkData2[, list(tree, ll, llOU)]
  
  benchmarkData2[, modelBM:=lapply(.I, function(i) {
    k <- PCMNumTraits(model[[i]])
    R <- PCMNumRegimes(tree[[i]])
    model <- PCM(
      model = PCMDefaultModelTypes()[[2]], 
      modelTypes = PCMDefaultModelTypes()[[2]], 
      k = k, 
      regimes = PCMRegimes(tree[[i]]))
    vec <- PCMParamRandomVecParams(model, k = PCMNumTraits(model), R = PCMNumRegimes(model))
    PCMParamLoadOrStore(model, vecParams = vec, k = PCMNumTraits(model), offset = 0, load = TRUE)
    model
  })]
  
  benchmarkData2[, llBM:=lapply(.I, function(i) {
    PCMLik(X[[i]], tree[[i]], modelBM[[i]])
  })]
  
  benchmarkData2[, list(tree, model, ll, llBM, llOU)]
  
  benchmarkData <- benchmarkData2
  usethis::use_data(benchmarkData, overwrite = TRUE)
}
