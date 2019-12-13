library(mvMORPH)
library(PCMBase)
library(PCMBaseCpp)
library(ape)
library(data.table)
library(xtable)

set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")

options(PCMBase.Threshold.EV = 1e-8, PCMBase.Threshold.SV = 1e-8, 
        PCMBase.Use1DClasses = TRUE)


dtTimes <- NULL

Ns <- c(25, 50, 100, 150, 200, 250)
ks <- c(1, 2, 4, 8)

matTimes_mvMORPH <- matTimes_PCMBaseCpp <- matrix(0, length(ks), length(Ns))

rownames(matTimes_PCMBaseCpp) <- rownames(matTimes_mvMORPH) <- as.character(ks)
colnames(matTimes_PCMBaseCpp) <- colnames(matTimes_mvMORPH) <- as.character(Ns)

for(N in Ns) {
  for(k in ks) {
    
    # Generating a random tree
    tree <- rtree(N)

    # Providing a tree whith the shift mapped on
    tot<-max(nodeHeights(tree))
    age=tot-3    # The shift occured 3 Ma ago
    tree<-make.era.map(tree,c(0,age))

    # Convert the tree with mapped regimes to a PCMTree object
    pcmTree <- PCMTree(map.to.singleton(tree))
    PCMTreeSetRegimesForEdges(pcmTree, names(pcmTree[["edge.length"]]))

    
    # Simulating the trait evolution.
    alpha <- t(benchmarkData$modelOU[[1]]$H[seq_len(k), seq_len(k),1])
    sigma <- t(benchmarkData$modelOU[[1]]$Sigma_x[seq_len(k), seq_len(k),1])
    theta<-benchmarkData$modelOU[[1]]$Theta[seq_len(k), 1]
    
    data<-mvSIM(tree, param=list(
      sigma=sigma, alpha=alpha, ntraits=k, theta=theta,
      names_traits=paste0("Trait.", seq_len(k))), model="ER", nsim=1)
    
    # Create a log-likelihood calculation function for an OUBM model using mvMORPH:
    llmvMORPH <- mvSHIFT(
      tree, data, model = "OUBM", optimization = "fixed")[["llik"]]
    
    vecParamsMVMORPH <- c(alpha[lower.tri(alpha, diag = TRUE)],
                          sigma[lower.tri(sigma, diag = TRUE)],
                          theta = theta)
    
    # Calculating the log-likelihood value of the parameters using mvMORPH:
    llmvMORPH(vecParamsMVMORPH, root.mle = FALSE)
    
    # Create a PCM model object using PCMBase. For simplicity, we use here a 2-regime
    # OU model. Alternatively, we could have used a mixed Gaussian model with an OU
    # and a BM regime.
    pcmOUBM <- PCM("OU", k = k, regimes = c("1", "2"))
    
    # Specify the parameter values for the model
    pcmOUBM[["H"]][,,1] <- alpha %*% t(alpha) 
    pcmOUBM[["Sigma_x"]][,,1] <- pcmOUBM[["Sigma_x"]][,,2] <- UpperTriFactor(sigma %*% t(sigma))
    pcmOUBM[["Theta"]][,1] <- theta
    pcmOUBM[["X0"]][] <- theta

    
    metaICpp <- PCMInfoCpp(t(data), pcmTree, pcmOUBM)
    
    nReps <- if(N < 100) {
      100
    } else {
      10
    }
    
    time_mvMORPH <- system.time(
      for(i in seq_len(nReps)) value_mvMORPH <- llmvMORPH(vecParamsMVMORPH, root.mle = FALSE)) / nReps
    
    time_PCMBaseCpp <- system.time(
      for(i in seq_len(nReps)) value_PCMBase <- PCMLik(t(data), pcmTree, pcmOUBM, metaI = metaICpp)) / nReps
    
    dtTimes <- rbindlist(
      list(dtTimes, 
           data.table(
             N = N, k = k, 
             mvMORPH=time_mvMORPH[3]*1000, PCMBaseCpp = time_PCMBaseCpp[3]*1000,
             ll_mvMORPH = value_mvMORPH,
             ll_PCMBaseCpp = value_PCMBase)))
    
    matTimes_mvMORPH[as.character(k), as.character(N)] <- time_mvMORPH[3]*1000
    matTimes_PCMBaseCpp[as.character(k), as.character(N)] <- time_PCMBaseCpp[3]*1000
    print(dtTimes)
  }
}

xtable(matTimes_mvMORPH)
xtable(matTimes_PCMBaseCpp)
