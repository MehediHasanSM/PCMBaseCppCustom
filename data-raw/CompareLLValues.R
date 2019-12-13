library(mvMORPH)
library(PCMBase)
library(PCMBaseCpp)

set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")

# Generating a random tree
tree <- ape::rtree(50)

# Providing a tree whith the shift mapped on
tot <- max(nodeHeights(tree))
age <- tot-3    # The shift occured 3 Ma ago
tree <- make.era.map(tree, c(0, age))

# Convert the tree with mapped regimes to a PCMTree object
pcmTree <- PCMTree(map.to.singleton(tree))
PCMTreeSetRegimesForEdges(pcmTree, names(pcmTree[["edge.length"]]))

# Plot of the phylogenies for illustration that they are the same and have
# the same regime assignment (uncomment the two lines below to see the plots).
#plotSimmap(tree,fsize=0.6,node.numbers=FALSE,lwd=3, pts=FALSE)
#PCMTreePlot(pcmTree)

# Parameters for simulate the trait evolution.
alpha<-matrix(c(1,0.1,0,2),2) 
sigma<-matrix(c(.1,.1,0,.1),2)
theta<-c(2,3)

data<-mvSIM(tree, param=list(
  sigma=sigma, alpha=alpha, ntraits=2, theta=theta,
  names_traits=c("head.size","mouth.size")), model="OUBM", nsim=1)

# Create a log-likelihood calculation function for an OUBM model using mvMORPH:
llmvMORPH <- mvSHIFT(
  tree, data, model = "OUBM", optimization = "fixed")[["llik"]]

# Calculating the log-likelihood value of the parameters using mvMORPH:
llmvMORPH(vecParams <- c(alpha[lower.tri(alpha, diag = TRUE)],
                         sigma[lower.tri(sigma, diag = TRUE)],
                         theta = theta), root.mle = FALSE)

# Create a PCM model object using PCMBase. For simplicity, we use here a 2-regime
# OU model. Alternatively, we could have used a mixed Gaussian model with an OU
# and a BM regime.
pcmOUBM <- PCM("OU", k = 2, regimes = c("1", "2"))

# Specify the parameter values for the model
pcmOUBM[["H"]][,,1] <- alpha %*% t(alpha) 
pcmOUBM[["Sigma_x"]][,,1] <- pcmOUBM[["Sigma_x"]][,,2] <- UpperTriFactor(sigma %*% t(sigma))
pcmOUBM[["Theta"]][,1] <- theta
pcmOUBM[["X0"]][] <- theta

# Calculate the log-likelihood value using PCMBase:
PCMLik(t(data), pcmTree, pcmOUBM)

# Calculate the log-likelihood value using PCMBaseCpp (faster execution):
metaICpp <- PCMInfoCpp(t(data), pcmTree, pcmOUBM)
PCMLik(t(data), pcmTree, pcmOUBM, metaI = metaICpp)

# The values are matchine up to numerical error
all.equal(target = llmvMORPH(vecParams, root.mle = FALSE), 
          current = PCMLik(t(data), pcmTree, pcmOUBM), check.attributes = FALSE)
all.equal(target = llmvMORPH(vecParams, root.mle = FALSE),
          current = PCMLik(t(data), pcmTree, pcmOUBM, metaI = metaICpp),
          check.attributes = FALSE)

