library(testthat)

context("Test singular matrix V")

library(PCMBase)
library(PCMBaseCpp)
library(ape)

modelBM <- PCM(model = "BM", k = 1)
phyltree <- rtree(5)
mTraits <- matrix(0,ncol=5,nrow=1)

metaI <- PCMInfo(mTraits, phyltree, modelBM)
metaICpp <- PCMInfoCpp(mTraits, phyltree, modelBM)


test_that(
  "Singular V matrix raises a warning (PCMBase.Errors.As.Warnings=TRUE):", {
    options(PCMBase.Errors.As.Warnings = TRUE)
    expect_warning(loglik_PCMBase <- PCMLik(mTraits, phyltree, modelBM, metaI = metaI))
    expect_warning(loglik_PCMBaseCpp <- PCMLik(mTraits, phyltree, modelBM, metaI = metaICpp))
    expect_true(is.na(loglik_PCMBase) && is.na(loglik_PCMBaseCpp))
  })

test_that(
  "Singular V matrix raises an error (PCMBase.Errors.As.Warnings=FALSE):", {
    options(PCMBase.Errors.As.Warnings = FALSE)
    expect_error(loglik_PCMBase <- PCMLik(mTraits, phyltree, modelBM, metaI = metaI))
    expect_error(loglik_PCMBaseCpp <- PCMLik(mTraits, phyltree, modelBM, metaI = metaICpp))
  })



