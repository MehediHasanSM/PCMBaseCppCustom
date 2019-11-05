library(PCMBaseCpp)

if(PCMBaseCppIsADevRelease()) {
  library(testthat)
  
  context("Test singular matrix V")  
  load("testobjects.RData")
  
  
  library(PCMBase)
  
  modelBM <- PCM(model = "BM", k = 1)
  
  mTraits <- matrix(0,ncol=PCMTreeNumTips(tree.a),nrow=1)
  
  metaI <- PCMInfo(mTraits, tree.a, modelBM)
  metaICpp <- PCMInfoCpp(mTraits, tree.a, modelBM)
  
  test_that(
    "Singular V matrix raises a warning (PCMBase.Errors.As.Warnings=TRUE):", {
      options(PCMBase.Errors.As.Warnings = TRUE)
      expect_warning(loglik_PCMBase <- PCMLik(mTraits, tree.a, modelBM, metaI = metaI))
      expect_warning(loglik_PCMBaseCpp <- PCMLik(mTraits, tree.a, modelBM, metaI = metaICpp))
      expect_true(is.na(loglik_PCMBase) && is.na(loglik_PCMBaseCpp))
    })
  
  test_that(
    "Singular V matrix raises an error (PCMBase.Errors.As.Warnings=FALSE):", {
      options(PCMBase.Errors.As.Warnings = FALSE)
      expect_error(loglik_PCMBase <- PCMLik(mTraits, tree.a, modelBM, metaI = metaI))
      expect_error(loglik_PCMBaseCpp <- PCMLik(mTraits, tree.a, modelBM, metaI = metaICpp))
    })
}