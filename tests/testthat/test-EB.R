library(testthat)
context("PCMLik, EB")

library(PCMBase)
library(PCMBaseCppCustom) 

if(PCMBaseCppIsADevRelease()) {
  
  load("testobjects.RData")
  
  set.seed(1)
  
  test_that("Equal R and Cpp likelihood values on a random model, single regime (a)", {
    expect_silent(model.a.123.EB <- PCM("EB", k = 3, regimes = "a"))
    expect_silent(PCMParamLoadOrStore(model.a.123.EB,
                                      PCMParamRandomVecParams(model.a.123.EB), 
                                      offset = 0, k = 3, load = TRUE))
    
    expect_silent(metaICpp <- PCMInfoCpp(X = traits.a.123[, 1:length(tree.a$tip.label)],
                                         tree = tree.a, model = model.a.123.EB))
    
    expect_equal(PCMLik(traits.a.123, tree.a, model.a.123.EB),
                 PCMLik(traits.a.123, tree.a, model.a.123.EB, metaI = metaICpp))
  })
  
  test_that("Equal R and Cpp likelihood values on a random model, single regime (a), SE>0", {
    expect_silent(model.a.123.EB <- PCM("EB", k = 3, regimes = "a"))
    expect_silent(PCMParamLoadOrStore(model.a.123.EB,
                                      PCMParamRandomVecParams(model.a.123.EB), 
                                      offset = 0, k = 3, load = TRUE))
    
    expect_silent(metaICpp <- PCMInfoCpp(
      X = traits.a.123[, 1:length(tree.a$tip.label)],
      tree = tree.a, 
      model = model.a.123.EB,
      SE = 0.01 * abs(traits.a.123[, 1:length(tree.a$tip.label)])))
    
    expect_equal(PCMLik(traits.a.123, tree.a, model.a.123.EB, SE = 0.01 * abs(traits.a.123[, 1:length(tree.a$tip.label)])),
                 PCMLik(traits.a.123, tree.a, model.a.123.EB, metaI = metaICpp))
  })
  
  test_that("Equal R and Cpp likelihood values on a random model, multiple regimes (ab)", {
    expect_silent(model.ab.123.EB <- PCM("EB", k = 3, regimes = c("a", "b")))
    expect_silent(PCMParamLoadOrStore(model.ab.123.EB,
                                      PCMParamRandomVecParams(model.ab.123.EB), 
                                      offset = 0, k = 3, load = TRUE))
    
    expect_silent(metaICpp <- PCMInfoCpp(X = traits.ab.123[, 1:length(tree.ab$tip.label)],
                                         tree = tree.ab, model = model.ab.123.EB))
    
    expect_equal(PCMLik(traits.ab.123, tree.ab, model.ab.123.EB),
                 PCMLik(traits.ab.123, tree.ab, model.ab.123.EB, metaI = metaICpp))
  })
}