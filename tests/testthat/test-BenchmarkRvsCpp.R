library(testthat)
context("BenchmarkRvsCpp")

library(PCMBase)
library(PCMBaseCppCustom)

test_that(
  "BenchmarkRvsCpp is passing.", 
  expect_error(
    expect_error(
      BenchmarkRvsCpp(ks = 2, includeParallelMode = FALSE, verbose = TRUE))))