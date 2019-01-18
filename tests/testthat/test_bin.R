context("test_bin")

library(ClusteredRanking)
setwd("/Users/cora/git_repos/ClusteredRanking/")
load(file = "data/binData.rda")

test_that("imports data", {
  expect_equal(5+1, 6)
  expect_equal(nrow(binData), 8)
})

test_that("ranking works", {
  lbw_rc_rnk <- rank_cluster.bin(binData$lbw,binData$births,row_names=binData$county,scale=rank)
  expect_equal(lbw_rc_rnk$theta, c(0.06100767, 0.07160035, 0.07573636, 0.08206402, 0.08547410))
})