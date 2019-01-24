context("test-pois")

test_that("imports data", {
  expect_equal(nrow(poisData), 8)
})

test_that("unweighted ranking works", {
  unweight <- rank_cluster.pois(poisData$lbw,poisData$births,row_names=poisData$county,weighted=FALSE)
  expect_equal(unweight$theta, c(0.06100767, 0.07160035, 0.07573636, 0.08206402, 0.08547410))
})

test_that("weighted ranking works", {
  weight <- rank_cluster.pois(poisData$lbw,poisData$births,row_names=poisData$county,weighted=TRUE)
  expect_equal(weight$theta, c(0.06100767, 0.07160035, 0.07573636, 0.08206402, 0.08547410))
})

test_that("unweighted ranking on rank scale works", {
  unweight_rank <- rank_cluster.pois(poisData$lbw,poisData$births,row_names=poisData$county, scale=rank, weighted=FALSE)
  expect_equal(unweight_rank$theta, c(0.06100767, 0.07160035, 0.07573636, 0.08206402, 0.08547410))
})

test_that("weighted ranking on rank scale works", {
  weight_rank <- rank_cluster.pois(poisData$lbw,poisData$births,row_names=poisData$county, scale=rank, weighted=TRUE)
    expect_equal(weight_rank$theta, c(0.06100767, 0.07160035, 0.07573636, 0.08206402, 0.08547410))
})

#test plotting later? plot_rt(rc = unweight)
