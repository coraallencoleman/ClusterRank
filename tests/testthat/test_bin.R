context("test_bin")
#getwd() #TODO
#load(file = "data/binData.rda")

test_that("imports data", {
  expect_equal(5+1, 6)
  expect_equal(nrow(binData), 8)
})

set.seed(123)

test_that("unweighted ranking works", {
  unweight <- ClusterRankBin(binData$lbw,binData$births,row_names=binData$county,weighted=FALSE)
  expect_equal(unweight$theta, c(0.06100767, 0.07160035, 0.07573636, 0.08206402, 0.08547410))
  expect_equal(as.vector(unweight$ranked_table$name), c("Middlesex", "Tolland", "Litchfield", "Windham", "New London", "Fairfield", "New Haven", "Hartford"))
})

test_that("weighted ranking works", {
  weight <- ClusterRankBin(binData$lbw,binData$births,row_names=binData$county,weighted=TRUE)
  expect_equal(weight$theta, c(0.06100767, 0.07160035, 0.07573636, 0.08206402, 0.08547410))
  expect_equal(as.vector(weight$ranked_table$name), c("Middlesex", "Tolland", "Litchfield", "Windham", "New London", "Fairfield", "New Haven", "Hartford"))
})

test_that("unweighted ranking on rank scale works", {
  unweight_rank <- ClusterRankBin(binData$lbw,binData$births,row_names=binData$county, scale=rank, weighted=FALSE)
  expect_equal(unweight_rank$theta, c(0.06100767, 0.07160035, 0.07573636, 0.08206402, 0.08547410))
})

test_that("weighted ranking on rank scale works", {
  weight_rank <- ClusterRankBin(binData$lbw,binData$births,row_names=binData$county, scale=rank, weighted=TRUE)
    expect_equal(weight_rank$theta, c(0.06100767, 0.07160035, 0.07573636, 0.08206402, 0.08547410))
})

#test plotting later? plot_rt(rc = unweight)
