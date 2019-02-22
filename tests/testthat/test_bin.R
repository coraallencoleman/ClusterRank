context("test_bin")

test_that("imports data", {
  expect_equal(5+1, 6)
  expect_equal(nrow(binData), 8)
})

test_that("unweighted ranking works", {
  set.seed(123)
  unweight <- ClusterRankBin(binData$lbw,binData$births,row.names=binData$county,weighted=FALSE)
  expect_equal(unweight$theta, c(0.06100767, 0.07160035, 0.07573636, 0.08206402, 0.08547410))
  expect_equal(as.vector(unweight$ranked_table$name), c("Middlesex", "Tolland", "Litchfield", "Windham", "New London", "Fairfield", "New Haven", "Hartford"))
})

test_that("weighted ranking works", {
  set.seed(123)
  weight <- ClusterRankBin(binData$lbw,binData$births,row.names=binData$county,weighted=TRUE)
  expect_equal(weight$theta, c(0.06100767, 0.07160035, 0.07573636, 0.08206402, 0.08547410))
  expect_equal(as.vector(weight$ranked_table$name), c("Middlesex", "Tolland", "Litchfield", "Windham", "New London", "Fairfield", "New Haven", "Hartford"))
})

test_that("unweighted ranking on rank scale works", {
  set.seed(123)
  unweight_rank <- ClusterRankBin(binData$lbw,binData$births,row.names=binData$county, scale=rank, weighted=FALSE)
  expect_equal(unweight_rank$theta, c(0.06100767, 0.07160035, 0.07573636, 0.08206402, 0.08547410))
})

test_that("weighted ranking on rank scale works", {
  set.seed(123)
  weight_rank <- ClusterRankBin(binData$lbw,binData$births,row.names=binData$county, scale=rank, weighted=TRUE)
    expect_equal(weight_rank$theta, c(0.06100767, 0.07160035, 0.07573636, 0.08206402, 0.08547410))
})

test_that("null row.names weighted ranking on rank scale works", {
  set.seed(123)
  weight_rank_NULL <- ClusterRankBin(binData$lbw,binData$births, scale=rank, weighted=TRUE)
  expect_equal(weight_rank_NULL$theta, c(0.06100767, 0.07160035, 0.07573636, 0.08206402, 0.08547410))
  expect_equal(as.vector(weight_rank_NULL$ranked_table$name), paste(c(4, 3, 7, 8, 6, 1, 5, 2)))
  })

#test plotting later? plot_rt(rc = unweight)
