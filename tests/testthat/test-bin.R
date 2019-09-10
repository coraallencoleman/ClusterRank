context("test_bin")

test_that("imports data", {
  expect_equal(nrow(binData), 8)
  expect_equal(nrow(AZbinData), 15)
})

test_that("unweighted ranking works", {
  set.seed(123)
  unweight <- ClusterRankBin(binData$lbw,binData$births,row.names=binData$county,weighted=FALSE)
  expect_equal(unweight$cluster.thetas, c(0.06100767, 0.07160035, 0.07573636, 0.08206402, 0.08547410))
  expect_equal(as.vector(unweight$ranked.table$name), c("Middlesex", "Tolland", "Litchfield", "Windham", "New London", "Fairfield", "New Haven", "Hartford"))
})

test_that("weighted ranking works", {
  set.seed(123)
  weight <- ClusterRankBin(binData$lbw,binData$births,row.names=binData$county,weighted=TRUE)
  expect_equal(weight$cluster.thetas, c(0.06100767, 0.07160035, 0.07573636, 0.08206402, 0.08547410))
  expect_equal(as.vector(weight$ranked.table$name), c("Middlesex", "Tolland", "Litchfield", "Windham", "New London", "Fairfield", "New Haven", "Hartford"))
})

test_that("unweighted ranking on rank scale works", {
  set.seed(123)
  unweight_rank <- ClusterRankBin(binData$lbw,binData$births,row.names=binData$county, scale=rank, weighted=FALSE)
  expect_equal(unweight_rank$cluster.thetas, c(0.06100767, 0.07160035, 0.07573636, 0.08206402, 0.08547410))
})

test_that("weighted ranking on rank scale works", {
  set.seed(123)
  weight_rank <- ClusterRankBin(binData$lbw,binData$births,row.names=binData$county, scale=rank, weighted=TRUE)
    expect_equal(weight_rank$cluster.thetas, c(0.06100767, 0.07160035, 0.07573636, 0.08206402, 0.08547410))
})

test_that("null row.names weighted ranking on rank scale works", {
  set.seed(123)
  weight_rank_NULL <- ClusterRankBin(binData$lbw,binData$births, scale=rank, weighted=TRUE)
  expect_equal(weight_rank_NULL$cluster.thetas, c(0.06100767, 0.07160035, 0.07573636, 0.08206402, 0.08547410))
  expect_equal(as.vector(weight_rank_NULL$ranked.table$name), paste(c(4, 3, 7, 8, 6, 1, 5, 2)))
  })

test_that("Errors with message if the NPMLE is point mass at a single value", { #GitHub issue #8
  y<-rep(10,20)
  n<-rep(50,20)
  expect_error(ClusterRankBin(y,n,row.names=paste("County",1:20),sig.digits=2), regexp = "^All units have identical point mass at")
})

##AZ case
test_that("Arizona clusters and ranks match Ron's code", {
  set.seed(123)
  weight_rank_NULL <- ClusterRankBin(AZbinData$lbw, AZbinData$births, scale=rank, weighted=TRUE, row.names = AZbinData$county)
  expect_equal(round(weight_rank_NULL$cluster.thetas, 3), c(0.059, 0.068, 0.070, 0.071, 0.078,0.084))
  expect_equal(as.vector(weight_rank_NULL$ranked.table$name), paste(c("Yuma", "La Paz", "Greenlee", "Pinal", "Maricopa", "Mohave", "Yavapai", "Pima", "Apache", "Graham", "Santa Cruz", "Coconino", "Cochise", "Gila", "Navajo")))
})

#possibel future TODO: this isnt caught as tie
# test_that("handles two-way tie case", { #GitHub issue #9
#   y <- c(8, 12, 13, 8)
#   n <- rep(100,4)
#   expect_error(ClusterRankBin(y,n,row.names=paste("County",1:4),sig.digits=2), regexp = "^Ties exist in cluster assignment")
# })
#
# test_that("handles multi-way tie case", { #GitHub issue #9
#   y <- c( 8, 12, 13, 8, 17, 12, 11, 14, 6, 17, 15, 9, 11, 6, 11, 7, 14, 14, 11, 8)
#   n <- rep(100,20)
#   tieMany <- ClusterRankBin(y,n,row.names=paste("County",1:20),sig.digits=2)
#   #"county 5" with the most (17) events is assigned rank 13,
#   #because rank positions 13-20 all belong with probability 1 to the cluster 2
# })

#test plotting later? plot_rt(rc = unweight)
