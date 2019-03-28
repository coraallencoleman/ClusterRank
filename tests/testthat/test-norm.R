context("test-test_norm")
set.seed(123)
test_that("imports data", {
  expect_equal(nrow(normData), 8)
})

test_that("unweighted ranking works", {
  set.seed(123)
  unweight <- ClusterRankNorm(normData$mean,se = normData$se, row.names=normData$county,weighted=FALSE)
  expect_equal(unweight$cluster.thetas, c(2.931671, 3.081624, 3.459036), tolerance = .001)
  expect_equal(as.vector(unweight$ranked.table$name), c("Middlesex", "Tolland", "Litchfield", "New London", "Fairfield",
                                                 "Hartford", "New Haven", "Windham"))
})

test_that("weighted ranking works", {
  set.seed(123)
  weight <- ClusterRankNorm(normData$mean,se = normData$se, row.names=normData$county,weighted=TRUE)
  expect_equal(weight$cluster.thetas, c(2.931671, 3.081624, 3.459036), tolerance = .001)
  expect_equal(as.vector(weight$ranked.table$name), c("Middlesex", "Tolland", "Litchfield", "New London", "Fairfield", "Hartford", "New Haven", "Windham"))
})

test_that("unweighted ranking on rank scale works", {
  set.seed(123)
  unweight_rank <- ClusterRankNorm(normData$mean,se = normData$se, row.names=normData$county, scale=rank, weighted=FALSE)
  expect_equal(unweight_rank$cluster.thetas, c(2.931671, 3.081624, 3.459036), tolerance = .001)
  expect_equal(as.vector(unweight_rank$pr_thetas), c(0.1454716, 0.4793399, 0.3751885), tolerance = .001)
  expect_equal(as.vector(unweight_rank$ranked.table$name), c("Middlesex", "Tolland", "Litchfield", "New London", "Fairfield",
                                                "Hartford", "New Haven", "Windham"))
})

test_that("weighted ranking on rank scale works", {
  set.seed(123)
  weight_rank <- ClusterRankNorm(normData$mean,se = normData$se, row.names=normData$county, scale=rank, weighted=TRUE)
    expect_equal(weight_rank$ranked.table, c(2.931671, 3.081624, 3.459036), tolerance = .001)
    expect_equal(as.vector(weight_rank$ranked.table$name), c("Middlesex", "Tolland", "Litchfield", "New London",
                                                  "Fairfield", "Hartford", "New Haven", "Windham"))
})

test_that("null row.names weighted ranking on rank scale works", {
  set.seed(123)
  weight_rank <- ClusterRankNorm(normData$mean,se = normData$se, scale=rank, weighted=TRUE)
  expect_equal(weight_rank$ranked.table, c(2.931671, 3.081624, 3.459036), tolerance = .001)
  expect_equal(as.vector(weight_rank$ranked.table$name), paste(c(4, 7, 3, 6, 1, 2, 5, 8)))
})

#TODO
test_that("Errors with message if the NPMLE is point mass at a single value", { #GitHub issue #8
  y<-rep(10,20)
  se<-rep(50,20)
  expect_error(ClusterRankBin(y,se,row.names=paste("County",1:20),sig.digits=2), regexp = "^All units have identical point mass at")
  })

test_that("handles two-way tie case", { #GitHub issue #9
  y <- c(8, 12, 13, 8)
  se <- rep(100,4)
})

test_that("handles multi-way tie case", { #GitHub issue #9
  y <- c( 8, 12, 13, 8, 17, 12, 11, 14, 6, 17, 15, 9, 11, 6, 11, 7, 14, 14, 11, 8)
  se <- rep(100,20)

  #"county 5" with the most (17) events is assigned rank 13,
  #because rank positions 13-20 all belong with probability 1 to the cluster 2
})
#test plotting later? plot_rt(rc = unweight)

