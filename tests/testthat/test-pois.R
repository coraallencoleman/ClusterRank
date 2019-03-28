context("test-pois")
test_that("imports data", {
  expect_equal(nrow(poisData), 8)
})

test_that("Wilson Hilferty CI works", {
  expect_equal(WilsonHilfertyPoiCI(33, 1, conf.level=0.95), c(33.00000, 22.71198, 46.34576), tolerance = .001)
})

test_that("unweighted ranking works", {
  set.seed(123)
  unweight <- ClusterRankPois(poisData$lbw, ti = poisData$births,row.names=poisData$county,weighted=FALSE)
  expect_equal(unweight$thetas, c(603.6664,744.9961, 1454.0000, 5306.1863, 5590.3703, 5802.2132))
  expect_equal(as.vector(unweight$ranked.table$name), c("Middlesex", "Tolland", "Windham", "Litchfield", "New London", "New Haven", "Fairfield", "Hartford"))
})

test_that("weighted ranking works", {
  set.seed(123)
  weight <- ClusterRankPois(poisData$lbw,  ti = poisData$births,row.names=poisData$county,weighted=TRUE)
  expect_equal(weight$thetas, c(603.6664,744.9961, 1454.0000, 5306.1863, 5590.3703, 5802.2132))
  expect_equal(as.vector(weight$ranked.table$name), c("Middlesex", "Tolland", "Windham", "Litchfield", "New London", "New Haven", "Fairfield", "Hartford"))

})

test_that("unweighted ranking on rank scale works", {
  set.seed(123)
  unweight_rank <- ClusterRankPois(poisData$lbw,  ti = poisData$births,row.names=poisData$county, scale=rank, weighted=FALSE)
  expect_equal(unweight_rank$thetas, c(603.6664,744.9961, 1454.0000, 5306.1863, 5590.3703, 5802.2132))
  expect_equal(as.vector(unweight_rank$ranked.table$name), c("Middlesex", "Tolland", "Windham", "Litchfield", "New London", "New Haven", "Fairfield", "Hartford"))
})

test_that("weighted ranking on rank scale works", {
  set.seed(123)
  weight_rank <- ClusterRankPois(poisData$lbw,  ti = poisData$births,row.names=poisData$county, scale=rank, weighted=TRUE)
  expect_equal(weight_rank$thetas, c(603.6664,744.9961, 1454.0000, 5306.1863, 5590.3703, 5802.2132))
  expect_equal(as.vector(weight_rank$ranked.table$name), c("Middlesex", "Tolland", "Windham", "Litchfield", "New London", "New Haven", "Fairfield", "Hartford"))
})

test_that("null row.names weighted ranking on rank scale works", {
  set.seed(123)
  weight_rank <- ClusterRankPois(poisData$lbw,  ti = poisData$births,scale=rank, weighted=TRUE)
  expect_equal(weight_rank$thetas, c(603.6664,744.9961, 1454.0000, 5306.1863, 5590.3703, 5802.2132))
  expect_equal(as.vector(weight_rank$ranked.table$name), paste(c(4, 7, 8, 3, 6, 5, 1, 2)))
})

test_that("Errors with message if the NPMLE is point mass at a single value", { #GitHub issue #8
  y<-rep(10,20)
  expect_error(ClusterRankPois(y,row.names=paste("County",1:20),sig.digits=2), regexp = "^All units have identical point mass at")
})

#TODO
test_that("handles two-way tie case", { #GitHub issue #9
  y <- c(8, 12, 13, 8)
})

test_that("handles multi-way tie case", { #GitHub issue #9
  y <- c( 8, 12, 13, 8, 17, 12, 11, 14, 6, 17, 15, 9, 11, 6, 11, 7, 14, 14, 11, 8)
  n <- rep(100,20)

  #"county 5" with the most (17) events is assigned rank 13,
  #because rank positions 13-20 all belong with probability 1 to the cluster 2
})

#test plotting later? plot_rt(rc = unweight)
