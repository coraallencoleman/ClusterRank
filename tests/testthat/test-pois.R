context("test-pois")

test_that("imports data", {
  expect_equal(nrow(poisData), 8)
})

test_that("unweighted ranking works", {
  unweight <- ClusterRankPois(poisData$lbw,poisData$births, row_names=poisData$county,weighted=FALSE)
  expect_equal(unweight$theta, c(603.6664,744.9961, 1454.0000, 5306.1863, 5590.3703, 5802.2132))
  expect_equal(as.vector(unweight$ranked_table$name), c("Middlesex", "Tolland", "Windham", "Litchfield", "New London", "New Haven", "Fairfield", "Hartford"))
})

test_that("weighted ranking works", {
  weight <- ClusterRankPois(poisData$lbw,poisData$births, row_names=poisData$county,weighted=TRUE)
  expect_equal(weight$theta, c(603.6664,744.9961, 1454.0000, 5306.1863, 5590.3703, 5802.2132))
  expect_equal(as.vector(weight$ranked_table$name), c("Middlesex", "Tolland", "Windham", "Litchfield", "New London", "New Haven", "Fairfield", "Hartford"))

})

test_that("unweighted ranking on rank scale works", {
  unweight_rank <- ClusterRankPois(poisData$lbw,poisData$births, row_names=poisData$county, scale=rank, weighted=FALSE)
  expect_equal(unweight_rank$theta, c(603.6664,744.9961, 1454.0000, 5306.1863, 5590.3703, 5802.2132))
  expect_equal(as.vector(unweight_rank$ranked_table$name), c("Middlesex", "Tolland", "Windham", "Litchfield", "New London", "New Haven", "Fairfield", "Hartford"))
})

test_that("weighted ranking on rank scale works", {
  weight_rank <- ClusterRankPois(poisData$lbw,poisData$births,  row_names=poisData$county, scale=rank, weighted=TRUE)
  expect_equal(weight_rank$theta, c(603.6664,744.9961, 1454.0000, 5306.1863, 5590.3703, 5802.2132))
  expect_equal(as.vector(weight_rank$ranked_table$name), c("Middlesex", "Tolland", "Windham", "Litchfield", "New London", "New Haven", "Fairfield", "Hartford"))
})

#test plotting later? plot_rt(rc = unweight)
