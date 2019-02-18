context("test-pois")
test_that("imports data", {
  expect_equal(nrow(poisData), 8)
})

test_that("Wilson Hilferty CI works", {
  expect_equal(WilsonHilfertyPoiCI(33, 1, conf.level=0.95), c(33.00000, 22.71198, 46.34576), tolerance = .001)
})

test_that("unweighted ranking works", {
  set.seed(123)
  unweight <- ClusterRankPois(poisData$lbw, row.names=poisData$county,weighted=FALSE)
  expect_equal(unweight$theta, c(603.6664,744.9961, 1454.0000, 5306.1863, 5590.3703, 5802.2132))
  expect_equal(as.vector(unweight$ranked_table$name), c("Middlesex", "Tolland", "Windham", "Litchfield", "New London", "New Haven", "Fairfield", "Hartford"))
})

test_that("weighted ranking works", {
  set.seed(123)
  weight <- ClusterRankPois(poisData$lbw, row.names=poisData$county,weighted=TRUE)
  expect_equal(weight$theta, c(603.6664,744.9961, 1454.0000, 5306.1863, 5590.3703, 5802.2132))
  expect_equal(as.vector(weight$ranked_table$name), c("Middlesex", "Tolland", "Windham", "Litchfield", "New London", "New Haven", "Fairfield", "Hartford"))

})

test_that("unweighted ranking on rank scale works", {
  set.seed(123)
  unweight_rank <- ClusterRankPois(poisData$lbw, row.names=poisData$county, scale=rank, weighted=FALSE)
  expect_equal(unweight_rank$theta, c(603.6664,744.9961, 1454.0000, 5306.1863, 5590.3703, 5802.2132))
  expect_equal(as.vector(unweight_rank$ranked_table$name), c("Middlesex", "Tolland", "Windham", "Litchfield", "New London", "New Haven", "Fairfield", "Hartford"))
})

test_that("weighted ranking on rank scale works", {
  set.seed(123)
  weight_rank <- ClusterRankPois(poisData$lbw, row.names=poisData$county, scale=rank, weighted=TRUE)
  expect_equal(weight_rank$theta, c(603.6664,744.9961, 1454.0000, 5306.1863, 5590.3703, 5802.2132))
  expect_equal(as.vector(weight_rank$ranked_table$name), c("Middlesex", "Tolland", "Windham", "Litchfield", "New London", "New Haven", "Fairfield", "Hartford"))
})

#test plotting later? plot_rt(rc = unweight)
