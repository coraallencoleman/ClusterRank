context("test-test_norm")

test_that("imports data", {
  expect_equal(nrow(normData), 8)
})

test_that("unweighted ranking works", {
  unweight <- ClusterRankNorm(normData$mean,se = normData$se, row_names=normData$county,weighted=FALSE)
  expect_equal(unweight$theta, c(2.931671, 3.081624, 3.459036), tolerance = .001)
  expect_equal(as.vector(unweight$ranked_table$name), c("Middlesex", "Tolland", "Litchfield", "Fairfield", "New London",
                                                 "Hartford", "New Haven", "Windham"))
})

test_that("weighted ranking works", {
  weight <- ClusterRankNorm(normData$mean,se = normData$se, row_names=normData$county,weighted=TRUE)
  expect_equal(weight$theta, c(2.931671, 3.081624, 3.459036), tolerance = .001)
  expect_equal(as.vector(weight$ranked_table$name), c("Middlesex", "Tolland", "Litchfield", "New London", "Fairfield", "Hartford", "New Haven", "Windham"))
})

test_that("unweighted ranking on rank scale works", {
  unweight_rank <- ClusterRankNorm(normData$mean,se = normData$se, row_names=normData$county, scale=rank, weighted=FALSE)
  expect_equal(unweight_rank$theta, c(2.931671, 3.081624, 3.459036), tolerance = .001)
  expect_equal(as.vector(unweight_rank$pr_theta), c(0.1454716, 0.4793399, 0.3751885), tolerance = .001)
  expect_equal(as.vector(unweight_rank$ranked_table$name), c("Middlesex", "Tolland", "Litchfield", "New London",
                                                "Fairfield", "Hartford", "New Haven", "Windham"))
})

test_that("weighted ranking on rank scale works", {
  weight_rank <- ClusterRankNorm(normData$mean,se = normData$se, row_names=normData$county, scale=rank, weighted=TRUE)
    expect_equal(weight_rank$theta, c(2.931671, 3.081624, 3.459036), tolerance = .001)
    expect_equal(as.vector(weight_rank$ranked_table$name), c("Middlesex", "Tolland", "Litchfield", "New London",
                                                  "Fairfield", "Hartford", "New Haven", "Windham"))
})

#test plotting later? plot_rt(rc = unweight)

