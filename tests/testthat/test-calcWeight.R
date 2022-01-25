test_that("calcWeight", {
  expect_equal(calcWeight(matrix(1:6, 2, 3)), c(0.5, 0.5))
})
