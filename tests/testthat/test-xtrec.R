test_that("grunfeld_data() returns correct structure", {
  dat <- grunfeld_data()
  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), 100L)
  expect_equal(ncol(dat), 4L)
  expect_named(dat, c("firm", "year", "invest", "mvalue"))
  expect_equal(length(unique(dat$firm)), 5L)
  expect_equal(length(unique(dat$year)), 20L)
})

test_that("xtrec t-REC runs without error on balanced panel", {
  dat <- grunfeld_data()
  res <- xtrec(dat, var = "invest", panel_id = "firm",
               time_id = "year", trend = 0L, robust = FALSE)
  expect_s3_class(res, "xtrec")
  expect_true(is.finite(res$statistic))
  expect_true(res$pvalue >= 0 && res$pvalue <= 1)
  expect_equal(res$N, 5L)
  expect_equal(res$TT, 20L)
  expect_false(res$robust)
  expect_equal(res$test_name, "t-REC")
})

test_that("xtrec t-REC with trend = 1 runs correctly", {
  dat <- grunfeld_data()
  res <- xtrec(dat, var = "invest", panel_id = "firm",
               time_id = "year", trend = 1L, robust = FALSE)
  expect_s3_class(res, "xtrec")
  expect_true(is.finite(res$statistic))
  expect_equal(res$kappa, 0.25)
})

test_that("xtrec stores correct bias coefficients for trend = 0", {
  dat <- grunfeld_data()
  res <- xtrec(dat, var = "invest", panel_id = "firm",
               time_id = "year", trend = 0L, robust = FALSE)
  expect_equal(res$a_p, 0.5, tolerance = 1e-10)
  expect_equal(res$b_p, 1/3, tolerance = 1e-10)
  expect_equal(res$kappa, 0.5)
})

test_that("print.xtrec produces output without error", {
  dat <- grunfeld_data()
  res <- xtrec(dat, var = "invest", panel_id = "firm",
               time_id = "year", trend = 0L, robust = FALSE)
  expect_invisible(print(res))
})

test_that("summary.xtrec produces output without error", {
  dat <- grunfeld_data()
  res <- xtrec(dat, var = "invest", panel_id = "firm",
               time_id = "year", trend = 0L, robust = FALSE)
  expect_invisible(summary(res))
})

test_that("xtrec errors on unbalanced panel", {
  dat <- grunfeld_data()
  dat_unbal <- dat[-1, ]
  expect_error(
    xtrec(dat_unbal, var = "invest", panel_id = "firm",
          time_id = "year", trend = 0L, robust = FALSE)
  )
})

test_that("xtrec errors when variable not found", {
  dat <- grunfeld_data()
  expect_error(
    xtrec(dat, var = "nonexistent", panel_id = "firm",
          time_id = "year", trend = 0L, robust = FALSE)
  )
})

test_that("xtrec errors on negative trend", {
  dat <- grunfeld_data()
  expect_error(
    xtrec(dat, var = "invest", panel_id = "firm",
          time_id = "year", trend = -1L, robust = FALSE)
  )
})
