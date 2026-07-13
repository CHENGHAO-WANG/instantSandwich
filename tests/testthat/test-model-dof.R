make_dof_gls_fit <- function(fixed_sigma = NULL) {
  ids <- c("M01", "M02", "M03", "F01", "F02", "F03")
  data <- droplevels(
    subset(as.data.frame(nlme::Orthodont), Subject %in% ids)
  )
  control <- if (is.null(fixed_sigma)) {
    nlme::glsControl()
  } else {
    nlme::glsControl(sigma = fixed_sigma)
  }

  nlme::gls(
    distance ~ age + Sex,
    data = data,
    correlation = nlme::corCompSymm(form = ~ 1 | Subject),
    method = "ML",
    control = control
  )
}

test_that("dof_gls_style counts all estimated gls parameters", {
  expect_equal(dof_gls_style(make_dof_gls_fit()), 5)
  expect_equal(dof_gls_style(make_dof_gls_fit(fixed_sigma = 1)), 4)
})

test_that("dof_gls_style rejects unsupported models", {
  expect_error(
    dof_gls_style(stats::lm(mpg ~ wt, data = mtcars)),
    "Unsupported model class"
  )
})

test_that("mmrm parameter counts include fixed and covariance parameters", {
  testthat::skip_if_not_installed("mmrm")
  data("fev_data", package = "mmrm")
  fit <- mmrm::mmrm(
    FEV1 ~ RACE + SEX + ARMCD * AVISIT + mmrm::us(AVISIT | USUBJID),
    data = fev_data
  )
  expected <- length(stats::coef(fit, complete = FALSE)) +
    length(mmrm::component(fit, "theta_est"))

  expect_equal(dof_mmrm(fit), expected)
  expect_equal(dof_gls_style(fit), expected)
})

test_that("dof_mmrm rejects unsupported models", {
  expect_error(
    dof_mmrm(stats::lm(mpg ~ wt, data = mtcars)),
    "requires an mmrm model"
  )
})
