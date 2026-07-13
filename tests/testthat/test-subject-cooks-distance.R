make_influence_data <- function(ids = c(
    "M01", "M02", "M03", "F01", "F02", "F03"
)) {
  droplevels(subset(as.data.frame(nlme::Orthodont), Subject %in% ids))
}

make_influence_lme <- function(data = make_influence_data()) {
  nlme::lme(
    distance ~ age + Sex,
    data = data,
    random = ~ 1 | Subject,
    method = "ML"
  )
}

make_influence_gls <- function(data = make_influence_data()) {
  nlme::gls(
    distance ~ age + Sex,
    data = data,
    correlation = nlme::corCompSymm(form = ~ 1 | Subject),
    method = "ML"
  )
}

manual_subject_cook <- function(fit, data, subject) {
  full_coef <- if (inherits(fit, "lme")) {
    nlme::fixef(fit)
  } else {
    stats::coef(fit)
  }
  reduced_data <- droplevels(data[as.character(data$Subject) != subject, ])
  reduced_coef <- if (inherits(fit, "lme")) {
    reduced_fit <- nlme::lme(
      distance ~ age + Sex,
      data = reduced_data,
      random = ~ 1 | Subject,
      method = "ML"
    )
    nlme::fixef(reduced_fit)
  } else {
    reduced_fit <- nlme::gls(
      distance ~ age + Sex,
      data = reduced_data,
      correlation = nlme::corCompSymm(form = ~ 1 | Subject),
      method = "ML"
    )
    stats::coef(reduced_fit)
  }
  delta <- reduced_coef - full_coef
  fixed_frame <- stats::model.frame(
    stats::formula(fit),
    data = data,
    na.action = stats::na.omit
  )
  rank <- qr(stats::model.matrix(stats::formula(fit), fixed_frame))$rank

  as.numeric(crossprod(delta, solve(stats::vcov(fit), delta)) / rank)
}

test_that("subject_cooks_distance handles lme fits and infers subjects", {
  data <- make_influence_data()
  fit <- make_influence_lme(data)
  result <- subject_cooks_distance(fit, data = data)

  expect_named(result, c("id", "cooks_distance", "converged", "message"))
  expect_setequal(result$id, unique(as.character(nlme::getGroups(fit))))
  expect_true(all(result$converged))
  first_id <- result$id[[1]]
  expect_equal(
    result$cooks_distance[result$id == first_id],
    manual_subject_cook(fit, data, first_id),
    tolerance = 1e-10
  )
})

test_that("subject_cooks_distance handles gls fits and explicit IDs", {
  data <- make_influence_data()
  fit <- make_influence_gls(data)
  result <- subject_cooks_distance(fit, data = data, id = "Subject")

  expect_true(all(result$converged))
  first_id <- result$id[[1]]
  expect_equal(
    result$cooks_distance[result$id == first_id],
    manual_subject_cook(fit, data, first_id),
    tolerance = 1e-10
  )
})

test_that("subject_cooks_distance validates model and ID inputs", {
  expect_error(
    subject_cooks_distance(stats::lm(mpg ~ wt, data = mtcars)),
    "requires an nlme::lme or nlme::gls model"
  )
  fit <- make_influence_gls()
  expect_error(
    subject_cooks_distance(fit, id = "missing_column"),
    "not a column"
  )
})

test_that("subject_cooks_distance retains failed subject refits", {
  data <- make_influence_data(c("M01", "M02", "M03", "F01"))
  fit <- make_influence_lme(data)
  result <- subject_cooks_distance(fit, data = data)
  failed <- result[result$id == "F01", ]

  expect_false(failed$converged)
  expect_true(is.na(failed$cooks_distance))
  expect_true(nzchar(failed$message))
  expect_true(any(result$converged))
})
