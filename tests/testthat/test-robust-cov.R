make_test_fit <- function() {
  nlme::gls(
    distance ~ age + Sex,
    data = nlme::Orthodont,
    correlation = nlme::corCompSymm(form = ~ 1 | Subject),
    method = "ML"
  )
}

manual_covariance <- function(fit) {
  form <- formula(fit)
  fit_data <- nlme::getData(fit)
  design <- model.matrix(form, model.frame(form, fit_data))
  residual <- residuals(fit, type = "response")
  group <- as.character(nlme::getGroups(fit))

  contributions <- lapply(unique(group), function(id) {
    index <- group == id
    design_i <- design[index, , drop = FALSE]
    residual_i <- residual[index]
    covariance_i <- as.matrix(
      nlme::getVarCov(fit, individual = id, type = "marginal")
    )
    precision_i <- solve(covariance_i)

    list(
      bread = t(design_i) %*% precision_i %*% design_i,
      meat = t(design_i) %*% precision_i %*%
        tcrossprod(residual_i) %*% precision_i %*% design_i
    )
  })

  bread_inverse <- Reduce(`+`, lapply(contributions, `[[`, "bread"))
  meat <- Reduce(`+`, lapply(contributions, `[[`, "meat"))
  covariance_model <- solve(bread_inverse)
  covariance_robust <- covariance_model %*% meat %*% covariance_model

  list(
    Sig.model = covariance_model,
    se.model = sqrt(diag(covariance_model)),
    Sig.robust = covariance_robust,
    se.robust = sqrt(diag(covariance_robust))
  )
}

test_that("robust.cov preserves the homework return interface", {
  fit <- make_test_fit()

  result <- robust.cov(fit)

  expect_named(
    result,
    c("Sig.model", "se.model", "Sig.robust", "se.robust"),
    ignore.order = FALSE
  )
  expect_equal(dim(result$Sig.model), c(length(coef(fit)), length(coef(fit))))
  expect_equal(dim(result$Sig.robust), c(length(coef(fit)), length(coef(fit))))
  expect_length(result$se.model, length(coef(fit)))
  expect_length(result$se.robust, length(coef(fit)))
})

test_that("robust.cov matches an analytical cluster calculation", {
  fit <- make_test_fit()

  result <- robust.cov(fit)
  expected <- manual_covariance(fit)

  expect_equal(result$Sig.model, expected$Sig.model, tolerance = 1e-10)
  expect_equal(result$se.model, expected$se.model, tolerance = 1e-10)
  expect_equal(result$Sig.robust, expected$Sig.robust, tolerance = 1e-10)
  expect_equal(result$se.robust, expected$se.robust, tolerance = 1e-10)
})

test_that("the placeholder hello function is not exported", {
  expect_false("hello" %in% getNamespaceExports("instantSandwich"))
})
