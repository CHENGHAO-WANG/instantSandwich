#' Compute Model-Based and Sandwich Covariance Matrices
#'
#' Computes model-based and cluster-robust sandwich covariance matrices and
#' standard errors from a fitted longitudinal model supported by [nlme].
#'
#' @param u A fitted `nlme` model that supports [nlme::getData()],
#'   [nlme::getGroups()], and [nlme::getVarCov()].
#'
#' @return A named list containing `Sig.model`, the model-based coefficient
#'   covariance matrix; `se.model`, its standard errors; `Sig.robust`, the
#'   sandwich-robust coefficient covariance matrix; and `se.robust`, its
#'   standard errors.
#' @export
#'
#' @examples
#' fit <- nlme::gls(
#'   distance ~ 0 + Sex + age + Sex:age,
#'   data = nlme::Orthodont,
#'   correlation = nlme::corCompSymm(form = ~ 1 | Subject),
#'   method = "REML"
#' )
#'
#' covariance <- robust.cov(fit)
#' estimate <- stats::coef(fit)
#' group <- nlme::getGroups(fit)
#' n_subjects <- length(unique(group))
#' n_observations <- length(stats::residuals(fit))
#'
#' # Between-within degrees of freedom from the BIOS 767 course example.
#' degrees_freedom <- c(
#'   rep(n_subjects - 2, 2),
#'   rep(n_observations - n_subjects - 2, 2)
#' )
#'
#' model_t <- estimate / covariance$se.model
#' robust_t <- estimate / covariance$se.robust
#' coefficient_tests <- data.frame(
#'   estimate = estimate,
#'   model_se = covariance$se.model,
#'   model_t = model_t,
#'   model_p = 2 * stats::pt(-abs(model_t), df = degrees_freedom),
#'   robust_se = covariance$se.robust,
#'   robust_t = robust_t,
#'   robust_p = 2 * stats::pt(-abs(robust_t), df = degrees_freedom)
#' )
#' coefficient_tests
#'
#' # Test parallel age trends by testing the Sex-by-age interaction.
#' contrast <- matrix(c(0, 0, 0, 1), nrow = 1)
#' contrast_estimate <- contrast %*% estimate
#' contrast_covariance <- contrast %*% covariance$Sig.model %*% t(contrast)
#' f_statistic <- as.numeric(
#'   contrast_estimate %*% solve(contrast_covariance) %*%
#'     t(contrast_estimate)
#' )
#' df_numerator <- 1
#' df_denominator <- n_observations - n_subjects - 2
#' parallelism_test <- data.frame(
#'   F = f_statistic,
#'   p_value = stats::pf(
#'     f_statistic,
#'     df1 = df_numerator,
#'     df2 = df_denominator,
#'     lower.tail = FALSE
#'   ),
#'   df_numerator = df_numerator,
#'   df_denominator = df_denominator
#' )
#' parallelism_test
robust.cov <- function(u) {
  form <- stats::formula(u)
  model_data <- stats::model.frame(form, nlme::getData(u))
  design <- stats::model.matrix(form, model_data)
  residual <- stats::residuals(u, type = "response")

  group <- as.character(nlme::getGroups(u))
  ids <- unique(group)
  n_coef <- ncol(design)
  bread_inverse <- matrix(0, nrow = n_coef, ncol = n_coef)
  meat <- matrix(0, nrow = n_coef, ncol = n_coef)

  # Accumulate bread and meat contributions independently for each subject.
  for (id in ids) {
    index <- which(group == id)
    design_i <- design[index, , drop = FALSE]
    residual_i <- residual[index]
    covariance_i <- as.matrix(
      nlme::getVarCov(u, individual = id, type = "marginal")
    )
    precision_i <- solve(covariance_i)

    bread_inverse <- bread_inverse +
      t(design_i) %*% precision_i %*% design_i
    meat <- meat +
      t(design_i) %*% precision_i %*%
      (residual_i %o% residual_i) %*% precision_i %*% design_i
  }

  covariance_model <- solve(bread_inverse)
  covariance_robust <- covariance_model %*% meat %*% covariance_model

  list(
    Sig.model = covariance_model,
    se.model = sqrt(diag(covariance_model)),
    Sig.robust = covariance_robust,
    se.robust = sqrt(diag(covariance_robust))
  )
}
