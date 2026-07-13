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
#'   distance ~ age + Sex,
#'   data = nlme::Orthodont,
#'   correlation = nlme::corCompSymm(form = ~ 1 | Subject),
#'   method = "ML"
#' )
#'
#' covariance <- robust.cov(fit)
#' estimate <- stats::coef(fit)
#'
#' # Illustrative between-within degrees of freedom, as used in the homework.
#' group <- nlme::getGroups(fit)
#' n_subjects <- length(unique(group))
#' n_observations <- length(stats::residuals(fit))
#' between_effect <- names(estimate) %in% c("(Intercept)", "SexFemale")
#' degrees_freedom <- ifelse(
#'   between_effect,
#'   n_subjects - sum(between_effect),
#'   n_observations - n_subjects - sum(!between_effect)
#' )
#'
#' model_t <- estimate / covariance$se.model
#' robust_t <- estimate / covariance$se.robust
#'
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
