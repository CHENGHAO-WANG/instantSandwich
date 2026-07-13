#' Count Parameters in an MMRM Fit
#'
#' Counts the estimated fixed-effect coefficients and covariance parameters in
#' a fitted mixed model for repeated measures.
#'
#' @param fit A fitted `mmrm` model.
#'
#' @return A single numeric value giving the total number of estimated model
#'   parameters.
#' @export
#'
#' @examples
#' if (requireNamespace("mmrm", quietly = TRUE)) {
#'   data("fev_data", package = "mmrm")
#'   fit <- mmrm::mmrm(
#'     FEV1 ~ RACE + SEX + ARMCD * AVISIT + mmrm::us(AVISIT | USUBJID),
#'     data = fev_data
#'   )
#'   dof_mmrm(fit)
#' }
dof_mmrm <- function(fit) {
  if (!inherits(fit, c("mmrm", "mmrm_fit", "mmrm_tmb"))) {
    stop("`dof_mmrm()` requires an mmrm model.", call. = FALSE)
  }
  if (!requireNamespace("mmrm", quietly = TRUE)) {
    stop("`dof_mmrm()` requires the mmrm package.", call. = FALSE)
  }

  n_fixed <- length(stats::coef(fit, complete = FALSE))
  n_covariance <- length(mmrm::component(fit, "theta_est"))
  n_fixed + n_covariance
}

#' Count Parameters Using GLS-Style Rules
#'
#' Counts all estimated parameters in an [nlme::gls()] or `mmrm` fit. The GLS
#' count includes regression coefficients, model-structure parameters, and the
#' residual standard deviation when it is estimated.
#'
#' @param fit A fitted [nlme::gls()] or `mmrm` model.
#'
#' @return A single numeric value giving the total number of estimated model
#'   parameters.
#' @export
#'
#' @examples
#' fit <- nlme::gls(
#'   distance ~ age + Sex,
#'   data = nlme::Orthodont,
#'   correlation = nlme::corCompSymm(form = ~ 1 | Subject),
#'   method = "ML"
#' )
#' dof_gls_style(fit)
dof_gls_style <- function(fit) {
  if (inherits(fit, "gls")) {
    fixed_sigma <- isTRUE(attr(fit[["modelStruct"]], "fixedSigma"))
    n_fixed <- fit$dims$p
    n_structure <- length(stats::coef(fit[["modelStruct"]]))
    return(n_fixed + n_structure + as.integer(!fixed_sigma))
  }

  if (inherits(fit, c("mmrm", "mmrm_fit", "mmrm_tmb"))) {
    return(dof_mmrm(fit))
  }

  stop(
    "Unsupported model class for degrees-of-freedom computation.",
    call. = FALSE
  )
}
