# Extract fixed-effect coefficients through the interface for each nlme class.
subject_coefficients <- function(fit) {
  if (inherits(fit, "lme")) nlme::fixef(fit) else stats::coef(fit)
}

# Resolve and validate one subject column, inferring it from the fit if omitted.
subject_id_column <- function(fit, data, id) {
  if (is.null(id)) {
    id_variables <- all.vars(nlme::getGroupsFormula(fit))
    if (length(id_variables) != 1L) {
      stop(
        "Supply `id` as one column name for this grouping structure.",
        call. = FALSE
      )
    }
    id <- id_variables
  }
  if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(id)) {
    stop("`id` must be one non-missing column name.", call. = FALSE)
  }
  if (!id %in% names(data)) {
    stop(sprintf("`%s` is not a column in `data`.", id), call. = FALSE)
  }
  if (anyNA(data[[id]])) {
    stop(
      sprintf("Subject column `%s` must not contain missing values.", id),
      call. = FALSE
    )
  }
  id
}

# Refit one model after removing a subject, optionally overriding fit control.
refit_without_subject <- function(fit, data, id, subject, control) {
  keep <- as.character(data[[id]]) != subject
  reduced_data <- droplevels(data[keep, , drop = FALSE])
  refit_call <- if (is.null(control)) {
    stats::update(fit, data = reduced_data, evaluate = FALSE)
  } else {
    stats::update(
      fit,
      data = reduced_data,
      control = control,
      evaluate = FALSE
    )
  }
  if (inherits(fit, "lme")) {
    refit_call[[1L]] <- quote(nlme::lme)
  }
  eval(refit_call, envir = environment())
}

#' Compute Leave-One-Subject-Out Cook's Distances
#'
#' Measures each subject's influence on the fixed-effect estimates by removing
#' that subject, refitting an [nlme::lme()] or [nlme::gls()] model, and computing
#' a coefficient-based Cook's distance.
#'
#' @param fit A fitted [nlme::lme()] or [nlme::gls()] model.
#' @param data Data used for the refits. By default, it is obtained with
#'   [nlme::getData()].
#' @param id An optional single character column name identifying subjects. If
#'   `NULL`, the function infers the single grouping variable from `fit`.
#' @param control An optional [nlme::lmeControl()] or [nlme::glsControl()] object
#'   passed to every refit.
#'
#' @return A data frame with `id`, `cooks_distance`, `converged`, and `message`.
#'   Failed refits have `NA` distance and retain their error message.
#' @export
#'
#' @examples
#' ids <- c("M01", "M02", "M03", "F01", "F02", "F03")
#' orthodont <- droplevels(subset(
#'   as.data.frame(nlme::Orthodont),
#'   Subject %in% ids
#' ))
#' lme_fit <- nlme::lme(
#'   distance ~ age + Sex,
#'   data = orthodont,
#'   random = ~ 1 | Subject,
#'   method = "ML"
#' )
#' subject_cooks_distance(lme_fit, data = orthodont)
#'
#' gls_fit <- nlme::gls(
#'   distance ~ age + Sex,
#'   data = orthodont,
#'   correlation = nlme::corCompSymm(form = ~ 1 | Subject),
#'   method = "ML"
#' )
#' subject_cooks_distance(gls_fit, data = orthodont, id = "Subject")
subject_cooks_distance <- function(
    fit,
    data = nlme::getData(fit),
    id = NULL,
    control = NULL
) {
  if (!inherits(fit, c("lme", "gls"))) {
    stop(
      "`subject_cooks_distance()` requires an nlme::lme or nlme::gls model.",
      call. = FALSE
    )
  }

  # Validate the refit data and recover the subjects used by the fitted model.
  id <- subject_id_column(fit, data, id)
  subjects <- unique(as.character(nlme::getGroups(fit)))
  if (!all(subjects %in% as.character(data[[id]]))) {
    stop("Not all fitted subjects occur in the supplied `data`.", call. = FALSE)
  }

  # Compute the full-fit ingredients shared by every subject deletion.
  full_coef <- subject_coefficients(fit)
  full_covariance <- stats::vcov(fit)
  fixed_formula <- stats::formula(fit)
  fixed_frame <- stats::model.frame(
    fixed_formula,
    data = data,
    na.action = stats::na.omit
  )
  fixed_rank <- qr(stats::model.matrix(fixed_formula, fixed_frame))$rank

  # Isolate errors to the affected subject so other distances remain available.
  results <- lapply(subjects, function(subject) {
    refit <- tryCatch(
      refit_without_subject(fit, data, id, subject, control),
      error = function(error) error
    )
    if (inherits(refit, "error")) {
      return(data.frame(
        id = subject,
        cooks_distance = NA_real_,
        converged = FALSE,
        message = conditionMessage(refit),
        stringsAsFactors = FALSE
      ))
    }

    reduced_coef <- subject_coefficients(refit)
    if (!identical(names(reduced_coef), names(full_coef))) {
      return(data.frame(
        id = subject,
        cooks_distance = NA_real_,
        converged = FALSE,
        message = "Refit changed the fixed-effect coefficient structure.",
        stringsAsFactors = FALSE
      ))
    }

    delta <- reduced_coef - full_coef
    distance <- as.numeric(
      crossprod(delta, solve(full_covariance, delta)) / fixed_rank
    )
    data.frame(
      id = subject,
      cooks_distance = distance,
      converged = TRUE,
      message = NA_character_,
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, results)
  rownames(result) <- NULL
  result
}
