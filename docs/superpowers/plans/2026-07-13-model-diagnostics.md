# Model Diagnostics Utilities Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add exported model-parameter counting helpers and a unified leave-one-subject-out Cook's-distance function for `nlme::lme` and `nlme::gls` fits.

**Architecture:** Keep parameter counting in a focused `R/model-dof.R` module and subject influence in `R/subject-cooks-distance.R`. The influence function infers or accepts one subject column, refits sequentially, isolates per-subject failures, and uses the analytical quadratic-form Cook's-distance calculation. Public behavior is specified first with `testthat`; roxygen2 generates exports and manuals.

**Tech Stack:** R, `nlme`, optional model support from `mmrm`, roxygen2, testthat edition 3, R CMD build/check.

---

## File map

- Create `R/model-dof.R`: exported `dof_mmrm()` and `dof_gls_style()`.
- Create `R/subject-cooks-distance.R`: exported Cook's-distance API and small internal helpers.
- Create `tests/testthat/test-model-dof.R`: parameter-count behavior and validation.
- Create `tests/testthat/test-subject-cooks-distance.R`: `lme`/`gls` behavior, analytical comparison, and validation.
- Modify `DESCRIPTION`: declare optional `mmrm` support.
- Regenerate `NAMESPACE`, `man/dof_mmrm.Rd`, `man/dof_gls_style.Rd`, and `man/subject_cooks_distance.Rd` from roxygen comments.

### Task 1: Model parameter-count helpers

**Files:**
- Create: `tests/testthat/test-model-dof.R`
- Create: `R/model-dof.R`
- Modify: `DESCRIPTION`
- Modify (generated): `NAMESPACE`
- Create (generated): `man/dof_mmrm.Rd`
- Create (generated): `man/dof_gls_style.Rd`

- [ ] **Step 1: Write failing tests for `gls` parameter counts and validation**

Create `tests/testthat/test-model-dof.R` with a compound-symmetry `gls` fit whose expected count is three regression coefficients, one correlation parameter, and one residual standard deviation. Also test fixed sigma and an unsupported model:

```r
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
  expect_error(dof_gls_style(stats::lm(mpg ~ wt, data = mtcars)),
               "Unsupported model class")
})
```

- [ ] **Step 2: Run the focused tests and verify the missing-function failure**

Run:

```powershell
$env:LC_ALL='C'; $env:LANG='C'; Rscript -e "testthat::test_file('tests/testthat/test-model-dof.R')"
```

Expected: FAIL because `dof_gls_style()` does not exist.

- [ ] **Step 3: Add optional `mmrm` tests before implementation**

Append tests that run when `mmrm` is installed and compare both public helpers with the independently extracted fixed-effect and covariance-parameter counts:

```r
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
  expect_error(dof_mmrm(stats::lm(mpg ~ wt, data = mtcars)),
               "requires an mmrm model")
})
```

Run the same focused test command. Expected: the new public functions are still missing; when `mmrm` is unavailable, its fitted-model test is explicitly skipped rather than causing an unrelated setup failure.

- [ ] **Step 4: Implement and document the minimal parameter-count functions**

Create `R/model-dof.R`. Use `inherits()` for class dispatch, `isTRUE()` for the fixed-sigma flag, and delegate the `mmrm` branch of `dof_gls_style()` to `dof_mmrm()`:

```r
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

  stop("Unsupported model class for degrees-of-freedom computation.",
       call. = FALSE)
}
```

Add `mmrm` under `Suggests` in `DESCRIPTION`; keeping it optional allows the
existing `nlme` functionality to install independently while the two helpers
give a clear runtime error if an `mmrm` object is used without its package.

- [ ] **Step 5: Run focused tests and verify green**

Run the focused test command from Step 2. Expected: the GLS tests pass, the
unsupported-model tests pass, and the real `mmrm` test either passes or is
reported as skipped because `mmrm` is not installed.

- [ ] **Step 6: Generate documentation and rerun package tests**

Run:

```powershell
$env:LC_ALL='C'; $env:LANG='C'; Rscript -e "roxygen2::roxygenise(); testthat::test_local('.')"
```

Expected: `NAMESPACE` exports both helpers, both `.Rd` files are generated, and
the package tests pass with only the conditional `mmrm` skip when applicable.

- [ ] **Step 7: Commit the completed parameter-count unit**

Review `git diff --check` and the focused diff. Stage only `DESCRIPTION`,
`NAMESPACE`, `R/model-dof.R`, the two generated manuals, and
`tests/testthat/test-model-dof.R`; do not stage `instantSandwich.Rproj`.

```powershell
git commit -m "Add model parameter count helpers"
```

### Task 2: Leave-one-subject-out Cook's distance

**Files:**
- Create: `tests/testthat/test-subject-cooks-distance.R`
- Create: `R/subject-cooks-distance.R`
- Modify (generated): `NAMESPACE`
- Create (generated): `man/subject_cooks_distance.Rd`

- [ ] **Step 1: Write failing tests and reusable fixtures**

Create `tests/testthat/test-subject-cooks-distance.R` with one small dataset
containing subjects from both sexes, helpers that fit `lme` and `gls` models,
and an independent one-subject calculation:

```r
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
  full_coef <- if (inherits(fit, "lme")) nlme::fixef(fit) else stats::coef(fit)
  reduced_data <- droplevels(data[as.character(data$Subject) != subject, ])
  reduced_fit <- update(fit, data = reduced_data)
  reduced_coef <- if (inherits(fit, "lme")) {
    nlme::fixef(reduced_fit)
  } else {
    stats::coef(reduced_fit)
  }
  delta <- reduced_coef - full_coef
  fixed_frame <- stats::model.frame(stats::formula(fit), data = data,
                                    na.action = stats::na.omit)
  rank <- qr(stats::model.matrix(stats::formula(fit), fixed_frame))$rank

  as.numeric(crossprod(delta, solve(stats::vcov(fit), delta)) / rank)
}
```

Add tests for the desired result schema, inferred IDs, both model classes, and
the independent analytical comparison:

```r
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
```

- [ ] **Step 2: Add validation and failed-refit tests**

Append tests for an unsupported `lm`, a missing ID column, and a dataset where
one subject is the only subject at a fixed-effect factor level. Removing that
subject must create one retained failure row while other subjects remain usable:

```r
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
```

- [ ] **Step 3: Run focused tests and verify the missing-function failure**

Run:

```powershell
$env:LC_ALL='C'; $env:LANG='C'; Rscript -e "testthat::test_file('tests/testthat/test-subject-cooks-distance.R')"
```

Expected: FAIL because `subject_cooks_distance()` does not exist.

- [ ] **Step 4: Implement internal extraction and validation helpers**

Create `R/subject-cooks-distance.R`. Add 1-5 line comments above each internal
function, per the repository instructions:

```r
# Extract fixed-effect coefficients through the interface for each nlme class.
subject_coefficients <- function(fit) {
  if (inherits(fit, "lme")) nlme::fixef(fit) else stats::coef(fit)
}

# Resolve and validate one subject column, inferring it from the fit if omitted.
subject_id_column <- function(fit, data, id) {
  if (is.null(id)) {
    id_variables <- all.vars(nlme::getGroupsFormula(fit))
    if (length(id_variables) != 1L) {
      stop("Supply `id` as one column name for this grouping structure.",
           call. = FALSE)
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
    stop(sprintf("Subject column `%s` must not contain missing values.", id),
         call. = FALSE)
  }
  id
}

# Refit one model after removing a subject, optionally overriding fit control.
refit_without_subject <- function(fit, data, id, subject, control) {
  keep <- as.character(data[[id]]) != subject
  reduced_data <- droplevels(data[keep, , drop = FALSE])
  if (is.null(control)) {
    update(fit, data = reduced_data)
  } else {
    update(fit, data = reduced_data, control = control)
  }
}
```

- [ ] **Step 5: Implement the exported Cook's-distance function and example**

Add the exported function below the helpers. Keep refit errors local to each
subject and use `solve(V, delta)` rather than forming a numerical inverse:

```r
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
    stop("`subject_cooks_distance()` requires an nlme::lme or nlme::gls model.",
         call. = FALSE)
  }
  id <- subject_id_column(fit, data, id)
  subjects <- unique(as.character(nlme::getGroups(fit)))
  if (!all(subjects %in% as.character(data[[id]]))) {
    stop("Not all fitted subjects occur in the supplied `data`.", call. = FALSE)
  }

  full_coef <- subject_coefficients(fit)
  full_covariance <- stats::vcov(fit)
  fixed_formula <- stats::formula(fit)
  fixed_frame <- stats::model.frame(
    fixed_formula,
    data = data,
    na.action = stats::na.omit
  )
  fixed_rank <- qr(stats::model.matrix(fixed_formula, fixed_frame))$rank

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
```

- [ ] **Step 6: Run focused tests and verify green**

Run the focused test command from Step 3. Expected: all `lme`, `gls`, analytical,
validation, and failed-refit tests pass without warnings.

- [ ] **Step 7: Generate documentation and run all package tests**

Run:

```powershell
$env:LC_ALL='C'; $env:LANG='C'; Rscript -e "roxygen2::roxygenise(); testthat::test_local('.')"
```

Expected: `NAMESPACE` exports `subject_cooks_distance`, its `.Rd` file is
generated, and the complete package test suite passes.

- [ ] **Step 8: Commit the completed influence unit**

Review `git diff --check` and the focused diff. Stage only `NAMESPACE`,
`R/subject-cooks-distance.R`, `man/subject_cooks_distance.Rd`, and
`tests/testthat/test-subject-cooks-distance.R`; leave `instantSandwich.Rproj`
untracked.

```powershell
git commit -m "Add subject-level Cook's distance diagnostics"
```

### Task 3: Package-level verification

**Files:**
- Review all files changed in Tasks 1 and 2.

- [ ] **Step 1: Run clean documentation and test verification**

Run:

```powershell
$env:LC_ALL='C'; $env:LANG='C'; Rscript -e "roxygen2::roxygenise(); testthat::test_local('.')"
```

Expected: all package tests pass; conditional tests are skipped only when the
optional `mmrm` package is unavailable.

- [ ] **Step 2: Build the source package**

From the parent directory, run:

```powershell
$env:LC_ALL='C'; $env:LANG='C'; R CMD build instantSandwich
```

Expected: exit code 0 and a newly built `instantSandwich_0.1.0.tar.gz`.

- [ ] **Step 3: Run R CMD check**

Run:

```powershell
$env:LC_ALL='C'; $env:LANG='C'; $env:_R_CHECK_FORCE_SUGGESTS_='false'; R CMD check --no-manual instantSandwich_0.1.0.tar.gz
```

Expected: exit code 0 with no errors or warnings. If `mmrm` is unavailable, R
may report only the standard note that a suggested package was not available;
install `mmrm` and rerun when a zero-note result is required.

- [ ] **Step 4: Review repository state and commits**

Run `git status --short --branch`, `git diff --check HEAD~2..HEAD`, and
`git show --stat --oneline HEAD~1..HEAD`. Confirm that only requested package
code, tests, dependency metadata, generated documentation, and planning records
were committed, and that no push occurred.
