# `robust.cov()` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an exported, documented, and tested `robust.cov(u)` function that is API-compatible with the function used in `767_homework`.

**Architecture:** Keep the estimator in one focused R source file. The function will recover model data from an `nlme` fit, accumulate model-based and cluster-level sandwich terms by subject, and return the exact four-element list expected by the homework. Tests will compare the result with an independently assembled analytical calculation, while the documentation example will demonstrate coefficient testing with two-sided p-values.

**Tech Stack:** R 4.5, `nlme`, roxygen2, testthat edition 3, base matrix algebra

---

### Task 1: Add package test infrastructure and a failing public-API test

**Files:**
- Modify: `DESCRIPTION`
- Create: `tests/testthat.R`
- Create: `tests/testthat/test-robust-cov.R`

- [ ] **Step 1: Declare runtime and test dependencies**

Update `DESCRIPTION`, preserving its existing version and `Authors@R`, with
these exact package fields:

```text
Title: Model-Based and Sandwich Covariance Estimation
Description: Computes model-based and cluster-robust sandwich covariance
    matrices and standard errors for fitted longitudinal models.
License: GPL-3
Imports:
    nlme
Suggests:
    testthat (>= 3.0.0)
Config/testthat/edition: 3
RoxygenNote: 7.3.3
```

- [ ] **Step 2: Add the standard testthat runner**

Create `tests/testthat.R`:

```r
library(testthat)
library(instantSandwich)

test_check("instantSandwich")
```

- [ ] **Step 3: Write the failing compatibility test**

Create `tests/testthat/test-robust-cov.R`:

```r
make_test_fit <- function() {
  nlme::gls(
    distance ~ age + Sex,
    data = nlme::Orthodont,
    correlation = nlme::corCompSymm(form = ~ 1 | Subject),
    method = "ML"
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
```

- [ ] **Step 4: Run the test and verify the RED state**

Run:

```powershell
Rscript -e "devtools::test()"
```

Expected: the test fails because `robust.cov` does not exist.

### Task 2: Implement the analytical covariance estimator

**Files:**
- Create: `R/robust-cov.R`
- Test: `tests/testthat/test-robust-cov.R`

- [ ] **Step 1: Add the minimal implementation**

Create `R/robust-cov.R` with this implementation beneath the complete roxygen2
documentation added in Task 3:

```r
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
```

- [ ] **Step 2: Run the compatibility test and verify the GREEN state**

Run:

```powershell
Rscript -e "devtools::test()"
```

Expected: the return-interface test passes.

- [ ] **Step 3: Add an independent numerical test**

Append this independent cluster-by-cluster calculation and numerical test:

```r
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

test_that("robust.cov matches an analytical cluster calculation", {
  fit <- make_test_fit()

  result <- robust.cov(fit)
  expected <- manual_covariance(fit)

  expect_equal(result$Sig.model, expected$Sig.model, tolerance = 1e-10)
  expect_equal(result$se.model, expected$se.model, tolerance = 1e-10)
  expect_equal(result$Sig.robust, expected$Sig.robust, tolerance = 1e-10)
  expect_equal(result$se.robust, expected$se.robust, tolerance = 1e-10)
})
```

- [ ] **Step 4: Run the numerical test and confirm it passes**

Run:

```powershell
Rscript -e "devtools::test()"
```

Expected: both tests pass with no failures or warnings.

### Task 3: Add exported documentation and the inference example

**Files:**
- Modify: `R/robust-cov.R`
- Regenerate: `NAMESPACE`
- Create: `man/robust.cov.Rd`

- [ ] **Step 1: Complete the roxygen2 documentation**

Add this complete roxygen2 block above `robust.cov()`:

```r
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
```

Immediately after `@examples`, include this runnable workflow:

```r
fit <- nlme::gls(
  distance ~ age + Sex,
  data = nlme::Orthodont,
  correlation = nlme::corCompSymm(form = ~ 1 | Subject),
  method = "ML"
)

covariance <- robust.cov(fit)
estimate <- stats::coef(fit)
degrees_freedom <- stats::df.residual(fit)

model_t <- estimate / covariance$se.model
robust_t <- estimate / covariance$se.robust

coefficient_tests <- data.frame(
  estimate = estimate,
  model_se = covariance$se.model,
  model_t = model_t,
  model_p = 2 * stats::pt(-abs(model_t), df = degrees_freedom),
  robust_se = covariance$se.robust,
  robust_t = robust_t,
  robust_p = 2 * stats::pt(-abs(robust_t), df = degrees_freedom)
)
coefficient_tests
```

- [ ] **Step 2: Generate package documentation**

Run:

```powershell
Rscript -e "roxygen2::roxygenise()"
```

Expected: `NAMESPACE` contains `export(robust.cov)` and
`man/robust.cov.Rd` contains the example with two-sided p-values.

- [ ] **Step 3: Run tests and documentation examples**

Run:

```powershell
Rscript -e "devtools::test(); devtools::run_examples()"
```

Expected: all tests and examples pass without errors.

### Task 4: Verify, review, and commit the package change

**Files:**
- Review: `DESCRIPTION`
- Review: `NAMESPACE`
- Review: `R/robust-cov.R`
- Review: `man/robust.cov.Rd`
- Review: `tests/testthat.R`
- Review: `tests/testthat/test-robust-cov.R`

- [ ] **Step 1: Run focused and full verification**

Run:

```powershell
Rscript -e "devtools::test()"
R CMD build .
R CMD check --no-manual instantSandwich_0.1.0.tar.gz
```

Expected: tests pass and `R CMD check` reports no errors or warnings. Any notes
caused by the pre-existing scaffold or local environment will be recorded
exactly rather than hidden.

- [ ] **Step 2: Confirm the homework repository is unchanged**

Run:

```powershell
git -C D:\RPackageDevelopment\homework\767_homework status --short
```

Expected: no changes attributable to this implementation.

- [ ] **Step 3: Review the focused diff**

Run:

```powershell
git diff --check
git diff -- DESCRIPTION NAMESPACE R/robust-cov.R man/robust.cov.Rd tests/testthat.R tests/testthat/test-robust-cov.R
```

Expected: no whitespace errors, no unrelated tracked changes, and exact API
compatibility with `robust.cov(u)` from the homework.

- [ ] **Step 4: Commit without pushing**

Run:

```powershell
git add DESCRIPTION NAMESPACE R/robust-cov.R man/robust.cov.Rd tests/testthat.R tests/testthat/test-robust-cov.R
git commit -m "Add homework-compatible sandwich covariance estimator"
```

Expected: one focused implementation commit; no push is performed.
