# F-Test Example and `hello()` Removal Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the placeholder `hello()` export and extend the `robust.cov()` example with course-matched coefficient and parallelism F tests.

**Architecture:** Keep `robust.cov()` unchanged and express the inferential workflow as a runnable roxygen2 example. Use the four-coefficient Orthodont model to reproduce the linked course model's coefficient ordering, between-within degrees of freedom, and fourth-coefficient parallelism contrast.

**Tech Stack:** R 4.5, `nlme`, roxygen2, testthat edition 3

---

### Task 1: Remove the placeholder `hello()` API

**Files:**
- Modify: `tests/testthat/test-robust-cov.R`
- Delete: `R/hello.R`
- Delete: `man/hello.Rd`
- Regenerate: `NAMESPACE`

- [ ] **Step 1: Add a failing namespace test**

Append this test to `tests/testthat/test-robust-cov.R`:

```r
test_that("the placeholder hello function is not exported", {
  expect_false("hello" %in% getNamespaceExports("instantSandwich"))
})
```

- [ ] **Step 2: Run the test and verify the RED state**

Run:

```powershell
$env:LC_ALL='C'; $env:LANG='C'; Rscript -e "devtools::test()"
```

Expected: one failure because `hello` remains exported.

- [ ] **Step 3: Delete the placeholder source and help files**

Delete `R/hello.R` and `man/hello.Rd`. Then run:

```powershell
$env:LC_ALL='C'; $env:LANG='C'; Rscript -e "roxygen2::roxygenise()"
```

Expected: `NAMESPACE` contains only `export(robust.cov)`.

- [ ] **Step 4: Run the test and verify the GREEN state**

Run:

```powershell
$env:LC_ALL='C'; $env:LANG='C'; Rscript -e "devtools::test()"
```

Expected: all estimator and namespace tests pass.

### Task 2: Add course-matched coefficient and F tests to the example

**Files:**
- Modify: `R/robust-cov.R`
- Regenerate: `man/robust.cov.Rd`

- [ ] **Step 1: Replace the runnable roxygen2 example**

Use this complete `@examples` body:

```r
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
```

- [ ] **Step 2: Regenerate documentation**

Run:

```powershell
$env:LC_ALL='C'; $env:LANG='C'; Rscript -e "roxygen2::roxygenise()"
```

Expected: `man/robust.cov.Rd` contains the four-coefficient model, two-sided
t-test p-values, and F test with denominator degrees of freedom
`n_observations - n_subjects - 2`.

- [ ] **Step 3: Run tests and examples**

Run:

```powershell
$env:LC_ALL='C'; $env:LANG='C'; Rscript -e "devtools::test(); devtools::run_examples()"
```

Expected: all tests and the updated example pass without errors.

### Task 3: Verify, review, and commit

**Files:**
- Review: `NAMESPACE`
- Review: `R/robust-cov.R`
- Review: `man/robust.cov.Rd`
- Review: `tests/testthat/test-robust-cov.R`
- Confirm deleted: `R/hello.R`
- Confirm deleted: `man/hello.Rd`

- [ ] **Step 1: Build and run the full package check**

Run:

```powershell
$env:LC_ALL='C'; $env:LANG='C'; R.exe CMD build .
$env:LC_ALL='C'; $env:LANG='C'; R.exe CMD check --no-manual instantSandwich_0.1.0.tar.gz
```

Expected: `R CMD check` ends with `Status: OK`.

- [ ] **Step 2: Review the complete diff**

Run:

```powershell
git diff --check
git diff -- NAMESPACE R/hello.R R/robust-cov.R man/hello.Rd man/robust.cov.Rd tests/testthat/test-robust-cov.R
```

Expected: only the approved API removal, example, generated documentation, and
namespace test are changed.

- [ ] **Step 3: Commit without pushing**

Run:

```powershell
git add NAMESPACE R/robust-cov.R man/robust.cov.Rd tests/testthat/test-robust-cov.R
git add -u R/hello.R man/hello.Rd
git commit -m "Document course-style sandwich F test"
```

Expected: one focused follow-up commit; no push is performed.
