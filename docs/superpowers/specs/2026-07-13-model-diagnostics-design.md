# Model Degrees of Freedom and Subject Influence Design

## Goal

Extend `instantSandwich` with exported helpers for counting fitted-model
parameters and computing leave-one-subject-out Cook's distances. The public
functions are `dof_mmrm()`, `dof_gls_style()`, and
`subject_cooks_distance()`.

## Public interfaces

`dof_mmrm(fit)` accepts an `mmrm` fit and returns the number of estimated fixed
effects plus the number of estimated covariance parameters.

`dof_gls_style(fit)` accepts either an `nlme::gls` fit or an `mmrm` fit. For a
`gls` fit, it counts regression coefficients, correlation or variance-function
parameters, and the residual standard deviation when it is estimated. For an
`mmrm` fit, it returns the same count as `dof_mmrm()`. Other model classes
produce an informative error.

`subject_cooks_distance(fit, data = nlme::getData(fit), id = NULL,
control = NULL)` accepts an `nlme::lme` or `nlme::gls` fit. `data` supplies the
rows used for leave-one-subject-out refits. `id` is an optional single character
column name; when omitted, the function infers the single grouping variable
from the fitted model. `control` optionally supplies an `nlme` control object
for every refit.

The Cook's-distance result is a data frame with one row per subject and four
columns: `id`, `cooks_distance`, `converged`, and `message`. A failed refit is
retained with `NA` distance, `FALSE` convergence, and the error message, so one
problematic subject does not discard all other results.

## Computation

For each subject, `subject_cooks_distance()` removes that subject's rows,
drops unused factor levels, and refits the original model. It uses
`nlme::fixef()` for an `lme` fit and `stats::coef()` for a `gls` fit. With
full-fit coefficient covariance matrix `V`, full coefficient vector `beta`,
leave-one-subject-out vector `beta[-i]`, and fixed-effect design rank `p`, it
computes the analytical measure

`(beta[-i] - beta)' solve(V) (beta[-i] - beta) / p`.

Refits run sequentially. This avoids adding parallel-processing dependencies
or introducing nondeterministic package examples. Nested or otherwise
ambiguous grouping structures require an explicit, single `id` column.

## Documentation and dependencies

Each exported function receives roxygen2 documentation and a runnable example.
Small built-in data subsets keep the repeated-refit example reasonably fast.
The `mmrm` package is declared as a package dependency because the public
degrees-of-freedom helpers call its model-component interface.

## Tests and verification

Tests will be written before implementation and will reuse the existing
`testthat` setup. They will cover parameter counts, supported and unsupported
model classes, subject-variable inference, both `lme` and `gls` refits, result
structure, and agreement with an independently calculated Cook's distance.
Roxygen output will be regenerated, package tests will be run, and a full
`R CMD check` will provide final verification. The implementation diff will be
reviewed and committed without pushing.
