# `robust.cov()` Design

## Goal

Add an exported `robust.cov(u)` function to `instantSandwich` that preserves the
interface and results used in `767_homework/hw3_code/part2and3.r`, without
modifying the homework repository.

## Public interface

`robust.cov(u)` accepts a fitted `nlme` model for which `formula()`,
`model.frame()`, `model.matrix()`, `residuals()`, `nlme::getData()`,
`nlme::getGroups()`, and `nlme::getVarCov()` provide the model components needed
for clustered sandwich covariance estimation.

It returns a named list with the same four elements as the homework function:

- `Sig.model`: model-based coefficient covariance matrix.
- `se.model`: model-based coefficient standard errors.
- `Sig.robust`: sandwich-robust coefficient covariance matrix.
- `se.robust`: sandwich-robust coefficient standard errors.

## Computation

For each subject, the function extracts the subject-specific design matrix,
response residuals, and marginal covariance matrix. It accumulates
`X' V^-1 X` and the cluster-level sandwich meat, then computes the analytical
model-based and robust covariance matrices. Character group identifiers are
used when requesting marginal covariance matrices to avoid factor-indexing
mismatches.

The implementation will preserve the homework calculation and use
namespace-qualified `nlme` calls where appropriate. Matrix singularities and
unsupported model objects will surface through informative errors from the
underlying R or `nlme` operations rather than introducing a new incompatible
validation layer.

## Documentation example

The roxygen2 example will fit a small `nlme` model, call `robust.cov()`, and show
how the homework consumes both `se.model` and `se.robust` for coefficient-level
inference. It will calculate conventional two-sided p-values as
`2 * pt(-abs(t_value), df)`, correcting the one-tail expression in the homework
script while leaving that script unchanged.

## Package integration

The package will export `robust.cov`, declare `nlme` as an imported dependency,
and replace the placeholder package metadata needed for a valid package check.
The default `hello()` scaffold is outside the feature and will be left intact
unless package checks require otherwise.

## Tests and verification

Tests will first establish the missing-function failure, then verify the public
return structure and compare the function's matrices and standard errors with
an independently assembled subject-level sandwich calculation on a small fitted
`nlme` model. Documentation generation, package tests, and `R CMD check` will be
run before completion. The final diff will be reviewed, and implementation
changes will be committed without pushing.
