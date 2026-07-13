# F-Test Example and `hello()` Removal Design

## Goal

Remove the placeholder `hello()` API and extend the `robust.cov()` documentation
example with the coefficient tests and parallelism F test demonstrated in the
BIOS 767 `GLM-MVN.Rmd` course example.

## Package API cleanup

Delete `R/hello.R` and its generated `man/hello.Rd` help page. Regenerating the
roxygen2 namespace will remove `export(hello)` while retaining
`export(robust.cov)`.

## Example model and coefficient tests

Fit the `nlme::Orthodont` data with the course-style model
`distance ~ 0 + Sex + age + Sex:age`. This produces two group intercepts, a
common age slope, and a sex-by-age interaction in the same structural order as
the linked course model.

For coefficient-level inference, use between-within degrees of freedom matching
the course calculation:

```r
c(
  rep(n_subjects - 2, 2),
  rep(n_observations - n_subjects - 2, 2)
)
```

Retain the requested conventional two-sided p-values for both model-based and
sandwich-robust t statistics.

## Parallelism F test

Use `L <- matrix(c(0, 0, 0, 1), nrow = 1)` to test whether the fourth
coefficient, the sex-by-age interaction, equals zero. Compute the Wald F
statistic from `covariance$Sig.model`, with numerator degrees of freedom 1 and
denominator degrees of freedom
`n_observations - n_subjects - 2`, matching `GLM-MVN.Rmd`. Report the statistic,
p-value, and both degrees of freedom in a data frame.

## Tests and verification

Add a namespace test confirming that `hello` is no longer exported. Generate
documentation, run package tests and examples, and run `R CMD check
--no-manual`. Review the diff and commit the follow-up without pushing. The
untracked RStudio project file remains outside the change.
