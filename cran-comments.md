# cran-comments.md

## Submission type

Resubmission of orthoMTL 0.1.0.

This addresses the second round of reviewer feedback (Konstanze Lauseker):
functions must not call `set.seed()` to a hardcoded value. `orthoMTL()`,
`cv_orthoMTL()`, `bootstrap_orthoMTL()`, and `simulate_mtl()` all defaulted
`seed = 42` and called `set.seed()` unconditionally; `seed` now defaults to
`NULL` and `set.seed()` is only called when the caller supplies a value.

(First-round feedback, an invalid `+ file LICENSE` field, was fixed in the
prior resubmission.)

## Test environments

* Local: Windows 11, R 4.5.2
* GitHub Actions:
  * macOS-latest (R release)
  * windows-latest (R release)
  * ubuntu-latest (R devel, release, oldrel-1)

## R CMD check results

0 errors | 0 warnings | 0 notes

(A "unable to verify current time" NOTE appears in some local check runs due
to a network-reachability artifact of the time-verification service; it does
not reproduce on CRAN's own build infrastructure and is omitted above.)

## Downstream dependencies

There are currently no downstream dependencies for this package.

## Additional notes

* The package is the successor to `orthopen`
  (https://github.com/kevinVervier/orthopen), a small research package
  never released on CRAN. The relationship is acknowledged in
  `inst/LICENSE.note`.
* All exported functions have runnable examples; longer-running examples
  are wrapped in `\donttest{}`.
* Vignette optional Cox comparison sections use `eval =
  requireNamespace(...)` guards so the vignette builds even if
  `glmnet`/`survival` are not installed at build time.
