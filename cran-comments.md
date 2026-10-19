# cran-comments.md

## Submission type

Patch release, orthoMTL 0.1.1 (update of the CRAN release 0.1.0).

* Fixes a correctness bug in `nnmaxheap_C()` (the projection onto
  non-negative, non-increasing vectors), reported externally in
  https://github.com/Novartis/orthoMTL/issues/12. In 0.1.0 it returned
  inflated values for some inputs, which affected `predict()` on survival
  fits and the scoring in `cv_orthoMTL(survival = TRUE)`. Regression and
  classification fits are not affected. The function is now tested against
  an isotonic-regression oracle (`stats::isoreg`). Details are in `NEWS.md`.
* Adds `URL` and `BugReports` fields to `DESCRIPTION`.

The maintainer address is unchanged.

## Test environments

* Local: Windows 11, R 4.5.3
* GitHub Actions:
  * macOS-latest (R release)
  * windows-latest (R release)
  * ubuntu-latest (R devel, release, oldrel-1)

## R CMD check results

0 errors | 0 warnings | 1 note

* "checking for future file timestamps: unable to verify current time".
  This comes from the local check machine's clock lookup, not from the
  package.

## Downstream dependencies

There are currently no downstream dependencies for this package.
