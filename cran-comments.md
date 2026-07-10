# cran-comments.md

## Submission type

Initial CRAN submission of orthoMTL 0.1.0.

## Test environments

* Local: Windows 11, R 4.5.2
* GitHub Actions:
  * macOS-latest (R release)
  * windows-latest (R release)
  * ubuntu-latest (R devel, release, oldrel-1)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

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
