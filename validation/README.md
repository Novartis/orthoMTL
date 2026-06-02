# validation/

Confidence harnesses for the `orthoMTL` solver. These are **developer
scripts, not part of the package** (excluded via `.Rbuildignore`); they are
not run by `R CMD check` or the `testthat` suite. They exist to justify the
solver-correctness changes on branch `fix/logistic-and-alpha-migration`
(REG-01 logistic gradient, REG-02 overflow-safe softplus, REG-03 `enet`→
`alpha` migration), which deliberately change numerical results versus the
original 2014 ECML-PKDD communication.

Each script prints a table and exits non-zero on failure, so they can be run
in a pre-release check. Run them against the **installed** package:

```
R CMD INSTALL .
"C:\Program Files\R\R-4.5.3\bin\Rscript.exe" validation/parity_orthopen.R
"C:\Program Files\R\R-4.5.3\bin\Rscript.exe" validation/simulation_recovery.R
"C:\Program Files\R\R-4.5.3\bin\Rscript.exe" validation/old_vs_fixed_logistic.R
```

## The three legs

| Script | Question it answers | Needs |
|---|---|---|
| `parity_orthopen.R` | Is the fixed solver a faithful port? | `orthopen` **v1.1.0** installed |
| `old_vs_fixed_logistic.R` | Did the bug actually matter? | legacy `orthopen` source on disk |
| `simulation_recovery.R` | Does the survival path extract real signal? | `orthoMTL` only |

### 1. `parity_orthopen.R` — bit-for-bit vs the validated reference
Feeds identical inputs to `orthoMTL()` and `orthopen::orthopen()`
(github.com/kevinVervier/orthopen, **v1.1.0** — the independently validated
solver these fixes were ported *from*). The non-survival path is line-for-line
identical, so all shared cases (regression/orthogonal, `alpha` 0.5 and 1,
logistic, disjoint) agree to **`max|B_ref - B_new| = 0`**. Proves the fixes
are a faithful port, not a re-derivation that could carry new bugs.
Skips gracefully if `orthopen` is not installed.

### 2. `old_vs_fixed_logistic.R` — the bug, made visible
Compares against the **buggy pre-v1.0.1 ancestor** orthoMTL was forked from
(`orthopen/new/legacy_orthopen.R`: correct `{-1,+1}` logistic *loss* but
`{0,1}` cross-entropy *gradient* — the REG-01 mismatch). The logistic
objective is invariant under `(Y, W) -> (-Y, -W)`, so a correct minimiser
must satisfy `B(-Y) = -B(Y)`. The fix holds this (`0.000`); the buggy ancestor
violates it (`~0.62`) and has higher held-out logistic loss. Point the script
at the legacy source with the `ORTHOPEN_LEGACY` env var; skips gracefully if
absent.

### 3. `simulation_recovery.R` — the survival path parity can't reach
`survival = TRUE` masking is orthoMTL-only, so parity cannot test it. Uses
`simulate_mtl()`'s known ground-truth coefficients: over 10 seeds, paired
real-vs-null, support-recovery AUC is ~0.66 vs ~0.50 null (real > null in
10/10 seeds). The absolute held-out C-index is modest **by design** — see the
header comment: `simulate_mtl` effects are on the hazard scale while labels
encode progression-free, and `cindex_mtl` row-sums across tasks, collapsing
the time-varying structure (the documented M-01/02/03 metric caveats). Even
the sign-corrected true coefficients only reach ~0.5, so this is a metric/sim
ceiling, **not** a solver defect.
