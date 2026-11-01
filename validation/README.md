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

## S-* solver A/B studies

Behaviour-comparison scripts for the `S-*` KANBAN cards (loss/penalty
normalisation, convergence budget, init asymmetry, gradient schedule). Each
compares the **shipped** behaviour (A) against a **candidate** (B) on the
same data and prints a verdict; none modifies the package. All four are
driven entirely through existing arguments (`ab_s02` via the `schedule`
argument that the S-02 study motivated and that now ships in `orthoMTL()`).

| Script | Card | Finding |
|---|---|---|
| `ab_s01_normalisation.R` | S-01 | **Keep.** The per-observation loss normalisation gives *active, n-stable* shrinkage; dividing the penalty by `n` instead makes regularisation wash out as `n` grows. Standard glmnet convention, not a defect. |
| `ab_s04_convergence.R` | S-04 | **Defaults adequate.** Extra patience leaves the objective unchanged and the iteration cap never binds, even at p>n, λ=1e-4, k=16. The lever (if any) is `max_iter`, not `stop_no_improve`. |
| `ab_s05_init.R` | S-05 | **Asymmetry necessary.** `disjoint`+zero init is a fixed point of `proj_disjoint()` (stays pinned at 0); non-disjoint is deterministic at zero and random init only adds seed noise. |
| `ab_s02_gradient_schedule.R` | S-02 | **Implemented as `schedule=`.** `sqrt(i)` is correct but ~12–17× slower than `log`/`const` to the *same* optimum (`linear`=1/i stalls). Now selectable via `orthoMTL(schedule=)`; default `"sqrt"` preserves published numerics. |

All four run against the **installed** package only:

```
"C:\Program Files\R\R-4.5.3\bin\Rscript.exe" validation/ab_s01_normalisation.R
"C:\Program Files\R\R-4.5.3\bin\Rscript.exe" validation/ab_s04_convergence.R
"C:\Program Files\R\R-4.5.3\bin\Rscript.exe" validation/ab_s05_init.R
"C:\Program Files\R\R-4.5.3\bin\Rscript.exe" validation/ab_s02_gradient_schedule.R
```

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

## O-01: the atomic norm `Omega_O` (`penalty = "omega_O"`) versus the pairwise `Omega_K`

`ab_omega_o_vs_omega_k.R` (run against a library holding this branch's build; parts and env vars in
its header). Background and proofs: ai-research `audit/ortho_penalty_screening/README.md` §3.4–3.5
and the draft note `LaTeX/omega_o_two_summing/` (branch `dr_k`). For `p >= T`,
`Omega_O(W)^2 = max over correlation matrices Z of tr(W Z W^T)` (a T x T max-cut SDP, solved each
iteration by a Burer-Monteiro mixing method with a dual certificate, `R/omega_o.R`); it is convex
for every T and needs no `K`; it equals `penalty = "K"` with an all-ones `K` exactly when `W^T W` is
balanced. The default path (`penalty = "K"`) is byte-identical to 0.1.1 (checked on survival,
elastic-net/log-schedule and disjoint fits).

**Run of 2026-11-01** (12 seeds, 3-fold inner CV for every arm, schedule `sqrt`; parts R and F paired
task-by-task with ab_k01c's stored K arms through a refit-and-compare check that passed on all 144 R
tasks). Paired test-MSE difference vs identity K (ridge), negative = better:

| scenario | best K practice (C) | `Omega_O` (O) | lasso (L) |
|---|---|---|---|
| orthogonal, low noise | -0.075 | -0.039 | +0.027 |
| orthogonal, mid/high noise | ~0 | +0.02 (n.s./worse) | +0.14 ... +0.24 |
| sliding | -0.02 ... -0.08 | +0.08 ... +0.23 | -0.27 ... -0.98 |
| blocks | -0.03 ... -0.06 | +0.22 ... +0.58 | -0.27 ... -0.98 |
| shared | ~0 (CV retreats to ridge) | +0.37 ... +1.16 | -0.26 ... -1.05 |
| frustrated (new) | ~0 | +0.01 ... +0.11 | -0.09 ... +0.08 |
| rotating (new) | ~0 | +0.03 ... +0.42 | -0.17 ... -0.03 |
| orthogonal20 (T = 20, new) | -0.069 (noise 1) | -0.043 (noise 1) | ~0 ... +0.18 |

Survival (vignette protocol, lambda grid widened to 1e-5..1e-1): all arms indistinguishable
(C-index differences within +-0.003, n.s.; support AUC near chance at T = 4); `omega_O` CV is ~2x
faster (no `diag_val` grid).

**Reading.** Plain `Omega_O` helps only on genuinely orthogonal truth and is beaten there by a tuned
K; it loses badly wherever tasks share features, because `Omega_O^2 = ||W||_F^2 + (sharing term)`
ties shrinkage to the price of sharing with a fixed ratio (it cannot shrink noise features without
crushing shared ones). Lasso dominates every sparse scenario, a baseline ab_k01c never had.
Caveats: the subgradient solver may under-converge `Omega_O` more than `Omega_K` (a
column-generation solver is planned); every truth here is 5-of-p sparse, so penalty + alpha
(elastic net) is the fair next comparison. **Next:** part OG, `Omega_O + gamma ||W||_F^2`
(`gamma =` argument, convex for every gamma, ridge as gamma -> infinity), (gamma, lambda) by CV;
a 1-seed smoke run brought the deltas to -0.06 ... +0.07; the 12-seed run is in progress.
