#!/usr/bin/env Rscript
# =====================================================================
# A/B (S-05): is the disjoint/non-disjoint init asymmetry a problem?
# ---------------------------------------------------------------------
# KANBAN S-05 (orthoMTL.R:58): W_0 is initialised differently by mode --
# zeros when disjoint = FALSE, random positive values when disjoint =
# TRUE -- which creates "asymmetric seed-dependence". This script asks
# whether that asymmetry is a defect or a necessity, by forcing each
# init into each mode (via the W_0 argument -- no solver edit) and
# measuring seed-to-seed spread of the solution.
#
# Four cells, each fitted over several seeds on identical data:
#   non-disjoint + zero  (shipped default)  -> expect deterministic
#   non-disjoint + random                   -> does randomness help or just add noise?
#   disjoint     + random (shipped default) -> seed-dependent by construction
#   disjoint     + zero                      -> KEY TEST
#
# The disjoint projection proj_disjoint() zeros every entry where v <= 0;
# with a zero start v = 0 everywhere, so the iterate is pinned at 0 -- a
# fixed point. We expect disjoint+zero to be DEGENERATE (no learning),
# which is precisely why the shipped code seeds disjoint with random
# positive values. The asymmetry would then be necessary, not a bug.
#
# Run (against the INSTALLED package):
#   "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" validation/ab_s05_init.R
# =====================================================================

suppressMessages(library(orthoMTL))

set.seed(31)
n <- 150L; p <- 12L; n_tasks <- 3L
lambda <- 1e-2
seeds  <- 1:6

# Fixed data with a roughly disjoint-support ground truth (each task uses a
# different feature block), so disjoint mode has real structure to find.
B_true <- matrix(0, p, n_tasks)
blk <- split(seq_len(p), cut(seq_len(p), n_tasks, labels = FALSE))
for (t in seq_len(n_tasks)) B_true[blk[[t]], t] <- 1
X <- matrix(rnorm(n * p), n, p)
Y <- X %*% B_true + matrix(rnorm(n * n_tasks), n, n_tasks) * 0.3

maxabs <- function(A) max(abs(A))

# Final objective + coefficient matrix for one (disjoint, init) over a seed.
run_cell <- function(disjoint, init, seed) {
  W0 <- switch(init,
    default = NULL,                                   # solver chooses by mode
    zero    = matrix(0, p, n_tasks),
    random  = { set.seed(seed); matrix(abs(rnorm(p * n_tasks)), p, n_tasks) }
  )
  fit <- orthoMTL(X, Y, lambda = lambda, disjoint = disjoint, W_0 = W0,
                  seed = seed, stop_no_improve = 300L, max_iter = 1e5,
                  verbose = 0)
  list(obj = fit$obj, B = fit$B)
}

# Spread of the solution across seeds: worst-case pairwise max|B_i - B_j|.
cell_summary <- function(disjoint, init) {
  fits <- lapply(seeds, function(s) run_cell(disjoint, init, s))
  objs <- vapply(fits, `[[`, numeric(1), "obj")
  Bs   <- lapply(fits, `[[`, "B")
  spread <- 0
  for (a in seq_along(Bs)) for (b in seq_along(Bs)) {
    spread <- max(spread, maxabs(Bs[[a]] - Bs[[b]]))
  }
  data.frame(
    mode    = if (disjoint) "disjoint" else "non-disjoint",
    init    = init,
    obj_mean = mean(objs),
    obj_sd   = stats::sd(objs),
    B_spread = spread,                  # 0 == identical across seeds
    Bnorm    = sqrt(mean(vapply(Bs, function(B) sum(B^2), numeric(1))))
  )
}

tab <- rbind(
  cell_summary(FALSE, "default"),   # non-disjoint zero (shipped)
  cell_summary(FALSE, "random"),    # non-disjoint random
  cell_summary(TRUE,  "default"),   # disjoint random (shipped)
  cell_summary(TRUE,  "zero")       # disjoint zero  (the key test)
)

cat("A/B S-05: initialisation asymmetry  (n =", n, ", p =", p,
    ", tasks =", n_tasks, ", lambda =", lambda, ", seeds =",
    length(seeds), ")\n\n")
print(format(tab, digits = 4), row.names = FALSE)

# Pull out the four cells for the verdict.
nd_zero  <- tab[tab$mode == "non-disjoint" & tab$init == "default", ]
nd_rand  <- tab[tab$mode == "non-disjoint" & tab$init == "random",  ]
dj_rand  <- tab[tab$mode == "disjoint"     & tab$init == "default", ]
dj_zero  <- tab[tab$mode == "disjoint"     & tab$init == "zero",    ]

# Is disjoint+zero degenerate? Compare its objective to the trivial all-zero
# objective 0.5*mean(Y^2) and check the coefficients never left ~0.
trivial_obj <- 0.5 * sum(Y^2) / nrow(X)
dj_zero_degenerate <- dj_zero$Bnorm < 1e-6 &&
  abs(dj_zero$obj_mean - trivial_obj) < 1e-6

TOL <- 1e-6
cat("\nReadout:\n")
cat(sprintf("  non-disjoint + zero (shipped)  : seed-spread = %.2e  %s\n",
            nd_zero$B_spread,
            if (nd_zero$B_spread < TOL) "[deterministic]" else "[seed-dependent]"))
cat(sprintf("  non-disjoint + random          : seed-spread = %.2e ; obj %+.2e vs zero\n",
            nd_rand$B_spread, nd_rand$obj_mean - nd_zero$obj_mean))
cat(sprintf("  disjoint     + random (shipped): seed-spread = %.2e  [seed-dependent]\n",
            dj_rand$B_spread))
cat(sprintf("  disjoint     + zero            : Bnorm = %.2e ; obj=%.4f vs trivial %.4f  %s\n",
            dj_zero$Bnorm, dj_zero$obj_mean, trivial_obj,
            if (dj_zero_degenerate) "[DEGENERATE: pinned at 0]" else "[learned]"))

cat("\nVERDICT:\n")
nd_random_no_help <- nd_rand$obj_mean >= nd_zero$obj_mean - TOL
if (nd_zero$B_spread < TOL && dj_zero_degenerate && nd_random_no_help) {
  cat("  The init asymmetry is NECESSARY, not a defect:\n",
      " - disjoint + zero is a FIXED POINT of proj_disjoint() (v<=0 zeros every",
      "\n   entry), so it never learns -- the solution stays pinned at 0. The",
      "random\n   positive seed exists precisely to escape this.\n",
      " - non-disjoint + zero is already DETERMINISTIC (zero seed-spread) and",
      "random\n   init only adds seed noise without lowering the objective.\n",
      " So each mode uses the only sensible init. CONCLUSION: keep the",
      "asymmetry;\n  the seed argument already lets users explore disjoint",
      "restarts if needed.\n")
} else {
  cat("  Unexpected pattern -- inspect the table. Either disjoint+zero learned",
      "\n  (projection fixed-point not hit) or random init changed non-disjoint",
      "\n  results materially. Re-examine S-05 before concluding.\n")
}

if (any(!is.finite(unlist(tab[, c("obj_mean","obj_sd","B_spread","Bnorm")])))) {
  cat("\nWARNING: non-finite values; inspect before trusting.\n")
  quit(status = 1L)
}
quit(status = 0L)
