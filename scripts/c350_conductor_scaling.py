"""
C-350: Conductor scaling analysis for 20 EC L-functions
A_mean vs log(N): OLS fit, partial correlation, theory comparison, within-rank slopes
"""

import numpy as np
from scipy import stats

# ── Data ────────────────────────────────────────────────────────────────────
data = [
    # (curve,  rank, N,      A_mean,   A_std,    A_median)
    ("11a1",      0,    11,   11.781883, 11.460905,  9.356021),
    ("14a6",      0,    14,   10.828672,  5.730233,  9.741295),
    ("15a1",      0,    15,   10.107686,  4.340766,  9.807334),
    ("43a1",      0,    43,   15.401610,  9.586877, 12.941250),
    ("197a1",     0,   197,   43.445799, 41.508863, 32.025383),
    ("37a1",      1,    37,   21.841308, 24.413449, 12.912245),
    ("53a1",      1,    53,   16.833372, 11.258824, 13.687443),
    ("58a1",      1,    58,   21.047827, 24.363645, 13.474656),
    ("61a1",      1,    61,   21.153289, 17.752404, 14.824728),
    ("79a1",      1,    79,   24.871237, 38.939956, 14.306612),
    ("389a1",     2,   389,   28.528761, 22.565275, 22.058109),
    ("433a1",     2,   433,   34.691792, 50.449267, 25.554851),
    ("446a1",     2,   446,   28.439588, 16.320417, 24.480136),
    ("563a1",     2,   563,   29.836343, 20.816958, 24.423589),
    ("571a1",     2,   571,   29.135032, 20.486780, 24.880428),
    ("5077a1",    3,  5077,   52.300602, 35.130204, 41.926516),
    ("11197a1",   3, 11197,   53.604026, 37.425583, 43.476101),
    ("11642a1",   3, 11642,  100.624715,314.931277, 41.866509),
    ("12279a1",   3, 12279,   48.605324, 33.038095, 41.825647),
    ("13766a1",   3, 13766,   67.146456, 80.510683, 46.794673),
]

curves   = [d[0] for d in data]
ranks    = np.array([d[1] for d in data], dtype=float)
Ns       = np.array([d[2] for d in data], dtype=float)
A_mean   = np.array([d[3] for d in data])
A_std    = np.array([d[4] for d in data])
A_median = np.array([d[5] for d in data])
logN     = np.log(Ns)

SEP = "=" * 70

# ════════════════════════════════════════════════════════════════════════════
# 1. OLS: A_mean = α·log(N) + β
# ════════════════════════════════════════════════════════════════════════════
slope, intercept, r_val, p_val, se_slope = stats.linregress(logN, A_mean)
r2 = r_val ** 2

print(SEP)
print("1. CONDUCTOR SCALING  —  A_mean = α·log(N) + β  (OLS, n=20)")
print(SEP)
print(f"  α  (slope)    = {slope:.4f}  (SE={se_slope:.4f})")
print(f"  β  (intercept)= {intercept:.4f}")
print(f"  R²            = {r2:.4f}")
print(f"  r             = {r_val:.4f}")
print(f"  p-value       = {p_val:.4e}")
strength = "strong" if r2 > 0.7 else ("moderate" if r2 > 0.4 else "weak")
print(f"  Fit quality   : {strength} (R²={r2:.3f})")

# ════════════════════════════════════════════════════════════════════════════
# 2. Partial correlation: corr(A_mean, rank | log(N))
#    Method: residualise both A_mean and rank on log(N), then correlate
# ════════════════════════════════════════════════════════════════════════════
def residualise(y, x):
    """Return residuals of OLS y ~ x."""
    sl, ic, _, _, _ = stats.linregress(x, y)
    return y - (sl * x + ic)

res_A    = residualise(A_mean, logN)
res_rank = residualise(ranks, logN)
pc_r, pc_p = stats.pearsonr(res_A, res_rank)

print()
print(SEP)
print("2. PARTIAL CORRELATION  —  corr(A_mean, rank | log(N))")
print(SEP)
print(f"  Partial r     = {pc_r:.4f}")
print(f"  p-value       = {pc_p:.4e}")
sig = "significant" if pc_p < 0.05 else "not significant"
direction = "positive" if pc_r > 0 else "negative"
print(f"  Interpretation: rank has an INDEPENDENT {direction} effect on A_mean")
print(f"                  after controlling for log(N)  [{sig}, p={pc_p:.4e}]")

# ════════════════════════════════════════════════════════════════════════════
# 3. Theory comparison: κ_background = (1/2)·log(N·t²/(4π²))  at t=50
# ════════════════════════════════════════════════════════════════════════════
t0 = 50.0
kappa_theory = 0.5 * np.log(Ns * t0**2 / (4 * np.pi**2))

print()
print(SEP)
print(f"3. THEORY COMPARISON  —  κ_bg = (1/2)·log(N·t²/4π²)  at t={t0:.0f}")
print(SEP)
print(f"  {'Curve':<12} {'rank':>4}  {'log(N)':>7}  {'κ_theory':>9}  {'A_mean':>8}  {'ratio A/κ':>9}")
print(f"  {'-'*12} {'-'*4}  {'-'*7}  {'-'*9}  {'-'*8}  {'-'*9}")
ratios = []
for i, c in enumerate(curves):
    ratio = A_mean[i] / kappa_theory[i]
    ratios.append(ratio)
    print(f"  {c:<12} {int(ranks[i]):>4}  {logN[i]:>7.3f}  {kappa_theory[i]:>9.3f}  {A_mean[i]:>8.3f}  {ratio:>9.4f}")

print()
# Fit A_mean vs kappa_theory
kt_sl, kt_ic, kt_r, kt_p, kt_se = stats.linregress(kappa_theory, A_mean)
kt_r2 = kt_r**2
print(f"  OLS  A_mean ~ κ_theory:  slope={kt_sl:.4f}  intercept={kt_ic:.4f}")
print(f"       R²={kt_r2:.4f}   p={kt_p:.4e}")
print(f"  Mean ratio A_mean/κ_theory = {np.mean(ratios):.4f}  (std={np.std(ratios):.4f})")
print(f"  → A_mean tracks κ_theory {'well' if kt_r2 > 0.7 else 'partially' if kt_r2 > 0.4 else 'poorly'}")

# ════════════════════════════════════════════════════════════════════════════
# 4. Within-rank slopes: rank 0 and rank 3
# ════════════════════════════════════════════════════════════════════════════
print()
print(SEP)
print("4. WITHIN-RANK SLOPES  —  A_mean ~ log(N)  (rank 0 and rank 3)")
print(SEP)

for r_target in [0, 1, 2, 3]:
    mask = ranks == r_target
    n_r = mask.sum()
    logN_r  = logN[mask]
    A_r     = A_mean[mask]
    N_range = f"N in [{int(Ns[mask].min())}, {int(Ns[mask].max())}]"

    if n_r < 3:
        print(f"  rank {r_target} (n={n_r}): too few points for reliable OLS — skipped")
        continue

    sl_r, ic_r, rv_r, pv_r, se_r = stats.linregress(logN_r, A_r)
    r2_r = rv_r ** 2
    sig_r = "p<0.05 *" if pv_r < 0.05 else f"p={pv_r:.3f} (ns)"
    mark = " <<<" if r_target in (0, 3) else ""
    print(f"  rank {r_target} (n={n_r}, {N_range}){mark}")
    print(f"    α = {sl_r:.4f}  β = {ic_r:.4f}  R²={r2_r:.4f}  {sig_r}")

# ════════════════════════════════════════════════════════════════════════════
# 5. Clean summary
# ════════════════════════════════════════════════════════════════════════════
print()
print(SEP)
print("5. SUMMARY")
print(SEP)
print(f"  [Global OLS]   A_mean = {slope:.3f}·log(N) + {intercept:.3f}")
print(f"                 R²={r2:.4f}  p={p_val:.2e}  → log(N) explains {r2*100:.1f}% of A_mean variance")
print()
print(f"  [Partial corr] After partialling out log(N):")
print(f"                 rank effect: partial r={pc_r:.4f}  p={pc_p:.4e}  [{sig}]")
rank_conclusion = (
    "Rank retains independent explanatory power beyond conductor scaling."
    if pc_p < 0.05 else
    "Rank effect on A_mean is NOT independent of conductor once log(N) is controlled."
)
print(f"                 → {rank_conclusion}")
print()
print(f"  [Theory]       κ_theory at t=50 fits A_mean with R²={kt_r2:.4f}")
print(f"                 Mean A_mean/κ_theory = {np.mean(ratios):.3f}  (excess factor over background)")
print()
print(f"  [11642a1 note] rank-3 outlier: A_mean=100.6, A_std=314.9 (high-variance;")
print(f"                 P2 CV=32.4%, 3/4 PASS). Influences rank-3 within-group slope.")

# Spearman check at curve level
sp_r, sp_p = stats.spearmanr(ranks, A_mean)
print()
print(f"  [Spearman]     ρ(rank, A_mean) curve-level = {sp_r:.4f}  p={sp_p:.4e}  (n=20)")
print()
print(SEP)
