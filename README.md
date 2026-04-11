# Resonant Deep Learning — Neural Detection of Riemann Zeros

A framework that detects Riemann ξ-function zeros via gauge ODE residuals, then systematically diagnoses *why* the detector fails in specific regions. The diagnostic cascade produced a **Nyquist-type bandwidth bound** connecting neural architecture constraints to analytic number theory.

## Key Results

### 1. Zero Detection: 50/50 on t ∈ [100, 200]

A 3-layer resonant network learns a gauge ODE whose residual F₂(t) = Im{e^{-iφ}(L̂ - ψ)} vanishes at every Riemann zero. Tested on 50 zeros with 100% detection rate.

### 2. Bandwidth–Density Bridge Inequality (Proposition 8.5)

When a cluster of zeros at t ∈ [165, 173] resisted detection, a 10-experiment cascade identified the cause: the Fourier feature embedding's frequency ceiling.

**Necessary condition for representability:**

```
K ≥ W / (2·Δt)
```

where K = number of Fourier frequencies, W = training window width, Δt = minimum zero spacing.

Substituting the Riemann–von Mangoldt zero density N'(t) ~ (2π)⁻¹ log(t/2π):

```
K(t, W) ≥ W · log(t/2π) / (4π · c_L)
```

This links a **neural architecture parameter** (embedding dimension) to an **analytic number theory quantity** (zero density). Validated with 5-seed ensemble: cluster |F₂| improves 25% when K crosses the threshold.

### 3. ±π Discontinuity Bottleneck → Smooth Target Resolution

The isolation experiment (p = 1.18 × 10⁻¹⁵) proved the supervisory chain arg ξ(½+it) → network is sound, but the {0, π} step-function encoding creates an unrepresentable target at zeros. Replacing it with a smooth surrogate ξ_real/max|ξ_real| reduces the zero/non-zero MSE ratio from **18.2× to 0.03×** (99.8% gap elimination).

### 4. Four-Failure Decomposition Methodology

The cos² PQO training appeared to "disprove" the phase-zero correspondence. A systematic isolation cascade revealed **four independent failure modes** stacking into one negative result:

1. Target decoupling (synthetic target carries no zero information)
2. Loss-weight imbalance (residual loss dominates supervision)  
3. Arithmetic-mean readout (destroys ±π wrap signal)
4. Smooth-field expressivity (network can't represent step functions)

Each was isolated and resolved independently. The methodology transfers to any ML experiment where "it doesn't work" needs decomposition.

## Architecture

```
Input t → Fourier embedding [t_norm, sin(ω₁t), cos(ω₁t), ...] 
       → ComplexLinearEmbedding
       → ResonantBlock × 3 (Lie exponential phase propagation)
       → GaugeStateBuffer
       → F₂ residual readout
```

## Project Structure

```
gdl/
├── core/           GDL Blueprint (Bronstein et al. 2021)
├── domains/        5G + Time equivariant implementations
├── rdl/            Resonant Deep Learning
│   ├── models/     MasterResonantNetwork
│   ├── losses/     F₂ discriminant + PQO + curvature + target
│   ├── dynamics/   Gauge ODE (separator, damping, Heun integrator)
│   └── pipeline/   XiFeatureDataset (Fourier features from ξ-cache)
├── bridge/         RDL ↔ GDL 5G Gauge adapter
└── viz/            Phase, curvature, dashboard visualisation
```

### Key Experiment Scripts

| Script | What it does |
|--------|-------------|
| `basis_bandwidth_scaling.py` | Tests in_features ∈ {32, 64, 128, 256}, identifies bandwidth ceiling |
| `multiseed_if128.py` | 5-seed ensemble validation of bandwidth resolution |
| `optim_budget_scaling.py` | Separates bandwidth effect from over-parameterisation |
| `phase_isolation_experiment.py` | Isolation experiment (p = 1.18e-15) |
| `smooth_zero_crossing.py` | Smooth target confirmation (ratio 18× → 0.03×) |
| `ensemble_hard_zeros.py` | Multi-seed F₂ ensemble analysis |
| `overnight_exploration.py` | Large-scale zero exploration (14,441 zeros) |

## Quick Start

```bash
# Requirements: Python 3.10+, PyTorch, mpmath, scipy, numpy
pip install torch numpy scipy mpmath einops tqdm

# Run tests (115 tests)
python -m pytest tests/ -q

# Train zero detector on t ∈ [100, 200]
python scripts/run_xi_pipeline.py --optimizer adam --pqo-mode cos2 \
  --epochs 200 --t-min 100 --t-max 200 --hidden 64

# Reproduce bandwidth scaling experiment
python scripts/basis_bandwidth_scaling.py

# Reproduce smooth target experiment
python scripts/smooth_zero_crossing.py
```

## Paper

The full unified paper (44 pages EN / 41 pages KO) is available in `paper/`:
- Three-tier classification: 🟢 Validated / 🟡 Ideas / 🔴 Failure diagnosis
- Complete experimental cascade with honest retraction records
- All failures preserved with diagnostic chains

## Detection Results

| Range | Zeros | Detection | cos² PQO | Optimiser |
|-------|-------|-----------|----------|-----------|
| t ∈ [10, 50] | 10 | **10/10** | < 0.0002 | Adam |
| t ∈ [100, 200] | 50 | **50/50** | < 0.0011 | Adam |
| t ∈ [10, 50] | 10 | 0/10 | — | PGGD (failed) |

## What This Is Not

- **Not a proof of the Riemann Hypothesis.** Structural parallels are not logical proofs.
- **Not a better zero-finding algorithm.** Odlyzko–Schönhage is faster. This is a diagnostic lens.
- **Not claiming per-zero F₂ values are observables.** They are model-dependent (CV 46–106% across seeds). Only population-level statistics are reproducible.

## Citation

```
Kim, Young-Hu. "Resonant Deep Learning: Neural Detection of Riemann Zeros 
via Gauge ODE Residuals." Preprint, 2026.
```

## License

GPL-3.0 License. See [LICENSE](LICENSE).
