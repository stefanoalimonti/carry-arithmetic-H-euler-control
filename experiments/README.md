# Experiments — Paper H

| Script | Description | Referenced in |
|--------|-------------|---------------|
| `H01_carry_zeta_control.py` | Full carry product Z_carry(s) on the critical line | §3.1 |
| `H02_carry_vs_euler.py` | Side-by-side: carry product vs Euler product at Riemann zeros | §3.2, §3.3 |
| `H03_corrected_carry_product.py` | Isolated carry correction ∏R(l, s) smooth envelope | §4.2 |
| `H04_convergence_rate.py` | Convergence rate of carry correction as L increases | App. B |
| `H05_shape_correlation.py` | KS test: ∏R extrema vs Riemann zero locations | §3.4 |
| `H06_complex_spectral.py` | Complex spectral determinant analysis | §4.1 |

## Parameter Scaling

By default, scripts use reduced parameters (D = 16–24, N = 100–2000, L ≤ 47)
for quick verification (~1 min).  The paper's quantitative results (Appendix A:
D = 40, N = 10 000, L = 97) require the full parameter set.

**H01** accepts command-line flags to reproduce paper-scale results:

```bash
python H01_carry_zeta_control.py              # quick mode (default)
python H01_carry_zeta_control.py --paper      # paper parameters (slow)
python H01_carry_zeta_control.py -D 40 -N 10000 -L 97   # equivalent
```

Other scripts (H02–H06) still use hard-coded reduced parameters.  To
reproduce paper-scale results for those, increase `d_test`, `n_semiprimes`,
and the prime list inside each script.

Qualitative conclusions (∏R is smooth, Euler product drives zero-detection)
hold at all tested scales; quantitative tables in the paper require
`--paper`-level runs.

## Requirements

- Python >= 3.8, NumPy, SymPy, mpmath
