# carry-arithmetic-H-euler-control

**Carry Polynomials and the Partial Euler Product: A Control Experiment**

*Author: Stefano Alimonti* · [ORCID 0009-0009-1183-1698](https://orcid.org/0009-0009-1183-1698)

## Main Result

The carry product $Z_{\text{carry}}(s) = \prod_l \langle |\det(I - M_l / l^s)| \rangle$ has minima near the first nine Riemann zeros. A control experiment decomposes it as $Z_{\text{carry}} = Z_{\text{Euler}} \cdot \prod R$, isolating the carry correction. **Result:** the isolated carry correction $\prod R(l,s)$ is a smooth function with no structure at the zeros. All zero-detection comes from the partial Euler product alone.

Clean negative result with methodological implications for empirical zeta function studies.

## Repository Structure

```
paper/
  carry_euler_control.md              The paper
  carry_euler_control_diagram.png     Control experiment figure
experiments/
  H01_carry_zeta_control.py           Full carry product computation
  H02_carry_vs_euler.py               Carry vs Euler product comparison
  H03_corrected_carry_product.py      Euler-corrected carry product (isolated R(l,s) factor)
  H04_convergence_rate.py             Convergence rate study
  H05_shape_correlation.py            Shape correlation analysis
  H06_complex_spectral.py             Complex spectral determinant
```

## Reproduction

```bash
pip install numpy sympy mpmath
python experiments/H01_carry_zeta_control.py
python experiments/H02_carry_vs_euler.py
python experiments/H03_corrected_carry_product.py
python experiments/H04_convergence_rate.py
python experiments/H05_shape_correlation.py
python experiments/H06_complex_spectral.py
```

## Dependencies

- Python >= 3.8, NumPy, SymPy, mpmath

## Companion Papers

| Label | Title | Repository |
|-------|-------|------------|
| [B] | Carry Polynomials and the Euler Product | [`carry-arithmetic-B-zeta-approximation`](https://github.com/stefanoalimonti/carry-arithmetic-B-zeta-approximation) |
| [C] | Eigenvalue Statistics of Carry Companion Matrices: Markov-Driven GOE↔GUE Transition | [`carry-arithmetic-C-matrix-statistics`](https://github.com/stefanoalimonti/carry-arithmetic-C-matrix-statistics) |
| [D] | The Carry-Zero Entropy Bound | [`carry-arithmetic-D-factorization-limits`](https://github.com/stefanoalimonti/carry-arithmetic-D-factorization-limits) |
| [P1] | Pi from Pure Arithmetic | [`carry-arithmetic-P1-pi-spectral`](https://github.com/stefanoalimonti/carry-arithmetic-P1-pi-spectral) |
| [P2] | The Sector Ratio in Binary Multiplication | [`carry-arithmetic-P2-sector-ratio`](https://github.com/stefanoalimonti/carry-arithmetic-P2-sector-ratio) |
| [L] | The Carry–Dirichlet Bridge | [`carry-arithmetic-L-dirichlet-bridge`](https://github.com/stefanoalimonti/carry-arithmetic-L-dirichlet-bridge) |
| [A] | Spectral Theory of Carries in Positional Multiplication | [`carry-arithmetic-A-spectral-theory`](https://github.com/stefanoalimonti/carry-arithmetic-A-spectral-theory) |

### Citation

```bibtex
@article{alimonti2026euler_control,
  author  = {Alimonti, Stefano},
  title   = {Carry Polynomials and the Partial Euler Product: A Control Experiment},
  year    = {2026},
  note    = {Preprint},
  url     = {https://github.com/stefanoalimonti/carry-arithmetic-H-euler-control}
}
```

## License

Paper: CC BY 4.0. Code: MIT License.
