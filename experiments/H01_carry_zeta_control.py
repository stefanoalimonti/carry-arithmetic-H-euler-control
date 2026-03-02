#!/usr/bin/env python3
"""
H01_carry_zeta_control.py — Control test for carry-zeta zeros

H01 found that |Z_carry(1/2+it)| has minima near Riemann zeros.
But Z_carry ≈ |ζ_partial|² · ∏R, and the minima come from the
partial Euler product, not from carry corrections.

This script:
1. Computes the pure partial Euler product |ζ_partial(1/2+it)|²
   (3 lines, no carry, no semiprimes)
2. Computes ∏R(l,s) = Z_carry / |ζ_partial|² to isolate the
   genuine carry contribution
3. Checks whether ∏R has ANY structure near Riemann zeros
"""

import argparse
import sys
import math
import numpy as np

PRIMES_25 = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
             37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
             79, 83, 89, 97]

PAPER_DEFAULTS = {"d_test": 40, "n_semiprimes": 10000, "l_bound": 97}
QUICK_DEFAULTS = {"d_test": 16, "n_semiprimes": 2000,  "l_bound": 47}


def pr(*a, **kw):
    print(*a, **kw)
    sys.stdout.flush()


def multiply_base2(p, q):
    """Return carry polynomial coefficients for p*q in base 2."""
    N = p * q
    p_bits = []
    pp = p
    while pp > 0:
        p_bits.append(pp & 1)
        pp >>= 1
    q_bits = []
    qq = q
    while qq > 0:
        q_bits.append(qq & 1)
        qq >>= 1
    n_bits = []
    nn = N
    while nn > 0:
        n_bits.append(nn & 1)
        nn >>= 1

    dp = len(p_bits)
    dq = len(q_bits)
    M = dp + dq - 2

    convolutions = [0] * (M + 1)
    for i in range(dp):
        for j in range(dq):
            if i + j <= M:
                convolutions[i + j] += p_bits[i] * q_bits[j]

    carries = [0] * (M + 2)
    for k in range(M + 1):
        f_k = n_bits[k] if k < len(n_bits) else 0
        s = convolutions[k] + carries[k]
        carries[k + 1] = (s - f_k) // 2

    return carries[:M + 2]


def companion_matrix(carries):
    """Build companion matrix from carry polynomial."""
    top = len(carries) - 1
    while top > 0 and carries[top] == 0:
        top -= 1
    if top <= 1:
        return None
    coeffs = carries[1:top + 1]
    D = len(coeffs)
    if D < 3:
        return None
    M = np.zeros((D, D))
    M[0, :] = [-c / coeffs[-1] for c in reversed(coeffs)]
    for i in range(1, D):
        M[i, i - 1] = 1.0
    return M


def parse_args():
    p = argparse.ArgumentParser(
        description="Carry-zeta control test (H01).",
        epilog="Paper parameters (Appendix A): --paper  or  -D 40 -N 10000 -L 97")
    p.add_argument("--paper", action="store_true",
                   help="Use full paper parameters (D=40, N=10000, L=97). Slow.")
    p.add_argument("-D", type=int, default=None, dest="d_test",
                   help=f"Semiprime bit length (default: {QUICK_DEFAULTS['d_test']})")
    p.add_argument("-N", type=int, default=None, dest="n_semiprimes",
                   help=f"Ensemble size (default: {QUICK_DEFAULTS['n_semiprimes']})")
    p.add_argument("-L", type=int, default=None, dest="l_bound",
                   help=f"Prime bound (default: {QUICK_DEFAULTS['l_bound']})")
    args = p.parse_args()
    defaults = PAPER_DEFAULTS if args.paper else QUICK_DEFAULTS
    d_test = args.d_test or defaults["d_test"]
    n_semiprimes = args.n_semiprimes or defaults["n_semiprimes"]
    l_bound = args.l_bound or defaults["l_bound"]
    return d_test, n_semiprimes, l_bound


def main():
    d_test, n_semiprimes, l_bound = parse_args()
    test_primes = [p for p in PRIMES_25 if p <= l_bound]

    pr("=" * 78)
    pr("  CARRY-ZETA CONTROL TEST")
    pr("  Is the structure from carry corrections or just the Euler product?")
    pr("=" * 78)
    pr(f"  Parameters: D={d_test}, N={n_semiprimes}, L={l_bound} ({len(test_primes)} primes)")
    if d_test < PAPER_DEFAULTS["d_test"]:
        pr(f"  NOTE: reduced parameters for quick run. "
           f"Paper uses D={PAPER_DEFAULTS['d_test']}, N={PAPER_DEFAULTS['n_semiprimes']}, "
           f"L={PAPER_DEFAULTS['l_bound']}.")
        pr(f"  Run with --paper to reproduce Appendix A results (slow).")
    pr()

    t_values = np.linspace(1, 50, 500)
    known_zeros = [14.1347, 21.022, 25.011, 30.425, 32.935, 37.586, 40.919, 43.327, 48.005]

    # ═══════════════════════════════════════════════════════════════
    pr("═" * 78)
    pr("  PART 1: PURE PARTIAL EULER PRODUCT (no carry, no semiprimes)")
    pr("═" * 78)
    pr()

    euler_product_sq = np.ones(len(t_values))
    for l in test_primes:
        for it_idx, t in enumerate(t_values):
            s = 0.5 + 1j * t
            euler_factor = abs(1 - l**(-s))
            euler_product_sq[it_idx] *= 1.0 / euler_factor**2

    pr("  |ζ_partial(1/2+it)|² computed for l ∈ {3,...,47}")
    pr(f"  Range: [{np.min(euler_product_sq):.4f}, {np.max(euler_product_sq):.4f}]")

    euler_mins = []
    for i in range(1, len(euler_product_sq) - 1):
        if euler_product_sq[i] < euler_product_sq[i-1] and euler_product_sq[i] < euler_product_sq[i+1]:
            euler_mins.append((t_values[i], euler_product_sq[i]))

    pr(f"\n  Local minima of |ζ_partial(1/2+it)|²:")
    for t, z in euler_mins[:15]:
        pr(f"    t = {t:8.4f}  |ζ_partial|² = {z:.6f}")

    pr(f"\n  Distance to known Riemann zeros:")
    for z in known_zeros:
        dists = [abs(z - t) for t, _ in euler_mins]
        nearest = min(dists) if dists else float('inf')
        pr(f"    ζ zero at t={z:.3f}: nearest Euler min at Δt = {nearest:.3f}")

    # ═══════════════════════════════════════════════════════════════
    pr()
    pr("═" * 78)
    pr("  PART 2: CARRY PRODUCT ⟨|det|⟩ · |1 - l^{-s}| = R(l,s)")
    pr("  (isolating the genuine carry correction)")
    pr("═" * 78)
    pr()

    np.random.seed(42)

    R_product = np.ones(len(t_values))

    for l in test_primes:
        det_avg = np.zeros(len(t_values))
        n_valid = 0

        for _ in range(n_semiprimes):
            p = (1 << (d_test - 1)) | np.random.randint(0, 1 << (d_test - 1)) | 1
            q = (1 << (d_test - 1)) | np.random.randint(0, 1 << (d_test - 1)) | 1
            carries = multiply_base2(int(p), int(q))
            M_c = companion_matrix(carries)
            if M_c is None:
                continue

            D = M_c.shape[0]
            I_D = np.eye(D)

            for it_idx, t in enumerate(t_values):
                s = 0.5 + 1j * t
                try:
                    det_val = np.linalg.det(I_D - M_c / l**s)
                    det_avg[it_idx] += abs(det_val)
                except:
                    pass

            n_valid += 1

        if n_valid > 0:
            det_avg /= n_valid

            R_l = np.zeros(len(t_values))
            for it_idx, t in enumerate(t_values):
                s = 0.5 + 1j * t
                euler_factor = abs(1 - l**(-s))
                R_l[it_idx] = det_avg[it_idx] * euler_factor

            R_product *= R_l
            pr(f"  l={l:2d}: ⟨R(l,s)⟩ mean={np.mean(R_l):.6f}, "
               f"std={np.std(R_l):.6f}, range=[{np.min(R_l):.4f}, {np.max(R_l):.4f}]")

    pr(f"\n  ∏R(l,s) product:")
    pr(f"  Range: [{np.min(R_product):.6f}, {np.max(R_product):.6f}]")
    pr(f"  Mean:  {np.mean(R_product):.6f}")
    pr(f"  Std:   {np.std(R_product):.6f}")
    pr(f"  CV:    {np.std(R_product)/np.mean(R_product):.4f}")

    R_mins = []
    for i in range(1, len(R_product) - 1):
        if R_product[i] < R_product[i-1] and R_product[i] < R_product[i+1]:
            R_mins.append((t_values[i], R_product[i]))

    pr(f"\n  Local minima of ∏R(l,1/2+it):")
    for t, z in R_mins[:15]:
        pr(f"    t = {t:8.4f}  ∏R = {z:.6f}")

    pr(f"\n  Distance from known Riemann zeros to nearest ∏R minimum:")
    for z in known_zeros:
        dists = [abs(z - t) for t, _ in R_mins]
        nearest = min(dists) if dists else float('inf')
        pr(f"    ζ zero at t={z:.3f}: nearest ∏R min at Δt = {nearest:.3f}")

    # ═══════════════════════════════════════════════════════════════
    pr()
    pr("═" * 78)
    pr("  PART 3: COMPARISON — EULER vs CARRY CORRECTIONS")
    pr("═" * 78)
    pr()

    pr("  Riemann zero distances:")
    pr(f"  {'ζ zero':>12s}  {'Δt(Euler)':>10s}  {'Δt(∏R)':>10s}  {'∏R smooth?':>12s}")
    pr(f"  {'-'*12}  {'-'*10}  {'-'*10}  {'-'*12}")
    for z in known_zeros:
        d_euler = min([abs(z - t) for t, _ in euler_mins]) if euler_mins else float('inf')
        d_R = min([abs(z - t) for t, _ in R_mins]) if R_mins else float('inf')
        smooth = "YES" if d_R > 1.0 else "no"
        pr(f"  {z:12.3f}  {d_euler:10.3f}  {d_R:10.3f}  {smooth:>12s}")

    # ═══════════════════════════════════════════════════════════════
    pr()
    pr("═" * 78)
    pr("  VERDICT")
    pr("═" * 78)
    pr()
    pr("  If Δt(Euler) ≈ Δt(H01 carry-zeta): the minima come from the")
    pr("  partial Euler product, NOT from carry corrections.")
    pr()
    pr("  If ∏R is smooth (no sharp minima near ζ zeros): the carry")
    pr("  correction adds no information about Riemann zeros beyond")
    pr("  what the Euler product already contains.")
    pr("=" * 78)


if __name__ == '__main__':
    main()
