"""
H02_carry_vs_euler.py — The Million Dollar Experiment

Question: Is the carry product BETTER, EQUAL, or WORSE than the raw
truncated Euler product at approximating |ζ(s)| on the critical line?

If BETTER: carries contain extra arithmetic information beyond Euler.
           The carry framework might be more fundamental than ζ(s).
If EQUAL:  carries are just a complicated way to compute the Euler product.
           The ~5% residual R is pure noise.
If WORSE:  the ensemble averaging adds noise that degrades the approximation.

Method:
1. Compute |ζ(1/2 + it)| exactly via mpmath
2. Compute truncated Euler product: E_L(t) = prod_{l<=L} |1 - l^{-1/2-it}|^{-1}
3. Compute carry product: C_L(t) = prod_{l<=L} <|det(I - M_l / l^{1/2+it})|>
   using ensemble of random semiprimes
4. Compare |E_L(t) - ζ| vs |C_L(t) - ζ| at Riemann zero locations
"""
import numpy as np
import math
import random
from collections import defaultdict

random.seed(42)

def to_digits(n, base):
    if n == 0: return [0]
    digits = []
    while n > 0:
        digits.append(n % base)
        n //= base
    return digits

def is_prime_miller_rabin(n, k=15):
    if n < 2: return False
    if n == 2 or n == 3: return True
    if n % 2 == 0: return False
    r, d = 0, n - 1
    while d % 2 == 0:
        r += 1
        d //= 2
    for _ in range(k):
        a = random.randrange(2, n - 1)
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def random_prime_range(lo, hi):
    while True:
        n = random.randrange(lo | 1, hi, 2)
        if is_prime_miller_rabin(n):
            return n

def get_carry_polynomial_roots(p, q, b=2):
    """Get eigenvalues of companion matrix from carry polynomial."""
    g = to_digits(p, b)
    h = to_digits(q, b)
    f = to_digits(p * q, b)
    
    max_len = len(g) + len(h)
    g += [0] * (max_len - len(g))
    h += [0] * (max_len - len(h))
    f += [0] * (max_len - len(f))
    
    carries = [0]
    for k in range(max_len):
        conv_k = sum(g[i] * h[k-i] for i in range(k+1) if 0 <= k-i < len(h))
        c_next = (conv_k + carries[-1] - f[k]) // b
        carries.append(c_next)
    
    while len(carries) > 1 and carries[-1] == 0:
        carries.pop()
    
    # Quotient polynomial coefficients: q_k = -carry_{k+1}
    coeffs = [-c for c in carries[1:]]
    
    if not coeffs or all(c == 0 for c in coeffs):
        return np.array([])
    
    # Remove trailing zeros
    while coeffs and coeffs[-1] == 0:
        coeffs.pop()
    
    if not coeffs:
        return np.array([])
    
    # Companion matrix eigenvalues = roots of Q(x)
    # numpy wants [leading, ..., constant] for np.roots
    poly_coeffs = [1] + [-c for c in reversed(coeffs)]
    
    try:
        roots = np.roots(poly_coeffs)
        return roots
    except:
        return np.array([])

def spectral_det(roots, l, s):
    """Compute |det(I - M/l^s)| = prod|1 - lambda_j/l^s|"""
    if len(roots) == 0:
        return 1.0
    
    ls = l ** s
    det = 1.0
    for lam in roots:
        det *= abs(1.0 - lam / ls)
    
    return det

def euler_factor(l, s):
    """Compute |1 - l^{-s}|^{-1}"""
    ls = l ** s
    return 1.0 / abs(1.0 - 1.0/ls)

def zeta_approx_hardy(t, num_terms=5000):
    """
    Approximate |ζ(1/2 + it)| using partial Dirichlet series.
    This is a rough approximation for comparison purposes.
    For t near a zero, |ζ| should be very small.
    """
    s = complex(0.5, t)
    total = complex(0, 0)
    for n in range(1, num_terms + 1):
        total += n ** (-s)
    return abs(total)

# Known Riemann zeros (imaginary parts)
RIEMANN_ZEROS = [
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
    37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
    52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
    67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
]

if __name__ == "__main__":
    print("=" * 80)
    print(" THE MILLION DOLLAR EXPERIMENT: Carry Product vs Euler Product")
    print("=" * 80)
    
    b = 2
    bit_size = 24
    N_ENSEMBLE = 300
    
    # Primes for the Euler product
    primes_list = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    L_max = max(primes_list)
    
    # Generate ensemble of semiprimes
    lo = 1 << (bit_size - 1)
    hi = 1 << bit_size
    half_lo = int(math.isqrt(lo))
    half_hi = int(math.isqrt(hi)) + 100
    
    print(f"\nGenerating {N_ENSEMBLE} semiprimes ({bit_size}-bit)...")
    semiprimes = []
    for _ in range(N_ENSEMBLE):
        p = random_prime_range(half_lo, half_hi)
        q = random_prime_range(half_lo, half_hi)
        semiprimes.append((p, q))
    
    # For each prime l, compute the ensemble-averaged carry determinant
    print(f"Computing carry determinants for {len(primes_list)} primes...")
    
    # Pre-compute all roots for the ensemble
    all_roots = []
    for p, q in semiprimes:
        roots = get_carry_polynomial_roots(p, q, b)
        all_roots.append(roots)
    
    # Now compare at each Riemann zero
    print(f"\n{'t':>10s} | {'|ζ| (Dirichlet)':>15s} | {'Euler prod':>12s} | {'Carry prod':>12s} | {'R = C/E':>10s} | {'Better?':>10s}")
    print("-" * 85)
    
    euler_errors = []
    carry_errors = []
    R_values = []
    
    for t in RIEMANN_ZEROS[:15]:
        s = complex(0.5, t)
        
        # 1. "True" |ζ(s)| from Dirichlet series (rough but unbiased)
        zeta_val = zeta_approx_hardy(t, num_terms=10000)
        
        # 2. Truncated Euler product
        euler_prod = 1.0
        for l in primes_list:
            euler_prod *= euler_factor(l, s)
        
        # 3. Carry product (ensemble averaged per-factor)
        carry_prod = 1.0
        for l in primes_list:
            det_sum = 0.0
            valid = 0
            for roots in all_roots:
                if len(roots) > 0:
                    d = spectral_det(roots, l, s)
                    if np.isfinite(d) and d > 0:
                        det_sum += d
                        valid += 1
            
            if valid > 0:
                avg_det = det_sum / valid
                # carry_factor = 1/avg_det (analogous to Euler factor)
                if avg_det > 1e-10:
                    carry_prod *= (1.0 / avg_det)
        
        # Ratio R = Carry/Euler
        if euler_prod > 0:
            R = carry_prod / euler_prod
        else:
            R = float('inf')
        
        R_values.append(R)
        
        # Which is closer to zeta?
        euler_err = abs(euler_prod - zeta_val)
        carry_err = abs(carry_prod - zeta_val)
        euler_errors.append(euler_err)
        carry_errors.append(carry_err)
        
        better = "CARRY" if carry_err < euler_err else "EULER" if euler_err < carry_err else "TIE"
        
        print(f"{t:10.4f} | {zeta_val:15.6f} | {euler_prod:12.6f} | {carry_prod:12.6f} | {R:10.4f} | {better:>10s}")
    
    # Summary
    print(f"\n{'='*80}")
    print(f" SUMMARY")
    print(f"{'='*80}")
    
    euler_wins = sum(1 for e, c in zip(euler_errors, carry_errors) if e < c)
    carry_wins = sum(1 for e, c in zip(euler_errors, carry_errors) if c < e)
    ties = sum(1 for e, c in zip(euler_errors, carry_errors) if abs(e-c) < 1e-10)
    
    print(f"  Euler wins: {euler_wins}/{len(euler_errors)}")
    print(f"  Carry wins: {carry_wins}/{len(euler_errors)}")
    print(f"  Ties:       {ties}/{len(euler_errors)}")
    
    mean_euler_err = sum(euler_errors) / len(euler_errors)
    mean_carry_err = sum(carry_errors) / len(carry_errors)
    
    print(f"\n  Mean Euler error:  {mean_euler_err:.6f}")
    print(f"  Mean Carry error:  {mean_carry_err:.6f}")
    print(f"  Carry/Euler ratio: {mean_carry_err/mean_euler_err:.4f}")
    
    mean_R = sum(R_values) / len(R_values)
    std_R = (sum((r - mean_R)**2 for r in R_values) / len(R_values)) ** 0.5
    
    print(f"\n  Mean R = Carry/Euler: {mean_R:.6f} ± {std_R:.6f}")
    print(f"  R = 1 would mean carries add nothing beyond Euler.")
    print(f"  R < 1 would mean carries are a WORSE approximation.")
    print(f"  R > 1 would mean carries contain EXTRA information.")
    
    if mean_R > 1.05:
        verdict = "CARRIES MAY CONTAIN EXTRA INFORMATION"
    elif mean_R < 0.95:
        verdict = "CARRIES ADD NOISE (WORSE THAN EULER)"
    else:
        verdict = "CARRIES ≈ EULER (NO EXTRA INFORMATION)"
    
    print(f"\n  VERDICT: {verdict}")
    
    # Detailed R analysis
    print(f"\n{'='*80}")
    print(f" R(t) PROFILE: Does R depend on t?")
    print(f"{'='*80}")
    for t, R in zip(RIEMANN_ZEROS[:15], R_values):
        bar = "#" * int(min(R, 3) * 20)
        print(f"  t={t:8.3f}: R={R:8.4f} {bar}")
