"""
H03_corrected_carry_product.py — Fix the sign: carry det IS the Euler factor

Key insight from prior experiments: r(carry_factor, euler_factor) ≈ -1 per prime.
This means 1/avg_det is anti-correlated with euler_factor = 1/|1-l^{-s}|.
Therefore avg_det IS correlated with euler_factor.

So the correct carry approximation is:
   |ζ(s)| ≈ ∏_l <|det(I - M_l/l^s)|>     (NOT 1/∏ <|det|>)

This script tests BOTH formulations rigorously.

This script evaluates carry products against known ζ(s) values as a baseline
comparison, not as an independent discovery of zeta structure.
"""
import numpy as np
import math
import random

try:
    import mpmath
    mpmath.mp.dps = 30
    HAS_MPMATH = True
except ImportError:
    HAS_MPMATH = False

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

def get_carry_roots(p, q, b=2):
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
    coeffs = [-c for c in carries[1:]]
    if not coeffs or all(c == 0 for c in coeffs):
        return np.array([])
    while coeffs and coeffs[-1] == 0:
        coeffs.pop()
    if not coeffs:
        return np.array([])
    poly_coeffs = [1] + [-c for c in reversed(coeffs)]
    try:
        return np.roots(poly_coeffs)
    except:
        return np.array([])

def spectral_det(roots, l, s):
    if len(roots) == 0:
        return 1.0
    ls = l ** s
    det = 1.0
    for lam in roots:
        det *= abs(1.0 - lam / ls)
    return det

def euler_factor(l, s):
    return 1.0 / abs(1.0 - l**(-s))

def normalize(vals):
    arr = np.array(vals, dtype=float)
    mu = np.mean(arr)
    sigma = np.std(arr)
    if sigma < 1e-15:
        return arr - mu
    return (arr - mu) / sigma

def pearson_r(x, y):
    xn = normalize(x)
    yn = normalize(y)
    return np.dot(xn, yn) / len(xn)

if __name__ == "__main__":
    if not HAS_MPMATH:
        print("ERROR: mpmath required")
        exit(1)
    
    print("=" * 80)
    print(" CORRECTED CARRY PRODUCT: Testing both ∏<det> and 1/∏<det>")
    print("=" * 80)
    
    b = 2
    bit_size = 24
    N_ENSEMBLE = 500
    primes_list = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    
    lo = 1 << (bit_size - 1)
    hi = 1 << bit_size
    half_lo = int(math.isqrt(lo))
    half_hi = int(math.isqrt(hi)) + 100
    
    print(f"\nGenerating {N_ENSEMBLE} semiprimes ({bit_size}-bit)...")
    all_roots = []
    for _ in range(N_ENSEMBLE):
        p = random_prime_range(half_lo, half_hi)
        q = random_prime_range(half_lo, half_hi)
        all_roots.append(get_carry_roots(p, q, b))
    
    t_values = np.linspace(10, 80, 300)
    
    zeta_vals = []
    euler_vals = []
    carry_direct = []    # ∏ <det>
    carry_inverse = []   # ∏ 1/<det>
    carry_geomean = []   # exp(<log det>) for each l, then ∏
    
    print(f"Computing for {len(t_values)} t-values (3 carry formulations)...")
    
    for idx, t in enumerate(t_values):
        s = complex(0.5, t)
        
        z = mpmath.zeta(mpmath.mpc(0.5, t))
        zeta_vals.append(float(abs(z)))
        
        ep = 1.0
        for l in primes_list:
            ep *= euler_factor(l, s)
        euler_vals.append(ep)
        
        cd = 1.0
        ci = 1.0
        cg = 1.0
        
        for l in primes_list:
            dets = []
            for roots in all_roots:
                if len(roots) > 0:
                    d = spectral_det(roots, l, s)
                    if np.isfinite(d) and d > 0:
                        dets.append(d)
            
            if dets:
                avg_det = np.mean(dets)
                geom_det = np.exp(np.mean(np.log(dets)))
                
                cd *= avg_det
                if avg_det > 1e-10:
                    ci *= (1.0 / avg_det)
                cg *= geom_det
        
        carry_direct.append(cd)
        carry_inverse.append(ci)
        carry_geomean.append(cg)
    
    zeta_arr = np.array(zeta_vals)
    euler_arr = np.array(euler_vals)
    cd_arr = np.array(carry_direct)
    ci_arr = np.array(carry_inverse)
    cg_arr = np.array(carry_geomean)
    
    # Also compute 1/carry_direct and 1/carry_geomean
    cd_inv_arr = 1.0 / np.where(cd_arr > 1e-15, cd_arr, 1e-15)
    cg_inv_arr = 1.0 / np.where(cg_arr > 1e-15, cg_arr, 1e-15)
    
    print(f"\n{'='*80}")
    print(f" PEARSON CORRELATIONS WITH |ζ(1/2+it)|")
    print(f"{'='*80}")
    
    formulas = [
        ("Euler prod (baseline)", euler_arr),
        ("∏ <det>  (carry direct)", cd_arr),
        ("∏ 1/<det>  (carry inverse)", ci_arr),
        ("∏ geom_mean(det)", cg_arr),
        ("1/∏ <det>", cd_inv_arr),
        ("1/∏ geom_mean(det)", cg_inv_arr),
    ]
    
    print(f"\n  {'Formula':>35s} | {'r (linear)':>12s} | {'r (log)':>12s}")
    print(f"  {'-'*65}")
    
    log_zeta = np.log(zeta_arr + 1e-15)
    
    for name, vals in formulas:
        r_lin = pearson_r(vals, zeta_arr)
        log_vals = np.log(np.abs(vals) + 1e-15)
        r_log = pearson_r(log_vals, log_zeta)
        marker = " <-- BEST" if r_lin > 0.95 else ""
        print(f"  {name:>35s} | {r_lin:12.6f} | {r_log:12.6f}{marker}")
    
    # ZERO DETECTION for each formula
    print(f"\n{'='*80}")
    print(f" ZERO DETECTION (minima near known Riemann zeros)")
    print(f"{'='*80}")
    
    known_zeros = [14.135, 21.022, 25.011, 30.425, 32.935, 37.586, 40.919,
                   43.327, 48.005, 49.774, 52.970, 56.446, 59.347, 60.832, 65.113,
                   67.080, 69.546, 72.067, 75.705, 77.145]
    
    def find_minima(vals, t_vals, q=0.3):
        minima = []
        for i in range(1, len(vals)-1):
            if vals[i] < vals[i-1] and vals[i] < vals[i+1]:
                if vals[i] < np.quantile(vals, q):
                    minima.append(t_vals[i])
        return minima
    
    def count_hits(minima, known, tol=1.5):
        hits = 0
        total = 0
        for z in known:
            if t_values[0] <= z <= t_values[-1]:
                total += 1
                for m in minima:
                    if abs(m - z) < tol:
                        hits += 1
                        break
        return hits, total
    
    for name, vals in formulas:
        minima = find_minima(vals, t_values)
        h, tot = count_hits(minima, known_zeros)
        print(f"  {name:>35s}: {h}/{tot} zeros detected")
    
    # Direct |ζ| check
    minima_z = find_minima(zeta_vals, t_values)
    hz, tz = count_hits(minima_z, known_zeros)
    print(f"  {'|ζ| (ground truth)':>35s}: {hz}/{tz} zeros detected")
    
    # INFORMATION CONTENT: Does carry product encode info beyond Euler?
    print(f"\n{'='*80}")
    print(f" INFORMATION CONTENT: Residual analysis")
    print(f" After removing Euler prediction, does carry add anything?")
    print(f"{'='*80}")
    
    # Fit: log(carry) = a * log(euler) + b + residual
    log_euler = np.log(euler_arr + 1e-15)
    log_cd = np.log(cd_arr + 1e-15)
    
    # Linear regression in log space
    A = np.column_stack([log_euler, np.ones(len(log_euler))])
    result = np.linalg.lstsq(A, log_cd, rcond=None)
    a, b_coeff = result[0]
    
    log_cd_predicted = a * log_euler + b_coeff
    residual = log_cd - log_cd_predicted
    
    # Does the residual correlate with log(ζ)?
    r_residual_zeta = pearson_r(residual, log_zeta)
    
    log_euler_residual = log_zeta - (np.linalg.lstsq(
        np.column_stack([log_euler, np.ones(len(log_euler))]),
        log_zeta, rcond=None
    )[0][0] * log_euler + np.linalg.lstsq(
        np.column_stack([log_euler, np.ones(len(log_euler))]),
        log_zeta, rcond=None
    )[0][1])
    
    r_carry_residual_vs_zeta_residual = pearson_r(residual, log_euler_residual)
    
    print(f"  log(carry) = {a:.4f} * log(euler) + {b_coeff:.4f}")
    print(f"  r(carry_residual, log ζ): {r_residual_zeta:.6f}")
    print(f"  r(carry_residual, ζ_residual after Euler): {r_carry_residual_vs_zeta_residual:.6f}")
    print()
    print(f"  If r(carry_residual, ζ_residual) >> 0:")
    print(f"    -> Carries capture ζ-information that Euler alone misses!")
    print(f"  If r(carry_residual, ζ_residual) ≈ 0:")
    print(f"    -> Carries add only noise beyond Euler.")
    
    # FINAL VERDICT
    best_r = max((pearson_r(v, zeta_arr), n) for n, v in formulas)
    
    print(f"\n{'='*80}")
    print(f" FINAL VERDICT")
    print(f"{'='*80}")
    print(f"  Best carry formula: {best_r[1]} (r = {best_r[0]:.4f})")
    print(f"  Euler baseline:     r = {pearson_r(euler_arr, zeta_arr):.4f}")
    print(f"  Residual info:      r = {r_carry_residual_vs_zeta_residual:.4f}")
    
    if abs(r_carry_residual_vs_zeta_residual) > 0.1:
        print(f"\n  ** CARRIES CONTAIN INFORMATION BEYOND THE EULER PRODUCT **")
        print(f"  The carry framework captures {abs(r_carry_residual_vs_zeta_residual)*100:.1f}% of")
        print(f"  the ζ variance that the truncated Euler product misses.")
    else:
        print(f"\n  Carries do NOT contain significant information beyond Euler.")
        print(f"  The carry product is approximately a rescaled Euler product.")
