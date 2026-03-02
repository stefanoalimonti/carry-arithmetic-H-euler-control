"""
H05_shape_correlation.py — Rigorous comparison: Carry vs Euler vs ζ(s)

Instead of asking "which product is closer to |ζ(s)|", we ask:
1. Does the carry product have a better SHAPE correlation with |ζ(s)|?
   (Remove multiplicative constants — compare normalized curves)
2. Is the carry factor for each prime just a rescaled Euler factor?
   (If carry_factor(l,s) ≈ C·euler_factor(l,s), carries add nothing)
3. Does the carry product correctly predict WHERE the zeros are?
   (Does it have minima near the Riemann zeros?)

This uses mpmath for exact ζ(s) computation.
"""
import numpy as np
import math
import random
from collections import defaultdict

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
    ls = l ** s
    return 1.0 / abs(1.0 - 1.0/ls)

def normalize(vals):
    """Normalize to zero mean, unit variance."""
    arr = np.array(vals)
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
        print("ERROR: mpmath required. Install with: pip install mpmath")
        exit(1)
    
    print("=" * 80)
    print(" SHAPE CORRELATION: Carry vs Euler vs |ζ(s)|")
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
    
    # Scan t from 10 to 80 in fine steps
    t_values = np.linspace(10, 80, 200)
    
    zeta_vals = []
    euler_vals = []
    carry_vals = []
    
    print(f"Computing |ζ(1/2+it)|, Euler product, Carry product for {len(t_values)} t-values...")
    
    for t in t_values:
        s = complex(0.5, t)
        
        # Exact ζ
        z = mpmath.zeta(mpmath.mpc(0.5, t))
        zeta_vals.append(float(abs(z)))
        
        # Euler product
        ep = 1.0
        for l in primes_list:
            ep *= euler_factor(l, s)
        euler_vals.append(ep)
        
        # Carry product
        cp = 1.0
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
                if avg_det > 1e-10:
                    cp *= (1.0 / avg_det)
        carry_vals.append(cp)
    
    zeta_arr = np.array(zeta_vals)
    euler_arr = np.array(euler_vals)
    carry_arr = np.array(carry_vals)
    
    # CORRELATION ANALYSIS
    print(f"\n{'='*80}")
    print(f" TEST 1: SHAPE CORRELATION (Pearson r)")
    print(f"{'='*80}")
    
    r_euler_zeta = pearson_r(euler_arr, zeta_arr)
    r_carry_zeta = pearson_r(carry_arr, zeta_arr)
    r_euler_carry = pearson_r(euler_arr, carry_arr)
    
    print(f"  Pearson r(Euler, ζ):  {r_euler_zeta:.6f}")
    print(f"  Pearson r(Carry, ζ):  {r_carry_zeta:.6f}")
    print(f"  Pearson r(Euler, Carry): {r_euler_carry:.6f}")
    
    if r_carry_zeta > r_euler_zeta + 0.01:
        print(f"\n  -> CARRY has BETTER shape correlation with ζ!")
    elif r_euler_zeta > r_carry_zeta + 0.01:
        print(f"\n  -> EULER has better shape correlation with ζ.")
    else:
        print(f"\n  -> Both have similar shape correlation with ζ.")
    
    # LOG-CORRELATION (for multiplicative structures)
    print(f"\n{'='*80}")
    print(f" TEST 2: LOG-SPACE CORRELATION")
    print(f"{'='*80}")
    
    log_zeta = np.log(zeta_arr + 1e-15)
    log_euler = np.log(euler_arr + 1e-15)
    log_carry = np.log(carry_arr + 1e-15)
    
    rl_euler_zeta = pearson_r(log_euler, log_zeta)
    rl_carry_zeta = pearson_r(log_carry, log_zeta)
    
    print(f"  Pearson r(log Euler, log ζ):  {rl_euler_zeta:.6f}")
    print(f"  Pearson r(log Carry, log ζ):  {rl_carry_zeta:.6f}")
    
    # ZERO DETECTION
    print(f"\n{'='*80}")
    print(f" TEST 3: ZERO DETECTION")
    print(f" Can the carry product locate Riemann zeros?")
    print(f"{'='*80}")
    
    known_zeros = [14.135, 21.022, 25.011, 30.425, 32.935, 37.586, 40.919, 
                   43.327, 48.005, 49.774, 52.970, 56.446, 59.347, 60.832, 65.113,
                   67.080, 69.546, 72.067, 75.705, 77.145]
    
    # Find local minima in each product
    def find_local_minima(vals, t_vals, threshold_quantile=0.3):
        """Find positions where the function is a local minimum."""
        minima = []
        for i in range(1, len(vals)-1):
            if vals[i] < vals[i-1] and vals[i] < vals[i+1]:
                if vals[i] < np.quantile(vals, threshold_quantile):
                    minima.append(t_vals[i])
        return minima
    
    zeta_minima = find_local_minima(zeta_vals, t_values, 0.25)
    euler_minima = find_local_minima(euler_vals, t_values, 0.25)
    carry_minima = find_local_minima(carry_vals, t_values, 0.25)
    
    def count_hits(minima, known, tol=1.5):
        hits = 0
        for z in known:
            if z < t_values[0] or z > t_values[-1]:
                continue
            for m in minima:
                if abs(m - z) < tol:
                    hits += 1
                    break
        total_in_range = sum(1 for z in known if t_values[0] <= z <= t_values[-1])
        return hits, total_in_range
    
    zh, zt = count_hits(zeta_minima, known_zeros)
    eh, et = count_hits(euler_minima, known_zeros)
    ch, ct = count_hits(carry_minima, known_zeros)
    
    print(f"  |ζ| minima near known zeros: {zh}/{zt} (sanity check)")
    print(f"  Euler product minima near known zeros: {eh}/{et}")
    print(f"  Carry product minima near known zeros: {ch}/{ct}")
    print(f"  (tolerance: 1.5)")
    
    # TEST 4: PER-FACTOR ANALYSIS
    print(f"\n{'='*80}")
    print(f" TEST 4: PER-PRIME FACTOR COMPARISON")
    print(f" Is carry_factor(l,s) ≈ C · euler_factor(l,s)?")
    print(f"{'='*80}")
    
    print(f"\n  {'l':>4s} | {'r(carry_f, euler_f)':>20s} | {'mean ratio C':>12s} | {'std ratio':>10s}")
    print(f"  {'-'*60}")
    
    for l in primes_list[:10]:
        euler_fs = []
        carry_fs = []
        for t in t_values[::5]:
            s = complex(0.5, t)
            ef = euler_factor(l, s)
            
            det_sum = 0.0
            valid = 0
            for roots in all_roots:
                if len(roots) > 0:
                    d = spectral_det(roots, l, s)
                    if np.isfinite(d) and d > 0:
                        det_sum += d
                        valid += 1
            if valid > 0:
                cf = 1.0 / (det_sum / valid)
            else:
                cf = 1.0
            
            euler_fs.append(ef)
            carry_fs.append(cf)
        
        r_val = pearson_r(euler_fs, carry_fs)
        ratios = [c/e if e > 1e-15 else 0 for c, e in zip(carry_fs, euler_fs)]
        mean_C = np.mean(ratios)
        std_C = np.std(ratios)
        
        print(f"  {l:4d} | {r_val:20.6f} | {mean_C:12.4f} | {std_C:10.4f}")
    
    print(f"\n  If r ≈ 1 and std(C) ≈ 0: carry factor is JUST a rescaled Euler factor.")
    print(f"  If r << 1 or std(C) >> 0: carry factor has DIFFERENT structure.")
    
    # FINAL VERDICT
    print(f"\n{'='*80}")
    print(f" FINAL VERDICT")
    print(f"{'='*80}")
    
    extra_info = r_carry_zeta > r_euler_zeta + 0.02
    different_shape = r_euler_carry < 0.9
    detects_zeros = ch > eh
    
    score = 0
    if extra_info:
        score += 1
        print(f"  [+] Carry has better shape correlation with ζ ({r_carry_zeta:.4f} vs {r_euler_zeta:.4f})")
    else:
        print(f"  [-] Carry does NOT have better shape correlation ({r_carry_zeta:.4f} vs {r_euler_zeta:.4f})")
    
    if different_shape:
        score += 1
        print(f"  [+] Carry product has different structure from Euler ({r_euler_carry:.4f})")
    else:
        print(f"  [-] Carry product tracks Euler closely ({r_euler_carry:.4f})")
    
    if detects_zeros:
        score += 1
        print(f"  [+] Carry detects MORE Riemann zeros ({ch}/{ct} vs {eh}/{et})")
    else:
        print(f"  [-] Carry does NOT detect more zeros ({ch}/{ct} vs {eh}/{et})")
    
    print(f"\n  Score: {score}/3")
    if score >= 2:
        print(f"  CONCLUSION: Evidence that carries contain information BEYOND Euler.")
        print(f"  The carry framework may capture arithmetic structure that the")
        print(f"  truncated Euler product misses.")
    elif score == 1:
        print(f"  CONCLUSION: Weak/mixed evidence. Needs further investigation.")
    else:
        print(f"  CONCLUSION: Carries appear to be a rescaled Euler product.")
        print(f"  No evidence of extra information.")
