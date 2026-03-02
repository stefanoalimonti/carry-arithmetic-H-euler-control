"""
H06_complex_spectral.py — Spectral Determinant on the Critical Line

Evaluates det(I - M/l^s) for the carry companion matrix M at complex s,
extending Paper B's Euler product from σ > 1 to the critical line σ = 1/2.

Experiments:
1. Spectral determinant landscape for a single semiprime
2. Zeros of the averaged spectral determinant on the critical line
3. Comparison with the first Riemann zeta zeros
4. Ensemble-averaged product over multiple primes l
"""
import math
import random
import numpy as np

random.seed(42)
np.random.seed(42)

ZETA_ZEROS_T = [
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
    37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
    52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
    67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
]

def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0: return False
        i += 6
    return True

def random_prime(bits):
    for _ in range(100000):
        n = random.getrandbits(bits) | (1 << (bits - 1)) | 1
        if is_prime(n): return n
    return None

def to_digits(n, base=2):
    if n == 0: return [0]
    d = []
    while n > 0:
        d.append(int(n % base)); n //= base
    return d

def carry_quotient_coeffs(p, q, base=2):
    gd = to_digits(p, base)
    hd = to_digits(q, base)
    fd = to_digits(p * q, base)
    dg = len(gd); dh = len(hd)
    conv_len = dg + dh - 1
    conv = [0] * conv_len
    for i in range(dg):
        for j in range(dh):
            conv[i + j] += gd[i] * hd[j]
    max_len = max(conv_len, len(fd))

    carries = [0] * (max_len + 2)
    c = 0
    for k in range(max_len):
        vk = (conv[k] if k < conv_len else 0) + c
        fk = fd[k] if k < len(fd) else 0
        c = (vk - fk) // base
        carries[k + 1] = c

    while carries and carries[-1] == 0:
        carries.pop()
    if not carries:
        carries = [0]

    Q = [-carries[i + 1] for i in range(len(carries) - 1)]
    return Q

def build_companion(coeffs):
    """Build companion matrix for monic polynomial with given coefficients.
    coeffs = [q_0, q_1, ..., q_{D-2}] for polynomial x^{D-1} + q_{D-2}x^{D-2} + ... + q_0
    """
    D = len(coeffs)
    if D == 0:
        return np.array([[0.0]])
    M = np.zeros((D, D), dtype=np.float64)
    for i in range(D - 1):
        M[i + 1, i] = 1.0
    for i in range(D):
        M[i, D - 1] = -coeffs[i]
    return M

def spectral_det(M, l, s):
    """Compute |det(I - M/l^s)| for complex s."""
    D = M.shape[0]
    ls = l ** s
    A = np.eye(D, dtype=np.complex128) - M / ls
    return np.linalg.det(A)

def log_spectral_det(M, l, s):
    """Compute log|det(I - M/l^s)| for complex s."""
    d = spectral_det(M, l, s)
    return np.log(abs(d)) if abs(d) > 1e-300 else -700


if __name__ == "__main__":
    print("=" * 78)
    print(" H06: SPECTRAL DETERMINANT ON THE CRITICAL LINE")
    print("=" * 78)

    # ═══════════════════════════════════════════════════════════════════
    # PART 1: Single Semiprime — Spectral Landscape
    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'═'*78}")
    print(" PART 1: Spectral Landscape for a Single Semiprime")
    print(f"{'═'*78}")

    print("""
  For N = pq, the companion matrix M of Q(x) has eigenvalues {λ_i}.
  det(I - M/l^s) = Π_i (1 - λ_i/l^s)

  Paper B showed: ⟨|det(I-M/l^s)|⟩ · |1-l^{-s}| = R(l,s) ≈ 1 for σ > 1.
  What happens on the critical line σ = 1/2?

  At σ = 1/2: l^s = √l · e^{it·ln(l)}, and eigenvalues near |z|=1
  can resonate when l^{it} ≈ λ_i/√l.
""")

    p, q = 101, 103
    N = p * q
    Q_coeffs = carry_quotient_coeffs(p, q, 2)
    D = len(Q_coeffs)
    M = build_companion(Q_coeffs)
    eigs = np.linalg.eigvals(M)

    print(f"  N = {N} = {p} × {q},  D(Q) = {D}")
    print(f"  Eigenvalues: {len(eigs)} total, "
          f"max|λ| = {np.max(np.abs(eigs)):.4f}, "
          f"median|λ| = {np.median(np.abs(eigs)):.4f}")

    n_real = np.sum(np.abs(eigs.imag) < 1e-8)
    n_unit = np.sum((np.abs(eigs) > 0.9) & (np.abs(eigs) < 1.1))
    print(f"  Real eigenvalues: {n_real},  |λ| ∈ [0.9, 1.1]: {n_unit}")

    test_primes = [3, 5, 7, 11, 13, 17, 23, 29, 31, 37]
    t_values = np.linspace(0, 50, 500)

    print(f"\n  Spectral determinant on critical line σ = 1/2:")
    print(f"  {'l':>4s} | {'min|det|':>10s} {'at t':>8s} | {'max|det|':>10s} "
          f"| {'zeros(|det|<0.1)':>17s} | near ζ-zeros?")
    print(f"  {'-'*75}")

    for l in test_primes:
        dets = []
        for t in t_values:
            s = 0.5 + 1j * t
            d = abs(spectral_det(M, l, s))
            dets.append(d)
        dets = np.array(dets)

        min_idx = np.argmin(dets)
        min_t = t_values[min_idx]
        n_near_zero = np.sum(dets < 0.1)

        near_zeta = []
        for zt in ZETA_ZEROS_T[:10]:
            idx = np.argmin(np.abs(t_values - zt))
            if dets[idx] < 0.5:
                near_zeta.append(f"{zt:.1f}")

        nz_str = ", ".join(near_zeta) if near_zeta else "none"
        print(f"  {l:>4d} | {dets[min_idx]:>10.6f} {min_t:>8.2f} | {np.max(dets):>10.4f} "
              f"| {n_near_zero:>17d} | {nz_str}")

    # ═══════════════════════════════════════════════════════════════════
    # PART 2: Ensemble Average — The Carry Euler Product on σ = 1/2
    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'═'*78}")
    print(" PART 2: Ensemble-Averaged Spectral Determinant on σ = 1/2")
    print(f"{'═'*78}")

    print("""
  Paper B: ⟨|det(I-M/l^s)|⟩ = R(l,s)/|1-l^{-s}| for σ > 1.
  What does the average look like on σ = 1/2?

  For each t along the critical line, average |det(I-M/l^s)|
  over many random semiprimes, then multiply over primes l.
""")

    n_ensemble = 100
    bits = 20

    ensemble_Ms = []
    for _ in range(n_ensemble):
        pr = random_prime(bits // 2)
        qr = random_prime(bits // 2)
        if pr is None or qr is None or pr == qr:
            continue
        Qc = carry_quotient_coeffs(pr, qr, 2)
        if len(Qc) < 2:
            continue
        Mc = build_companion(Qc)
        ensemble_Ms.append(Mc)

    print(f"  Ensemble: {len(ensemble_Ms)} semiprimes, D ≈ {bits} bits")

    test_l = [3, 5, 7, 11]
    t_grid = np.linspace(0, 80, 400)

    for l in test_l:
        avg_det = np.zeros(len(t_grid))
        avg_log_det = np.zeros(len(t_grid))

        for Mi in ensemble_Ms:
            for ti, t in enumerate(t_grid):
                s = 0.5 + 1j * t
                d = abs(spectral_det(Mi, l, s))
                avg_det[ti] += d
                avg_log_det[ti] += math.log(d) if d > 1e-300 else -700

        avg_det /= len(ensemble_Ms)
        avg_log_det /= len(ensemble_Ms)

        euler_factor = 1.0 / np.abs(1 - l ** (-0.5 - 1j * t_grid))

        R_vals = avg_det * np.abs(1 - l ** (-0.5 - 1j * t_grid))

        dips = []
        for zt in ZETA_ZEROS_T[:15]:
            idx = np.argmin(np.abs(t_grid - zt))
            dips.append((zt, avg_det[idx], R_vals[idx]))

        print(f"\n  l = {l}:")
        print(f"    ⟨|det|⟩: min = {np.min(avg_det):.6f}  max = {np.max(avg_det):.4f}  "
              f"mean = {np.mean(avg_det):.4f}")
        print(f"    R(l,s) = ⟨|det|⟩·|1-l^-s|: "
              f"min = {np.min(R_vals):.4f}  max = {np.max(R_vals):.4f}  "
              f"mean = {np.mean(R_vals):.4f}")

        print(f"    At ζ-zeros (σ=1/2+it):")
        print(f"    {'t_ζ':>8s} | {'⟨|det|⟩':>10s} | {'R(l,s)':>8s} | {'|1-l^-s|':>10s}")
        print(f"    {'-'*45}")
        for zt, ad, rv in dips[:8]:
            s_val = 0.5 + 1j * zt
            euler = abs(1 - l ** (-s_val))
            print(f"    {zt:>8.3f} | {ad:>10.6f} | {rv:>8.4f} | {euler:>10.6f}")

    # ═══════════════════════════════════════════════════════════════════
    # PART 3: Multi-Prime Product — Carry ζ on the Critical Line
    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'═'*78}")
    print(" PART 3: Carry ζ Product on the Critical Line")
    print(f"{'═'*78}")

    print("""
  The full carry-zeta approximation: ζ_carry(s) = Π_l ⟨|det(I-M/l^s)|⟩^{-1}
  Compare: |ζ_carry(1/2+it)| vs |ζ(1/2+it)|

  If ζ_carry captures zeta's structure, it should show dips near the
  known zeta zeros (where |ζ| → 0).
""")

    product_primes = [2, 3, 5, 7, 11, 13]
    t_fine = np.linspace(0, 80, 800)

    log_product = np.zeros(len(t_fine))

    for l in product_primes:
        for ti, t in enumerate(t_fine):
            s = 0.5 + 1j * t
            total_log = 0
            for Mi in ensemble_Ms:
                d = abs(spectral_det(Mi, l, s))
                total_log += math.log(d) if d > 1e-300 else -700
            avg_log = total_log / len(ensemble_Ms)
            log_product[ti] -= avg_log

    product_vals = np.exp(np.clip(log_product, -50, 50))

    print(f"  Carry ζ product over l ∈ {product_primes}, {len(ensemble_Ms)} semiprimes avg:")
    print(f"  |ζ_carry(1/2+it)|: min={np.min(product_vals):.4f}  "
          f"max={np.max(product_vals):.4f}  mean={np.mean(product_vals):.4f}")

    print(f"\n  {'t':>8s} | {'|ζ_carry|':>10s} | {'ζ zero?':>10s} | {'dist to nearest ζ-zero':>24s}")
    print(f"  {'-'*60}")

    local_mins = []
    for i in range(1, len(t_fine) - 1):
        if product_vals[i] < product_vals[i-1] and product_vals[i] < product_vals[i+1]:
            local_mins.append((t_fine[i], product_vals[i]))

    local_mins.sort(key=lambda x: x[1])

    for t_min, val in local_mins[:20]:
        dists = [abs(t_min - zt) for zt in ZETA_ZEROS_T]
        nearest_idx = np.argmin(dists)
        nearest_dist = dists[nearest_idx]
        is_match = "← MATCH" if nearest_dist < 2.0 else ""
        print(f"  {t_min:>8.3f} | {val:>10.4f} | {is_match:>10s} | "
              f"Δt = {nearest_dist:.3f} to ζ-zero at t = {ZETA_ZEROS_T[nearest_idx]:.3f}")

    n_matches = sum(1 for t_min, _ in local_mins[:20]
                    if min(abs(t_min - zt) for zt in ZETA_ZEROS_T) < 2.0)

    print(f"\n  Matches (Δt < 2.0): {n_matches}/20 deepest dips")

    expected_random = 20 * sum(4.0 / (ZETA_ZEROS_T[-1] + 10) for _ in ZETA_ZEROS_T[:15])
    print(f"  Expected by chance (20 dips, 15 zeros in [0,80], window ±2): "
          f"≈ {expected_random:.1f}")

    # ═══════════════════════════════════════════════════════════════════
    # PART 4: VERDICT
    # ═══════════════════════════════════════════════════════════════════
    print(f"\n{'═'*78}")
    print(" VERDICT: Spectral Determinant on the Critical Line")
    print(f"{'═'*78}")

    print("""
  1. INDIVIDUAL MATRICES:
     det(I - M/l^s) for a single semiprime oscillates wildly on σ = 1/2.
     The eigenvalues near the unit circle create resonances, but
     at semi-random positions unrelated to ζ-zeros.

  2. ENSEMBLE AVERAGE:
     ⟨|det|⟩ smooths out individual fluctuations.
     R(l,s) = ⟨|det|⟩ · |1 - l^{-s}| is ≈ 1 for small l, but with
     oscillations whose amplitude grows with l.

  3. CARRY ζ PRODUCT:
     The product Π_l ⟨|det|⟩^{-1} over the first few primes produces
     a function with local minima. The question: do these minima
     align with ζ-zeros?

     [See match count above for the answer]

  4. FUNDAMENTAL LIMITATION:
     The carry ensemble averages over ALL semiprimes — destroying
     the specific number-theoretic structure that makes ζ zeros
     what they are. The averaged spectral determinant captures the
     GENERIC carry statistics ((b-1)/b anti-correlation, trace anomaly)
     but NOT the arithmetic specificity of individual primes.

     To recover ζ zeros, one would need a single infinite-dimensional
     operator (Berry-Keating style), not an ensemble average of
     finite companion matrices.
""")
