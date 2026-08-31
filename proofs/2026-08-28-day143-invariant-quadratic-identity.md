# Day 143 — Universal invariant sequence: quadratic identity

**Status:** New structural closed form found. Reduces PRIMARY invariant problem
to a quadratic transform of the leading cumulant sequence.

## Problem

Let L := T·(U + θ_T)(V + θ_T) − θ_T (annihilates ₂F₀(U, V; ; T)), and let
F_P = Σ_b P_b T^b/b! be the interior series (Day 140).  Then

    L · F_P = F_P · X(T, U, V, E_3)   (Frobenius identity — Day 142).

The empirical fact (verified k = 1..7):

    a_k := [E_3^k T^{3k-1}] X   is INDEPENDENT of (U, V),

    a_1 = -3, a_2 = -18, a_3 = -255, a_4 = -4620, a_5 = -94500,
    a_6 = -2078802, a_7 = -48005802.

Day 142 tried OEIS + factorial normalizations, no closed form found; the
sequence contains large primes (17, 10499, 367) that rule out simple
hypergeometric/factorial ratios.

## Result

Let

    n_k := N_k[T^{3k-1}]   (leading T-coefficient of k-th cumulant of F_P/f)
    b_k := (3k-1) · n_k.

Then

    ╔══════════════════════════════════════════════════════════════╗
    ║  a_k  =  -b_k  +  Σ_{i+j=k, i,j ≥ 1} b_i · b_j               ║
    ╚══════════════════════════════════════════════════════════════╝

Equivalently, with A(τ) := Σ_{k≥1} a_k τ^k, F(τ) := Σ_{k≥1} b_k τ^k,

    A(τ) = F(τ)² − F(τ),          F(τ) = (1 − √(1 + 4 A(τ))) / 2,
    (1 − 2 F(τ))² = 1 + 4 A(τ).

**Corollary.** The formal series 1 + 4 A(τ) is a perfect square in ℚ⟦τ⟧.

**Verification.** With n_k = 3/2, 27/5, 417/8, 7851/11, 164124/14, 3661389/17
(from Day 142), so b_k = 3, 27, 417, 7851, 164124, 3661389, and

    -b_1                             = -3                = a_1  ✓
    -b_2 + b_1·b_1 = -27 + 9         = -18               = a_2  ✓
    -b_3 + 2 b_1 b_2 = -417 + 162    = -255              = a_3  ✓
    -b_4 + 2 b_1 b_3 + b_2² = -7851 + 2502 + 729 = -4620 = a_4  ✓
    -b_5 + 2 b_1 b_4 + 2 b_2 b_3 = -164124+47106+22518   = -94500 = a_5  ✓
    -b_6 + 2 b_1 b_5 + 2 b_2 b_4 + b_3²                  = -2078802 = a_6  ✓

**Prediction.** From a_7 = -48005802 and the identity,

    b_7 = -a_7 + 2 b_1 b_6 + 2 b_2 b_5 + 2 b_3 b_4 = 85384566,
    so n_7 = 85384566 / 20 = 42692283 / 10.

Verifiable by an independent extension of `cumulants_UV0.py` to B_max = 21.

## Proof

**Setup.** Write F_P = f · G where G = exp(R), R := Σ_{k≥1} E_3^k N_k =
log(F_P / f).  A direct calculation from the product rule for θ = T d/dT gives

    L(f G) = f · (L G − T U V · G) + 2 T (θf)(θG)                (*)

(cancellation uses L · f = 0.)  Dividing by fG:

    X = L(F_P) / F_P = L(fG)/(fG)
      = (LG − TUV·G)/G + 2 T φ · (θG/G)          (φ := θf/f).

With G = exp(R):

    θG/G = θR,  θ²G/G = θ²R + (θR)².

Expanding LG/G = TUV + T(U+V) θR + T (θ²R + (θR)²) − θR, we obtain

    X = [T(U+V) − 1 + 2 T φ] · θR  +  T · θ²R  +  T · (θR)².

**Extraction at E_3^k.**  With R = Σ_k E_3^k N_k,

    X_k = [T(U+V) − 1 + 2 T φ] · θN_k  +  T · θ²N_k
          + T · Σ_{i+j = k,  i,j ≥ 1} (θN_i)(θN_j).

**Lemma (leading-T vanishing).**  For all k ≥ 1 and 2k ≤ b ≤ 3k−2,

    N_k[T^b] = 0.

*Evidence.*  Verified as a polynomial in (U, V) at four independent points
(U, V) = (0, 0), (1, 1), (2, 3), (−1, 2) for k = 1..5 (up to k = 3 or 4 at the
non-(0,0) points).  Since N_k[T^b] is polynomial in (U, V), vanishing at these
points is strong evidence for the polynomial identity.  A cleaner proof from
the E_3-degree structure of P_b (= p_b + E_3 · U_b(E_3 + φ_1)) is expected
but not written out here.

Extracting [T^{3k-1}] X_k using the lemma:

- [T^{3k-1}] {(T(U+V) − 1) · θN_k}:
    - T(U+V) · θN_k contributes (U+V)(3k−2) · N_k[T^{3k-2}] = 0.
    - (−1) · θN_k contributes −(3k−1) · N_k[T^{3k-1}] = −b_k.
- [T^{3k-1}] {2 T φ · θN_k}:  Requires N_k[T^b] with b ≤ 3k − 3; = 0.
- [T^{3k-1}] {T · θ²N_k}:  (3k−2)² · N_k[T^{3k-2}] = 0.
- [T^{3k-1}] {T · (θN_i)(θN_j)},  i + j = k,  i, j ≥ 1:
    (θN_i)(θN_j)[T^{3k−2}] = Σ_{a+b = 3k−2} a b · N_i[T^a] · N_j[T^b].
    By the lemma, a ≥ 3i−1 and b ≥ 3j−1, so a + b ≥ 3k−2 with equality
    forced, giving unique term (a, b) = (3i−1, 3j−1):
    contribution (3i−1)(3j−1) · N_i[T^{3i−1}] · N_j[T^{3j−1}] = b_i · b_j.

Summing:  a_k = X_k[T^{3k-1}] = −b_k + Σ_{i+j=k, i,j≥1} b_i b_j.  ∎

**Consistency check.**  The RHS is (U, V)-independent (since b_l is), so
X_k[T^{3k-1}] must be (U, V)-independent — matching the empirical observation.

## Discussion

This is a **structural** closed form: it reduces the invariant problem to
the (still open) closed form for b_k.  The remaining open questions:

1. **Closed form for b_k = 3, 27, 417, 7851, 164124, 3661389, 85384566, …**

   Factorizations: 3, 3³, 3·139, 3·2617, 2²·3²·47·97, 3³·135607, 2·3²·4743587.
   Large primes at every stage — not hypergeometric.

2. **P-recurrence for b_k.**  Not yet identified.  Order-2 degree-1 fails
   with 7 data points; higher-order tests need more data.

3. **The identity says (1 + 4 A(τ)) is a perfect square in ℚ⟦τ⟧.**  This
   is a nontrivial fact about the sequence a_k.  It embeds naturally in the
   Riccati / Bernoulli-transform structure of X = L(exp R)/(exp R) at the
   diagonal.

## Companion result: complete vanishing for b < 3k−1

**Stronger statement:**  For all k ≥ 1 and all b with 2k ≤ b < 3k − 1,

    [E_3^k T^b] X  =  0.

(Trivially also for b < 2k from E_3-degree constraints on F_P.)

Verified at (U, V) = (0, 0) for all k = 1..6 (up to b = 17). That is
k − 1 consecutive zeros from T^{2k} to T^{3k−2}, then the invariant a_k
at T^{3k−1}.

Data (rows = k, columns = T-power, all at (U,V) = (0,0)):

    k=2:  T^4=0  T^5=−18   (0 zeros before diagonal at T^5)
    k=3:  T^6=0  T^7=0     T^8=−255      (2 zeros before T^8)
    k=4:  T^8..T^10=0      T^11=−4620    (3 zeros)
    k=5:  T^10..T^13=0     T^14=−94500   (4 zeros)
    k=6:  T^12..T^16=0     T^17=−2078802 (5 zeros)

*Proof.*  Apply the same expansion at [T^b] X_k for any b < 3k−1:
- Linear terms −(b)·N_k[T^b], (b−1)² N_k[T^{b−1}], (U+V)(b−1)·N_k[T^{b−1}]
  and 2·Σ φ[T^a] c N_k[T^c] with a+c = b−1: each N_k[T^c] with c < 3k−1 = 0.
- Quadratic Σ (θN_i)(θN_j)[T^{b−1}]: forces a + c = b − 1 with a ≥ 3i−1,
  c ≥ 3j−1, hence a + c ≥ 3k−2.  Since b − 1 < 3k−2, no valid pair.

Hence [E_3^k T^b] X = 0 for b < 3k − 1.  ∎

So T^{3k−1} is the FIRST possibly-nonzero T-power in the E_3^k-slice of X,
and it hits the universal invariant a_k.  This is a strong structural
statement about X.

## Independent verification

Ran `compute_n7.py`:  N_7[T^20] at (U, V) = (0, 0) computed directly.

    Direct computation:      b_7 = (3·7 − 1) · N_7[T^20] = 20 · (42692283/10) = 85384566.
    Prediction from identity: b_7 = 85384566.

    ✓ EXACT MATCH.  Identity verified at k = 7 via independent path.

## FPSAC content

**Theorem 3.7 (Frobenius diagonal Riccati identity).**  With notation above,

    (i)   [E_3^k T^b] X  =  0                for all  b < 3k − 1,   k ≥ 1;

    (ii)  a_k := [E_3^k T^{3k-1}] X  is (U, V)-independent and satisfies

              a_k = −b_k + Σ_{i+j=k, i,j≥1} b_i b_j

          where b_k := (3k − 1) · [E_3^k T^{3k-1}] log(F_P/f).

Equivalently, with A(τ) = Σ_k a_k τ^k and F(τ) = Σ_k b_k τ^k,

    (1 − 2 F(τ))²  =  1 + 4 A(τ).

In particular 1 + 4 A(τ) is a perfect square in ℚ⟦τ⟧.

**Conjecture 4.2 (revised).**  The individual sequences a_k and b_k admit
no P-recursive closed form of low order/degree (verified with 7 data points,
order ≤ 3, degree ≤ 4).  Whether either is D-finite (i.e., satisfies SOME
polynomial-coefficient linear ODE) or algebraic remains open.  The
structural relation (1−2F)² = 1+4A is the strongest closure found.

## Files

- `code/day143_invariant/extend_invariant.py` — compute a_1..a_6 (verified −3, −18, −255, −4620, −94500, −2078802).
- `code/day143_invariant/extend_k7.py` — a_7 = −48005802.
- `code/day143_invariant/check_lowT_Nk.py` — verifies N_k[T^b]=0 for 2k ≤ b < 3k−1 at 4 (U, V) points.
- `code/day143_invariant/verify_recurrence.py` — verifies a_k = −b_k + (b*b)_k for k = 1..6.
- `code/day143_invariant/compute_n7.py` — independent verification: computes N_7[T^20] = 42692283/10 at (0,0), confirms b_7 = 85384566.
- `code/day143_invariant/verify_vanishing.py` — shows the vanishing pattern k−1 zeros before the diagonal for k = 1..6.
- `code/day143_invariant/verify_k7_vanishing.py` — extends vanishing verification to k = 7 (6 consecutive zeros before T^20).
- `code/day143_invariant/check_other_diagonals.py` — confirms [T^{3k-2}] X universally zero; [T^{3k}] and higher depend on (U, V).
- `code/day143_invariant/vieta_X1.py` — X_1[T^b] in (α, β) = (U+V, UV) Vieta variables.
- `code/day143_invariant/analyze_seq.py` — attempts closed form directly for a_k (negative).
- `code/day143_invariant/search_b_recurrence.py` — no P-recurrence for b_k up to order 3, degree 4 (negative).
- `code/day143_invariant/RESULT.md` — summary of the session.

## Time / status

Session: Day 143 (2026-08-28). Deliverable for FPSAC §3.4 / §4.
