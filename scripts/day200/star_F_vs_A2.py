"""
Day 200 sharpened: does F⋆ = c(q,t) A^{(2)} on Λ_{q,t}?

F := p_2(Y) . 1  =  (1/q) m_(2) - (q+1)(t-1)/q m_(1,1)
                 =  (1/q) P_(2) - t(q+1)(t-1)/(qt-1) P_(1,1).

We test Option 1 (spectral): compute F ⋆ P_λ for small λ via Hikita's identity
    e_a(X) ⋆ G(X) = t^{-a(a-1)/2} e_a(Y) . G(X)   (bullet action of Y_i's)
and its consequence for products in the star algebra:
    e_1 ⋆ e_1 ⋆ G = e_1(Y) . (e_1(Y) . G)
    e_2 ⋆ G       = t^{-1} e_2(Y) . G

Expanding F in the e-basis (deg 2):
    F(X) = (1/q) e_1^2 - C e_2,   where
    m_(2) = e_1^2 - 2 e_2   (in ≥2 vars),   m_(1,1) = e_2
    ⇒  F = (1/q)(e_1^2 - 2 e_2) - (q+1)(t-1)/q · e_2
         = (1/q) e_1^2 - [2 + (q+1)(t-1)]/q · e_2.

So
    F ⋆ G = (1/q) e_1(Y)^2 . G   -  C/t · e_2(Y) . G
with C = [2 + (q+1)(t-1)]/q.

We check whether F ⋆ P_λ is a scalar (q,t-rational) multiple of P_λ for
each |λ| ≤ 3 (using m = 4 variables, so ≥ n+2 for the |λ|=2 test).
If yes, extract μ_F(λ).  Compare against a directly-computed A^{(2)}
eigenvalue: A^{(2)} P_λ = t^{-1} P^*_{(1,1)}(q^{-λ}; q^{-1}, t^{-1}) P_λ.

We compute P^*_{(1,1)} evaluations via the small-|μ|=2 shifted Macdonald
polynomial in TWO variables (since only λ_1, λ_2 appear in the (1,1)
shifted polynomial evaluated at (q^{-λ}; q^{-1}, t^{-1}) — the higher
parts contribute trivially only if μ has parts limited to length ≤ 2).
We use Okounkov's construction: P^*_{(1,1)} is the unique symmetric
polynomial (in the shifted sense) of degree ≤ 2 with:
   P^*_(1,1)(∅) = 0, P^*_(1,1)(μ) = 0 for μ = (1), (2), (3), ...
   P^*_(1,1)((1,1)) = 1 (or ≠ 0, we extract ratios so leading normalization
   scales into c(q,t) if any).

We verify the identity F⋆ = c(q,t) A^{(2)} by computing the ratio
    r(λ) := μ_F(λ) / eig_A2(λ)
and testing whether r is λ-independent.  Two λ suffice; three provide a
tie-breaker.

If Option 1 fails (F⋆ P_λ not scalar-multiple), we report which λ breaks
diagonality and fall back to Option 2 (e-basis coefficients of F ⋆ e_r).
"""
import sys
import os
import time
sys.path.insert(0, '/home/agent/projects/proofs/scripts/day192')
sys.path.insert(0, '/home/agent/projects/proofs/scripts/day200')

import sympy as sp
from hikita_star import (
    build_action, e_r_X, e_lambda_X, partitions_of,
    expand_symmetric_in_e_basis,
)
from A2_vs_p2Y import (
    macdonald_P_in_powersum, eval_symfun_from_p,
    sym_poly_to_monomial_basis, dominance_ge,
    p_lambda_expand_monoms, monomial_symmetric_m,
    is_scalar_multiple, degree_of_sym_poly,
    expand_in_macdonald_basis,
)

q, t = sp.symbols('q t')

T0 = time.time()
def _t():
    return f"[t={time.time()-T0:6.1f}s]"

def _pr(*args):
    sys.stdout.write(" ".join(str(a) for a in args) + "\n")
    sys.stdout.flush()


# ---- F ⋆ G via Hikita/bullet ----
def F_star_G(G, m, Y_apply):
    """Compute F ⋆ G where F = p_2(Y).1.

    By Hikita's identification (bullet action of Y-symmetric-polynomial
    corresponds to ⋆-multiplication by (that polynomial)·1 in Λ_{q,t}),
    we have F ⋆ G = p_2(Y) · G in the bullet action, i.e.,
        F ⋆ G = sum_{i=1}^m Y_i^2 · G.
    """
    total = sp.Integer(0)
    for i in range(1, m+1):
        Yi_G = Y_apply(G, i)
        Yi2_G = Y_apply(Yi_G, i)
        total = sp.expand(total + Yi2_G)
    return total


# ---- Shifted Macdonald P^*_(1,1) for evaluation at spectral point (q^{-λ}; q^{-1}, t^{-1}) ----
# We use Okounkov's inhomogeneous shifted Macdonald polynomials.
# For μ = (1,1), P^*_(1,1)(x_1,...,x_N; q,t) is a symmetric polynomial in
# shifted variables such that its evaluation at the spectral point of λ is
# the standard Macdonald "spectral function" e^*_(1,1)(λ). For small λ we
# can just directly compute the A^{(2)} eigenvalue from the known
# Nazarov-Sklyanin formula in the p-basis form.
#
# Thibon 2608.30791 Thm 2.3 (as quoted in the task) says:
#   A^{(2)} P_λ = t^{-1} P^*_(1,1)(q^{-λ}; q^{-1}, t^{-1}) P_λ.
#
# We can compute the "spectral function" of P^*_(1,1) at (q^{-λ_1}, q^{-λ_2}, ...; q^{-1}, t^{-1})
# by writing P^*_(1,1) explicitly in shifted-power-sum form. But an
# equivalent (and simpler) recipe: since P^*_(1,1) is degree 2 and
# symmetric, one form in TWO variables (x_1, x_2) is:
#   P^*_(1,1)(x_1, x_2; q, t) = (x_1 - 1)(x_2 - t^{-1})  + (x_2 - 1)(x_1 - t^{-1}) - const,
# adjusted so that it vanishes on μ = (), (1), (2), (3), ...
#
# Actually, the clean formulation:
#   For μ ⊂ (single row), the spectral point ρ(μ) = (q^{μ_1}, t^{-1}, t^{-2}, ...).
#   P^*_(1,1) vanishes iff (1,1) ⊄ μ, i.e., μ has at most one part.
#   P^*_(1,1) at ρ((1,1)) = (q, q · t^{-1}) or similar — depends on convention.
#
# GIVEN the CONVENTION mismatch is a nightmare in shifted-Macdonald-land, we
# ALTERNATIVELY compute the A^{(2)} eigenvalue via the Nazarov-Sklyanin
# construction directly.  A^{(2)} on Λ acts by known formula in the p-basis
# (Thibon eq. (10) or so):
#    A^{(2)} = (1/2) [p_1^2 * partial_{p_1}^2 ... ]  -- but this is level-(2)
# For the purpose of THIS task we shortcut: rather than compute A^{(2)} eig
# from scratch, we exploit the fact that
#
#   IF F⋆ = c(q,t) A^{(2)}, THEN the ratio μ_F(λ_1) / μ_F(λ_2) for two
#   different λ_1, λ_2 must equal eig_A2(λ_1) / eig_A2(λ_2).
#
# So we compute μ_F for several λ and then check if the ratios pattern-match
# a Macdonald spectral function of type (1,1) — i.e., a spectral function
# whose only vanishings are at partitions with ≤ 1 row.
#
# Concretely: eig_A2(∅) = 0, eig_A2((1)) = 0, eig_A2((r)) = 0 for all r,
# and eig_A2(λ) ≠ 0 iff ℓ(λ) ≥ 2.
#
# So the CORE spectral diagnostic is:
#   Does μ_F vanish exactly on single-row partitions?
#
# Then compare μ_F((1,1)) to a directly-computed A^{(2)} eigenvalue on P_(1,1)
# to extract c(q,t).

# ---- Direct A^{(2)} eigenvalue on P_λ via known Thibon Thm 2.3 form ----
#
# For μ = (1,1), the shifted Macdonald P^*_(1,1)(x; q,t) can be normalized as
#   P^*_(1,1)(x_1, x_2, ...; q, t)  =  \sum_{i<j} (x_i - 1)(x_j - t^{-1}) + ...
# One convenient form (Okounkov-Olshanski, symmetric functions of shifted vars):
#   P^*_(1,1)(x; q, t) = e_2(x - 1) - some correction
# For SMALL partitions we compute the spectral function directly by:
#   A^{(2)} eig on P_λ  =  \sum_{(i,j): i<j, λ_j ≥ 1}  q^{λ_i - λ_j}·(t^{i - j}-terms)
# This is complicated; instead we use the following EMPIRICAL BENCHMARK:
#
#   The Macdonald DELTA_(1,1) operator on Λ_{q,t} has known eigenvalues:
#      Δ_{(1,1)} P_λ = e_2(q^{λ_i} t^{-i+1}; i=1..∞ / suitably regularized) · P_λ
#   with some subtraction to make it finite. In HIKITA's ⋆-algebra, the
#   ⋆-multiplication by P_(1,1) gives eigenvalue = P_(1,1) evaluated at the
#   spectral point of λ.
#
# The clean, tool-free identification: F ⋆ P_λ = μ_F(λ) P_λ, and since we
# EXPLICITLY have F = (1/q) P_(2) - t(q+1)(t-1)/(qt-1) P_(1,1), the
# ⋆-multiplication by F is
#     μ_F(λ) = (1/q) · e_(2)(λ) - t(q+1)(t-1)/(qt-1) · e_(1,1)(λ)
# where e_μ(λ) is the ⋆-eigenvalue of P_μ on P_λ, i.e., P_μ evaluated at
# the SPECTRAL POINT of λ.
#
# We could then compare to A^{(2)} eigenvalues, but the question actually
# collapses to: IS F ⋆ P_λ scalar-multiple of P_λ? (i.e., is the ⋆-product
# diagonal on Macdonald in Hikita's setup at level 1 with m variables?)
#
# THIS is a first-principles computational test we CAN run.

# ---- Main test ----

def analyze():
    m = 4  # Use 4 variables so partitions up to size 3 have room to breathe
    _pr(f"{_t()} Building Y-action for m = {m} ...")
    X, Ti_apply, Ti_inv_apply, Pi_apply, Y_apply = build_action(m)
    _pr(f"{_t()} Y-action built.")

    # Sanity: verify F ⋆ 1 = F.
    _pr(f"\n{_t()} Sanity: computing F ⋆ 1 (should equal F).")
    F_star_1 = F_star_G(sp.Integer(1), m, Y_apply)
    _pr(f"{_t()}   F ⋆ 1 = {sp.factor(F_star_1)}")
    F1_m = sym_poly_to_monomial_basis(F_star_1, 2, X)
    _pr(f"    m_(2) coeff:  {sp.factor(F1_m.get((2,), 0))}")
    _pr(f"    m_(1,1) coeff:  {sp.factor(F1_m.get((1,1), 0))}")
    # Expected: (1/q, -(q+1)(t-1)/q)
    exp_c2 = sp.Rational(1)/q
    exp_c11 = -(q+1)*(t-1)/q
    match_c2 = sp.simplify(F1_m.get((2,), 0) - exp_c2) == 0
    match_c11 = sp.simplify(F1_m.get((1,1), 0) - exp_c11) == 0
    _pr(f"    match m_(2): {match_c2} ; match m_(1,1): {match_c11}")

    # ---- Compute F ⋆ P_λ for each small λ ----
    _pr(f"\n{_t()} === Spectral diagnostic: F ⋆ P_λ vs P_λ ===")

    # For n = |λ|, we build the Macdonald basis in degree n (need m >= n),
    # compute F ⋆ P_λ (which has X-degree n+2 since F is degree 2 and ⋆ preserves
    # degree in the *derived* sense — wait: at finite m, we should check).
    #
    # Actually in Hikita's setup ⋆ is degree-PRESERVING: F ⋆ G lives in same
    # degree as ... hmm, no: F ⋆ G has X-degree = deg F + deg G in general.
    # Let me reconsider: ⋆ is a NEW multiplication on Λ_{q,t} that DOES
    # preserve grading (Λ_{q,t} = ⊕ Λ^d).  For deg(F)=2, deg(G)=n, F⋆G has
    # degree n+2.  But then F⋆P_λ CANNOT be a scalar times P_λ (they have
    # different degrees).  So the eigenvalue interpretation is wrong at
    # finite m.
    #
    # UNLESS Hikita's setup identifies P_λ ∈ Λ^n with something in Λ^{n+2}
    # via ⋆-multiplication.  Let me re-read the task statement...
    #
    # Task: "If it's a scalar multiple of P_λ (i.e., F is central-like),
    # extract the eigenvalue μ_F(λ). Compare against ..."
    #
    # So the task ASSUMES F ⋆ P_λ = μ_F(λ) P_λ (same degree). This requires
    # a degree-shift somewhere. Let's just compute and see what happens.

    # We'll compute F ⋆ P_λ for λ ∈ {∅, (1), (2), (1,1), (3), (2,1), (1,1,1)}
    # and try to express in Macdonald basis of degree |λ|+2 (which is the correct
    # target degree at finite m). Then check whether the coefficient supported on
    # P_λ (padded to degree |λ|+2)... but the "same-degree" idea doesn't apply
    # unless we have some way to promote P_λ to degree |λ|+2.
    #
    # SO the honest test is: expand F ⋆ P_λ in Macdonald basis of degree |λ|+2
    # and see if only ONE Macdonald poly appears (indicating a very special
    # structure — but that would mean F ⋆ P_λ = (constant) · P_ν for some
    # ν of size |λ|+2, not P_λ itself).

    max_n = 2  # partitions of size 0,1,2 for eigenvalue test
    for n in range(0, max_n + 1):
        parts_n = partitions_of(n)
        _pr(f"\n{_t()} --- Testing |λ| = {n}, partitions = {parts_n} ---")

        # Build Macdonald basis in degree n and in degree n+2
        _pr(f"{_t()}   Building Macdonald basis in degrees {n} and {n+2} ...")
        Pdict_n = macdonald_P_in_powersum(n, X) if n > 0 else {(): {(): sp.Integer(1)}}
        Pdict_n2 = macdonald_P_in_powersum(n+2, X)
        _pr(f"{_t()}   Macdonald basis built.")

        for lam in parts_n:
            if len(lam) > m:
                continue
            _pr(f"\n{_t()}   λ = {lam}")
            if n == 0:
                Pl_poly = sp.Integer(1)
            else:
                Pl_poly = eval_symfun_from_p(Pdict_n[lam], X)
            _pr(f"{_t()}     P_λ = {sp.factor(Pl_poly)}")

            _pr(f"{_t()}     Computing F ⋆ P_λ ...")
            FstarP = F_star_G(Pl_poly, m, Y_apply)
            deg = degree_of_sym_poly(FstarP, X)
            _pr(f"{_t()}     deg(F ⋆ P_λ) = {deg}  (expected {n+2})")

            if deg == n + 2:
                _pr(f"{_t()}     Expanding F ⋆ P_λ in Macdonald basis of deg {n+2} ...")
                coeffs, residual = expand_in_macdonald_basis(FstarP, n+2, X, Pdict_n2)
                nonzero = [(mu, c) for mu, c in coeffs.items() if sp.simplify(c) != 0]
                _pr(f"{_t()}     # nonzero P_μ coefficients: {len(nonzero)}")
                for mu, c in nonzero:
                    _pr(f"       P_{mu}:  {sp.factor(sp.cancel(c))}")

    # ---- Section: e-basis diagnostic (Option 2 fallback) ----
    _pr(f"\n{_t()} === Option-2 diagnostic: F ⋆ e_r in e-basis ===")
    for r in [0, 1, 2]:
        _pr(f"\n{_t()}   r = {r}")
        if r == 0:
            G = sp.Integer(1)
        else:
            G = e_r_X(m, r)
        FstarG = F_star_G(G, m, Y_apply)
        n = r + 2
        _pr(f"{_t()}     Expanding F ⋆ e_{r} in e-basis of degree {n} ...")
        expansion = expand_symmetric_in_e_basis(FstarG, m, n)
        for lam in partitions_of(n):
            c = expansion.get(lam, sp.Integer(0))
            cs = sp.simplify(c)
            if cs != 0:
                _pr(f"       e_{lam}: {sp.factor(cs)}")


def main():
    _pr("=" * 78)
    _pr(f"Day 200 sharpened: does F⋆ = c(q,t) A^{{(2)}} on Λ_{{q,t}}?")
    _pr(f"  F = (1/q) e_1^2  -  [2+(q+1)(t-1)]/q · e_2  (= p_2(Y).1)")
    _pr("=" * 78)
    analyze()


if __name__ == "__main__":
    main()
