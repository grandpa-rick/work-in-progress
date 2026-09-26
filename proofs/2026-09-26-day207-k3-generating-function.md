# Day 207 (Test 3): k=3 generating function for e_3(Y)•e_r, which telescopes completely and matches Day 193

**Date:** 2026-09-26. **Status:** a proof SKETCH. Each step carries its own grade.
- The whole chain is proved-or-symbolic **except** steps A3 and K3, which are `computed` only (Day 207 Test 1).
- The final e-basis result is `computed`, and **symbolic in r** (u = t^r).

**Scripts** (all in `scripts/day207/`):
- `k3_gf_engine.py` + `.log`: iterated Ω, run symbolically.
- `k3_gf_extract.py` + `.log`: e-basis extraction and the Day 193 comparison.
- `k3_gf_verify.py` + `.log`: check against the direct AHA pipeline.

**Conventions.** As in Day 206b §0: s = q^{-1}, a_ij = (X_i − tX_j)/(X_i − X_j), E(z) = ∏(1+X_iz), Q(v) = ∏(1−tX_iv)/(1−X_iv) = E(−tv)/E(−v), and [n] = [n]_t.

## Steps

**A3** [`computed`, m ≤ 7]. On symmetric F, e_3(Y)F = t³ σ^{(3)} π³ F.
- The Day 206b §1 per-tuple shift argument should extend verbatim. That extension is sketched, not written.

**K3** [`computed`, m ≤ 7]. For G = π³e_r, σ^{(3)}G = Σ_{|A|=3} G^{(A)} ∏^×_A.
- Expected proof: the coset identity σ_mσ'σ'' = [3]!·σ^{(3)} on G, plus Lemma 1 applied three times, mirroring Day 206b §2. Not written.

**S** [proved, classical]. Σ_{w∈S_3} ∏_{i<j} a_{w_i w_j} = [3]!.
- This is the Macdonald III.1 Poincaré identity.
- Hence Σ_{|A|=3} = [3]!^{-1} Σ_{ordered (a,b,c)} with the weight a_ab a_ac a_bc ∏^×.

**H3** [proved-by-argument]. The ordered weight factors as ∏_{j∉{a,b,c}}a_cj · ∏_{j∉{a,b}}a_bj · ∏_{j≠a}a_aj. Apply Lemma 2 three times: in w over X̂_{ab}, then in y over X̂_a, then in x over all variables.
- Each application gives Ω[u^n] = [v^n]P(v) for n ≥ 1.
- The ambient series are P = Q·h(x)h(y), then Q·h(x), then Q, where h(u,v) = (1−uv)/(1−tuv).
- The engine asserts at every stage that the input has f(0) = 0 (exponents ≥ 1), so Lemma 2 applies. Test 2 corroborates this numerically.

**P** [symbolic computation, free of r and m]. G^{(abc)} = [z^r]E(z)φ(x)φ(y)φ(w), where φ(u) = u(1+suz)/(1+uz). Apply Ω by partial fractions in u, using A/(1−cu) ↦ A(P(c)−1).
- All poles are simple, at c ∈ {−z, −tz, −t²z}.
- The mechanism is the k=2 one: the numerator (1+uz) of h(u,−z) cancels the denominator of φ. So Q(−z) only ever meets the pole at −tz, and Q(−tz) only the pole at −t²z.

**T** [proved + symbolic]. After the three stages, **the only Q-monomials that occur are the chain prefixes** 1, Q(−z), Q(−z)Q(−tz), Q(−z)Q(−tz)Q(−t²z).
- So E(z)·Q(−z)⋯Q(−t^{j−1}z) = E(t^jz) telescopes every term.
- **There is no non-telescoping remainder.**

## The generating function (T + extraction; q_n → e_n substituted)

  (1−t)(1−t²)(1−t³) · Σ_{r≥0} z^r t^{−3}e_3(Y)•e_r = Σ_{j=0}^{3} Σ_{b=0}^{3−j} M_{jb} e_b z^{b−3} E(t^j z)

Here the entry at row j, column b is M_{jb}, with [3] = 1+t+t²:

| j \ b | b=0 | b=1 | b=2 | b=3 |
|---|---|---|---|---|
| 0 | −(s−1)(st−1)(st²−1) | −s(s−1)(t−1)(st−1)[3] | −s²(s−1)(t−1)²(t+1)[3] | −s³(t−1)³(t+1)[3] |
| 1 | (s−1)(st−1)[3](st³−st²+s−t)/t³ | s(s−1)(t−1)(t+1)[3](st²−st+s−t)/t³ | s²(s−1)(t−1)²(t+1)[3]/t² | — |
| 2 | −(s−1)(s−t)[3](st³−st+s−t²)/t⁶ | −s(s−1)(s−t)(t−1)[3]/t⁴ | — | — |
| 3 | (s−1)(s−t)(s−t²)/t⁶ | — | — | — |

**Structural observations.**
- **No straightening is needed.** The q-polynomial coefficients collapse to terms *linear* in e_b.
- **The b + j ≤ 3 support is a proof of the support shape.** [z^r] gives Σ M_{jb} t^{j(r+3−b)} e_b e_{r+3−b}, so the support lies in {e_{r+3−b, b}}. This is the min(3, r)+1 meta-conjecture support at a = 3, modulo A3/K3.
- The negative z-powers cancel. This was checked numerically.

## Verification
- **(a) Against the direct pipeline** [`computed`]. The GF (with q_n computed directly from X) equals t^{−3}·e_3(Y)•e_r from the direct pipeline (`ekY_er`) at 3 random rational points for m = 6, 7 and r = 0..5: 12/12 OK. The [z^{−1,−2,−3}] parts vanish.
- **(b) Against Day 193** [symbolic in u = t^r]. The extracted coefficients of e_{r,3}, e_{r+1,2}, e_{r+2,1} and e_{r+3} equal Day 193's c_3, c_2, c_1, c_0 with s = 1/q, **for generic r**. All spurious monomials (e_re_1³, e_re_1e_2, …) cancel.
  - Since the GF holds for all r ≥ 0, this upgrades Day 193 from "r ≤ 6 fit" to "all r, conditional on A3 + K3".

## What remains for a proof of e_3⋆e_r
1. Write out A3 (the per-triple braid shift).
2. Write out K3 (the coset identity with the [3]! factor, and Lemma 1 applied three times).
3. Replace the sympy partial-fraction stage by a hand computation, or accept it as a certified symbolic identity; it is finite and free of r and m.

Watch item for k = 4: whether higher-order poles or non-chain Q-monomials appear. At k = 3 neither occurs.
