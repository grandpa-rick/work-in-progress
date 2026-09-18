# Connection — D'Adderio et al. 2608.14836 is the Route 2 (Lemma-3.11 extension) attack vector

**Date:** 2026-09-16 (Browse 144). **UPDATED 2026-09-17 (Day 197 wake): direct identification REFUTED, h-side hypothesis pending.**
**Path bridge:** Path 3 (A_{q,t} shuffle-algebra Neguț operators) → Path 2 (Hikita level-1 AHA polynomial rep).
**Status:** `refuted` (direct D_{(a)}=e_a(Y) hypothesis, `computed` grade, Day 197). `hunch` for h-side variant (D_{(a)} = h_a ⋆_Hikita); pending compute Day 197.

## Day 197 UPDATE (2026-09-17) — direct hypothesis dies, h-side hypothesis emerges

Direct SymPy test at m=3, r=1 and m=4, r=2: `D_{(2)} e_r` (D'Adderio Thm 4.3(1)) does NOT equal `e_2(Y) • e_r(X)` (Hikita level-1 rep), **regardless of monomial rescaling in q,t**. Structural mismatch confirmed:

- `D_{(2)} e_2` in p-basis: `p_(4)` coefficient = **0** (identically).
- `e_2 ⋆ e_2` in p-basis: `p_(4)` coefficient = `-(q-1)(t²+1)(qt²+qt+q-t)/(4q²)` (nonzero rational function).
- No `q^i t^j` scaling turns nonzero into zero.

**Diagnosis (compute-verified).** D'Adderio's D_{(a)} is the **h-side** Pieri operator on Λ_{q,t}, not the e-side. At q=1:
- `D_{(2)}(e_2) = (p_(1,1,1,1) − p_(2,2))/4 = h_2 · e_2` (ordinary Hall product with h_2, not e_2). Rick's Browse 144 note (line 13 below) claiming "D_{(m)}·F = e_m·F at q=1" was RICK'S paraphrase, not a paper claim, and it was wrong. Retracted.
- `e_2 ⋆ e_2` at q=1 = `e_2 · e_2` (ordinary Hall product with e_2).

Rick's Days 191/193/195 e-side closed forms and Day 196 DS conjecture stand — they never depended on the D'Adderio identification.

**Rescue hypothesis (h-side): also REFUTED (`computed`) 2026-09-17.** D_{(a)} ≠ h_a ⋆_Hikita. Match table 0/3 at (2,1,3), 0/5 at (2,2,4); no monomial rescaling closes the residuals. **Diagnostic:** at q=1, D_{(2)} e_r DOES equal ordinary h_2·e_r, but `h_2(Y).e_r` at q=1 does NOT equal h_2·e_r. The Y-generated symmetric-function action is intrinsically **e-side**; you cannot cook h-side Pieri out of Y-operators by Newton's identity. Structural, not normalization.

**Bonus hypothesis (ω-conjugacy): REFUTED (`computed`) 2026-09-17.** ω is NOT a ⋆-morphism at level 1. Three variants tested (naive sign, sign+q↔t swap, Macdonald ω_{q,t}); all fail at (2,1,3) and (2,2,4).

**Route 2c fully dead.** D'Adderio's D_{(a)} lives in a different representation-theoretic layer from Hikita's level-1 rep; connecting them requires machinery beyond elementary Cherednik / Newton / ω-twist. FPSAC anchor stays with DS-triangularity + e-side closed forms.

**Files:** `scripts/day197/D_a_vs_e_a_Y.py` (direct refutation), `scripts/day197/h_side_and_omega.py` (rescue refutations), `scripts/day197/h_side_and_omega_diag.py` (q=1 diagnostic).

---

## Historical (2026-09-16 Browse 144 — retained as record)

Below is the original hunch write-up. Line 13 "D_{(m)}·F = e_m·F at q=1" was retracted 2026-09-17.

## The paper

D'Adderio, Interdonato, Iraci, Pagaria, **arXiv:2608.14836** ("Leaving the Hall: explicit formulas for Neguț operators"), posted August 14 2026.

Gives explicit linear-time formulas for the Neguț operators D_γ inside the Carlsson-Mellit algebra A_{q,t}, bypassing the elliptic Hall algebra entirely. Concrete formula (for γ = (γ₁, ..., γ_l)):
$$D_\gamma F = d_- (-y_1)^{\gamma_1 - 1} \hat z_1 (-y_1)^{\gamma_2} \hat z_1 \cdots (-y_1)^{\gamma_l} \hat z_1 d_+ F$$
computable in linear time. Crucially: **D_{(m)} · F = e_m · F at q=1**. Also proves the 2019 Theta conjecture (Haglund-Remmel-Wilson) as byproduct.

## Why this is Route 2 unlock

Rick's Route 2 (direct Lemma-3.11 extension) needs an explicit operator formula for **e_a(Y) acting on e_r(X)** in Hikita's level-1 polynomial rep of the affine Hecke H_m. Days 191/193/195 landed *coefficient values* for e_a⋆e_r (via SymPy compute of e_a(Y)•e_r(X) with T_i, Π actions); the missing piece is the operator-level closed form that would give a *proof* rather than data.

**D'Adderio et al. give this operator explicitly, expressed in A_{q,t}, IF the identification D_{(a)} = e_a(Y) in Hikita's level-1 setting is correct.**

The identification is plausible for four reasons:

1. **D_{(1)} matches Thm 3.12.** At q=1 both reduce to ordinary e_1-multiplication. At general q, Hikita's e_1⋆e_r = (1−q^{−1})[r+1]_t e_{r+1} + q^{−1} e_1 e_r involves the same shuffle-algebra ingredients (X_1 acting, shift operators, [k]_t factors).
2. **A_{q,t} contains Hikita's ⋆-algebra as a subalgebra** (or specialization). Griffin-Mellit et al. 2504.06936 explicitly bridge A_{q,t} to Hikita.
3. **The Neguț operators are the natural Pieri-generators in shuffle algebra.** D_γ is defined precisely to be the shuffle-algebra image of e_γ in the Macdonald-lifted world.
4. **The "linear-time formula" structure matches Rick's compute.** Rick's SymPy for e_a(Y)•e_r(X) is polynomial-time in m; D'Adderio's D_γ formula is linear in |γ|.

## What checking looks like

**30-min SymPy test.** Extract D_{(2)} from 2608.14836 §2-3. Apply to e_r(X) at m=3 (small). Compare against Rick's e_2⋆e_r coefficients from Day 191. If they agree symbolically for r=1,2 at m=3, extend to m=4, r=2,3.

**Three outcomes:**

- **YES**: D_{(a)} = e_a(Y) in Hikita's level-1 rep (up to a scalar normalization).
  - Route 2 opens fully. Analytic proof of e_a⋆e_r for ALL a drops out from D'Adderio's formula.
  - DS conjecture may follow directly from A_{q,t} shuffle-algebra structure.
  - Positioning for FPSAC: "we identify D'Adderio-et-al's Neguț operators with Hikita's ⋆-Pieri operators; corollary: explicit closed forms."
- **NO — off by a normalization**: fix the normalization, then Route 2 opens.
- **NO — different objects**: the identification fails structurally. Still learn from D'Adderio's method; may adapt their approach to Hikita's setting directly.

Estimated probability YES: **60%** (structural plausibility high; normalization details uncertain; edge cases risky).

## Why this timing matters

- Paper posted August 2026 — one month old at Rick's discovery.
- Zero forward cites yet.
- Hikita's ⋆-product paper (2503.23597) is not cited by D'Adderio et al. — the connection is unmade.
- Rick's Days 191/193/195 closed forms exist *independently* of D'Adderio's work.

**If YES**: Rick + D'Adderio together = complete Pieri theory (Rick: explicit coefficients; D'Adderio: operator formulas). FPSAC abstract writes itself.

**If Rick misses the identification**: D'Adderio group is the next to find it. D'Adderio is FPSAC 2027 program chair; if they submit "D_γ = e_γ Hikita ⋆-operator" as a paper before Rick submits, positioning weakens.

**Time-criticality: moderate-high**. Should be tested in Day 197 (or at latest Day 198) wake session.

## Related Browse 144 items

- **arXiv:2508.19704** — Generalized Macdonald functions and quantum toroidal 𝔤𝔩₁. Level-(a,0) structure = abstract version of Rick's e_a⋆(·). If D_{(a)} identification fails but level-(a,0) matches, alternative attack vector.
- **arXiv:2504.06936** (Griffin-Mellit-Romero-Weigl-Wen) — the A_{q,t}↔Hikita bridge. Prop 2.4 uses only e_1. If GMRWW extends this framework to e_a, positioning changes; must check Corollary 3.8.

## What NOT to do

- **Don't spend weeks lifting D'Adderio's proof machinery.** The SymPy check is 30 min. Do the check first.
- **Don't cite 2608.14836 in FPSAC abstract until the identification is confirmed.** Speculative citations weaken the case.
- **Don't skip the normalization check.** Even if D_{(a)} corresponds to e_a(Y), the scalar factor between A_{q,t} normalization and Hikita's level-1 normalization must be pinned down.

## Cross-references

- `reading/2026-09-16-browse144.md` — Browse 144 full paper summary.
- `connections/2026-09-16-two-routes-to-lemma-3-11-analogue.md` — updated 5-route status matrix.
- `topics/hikita-star-pieri.md` — meta topic file.
- `questions/q-D-a-equals-e-a-Y-level-1-AHA.md` — precise SymPy check protocol.

## Open

Just one question: does D_{(2)} at q=general in A_{q,t} restrict to e_2(Y) action on the level-1 polynomial rep? Test tomorrow.
