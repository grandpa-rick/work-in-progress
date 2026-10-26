# PROVE 258 — (★) on LG without Buch–Ravikumar

Rick, 2026-10-26 (file dated for the PROVE slot 2026-10-27, as PROVE.md asked).
Draft: `writing/lg-og-negative-degree/main.tex`. Registry: `proofs/registry/q488-converse.json`, node `q488-oc-star-BR-free`.

## Problem

On LG(m,2m) put s_x(κ) = (-1)^{|κ|} H_κ(-x), x ∈ ℂ. Equation (★) (`prop:star`, `eq:star`):

    s_{x+1}(κ) = s_x(κ) + Σ_{ρ ⊊ κ} V(D(κ)∖D(ρ)) s_x(ρ)          for all x ∈ ℂ,

V(θ) = #{strictly order-preserving θ → {e<o}, diagonal boxes ↦ e}. Equivalently (`lem:layers`(a), pure combinatorics)
s_{x+1} = A_oA_e s_x.

As printed, (★) is proved from `eq:hyperplane` + `thm:Q30`, and `thm:Q30` stands on `prop:BRp1` (BR Cor 5.6 at p=1,
derived from BR Thm 5.5) + BR Cor 5.4. Goal: prove (★) from (G1) + (LS) + `lem:peel` only, then fix the dependency
graph of everything downstream.

## Solution

**Inputs.**
1. H_κ is a polynomial in j (Snapper/Kleiman; also `lem:P` given (G1)). H_∅ ≡ 1.
2. `thm:main`(b), proved in §3 (`sec:ls`) from (G1) + (LS) + Ehrhart–Macdonald reciprocity, with no (BR) and no (G0):
   s_j(κ) = N_κ(j) for **every** integer j ≥ 1 and **every** strict κ with parts ≤ m. For κ = ∅ this reads 1 = 1
   (N_∅ ≡ 1 by definition).
3. `lem:peel` (pure counting): for j ≥ 1, N_κ(j+1) = N_κ(j) + Σ_{ρ⊊κ} V(θ_ρ) N_ρ(j).

**Proposition.** On LG, (★) holds for all x ∈ ℂ, using only inputs 1–3.

*Proof.* Fix κ and put

    Δ_κ(x) = s_{x+1}(κ) − s_x(κ) − Σ_{ρ⊊κ} V(θ_ρ) s_x(ρ).

The sum is finite (at most 2^m terms), the V(θ_ρ) are integers independent of x, and each s_·(ρ) is a polynomial by input 1.
So Δ_κ is a polynomial in x. For an integer j ≥ 1 we have j+1 ≥ 1 and j ≥ 1, so input 2 applies to every term:
Δ_κ(j) = N_κ(j+1) − N_κ(j) − Σ V(θ_ρ) N_ρ(j), which is 0 by input 3. A polynomial with infinitely many zeros is 0, so
Δ_κ ≡ 0. ∎

**The two caveats from PROVE.md, checked.**
- (i) Polynomiality on the whole range: s_x is *defined* as (-1)^{|κ|}H_κ(-x) for all x ∈ ℂ, including half-integers
  (sec:half opening line uses the same definition). So it is a polynomial everywhere; no separate half-integer definition
  exists to worry about.
- (ii) The integer range: `thm:main`(b) is stated and proved for all j ≥ 1 with no κ-dependent lower bound — the window
  {1,…,⌈r/2⌉} where N_κ(j) = 0 is *included* (the count is 0 and so is s_j, by the same theorem). So for each κ the
  zero set of Δ_κ contains all of ℤ_{≥1}. Note `lem:peel` is also stated for all j ≥ 1. No gap.

**Non-circularity.** The draft's "second proof of `thm:main`" uses (★) + `lem:peel` to prove `thm:main`(b). That proof must
keep using the K-theoretic proof of (★). The new proof of (★) uses `thm:main`(b) from §3. Two proofs of (★), each
non-circular; they are just run in opposite directions, as with (a)⇔(b) already.

## Downstream: what is now BR-free on LG

Trace (every \ref read in the draft):

| result | uses | BR-free now? |
|---|---|---|
| `lem:layers`(a) | (★) + combinatorics | yes |
| `lem:layers`(b), `lem:eps` | combinatorics + Thm A (Demazure, (G1)) | yes |
| `cor:row` | (★) at half-integers, layers, eps | yes |
| eq:fibre, `lem:shift`, `lem:onerow` | combinatorics | yes |
| `thm:HZ` (i)–(iii) | cor:row, onerow, `lem:Npos`, Thm A | **yes** |
| `cor:OC` (= Thm B) | HZ + lem:chain | **yes** |
| `prop:T`(a) | (★) via A_oA_e, `lem:anchor`, Thm A, `lem:mobius`, `lem:XP` | **yes** |
| `prop:T`(b) | matrix combinatorics | yes (was already) |
| `thm:Z` as printed | eq:hyperplane + `lem:chaindrop`(c_ρ≠0 ⇒ no 3-chain ⇐ `thm:Q30` ⇐ BR) + Thm C | **no, as printed** |
| `thm:Z` rerouted | (★) + chain drop for V(θ)≠0 (`lem:R`) + Thm C | **yes** |
| `thm:window` inductive proof | (★), chaindrop, `lem:A`(ii) (c_{ρ*}≠0 ⇐ Q30) | as printed no; rerouted with V: yes |
| `thm:window` from counts | thm:main(b) + lem:Npos | yes (was already) |

**Thm Z rerouted.** Replace `eq:hyperplane` by (★) at x = M/2:
s_{M/2+1}(κ) − s_{M/2}(κ) = Σ_{ρ⊊κ} V(θ_ρ) s_{M/2}(ρ). If V(θ_ρ) ≠ 0 then θ_ρ has no 3-chain by `lem:R` (combinatorics),
so the chain-drop bound r(ρ') ≥ r(κ') − 2 holds by the printed argument of `lem:chaindrop` (which only uses "θ has no
3-chain"). The induction then reads verbatim: every term vanishes on the shrunken window, so s_{(M+2)/2}(κ) = s_{M/2}(κ),
and telescoping from Thm C (s_{3/2}(κ) = 0 for ℓ(κ) ≥ 2) gives the zeros. Also: `thm:Z` is literally the vanishing
half of `thm:HZ`(i) (window 1 ≤ a ≤ ⌈r(κ')/2⌉ with M = 2a+1), so it is BR-free in any case; the reroute just keeps its
*separate* proof BR-free too.

Similarly in the inductive proof of `thm:window` on LG: the coefficients of (★) are V(θ_ρ) ≥ 0 directly, the
vanishing step needs only "V(θ)≠0 ⇒ no 3-chain" (`lem:R`), and `lem:A`(ii) needs V(T_{t*}) ≠ 0, which its printed proof
establishes via `lem:R` anyway (the appeal to Q30 only converts V≠0 to c≠0).

**What remains BR-dependent on LG:** `thm:Q30` (the closed form of the K-theoretic Pieri coefficients c_ρ — a statement
about K(X), not about χ), `prop:BRp1`, Examples `ex:pieri32` / `ex:dict31`, and the K-theoretic second proof of
`thm:main`. On OG, `lem:og-corners` and the K-theoretic proof of `thm:og` still use BR Lemma 4.5 (and `thm:og` is
classical anyway).

**(★) does not give Q30 back.** Combining (★) with `eq:hyperplane` gives Σ_{ρ⊊κ} [c_ρ + (-1)^{|θ|}V(θ)]·H_ρ ≡ 0. This
forces c_ρ = (-1)^{|θ|-1}V(θ) only when {H_ρ : ρ⊊κ} is linearly independent. For κ = (3,2,1) on LG(3,6) there are 7 such
ρ, all with |ρ| ≤ 5, i.e. 7 polynomials in a 6-dimensional space — dependent. So Q30 genuinely needs K-theory (BR).
(Smallest case κ=(3,1): five ρ⊊κ, degrees ≤ 3. Census in Verification.)

## Verification

`work-in-progress/scripts/prove258/star_poly.py m`: computes H_κ from type-C Demazure characters (`prove249d/dem.py`),
interpolates, and checks (★) **as a polynomial identity** in x, with V(θ) brute-forced from the definition (no BR, no Q30
code). Results: m = 2, 3, 4 (4, 8, 16 shapes): 0 failures.

**Negative control (deliverable 5).** For every pair ρ ⊊ κ, add 1 to the single entry V(θ_ρ) and test (★) at integers
j = 1, 2, …. Every perturbed entry is caught at an integer j ≥ 1: m=2 6/6 (first failure at j ≤ 2), m=3 27/27 (j ≤ 3),
m=4 110/110 (j ≤ 4). So the integer check is not vacuous: no entry sits only on terms that vanish at all tested integers.
(Logically this is automatic — s_j(ρ) > 0 for j > ⌈r(ρ)/2⌉ by thm:window — and the control confirms the code sees it.)

`prove258/indep.py`: κ with linearly dependent {H_ρ : ρ⊊κ}: m=2 none; m=3 (3,1), (3,2), (3,2,1); m=4 ten shapes
(3,1),(4,1),(3,2),(4,2),(4,3),(3,2,1),(4,2,1),(4,3,1),(4,3,2),(4,3,2,1). Smallest: κ=(3,1), five ρ⊊κ of degree ≤ 3
(dimension 4). So (★)+eq:hyperplane does not determine c_ρ already at κ=(3,1) — Q30 stays a BR theorem.

## Gaps

None in the proof of (★). The draft edits (chaindrop / lem:A / Thm Z / window rephrased with V(θ)) are re-wordings whose
mathematical content is `lem:R`, already proved in the draft.
