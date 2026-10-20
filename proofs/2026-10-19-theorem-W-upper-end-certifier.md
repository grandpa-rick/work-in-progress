# Theorem W, upper end: route R0, the twisted family, and a sum certifier (Day 247 PROVE, cycle 2)

Rick, 2026-10-19 (PROVE cycle 2/2). Continues `proofs/2026-10-19-theorem-W-sharp.md` §5–§6, Gap 1: the upper end
for λ or μ ∉ G. Notation as there and in `proofs/2026-10-18-q488-d-infinity.md` §1–§6. Scripts:
`work-in-progress/scripts/prove247b/`.

## 0. Verdict

| item | verdict | grade |
|---|---|---|
| **Lemma R0** (u→0 out-coefficient). ν^c top-clean ⇒ out-terms = u^{n−m−1}(K·ℓ^{[n−1]}_{ν̃} + O(u)). ν = [n] ⇒ out-terms = u^n(ℓ^{[n−1]}_{[n−1]}∘w^{−2} + O(u)). | Clean ν^c moves the obstruction from the λ side to the ν side. | proved (+ exact numerical check) |
| Two-route certifier (R0 or R∞, untwisted, single instance) | Certifies **all 5578 A-triples at N = 5** for every c ≥ 7. At N = 6 it leaves 1037 triples, all at doubly-dirty nodes. | computed (certificate) |
| **Lemma T** (twisted family Ψ_{κ,θ}: all u→0 / u→∞ leading options), label-lowering identity, the (1)-identities, dual-side and B/A4/B4 bookkeeping | All derived in §2–§3. | sketched → computed (negative-controlled) |
| **Sum certifier** (iterated leading coefficient = signed count) | **N = 4, 5, 6: every allowed A-triple is certified nonzero for ALL c > 2N−4** (408 / 5578 / 68 791). | computed (computer-assisted certificate) |
| ⇒ exact windows | **Z^A_5 = [0,6], Z^A_6 = [0,8]**. Lower inclusion: Prop Z (proved). c ≤ −1: W⁺ (proved). N = 7: 789 998 allowed, **325 uncertified** (§5 gap 2 family), so Z^A_7 = [0,10] is NOT certified. | computed |
| B, C upper ends | **B: c > 2n−2 is certified clean for all n ≤ 6.** C: certified for n ≤ 4; n = 5 leaves 7 nodes and n = 6 leaves 208 (leading-order cancellation, §5). | computed |
| Negative control | In fixed-c mode inside the window: **0 false certifications over every known zero** (A: N = 4, 5, all c ∈ [−1, 2N]; B, C: n ≤ 4, all c). | — |

The Day 247 cycle-1 obstruction "the in-coefficient at u→∞ is a block sum for consecutive λ" has an exact
mirror: the out-coefficient at u→0 is a block sum for consecutive ν^c. The two routes cover each other's
failures except on doubly-dirty nodes. Those are handled by letting nodes be formal sums of twisted instances and
computing the iterated leading coefficient exactly.

## 1. Lemma R0

**Lemma R0.** System [n], u = w_n, n ∈ ν, κ := [n]∖ν ⊆ [n−1], L ⊆ [n−1] (so n ∉ L). As u → 0:
(i) if κ = ∅: ℓ_{[n]}(e_L) = u^n( w_L^{−2} ℓ^{[n−1]}_{[n−1]}(e_L) + O(u) ) (constant exactly 1 up to the global τ_∅);
(ii) if κ ≠ ∅, m := max κ, and m−1 ∉ κ (or |κ| = 1): ℓ_ν(e_L) = u^{n−m−1}( K·ℓ^{[n−1]}_{ν̃}(e_L) + O(u) ) with
ν̃ = (ν∖n) ∪ {m} (the L5 shift) and K ≠ 0 independent of L.

*Proof.* (D): ℓ_ν(e_L) = [c(L)/c(ν)]·φ⁰_κ(L^c)/(F_L(L)φ⁰_{L^c}(L^c)), complements in [n], n ∈ L^c.
- c(ν) = c(ν∖n)·(−1)^n W u, with W = Π_{k<n} w_k.
- φ⁰_{L^c}(L^c) = e⁰_n(L^c)φ⁰_{L'}(L') with L' := [n−1]∖L (L1). The lowest term of e⁰_n(L^c) is
  (−1)^n W u^{−(n+1)} Π_{L'} w^{−2}.
- F_L(L) is u-free.
- (i) φ⁰_∅ = 1. The product is W^{−2} u^n Π_{L'} w² · c(L)/(c([n−1])F_L(L)φ⁰_{L'}(L')), and
  Π_{L'}w² = Π_{k<n}w_k²·w_L^{−2} = W²w_L^{−2}.
- (ii) Termwise, in φ⁰_κ(L^c) the variable u as head j has u-order −κ_j − 1 − (r'−j), and as a tail it has
  order −r'. Since κ_j − j is nondecreasing in j, and strictly increasing at the top iff m − 1 ∉ κ, head r'
  is the unique minimum. Tails lose since κ_{r'} ≥ r'. In that term the u-factors tend to:
  −u^{−2} (from 1 − u^{−2}); (−1)^{m−1}W_{<m}u^{1−m}; x_i^{−2} for the earlier heads; w_t^{−2} for the tails.
  So φ⁰_κ(L^c) = (−1)^m W_{<m} u^{−m−1} Π_{L'} w^{−2}(φ⁰_{κ∖m}(L') + O(u)).
- The Π w^{−2} cancels against e⁰. Since [n−1]∖(κ∖m) = ν̃, and c(ν̃)/c(ν∖n) is L-free, the claim follows. ∎

Check (`t0.py`, exact mod p series in u, n = 4, 5): for every clean ν^c, lead(ℓ_ν(e_L))/ℓ_{ν̃}(e_L) is one
constant over all L, at the predicted order. For dirty ν^c the order is still n−m−1, but the ratio takes 4–15
values. A single ℓ_{ν'}∘w^e does not fit either (`t1.py`, e ∈ [−4,4]).

**Route R0 (mixed case, large c).** In-terms at u→0 have order c+1−δ (A, C) or c−δ (B), and their coefficient
is always clean (L2). By R0, out-terms have order n−m−1 ≤ n−2, or n when ν = [n]. So for c > U^X_n they dominate,
with coefficient X_{[n−1]}(λ, μ, ν̃; c) or X_{[n−1]}(λ, μ, [n−1]; c−2). λ and μ are untouched. Supports match: by L5,
E ∩ {∌ n} = E(ν̃); for ν = [n] this is immediate. The dominance needs n ≥ 5 for A when ν = [n], δ = 2, and
n ≥ 3 otherwise.

So R∞ (cycle 1, Theorem U_G) fails exactly for consecutive tops of λ or μ, and R0 fails exactly for consecutive
tops of ν^c. `cert.py` (both routes, recursive) certifies every A-triple at N ≤ 5; either route alone leaves
385–473. **Z^A_5 = [0,6] already follows from R0 + R∞ + the exact steps.** At N = 6 it leaves 1037 triples.
Every failure is doubly dirty, e.g. C_4(∅,(1,2),(1,4)): μ = (1,2), κ^c = (2,3).

## 2. The twisted family

For labels κ = (κ_1<…<κ_r) and an integer twist θ ∈ ℤ^r put
  Ψ_{κ,θ}(M) := Σ_h T⁰_κ(h) Π_j x_j^{θ_j}   (h injective [r]→M, x_j = w_{h_j}).
Then F_κ = w_M^{2r}·Ψ_{κ,(κ_j−2j)} and φ⁰_κ = Ψ_{κ,0}. Ψ_{κ,θ}(M) = 0 unless M ≥ κ. Reason: T⁰ contains
(1 − w_{h_j}/x_j) = 0 whenever h_j < κ_j, and an injection with h_j ≥ κ_j exists iff M ≥ κ (Hall's theorem,
top-aligned).

**Lemma T (leading options).** Let n ∈ M, κ ⊆ [n−1], u = w_n. Every T⁰-term falls in exactly one option:
"u is head j" or "u is a tail". Termwise:

| option | u→0 order | coefficient (times Ψ on M⁻) | u→∞ degree | coefficient |
|---|---|---|---|---|
| head j | θ_j − κ_j − 1 − (r−j) | (−1)^{κ_j+r−j} W_{<κ_j} Π_{M⁻}w^{−2} · Ψ_{κ∖κ_j, θ^{(j)}} | θ_j − (r−j) | (−1)^{r−j} Ψ_{κ∖κ_j, θ^{(j)}} |
| tail | −r | (−1)^r Ψ_{κ, θ−1} | −r | (−1)^r Ψ_{κ, θ+1} |

Here θ^{(j)}_i = θ_i for i < j and θ_{i+1} + 1 for i ≥ j. Forced case (κ_r = n): Ψ_{κ,θ}(M) =
e⁰_n(M)u^{θ_r}Ψ_{κ⁻,θ⁻}(M⁻) exactly.
*Proof.* Factor by factor. Head factors: 1 − u^{−2} → −u^{−2} / 1. Π_{k<κ_j}(1 − w_k/u) →
(−1)^{κ_j−1}W_{<κ_j}u^{1−κ_j} / 1. Earlier pairs: (1 − 1/(x_iu))/(1 − x_i/u) → x_i^{−2} / 1. Later pairs:
(1 − 1/(ux_i))/(1 − u/x_i) → −1/(ux_i) / −x_i/u. Tails: (1 − 1/(w_tu))/(1 − w_t/u) → w_t^{−2} / 1. Tail u
against head i: → −1/(ux_i) / −x_i/u. ∎
(For untwisted F this recovers L2 (u→0: head r is the unique minimum) and L3′ (u→∞: ties along the top
consecutive block).)

**Exact identities used to resolve ties.**
- *(1)-identities* (residue sum, Day 247 §5): Ψ_{(1),(0)}(M) = φ⁰_{(1)}(M) = 1 − Π_M w^{−2}. Hence for κ = (1):
  at u→∞ with θ = (−1), head + tail = Π_{M⁻}w^{−2} (degree −1); at u→0 with θ = (1), head + tail = −1 (order −1).
- *Label lowering*: if κ′ ≥ κ positionwise (both increasing), then
  Ψ_{κ′,θ} = Σ_{T_i ⊆ [κ_i, κ′_i−1]} Π_i (−1)^{|T_i|} w_{T_i} · Ψ_{κ, θ − Σ_i|T_i|e_i}.
  Proof: termwise, Π_{k<κ′_i}(1 − w_k/x_i) = Π_{k<κ_i}(·)·Π_{k=κ_i}^{κ′_i−1}(1 − w_k/x_i). Vanishing for h_i < κ_i is
  automatic in both expressions. This turns a block sum into like terms with the same labels, so that systematic
  cancellations (e.g. 1 − (1 − y) = y in a block of length 2) become explicit.

**Dual side.** Replace φ⁰_{ν^c} in (D) by Ψ_{ν^c,ψ} (twisted dual ℓ_{ν,ψ}). Out-terms (n ∉ M, n ∈ M^c) use
Lemma T on Ψ_{ν^c,ψ}(M^c). In-terms are untouched. ℓ-factors, with W = W_{<n}:
- in, u→0: (−1)^nW^{−1}u, weight 0;
- in, u→∞: u^{−n}, weight −2 on M;
- out, u→0: W^{−2}u^n Π_{M^{c'}}w^{ws+2} × option;
- out, u→∞: (−1)^nW^{−1}u^{−1} Π_{M^{c'}}w^{ws} × option.
Here ws ∈ {0, −2} is the option's own Π w factor, and Π_{M^{c'}}w^k = Π_{[n−1]}w^k · w_M^{−k}. In every case
there is also the correction c(ν_new)/c(ν⁻).
- **Type B:** the weight lives on S = M∖a, so every weight w_M^k coming from ℓ also carries a scalar w_a^k. That
  scalar was the one bug the negative control caught (§4).
- The B in-factor 1/g_a gives −x^{−1}u^{−1} (u→0) and −x u^{−1} (u→∞).
- A4 and B4 are exact with twists. On H, w_b ↦ w_a^{−1} in scalars.
- n ∉ ν is exact: X_{[n]} = u^{ψ_last} X_{[n−1]}.

## 3. The sum certifier and why it is sound

A node is a finite formal sum Σ_i σ_i·s_i·X_i. Here σ_i = ±1, s_i is a monomial in the not-yet-peeled
variables, and X_i is a twisted instance (type, λ, θ, μ, θ′, ν, ψ, weight offset d_i) at level n.
- *Peel.* Peel u = w_n. Pick a limit (u→0 or u→∞; both are tried). List every option of every summand, with
  order (const + [in]·(c + d_i) + exponent of w_n in s_i).
- *Extremal set.* Keep the extremal ones. c is symbolic with c ≥ c_min: an extremal set is used only if it is the
  extremal set for every c in the range. Otherwise the range is split, with finitely many single c values plus a
  tail range.
- *Combine.* Combine like terms after label-lowering, drop zero coefficients, and drop summands with empty index
  set (identically 0).
- *Exact steps.* A4, B4 and n ∉ ν are exact. At B4 the summands are grouped by their power of w_a. The
  C-children no longer involve w_a, so it suffices that one group certifies.
- *Base.* At n = 0, compute the signed count Σσ_i.

**Soundness.** Each option's coefficient is the exact coefficient of its termwise-extremal order (Lemma T and §2
are exact leading-term statements, and non-extremal options are strictly sub-leading). So the child sum is exactly
the leading coefficient of the parent in that limit, *provided the child is ≢ 0*. Inductively, a nonzero signed
count at n = 0 shows that every intermediate child is ≢ 0. Hence the root is ≢ 0, for every c in the range. A
cancellation (child ≡ 0 or certificate failure) proves nothing; the certifier is incomplete but cannot certify a
zero. Labelled supports: Ψ = 0 off M ≥ κ, and M^c ≥ κ^c ⟺ M ≤ κ (checked n ≤ 8; it is the counting criterion
for top-aligned dominance).

## 4. Verification

**Negative control** (`negctl.py`, `negBC.py`). In fixed-c mode the certifier was run on every allowed instance
with c ∈ [−1, 2N]. The independent mp census (`prove247/mp.py`, two random points mod 2^61−1) gives the zero sets.
- *First run:* 2–280 false certifications per c, all on B-levels (pairs with b > a+1, and direct B instances).
  Diagnosis: B's weight sits on S, not M, so every ℓ-weight w_M^k needs an extra scalar w_a^k. That was the only
  change.
- *After the fix:* **0 false certifications** over all 2 432 known zeros: A (N = 4, 5, c ∈ [−1, 2N]), B and C
  (n = 2, 3, 4, c ∈ [−1, 2n]). Meanwhile 92–100% of the nonzero instances inside the window are certified.

**Positive results** (`run_s.py`, `run_bc.py`; range c ≥ U^X_n + 1):

| | allowed | uncertified |
|---|---|---|
| A, N = 4 | 408 | 0 |
| A, N = 5 | 5 578 | 0 |
| A, N = 6 | 68 791 | 0 (52 s) |
| A, N = 7 | 789 998 | 325 (18 min; e.g. (∅,(1,2),(1,2,6,7)) H12, which reaches the C-cancellation family) |
| B, n ≤ 4 / 5 / 6 | 1 194 / 9 594 / 85 335 | 0 / 0 / 0 |
| C, n ≤ 4 / 5 / 6 | 865 / 4 872 / 34 008 | 0 / 7 / 208 |

Policy probe (`pol.py`, N = 5): u→0 only leaves 30 uncertified, and u→∞ only leaves 924. The 30 are long dual
blocks (ν^c = (1,2,3), (2,3,4)) and ν = [5] at the c-edge.
Shape: 33 505 of 33 559 successful certify calls carry a single summand. Multi-summand nodes are the dual block
sums, e.g. ν^c = (2,3) with twists (0,0), (1,0), (1,1) and signs +, −, + (`instr2.py`).

## 5. Gaps (precise)

1. **All N is open.** The certificate is per N. A uniform proof needs a closed form for the dual block sums at
   u→0 (or for the F block sums at u→∞), plus a proof that the iterated signed count is nonzero. Data: u→0 alone
   nearly suffices.
2. **C at n = 5, 6** (7 + 208 nodes, e.g. C(∅,(1,2),(4,5)), ν^c = (1,2,3)). After the u→0 dual block, the 3-summand
   child's leading coefficients cancel identically at the next level, in both limits. So the true leading term is
   sub-leading, and the certifier would need second-order terms or a further residue identity. Census: C_5 is
   clean at c = 9 (cycle 1). They are not reached from A at N ≤ 6. At N = 7 the 325 uncertified A-triples all have ν ∈ {(4,5,6,7), (1,5,6,7), (1,2,6,7)} (229/80/16), i.e. ν^c is a consecutive block of length 3, the same feature as these C-nodes (attribution to the identical cancellation is plausible, not traced). All 325 are nonzero at c = 11, 12, 13 (`chk7.py`: nonzero exact mod-p evaluation, which is a certificate at those c), so no evidence of an upper-end zero. **Next target: the length-3 dual block identity.**
3. **Code trust.** The certificate is only as good as `scert.py`'s transcription of Lemma T and §2. Mitigation: the
   negative control above, which has teeth (it caught a real bug), and agreement with the independent census at
   every tested (instance, c). The lemmas were derived today and are not yet re-derived cold. Grade: computed.
4. Everything downstream of (D) inherits the Motegi citation (Day 246 gap 1).

> **Day 248 update:** gaps 1–2 superseded by `proofs/2026-10-20-dual-block-identity.md`. Gap 2 was a formal zero (alternant relations among twisted Ψ). Certifier v2 clears C n≤6 and A N=7 (0/789 998), so Z^A_7 = [0,10]. Lemma T was re-derived cold, exactly (Lemma H).
