# Q488 converse at first order in q — what's true, what's false, what's proved

Rick, Day 243 PROVE (session dated 2026-10-17; file name per PROVE.md).
Scripts: `work-in-progress/scripts/day243/q488c/`.
Source: Motegi, arXiv:2610.15344, e-print source `typebcweightfunctions_arxivversionfinal-3.tex`. I deep-read
L474–487 (Thm 3.3, the symmetrisation formula), L1057–1062 (E^C), L603–620 (Thm 2.6(i), triangularity and the
diagonal) and L2641–2848 (stability) first-hand today. Clio's Q488 tex (0e59afb) supplies notation.

## 0. Verdict in one screen

PROVE asked for (a)–(d). **Two of the three main targets are FALSE as stated, and both fail first at N = 4.**

| item | statement | status |
|---|---|---|
| (a) | c ∈ C(w)[[q]] | **PROVED**, all N (§2) |
| (b) | jump ≥ 2 ⇒ c₀ = 0 | **FALSE** at N = 4: c₀^{1234}_{12,34} = 1/(w₃²w₄²). 17 jump-≥2 triples at N=4 have c₀ ≠ 0. True for N ≤ 3 (exhaustive) |
| (c) | jump ≥ 2 ⇒ c₁ has a long-root pole | **FALSE** at N = 4: (λ,μ,ν) = (∅,{3},{1,2,4}) has c₁ Laurent on every H_ab. Its full-q residue on H₁₂ is ≠ 0 but has q-valuation **2**. Holds at N ≤ 3 (13/13), and for 126/127 jump-≥2 triples at N=4 |
| (c-top) | if ν minus its two largest elements is ≥ λ, μ, then c₁^ν has a simple pole on H_{top two} | **PROVED**, all N, all jumps (§6, Thm 9) |
| (c-agg) | for each (λ,μ): some c^ν_{λμ} (some ν) has a pole on H_ab ⇔ ∃ M'' ⊆ [N]∖{a,b} with M'' ≥ λ, μ | **PROVED**, all N, full q (§4, Thm 6) |
| (d) | jump ≤ 1 ⇒ c₁ Laurent | **PROVED**, all N (§6bis, Thm 10). Also c₀ is Laurent for every triple. Computed: N ≤ 3 full; N = 4, 0 long-root poles in 599 jump-≤1 triples on all six pairs |

So **the order-q¹ route cannot settle Clio's converse (Gap 2) in general**. It settles the converse exactly
on the "top-two" family, which includes every (∅,∅,ν) with |ν| ≥ 2, the N=2 witness, and 6/13 of the N=3 jump-≥2
triples; with the single-defect theorem (§6, Thm 9') all 13 at N=3 and 111/127 at N=4. The missing triples are those whose created pole sits in a *high block* (§7). There the q¹ coefficient of
the residue can vanish identically while the full-q residue survives.

Engine of everything: **Lemma QC (§3)**, a q-twisted collision on the long-root divisors, proved from Motegi's
formula:
> a<b ∈ M, on H_ab = {w_a w_b = 1}:  Ã_{M,I} = g(q,w_a) · Ã_{M∖{a,b},I} for every I,
> g(q,x) = (1−q)² / ((1−qx²)(1−qx⁻²)).

At q=0 we have g=1: that's GKM, and the q=0 basis is a K_T(LG(N,2N))-type Schubert basis. At first order,
g' = (x − x⁻¹)² =: γ. Clio's long-root poles are exactly the failure of the twisted GKM condition to be untwisted.

## 1. Setup and the facts I take from Motegi

φ_I(M) := Ã_{M,I} := W̃C_{N,n,|I|}(w_M|w|I), for any n ≥ |M|, |I| (n-independent by stability, Motegi L2823–2848).
The c^ν_{λμ} are defined by φ_λφ_μ = Σ_ν c^ν_{λμ}φ_ν (pointwise on the 2^N points M).
M ≥ I is Motegi's order (|M| ≥ |I|, and the top |I| elements of M dominate I entrywise).
Notation: M'' := M∖{a,b}; h := 1 − 1/(w_a w_b) is the local parameter of H = H_ab; O_H is the local ring.

Used from Motegi (all deep-read today):
- (M1) the symmetrisation formula (Thm 3.3, L474–487). WC = F at β=1, γ=0. F is a Laurent polynomial in u, w and a
  polynomial in q (a symmetrised rational function with at most simple poles on u_j = u_k, which cancel).
- (M2) E^C (L1057–1062); Ã = WC/E^C.
- (M3) Triangularity: Ã_{M,I} = 0 unless M ≥ I (Thm 2.6(i)).
- (M4) The diagonal Ã_{I,I} is Clio's eq:diag: a product over all positive roots of C_N. Long-root factors
  (1 − w_{i}^{-1}w_{i'}^{-1}) for i<i' in I each occur once. (Q)-factors (1 − q·monomial). No other factor vanishes
  identically on any H_ab for generic q. At q=0 the diagonal is the product of the L, S, A factors, and it is ≠ 0.
- (M5) Stability (L2641–2647, L2823–2848), so we may take n = |M| ("unpadded", all u-values distinct).

## 2. (a): c ∈ C(w)[[q]]  — PROVED

R := C[q, w^±]_{(q)} is a DVR with residue field C(w), and R ⊂ C(w)[[q]].
By (M1) WC ∈ C[q][w^±] after specialisation. Every factor of E^C is (1−q)^{±1} or (1 − q·monomial)^{±1}, a unit in R.
So Ã_{M,I} ∈ R. By (M3), det(Ã) = ± Π_I Ã_{I,I}. Each Ã_{I,I} mod q is the nonzero product of its L, S, A factors
(M4), so it is a unit of R, and hence Ã^{-1} has entries in R. Then c^•_{λμ} = Ã^{-1}(φ_λφ_μ) ∈ R^{2^N}. ∎

Consequence. c₀ := c|_{q=0} is the structure-constant tensor of the basis φ⁰ := Ã|_{q=0}, and
c₁ = ∂_q c|_{q=0}. Write δ_I := ∂_qφ_I|_{q=0}. Differentiating φ_λφ_μ = Σ c^ν φ_ν gives

  (1.1)  Σ_ν c₁^ν φ⁰_ν = Φ := φ⁰_λ δ_μ + δ_λ φ⁰_μ − Σ_κ c₀^κ δ_κ .

**Laurent at every order is necessary.** If c is Laurent in w over C(q), then since c ∈ R its coefficients are
regular at q=0 (a negative q-power would leave a nonzero Laurent polynomial in w as its coefficient). So every c_k
is Laurent. Hence: **a pole of c₁ certifies that c is not Laurent.**

## 2bis. The explicit q = 0 formula — PROVED

Put q=0 in (M1) with n=|M| (M5), u = w_M. Then E^C|_{q=0} = 1. Per head h (label i): (1−h⁻²)h^{N+1}Π_{k<i}(1 − w_k/h).
Head pair: h_jh_k(1 − 1/(h_jh_k))/(1 − h_j/h_k). Tail t: t^{N+1}. Tail–head: th(1 − 1/(th))/(1 − t/h).
Tail pair: tt'. Every value collects exactly u^{N+n}, which cancels the prefactor. The (n−r)! tail orderings
cancel the 1/(n−r)!. Hence, summing over injective h: [r] → M (h_j is the head carrying label κ_j):

  (2.1)  φ⁰_κ(M) = Σ_h Π_j (1−h_j⁻²) Π_{k<κ_j}(1 − w_k/h_j) · Π_{j<k} (1 − 1/(h_jh_k))/(1 − h_j/h_k)
                    · Π_{t ∈ M∖h} Π_j (1 − 1/(t h_j))/(1 − t/h_j).

(Checked against the engine on all 64 (M,I) pairs at N=3: `testq0.py`.) Immediate consequences:
- **φ⁰_∅ ≡ 1** (r = 0: the empty product). Hence c₀^ν_{∅∅} = δ_{ν∅}.
- **N-stability at q=0.** (2.1) only involves w_k for k < max κ, together with the values w_M. So for M, κ ⊆ [N'],
  φ⁰_κ(M) does not depend on N ≥ N'. In particular c₀^ν_{λμ}, for λ, μ, ν ⊆ [N'], does not depend on N.
- **Vanishing terms.** If h_j = w_m with m < κ_j, the factor k = m is (1 − 1) = 0, identically in w (and at all q:
  the factor (u_j − w_k), k < i_j, carries no q).

**Lemma L (limit lemma).** Let a, b ∉ M, and let M ≥ κ. Restrict to H_ab: w_a = x, w_b = 1/x, with x fixed and
generic. Take w_c = t^{−e_c} for c ∉ {a,b}, with 0 < e_1 < e_2 < … increasing in c, and let t → 0⁺. Then
φ⁰_κ(M) → 1. Hence φ⁰_κ(M)|_H ≢ 0. (Also without H: φ⁰_κ(M) → [M ≥ κ]. This is the q=0 converse of (M3).)

*Proof.* Every surviving term has h_j = w_m with m ≥ κ_j. In such a term every factor (1 − X) whose X → 0 tends to 1.
That covers (1 − h⁻²), (1 − 1/(h h')), (1 − 1/(t h)), and the factors (1 − w_k/h_j) with k < κ_j ≤ m. (For
k ∈ {a,b} these are 1 − x^{±1}/h_j → 1, since h_j → ∞.) The denominators (1 − h_j/h_k) and (1 − t/h_j) tend to 1
when the numerator variable has the smaller index, and to ∞ otherwise. So a term tends to 1 if its heads increase
and every tail lies below every head, and to 0 otherwise. Exactly one such h exists: the top r elements of M, in
increasing order. It is a surviving term iff m_{s−r+j} ≥ κ_j for all j, i.e. iff M ≥ κ. ∎

## 3. Lemma QC (q-twisted long-root collision), all q — PROVED

**Lemma QC.** Let a < b with a, b ∈ M. On H_ab, for every I: Ã_{M,I} = g(q, w_a) · Ã_{M'',I}, where
g = (1−q)²/((1−qw_a²)(1−qw_a^{-2})).

*Proof.* Use (M1) unpadded: n = s = |M| for Ã_{M,I}, and n = s−2 for Ã_{M'',I} (M5). Write x := w_a, so w_b = 1/x.
1. *Special values in the head die.* If both x and 1/x are heads, the head-pair factor (u_ju_k − 1) = 0. If one
   is a head and the other a tail, the tail–head factor (u_ju_k − 1) = 0. So only permutations that put both special
   values in the tail survive. Tail positions are symmetric, so after cancelling (n−r)! the sum runs over ordered
   head selections h from the remaining s−2 values. These are exactly the head selections of the (s−2)-formula.
2. *Compare term by term* (same h, remaining tails T). The s-term is the (s−2)-term times the extra factors:
   - the prefactor (x·x⁻¹)^{−(N+s)} · Π_{v∈rest} v^{−2} (the change N+s−2 → N+s);
   - (βx − γ)(β/x − γ) = 1;
   - Π_k (x − q w_k)(1/x − q w_k);
   - for each head h: (qx − h)(xh − 1)/(x − h) · (q/x − h)(h/x − 1)/(1/x − h) = (qx − h)(q/x − h)
     = h²(1 − qx/h)(1 − q/(xh));
   - for each other tail t: (xt − q)(t/x − q) = t²(1 − q/(xt))(1 − qx/t);
   - the pair itself: (x·x⁻¹ − q) = (1 − q).
   The x-dependence is the same for heads and tails. So the bracket equals
   (1−q) Π_k(x−qw_k)(1/x−qw_k) Π_{v∈rest} v²(1 − qx/v)(1 − q/(xv)), independent of h, and the v² cancels the
   prefactor. Hence F_s(rest, x, 1/x) = (1−q) Π_k(x−qw_k)(x⁻¹−qw_k) Π_v (1−qx/v)(1−q/(xv)) · F_{s−2}(rest).
   (If s−2 < r, no head selection survives and both sides are 0.)
3. *Normalisation.* By (M2), E_s/E_{s−2} = (1−qx²)(1−qx⁻²) Π_v(1−qx/v)(1−q/(xv)) Π_k(1−qxw_k)(1−qw_k/x) / (1−q).
   (The (1−q)^{−2} offset against the diagonal pair factors (1−q)², and the Π_v(1 − q x v)(1 − q v/x) cancel between
   numerator and denominator.) Since (x − qw)(x⁻¹ − qw) = (1 − qxw)(1 − qw/x), the ratio is (1−q)²/((1−qx²)(1−qx⁻²)). ∎

Computed on all (M, pair, I) at N=2, 3 (`ratio.py`).

**Corollary 4.** (i) At q=0: φ⁰_I(M) = φ⁰_I(M'') on H_ab (long-root GKM). (ii) At first order, since g(0) = 1 and
∂_q g(0) = x² + x⁻² − 2 = γ := (w_a − w_a⁻¹)²: δ_I(M) − δ_I(M'') = γ φ⁰_I(M'') on H_ab.

## 4. Lattice lemma and the aggregate converse, full q — PROVED

Fix a<b. For g ∈ k(H)ˣ let L_g := {f ∈ O_H^{2^N} : f(M) − g f(M'') ∈ hO_H for all M ⊇ {a,b}}.

**Prop 5.** span_{O_H}{φ_I} = L_g (generic q, with g from QC), and span_{O_H}{φ⁰_I} = L_1. Hence, for any
f ∈ O_H^{2^N}: every coordinate of f in the basis φ is regular on H iff f ∈ L_g. Poles are at most simple, and the
polar part Res_H ℓ_ν(f) depends only on the defects (f(M) − g f(M''))|_H, linearly.

*Proof.* (i) ⊆: the entries lie in O_H (they are Laurent in w over C(q), or over C at q=0), and they satisfy the
congruence by QC. (ii) The colength of L_g is #{M ⊇ ab} = 2^{N−2}, since the defect map onto k(H)^{2^{N−2}} is
surjective (take f = c·e_M) and has kernel L_g. (iii) The colength of the span is v_H(det Ã) = Σ_I v_H(Ã_{I,I}),
and that is 2^{N−2} by (M4): exactly the I ⊇ {a,b}, each with multiplicity one. Equal colengths and an inclusion
give equality. The rest follows because h·O^{2^N} ⊆ L_g. ∎

**Theorem 6 (aggregate converse, all N, generic q).** For fixed λ, μ and a<b, the following are equivalent:
(i) every c^ν_{λμ} (ν ranging over all subsets) is regular on H_ab;
(ii) φ_λφ_μ(M'')|_H = 0 for every M'' ⊆ [N]∖{a,b};
(iii) no M'' ⊆ [N]∖{a,b} satisfies both M'' ≥ λ and M'' ≥ μ.

*Proof.* f := φ_λφ_μ satisfies f(M) = g² f(M'') on H by QC. So f ∈ L_g iff (g² − g) f(M'')|_H = 0. Since g ≠ 0, 1
on H, Prop 5 gives (i)⇔(ii). (iii)⇒(ii) is (M3). (ii)⇒(iii): if M'' ≥ λ, μ, then by Lemma L
φ⁰_λ(M'')φ⁰_μ(M'')|_H ≠ 0, and so its deformation φ_λφ_μ(M'')|_H is ≠ 0 at generic q. ∎

So, for every (λ, μ), *some* structure constant c^•_{λμ} acquires a long-root pole as soon as N ≥ max(|λ|,|μ|) + 2,
which is exactly when a jump-≥2 ν exists. The per-ν converse is the hard part.

## 5. The first-order residue formula — PROVED

**Theorem 7.** For every N, every triple, and every a<b:
  Res_{H_ab} c₁^ν_{λμ} = γ · Res_{H_ab} ℓ⁰_ν(1_{⊇ab} · P),  where P := φ⁰_λφ⁰_μ.
Here Res f := (h f)|_H, ℓ⁰_ν is the ν-coordinate in the basis φ⁰, and 1_{⊇ab} is the indicator of {M ⊇ {a,b}}.

*Proof.* By (1.1), c₁ = ℓ⁰(Φ). Φ ∈ O_H^{2^N}: δ is Laurent in w, and c₀ is regular on H by Prop 5 at q=0, because
P ∈ L_1 (Corollary 4(i)). By Prop 5 the residue depends only on the defects of Φ. By Corollary 4:
Φ(M) − Φ(M'') = γP(M'') + γP(M'') − γΣ_κ c₀^κ φ⁰_κ(M'') = γP(M'') on H. The function 1_{⊇ab}P has defect
P(M) − 0 = P(M'') on H at every M ⊇ ab. ∎

Checked against direct c₁ on all 167 (triple, pair ⊆ ν) at N=3 (`testresq0.py`).
Note: ℓ⁰_ν only sees values on ↓ν, because the inverse of a matrix triangular for the order is again triangular.

## 6. The top-two theorem — PROVED

**Prop 8 (forced heads).** Let a<b, let κ = κ''∪{a,b} and M = M''∪{a,b} with κ'', M'' ⊆ [a−1]. Then
φ⁰_κ(M) = e(M'') φ⁰_{κ''}(M''), where φ⁰_{κ''}(M'') is formula (2.1) for subsets of [a−1], and
  e(M'') = (1−w_a⁻²)(1−w_b⁻²)(1−1/(w_aw_b)) Π_{c∈M''}(1−1/(w_cw_a))(1−1/(w_cw_b))
           · Π_{c∈[a−1]∖M''}(1−w_c/w_a)(1−w_c/w_b) · Π_{a<c<b}(1−w_c/w_b).

*Proof.* In (2.1), label b needs a head w_m with m ≥ b, otherwise the term vanishes identically. The only such
element of M is b. Label a then needs m ≥ a with m ≠ b, which forces m = a. In the cross factors, every c ∈ M''
(head or tail) meets the heads w_a, w_b through the same label-free expression
(1−1/(w_cw_a))(1−1/(w_cw_b))/((1−w_c/w_a)(1−w_c/w_b)). Its denominators cancel the factors k = c of
Π_{k<a}, Π_{k<b}, and (1 − w_a/w_b) cancels the pair denominator. What remains is exactly the [a−1]-formula
for κ'' on M''. ∎

On H: e = h·e', and the identity (1−1/(w_cw_a))(1−w_a/w_c) = w_c⁻²(1−w_c/w_a)(1−w_aw_c) gives
  (6.1)  e'(M'')|_H = −γ·K·Π_{c∈M''} w_c^{−2},  K := Π_{c<a}(1−w_c/w_a)(1−w_aw_c) · Π_{a<c<b}(1−w_aw_c) ≠ 0.

**Theorem 9 (top-two).** Let a<b be the two largest elements of ν, put ν'' := ν∖{a,b}, and suppose ν'' ≥ λ and
ν'' ≥ μ. Then
  Res_{H_ab} c₁^ν_{λμ} = −K^{−1} · Ξ,  where Ξ := [φ⁰_{ν''}]( M'' ↦ L(M'') φ⁰_λ(M'')φ⁰_μ(M'') ) ≠ 0,
with L(M'') := Π_{c∈M''} w_c², and the coordinate taken in the q=0 system on subsets of [a−1]. In particular
c₁^ν has a simple pole on H_ab, so **c^ν_{λμ} is not Laurent** (§2).

*Proof.*
1. U := {M ≤ ν : a, b ∈ M} = {M''∪ab : M'' ⊆ [a−1], M'' ≤ ν''}, and U is an up-set of ↓ν. Reason: M ≤ ν forces
   the top two elements of M to be ≤ (a, b), and a, b ∈ M forces them to equal (a, b). For M ≤ M̃ ≤ ν the top two
   of M̃ are squeezed between those of M and those of ν.
2. On ↓ν we have 1_{⊇ab} = 1_U. Functions on ↓ν supported in U are spanned by {φ⁰_κ : κ ∈ U}, which is
   triangular on U. So, by Theorem 7, Res c₁^ν = γ Res ℓ^U_ν(P|_U), with ℓ^U the coordinates in {φ⁰_κ|_U}.
3. By Prop 8, ℓ^U_ν(F) = ℓ^{[a−1]}_{ν''}(F(·∪ab)/e). The [a−1]-basis involves only w_1, …, w_{a−1}, and its
   diagonal factors are nonzero, so it is regular and invertible on H. Hence
   Res ℓ^U_ν(P|_U) = ℓ^{[a−1]}_{ν''}(P(·∪ab)/e')|_H. Since P(M''∪ab)|_H = P(M'')|_H (Corollary 4), (6.1) gives
   Res c₁^ν = γ·(−1/(γK))·Ξ.
4. Ξ ≠ 0. Take w_c = t^{−e_c} with e_1 < … < e_{a−1}. By Lemma L the [a−1]-matrix tends to the zeta matrix, so its
   inverse stays bounded and [ν'', ν'']-entry → 1. Also P(M'') → [M'' ≥ λ][M'' ≥ μ], and for M'' < ν'' we have
   L(M'')/L(ν'') → 0, because Σ_{c∈M''} e_c < Σ_{c∈ν''} e_c under top-aligned domination. So Ξ/L(ν'') → P-limit
   at ν'' = 1. ∎

Checked: the formula −K⁻¹Ξ matches direct c₁ on all 105 triples at N = 2, 3 (`toptwo.py`), zero cases included.

Coverage: every (∅, ∅, ν) with |ν| ≥ 2; the N=2 witness (residue −1, computed by hand); 6 of 13 jump-≥2 triples
at N=3. The other seven (λ or μ ∈ {{2},{3}}, ν = {1,2,3}) have residue 0 on the top pair (2,3), as the formula
predicts. They fall to the next theorem, at pair (1,2).

**Theorem 9' (single defect).** Suppose some a<b in ν has ν'' := ν∖{a,b} ≥ λ, μ, and no M < ν has a, b ∈ M and
M'' ≥ λ, μ. Then Res_{H_ab} c₁^ν = γ·ρ_{νν}·P(ν'')|_H ≠ 0, where ρ_{νν} = Res_H(1/φ⁰_ν(ν)).
*Proof.* By Prop 5 and Theorem 7, Res = γ Σ_{M⊇ab} ρ_{νM}P(M'')|_H with ρ_{νM} = Res ℓ⁰_ν(e_M). That vanishes unless
M ≤ ν, and P(M'') = 0 unless M'' ≥ λ, μ (M3). Only M = ν survives. ρ_{νν} ≠ 0 because φ⁰_ν(ν) has a simple zero on
H (M4), and P(ν'')|_H ≠ 0 by Lemma L, since ν'' ∌ a, b. ∎

Coverage of Theorems 9 ∪ 9' (`coverage.py`, purely combinatorial): N=2: 1/1; **N=3: 13/13**; N=4: 111/127;
N=5: 792/1136; N=6: 5117/9748. The 16 uncovered N=4 triples are (λ,μ) ⊆ {∅,{1},{2},{3},{4}}, ν ∈ {134, 124, 1234}.
15 of them have a c₁ pole numerically; the 16th is the counterexample of §7.

**Full-q version.** The forced-head argument uses only the factor Π_{k<i}(u − w_k), which carries no q. Together
with QC this gives Res_{H_ab} c^ν = (g² − g)·[ψ_{ν''}](φ_λφ_μ(M'')/e'^{(q)}(M''))|_H at full q, where ψ is the
reduced system. This is sketched, not written out. Its q → 0 limit is Theorem 9.

## 6bis. (d) and Laurentness of c₀ — PROVED, all N

**Lemma C (collapse; Clio's Lemma 3, re-derived).** For all q, with τ as in Clio:
(i) on {w_a = 1}, for a ∈ M: Ã_{M,I} = Ã_{M∖a,I};
(ii) on {w_a = w_b}, a<b, for b ∈ M, a ∉ M: Ã_{M,I} = Ã_{(M∖b)∪a,I};
(iii) Ã(−w) = Ã(w), and hence (i) also holds on {w_a = −1}.
*Proof.* (i), (ii): pad to a common n. The two specialisation tuples then have the same multiset, and F and E^C are
symmetric in u (Motegi L455–458, L1057). (iii): under (u,w) ↦ (−u,−w), the factors of F change sign with total
exponent r(N−1) + (n−r)N + (n−r) + n(N+n) ≡ n(n+1) ≡ 0 (mod 2). The heads give (N−i)+(i−1), the tails N, (βu−γ) = u
gives 1, and the pair factors are invariant. Every factor of E^C is invariant. Then apply (i) to −w. ∎

**Lemma C' (collision regularity at q=0).** Let D be one of {w_a = ±1}, {w_a = w_b}, with collapse map τ on 𝒟.
If f is a vector of Laurent polynomials with f(M)|_D = f(τM)|_D for all M ∈ 𝒟, then ℓ⁰(f) is regular on D.
*Proof.* This is Clio's thm:collision induction, verbatim, with f in place of φ_λφ_μ. It needs: τM < M (Clio);
φ⁰ has the collapse property (Lemma C at q=0); and A₀_{ν,ν} vanishes on D to order exactly 1 for ν ∈ 𝒟 and to
order 0 otherwise. The last holds because the S-factor (1 − w_a⁻²) is present iff a ∈ ν, and the A-factor
(1 − w_a/w_b) is present iff a ∉ ν ∋ b, each once. ∎

**Theorem 10.** For all N: (i) c₀^ν_{λμ} is a Laurent polynomial in w for every triple; (ii) if jump ≤ 1, then
c₁^ν_{λμ} is a Laurent polynomial in w. (iii) For every triple, every non-monomial pole of c₁ lies on a long-root
divisor, and it is simple.

*Proof.* ℓ⁰ = A₀⁻¹ = adj/det, and det A₀ = Π_I A₀_{I,I} is a product of L, S, A factors and monomials. So for a
Laurent f, the only possible non-monomial poles of ℓ⁰(f) lie on L, S, A divisors. The vectors P = φ⁰_λφ⁰_μ and Φ
(1.1) are Laurent. (A₁ = ∂_qÃ|₀ is Laurent, because Ã ∈ C[q][w^±] localised at units 1 − q·monomial. And c₀ is
Laurent by (i), which is proved first.)
- S, A divisors: P and Φ have the collapse property. φ⁰ and δ have it by Lemma C (differentiate in q), and c₀ is
  regular on D. Lemma C' applies.
- Long roots: c₀ is regular, because P ∈ L₁ (Corollary 4, Prop 5). For c₁, Theorem 7 gives Res = γ Res ℓ_ν(1_{⊇ab}P).
  The residue is Σ_{M⊇ab, M≤ν} ρ_{νM}·P(M'')|_H. When jump ≤ 1, every such M has |M''| ≤ |ν| − 2 < max(|λ|,|μ|),
  so P(M'') = 0 by (M3). Simplicity of the poles comes from Prop 5. ∎

So (d) is the first-order shadow of Clio's thm:main, and **at first order it needs no Gap q**: the (Q)-family is
not a divisor at q = 0. Theorem 10(i) also turns Rick's Day-242 census "195/195 Laurent at q = 0" into a theorem.

## 7. Why (c) fails: the high-block cancellation (N = 4)

Triple (λ, μ, ν) = (∅, {3}, {1,2,4}): jump 2, and Lemma-comb pair (1,2).
- Numeric engine (`num4.py`, independent of the §3–6 theory, Motegi summand with q symbolic): c₁ has no pole on any
  of the six H_ab. Confirmed by the Theorem 7 route (`resnum.py`) and on two further seeds.
- Full q (`fullc.py`): c^{124}_{∅,3} has a simple pole on H₁₂ only. Its residue is
  7q²(q−1)⁶·(quintic in q)/const: **q-valuation 2**. So Clio's converse (and her Conjecture "exact pole divisor")
  survives for this triple, but the pole is invisible at q¹.
- Mechanism. The defect set is D = {123, 124}, an interval and an up-set of ↓ν. The pair (1,2) has an empty low
  block and a high block {3} ↔ {4}. The ν-coordinate is a divided difference over the high block, and at q=0 the
  restricted function φ⁰₃|_D is proportional (through the label-free factor Γ) to the basis element φ⁰_{123}|_D.
  So the q¹ residue is 0. In Theorem 9 it is the *low* block's weight L(M'') that forces non-vanishing. A
  high block alone gives a divided difference that can vanish at q = 0.

Also: (b) fails because φ⁰ is (by every check: diagonal = the codim-many tangent weights, support = containment of
strict partitions, GKM) the K_T(LG(N,2N)) Schubert basis. In cohomology σ₂₁·σ₄₃ = σ₄₃₂₁ in LG(4,8) is a jump-2
product, and c₀^{1234}_{12,34} = 1/(w₃²w₄²) is its equivariant K-lift (`c0_C4.txt`). The identification with K_T(LG)
is a **hunch** (strong): nothing above uses it.

## 8. Computational record

- q=0, N=3: no jump-≥2 c₀ ≠ 0; all c₀ Laurent. N=4: 17 jump-≥2 c₀ ≠ 0 (list in `c0_C4.txt`), all c₀ Laurent.
- c₁, N=2,3 (symbolic, `c1N.py` with engine3, validated vs engine2 vs the Day-242 engine): every jump-≥2 triple has a
  c₁ pole (H₁₂ always among them); every jump-≤1 c₁ is Laurent.
- N=4 numeric (`num4.py`), H₁₂: jump 2: 107 poles, 1 none; jump 3: 18/18; jump 4: 1/1; jump ≤ 1: 0 poles in 599.
  Remaining pairs: `num_N4_{13,14,23,24,34}.txt`, see the addendum below.
- Controls: QC fails off H. The ratio is I-independent only on H (a single-I comparison would have been blind,
  which is why every I was compared). The top-two formula predicts zeros on the non-covered triples, and they are
  zero.

## 9. Gaps (precise)

1. Per-ν converse for triples not covered by Theorem 9. Needed: for some pair (a,b), the full-q residue
   (g² − g)·Res ℓ_ν((1−1/g)… ) ≠ 0. That is a statement about the *q-valuation* of a high-block divided
   difference, and order q¹ is provably insufficient. Candidate tools: (i) the q ↔ 1/q symmetry of g, which
   suggests the q→∞ expansion, where Ã·q^{−(dim − codim I)} has a finite limit (probed at N=2), may cover high
   blocks; (ii) the block factorisation A|_V = diag(Γ)(A^{high} ⊗ A^{low}), valid when the pair is consecutive in
   ν and V is an interval with fixed block sizes (sketched, verified on the N=6 example (5,5,{1,3,4,6}), where the
   poles are exactly Clio's predicted pairs).
2. (d) is closed (Thm 10). Not closed: (d) at full q, which is Clio's thm:main and still needs her Gap q.
3. Theorem 6 and Prop 5 need "no factor of Ã_{I,I} other than the long root vanishes on H_ab", which I read off
   (M4) and Clio's eq:diag (and checked at N ≤ 3). Fine, but noted.
