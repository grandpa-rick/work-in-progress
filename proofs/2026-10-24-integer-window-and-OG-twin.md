# The integer window, exactly, on LG and OG: monotone hyperplane induction (PROVE 253)

Rick, 2026-10-23 (PROVE 253; the file name is the date PROVE.md asked for). Conventions follow
`proofs/2026-10-24-Z-by-hyperplane-induction.md` (Theorem Z, "the Z file") and `proofs/2026-10-23-odd-zero-law.md`.
Scripts are in `work-in-progress/scripts/prove253/` (`check.py`, `mono.py`, `drop.py`, `strip.py`, `negctl.py`).

## 0. Verdict

PROVE.md asked for the vanishing halves (Z_int) and (Z_OG). The session got **more**: the exact integer zero sets, with signs,
on both LG and OG, by **one** induction. The integer case has a structural advantage the half-integer case lacks. The anchor
value is **0**, not a nonzero number with the "wrong" sign, so BR's alternating Pieri signs make (−1)^{|κ|}H_κ(−j)
**monotone** in j. That is exactly the sign-coherence that died at half-integers (Z file §5(f)).

| item | statement | grade |
|---|---|---|
| Lemma 0 (chain formula) | r(κ) = max_i (κ_i + 2i − 2) = 1 + max{a+b : (a,b) ∈ D(κ)} | proved (§2), checked m ≤ 7 |
| Lemma 1 (anchor) | H_κ(−1) = 0 for every κ ≠ ∅, on LG and on OG | proved given (G1) (§3) |
| Lemma 2 (Pieri support + sign) | dimension-indexed K-Pieri for O^1 on LG and OG: support, sign, and OG coefficient exactly | proved given (BR), first-hand (§4) |
| Lemma Ch′ (chain drop, full diagram) | LG: r(ρ) ≥ r(κ)−2 on the support; OG: r(ρ) ≥ r(κ)−1 | proved (§5), checked m ≤ 7 |
| Lemma A (antidiagonal strip) | stripping the top 1 (r odd) or 2 (r even) antidiagonals gives c_ρ ≠ 0 on LG; stripping 1 gives c_ρ = (−1)^{|θ|−1} on OG | proved (§6), checked m ≤ 7 against BR's h-sum formula; negative control §8 |
| **Theorem I (LG)** | s_j := (−1)^{|κ|}H_κ(−j) satisfies 0 = s_1 = ⋯ = s_W < s_{W+1} ≤ s_{W+2} ≤ ⋯, W = ⌈r(κ)/2⌉. So (Z_int) and (NZ_int) both hold, with sign (−1)^{|κ|} | **proved given (G0), (G1), (BR)** (§7) |
| **Theorem I_OG** | the same with W = r(κ) | **proved given (G0), (G1), (BR)** (§7) |

Theorem I discharges **gap 4** of `2026-10-23-odd-zero-law.md` §7 (the integer window used in Prop. R's application to (Z)), and
it does so without the LS-parity model. It also proves (NZ_int), which PROVE.md listed as a stretch goal. Theorem I_OG is
Crown E (Day 252 dream) made a theorem: zero-window length = Pieri support height (2 on LG, 1 on OG).

## 1. Setup

X is LG(m,2m) (G = Sp_{2m}) or OG(m,2m+1) (G = SO_{2m+1}). L = L_{ω_m} is the ample generator of Pic X: the Plücker O(1) on LG,
the spinor line bundle on OG (`wake253/og.py` convention). Schubert varieties follow BR 1005.2605 §4 (OG) / §5 (LG):
X^λ = X^λ(E_•) has codimension |λ|, for λ strict with λ₁ ≤ m, and λ^∨ = [m]∖λ. **X_κ := X^{κ^∨}**, of dimension |κ|.
H_κ(j) := χ(X_κ, L^j), a polynomial in j of degree |κ| (Snapper). D(κ) is the shifted diagram (Z file §1). r(κ) is the
longest chain in D(κ) under the product order. W(κ) := ⌈r(κ)/2⌉ on LG and r(κ) on OG, with W(∅) = 0.

Inputs (standard, cited, not deep-read; same family as the Z file):
- **(G0)** Pic X = ℤ·L, and the Schubert divisor X^{(1)} is the zero scheme of a section of L. Hence O^1 := [O_{X^{(1)}}] = 1 − [L^{−1}]
  in K(X). (Chevalley; X smooth, so the divisor is Cartier with ideal sheaf L^{−1}.)
- **(G1)** Demazure character formula for line bundles on Schubert varieties (any G).
- **(BR)** Buch–Ravikumar, arXiv:1005.2605, read first-hand from the v1 LaTeX source (§4 quoted below).

## 2. Lemma 0 (the chain formula)

**Lemma 0.** r(κ) = 1 + max{a+b : (a,b) ∈ D(κ)} = max_{1≤i≤ℓ(κ)} (κ_i + 2i − 2), and r(∅) = 0.

*Proof.* Along a chain x₁ < x₂ < ⋯ < x_r in the product order, the coordinate sum a+b strictly increases. It starts at ≥ 0
and ends at ≤ max(a+b), so r ≤ 1 + max(a+b). Conversely, take (a,b) ∈ D(κ). The staircase path (0,0),(0,1),(1,1),(1,2),…,
(a,a),(a,a+1),…,(a,b) has a+b+1 boxes. Each satisfies i ≤ j and lies below (a,b) in the product order. D(κ) is an order
ideal of Δ_m, so every box on the path is in D(κ), and the path is a chain. The last box of row i−1 (0-indexed) is
(i−1, i+κ_i−2), with coordinate sum κ_i + 2i − 3. ∎

Example: r(4,3,2,1) = max(4,5,6,7) = 7. r(3,2) = max(3,4) = 4.

**Fact (order ideals).** Every order ideal I of Δ_m equals D(ρ) for a unique strict ρ.
*Proof.* (r,c) ∈ I with c > r implies (r,c−1) ∈ I, so each row of I is an interval starting on the diagonal.
(r,c) ∈ I with r ≥ 1 implies (r−1,c) ∈ I, so row r−1 reaches column ≥ c. Therefore each nonempty row r ≥ 1 forces row r−1
to be nonempty and strictly longer. ∎

## 3. Lemma 1 (the anchor H_κ(−1) = 0)

**Lemma 1.** If κ ≠ ∅ then H_κ(−1) = 0, on LG and on OG.

*Proof.* X_κ is a G-translate of a B-stable Schubert variety \overline{BwP/P}, where w ∈ W^P is a minimal coset representative
of length |κ| ≥ 1. L is G-linearised, so χ(X_κ, L^j) is unchanged by translation. By (G1) and Lemma P of the odd-zero-law
file (stated and proved there for any root system), H_κ(j) = [D_w e^{jω_m}](h=0) as polynomials in j. Since w is minimal in
wW_P, ws_i > w for every s_i ∈ W_P, i.e. every i ≠ m. Since w ≠ e, w has some right descent, so ws_m < w. Write
w = w′s_m (lengths add). Then D_w = D_{w′}D_m, and D_m acts first. Now ⟨−ω_m, α_m^∨⟩ = −1 by the definition of the fundamental
weight. With s_m e^λ = e^{λ+α_m} this gives D_m e^{−ω_m} = (e^{−ω_m} − e^{−α_m}e^{−ω_m+α_m})/(1 − e^{−α_m}) = 0. Hence
D_w e^{−ω_m} = 0 as a Laurent polynomial, and H_κ(−1) = 0. ∎

(On LG this is the remark "D_m e^{−ω} = 0" already used in Thm C. Lemma W of the odd-zero-law file exhibits the reduced word
ending in s_m explicitly. The coset argument does both types at once.)

## 4. Lemma 2 (the Pieri input, first-hand)

BR define c^ν_{λ,p} by O^λ·O^p = Σ_ν c^ν_{λ,p} O^ν (eq. (kcoef) of §1). Quotes from the v1 source:

- **OG, §4.** "c^ν_{λ,p} = 𝓑(ν/λ,p)" for λ ⊂ ν (unlabelled Cor after eq. (hsumB)). **Lemma vanishB**: "(a) If θ is not a rim
  then 𝓑(θ,p)=0 for all p. (b) If θ contains a row or column with a boxes, then 𝓑(θ,p)=0 for p<a. (c) If θ is a rook strip,
  then 𝓑(θ,1)=(−1)^{|θ|−1}." A rook strip has at most one box in each row and each column (§3).
- **LG, §5, Cor pieriCtab.** "The constant c^ν_{λ,p} for K(LG(n,2n)) is equal to (−1)^{|ν/λ|−p} times the number of
  KLG-tableaux of shape ν/λ with content {1,2,…,p}." KLG-tableau (verbatim conditions): labels from {1′<1<2′<2<⋯}, "(i) each
  row of θ is strictly increasing from left to right; (ii) each column of θ is strictly increasing from top to bottom;
  (iii) each box containing an unprimed integer must be larger than or equal to all boxes southwest of it, (iv) each box
  containing a primed integer must be smaller than or equal to all boxes southwest of it, and (v) no diagonal box contains a
  primed integer. If θ is not a rim, then there are no KLG-tableaux of shape θ."
- **λ ⊄ ν.** c^ν_{λ,p} = χ(O^λ·O^p·O^*_ν), where O^*_ν is a ℤ-combination of Richardson classes O_τ = O^τ·O^{… }_{opp} with
  τ ⊆ ν (eqs. (dualA), (dualB), and the LG analogue at line ~1426). By BR §4, X^λ(E_•) ∩ X^μ(E^op_•) ≠ ∅ iff λ ⊆ μ^∨. So
  O^λ·O_τ = 0 unless λ ⊆ τ ⊆ ν. Hence c^ν_{λ,p} = 0 unless λ ⊆ ν. BR state this in type A ("In particular, we have
  c^ν_{λμ}=0 whenever λ ⊄ ν") with the same proof. Also c^λ_{λ,1} = 0 (codimension).

**Lemma 2 (dimension indexing).** Write [O_{X_κ}]·O^1 = Σ_ρ c_ρ [O_{X_ρ}] and θ := D(κ)∖D(ρ).
- (LG) c_ρ ≠ 0 ⇒ ρ ⊊ κ, θ is a rim (no two boxes strictly SE of each other), and θ has ≤ 2 boxes in every row and column.
  Moreover c_ρ = (−1)^{|θ|−1} N(θ), where N(θ) ≥ 0 is the number of KLG-tableaux of shape R(θ) with content {1}.
- (OG) c_ρ ≠ 0 ⇒ ρ ⊊ κ and θ is a rook strip. Conversely c_ρ = (−1)^{|θ|−1} for **every** ρ ⊊ κ with θ a rook strip.
- In both cases a rook strip θ = D(κ)∖D(ρ) is an antichain (no two boxes comparable).

*Proof.* c_ρ = c^{ρ^∨}_{κ^∨,1}. By Lemma S of the Z file, θ = R(D(ρ^∨)∖D(κ^∨)) with R(r,c) = (m−1−c, m−1−r), an order-reversing
involution of Δ_m that swaps rows and columns and preserves the diagonal. So R preserves "rim", "≤ k per row/column", and
"rook strip". It also preserves |θ|, and |θ| = |κ| − |ρ|.
- LG: Cor pieriCtab with p = 1. Labels lie in {1′,1}, and strict rows/columns allow ≤ 2 boxes per row/column. Non-rims have no
  KLG-tableaux.
- OG: Lemma vanishB (a),(b) with p = 1 < a = 2 give the support statement. (c) gives the value.
- Antichain: suppose x = (i,j) and y = (k,l) in θ with x < y and i < k, j < l. Then (i,l) is in the shifted plane
  (i ≤ j < l). It is ≤ y ∈ D(κ), so it lies in D(κ). It is ≥ x ∉ D(ρ), so it is not in D(ρ). Hence (i,l) ∈ θ, in the same
  row as x, and θ is not a rook strip. If instead x < y share a row or a column, θ is not a rook strip either. ∎

## 5. Lemma Ch′ (chain drop on the full diagram)

**Lemma Ch′.** If c_ρ ≠ 0, then every chain of D(κ) meets θ in ≤ 2 boxes (LG) or ≤ 1 box (OG). Hence r(ρ) ≥ r(κ) − 2 (LG),
resp. r(κ) − 1 (OG).

*Proof.* LG: the proof of Lemma Ch in the Z file §3 is a statement about any chain C of D(κ) and any θ that is a rim with
≤ 2 boxes per row/column. It never uses row 0. OG: θ is an antichain (Lemma 2). Now take a longest chain C of D(κ). Then
C∖θ ⊆ D(ρ) is a chain of the claimed size. ∎

(PROVE.md step 2 asked whether Lemma Ch is about P_κ in general or only about P_κ minus row 1. It is general. `drop.py`:
0 violations over 1337 LG-admissible and 449 OG-admissible pairs at m = 7. The maximum drop is exactly 2 on LG and 1 on OG,
so both bounds are sharp.)

## 6. Lemma A (a support term that drops the window)

For κ ≠ ∅ and t ∈ {1,2}, put T_t := {(a,b) ∈ D(κ) : a+b ≥ r(κ)−t} and D(ρ_t) := D(κ)∖T_t.

**Lemma A.**
- (i) D(κ)∖T_t is an order ideal of Δ_m, so ρ_t is a strict partition, ρ_t ⊊ κ, and r(ρ_t) ≤ r(κ)−t.
- (ii) T_t is a rim with ≤ 2 boxes per row/column. T_1 is an antichain.
- (iii) LG: put t* = 1 if r(κ) is odd and t* = 2 if r(κ) is even. Then N(T_{t*}) ≥ 1, so c_{ρ_{t*}} ≠ 0.
- (iv) OG: c_{ρ_1} = (−1)^{|T_1|−1}.
- (v) W(ρ_{t*}) ≤ W(κ) − 1 on LG, and W(ρ_1) ≤ W(κ) − 1 on OG.

*Proof.*
- (i) T_t is an up-set of D(κ), so its complement is an order ideal of D(κ), hence of Δ_m. Apply the Fact of §2.
  T_t ∋ the box realising max(a+b) (Lemma 0), so ρ_t ⊊ κ. Every box of ρ_t has a+b ≤ r−t−1, so r(ρ_t) ≤ r−t by Lemma 0.
- (ii) By Lemma 0, T_t lies on the antidiagonals a+b ∈ {r−t, …, r−1}, at most two of them. Two boxes in strictly-SE
  position differ in a+b by ≥ 2, so T_t is a rim. A row or column meets each antidiagonal at most once, which gives
  ≤ 2 boxes per row/column. T_1 lies on one antidiagonal, so it is an antichain.
- (iii) We exhibit a KLG-tableau of shape R(T_{t*}) with content {1}. With labels from {1′,1}, conditions (iii)/(iv) always
  hold: 1 is ≥ every label and 1′ is ≤ every label. Content {1} just needs T ≠ ∅.
  - *r odd, t* = 1.* T_1 is an antichain, so R(T_1) is one too: no two boxes share a row or column. Label every box 1. Rows
    and columns are trivially strict, and there are no primes, so (v) holds.
  - *r even, t* = 2.* Label a box of R(T_2) with 1′ if it is the image of an *outer* box (a+b = r−1), and with 1 if it is the
    image of an *inner* box (a+b = r−2). A row of R(T_2) with two boxes is the image of a column pair (a,b),(a+1,b) of T_2,
    inner over outer. R maps these to (m−1−b, m−1−a) and (m−1−b, m−2−a), so the outer image is on the left: 1′ then 1,
    strictly increasing. A column of R(T_2) with two boxes is the image of a row pair (a,b),(a,b+1), inner then outer. Their
    images are (m−1−b, ·) and (m−2−b, ·), so the outer image is on top: 1′ over 1, strictly increasing. For (v): primed boxes
    are images of outer boxes. Those have a+b = r−1, which is **odd**, so they are off the diagonal, and R preserves the
    diagonal.
- (iv) T_1 is an antichain, hence a rook strip. Apply Lemma 2 (OG).
- (v) LG, r odd: r(ρ₁) ≤ r−1, which is even, so ⌈r(ρ₁)/2⌉ ≤ (r−1)/2 = W−1. LG, r even: r(ρ₂) ≤ r−2, so
  ⌈r(ρ₂)/2⌉ ≤ r/2 − 1 = W−1. OG: r(ρ₁) ≤ r−1 = W−1. ∎

*Why t* depends on parity (negative control).* For r odd, the 2-strip T_2 puts the diagonal box (a,a) with 2a = r−1 on the
**outer** antidiagonal. If (a−1,a) ∈ T_2 as well, (v) forbids the forced prime and N = 0. `negctl.py`: on the r-odd 2-strips,
"c = 0 ⟺ such a conflict" matches BR's h-sum formula 67/67 (m ≤ 6). 42 of the 67 are zero, so the parity choice is
load-bearing and the test separates.

## 7. The theorems

**Theorem I.** Assume (G0), (G1), (BR). On X = LG(m,2m) (W = ⌈r/2⌉) or X = OG(m,2m+1) (W = r), for every strict κ with parts in
[m], the numbers s_j(κ) := (−1)^{|κ|}H_κ(−j) satisfy

  s_j(κ) = 0 for 1 ≤ j ≤ W(κ), and 0 < s_{W+1}(κ) ≤ s_{W+2}(κ) ≤ s_{W+3}(κ) ≤ ⋯.

In particular, the zero set of H_κ on the negative integers is exactly {−1, …, −W(κ)}. For LG this is (Z_int) plus (NZ_int).
For OG it is (Z_OG) plus nonvanishing. Beyond the window, sgn H_κ(−j) = (−1)^{|κ|}.

*Proof.* We use strong induction on |κ|. The claim P(κ) is: s_j(κ) ≥ 0 for all integers j ≥ 1, with equality iff j ≤ W(κ).
For κ = ∅, X_∅ is a point, H ≡ 1, W = 0, and P(∅) holds. Let κ ≠ ∅, so r ≥ 1 and W ≥ 1.

*Recursion.* By (G0), for every integer j,
  b_κ(j) := H_κ(j) − H_κ(j−1) = χ([O_{X_κ}]·(1−[L^{−1}])·[L^j]) = Σ_ρ c_ρ H_ρ(j).
(No section or Tor is needed, since L^{−1} is a line bundle.) Put j ↦ −j, multiply by (−1)^{|κ|}, and use
c_ρ = (−1)^{|θ|−1}|c_ρ| (Lemma 2) together with |κ| = |ρ| + |θ|:

  **s_{j+1}(κ) = s_j(κ) + Σ_ρ |c_ρ|·s_j(ρ)**  for every integer j.   (★)

*Sign of the increments.* For j ≥ 1, every ρ in the sum has |ρ| < |κ|, so s_j(ρ) ≥ 0 by P(ρ). Hence s is nondecreasing on
j ≥ 1.

*Anchor.* s_1(κ) = ±H_κ(−1) = 0 by Lemma 1. So s_j(κ) ≥ 0 for all j ≥ 1.

*Vanishing.* Let 1 ≤ j ≤ W−1. For every ρ with c_ρ ≠ 0, Lemma Ch′ gives W(ρ) ≥ W−1 ≥ j. On LG this is
⌈r(ρ)/2⌉ ≥ ⌈(r−2)/2⌉ = W−1. On OG it is r(ρ) ≥ r−1. So s_j(ρ) = 0 by P(ρ), and (★) gives s_{j+1} = s_j. From s_1 = 0 this
reaches s_W = 0.

*Nonvanishing.* At j = W, (★) gives s_{W+1}(κ) ≥ |c_{ρ*}|·s_W(ρ*), where ρ* = ρ_{t*} (LG) or ρ₁ (OG). The other terms are
≥ 0. By Lemma A, c_{ρ*} ≠ 0 and W(ρ*) ≤ W−1 < W, so s_W(ρ*) > 0 by P(ρ*). Hence s_{W+1} > 0, and by monotonicity s_j > 0 for
all j > W. ∎

**Worked example (LG, κ = (3,2), m = 3).** D(κ) = {(0,0),(0,1),(0,2),(1,1),(1,2)} has coordinate sums 0,1,2,2,3. So r = 4
(Lemma 0: max(3, 2+2)) and W = 2. BR's formula (`brpieri.pieri_dim`) gives

  [O_{X_{(3,2)}}]·O^1 = 2[O_{X_{(3,1)}}] − [O_{X_{(2,1)}}] − [O_{X_{(3)}}] + [O_{X_{(2)}}].

The θ's are {(1,2)}, {(0,2),(1,2)}, {(1,1),(1,2)}, {(0,2),(1,1),(1,2)}, with signs (−1)^{|θ|−1} = +,−,−,+ as Lemma 2 says.
The r's are 3, 3, 3, 2, so every drop is ≤ 2 (Lemma Ch′). r is even, so t* = 2 and T_2 = {a+b ≥ 2} = {(0,2),(1,1),(1,2)},
ρ* = (2). Lemma A's tableau labels the outer box (1,2) with 1′ and the inner boxes (0,2), (1,1) with 1. The inner diagonal
box (1,1) is unprimed, as required. It is the unique labelling, matching c_{(2)} = +1. Then (★) at j = 2 reads
s_3(3,2) = s_2(3,2) + 2s_2(3,1) + s_2(2,1) + s_2(3) + s_2(2). The first four terms are 0, since the window of each is
W = 2 ≥ 2, and s_2((2)) > 0, since W((2)) = 1. So H_{(3,2)} vanishes at −1, −2 and not at −3.

**Example (κ = (4,3,2,1), m = 4).** r = 7 (odd), W = 4. The support is the single term (4,3,2), with θ = {(3,3)}, which is T_1
(the only box with a+b = 6). It is diagonal and gets label 1, and c = +1. r(4,3,2) = 6, so W(ρ*) = 3 < 4. Theorem I: zeros at
−1..−4, and (−1)^{10}H(−5) > 0.

## 8. Verification

All checks use **Demazure characters directly** (`check.py`, `mono.py`: types C and B, interpolated to negative j). They are
independent of the LS-parity model and of the Ω/SMT model.
- `check.py m`, m = 2,3,4: LG integer zero set = {1..⌈r/2⌉}: 4/4, 8/8, 16/16. OG zero set = {1..r} and
  H^OG(−b) = (−1)^{|κ|}·#{strict order-preserving P_κ → [b−1]}: 4/4, 8/8, 16/16. m = 5 result in §8a.
- `mono.py m`, m = 2,3,4: s_j ≥ 0, nondecreasing, and = 0 iff j ≤ W, for j ≤ |κ|+7. LG 3/3, 7/7, 15/15; OG 3/3, 7/7, 15/15.
  m = 5 in §8a.
- `drop.py 7`: Lemma Ch′, 0 violations (§5).
- `strip.py 7`: Lemma 0 127/127. c(ρ*) ≠ 0 with sign (−1)^{|θ|−1}: 127/127. All c sign-alternating: 127/127 (BR's h-sum
  formula `brpieri.C`, which is independent of the KLG count used in the proof).
- `negctl.py 6`: 67/67 (§6).
- Wake 253, earlier: LG zero set exactly {−1..−⌈r/2⌉} 127/127 at m = 7 (LS model). OG zero set exactly {−1..−r} 63/63 at
  m = 6 (Ω model, checked against Demazure).

## 9. Remarks

- **OG via SMT, a second proof of Theorem I_OG.** X^OG is minuscule. Standard monomial theory gives H^OG_κ(j) = Ω(P_κ, j+1),
  the order polynomial (`og.py` checks this against Demazure). Stanley reciprocity then gives
  H^OG_κ(−b) = (−1)^{|κ|}Ω̄(P_κ, b−1), the count of strict order-preserving maps, which is nonzero iff b−1 ≥ r. That is
  Theorem I_OG in one line, but it needs SMT (Lakshmibai–Seshadri / Proctor; not deep-read). The induction proof needs only
  (G0), (G1), (BR). In type A the same mechanism (Lenart rook strips, height 1) is Stanley reciprocity for Grassmannian
  Schubert varieties (Bruns–Herzog; Fong–Miyazaki 2502.09559, owner search Wake 253).
- **LG strict-map formula.** (−1)^{|κ|}H_κ(−j) = #{strict maps P_κ → [2j−1], even diagonal values} (computed 127/127). It is
  consistent with Theorem I, but it stays computed: it rests on the LS-parity model. Theorem I does not use it.
- **Why integers win and half-integers lose.** (★) holds at every j ∈ ½ℤ. At half-integers the anchor H_κ(−½) = 2^{−κ₁} is
  nonzero, and the Prop-R regime flips sgn H_ρ to (−1)^{|ρ′|}. So the increments of (★) have mixed signs: the Z file §5(f),
  κ = (3,2), −1/8 vs +1/4. At integers, everything starts at 0 and P(ρ) keeps every increment ≥ 0.

## 10. Gaps (precise)

1. (G0), (G1): standard, cited, not deep-read. These are the same debt as Theorem Z and Thm C.
2. (BR): read first-hand (v1 LaTeX source). BR leave the proof of Cor pieriCtab (the KLG count, which supplies the LG sign
   and Lemma A(iii)) "to the reader". We cite the stated corollary. Corroboration: BR's h-sum formula 𝓒(θ,1) (BR Cor 5.4,
   proved in full there) gives the same sign and the same nonvanishing of c_{ρ*} for every κ at m ≤ 7 (`strip.py`). The
   true K-class computed from (G1) matches at m ≤ 5 (Wake 252 `kclass.py`). The OG input (Lemma vanishB) is proved in full
   in BR.
3. Nothing else. Lemmas 0, 1, Ch′, A and the induction are self-contained.

## 11. Lesson

The same (★) that has no sign-coherence at half-integers has it for free at integers, because the anchor is **zero**. When an
alternating-sign recursion refuses to be coherent, look for the specialisation where the base case vanishes: there the signs
have nothing to fight.
