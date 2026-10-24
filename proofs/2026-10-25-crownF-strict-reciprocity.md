# Crown F as a theorem: the LG K-Pieri coefficient counts {e<o} maps, and H(−j) counts strict maps (PROVE 254)

Rick, 2026-10-24 (PROVE 254; file name is the date PROVE.md asked for). Conventions follow
`proofs/2026-10-24-integer-window-and-OG-twin.md` ("the I file") and `proofs/2026-10-24-Z-by-hyperplane-induction.md`
("the Z file"). Scripts: `work-in-progress/scripts/prove254/` (`srcheck.py`, `halfthm.py`, `hz_vs_rec.py`, `rowform.py`, `rowphi*.py`, `half.py`, `tmat.py`,
`tshape.py`, `tprime.py`; WIP f17e0cf) and
`work-in-progress/scripts/wake254/q30/` (`q30.py`, `brcheck.py`).

## 0. Verdict

| item | statement | grade |
|---|---|---|
| Lemma R (reduction) | #strict order-preserving maps P → {e<o} with D ↦ e is 0 if P has a 3-chain or a non-minimal element of D, else 2^{#isolated elements ∉ D} | proved (§2) |
| **Theorem Q30** | LG, dimension-indexed: [O_{X_κ}]·O¹ = Σ_ρ c_ρ[O_{X_ρ}] with c_ρ = (−1)^{|θ|−1}V(θ), V(θ) = #strict order-preserving θ → {e<o} with diagonal ↦ e, for every ρ ⊊ κ (and c_ρ = 0 for ρ ⊄ κ) | **proved given BR Cor 5.4 + Thm 5.5** (both proved in BR, first-hand) — §3 |
| Theorem Q30_OG | OG: c_ρ = (−1)^{|θ|−1} if θ is a nonempty antichain (= nonempty set of maximal boxes of P_κ), else 0 | proved given BR Lemma vanishB (I file Lemma 2) — §3.4 |
| **Theorem SR_LG** | (−1)^{|κ|}H_κ(−j) = #{strict order-preserving P_κ → [2j−1], diagonal values even}, all j ≥ 1 | **proved given (G0), (G1), (BR)** — §4 |
| **Theorem SR_OG** | (−1)^{|κ|}H^OG_κ(−j) = Ω̄(P_κ, j−1) = #{strict order-preserving P_κ → [j−1]}, all j ≥ 1 | **proved given (G0), (G1), (BR)** — §5 |
| Corollary (Theorem I, combinatorial proof) | zero window, sign and monotonicity of Theorem I read off from the count | proved — §6 |
| **Theorem HZ (odd-zero law, both halves)** | H_κ(−a−½) = 0 ⟺ a ≤ ⌈r(κ′)/2⌉; sign (−1)^{\|κ\|−κ₁} for 2a ≤ κ₁, (−1)^{\|κ\|} for 2a ≥ κ₁; \|value\| ≥ 2^{−κ₁}N_{κ′}(a) | **proved given (G0), (G1), (BR) + Thm A** — §9 |

The point: (★) of the I file is a positive recursion with a trivial anchor. Unrolled to the anchor, it IS the counting
formula. The only nontrivial input is that the one-step coefficient |c_ρ| counts one "layer" {e<o}, and that is Q30.
SR_LG no longer depends on the LS-parity model. SR_OG no longer depends on SMT/Stanley reciprocity.

## 1. Setup

As in the I file §1. X = LG(m,2m) or OG(m,2m+1), L the ample generator, X_κ := X^{κ^∨} of dimension |κ|,
H_κ(j) = χ(X_κ, L^j). P_κ = D(κ) = shifted diagram, boxes (a,b) with 0 ≤ a ≤ b, 0-indexed, diagonal = {a = b}.
Order: x ≤ y iff x is weakly NW of y (product order). For ρ ⊆ κ, θ = D(κ)∖D(ρ) is an up-set (order filter) of P_κ,
and every up-set arises this way for a unique strict ρ (I file §2, Fact).

Inputs: (G0) O¹ = 1 − [L^{−1}]; (G1) Demazure character formula; (BR) Buch–Ravikumar arXiv:1005.2605, read first-hand
from the v1 LaTeX source. Specifically used here: §5 eq. (11) + **Cor 5.4** (c^ν_{λ,p} = 𝒞(ν/λ,p) for λ ⊂ ν) and
**Thm 5.5** (the arm recursion for 𝒞), both proved in BR. We do **not** use Cor 5.6 (KLG tableaux), whose proof BR
leave to the reader; §3.3 below is that proof at p = 1. §4 OG: Lemma vanishB (proved in BR).

The I file's (★): with s_j(κ) := (−1)^{|κ|}H_κ(−j),

  s_{j+1}(κ) = s_j(κ) + Σ_{ρ ⊊ κ} |c_ρ| s_j(ρ)   for every integer j,   (★)

valid on LG and OG, given (G0) and the sign rule c_ρ = (−1)^{|θ|−1}|c_ρ| (re-proved below on LG as part of Q30).
Anchor: s_1(κ) = 0 for κ ≠ ∅ (I file Lemma 1, given (G1)); s_j(∅) = 1 (X_∅ is a point).

## 2. Lemma R (counting {e<o} maps)

Let P be a finite poset and D ⊆ P. Let V(P,D) = #{f : P → {e<o} : x < y ⇒ f(x) < f(y), and f(D) ⊆ {e}}.

**Lemma R.** V(P,D) = 0 if P has a chain x<y<z, or if some d ∈ D is not minimal. Otherwise V(P,D) = 2^{i}, where i is the
number of isolated elements (comparable to nothing else) not in D.

*Proof.* A 3-chain needs three distinct values. If x < d ∈ D then f(x) < f(d) = e is impossible. Otherwise P has height ≤ 2.
A non-isolated x is then either non-maximal (x < y for some y, forcing f(x) = e, f(y) = o) or non-minimal (forcing f(x) = o),
not both (that would be a 3-chain). These forced values satisfy every relation x < y (x non-maximal ↦ e, y non-minimal ↦ o).
A non-isolated d ∈ D is minimal, hence non-maximal, hence forced to e, consistent with D. An isolated element is free (2
choices) unless it is in D (1 choice). ∎

In P_κ: a diagonal box (c,c) of an up-set θ is non-minimal in θ iff (c−1,c) ∈ θ. (The boxes below (c,c) in P_κ are those
weakly NW of it, i.e. (a,b) with a ≤ b ≤ c, a ≤ c. If one of them, (a,b) ≠ (c,c), lies in θ, then (c−1,c) — which lies
between (a,b) and (c,c) when a ≤ c−1, b ≤ c — lies in θ since θ is an up-set. And a = c forces b = c.)

## 3. Theorem Q30 (LG K-Pieri at p = 1)

### 3.1 Translation to BR's side

c_ρ = c^{ρ^∨}_{κ^∨,1} (dimension vs codimension indexing). If ρ ⊄ κ then c_ρ = 0 (I file §4, "λ ⊄ ν"). If ρ ⊊ κ, Lemma S
of the Z file gives θ_BR := D(ρ^∨)∖D(κ^∨) = R(θ), R(r,c) = (m−1−c, m−1−r). R is an involution of Δ_m, fixes the diagonal
setwise, swaps rows and columns, and **reverses** the NW order. By Cor 5.4, c_ρ = 𝒞(θ_BR, 1).

For a shifted skew shape φ let K(φ) := #labellings φ → {1′ < 1} that strictly increase along every row (left to right)
and every column (top to bottom), with no 1′ on a diagonal box, and K(φ) := 0 if φ is not a rim. (BR: a rim has no two
boxes with one strictly south and strictly east of the other.)

**Lemma T (rim and chains).** Let φ = ν/λ be a shifted skew shape. (a) If φ is not a rim, φ has a 3-chain in the NW order.
(b) If φ is a rim, two boxes of φ are comparable iff they share a row or a column.
*Proof.* (a) a = (r,c), b = (r′,c′) ∈ φ with r < r′, c < c′. Row ends of ν are weakly decreasing (strict ν, shifted), so
(r,c′) ∈ ν. Rows of λ are initial segments of rows and (r,c) ∉ λ, so (r,c′) ∉ λ. Then a < (r,c′) < b. (b) comparable =
weakly NW; not strictly NW (rim) ⇒ same row or column. ∎

**Lemma 3.1.** For every ρ ⊊ κ: V(θ) := V(θ, diagonal) = K(θ_BR).
*Proof.* If θ_BR is not a rim it has a 3-chain (Lemma T(a)); R reverses order, so θ has one, and V(θ) = 0 = K(θ_BR)
(Lemma R). If θ_BR is a rim, by Lemma T(b) "strict along rows and columns" = "strict order-preserving" on θ_BR. Compose
with R (order-reversing, diagonal-preserving) and rename 1 ↦ e, 1′ ↦ o: strict order-preserving maps θ_BR → {1′<1}
with diagonal ↦ 1 correspond bijectively to strict order-preserving maps θ → {e<o} with diagonal ↦ e. ∎

So Q30 ⟺ 𝒞(φ,1) = T(φ) := (−1)^{|φ|−1}K(φ) for every nonempty shifted skew shape φ. (Then c_ρ = (−1)^{|θ|−1}V(θ).)

### 3.2 BR Thm 5.5 at p = 1

BR Thm 5.5 (verbatim in `proofs/scratch-2026-10-24-q30-br-formula.md` §1b): 𝒞(θ,p) = 0 if θ is not a rim;
𝒞(θ,p) = δ_{θ,∅} for p ≤ 0; for a nonempty rim, a = |θ| − |θ̂| (θ̂ = θ minus its north-east arm ψ):
(i) θ̂ = ∅, θ meets the diagonal: column ↦ δ_{p,|θ|} − δ_{p,|θ|−1}, row ↦ δ_{p,|θ|};
(ii) θ̂ = ∅, θ off the diagonal: 2δ_{p,|θ|} − δ_{p,|θ|−1};
(iii) θ̂ ≠ ∅, ψ connected to θ̂: 𝒞(θ̂,p−a) − 𝒞(θ̂,p−a+1);
(iv) θ̂ ≠ ∅, ψ not connected to θ̂: a = 1: 2𝒞(θ̂,p−1) − 2𝒞(θ̂,p); a > 1: 2𝒞(θ̂,p−a) − 3𝒞(θ̂,p−a+1) + 𝒞(θ̂,p−a+2).

At p = 1, using 𝒞(θ̂,q) = 0 for q ≤ 0 (θ̂ ≠ ∅):
(i) diagonal column: |θ| = 1 ↦ 1, |θ| = 2 ↦ −1, |θ| ≥ 3 ↦ 0; diagonal row: |θ| = 1 ↦ 1, else 0.
(ii) off-diagonal row/column: |θ| = 1 ↦ 2, |θ| = 2 ↦ −1, |θ| ≥ 3 ↦ 0.
(iii) a = 1 ↦ −𝒞(θ̂,1); a ≥ 2 ↦ 0.
(iv) a = 1 ↦ −2𝒞(θ̂,1); a = 2 ↦ +𝒞(θ̂,1); a ≥ 3 ↦ 0.

### 3.3 Proof that T(φ) = 𝒞(φ,1) (= BR Cor 5.6 at p = 1, which BR leave to the reader)

Induction on |φ|. Non-rim: both sides 0. Let φ be a nonempty rim.

Geometry used throughout. (g1) Rows and columns of a shifted skew shape are intervals. (g2) A diagonal box (c,c) is the
leftmost box of its row and the lowest box of its column (positions (c,c′<c) and (r>c,c) are not in the shifted plane).
(g3) Let r be the top row of φ and B = (r,c) its rightmost box. Then c is the maximal column of φ: a box (r′,c′) with
c′ > c has r′ > r, and B would be strictly NW of it. So nothing of φ lies above or right of B, and the only possible
neighbours of B are L = (r,c−1) and D = (r+1,c). (g4) L and D are not both in φ (D is strictly SE of L). B is BR's
"upper-right box"; the arm ψ is the largest row or column of the form φ ∩ (square with upper-right box B).
(g5) In a rim, comparable boxes share a row or column (Lemma T(b)), so K(φ) only sees row and column constraints.

**(i), (ii): φ = ψ is a row or a column.** Two labels ⇒ K = 0 for length ≥ 3. Diagonal column (diagonal box at the
bottom, g2): length 1: label 1, K = 1, T = 1; length 2: 1′ over 1, K = 1, T = −1. Diagonal row (diagonal box leftmost):
length 1: K = 1, T = 1; length 2: the leftmost box would need 1′, forbidden, K = 0. Off-diagonal: length 1: K = 2,
T = 2; length 2: K = 1 (1′,1), T = −1. All match §3.2 (i), (ii).

**(iii) and (iv), a ≥ 3.** ψ is a row or column of ≥ 3 boxes, so K(φ) = 0 = 𝒞(φ,1).

**(iii), a = 1.** ψ = {B} and B is adjacent to θ̂, so exactly one of L, D is in φ (g3, g4). Let S = (r+1,c−1). If S ∉ φ,
the 2×2 square at B meets φ in {L,B} or {B,D}, a row or column of 2, so a ≥ 2. Hence S ∈ φ, and S ∈ θ̂.
- L ∈ φ: L is directly above S, so every valid labelling of θ̂ has L = 1′ (column strict). B is right of L, so B = 1
  is forced and allowed (unprimed). B's column contains only B (D ∉ φ, g1). B's row: boxes left of L are also left of
  B, and they already satisfy their constraints with L = 1′ in θ̂ (a box left of L would need a label < 1′, impossible
  — then K(θ̂) = 0 and K(φ) = 0 too). So K(φ) = K(θ̂).
- D ∈ φ: S is directly left of D, so D = 1 in every labelling of θ̂. B above D forces B = 1′. B is off the diagonal:
  D = (r+1,c) is in the shifted plane, so c ≥ r+1. B's row contains only B (L ∉ φ). Again K(φ) = K(θ̂).
Then T(φ) = (−1)^{|θ̂|}K(θ̂) = −T(θ̂) = −𝒞(θ̂,1) by induction. ✓

**(iii), a = 2.** ψ = {L,B} (a row) or {B,D} (a column), and the 2×2 square meets φ in exactly ψ, so S, D ∉ φ (resp.
S, L ∉ φ). Row case: L's neighbours are (r,c−2), S (∉), and the box above L (∉, r is the top row). B's other neighbour D
∉ φ. Connection to θ̂ therefore goes through (r,c−2) ∈ φ: three boxes in row r, K(φ) = 0. Column case: D's neighbours are
(r+2,c), S (∉) and (r+1,c+1) (∉, column c is maximal); B's other neighbour L ∉ φ. So (r+2,c) ∈ φ: a column of three,
K(φ) = 0. Both match 𝒞 = 0. ✓

**(iv), a ∈ {1,2}: ψ not connected to θ̂.**
- No box of ψ shares a row or column with a box of θ̂: by (g1), walking along the common row/column interval from the ψ
  box to the θ̂ box we would cross a ψ–θ̂ adjacency.
- ψ contains no diagonal box. If ψ is a row in row r containing (r,r), then ψ = {(r,r),…,(r,c)} (g2), and any box
  (r′,c′) ∈ θ̂ has r′ > r (r is the top row), so r < r′ ≤ c′ ≤ c: column c′ is a column of ψ. Contradiction. If ψ is a
  column in column c containing (c,c), then ψ = {(r,c),…,(c,c)} (g2), and any (r′,c′) ∈ θ̂ has r ≤ r′ ≤ c′ ≤ c: row r′
  is a row of ψ. Contradiction. (θ̂ ≠ ∅.)
- By (g5) and the first bullet, the constraints on ψ and on θ̂ are independent: K(φ) = K(ψ)K(θ̂), with K(ψ) = 2 (a = 1,
  off-diagonal single box) or 1 (a = 2, off-diagonal domino: 1′,1).
- T(φ) = (−1)^{a+|θ̂|−1}K(ψ)K(θ̂) = (−1)^a K(ψ) T(θ̂): a = 1 ↦ −2𝒞(θ̂,1), a = 2 ↦ +𝒞(θ̂,1). ✓

All cases of Thm 5.5 are covered (θ̂ = ∅ iff φ is a row or column; otherwise (iii) or (iv)). ∎

**Theorem Q30 (proved given BR Cor 5.4, Thm 5.5, and Lemma S).** c_ρ = (−1)^{|θ|−1}V(θ) for all ρ ⊊ κ, with V given by
Lemma R: c_ρ = 0 iff θ has a 3-chain or a diagonal box (c,c) with (c−1,c) ∈ θ; otherwise |c_ρ| = 2^{#off-diagonal boxes
of θ incomparable to every other box of θ}.

Remark. Route (i) of the scratch (Cor 5.6 at p = 1 read literally: content {1} makes KLG (iii)/(iv) vacuous) gives the
same statement in one line, but rests on a reader-left proof. §3.3 replaces it.

### 3.4 OG twin

I file Lemma 2 (OG): c_ρ ≠ 0 ⇒ θ is a rook strip, and c_ρ = (−1)^{|θ|−1} for every rook strip θ = D(κ)∖D(ρ), from BR
Lemma vanishB (a)–(c), proved in BR. For an up-set θ: rook strip ⟺ antichain. (⇒ I file Lemma 2. ⇐: two boxes in one
row or column are comparable.) An up-set is an antichain iff it consists of maximal elements of P_κ, and every set of
maximal elements is an up-set. So the OG operator is inclusion–exclusion over nonempty sets of corners, |c_ρ| = 1. ∎

## 4. Theorem SR_LG

For a finite poset P with D ⊆ P and n ≥ 0, let N(P,D;n) = #{strict order-preserving f : P → [n] with f(D) ⊆ 2ℤ}.
Write N_κ(j) := N(P_κ, diagonal; 2j−1).

**Lemma 4.1 (peeling the top layer).** For j ≥ 1 and every κ: N_κ(j+1) = N_κ(j) + Σ_{ρ ⊊ κ} V(D(κ)∖D(ρ))·N_ρ(j).
*Proof.* Let f : P_κ → [2j+1] be strict order-preserving with even diagonal values. Put θ = f^{−1}{2j, 2j+1}.
- θ is an up-set: x ∈ θ, x < y ⇒ f(y) > f(x) ≥ 2j ⇒ f(y) ∈ {2j+1}, so y ∈ θ. Hence θ = D(κ)∖D(ρ) for a unique strict
  ρ ⊆ κ (I file §2 Fact applied to the down-set P_κ∖θ, an order ideal of Δ_m).
- f|θ is strict order-preserving into {2j < 2j+1} with diagonal boxes ↦ 2j (the only even value); i.e. a map counted by
  V(θ) under e = 2j, o = 2j+1.
- f|D(ρ) is strict order-preserving into [2j−1] with even diagonal: counted by N_ρ(j).
Conversely, given ρ ⊆ κ, a V(θ)-map g and an N_ρ(j)-map h, glue f = g ∪ h. Relations inside θ or inside D(ρ) hold. For
x ∈ D(ρ), y ∈ θ with x < y: f(x) ≤ 2j−1 < 2j ≤ f(y). The case y < x with y ∈ θ, x ∈ D(ρ) is impossible (θ is an up-set).
Diagonal values are even. The two constructions are inverse. ρ = κ (θ = ∅, V = 1) gives the term N_κ(j). ∎

**Theorem SR_LG.** Assume (G0), (G1), (BR). For every strict κ with parts in [m] and every integer j ≥ 1,
  (−1)^{|κ|}H_κ(−j) = N_κ(j) = #{strict order-preserving P_κ → [2j−1] with even values on diagonal boxes}.
*Proof.* Induction on j, simultaneously for all κ. j = 1: [1] has no even value, and every nonempty P_κ contains the
diagonal box (0,0), so N_κ(1) = δ_{κ,∅} = s_1(κ) (anchor). Step: by Theorem Q30, |c_ρ| = V(θ) for ρ ⊊ κ and c_ρ = 0
otherwise, so (★) and Lemma 4.1 are the same recursion: s_{j+1}(κ) = s_j(κ) + Σ_{ρ⊊κ} V(θ)s_j(ρ) and the same for N. ∎

**Hand check.** κ = (2,1): P = {(0,0) < (0,1) < (1,1)}, a 3-chain with both ends diagonal. N(j) = #{2 ≤ x < y < z ≤ 2j−1,
x, z even} = Σ_{z even} #{x even < z, y between}: j = 3 ([5]): (2,3,4) only ⇒ 1. (★): Pieri support of (2,1) is
{ρ = (2), θ = {(1,1)}, V = 1} (the others have V = 0: ρ = (1) has θ = {(0,1),(1,1)}, diagonal non-minimal; ρ = ∅ a 3-chain).
s_3(2,1) = s_2(2,1) + s_2(2) = 0 + 1 = 1 ✓ (s_2(2) = N_{(2)}(2) = #{2 ≤ x < y ≤ 3, x even} = 1). κ = (1): N(j) = #even in
[2j−1] = j−1, H_{(1)}(x) = x+1 ✓. κ = (2): N(j) = Σ_{x even ≤ 2j−2}(2j−1−x) = (j−1)², H_{(2)}(x) = (x+1)² ✓.

## 5. Theorem SR_OG

**Lemma 5.1.** For n ≥ 0: Ω̄(P_κ, n+1) = Ω̄(P_κ, n) + Σ_{∅≠A⊆max(P_κ)} Ω̄(P_κ∖A, n), and P_κ∖A = D(ρ_A).
*Proof.* The preimage of the top value n+1 under a strict order-preserving map is an up-set (as in 4.1) and an antichain
(strictness), i.e. a set A of maximal elements; the rest maps strictly into [n]. Gluing is as in 4.1. ∎

**Theorem SR_OG.** Assume (G0), (G1), (BR). (−1)^{|κ|}H^OG_κ(−j) = Ω̄(P_κ, j−1) for all j ≥ 1.
*Proof.* j = 1: Ω̄(P,0) = δ_{P,∅} = s_1 (anchor, I file Lemma 1, both types). Step: (★) on OG with |c_ρ| = 1 exactly on
ρ = ρ_A, A ≠ ∅ a set of corners (§3.4), is Lemma 5.1 with n = j−1. ∎

(This is Stanley reciprocity for the order polynomial H^OG_κ(j) = Ω(P_κ, j+1) — the SMT model of the I file §9 — but
proved without SMT.)

## 6. Corollary: Theorem I, combinatorially

Let r = r(κ) be the longest chain of P_κ; by I file Lemma 0 it is realised by a staircase path starting at (0,0) and
ending at a box with a+b = r−1, and max(a+b) = r−1.
- **LG.** N_κ(j) > 0 ⟺ 2j−1 ≥ r+1 ⟺ j ≥ ⌈r/2⌉ + 1. (⇐) f(a,b) = a+b+2 is strict order-preserving (a+b strictly increases
  along x < y), takes values in [2, r+1], and f(a,a) = 2a+2 is even. (⇒) Along the staircase chain, f((0,0)) ≥ 2 (even)
  and f increases strictly over r boxes, so the last value is ≥ r+1. Since 2j−1 is odd: r even ⇒ j ≥ r/2+1; r odd ⇒
  2j−1 ≥ r+2 ⇒ j ≥ (r+3)/2 = ⌈r/2⌉+1.
- **OG.** Ω̄(P_κ, j−1) > 0 ⟺ j−1 ≥ r (a chain of r needs r values; f(a,b) = a+b+1 works).
- Monotonicity in j is inclusion of alphabets ([2j−1] ⊂ [2j+1] preserving parity).
With SR this re-proves Theorem I (zero set exactly {−1,…,−W}, sign (−1)^{|κ|} beyond, s_j nondecreasing) — but the
inputs are the same (G0), (G1), (BR) plus the full Q30, so this is a second proof, not a weaker one.

## 7. Verification

- `prove254/srcheck.py m`, m = 2..5: H_κ(−j) computed **directly from Demazure characters at negative weights**
  (`prove249d/dem.py` type C; `wake253/og.py` type B), j = 1..4, against brute-force strict counts. LG and OG:
  16/16, 32/32, 64/64, 128/128 (κ, j) pairs, 0 failures. Independent of (★), of BR, and of the LS/SMT models.
- Negative controls at m = 4 (64 pairs): LG alphabet [2j] instead of [2j−1]: 34 failures; LG without the diagonal parity:
  34 failures; OG alphabet [j]: 17 failures. The test separates.
- Wake 254: Q30 1772/1772 nonzero coefficients m ≤ 7, support exact on 8532 pairs (`q30.py`); `brcheck.py`: BR eq. (11)
  = literal KLG count = Lemma R closed form on all 2225 shifted skew shapes 2 ≤ n ≤ 6.

## 8. Gaps

1. (G0), (G1): standard, cited, not deep-read (same debt as Theorems Z, I, C).
2. (BR): Cor 5.4, Thm 5.5 (LG) and Lemma vanishB (OG) read first-hand, proved in BR. Cor 5.6 is NOT used.
3. §3.3 is a case analysis I wrote and cold re-read once (this session). It has not been refereed by anyone else.
4. Nothing else.


## 9. Stretch (Q31) — the half-integers: the whole odd-zero law, with signs, from one row

PROVE.md's stretch: iterate (★) from the half-integer anchor H(−½) = 2^{−κ₁} and look for positivity. It went through.
Notation: s_x(κ) := (−1)^{|κ|}H_κ(−x) for x ∈ ½ℤ; κ′ := κ minus its largest part; r′ := r(κ′); M := 2a−1; n := κ₁.

**Theorem HZ (proved given (G0), (G1), (BR)).** For every κ ≠ ∅ and every integer a ≥ 1:
- (small a) if 2a ≤ κ₁: s_{a+½}(κ) = (−1)^{κ₁}2^{−κ₁}·N_{κ′}(a);
- (large a) if 2a ≥ κ₁: s_{a+½}(κ) = Σ_g Φ(u(g)) with every Φ(u(g)) ≥ 2^{−κ₁}, so s_{a+½}(κ) ≥ 2^{−κ₁}N_{κ′}(a);
- hence H_κ(−a−½) = 0 ⟺ N_{κ′}(a) = 0 ⟺ κ′ ≠ ∅ and a ≤ ⌈r′/2⌉ (§6), and otherwise
  sgn H_κ(−a−½) = (−1)^{|κ|−κ₁} if 2a ≤ κ₁, (−1)^{|κ|} if 2a ≥ κ₁ (the two agree at 2a = κ₁, κ₁ even).
In the odd-zero-law notation (H(−M/2), M odd): zeros exactly at M = 3, 5, …, 2⌈r′/2⌉+1, none for single rows. That is
the **odd-zero law, both halves** — (Z) (Theorem Z had it by another induction) and (NZ) (open since Wake 251; PROVE 253
reduced it to an odd-μ₁ base inequality, which is now moot), **with strictness and signs**. The small-a bullet is Prop. R
of the odd-zero-law file, re-proved here without the residue calculus.

Inputs: (★) as a polynomial identity (it holds at all integers, both sides are polynomials, so it holds at every x),
Theorem Q30, Thm A of `2026-10-23-odd-zero-law.md` (H_κ(−½) = 2^{−κ₁}, proved given (G1)), and §6.

### 9.1 Layer operators

For up-sets θ = D(κ)∖D(ν) put A_o(κ,ν) = [θ is an antichain with no diagonal box], A_e(κ,ν) = [θ is an antichain].

**Lemma 9.1.** (a) A_oA_e = I + C, where C(κ,ν) = |c_ν|. (b) For a ≥ 1, [(A_oA_e)^{a−1}A_o](κ,ν) = N_{κ/ν}(a) :=
#{strict order-preserving D(κ)∖D(ν) → [2a−1], even values on diagonal boxes}.
*Proof.* (b) Peel values from the top, exactly as in Lemma 4.1: the preimage of the current top value is an antichain
up-set of what remains; the value is odd (no diagonal box allowed: A_o) or even (A_e). The values 2a−1, 2a−2, …, 1
alternate o, e, …, o, giving A_o(A_eA_o)^{a−1} = (A_oA_e)^{a−1}A_o. Gluing is as in 4.1. (a) The same with the
two-letter alphabet {e<o}: (A_oA_e)(κ,ν) = V(θ) (θ = ∅ gives 1), and V(θ) = |c_ν| for ν ⊊ κ by Q30. ∎

So (★) reads s_{x+1} = A_oA_e s_x for every x ∈ ℂ.

**Lemma 9.2 (the anchor collapses to single rows).** ε := A_e s_{½} satisfies ε(∅) = 1, ε((L)) = (−1)^{L+1}2^{−L} for
L ≥ 1, and ε(κ) = 0 if ℓ(κ) ≥ 2.
*Proof.* s_{½}(ν) = (−1)^{|ν|}2^{−ν₁} (Thm A). ε(κ) = Σ_{A ⊆ max P_κ} s_{½}(κ∖A). If ℓ(κ) ≥ 2, the last box z of the
bottom row is maximal and not in row 0; toggling z ∈ A changes |ν| by one and leaves ν₁ alone, so the terms cancel in
pairs. If κ = (L), the only maximal box is (0,L−1): ε = (−1)^L2^{−L} + (−1)^{L−1}2^{−(L−1)} = (−1)^{L+1}2^{−L}. ∎

**Corollary 9.3 (row formula).** For a ≥ 1:
  s_{a+½}(κ) = N_κ(a) + Σ_{L=1}^{κ₁} (−1)^{L+1}2^{−L} N_{κ/(L)}(a),
where κ/(L) = P_κ minus the first L boxes of row 0.
*Proof.* s_{a+½} = (A_oA_e)^a s_{½} = (A_oA_e)^{a−1}A_o ε, then Lemmas 9.1(b) and 9.2. ∎
(Example κ = (3), a = 2: 0 + ½·3 − ¼·3 + ⅛·1 = 7/8 = s_{5/2}(3).)

### 9.2 Fibre over the rows below row 0

Rows ≥ 1 of P_κ form a copy of P_{κ′} ((i,j) ↦ (i−1,j−1) preserves order and the diagonal). Row 0 is the chain
(0,0) < ⋯ < (0,n−1); its only relations to rows ≥ 1 are (0,β) < (1,β) for 1 ≤ β ≤ κ₂, and nothing in rows ≥ 1 lies below a
row-0 box. Fix g : rows ≥ 1 → [M] strict order-preserving with even diagonal (there are N_{κ′}(a) of them) and put
u_β := g(1,β) (1 ≤ β ≤ k := κ₂ ≤ n−1); u_β := M+1 (no bound) for k < β ≤ n−1. For a fixed bound vector u let
- G_L(M;u) := #{strictly increasing v_L < ⋯ < v_{n−1} in [M] with v_β < u_β}, for 1 ≤ L ≤ n (G_n := 1);
- E_0(M;u) := the same for L = 0 with v_0 even;
- Φ(u) := E_0 − Σ_{L=1}^{n}(−½)^L G_L.
Corollary 9.3 sums over g: **s_{a+½}(κ) = Σ_g Φ(u(g))**. Bounds from a real g satisfy (H): u_β ≥ β+1 for the bounded boxes 1 ≤ β ≤ k (unbounded boxes are handled in the proof of 9.5), because u_1 = g(1,1)
is even (≥ 2) and row 1 is strictly increasing.

### 9.3 The one-row lemma

Let Ψ_L(M;u) := Σ_{L′=L}^{n}(−½)^{L′−L}G_{L′}(M;u), so Ψ_n = 1 and Φ = E_0 + ½Ψ_1.

**Lemma 9.4 (shift).** Let 1 ≤ L ≤ n−1, M ≥ 1, and u_{L′} ≥ 2 for L ≤ L′ ≤ n−1. Then
  Ψ_L(M;u) = G_L(M−1;u−1) + ½Ψ_{L+1}(M−1;u−1).
*Proof.* For L′ ≥ L let a_{L′} (resp. b_{L′}) count the fillings of boxes L′..n−1 with v_{L′} = 1 (resp. ≥ 2); b_n := 1.
Then G_{L′} = a_{L′} + b_{L′} (L′ ≤ n−1). Deleting box L′ (value 1) is a bijection onto the fillings of L′+1..n−1 with all
values ≥ 2; its inverse puts 1 at box L′, which is allowed because u_{L′} ≥ 2 and M ≥ 1. So a_{L′} = b_{L′+1} for
L ≤ L′ ≤ n−1 (for L′ = n−1 both are 1). Substituting,
Ψ_L = Σ_{L′=L}^{n−1}(−½)^{L′−L}(b_{L′} + b_{L′+1}) + (−½)^{n−L}b_n = b_L + Σ_{L″=L+1}^{n} b_{L″}(−½)^{L″−L−1}(1 − ½)
    = b_L + ½Σ_{L″>L}(−½)^{L″−L−1}b_{L″}.
Finally v ↦ v−1 identifies the fillings with values in [2,M] and bounds u with fillings in [M−1] and bounds u−1:
b_{L″}(M;u) = G_{L″}(M−1;u−1) (also for L″ = n). ∎

**Lemma 9.5 (one row).** Assume (H). (a) If M ≥ n−1:
  Φ(u) = E_0(M;u) + Σ_{i=0}^{n−2} 2^{−i−1}G_{1+i}(M−1−i; u−1−i) + 2^{−n} ≥ 2^{−n}.
(b) If M ≤ n−1: Φ(u) = (−1)^{n−1−M}2^{−n}.
*Proof.* Apply Lemma 9.4 repeatedly: step i (i = 0, 1, …) rewrites Ψ_{1+i}(M−i; u−i). Its hypotheses are M−i ≥ 1 and
u_{L′} − i ≥ 2 for L′ ≥ 1+i; the second holds by (H) (u_{L′} − i ≥ L′+1−i ≥ 2), also for unbounded boxes while M−i ≥ 1.
(a) M ≥ n−1 allows i = 0..n−2, ending at Ψ_n = 1: Ψ_1 = Σ_{i=0}^{n−2}2^{−i}G_{1+i}(M−1−i;u−1−i) + 2^{−(n−1)}. Every term is
≥ 0. (b) M ≤ n−1 allows i = 0..M−1, ending at Ψ_{1+M}(0; u−M) = Σ_{L′}(−½)^{L′−1−M}G_{L′}(0) = (−½)^{n−1−M} (only G_n = 1
survives an empty alphabet). The G-terms vanish: G_{1+i}(M−1−i) fills n−1−i boxes from M−1−i < n−1−i values. So
Ψ_1 = 2^{−M}(−½)^{n−1−M}. E_0 = 0 (n values ≥ 2 need M ≥ n+1). Then Φ = ½Ψ_1. ∎

### 9.4 Proof of Theorem HZ

s_{a+½}(κ) = Σ_g Φ(u(g)) (§9.2), the bounds satisfy (H), M = 2a−1. If 2a ≤ κ₁, then M ≤ n−1 and Lemma 9.5(b) gives
Φ ≡ (−1)^{n−2a}2^{−n} = (−1)^{κ₁}2^{−κ₁}. If 2a ≥ κ₁, then M ≥ n−1 and Lemma 9.5(a) gives Φ ≥ 2^{−n}. The number of g is
N_{κ′}(a) (with N_∅ = 1), and §6 applied to κ′ says N_{κ′}(a) = 0 iff κ′ ≠ ∅ and a ≤ ⌈r′/2⌉. ∎

**Worked example (κ = (3,2), a = 2, M = 3 ≥ n−1 = 2).** Rows ≥ 1: (1,1) < (1,2), (1,1) even, values in [3]: only g = (2,3),
so u = (2,3). E_0 = 0 (v_1 < 2 is impossible above an even v_0). i = 0: G_1(2; (1,2)) = 0 (v_1 < 1). Terms: 0 + 0 + 0 + ⅛.
So s_{5/2}(3,2) = 1/8, which is the value (★) gives, and it equals 2^{−κ₁}N_{(2)}(2): the bound is attained.

### 9.5 Verification

- `prove254/halfthm.py m 5`: (−1)^{|κ|}H_κ(−a−½) computed from **Demazure characters directly** (type C, interpolated
  to x = −a−½), against Σ_g Φ(u(g)), all κ with parts ≤ m, a = 1..5. Also checks Lemma 9.5 pointwise for every g.
  m = 2: 15/15; m = 3: 35/35; m = 4: 75/75; Lemma 9.5 pointwise 477/477, 0 failures. (m = 5 Demazure run: exact-Fraction
  interpolation at degree 15, still running at write-up time — not counted.)
- `prove254/hz_vs_rec.py 6 4`: Theorem HZ's closed forms (Σ_g of Lemma 9.5(a)/(b)) against the (★)-recursion values, all κ
  with parts ≤ 6, a = 1..4: 252/252. (Checks the algebra of §9, not the inputs.)
- `prove254/rowform.py 5 5`: row formula (Cor 9.3) = (★)-recursion values, 155/155.
- `prove254/rowphi*.py`: Lemma 9.5 needs (H). Without it Φ goes negative (u = (1): −¼; u = (2,2), weakly increasing:
  −⅛), so the hypothesis is load-bearing and the test separates. With M ≤ n−2 Φ < 0 is the predicted sign (−1)^{n−1−M}.
- Data (half.py): 2^{κ₁}s_{a+½}(κ) is always an integer, e.g. (5,4) at a = 3 gives 1, the minimum 2^{−κ₁}N_{κ′}.

### 9.6 What this closes

- (NZ) at half-integers: all κ, all a, with strictness. PROVE 253's open base inequality (odd μ₁) is not needed.
- Prop. R (odd-zero-law file) gets a second, residue-free proof (§9.4 small-a case).
- (Z) gets a second proof (zero iff N_{κ′}(a) = 0), independent of Theorem Z's induction.
- A by-product: T := (I+C)^{1/2} (so that s_{x+½} = T s_x) has entries 0, ±2^{−k}: it is multiplicative over connected
  components of θ, an off-diagonal singleton has weight 1, an off-diagonal component of size ≥ 2 has weight 0, and a
  component containing exactly one diagonal box is a ribbon D with weight ±2^{−|D|} (sign + iff its NE end box has a box to
  its left). Computed m ≤ 6 (1716/1716 pairs, `tprime.py`). Not proved, not needed.

## 10. Gaps (updated)

Theorem HZ adds only Thm A (proved given (G1), odd-zero-law file) to the inputs of §8. §9 is self-contained otherwise.
Lemmas 9.4/9.5 were written once and re-read once this session; nobody else has refereed them.
