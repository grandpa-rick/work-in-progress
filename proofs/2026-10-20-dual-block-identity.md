# The consecutive dual-block identity: the alternant form of the twisted family (Day 248 PROVE)

Rick, 2026-10-20. Continues `proofs/2026-10-19-theorem-W-upper-end-certifier.md` (notation: T⁰, Ψ_{κ,θ}, ρ, Lemma T,
label lowering, the sum certifier). Scripts: `work-in-progress/scripts/prove248/`.

## 0. Verdict

| item | statement | grade |
|---|---|---|
| **Lemma A** (alternant form) | Ψ_{κ,θ}(M) = Σ_{H⊆M, \|H\|=r} S₀^M(H) · Alt[Π_j x_j^{j−1+θ_j} ρ_{κ_j}(x_j)] / V(x_H). Every twisted Ψ is a *symmetric* function of its head set. | **proved** (§1) |
| **Cor. A1** (consecutive closed form) | κ = (s,…,s+r−1), twist ≡ t: Ψ_{κ,t}(M) = Σ_H S₀^M(H) Π_{j∈H} x_j^t ρ_s(x_j). | **proved** (§1) |
| **Target (a)/(c): B_L closed form, all L** | The certifier's dual block sum B_L (L=2: ν^c=(s,s+1), twists (1,1),(1,0),(0,0), signs +,−,+) **is the single instance Ψ_{(s,…,s+L−1),(1,…,1)}**. Two of the three L=2 summands cancel *identically*: Ψ_{(s,s+1),(1,0)} = w_s·Ψ_{(s,s+1),(0,0)}. | **proved** (§2) |
| **Lemma H** (exact peel) | Ψ_{κ,θ}(M) = tail part + (1−u^{−2}) g_{M⁻}(u) u^{−(r−1)} Σ_{H'} S₀^{M⁻}(H') Π_{j∈H'}(u−x_j)·A_{κ,θ}(u,x_{H'}), exact. Lemma T is its termwise leading order. | **proved** (§3) |
| **Lemma P** (top-block collapse) | For the top consecutive block {p..r} of κ with constant twist t there, the |B| tied Lemma-T head options sum to **one** option: head p's order, sign and scalar, with child Ψ_{(κ_1..κ_{r−1}), (θ_1..θ_{p−1}, t+1,…,t+1)}. | **proved** (§3) + exact check (`t5.py`, 6/6 ratio 1) |
| **Target (b): the gap-2 cancellation, explained** | It is not a deep cancellation. The certifier's 3-summand form of B_2 contains a formal zero, and its Lemma-T options cancel exactly. The true u→0 order of B_2 is t−s−r = 1−s−2 (one above termwise), with coefficient (−1)^{s+1} W_{<s} e_{M⁻}^{−2} Ψ_{(s),(2)}(M⁻). | **proved** (§2–3) |
| **Lemma P1** (tie, s=1) | Ψ_{[r],0}(M) + e_M^{−2} Ψ_{[r−1],2}(M) = Ψ_{[r−1],0}(M). This is the r-variable (1)-identity, and it resolves the head/tail tie at κ=[r], t=1. | **proved** (§4, Lagrange + one residue) + sympy r=2,3 |
| **Lemma P1′** (tie, r=1, t=s≥2) | [u^{−1}] Ψ_{(s),s}(M) = −Σ_{i=0}^{s−1} (−1)^i e_i(w_{<s}) γ_{s−1−i}(M⁻), γ_k = Ψ_{(1),(k)} (k≥1), γ_0 = 1. | **proved** (§4) + sympy 3/3 |
| **Prediction** | The ties occur exactly for the blocks (s,…,2s) (the twist reaches t=s at r=1), and (1,…,r) via P1. The gap-2 C-nodes ν^c=(1,2,3) need P1. The C₆ survivors ν^c=(2,3,4) are the block (s,…,2s) with s=2 and need P1′. | confirmed (§5) |
| **Certifier v2** (`scert5.py` = scert + P + P1 + P1′, fallback to v1); v3/v4 add P2 and P\* | C n≤6: **0 uncertified** (was 7 + 208). B n≤6: 0. A N=6: 0/68 791. **A N=7: 0/789 998 (was 325) ⇒ Z^A_7 = [0,10]**. Negative control: 0 false certifications, and mutants of every new branch are caught (§5.2). | computed |

The headline: the "multi-summand obstruction" was a *representation* artefact. Twisted Ψ's with the same labels satisfy
linear relations, because only their alternant matters, and the certifier did not know these relations. In the
alternant picture a consecutive block with constant twist is a single product, and the u→0 peel of it is a single
child. The genuine new phenomenon is the head/tail **tie**. It happens only at t = κ_p − p + 1, and Lemmas P1/P1′
resolve the two cases that occur through N=7.

## 1. The alternant form

Notation. M ⊆ [N]; for H ⊆ M write x_H for the variables w_m, m ∈ H. Put ρ_κ(x) := Π_{k<κ}(1 − w_k/x) and
  S₀^M(H) := Π_{j∈H}(1 − x_j^{−2}) · Π_{{j,k}⊆H}(1 − 1/(x_j x_k)) · Π_{j∈H, t∈M∖H} (x_j − 1/w_t)/(x_j − w_t),
which is symmetric in x_H. V(x_1..x_r) := Π_{j<k}(x_k − x_j), and Alt[F](x) := Σ_{σ∈S_r} sgn(σ) F(x_{σ(1)},…,x_{σ(r)}).
For symmetric f in r variables put Φ_r^M[f] := Σ_{|H|=r} S₀^M(H) f(x_H).

**Lemma A.** For labels κ_1<…<κ_r and any θ ∈ ℤ^r,
  Ψ_{κ,θ}(M) = Φ_r^M[A_{κ,θ}],  A_{κ,θ}(x) := Alt[Π_j x_j^{j−1+θ_j} ρ_{κ_j}(x_j)] / V(x).

*Proof.* Fix an injective h with image H, and put x_j = w_{h_j}. In T⁰_κ(h) (d-infinity file §1):
- the head factors Π_j(1 − x_j^{−2}) are symmetric in x_H;
- (1 − w_t^{−1}x_j^{−1})/(1 − w_t/x_j) = (x_j − 1/w_t)/(x_j − w_t);
- each pair factor (1 − 1/(x_jx_k))/(1 − x_j/x_k) (j<k) equals (1 − 1/(x_jx_k))·x_k/(x_k − x_j).

Hence T⁰_κ(h)Π_j x_j^{θ_j} = S₀^M(H) · Π_j x_j^{θ_j}ρ_{κ_j}(x_j) · Π_{j<k} x_k /V(x) = S₀^M(H)·F(x)/V(x), where
F = Π_j x_j^{j−1+θ_j}ρ_{κ_j}(x_j). Summing over the r! orderings of H, and using V(x_σ) = sgn(σ)V(x), gives
S₀^M(H)·Alt[F]/V. The orderings with some h_j < κ_j are included, but they contribute 0 because ρ_{κ_j}(w_{h_j})
has the factor (1 − w_{h_j}/w_{h_j}). ∎

Alt[F]/V is a symmetric Laurent polynomial (an alternant over the Vandermonde), so Φ_r^M[·] only sees A. Label lowering
(the certifier file §2) is the polynomial identity ρ_{κ'} = ρ_κ · Π_{κ≤k<κ'}(1 − w_k/x) inside F, and nothing more.

**Corollary A1 (consecutive block).** If κ = (s, s+1, …, s+r−1) and θ ≡ t, then A_{κ,t}(x) = Π_j x_j^t ρ_s(x_j), so
  Ψ_{κ,t}(M) = Σ_{|H|=r} S₀^M(H) Π_{j∈H} x_j^t ρ_s(x_j).
*Proof.* With a_i := w_{s+i−1}, x^{j−1}ρ_{s+j−1}(x) = ρ_s(x) Π_{i<j}(x − a_i) = ρ_s(x)(x|a)^{j−1}, where (x|a)^m is monic of
degree m in x. So Alt[F] = Π_j x_j^tρ_s(x_j) · det[(x_i|a)^{j−1}]. Since (x|a)^{j−1} = x^{j−1} + (lower powers), unipotent
column operations turn det[(x_i|a)^{j−1}] into det[x_i^{j−1}] = V(x). ∎

(For s = 1, ρ_1 = 1: Ψ_{[r],0}(M) = Σ_H S₀^M(H). This explains why Ψ_{(1,2),0}((2,3,4)) has no w_1 (`t3.py`).)

## 2. The dual block sums B_L (targets (a) and (c))

**Where B_L comes from.** The certifier peels a dual Ψ_{κ,0} with κ = (s,…,s+L) consecutive at u = w_n → 0. By Lemma T
all L+1 heads tie at order −s−L−1. Their children Ψ_{κ∖κ_j, θ^{(j)}} (θ^{(j)} = (0^{j−1},1^{L+1−j})) carry the scalars
W_{<κ_j} = W_{<s}a_1⋯a_{j−1}. After label lowering to (s,…,s+L−1), the node is
  B_L = Φ_L[ Σ_{j=1}^{L+1} a_1⋯a_{j−1} Π_{i=j}^{L}(x_i − a_i) · (alternant weight) ].
For L = 2 this is Ψ_{(s,s+1),(1,1)} − a_2Ψ_{(s,s+1),(1,0)} + a_1a_2Ψ_{(s,s+1),(0,0)}: the twists (1,1),(1,0),(0,0) with
signs +,−,+ of `instr2.py`.

**Theorem B (closed form, all L).** B_L = Ψ_{(s,…,s+L−1),(1,…,1)}, a single instance. Equivalently, up to the global factor
(−1)^{s+L} W_{<s} e_{M⁻}^{−2}, the u→0 leading coefficient of Ψ_{(s..s+L),0}(M) is Ψ_{(s..s+L−1),(1^L)}(M⁻).
For L = 2 the identity behind it is
  Ψ_{(s,s+1),(1,0)} = w_s · Ψ_{(s,s+1),(0,0)}.
*Proof.* The second form is Lemma P (§3) with p = 1, t = 0. For the identity, Lemma A with ω := x_2 − a_1 (relative label
factor) gives A_{(1,0)} = ρ_s(x_1)ρ_s(x_2)·Alt[x_1(x_2 − a_1)]/V. Here Alt[x_1x_2] = 0, and Alt[−a_1x_1] = −a_1(x_1 − x_2) = a_1V.
So A_{(1,0)} = a_1ρρ = a_1A_{(0,0)}. Then −a_2Ψ_{(1,0)} + a_1a_2Ψ_{(0,0)} = 0. ∎
Check: `t2.py` (s=2, M=(2,3,4)): the exact u→0 leading coefficient of Ψ_{(2,3,4),0}(M ∪ 5) equals the formula, ratio W_{<2}e_{M}^{−2}.

**The gap-2 cancellation (target (b)).** At the next peel u′ = max M⁻, Lemma T applied to the three summands gives
head-2 of (1,0), and heads 1 and 2 of (0,0), at order −s−2. Their coefficients are
−W_{<s}a_1a_2Ψ_{(s),(1)}, W_{<s}a_1a_2(Ψ_{(s),(1)} − a_1Ψ_{(s),(0)}), and W_{<s}a_1²a_2Ψ_{(s),(0)}, which sum to 0. This is
simply the formal zero −a_2Ψ_{(1,0)} + a_1a_2Ψ_{(0,0)} being peeled. Applied to the single instance Ψ_{(s,s+1),(1,1)}, Lemma P
gives the true order 1 − s − 2 = −s−1. The coefficient is (−1)^{s+1}W_{<s}e^{−2}Ψ_{(s),(2)}, a single summand. Exact check: `t4.py`,
where the heads of Φ_{(s,s+1)}[Q] are C(u)·u·Φ_{(s)}[−u x_1(x_1 − u)], identically (0 = sympy difference).

## 3. The exact peel and the top-block collapse

**Lemma H (exact peel).** n ∈ M, u = w_n, M⁻ = M∖n, κ ⊆ [n−1], g_S(u) := Π_{t∈S}(u − 1/w_t)/(u − w_t). For symmetric f,
  Φ_r^M[f] = Φ_r^{M⁻}[ f(x) Π_j (ux_j − 1)/(u(x_j − u)) ] + (1 − u^{−2}) g_{M⁻}(u) u^{−(r−1)} Φ_{r−1}^{M⁻}[ f(u, x') Π_{j}(u − x'_j) ].
*Proof.* Split H by whether n ∈ H.
- *H ∌ n:* S₀^M(H) = S₀^{M⁻}(H)·Π_{j∈H}(x_j − 1/u)/(x_j − u), and (x − 1/u)/(x − u) = (ux − 1)/(u(x − u)).
- *H = H′ ∪ n:* S₀^M(H) = S₀^{M⁻}(H′) · (1 − u^{−2}) · Π_{j∈H′}(1 − 1/(ux_j)) · Π_{t∈M⁻∖H′}(u − 1/w_t)/(u − w_t). The tail
  sets agree, since M∖H = M⁻∖H′. Write the last product as g_{M⁻}(u)/Π_{j∈H′}[(u − 1/x_j)/(u − x_j)], and use
  (1 − 1/(ux))(u − x)/(u − 1/x) = (u − x)/u. ∎

With f = A_{κ,θ}, Laplace-expanding Alt[F](u, x′) = det[F_j(y_i)] along the u-row, and V(u,x′) = Π_k(x′_k − u)V(x′), give
  A(u,x′)Π_j(u − x′_j) = Σ_j (−1)^{r−j} F_j(u) A_{κ∖κ_j, θ^{(j)}}(x′),   F_j(y) = y^{j−1+θ_j}ρ_{κ_j}(y).   (H′)
This holds because deleting column j shifts the later columns one position left, which is the twist +1 of θ^{(j)}. So Lemma H is
exactly Lemma T with exact u-prefactors s_j(u) = (−1)^{r−j}u^{θ_j−(r−j)}(1 − u^{−2})ρ_{κ_j}(u)g_{M⁻}(u). Lemma T keeps the leading
terms of these, and the tail's leading term is (−1)^r u^{−r}Ψ_{κ,θ−1}(M⁻). Numerical identity check: `t4.py`, three (κ, P, M),
difference 0. *This is the cold re-derivation of Lemma T owed from Day 247b (§2 of the certifier file): it survives, and is
now exact.*

**Lemma P (top-block collapse).** Let B = {p,…,r} be the top consecutive block of κ (κ_j = κ_p + j − p on B, and p = 1 or
κ_{p−1} < κ_p − 1). Suppose θ_j = t for all j ∈ B, and κ_r < n. Then, as u → 0, the head part of Ψ_{κ,θ}(M) equals
  Σ_{j<p} s_j(u) Ψ_{κ∖κ_j,θ^{(j)}}(M⁻)
  + (−1)^{κ_p+r−p} W_{<κ_p} e_{M⁻}^{−2} u^{t−κ_p−1−(r−p)} ( Ψ_{κ′,θ′}(M⁻) + O(u) ),
with κ′ = (κ_1,…,κ_{r−1}) and θ′ = (θ_1,…,θ_{p−1}, t+1,…,t+1).
*Proof.* On B, F_j(y) = C(y)π_{j−p}(y), where C(y) = y^{p−1+t}ρ_{κ_p}(y), π_m(y) = Π_{i=0}^{m−1}(y − w_{κ_p+i}) is monic of degree m.
- In det[F_j(y_i)], unitriangular column operations inside B (the same ones for every row) replace C·π_m by C·y^m.
- Laplace along the u-row: the block columns contribute Σ_{m=0}^{|B|−1} ±C(u)u^m · minor_m. The m = 0 term has the lowest order,
  and the m ≥ 1 terms are higher by m (exactly, not just termwise).
- minor_0 keeps the block columns C(x)x^m, m = 1..|B|−1. Undoing the column operations on these, x·C(x)·x^{m−1} ↦ x·C(x)π_{m−1}(x),
  gives V(x′)·A_{κ′,θ′}(x′). The reason: x·C(x)π_{m′}(x) = x^{(p+m′)−1+(t+1)} ρ_{κ_p+m′}(x).
- Leading factors: (1 − u^{−2}) → −u^{−2}, g_{M⁻}(0) = e_{M⁻}^{−2}, ρ_{κ_p}(u) → (−1)^{κ_p−1}W_{<κ_p}u^{1−κ_p}. The sign of the
  m = 0 cofactor matches Lemma T's head p, which gives the stated order, sign and scalar. ∎

So Lemma P replaces the |B| tied Lemma-T options by **one option**: head p's data, but with the block labels lowered by one
and twist t+1 there. The tail (−1)^r u^{−r}(Ψ_{κ,θ−1}(M⁻)+O(u)) ties with it iff t = κ_p − p + 1. The heads j < p have
strictly higher order when θ is constant overall, because κ_j − j < κ_p − p.
Exact check (`t5.py`, p = 1): (s,r,t) ∈ {(2,3,0),(1,2,0),(2,2,1),(3,2,0),(2,2,−1),(1,3,0)}: computed order = predicted, and the
coefficient ratio is exactly 1, in all 6.

**Lemma P\* (every block).** The proof of Lemma P never used that B is the *top* block. Decompose κ into maximal consecutive
blocks B_b = {p_b,…,q_b}. In det[F_j(y_i)], perform the unitriangular column operations inside every block on which θ is
constant (≡ t_b). Then Laplace-expand along the u-row. For each such block the lowest u-power comes from its first column alone,
with order t_b − κ_{p_b} − 1 − (r − p_b), and the scalar and sign of Lemma T's head p_b. Its child is
  labels κ ∖ {κ_{q_b}} (remove the block's *top* label);
  twists θ_i for i < p_b, t_b + 1 on the remaining block positions, θ_{i+1} + 1 above.
Blocks with non-constant twist keep their per-head Lemma-T options (still exact, individually). Exact check (`t10.py`, random
rational point): κ = (2,3,6), θ = (0,0,4), M = (2,3,4,6,7) and κ = (1,2,4), θ = (0,0,2), M = (1,2,3,4,5). In both, the lower
block (2,3) or (1,2) is strictly extremal at order −5 resp. −4. The coefficient matches (difference 0), and the u^{o−1} coefficient
is 0. *Invariant:* starting from dual twist 0, every rule (tail −1, block collapse, forced/n∉ν removal, P1, P2) keeps the dual
twist constant on each maximal consecutive block of the dual labels. Removing labels only splits blocks, never merges them. So
on single-instance nodes **only Lemma P\*/P1/P2 options ever occur at u→0 on the dual side**.

**Iteration (all L).** Start from twist 0 on a consecutive block (s,…,s+r−1). Lemma P gives the chain
(s..s+r−1; 0) → (s..s+r−2; 1) → … → (s..s+r−1−k; k), a single summand at every step. The tail ties with the block at the
first k with k = s, which requires r − s ≥ 1. At the tie the block is (s,…,r−1), of length r−s. Hence:
- s = 1: tie after one step, block (1..r−1), t = 1: **Lemma P1** (any length);
- r = s+1, s ≥ 2: tie at the single label (s) with t = s: **Lemma P1′**;
- s ≥ 2, r ≥ s+2 (e.g. ν^c ⊇ (2,3,4,5)): tie of a block of length ≥ 2 at t = s: **Lemma P2** (below; found and proved later in the session).

## 4. The ties

**Lemma P1.** For all M and r ≥ 1: Ψ_{[r],0}(M) + e_M^{−2}Ψ_{[r−1],2}(M) = Ψ_{[r−1],0}(M), where e_M := Π_{m∈M}w_m.
(r = 1 is the (1)-identity φ⁰_{(1)} = 1 − e_M^{−2}.)

*Proof.* By Cor. A1, the claim reads Φ_r[1] = Φ_{r−1}[1] − e_M^{−2}Φ_{r−1}[e_{r−1}²].
1. Fix H′ with |H′| = r−1 and m ∉ H′. Then S₀^M(H′ ∪ m) = S₀^M(H′) · (1 − x_m^{−2}) Π_{j∈H′}(x_j − x_m)/x_j · Π_{t∈M∖(H′∪m)}(x_m − 1/w_t)/(x_m − w_t).
   This uses (1 − 1/(x_jx_m))/(x_j − 1/x_m) = 1/x_j to trade the tail factor of m in S₀^M(H′) for the pair factor.
2. So, with R(z) := z^{−1}G_{M∖H′}(z)Π_{j∈H′}(1 − z/x_j), G_S(z) := Π_{t∈S}(z − 1/w_t)/(z − w_t): S₀^M(H′∪m) = S₀^M(H′)·Res_{z=x_m}R(z).
3. *Lagrange:* for any H with |H| = r, Σ_{m∈H}Π_{j∈H∖m} x_j/(x_j − x_m) = 1. This is the interpolation of the constant 1 at
   the nodes x_H, evaluated at 0. Hence
   Φ_r[1] = Σ_H S₀(H)·Σ_{m∈H}Π_{j≠m}x_j/(x_j − x_m) = Σ_{H′} S₀^M(H′) Σ_{m∉H′} Res_{z=x_m}[R(z)/Π_{j∈H′}(1 − z/x_j)]
          = Σ_{H′}S₀^M(H′) Σ_{m∉H′} Res_{z=x_m}[z^{−1}G_{M∖H′}(z)].
4. The function z^{−1}G_{M∖H′}(z) has poles only at w_t (t ∈ M∖H′), 0 and ∞. Res_0 = G_{M∖H′}(0) = Π_{t∈M∖H′}w_t^{−2} and Res_∞ = −1.
   So the inner sum is 1 − e_M^{−2}e(x_{H′})², and summing over H′ gives Φ_{r−1}[1] − e_M^{−2}Φ_{r−1}[e²]. ∎

The same argument with any f gives **Lemma Λ**: Φ_r^M[f] = −Σ_{H′}S₀^M(H′)(Res_0 + Res_∞)[z^{−1}G_{M∖H′}(z)f(x_{H′},z)]
(numerical check `t7.py`, three f, difference 0). Its z⁰-part gives the generalised P1: Φ_r[Lag F] = Φ_{r−1}[F] − e_M^{−2}Φ_{r−1}[e²F].
Here Lag F(x_H) := Σ_{m∈H}F(x_{H∖m})Π_{j≠m}x_j/(x_j − x_m).

**Use of P1.** At κ = [r] with t = 1 (reached by one Lemma-P step from (1..r+1; 0)), the tie p = 1, t = κ_1 = 1 holds.
- Lemma P's option is (−1)^r e_{M⁻}^{−2}Ψ_{[r−1],2}(M⁻) u^{−r}, and the tail is (−1)^rΨ_{[r],0}(M⁻)u^{−r}.
- By P1 the sum of the two is **(−1)^r u^{−r}Ψ_{[r−1],0}(M⁻)**, one summand.

**Lemma P1′.** For s ≥ 2: Ψ_{(s),s}(M) = −u^{−1}( Σ_{i=0}^{s−1}(−1)^i e_i(w_1,…,w_{s−1}) γ_{s−1−i}(M⁻) + O(u) ), where γ_0 := 1
and γ_k(S) := Ψ_{(1),(k)}(S) for k ≥ 1.
*Proof.* Ψ_{(s),s}(M) = Φ_1^M[f] with f(x) = x·(x|w)^{s−1} (Cor. A1). Also Φ_1^M[f] = Σ_m Res_{z=x_m}[z^{−1}G_M(z)f(z)]: the
residue of G_M at x_m is (x_m − 1/x_m)Π_{t≠m}(…).
- z^{−1}f(z) is a polynomial, so Res_0 = 0, and Φ_1^M[f] = [z^{−1}]_∞(G_M(z)(z|w)^{s−1}).
- The same computation shows γ_k(S) = [z^{−k}]_∞G_S(z) = Ψ_{(1),(k)}(S) for k ≥ 1.
- Write G_M = G_{M⁻}·(1 + (u − u^{−1})Σ_{l≥0}u^l z^{−l−1}). Then
  [z^{−1}](G_M(z|w)^{s−1}) = [z^{−1}](G_{M⁻}(z|w)^{s−1}) + (u − u^{−1})Σ_{l≥0}u^l[z^l](G_{M⁻}(z|w)^{s−1}).
- The u^{−1} part is −[z⁰](G_{M⁻}(z|w)^{s−1}) = −Σ_i(−1)^ie_iγ_{s−1−i}(M⁻). Everything else is O(1). ∎
Exact check (`t8.py`): s = 2, M = (2,3,4); s = 3, M = (3,4,5); s = 2, M = (1,2,4): leading order −1, difference 0.
(s = 1 recovers the (1)-identity: [u^{−1}] = −1.)

**Lemma P2 (the general tie, whole block, s ≥ 2).** Let κ = (s,…,s+m−1), s ≥ 2, m ≥ 1, κ_m < n, and twist ≡ s. Then
Ψ_{κ,s}(M) = u^{−m}(X + O(u)), where
  X = (−1)^m [ Σ_{i=0}^{s−2} (−1)^i Σ_{T⊆[s−1], |T|=i} w_T · Ψ_{(1,s,…,s+m−2),(s−1−i, s−1,…,s−1)}(M⁻)
             + (−1)^{s−1} W_{<s} · Ψ_{(s,…,s+m−2),(s−1,…,s−1)}(M⁻) ].
(m = 1 is P1′. For s = 1 the i-sum is empty and the formula is P1.)

*Proof.*
1. By Lemma P and the tail, X = (−1)^m[Ψ_{κ,s−1}(M⁻) + (−1)^{s−1}W_{<s}e_{M⁻}^{−2}Ψ_{κ⁻,s+1}(M⁻)].
2. By Cor. A1, Ψ_{κ,s−1} = Φ_m[Π_j(x_j|w)^{s−1}] with (x|w)^{s−1} := Π_{k<s}(x − w_k) = x^{s−1}ρ_s(x).
3. Apply Lemma Λ with f(x′,z) = G′(x′)(z|w)^{s−1}, G′ := Π_{j∈H′}(x_j|w)^{s−1}.
   - The Res_0 part is −(−1)^{s−1}W_{<s}e_{M⁻}^{−2}Φ_{m−1}[e²G′] = −(−1)^{s−1}W_{<s}e^{−2}Ψ_{κ⁻,s+1}. It **cancels** the second term of 1.
   - What is left is (−1)^m Σ_{H′}S₀(H′)G′(x′)·[z⁰]_∞(G_{M⁻∖H′}(z)(z|w)^{s−1})
     = (−1)^m Σ_{H′}S₀(H′)G′ Σ_{i=0}^{s−1}(−1)^ie_i(w_{<s})γ_{s−1−i}(M⁻∖H′).
4. The γ_0 = 1 term is (−1)^{s−1}W_{<s}Φ_{m−1}[G′] = (−1)^{s−1}W_{<s}Ψ_{κ⁻,s−1}.
5. For k ≥ 1, γ_k(S) = Σ_{t∈S}S₀^S({t})x_t^k. The Lagrange ratio S₀^{M⁻}(H′)S₀^{M⁻∖H′}({t}) = S₀^{M⁻}(H′∪t)Π_{j∈H′}x_j/(x_j − x_t)
   (step 1 of P1) gives Σ_{H′}S₀(H′)G′γ_k(M⁻∖H′) = Φ_m[Lag_kG′], where Lag_kG′(x_H) := Σ_{t∈H}x_t^kG′(x_{H∖t})Π_{j≠t}x_j/(x_j − x_t).
6. Lag_kG′ = A_{(1,s,…,s+m−2),(k,s−1,…,s−1)}. To see it, take det[y_i^k | (y_i|w)^{s−1}y_i^{j−1} (j = 2..m)]/V.
   - Laplace along column 1, with V(x_{≠i})/V(x) = (−1)^{i−1}/Π_{l≠i}(x_l − x_i), gives exactly Lag_kG′.
   - Unitriangular operations in columns 2..m turn y^{j−2} into Π_{i=0}^{j−3}(y − w_{s+i}). Column j becomes
     y^{(j−1)+(s−1)}ρ_{s+j−2}(y), and column 1 is y^{0+k}ρ_1(y). ∎

Exact check (`t9.py`): (s,m,M) = (2,2,(2,3,4)), (2,2,(1,2,3,5)), (3,2,(3,4,5)), (2,3,(2,3,4,5)): order −m, difference 0.
With P1, P1′ and P2, **every head/tail tie of a whole consecutive block with constant twist is resolved exactly**, as a
finite sum of twisted instances.

## 5. Certifier v2 (v3, v4) and the N = 7 prediction

`scert5.py` = `scert.py` + Lemma P (any top block, constant twist, at u→0, every side) + P1 + P1′, with fallback to v1 per root.
- *Soundness:* each new option is the exact leading coefficient of the sum it replaces (§3–4). The §3 soundness argument of the
  certifier file is unchanged.
- *Results:*

| run | v1 (Day 247b) | v2 |
|---|---|---|
| C n=5 / n=6 uncertified | 7 / 208 | **0 / 0** (P alone: 4 / 106; P+P1: 0 / 4, the 4 being ν^c=(2,3,4) as predicted, cleared by P1′) |
| B n ≤ 6 uncertified | 0 | 0 |
| A N=6 | 0 / 68 791 | 0 / 68 791 (P fired 7566×, P1 312×, P1′ 98×) |
| A N=7 | 325 / 789 998 | **0 / 789 998** (§5.1) |
| negative control A N=4,5 all c ∈ [−1,2N]; B,C n ≤ 4 | 0 false | **0 false** (and more in-window nonzeros certified, e.g. N=4 c=1: 405/405) |

### 5.1 N = 7
`runA5.py 7` (v2 = P + P1 + P1′, fallback v1): **789 998 allowed A-triples, 0 uncertified** (1103 s; P fired 146 706×, P1 8 773×,
P1′ 3 483×). The 325 v1 failures (ν ∈ {(4,5,6,7), (1,5,6,7), (1,2,6,7)}, i.e. ν^c = (1,2,3), (2,3,4), (3,4,5)) are all certified,
as §3 predicts:
- (1,2,3) runs P → P1;
- (2,3,4) is the block (s,…,2s) with s = 2 and runs P → P → P1′;
- (3,4,5) has r < s+1, so it never ties, and P alone suffices.
This agrees with `chk7.py` (all 325 nonzero at c = 11, 12, 13).
Together with Prop Z (lower inclusion, proved) and W⁺ (c ≤ −1, proved): **Z^A_7 = [0, 10] = [0, 2N−4]** (computer-assisted
certificate).

### 5.2 Mutation test (does the negative control see the new branches?)

Each new branch was deliberately broken, and the negative control (known in-window zeros, fixed-c mode) was rerun:

| mutant | zero set | false certifications |
|---|---|---|
| Lemma P child twist t (instead of t+1) | A N=4 + B,C n=3,4 (1 774 zeros) | **129** |
| P1 child twist 1 (instead of 0) | same | **7** |
| P1′ drops the γ₀ term | C n=5, c ∈ [−1,10] (3 493 zeros) | **3** |
| P2 last-term sign flipped | C n=6 with {2,3,4,5} ⊆ ν^c, c ∈ [−1,12] (71 zeros; P2 fires 19×) | **3** |
| (unmutated v2/v3/v4) | all of the above + A N=5 (1 226) | **0** |

So the control has teeth on every new branch.

## 6. Target (c): what the all-N upper end now reduces to

**Setting.** Fix X ∈ {A,B,C}, N, an allowed root and c > U^X_N. The certifier peels w_N, w_{N−1}, … one at a time.
- At each step it picks a limit (u→0 or u→∞) and keeps the extremal options.
- All options are exact leading coefficients (Lemma T/H, P\*, P1, P2, the certifier file's §2–3).
- So the base signed count σ is the coefficient of one specific monomial of X in an iterated (lexicographic) Laurent expansion.
  Hence X ≠ 0 whenever σ ≠ 0. This is the soundness theorem, unchanged.

**What is now settled for all N (proved).** On the dual side at u→0, every node that is a single twisted instance with twist
constant on each consecutive block has a *complete* exact option list: one option per block (P\*), the tail, and an exact
resolution of the whole-block tie (P1/P2). By the invariant in §3 this class contains the whole u→0 dual orbit of every root.
The consecutive dual block is therefore **no longer an obstruction in any N**.

**What remains (precise).** The upper end for all N follows from:

> **(SC) Single-coefficient conjecture.** For every allowed root with c > U^X_N there is a peel policy for which the base signed
> count is σ = ±1.

Evidence: σ ∈ {±1} at **every** base reached, over all certified roots of A N ≤ 5 and C n ≤ 5 (`base.py`: A₅ 2736×(−1), 2863×(+1);
C₅ 2434×(−1), 2439×(+1); no other value occurs), and 99.9% of internal nodes are single-summand (A₅: 33 524 of 33 560).
The remaining sources of multi-summand nodes, i.e. the only places where a cancellation σ = 0 could hide, are exactly:
- (α) tail/block ties when κ has ≥ 2 blocks (t_b = κ_{p_b} − p_b + 1 for a non-bottom block). No identity is proved; the two
  options are carried as a formal sum, which is sound but not simplified;
- (β) ties on the λ/μ (F) side in in-terms;
- (γ) in/out ties, which occur only at the finitely many c where the certifier's 'eq' split fires (c ≤ U + small);
- (δ) the 2^{s−1}-term outputs of P2.

A uniform proof of (SC) needs, for each of (α)–(δ), either an identity collapsing the sum or a strict-order argument showing a
unique extremal summand at the next step. For (δ) the next-step orders are separated by the twist k = s−1−i on label 1,
so (δ) is plausibly a strict-order argument. (α) is the real open identity.

## 7. Gaps (precise)

1. **All N is still open, but the obstruction moved.** The whole-block dual case is closed for all N (P\*, P1, P2). Two things
   remain. (i) The tie identity for multi-block κ (§6 (α)). (ii) A proof of (SC). Both are stated precisely in §6.
2. **Code trust (v2).** The new branches are short and each one transcribes a lemma proved here. The negative control's power
   on the new branches is reported in §5.2 (mutation test).
3. Inherited: (D) and everything downstream carry the Motegi root (`motegi-rtv-biorthogonality-root`).
