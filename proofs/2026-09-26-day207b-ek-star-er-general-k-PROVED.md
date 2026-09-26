# Day 207b: e_k⋆e_r Pieri rule PROVED for all k ≥ 1, r ≥ 0

**Date:** 2026-09-26 (deep-work session)
**Author:** Rick
**Status:** PROVED as an operator identity on Λ_m for every m ≥ 0, k ≥ 1, r ≥ 0. It covers:
- **(A_k)**, the parabolic kernel;
- **(K_k)**, the k-set kernel;
- **(E_k)**, the generating function, together with its e-basis Pieri corollary.

The ⋆-reading uses Hikita 2503.23597 Def 3.4 / Lemma 3.3, exactly as in Day 206b (R0, verified-quote). No other literature is load-bearing. The standard Coxeter/Hecke facts in §1 are textbook (Humphreys, *Reflection groups and Coxeter groups* §§1.10 and 7.1–7.3; Björner–Brenti Prop. 2.4.4).

This supersedes the `computed` grades of Day 207, namely the k = 3, 4, 5 kernels, the telescoping GF, and the general M-formula. The PROVE.md "minimum deliverable" (A_3 + K_3) and the stretch goal (A_k, K_k, E_k for all k) are both done.

**Scripts** (`scripts/day207b/`, all logs ALL OK):
- `check_general_k.py`: checks KL, A_k per tuple, OL, (★), T^{(0)}_k = 0, M ≡ c, and the Pieri formula vs direct AHA.
  - KL and A_k: exact in ℤ[X,s,t] for m = 4, 5, 6 and k ≤ 4.
  - OL: 5 (m,k) cases.
  - (★): symbolic for n ≤ 9.
  - Pieri vs AHA: 216 cases, m ≤ 6, k ≤ 4, all r ≤ m, *including m < k*.
- `check_recursion_step.py`: the E_k induction step on closed forms, exact, z symbolic, m ≤ 5, k ≤ 4.
- `check_specializations.py`:
  - k = 1 is Thm 3.12, k = 2 is W_r, and k = 3, 4 are Day 193/195, all symbolically in u = t^r.
  - The r = 0 normalization e_k(Y)•1 = t^{C(k,2)}e_k holds for k ≤ 7.

## 0. Conventions and statement

Conventions as in Day 206b §0 (`scripts/day198/p2Y_er.py`, `scripts/day205/k3_fast_pipeline.py`):
- T_iF = t s_iF + (t−1)X_{i+1}(s_iF − F)/(X_i − X_{i+1}), which satisfies (T_i − t)(T_i + 1) = 0;
- πF = X_1F(X_2..X_m, sX_1), with s = q^{−1};
- Y_i = t^{m−i}T_{i−1}⋯T_1 π T_{m−1}^{−1}⋯T_i^{−1};
- e_k(Y) = Σ_{b_1<⋯<b_k} Y_{b_1}⋯Y_{b_k}. The Y_i commute, so the order is irrelevant; the proof uses only this order.

Further notation:
- a_{ij} = (X_i − tX_j)/(X_i − X_j);
- for A ⊆ [m], ∏^×_A := ∏_{i∈A, j∉A} a_{ij};
- E(z) = ∏(1 + X_iz) = Σe_nz^n;
- Q(y) = Σq_ny^n = E(−ty)/E(−y);
- [n] = [n]_t, and (a;t)_n = ∏_{i<n}(1 − at^i).

Define, for 0 ≤ j ≤ n (and c(n,j) := 0 otherwise),

  α_j := ∏_{i=1}^{j}(s − t^i)/(t;t)_j  (α_{−1} := 0),
  **c(n,j) := (s;t)_{n−j}/(t;t)_{n−j} · (α_j − s t^{n−j} α_{j−1}),   F_n(w) := Σ_{j=0}^{n} c(n,j) w^j.**

This is PROVE.md's F_n. The identity M^{(k)}_{jb} = (t;t)_k s^b t^{−kj} c(k−b, j) is checked symbolically for k ≤ 6, and by hand it is just R_{jb} = (s − t^j) − s t^{n−j}(1 − t^j).

**Theorem.** For all m ≥ 0, k ≥ 1 and r ≥ 0, in Λ_m:

  **t^{−C(k,2)} e_k(Y)•e_r = Σ_{b=0}^{k} s^b F_{k−b}(t^{r−b}) e_b e_{r+k−b}.**

Equivalently, as a generating function (E_k):

  Σ_r z^r t^{−C(k,2)} e_k(Y)•e_r = Σ_{j+b≤k} s^b t^{−kj} c(k−b,j) · e_b z^{b−k} E(t^j z).

Via Hikita Def 3.4 / Lemma 3.3, the left side is e_k⋆e_r.

For r < k, distinct b can give the same monomial e_be_{r+k−b}; the identity is between sums. For example, r = 0 gives e_k(Y)•1 = t^{C(k,2)}e_k, which is Lemma 3.3, re-derived here.

## 1. Hecke facts used (textbook)

- **(H1)** The T_i satisfy the braid relations and (T_i − t)(T_i + 1) = 0. So T_w := T_{i_1}⋯T_{i_ℓ} is well defined for any reduced word of w ∈ S_m (Matsumoto), and T_uT_v = T_{uv} whenever ℓ(uv) = ℓ(u) + ℓ(v).
  - The relations are classical for Demazure–Lusztig operators and are machine-checked in every pipeline since Day 198.
  - The words in play are s_{i_1}⋯s_{i_ℓ}, read as function composition, rightmost factor first.
- **(H2)** For J ⊆ S, let W^J = {w : ℓ(ws) > ℓ(w) ∀s ∈ J}. Each w factors uniquely as w = w^J w_J with w^J ∈ W^J, w_J ∈ W_J and ℓ(w) = ℓ(w^J) + ℓ(w_J). In particular w^J is the unique minimal element of wW_J.
- **(H3) Chain.** For K ⊆ J, the map (u, v) ↦ uv is a length-additive bijection W^J × (W_J)^K → W^K.
  - *Proof.* Let u ∈ W^J, v ∈ (W_J)^K and s ∈ K. Then ℓ(uvs) = ℓ(u) + ℓ(vs) = ℓ(u) + ℓ(v) + 1, since vs ∈ W_J. So uv ∈ W^K.
  - Conversely, let x ∈ W^K and write x = x^J x_J. If x_J s < x_J for some s ∈ K, then ℓ(xs) = ℓ(x^J) + ℓ(x_Js) < ℓ(x), which is a contradiction. So x_J ∈ (W_J)^K.
  - Uniqueness follows from (H2). ∎
- **(H4)** If T_kH = tH for all s_k ∈ J, then T_vH = t^{ℓ(v)}H for v ∈ W_J, and T_k^{−1}H = t^{−1}H.
  - A function symmetric under s_k satisfies T_kF = tF.
- **(H5)** πT_k = T_{k+1}π for 1 ≤ k ≤ m−2 (Day 206b (iv)).

**Minimal coset reps for J_n := ⟨s_i : i ≠ n⟩ = S_n × S_{m−n}.**
- wW_{J_n} = {u : u([n]) = w([n])}. The minimal element w^D, with D = w([n]), is increasing on [n] and on [n+1, m], and ℓ(w^D) = Σ_l (d_l − l).
- For D = {d_1 < ⋯ < d_n}, the word

  **w(D) := (s_{d_1−1}⋯s_1)(s_{d_2−1}⋯s_2)⋯(s_{d_n−1}⋯s_n)**

  sends l ↦ d_l for every l ≤ n. To see this:
  - the blocks with index l' > l fix l;
  - block l sends l to d_l;
  - the blocks l' < l only move points ≤ d_{l'} < d_l.
- So w(D) lies in the coset and has length ≤ Σ(d_l − l) = ℓ(w^D). Hence w(D) = w^D, and the word is reduced.
- Write T_{w(D)} for the corresponding operator; σ^{(k)} := Σ_{|D|=k} T_{w(D)} is PROVE.md's coset sum.

## 2. (A_k): the parabolic kernel

**Key Lemma (KL).** Let 0 ≤ n ≤ m and 1 ≤ c < d_1 < ⋯ < d_n ≤ m. Let H satisfy T_kH = tH for all k ≠ n (1 ≤ k ≤ m−1). Then

  (T_cT_{c+1}⋯T_{m−1})^{−1} T_{w(D)} H = t^{n−(m−c)} T_{w(D−1)} H,  D − 1 := {d_1−1, …, d_n−1}.

*Proof.* Let x := s_cs_{c+1}⋯s_{m−1}. It is reduced of length m−c, and it acts as m ↦ c and i ↦ i+1 for c ≤ i ≤ m−1, fixing i < c.

Let y := w(D−1). Its word has length Σ(d_l − 1 − l) = ℓ(w^D) − n, and it is valid because d_l − 1 ≥ c + l − 1 ≥ l.

- **y(m) = m.** The value m is not in D−1, so it is the largest value on the increasing block [n+1, m], and it sits at position m.
- **ℓ(xy) = ℓ(y) + (m−c).** x is increasing on [1, m−1], so the inversions of y among positions < m survive. At position m, xy(m) = c, and it is inverted with exactly those i < m with y(i) ≥ c. There are m − c of these, namely all the values c..m−1.
- **xy ∈ w^D W_{J_n}.** We have xy([n]) = x(D−1) = D, because every d_l − 1 ≥ c. So xy = w^D v with v ∈ W_{J_n}, and ℓ(v) = ℓ(xy) − ℓ(w^D) = m − c − n.

By (H1) and (H4), T_xT_yH = T_{xy}H = T_{w^D}T_vH = t^{m−c−n}T_{w^D}H. Now apply T_x^{−1}. ∎

**Proposition (A_k).** For symmetric F and 1 ≤ b_1 < ⋯ < b_k ≤ m:

  Y_{b_1}Y_{b_2}⋯Y_{b_k} F = t^{C(k,2)} T_{w(b)} π^k F.

Hence **e_k(Y)F = t^{C(k,2)} σ^{(k)} π^k F**.

*Proof.* Induct on the number of factors; the case of zero factors is trivial.

Let c < d_1 < ⋯ < d_n and assume Y_{d_1}⋯Y_{d_n}F = t^{C(n,2)}T_{w(D)}H, where H := π^nF.
- H = X_1⋯X_n F(X_{n+1}..X_m, sX_1..sX_n) is symmetric in the head and in the tail, so T_kH = tH for k ≠ n.

By KL,

  Y_c(t^{C(n,2)}T_{w(D)}H) = t^{C(n,2)} t^{m−c} T_{c−1}⋯T_1 π (T_c⋯T_{m−1})^{−1}T_{w(D)}H = t^{C(n,2)+n} T_{c−1}⋯T_1 π T_{w(D−1)} H.

Then use (H5). The indices in w(D−1) lie in [1, m−2], so

  π ∏_l (T_{d_l−2}⋯T_l) = ∏_l (T_{d_l−1}⋯T_{l+1}) π.

Finally, T_{c−1}⋯T_1 · ∏_l(T_{d_l−1}⋯T_{l+1}) is exactly the word w({c} ∪ D), and πH = π^{n+1}F. Since C(n,2) + n = C(n+1,2), the induction closes. ∎

(At k = 2 this is Day 206b A2. There, the braid-shift bookkeeping was done by hand. Here it is a single length count.)

## 3. (K_k): the k-set kernel

Let σ_{[l,m]} := Σ_{a=l}^{m} T_{a−1}⋯T_l, the partial symmetrizer on X_l..X_m. Let K_l := ⟨s_{l+1}, …, s_{m−1}⟩.

**(C1) Chain factorization.**
- The minimal reps of W_{K_{l−1}}/W_{K_l} are s_{a−1}⋯s_l for l ≤ a ≤ m. The coset is determined by the image of l, and the stated rep has the minimal a − l inversions.
- By (H3) along K_k ⊂ K_{k−1} ⊂ ⋯ ⊂ K_0:

  Σ_{u∈W^{K_k}} T_u = σ_{[1,m]} σ_{[2,m]} ⋯ σ_{[k,m]}.

- By (H3) along K_k ⊂ J_k, with (W_{J_k})^{K_k} = S_k = ⟨s_1..s_{k−1}⟩:

  Σ_{u∈W^{K_k}} T_u = σ^{(k)} · Σ_{v∈S_k} T_v.

- So if G is symmetric in the head X_1..X_k and in the tail, then (H4) gives

  **σ_{[1,m]}⋯σ_{[k,m]} G = [k]! σ^{(k)} G.**

(At k = 2 this is Day 206b's coset identity (a).)

**(C2) Ordered formula (iterated Lemma 1).** Let G be symmetric in the tail X_{k+1}..X_m only. Then

  σ_{[1,m]}⋯σ_{[k,m]} G = Σ_{(a_1,…,a_k) distinct} G(X_{a_1},…,X_{a_k}; rest) · ∏_{l=1}^{k} ∏_{j∉{a_1..a_l}} a_{a_l j}.

*Proof.* Induct on k; k = 1 is Day 205b Lemma 1.

- **Inner step.** Apply the induction hypothesis to H := σ_{[2,m]}⋯σ_{[k,m]}G in the variables X_2..X_m, over ℚ(s,t)(X_1). The head is at positions 2..k. H is the displayed sum over distinct tuples in [2,m], and it is symmetric in X_2..X_m, since a permutation of those variables permutes the terms.
- **Outer step.** Lemma 1 gives σ_{[1,m]}H = Σ_{a_1} H^{(1↔a_1)} ∏_{j≠a_1} a_{a_1 j}.
- **The swap.** Under X_1 ↔ X_{a_1}, the tuples in [2,m] become tuples in [m]∖{a_1}. Each index set [2,m]∖{…} becomes [m]∖{a_1, …}. ∎

**(C3) Symmetrization.** Σ_{w∈S_k} ∏_{l<l'} a_{w_l w_{l'}} = [k]!. This is Macdonald III (1.4); here it is re-derived.
- Take m = k and G = 1 in (C1) and (C2). The left side is Σ_{u∈S_k} T_u 1 = Σ t^{ℓ(u)} = [k]!.
- The right side is the displayed sum.

**Proposition (K_k).** For G symmetric in head and tail,

  σ^{(k)}G = Σ_{|A|=k} G^{(A)} ∏^×_A.

*Proof.* In (C2), G(X_{a_1},…) = G^{(A)} depends only on A = {a_1..a_k}. The weight splits as

  ∏_{i∈A, j∉A} a_{ij} · ∏_{l<l'} a_{a_l a_{l'}}.

Summing over the orderings of A gives [k]! ∏^×_A by (C3). Compare with (C1) and divide by [k]!. ∎

## 4. (E_k): the generating function, by induction on k

Put φ(u) := u(1+suz)/(1+uz). Since (π^ke_r)^{(A)} = X_A e_r(X_{A^c}, sX_A) = [z^r] E(z)∏_{a∈A}φ(X_a), §§2–3 give

  Γ_k(X) := Σ_r z^r t^{−C(k,2)} e_k(Y)•e_r = E(z) Σ_{|A|=k} ∏^×_A ∏_{a∈A} φ(X_a).

This holds for every m ≥ 0. For m < k both sides are 0, since the product of k distinct Y's is empty.

Let T_k(X) := Σ_{j+b≤k} N^{(k)}_{jb} e_b z^{b−k}E(t^jz), with N^{(k)}_{jb} := s^b t^{−kj} c(k−b, j). Both are Laurent polynomials in z over ℚ(s,t)[X]^{S_m}. We show Γ_k = T_k.

**(R) Recursion.** For k ≥ 1 and m ≥ 1, with X̂_i the other m−1 variables:

  [k] Γ_k(X) = Σ_i X_i(1 + sX_iz) ∏_{j≠i} a_{ij} · Γ_{k−1}(X̂_i).

*Proof.*
1. By (C2) and (C3), [k]! Σ_A ∏^×_A Φ(X_A) = Σ_{ordered} ∏_l ∏_{j∉{a_1..a_l}} a_{a_l j} Φ. This holds for any symmetric Φ, coefficientwise in z.
2. Separate a_1 =: i. The remaining ordered sum is [k−1]! times the (k−1)-set sum in X̂_i.
3. Use φ(X_i)E(X;z) = X_i(1+sX_iz)E(X̂_i;z). ∎

**(L) Step lemma.** For every m ≥ 0 and k ≥ 1,

  Σ_i X_i(1+sX_iz) ∏_{j≠i}a_{ij} · T_{k−1}(X̂_i) = [k] T_k(X).

For m = 0 the left side is the empty sum.

*Proof.* Fix a term N'_{jb} e_b(X̂_i) z^{b−k+1} E(X̂_i; γ) of T_{k−1}(X̂_i), where γ := t^jz and N' := N^{(k−1)}. Use e_b(X̂_i) = Σ_{c=0}^b (−X_i)^c e_{b−c} and E(X̂_i;γ) = E(γ)/(1+γX_i). The term times X_i(1+sX_iz) becomes

  N' z^{b−k+1} E(γ) Σ_c (−1)^c e_{b−c} g_c(X_i),  g_c(u) = u^{1+c}(1+suz)/(1+γu) = s t^{−j}u^{1+c} + (1 − s t^{−j}) u^{1+c}/(1+γu).

The coefficients are now fully symmetric and g_c(0) = 0. So Day 205b Lemma 2 applies, coefficientwise in z: Σ_i g(X_i)∏_{j≠i}a_{ij} = (1−t)^{−1}Ω[g], with Ω[u^n] := q_n for n ≥ 1. For m = 0 both sides are 0.

- **(a) The s-part.** Σ_c (−1)^c e_{b−c} q_{1+c} = (1 − t^{b+1}) e_{b+1}. This is [y^{b+1}] of Q(y)E(−y) = E(−ty).
- **(b) The pole part.** We have Σ_c(−1)^c e_{b−c}u^c = [w^b]E(w)/(1+uw), and the exact power-series identity

  Ω[u/((1+wu)(1+γu))] = Σ_{a,p≥0}(−w)^a(−γ)^p q_{1+a+p} = (Q(−γ) − Q(−w))/(w − γ).

  So E(γ) × (part b) = [w^b] P(w,γ), where

  P(w,γ) := (E(w)E(tγ) − E(γ)E(tw))/(w − γ).

  Here E(x)Q(−x) = E(tx) has been used twice. P is a polynomial, since the numerator vanishes at w = γ. Writing the numerator as Σ_n e_n w^n (E(tγ) − t^nE(γ)) gives

  [w^b]P = Σ_{n>b} e_n γ^{n−1−b}(E(tγ) − t^nE(γ)) = −Σ_{n≤b} e_n γ^{n−1−b}(E(tγ) − t^nE(γ)).

  The two forms agree because the full sum over n is γ^{−1−b}(E(γ)E(tγ) − E(tγ)E(γ)) = 0.

Collect the coefficient of e_{b'} z^{b'−k} E(t^jz) and multiply by (1−t). It is

  s t^{−j}(1−t^{b'}) N'_{j,b'−1} + (1−st^{−j}) Σ_{b≥b'} t^{b'} t^{j(b'−1−b)} N'_{jb} − (1 − st^{1−j}) Σ_{b≥b'} t^{(j−1)(b'−1−b)} N'_{j−1,b}.

The three terms come from, respectively:
- (a) at (j, b'−1);
- the E(γ)-part of (b) at (j, b);
- the E(tγ)-part of (b) at (j−1, b), because t·t^{j−1}z = t^jz.

Now insert N'_{jb} = s^b t^{−(k−1)j}c(k−1−b, j), and put n := k − b' and p := b − b'. The coefficient becomes s^{b'}t^{−kj} times

  (1 − t^{b'}) c(n,j) + t^{b'} [ (1−st^{−j}) Σ_{p≥0} (st^{−j})^p c(n−1−p, j) − t^n (1−st^{1−j}) Σ_{p≥0} (st^{1−j})^p c(n−1−p, j−1) ].

By (★) below, the bracket equals (1 − t^n)c(n,j). So the total is (1 − t^{b'} + t^{b'} − t^{b'+n}) c(n,j) = (1 − t^k)c(n,j). Dividing by (1−t) gives [k]N^{(k)}_{jb'}. ∎

**(★) The one-index identity.** For n ≥ 1 and all j ≥ 0:

  (1 − t^n) c(n,j) = (1 − st^{−j}) Σ_{p≥0} (st^{−j})^p c(n−1−p, j) − t^n (1 − st^{1−j}) Σ_{p≥0} (st^{1−j})^p c(n−1−p, j−1).

*Proof.* Let G(x) := Σ_{n≥0} (s;t)_n/(t;t)_n x^n. Comparing coefficients gives (1−x)G(x) = (1−sx)G(tx). Let C_j(x) := Σ_n c(n,j)x^n. Then

  C_j(x) = x^j [α_j G(x) − sα_{j−1} G(tx)] = x^j G(tx) [α_j(1−sx)/(1−x) − sα_{j−1}].

In generating functions, (★) for all n ≥ 0 (n = 0 is 0 = 0) reads

  C_j(x) − C_j(tx) = (1−st^{−j}) x C_j(x)/(1 − st^{−j}x) − (1−st^{1−j}) tx C_{j−1}(tx)/(1 − st^{2−j}x).

Move the first right-hand term to the left. Using C_j(x)[1 − (1−st^{−j})x/(1−st^{−j}x)] = C_j(x)(1−x)/(1−st^{−j}x), this becomes

  C_j(x)(1−x)/(1−st^{−j}x) − C_j(tx) = −(1−st^{1−j}) tx C_{j−1}(tx)/(1−st^{2−j}x).   (★★)

Put D_j := α_j − sα_{j−1}. Since α_j = α_{j−1}(s − t^j)/(1 − t^j) for j ≥ 1, we get the key linear factorization

  α_j(1 − sy) − sα_{j−1}(1 − y) = D_j (1 − st^{−j}y)   (for all y; both sides are linear in y and agree at y = 0 and y = t^j/s).

It holds trivially for j = 0.

Apply it with y = x and y = tx, and G(tx) = G(t²x)(1−stx)/(1−tx):
- the left side of (★★) is x^jG(t²x) D_j [(1−stx) − t^j(1 − st^{1−j}x)]/(1−tx) = x^jG(t²x) D_j(1 − t^j)/(1−tx);
- the right side is −x^jG(t²x) t^j(1−st^{1−j}) D_{j−1}/(1−tx). Here the factor (1 − st^{2−j}x) cancels by the same factorization at j−1.

So (★★) is equivalent to (1 − t^j) D_j = t(s − t^{j−1}) D_{j−1}. For j ≥ 1, D_j = α_{j−1} t^j (s−1)/(1−t^j), and both sides equal t^j(s−1)α_{j−1}:
- for j = 1, D_0 = 1;
- for j ≥ 2, use α_{j−1} = α_{j−2}(s−t^{j−1})/(1−t^{j−1}).

For j = 0 both sides of (★★) vanish, because D_0(1 − t^0) = 0 and C_{−1} = 0. ∎

**Proof of (E_k).** Induct on k, for all m at once. For k = 0, Γ_0 = E(z) = T_0, since c(0,0) = 1. Let k ≥ 1.
- **m = 0.** Γ_k = 0, and (L) gives [k]T_k = 0. As a by-product, Σ_j t^{−kj}c(k,j) = 0.
- **m ≥ 1.** By (R), the induction hypothesis in m−1 variables, and (L):

  [k]Γ_k = Σ_i(…)Γ_{k−1}(X̂_i) = Σ_i(…)T_{k−1}(X̂_i) = [k]T_k. ∎

**Corollary (the Pieri rule).** Since [z^r] z^{b−k}E(t^jz) = t^{j(r+k−b)}e_{r+k−b},

  t^{−C(k,2)}e_k(Y)•e_r = Σ_b s^b e_b e_{r+k−b} Σ_j c(k−b,j) t^{j(r−b)} = Σ_b s^b F_{k−b}(t^{r−b}) e_b e_{r+k−b}. ∎

## 5. What was actually new today, and why it was short

- **(A_k).** Day 206b's hand braid-shift is replaced by one Coxeter length count (KL): T_x^{−1}T_{w^D} = t^{−ℓ(v)}T_{w^{D−1}} on J-spherical vectors.
- **(K_k).** Parabolic chain factorization plus iterated Lemma 1. The Poincaré identity (C3) falls out of the same two facts at m = k.
- **(E_k).** PROVE.md's advice held: "do NOT fall back on symbolic partial fractions; reuse the proved atoms."
  - The induction peels off the **outer** variable, with Lemma 1 at the top level, rather than stacking inner Ω-stages.
  - The induction hypothesis is then only ever evaluated at a single variable X_i. That is a one-pole Lemma 2 computation, with poles at u = −1/γ and u = −1/w only.
  - The e_b(X̂_i) expansion is absorbed by the (w,γ) difference quotient.
  - The M-recursion loses both k and b and becomes the single-index identity (★). The two-term shape of R_{jb} is exactly the pair of (a)/(b) contributions.
- **The chain telescoping** is E(γ)Q(−γ) = E(tγ), applied once per step. The "only chain prefixes occur" phenomenon from Day 207 is simply that each step raises j by at most 1.

## 6. Verification

Everything is in `scripts/day207b/`; see the header. KL, A_k, (C2), (★), (L) and the final formula were each checked independently. In particular:
- (L) is checked with z symbolic, which catches algebra slips in the proof itself and not just in the final statement;
- the final formula is checked against direct AHA computation for m < k, which confirms the m = 0 base-case mechanism.

## 7. Scope, gaps, credits

- **Scope.** A polynomial identity in Λ_m ⊗ ℚ(s,t), for every m ≥ 0, k ≥ 1, r ≥ 0. No division by anything depending on m.
- **Gaps.**
  - None in the operator identity.
  - The ⋆-reading uses Hikita Def 3.4 / Lemma 3.3 (verified-quote), as for W_r.
  - The Hecke facts (H1)–(H2) are textbook. The braid relations for these particular T_i are classical, not re-proved here, and are machine-checked.
- **Not proved.** The ₂φ₁ form F_n(w) = (1 − t^nw)(s;t)_n/(t;t)_n · ₂φ₁(t^{1−n}, t/s; t^{1−n}/s; t; w) remains `computed`, n ≤ 6. It is not needed; C_j gives an equally explicit form. It is a plausible one-afternoon check, via contiguity, later.
- **Credits.**
  - The kernel form of e_k(Y) on symmetric functions is classical in the q-shift DAHA (Macdonald D-operators; Concha–Lapointe 2307.02385 Lemmas 8 and 10 give the template).
  - The level-1 π^k transfer, the per-tuple KL argument in this convention, and the e-basis Pieri rule for all k are new per `memory/reading/2026-09-26-novelty-ek-star-er.md`. That audit is not re-run today, since there was no browsing.
- **Consequences.**
  - Day 193 (e_3⋆e_r) and Day 195 (e_4⋆e_r) closed forms: PROVED for all r (`check_specializations.py`, symbolic in u).
  - The support claim of the min(k,r)+1 meta-conjecture: PROVED. The support is {e_b e_{r+k−b} : 0 ≤ b ≤ k}.
