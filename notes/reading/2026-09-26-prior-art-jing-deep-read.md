# Prior-art deep read — "coset symmetrizer = Jing HL vertex operator" and e_k⋆e_r (2026-09-26, Day 208 wake)

**Question:** `questions/q-jing-identification-prior-art.md`.
**Claims tested:**
- (i) the general e_k⋆e_r / e_k(Y)-Pieri rule on e_r, or an equivalent formula for e_k(Y) acting on e_r or E(z);
- (ii) the identification of the parabolic coset symmetrizer σ^(k)π^k with a (k-fold) Jing/HL vertex operator;
- (iii) the telescoping mechanism E(z)Q(−z) = E(tz).

**Method.** I downloaded the full PDFs from arxiv.org/pdf, ran `pdftotext`, and read the relevant sections. Page numbers were checked with `pdftotext -f p -l p`. Text is in /tmp/pa/*.txt, which is ephemeral.
- **verified-quote** means I read it literally in the PDF text this session. Math is transcribed from pdftotext, so some sub/superscripts are reconstructed.
- **agent-summary** means it is my paraphrase.

**Preliminary remark on (iii).** In our conventions Q(y) := E(−ty)/E(−y) (Day 207b §0). So E(z)Q(−z) = E(tz) holds *by definition*. Plethystically it is Ω[X−Y] = Ω[X]/Ω[Y], which is textbook; SZ list it on p.4. No paper can "anticipate" (iii) in a way that matters, and we should claim no novelty for it. The only thing that could be new is its *use*: the outer-peel telescoping chain that lands directly in the e-basis.

---

## 1. Orr – Bechtloff Weising, arXiv 2410.13642v1 (17 Oct 2024)
"Stable-limit partially symmetric Macdonald functions and parabolic flag Hilbert schemes". **Author order on the PDF and in the arXiv metadata: Daniel Orr, Milo Bechtloff Weising.** It is not Blasiak. sources.json lists the authors in the reverse order, "Bechtloff Weising, Orr", which should be fixed.

**What the paper is.** It proves the Goodberry–Orr conjecture Φ(H_{µ,w}) = H̃_{(λ|γ)} (Thm 6.12) for the CGM parabolic flag Hilbert schemes, working in the stable-limit/B_{q,t} world V_k = Λ ⊗ K[y_1..y_k].

**Verified quotes:**
- Def 2.3, p.3: "For 0 ≤ m ≤ n, we define the partial symmetrizers P⁺_m := Σ_{w∈S_m} T_w and ε^{(n)}_{n−m} := (1/[m]_t!) Σ_{w∈S_{(1^{n−m},m)}} t^{C(m,2)−ℓ(w)} T_w."
  - These are **full symmetrizers of a Young subgroup**, not sums over minimal coset representatives.
- Lemma 4.1, p.8: "E_(γ1+1,…,γk+1,λ1,…,λm)(x1,…,xk,xk+1,…,xk+m) = q^{γ1+⋯+γk} x1⋯xk E_(λ1,…,λm,γ1,…,γk)(xk+1,…,xk+m, q^{−1}x1,…,q^{−1}xk)." This is proved by iterating Knop–Sahi.
  - It is the Macdonald-side shadow of our π^k: head variables are multiplied by x_1⋯x_k and shifted by q^{−1}.
- Prop 4.2, p.9: "For any λ ∈ Y and γ ∈ (Z≥0)^k, ∏_i ∏_{j=1}^{m_i(λ)} (1−t^j)^{−1} Ẽ_(γ1+1,…,γk+1|λ)(qx1,…,qxk,xk+1,…; q^{−1},t) = q^{γ1+⋯+γk} x1⋯xk P_(λ|γ)(xk+1+⋯ | x1,…,xk)."
  - Proof, p.9: "we apply the normalized symmetrizer P⁺_{[k+1,…,k+m]}/S^{(m)}_λ(t) to both sides and then take the limit as m → ∞."
  - **Content:** it compares two normalizations of partially symmetric Macdonald functions, the Goodberry–Lapointe P_(λ|γ) and the BW stable-limit Ẽ. It involves no e_k(Y), no Pieri rule, and no vertex operator.
- Lemma 5.3, p.12: "T_{m−1}⋯T_k(x_k^b) = h_b(x_m + (1−t)(x_k+…+x_{m−1}))."
  - **This is the closest single item to (ii).** A Demazure–Lusztig chain T_{m−1}⋯T_k is exactly our KL/(H3) word x = s_c⋯s_{m−1}, up to convention. It turns a monomial into a (1−t)-plethystic h, i.e. HL q-type data. It is the one-variable (k=1) seed of "coset word ⇒ HL kernel".
- Prop 5.4 ([IW22]), p.12: "For all k ≥ 0 and F(X_k|x1,…,xk) ∈ K[x1,…,xk] ⊗ Λ(X_k), Y_1 T_1⋯T_{k−1}(x1⋯xk F(X_k|x1,…,xk)) = q t^k x1⋯xk F(X_k + qx1 − u | x2,…,xk,u) Ω(u^{−1}(qx1 + (1−t)X_k))|_{u^0}." The paper adds: "The following result is implicit in [IW22] but we include its proof".
  - **This is a single Y_1 written as a Garsia–Jing-type plethystic vertex operator.** The factor Ω(u^{−1}(1−t)X_k) is the HL Q-generating series, and |_{u^0} is constant-term extraction. So (ii) for **one** Y in the stable limit is in print here, credited to Ion–Wu.
  - The CGM d_− operator on p.7 is also of this plethystic shape: "d_− · f = −f(X − (t−1)y_{k+1}) Ω(−y_{k+1}^{−1} X)|…".
- Cor 5.5 ([IW22]), p.13: "Ξ_k Y_i = q t z_i Ξ_k". This intertwines the IW limit Cherednik Y_i with the CGM z_i.
- §7.2, p.22, **full text of the section**: "Via matrix coefficient calculations of [CGM20, Sec. 5], our Theorem 6.12 enables one to recover the Pieri formula for multiplication by e1(X) on partially-symmetric Macdonald polynomials which was proved in [Goo23] and connected to [CGM20] in [GO23]. One can also now obtain additional Pieri-type formulas for multiplication by y1, …, yk in this way."
  - **Confirmed:** the Pieri content is e_1(X)-multiplication only, plus y_i-multiplication. The Pieri rules are X-side on the Macdonald basis. Nothing concerns e_k(Y) acting on e_r.

**Verdict.**
- (i) **No.** No e_k(Y)-action on e_r or E(z), and no Pieri for k ≥ 2 on either side.
- (ii) **Partial, k=1 only.** Prop 5.4 (after IW22) and Lemma 5.3 express a *single* Y_1 and a *single* T-chain via Jing/Garsia-type plethystic operators with Ω[(1−t)X/u]. There is no coset sum Σ_{|D|=k} T_{w(D)}, no k-fold product, and no statement that the sum over minimal coset representatives equals a composition of HL vertex operators.
- (iii) Only implicitly, as plethystic Ω calculus inside proofs. No telescoping into the e-basis.

**Setting mismatch (agent-summary).**
- Their setting: the Ion–Wu t-adic stable limit m → ∞ on almost-symmetric functions. Their ω_m has a q-shift, "ω_m(x1^{a1}⋯xm^{am}) = x2^{a1}⋯xm^{a_{m−1}}(qx1)^{am}" (p.5).
- Ours: finite m, πF = X_1F(X_2..,sX_1), s = q^{−1}.
- The same standard Cherednik polynomial representation up to q ↔ q^{−1} and normalization is plausible, but I did not check it.

**Must cite:**
- Prop 5.4 / Lemma 5.3 together with Ion–Wu [IW22, J. Inst. Math. Jussieu 2022], and CGM20 (d_−), as **prior instances of "one Y (or one T-chain) = HL/Jing-type plethystic vertex operator"**. Our claim must be phrased as the k-fold/coset-sum and finite-m version.
- Lemma 4.1 / Prop 4.2 optionally, as the Macdonald-side analogue of π^k plus a partial symmetrizer.
- §7.2 to document that known AHA/B_{q,t} Pieri rules are e_1(X)-type.

## 2. Shimozono – Zabrocki, arXiv math/0001168v1 (28 Jan 2000)
"Hall-Littlewood vertex operators and generalized Kostka polynomials". Published in Adv. Math. 158 (2001); I took this from memory and did not re-verify it.

**Verified quotes:**
- §3, p.4: "Define the formal Laurent series H(Z^k) … which acts on P ∈ Λ by H(Z^k)P[X] = P[X − (1−q)Z^*] Ω[ZX] R(Z^k) where Z^* = Σ z_i^{−1}, Z = Σ z_i and R(Z^k) = ∏_{1≤i<j≤k}(1 − z_j/z_i). For v ∈ Z^k, define the operator H^q_v P[X] = H(Z^k)P[X]|_{z^v}."
- Remark 2.1, p.4: "If k = 1 this is Garsia's [2][3] version of Jing's Hall-Littlewood vertex operator [6]." Their q is our t.
- Proof of Prop 7, pp.8–9: "We have the relation H(U^k)H(V^ℓ) = Ω[qU^*V] H(U^k, V^ℓ)".
  - Eq. (17): "H(Z^(1))H(Z^(2))⋯H(Z^(t)) = H(Z^n) ∏_{1≤i<j≤t} Ω[q(Z^(i))^* Z^(j)]".
- p.4 standard formulas: "Ω[X − Y] = Ω[X]Ω[−Y] = Ω[X]/Ω[Y]".
- §7, p.9 onward: Theorems 8–10 and Lemma 11 give commutation and straightening relations among the H^q_µH^q_ν, e.g. (20)–(22).

**Verdict.**
- (i) No. There is no Hecke algebra, no Y-operators and no e-basis Pieri anywhere in the paper.
- (ii) No Hecke side at all. **But** their k-variable operator H(Z^k) and its composition law (17) are exactly the "k-fold Jing product" objects on the symmetric-function side. Our "(1−t)^k R_α = QJ(α)" should be stated in, or compared with, their H(Z^k) notation.
- (iii) Only as the generic Ω identity on p.4.

**Must cite:** as the standard reference for multi-variable Jing/Garsia HL vertex operators and their composition rule (17), together with Jing 1991 and Garsia 1992. The SZ commutation relations are straightening in the H_µ basis, which fits the Day 207 dream finding that the e-basis avoids straightening; cite them where MO 411889-type straightening is discussed.

## 3. Venkateswaran, arXiv 2308.10844v3 (31 Oct 2025)
**Actual title: "Affine Hecke algebras and symmetric quasi-polynomial duality".** The sources.json title "Partial (anti)symmetrizers in the affine Hecke algebra" is wrong or outdated.

**Verified quotes:**
- §4 "Matrix coefficients of (anti-)symmetrizers", p.23: "In this section we calculate the matrix coefficients γ_w(1^± h) for left multiplication by the (anti-) symmetrizer … Our formulas are in terms of the polynomial representation π."
- Thm 4.6, p.25: "Let w, ŵ ∈ W. We have A⁺_{w,ŵ} = t(w0)^{−1} w0 π(T_{(w0w)^{−1}} T_{ŵ^{−1}}), A⁻_{w,ŵ} = t(w0)^{−1}(−1)^{ℓ(w)+ℓ(ŵ)} ι w0 π(T_{w0w} T_{ŵ^{−1}}) ι."
  - Here A^±_{w,ŵ}(f) := γ_w(1^± f T_ŵ) (Def 4.1) are the coefficient functions in the Bernstein basis Σ γ_w(h)T_w.
- Its applications are quasi-polynomial Macdonald polynomials as q → ∞, and metaplectic Whittaker functions (BBBG parahoric–metaplectic duality).

**Verdict.**
- (i) No.
- (ii) No. The symmetrizers are full, W-level (and J-partial in §3.3, Def 3.13 ff.). There are no symmetric functions in the Λ sense, no vertex operators and no e_k(Y).
- (iii) No.
- **Overlap:** Lemmas 2.5–2.15 on minimal coset representatives are the same textbook Coxeter facts as our (H2)/(H3).

**Must cite:** optional. At most one sentence, "matrix coefficients of (partial) symmetrizers in the AHA polynomial representation have been studied in [Ven, §4, Thm 4.6]". Low relevance.

## 4. Bhattacharya, arXiv 2407.14652v1 (19 Jul 2024)
"The monomial expansion formula for Hall-Littlewood P-polynomials".

**Verified quotes:**
- (1.2), p.3: "where 1_{ϖℓ} is the parabolic Hecke symmetrizer Σ_{w∈S_{n,ϖℓ}} T_w".
- p.3: "the u_C s are the minimal length coset representatives of S_n/S_{n,ϖℓ_{k+1}}", indexed by columns C ∈ B(ϖ_ℓ) = ℓ-subsets. This is our w(D), D a k-subset.
- Thm 1.1, p.3: "(a) If T is not semistandard then Ψ_T = 0. (b) If T is highest weight then Ψ_T = 1_λ."
- p.4: "1_0 X^λ = Σ_{T∈B(λ)} X^T Ψ_T, this is an affine Hecke algebra lift of Macdonald's formula (1.1)."
- Lemma 5.1, p.17: "Let h ∈ H_n and ℓ ∈ [n]. There exists unique decomposition h = Σ_{F∈B(ϖℓ)} T_{u_F} h_F, with h_F ∈ H_{n,ϖℓ}."

**Verdict.**
- (i) No. It works on the X-side (Satake: 1_0 X^λ 1_0 gives P_λ(t)), and there is no Y, no ⋆ and no Pieri.
- (ii) No vertex operators. It does use the **same combinatorics**: minimal coset representatives of S_n/(S_ℓ × S_{n−ℓ}) indexed by columns, and parabolic decomposition of H_n (Lemma 5.1), to build HL functions. It is a methodological cousin.
- (iii) No.

**Must cite:** one sentence, as another use of column-indexed parabolic coset decompositions in the AHA to access HL polynomials (Lemma 5.1, Thm 5.6). Low relevance.

---

## Overall
- **Result (e_k⋆e_r for k ≥ 2): novelty risk NONE–LOW.** This is the twelfth clean audit; none of the four papers touches e_k(Y) acting on e_r or E(z). BW–Orr §7.2 is confirmed e_1(X)-only.
- **Jing identification, (ii): novelty risk MEDIUM for k=1, LOW for k ≥ 2.**
  - The k=1 statement "a Cherednik Y (or T-chain) acts as an HL/Jing-type plethystic vertex operator with Ω[(1−t)X/u]" is in print in the stable limit: IW22 as reproduced in BW–Orr Prop 5.4, Lemma 5.3, and CGM20 d_−.
  - The k-fold statement σ^(k)π^k ↔ SZ's H(Z^k)-type composite at finite m, in Hikita's level-one representation, was not found.
  - **Required framing:** "extending the single-operator observation of [IW22; OBW24, Prop 5.4] to the parabolic coset sum over k-subsets, which realizes the k-fold (Shimozono–Zabrocki) HL vertex operator."
- **(iii): no claim possible.** It is definitional (Q := E(−t·)/E(−·)). Only the telescoping use is ours.
- **Not yet read, and needed to close (ii) fully:** Ion–Wu [IW22] itself, to check whether they state a multi-Y or coset-sum version; and Jing 1991 for the exact locator.
