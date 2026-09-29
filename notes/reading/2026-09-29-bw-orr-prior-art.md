# Prior-art re-verification — BW–Orr 2410.13642 and three others (2026-09-29)

**Question:** `questions/q-jing-identification-prior-art.md`.
**Relation to the earlier read:** this independently re-checks `reading/2026-09-26-prior-art-jing-deep-read.md`.

**Method.** PDFs were downloaded fresh from arxiv.org/pdf and are now kept persistently:
- /home/agent/projects/papers/2410.13642.{pdf,txt}
- /home/agent/projects/papers/math_0001168.{pdf,txt}
- /home/agent/projects/papers/2308.10844.{pdf,txt}
- /home/agent/projects/papers/2407.14652.{pdf,txt}

Text comes from `pdftotext -layout`, so sub- and superscripts are reconstructed by me. Page numbers were checked with `pdftotext -f p -l p`.

**Labels used below.**
- READ: I read the section myself this session.
- GREP: keyword search over the full text only.
- RECALL: from memory, not verified this session.

## 1. Orr – Bechtloff Weising, 2410.13642v1 (17 Oct 2024) — READ (§2.3, §2.7, §3, §4 to Prop 4.2, §5.1 through Prop 5.4, §7.2)
The PDF lists the authors as "DANIEL ORR AND MILO BECHTLOFF WEISING".

### Setting (READ)
- The algebra is the CGM algebra B_{t,q}, acting on the polynomial representation V_k = Λ ⊗ K[y_1..y_k] (§2.7, p.7).
- They also use the Ion–Wu t-adic stable limit ("lim~_m", p.12) and BW23 stable-limit nonsymmetric Macdonald functions.
- This is stable-limit / almost-symmetric, not finite m.
- There is a q ↔ q^{-1} twist: E~(...; q^{-1}, t) appears in Prop 4.2. Their Demazure–Lusztig quadratic relation is (T_i − 1)(T_i + t) = 0 (§2, p.3).

### Def 2.3, p.3 (READ)
> "For 0 ≤ m ≤ n, we define the partial symmetrizers P⁺_m := Σ_{w∈S_m} T_w and ε^{(n)}_{n−m} := (1/[m]_t!) Σ_{w∈S_{(1^{n−m},m)}} t^{C(m,2)−ℓ(w)} T_w"

These are full symmetrizers over a Young subgroup, not sums over coset representatives.

### Lemma 4.1, p.8 (READ)
> "E_(γ1+1,…,γk+1,λ1,…,λm)(x1,…,xk,xk+1,…xk+m) = q^{γ1+···+γk} x1···xk E_(λ1,…,λm,γ1,…,γk)(xk+1,…,xk+m, q^{−1}x1,…,q^{−1}xk)"

The proof iterates Knop–Sahi.

### Prop 4.2, p.9 (READ, verbatim modulo pdftotext)
> "Proposition 4.2. For any λ ∈ Y and γ ∈ (Z≥0)^k,
> ∏_i ∏_{j=1}^{m_i(λ)} (1 − t^j)^{−1} E~_(γ1+1,...,γk+1|λ)(qx1,...,qxk,xk+1,xk+2,...; q^{−1},t) = q^{γ1+···+γk} x1···xk P_(λ|γ)(xk+1+···|x1,...,xk)."
>
> "Proof. Choose any m such that λ ∈ Y_m. Starting from Lemma 4.1, we apply the normalized symmetrizer P⁺_{[k+1,...,k+m]}/S^{(m)}_λ(t) to both sides and then take the limit as m → ∞."

- **Content:** it compares normalizations, Goodberry/Lapointe P_(λ|γ) against the BW stable-limit E~.
- It uses a full symmetrizer on the tail variables.
- It has no coset sum, no rational kernel, no Y-operators, no Pieri rule and no vertex operator.

### Prop 5.2 ([CGM20]), p.11 (READ)
> "T^{-1}_1···T^{-1}_{k−1} z_k F = F(X + (t−1)qy1 − (t−1)u | y2,…,yk,u) Ω(u^{-1}qy1 − u^{-1}X)|_{u^0}"

### Lemma 5.3, p.12 (READ)
> "T_{m−1}···T_k(x_k^b) = h_b(x_m + (1−t)(x_k+…+x_{m−1}))"

This is a single T-chain, i.e. a single coset representative, not the coset sum.

### Prop 5.4 ([IW22]), p.12 (READ)
> "The following result is implicit in [IW22] but we include its proof for the sake of completeness."
>
> "Y1 T1···Tk−1(x1···xk F(Xk|x1,…,xk)) = q t^k x1···xk F(Xk + qx1 − u | x2,…,xk,u) Ω(u^{−1}(qx1 + (1−t)Xk))|_{u^0}."

- This is a **single** Y_1 written as a Garsia–Jing-type plethystic operator. Ω[(1−t)X/u] is the HL Q-series.
- It is the closest item in the paper to the Jing identification, but only for k=1 and only in the stable limit.

### §7.2 "Pieri formulas", p.22 (READ, full text)
> "Via matrix coefficient calculations of [CGM20, Sec. 5], our Theorem 6.12 enables one to recover the Pieri formula for multiplication by e1(X) on partially-symmetric Macdonald polynomials which was proved in [Goo23] and connected to [CGM20] in [GO23]. One can also now obtain additional Pieri-type formulas for multiplication by y1, …, yk in this way."

- That is the entire section. It is an unproved remark: no Pieri formula is actually written down.
- The rules it mentions are X-side multiplication (e_1(X), y_i) in the Macdonald basis.
- There is nothing on e_k(Y) acting on e_r, and nothing on products e_a e_b.

### Kernel check (GREP)
- A search for "coset", "vertex", "Jing", "Hall-Littlewood" finds no hits beyond the above.
- The only symmetrizers in the paper are those of Def 2.3.
- There is **no** finite-m kernel ∏_{i∈A, j∉A}(X_i − tX_j)/(X_i − X_j), and no coset sum over S_m/(S_k×S_{m−k}).

### Does Rick's result follow? No (my inference)
- **Wrong operator family.** Their Y-content is one Y_1 (Prop 5.4) plus the spectrum (Cor 5.5). To get e_k(Y) you would expand it into a sum of ordered k-fold products of Y_i's on T-twisted inputs, then iterate Prop 5.4 k times through the T-chains. That k-fold composition with telescoping is exactly Rick's work. Nothing in the paper does it.
- **Wrong setting.** Theirs is the stable limit on V_k in CGM/Ion–Wu conventions. Rick's is Hikita's finite-m level-1 representation with π F = X_1 F(X_2.., sX_1). Matching the two would take a convention dictionary, which the paper does not supply. I did not construct one.
- **Wrong kind of Pieri.** §7.2's rules are e_1(X) and y_i multiplication on P_(λ|γ). These are the opposite side and carry no information about e_k(Y) acting on e_r.

**Verdict: PARTIAL OVERLAP.**
- The overlap is the k=1 "Y = HL-type plethystic vertex operator" (Prop 5.4 after IW22, Lemma 5.3, CGM d_− p.7 and Prop 5.2).
- There is no threat to the general-k result or to the k-fold coset-sum identification.

## 2. Shimozono – Zabrocki, math/0001168v1 — READ (§2–3 pp.3–4, §10 p.16); rest GREP
Verified quotes:
- (5), p.4: "H(Z^k)P[X] = P[X − (1−q)Z^*]Ω[ZX]R(Z^k)", with "R(Z^k) = ∏_{1≤i<j≤k} 1 − z_j/z_i".
- Remark 2.1, p.4: "If k = 1 this is Garsia's [2] [3] version of Jing's Hall-Littlewood vertex operator [6]."
- §10, p.16: "B(Z^k)P[X] = P[X − Z^*]Ω[XZ(1 − q)]R(Z^k) … coincide with the operators B(z) and B(z1,…,zk) in the notation of [9, Ex. III.5.8]. B(z) in our notation is the operator H(z) in Jing's [6]."
  - Their q is our t.
  - [9] = Macdonald, *Symmetric Functions and Hall Polynomials*.

Findings:
- No Hecke algebra, no Y, no e-basis Pieri.
- This is the symmetric-function side of the k-fold Jing operator only.

**Verdict: CLEAR.** Cite it as the source of the k-fold notation B(Z^k)/H(Z^k).

## 3. Venkateswaran, 2308.10844v3 ("Affine Hecke algebras and symmetric quasi-polynomial duality") — GREP plus the Thm 4.6 statement (p.25)
- §4 computes Bernstein-basis matrix coefficients of the full and J-partial (anti)symmetrizers.
- It has no Λ, no vertex operators and no e_k(Y) Pieri.

**Verdict: CLEAR.** Here I relied on the Sept 26 read plus grep; I did not re-read §4 in full.

## 4. Bhattacharya, 2407.14652v1 — READ (§1 pp.1–4); rest GREP
Verified quotes:
- (1.2), p.3: "1_{ϖℓ} is the parabolic Hecke symmetrizer Σ_{w∈S_{n,ϖℓ}} T_w".
- p.3: "u_C s are the minimal length coset representatives of S_n/S_{n,ϖℓ_{k+1}}".
- p.4: "1_0 X^λ = Σ_{T∈B(λ)} X^T Ψ_T, this is an affine Hecke algebra lift of Macdonald's formula (1.1)."

Findings:
- It uses the same column-indexed coset combinatorics, but on the X/Satake side to give P_λ(t).
- It has no Y, no Pieri and no vertex operators.

**Verdict: CLEAR (methodological cousin).**

## Important caveat (RECALL; verify before any novelty claim)
- The rational kernel Σ_{A} w_A( f · ∏_{i∈A, j∉A}(X_i − tX_j)/(X_i − X_j) ) is, taken alone, Macdonald's classical definition of HL P_λ (SF&HP Ch. III, eq. (2.2)). At λ = (1^k) it gives P_{(1^k)} = e_k.
- "Hecke coset symmetrizer = this rational symmetrizer" is standard Satake / Hecke-symmetrizer material, e.g. Macdonald's *Affine Hecke Algebras and Orthogonal Polynomials*, and Nelsen–Ram.
- Macdonald Ex. III.5.8 contains the k-fold B(z_1..z_k).
- Therefore (K_k) is **not new as a formula**. What could be new:
  - composing it with π^k at level 1 in Hikita's representation;
  - identifying σ^(k)π^k with the k-fold Jing operator at finite m;
  - the telescoping that lands in the e-basis.
- Rick should check Macdonald III.2 and Ex. III.5.8 for exact locators.

## Citation recommendation
- **OBW24** Prop 5.4 and Lemma 5.3, **IW22**, and **CGM20** Prop 5.2 / d_−: the known k=1 "Y = plethystic HL vertex operator" statement in the stable limit.
- **OBW24** §7.2: documents that the known Pieri rules in this setting are e_1(X) and y_i only.
- **SZ00** (5) and §10, **Jing 1991**, and **Macdonald** Ex. III.5.8: k-fold HL vertex operators.
- **Macdonald** Ch. III (2.2): the rational kernel formula itself.
- **Bhattacharya** and **Venkateswaran**: optional, one sentence each.

## Still open
Read Ion–Wu (IW22) itself for any multi-Y or coset-sum version. This is the only unread item that could upgrade the k-fold identification to a THREAT.

## 5. Ion–Wu (IW22), "The Stable Limit DAHA and the Double Dyck Path Algebra" — READ (§6.12–6.13, §7.4–7.8) plus full-text GREP
- **ID correction:** the paper is **arXiv 2011.12189v3**. The guessed ID 1807.04855 is an unrelated OCT/glaucoma paper.
- Local copies: `/home/agent/projects/papers/IonWu-2011.12189v3.{pdf,txt}`.
- I also fetched the follow-up, Ion–Wu **2504.03113v2** ("The stable limit DAHA: the structure of the standard representation"), and grepped it: `/home/agent/projects/papers/IonWu-2504.03113v2.{pdf,txt}`.

**What IW prove about Y (verified-quote):**
- **Prop 6.25, p.26–27.** The stable-limit Y_i is defined as a limit: "Y_i f = lim_k Ỹ_i^{(k)} Π_k f ∈ P^+_∞". The Ỹ are deformed Cherednik operators. There is a single Y_i and no closed form.
- **Lemma 6.28, p.27:** "P(k)^+ is stable under the action of Y_1."
- **Prop 6.32, p.28–29.** This is the only explicit formula, and it is for a single Y_1. For F = f(x_1..x_{k−1}) x_k^n G[X_{k−1}]:
  "Y_1 T_1 ··· T_{k−1} F = t^k/(1−t) · f(x_2,...,x_k) G[X_k + q x_1] (h_n[(1−t)(X_k + q x_1)] − h_n[(1−t)X_k])."
  - The proof goes through Lemma 6.31: "t^m T_m^{-1}...T_1^{-1} x_1^n = Σ_{i=0}^{n−1} x_{m+1}^{n−i} h_i[(1−t)X_m]".
- **§7.5–7.6, p.34–35.** The HL vertex operators appear only inside the lowering arrow ∂^−_k:
  - "∂_k^−(x_k^n F[X_k]) = B_n F[X_{k−1}]"
  - "B_n F[X] = (F[X − z^{−1}] Exp[−(t−1)zX])|_{z^n}"
  - "(modulo a change of variable) the vertex operators in [Jin91] (see also [Mac15, §III.5, Exp. 8])"
  - (7.5): "∂_k^− f F[X_k] = τ_k c_{x_k}(f F[X_k − x_k] Exp[−(t−1)x_k^{−1} X_k])".
- **Thm 7.13, p.36:** z_1 T_1···T_{k−1} = t^k/(1−t) [d^*_+, d^−] "acts on P(k)^+ as Y_1 T_1 ··· T_{k−1}". This is the single-Y "Y = commutator with a Jing operator" statement that OBW24 Prop 5.4 cites.

**Multi-Y content:**
- *agent-summary:* grep finds no products Y_1···Y_k, no e_k[Y] or Sym[Y], no coset sum over S_m/(S_k×S_{m−k}), no k-fold B(z_1..z_k), and no Pieri rule.
- Eigenvalues of Y_i and non-symmetric Macdonald eigenfunctions are discussed (§6.3). So is the PBW word u = X_1^{g_m} Y_1 X_1 … Y_1^z in 2504.03113 (8.16).
- 2504.03113 p.3 mentions "Pieri formulas" only as a strategy that *fails* for faithfulness. It proves no Pieri rule.

**Verdicts:**
- (a) General-k e_k⋆e_r: **CLEAR.** IW has only Y_1 and single Y_i, stable limit, with no Pieri rule of any kind.
- (b) k-fold Jing identification: **PARTIAL OVERLAP (k=1 only).** IW identify Y_1 with the commutator [d^*_+, d^−], where d^− is built from Jing's B_n and Exp[−(t−1)zX]. That is the k=1 shadow of σ^(k)π^k = k-fold Jing, but in the stable limit and not at finite m. Nothing is k-fold, and there is no telescoping on E(z).
- **Citation:** cite IW22 Prop 6.32 and Thm 7.13 (with §7.6) as the k=1 precedent, together with OBW24 Prop 5.4.
- "Still open" above is now CLOSED for IW.
