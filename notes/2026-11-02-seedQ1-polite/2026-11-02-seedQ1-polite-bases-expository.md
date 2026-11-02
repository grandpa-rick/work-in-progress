# Seed Q1/Q4, "Path 4": polite bases, coproducts, and the tensor-product rule

2026-11-02. Expository write-up of arXiv:2401.10340, read first-hand (pdftotext of v1, the only version), plus a test case.
Grades: hunch < sketched < computed < checked-sober < proved. **[P: locator]** means read first-hand in the paper today.
**[mem]** means a literature locator recalled from memory and NOT checked first-hand this session.

Files:
- PDF `reading/pdfs/2401.10340.pdf`; also `reading/pdfs/1905.08460.pdf` (BKK, read §2–3.1 first-hand for the biperfect material)
- scripts `proofs/seedQ1/sym_vs_crystal.py` (+ `.out`), `proofs/seedQ1/a2_polite_factorisation.py` (+ `.out`)

---

## 0. Bibliographic correction (read this first)

Copied from arxiv.org/abs/2401.10340 (citation_title / citation_author meta tags):

> **Title:** On Mirković–Vilonen polytopes
> **Author:** Baumann, Pierre (sole author). Submitted 2024-01-18, v1 only, 30 pages, math.RT.

Three things in the brief (and in `memory/questions/q-zero-CHA.md` l.76–83) are wrong when checked against the paper:
1. It is **not** "Baumann et al."; the only author is Pierre Baumann. "Polite bases" is a notion introduced in the paper, not its title.
2. The definition of polite basis is **Definition 3.1** (§3.3). There is no Def 3.3.
3. The paper works with the **classical** Hopf algebra O(N) (q = 1), not with U_q(n). And polite bases do **not** take crystal structure as an assumption: the crystal is an output (§2(b) below).

---

## 1. Setting (a)

**[P: §3.1]** G is split connected reductive over a field k of characteristic 0, with B ⊃ T, N = unipotent radical of B, n = Lie N, Chevalley generators e_i (i ∈ I), and U(n) graded by Q_+.
O = O(N) is the commutative Hopf algebra of regular functions on N, graded by Q_− through T-conjugation. Its coproduct is the **ordinary** one, dual to group multiplication. It is not twisted or braided: O(N) is commutative, and the coproduct is not cocommutative.
N acts on O from both sides: "n · f = f(?n) and f · n = f(n?)". The pairing is ⟨x, f⟩ = (x·f)(1_N), and "This pairing is perfect and identifies O(N) as a bialgebra to the graded dual of U(n)."
ζ_i ∈ O_{−α_i} is defined by ⟨e_i, ζ_i⟩ = 1. The involution ∗ is f*(n) = (−1)^{ht λ} f(n^{−1}).

*One-line lemma (proved).* With these conventions, x·f = (id ⊗ ⟨x,·⟩)∘Δ(f) and f·x = (⟨x,·⟩ ⊗ id)∘Δ(f) for x ∈ U(n).
Proof: (x·f)(g) = d/dt f(g e^{tx})|_{t=0} = Σ f_(1)(g)⟨x, f_(2)⟩, and higher products are handled by induction. So **the left and right U(n)-actions on O are components of Δ**. This is the classical shadow of Kashiwara's e_i′ (§4).

**Stability machinery [P: §2].** A is a connected bialgebra graded by a discrete submonoid Γ of the right half-plane.
- L(a) is the set of left degrees β with a nonzero component of Δ(a) in A_β ⊗ A_γ. Lemma 2.1: if Δa = Σ b_i ⊗ c_i is a minimal writing, then L(b_i) ⊂ L(a) and R(c_i) ⊂ R(a). The proof is coassociativity.
- Semistable of slope μ means L(a) ⊂ Γ_{≤μ}.
- Theorem 2.6: multiplication ⊗_{μ∈R̄} A_[μ] → A is bijective. This is a Harder–Narasimhan factorization.
- Prop 2.9: there is a filtration F′ (via the total order ≤′) and an algebra isomorphism
  Δ′ : gr′ A → A_[>0] ⊗ A_[≤0], which "sends a homogeneous element a ∈ gr′_α A to the homogeneous component of degree (α, |a| − α) of Δ(a)".
- Prop 2.10 is the mirror version Δ″ : gr″ A → A_[≥0] ⊗ A_[<0]. Prop 2.11 is the three-factor version Δ_2.

**[P: §3.2]** For θ ∈ Q_R^*, regrade by λ ↦ (ht λ, θ(λ)). Then U_[μ] = U(n_μ), where n_μ is spanned by root vectors of slope θ(α)/ht α = μ. N factors as the Zappa–Szép product N_{>0} ⋈ N_0 ⋈ N_{<0}, and
O_[≥0] = ^{N<0}O, O_[>0] = ^{N0 N<0}O, O_[≤0] = O^{N>0}, O_[<0] = O^{N>0 N0}, O_[0] = ^{N<0}O^{N>0}.

**Definition 3.1 (polite) [P: §3.3], verbatim:** "A basis B of O is said to be polite if it satisfies the following four conditions:
(i) The elements of B are weight vectors.
(ii) For each (i, n) ∈ I × N, the function ζ_i^n belongs to B.
(iii) For each θ ∈ Q_R^*, the basis B is compatible with the filtration ^θF′_• O and with the isomorphism defined in Proposition 2.9.
(iv) For each θ ∈ Q_R^*, the basis B is compatible with the filtration ^θF″_• O and with the isomorphism defined in Proposition 2.10."

Compatibility in (iii) means: "for each degree γ ∈ Γ_θ, the subspace ^θF′_γ O should be spanned by B ∩ ^θF′_γ O … Then Δ′ should map this basis of gr′_• O onto the tensor product basis of O_[≤0] ⊗ O_[>0]." The author adds: "one can show that conditions (iii) and (iv) are in fact equivalent".
Intro paraphrase [P: p.3]: B is polite iff each S(K) is spanned by B ∩ S(K) and, for each f ∈ B and each vertex μ of Pol(f), Δ(f) has exactly one term b_i ⊗ c_i with wt b_i = μ, with b_i, c_i ∈ B.

**Perfect / biperfect [P: §3.5], verbatim:** "a basis B of O consisting of weight vectors is said to be perfect if, for each (i, n) ∈ I × N, the rule b ↦ e_i^{(n)}·b defines an injective map {b ∈ B | ℓ_i(b) = n} → {b ∈ B | e_i·b = 0}." Here ℓ_i(f) is the smallest n with e_i^{n+1}·f = 0.
"We will say that a basis B of O is biperfect if both B and its image by the involution ∗ are perfect."
Baumann calls this "a strengthened version" of Berenstein–Kazhdan's definition [BeKa, Contemp. Math. 433 (2007)].
Compare BKK Def 2.1 [P: 1905.08460 §2.3]: there a basis is perfect if it is "endowed with an upper semi-normal crystal structure" satisfying e_i·b = ε_i(b) ẽ_i b + Σ_{ε_i(b′)<ε_i(b)−1} a_{b′} b′. Even there the crystal is determined by the basis, since ẽ_i b is the leading term of e_i·b.

Main results [P]:
- Prop 3.3: polite ⇒ every S(K) is spanned by B ∩ S(K).
- Cor 3.5: Pol(f) is GGMS.
- Cor 3.7: Pol is injective on B.
- Prop 3.9: transition matrices between polite bases are unitriangular.
- Prop 3.11: the dual semicanonical basis is polite.
- **Prop 3.12: "A polite basis of O is biperfect."**
- Example 3.2: the dual canonical basis is polite (G simply laced), via Varagnolo–Vasserot + Tingley–Webster.
- Thm 4.8: the MV basis is polite.
- Rem 3.13: biperfect does not imply polite (type A_4 example, after Fei).

## 2. Crystal structure: input or output? (b)

**Verdict (checked first-hand): OUTPUT.** Def 3.1 has no crystal datum. Its ingredients are weights, the elements ζ_i^n, and the filtrations and isomorphisms Δ′, Δ″, all built from Δ alone.
Crystal structure first appears after Prop 3.12. Proof of 3.12 [P: §3.5]: take θ with θ(α_i) < 0 and θ > 0 on the other positive roots. Then N_{<0} = G_a^{(α_i)}, O_[<0] = k[ζ_i], O_[≥0] = ker e_i, and for n = ℓ_i(f)

> "∆″ f̄ = e_i^{(n)}·f ⊗ ζ_i^n … the axioms of a polite basis require that B_0 is a basis of F″_0 O = O_[≥0], that {b̄ | b ∈ B_n} is a basis of gr″_{nγ} O, and that b̄ ↦ e_i^{(n)}·b maps bijectively this basis onto B_0."

So the polite axiom at this θ directly produces, out of the coproduct, the map b ↦ (e_i^{(ε_i(b))}·b, ε_i(b)) ∈ B_0 × N, and this map is a **bijection**.
That is the i-string decomposition, and it gives ẽ_i^{max} and hence ẽ_i, f̃_i: move n by ±1 in the fibre. Running the same argument with −θ gives the ∗-crystal.
Then [P: p.18–19]: "Perfect bases are automatically endowed with a structure of crystal in the sense of Kashiwara. We will not recall the definition … One important result from [8] is that the crystals of the perfect bases of O are all isomorphic. The abstract crystal they share is denoted by B(∞)." Also: "the set of Mirković–Vilonen polytopes has a natural crystal structure, isomorphic to B(∞). This result goes back to Kamnitzer".
The abstract says the extremal-point behaviour "extends the crystal operations".

Where outside input does enter:
- **Identifying** the crystal with Kashiwara's B(∞) uses BK [8] Thm 5.37 and BKK Thm 2.4 (biperfect bases have a unique bicrystal, proved by comparing transition matrices with a known biperfect basis). In the end this rests on the dual canonical basis, i.e. on Kashiwara's/Lusztig's U_q theory.
- **Existence** of polite bases uses KLR/Tingley–Webster (Ex 3.2), preprojective modules (Prop 3.11), or geometric Satake (Thm 4.8).

The memory note "It *assumes* Kashiwara's B(∞) operators" is therefore wrong. The single-crystal operators ẽ_i, ε_i come **out of** Δ.

## 3. Does the paper get the tensor-product rule from Δ? (c)

**No.**
- `grep -i tensor` on the full text gives only vector-space tensor products: "Δ′ should map this basis of gr′_• O onto the tensor product basis of O_[≤0] ⊗ O_[>0]" (§3.3), and "the tensor product is computed in the decreasing order" (§2.2).
- The signature rule, B ⊗ B′ as a crystal, and Kashiwara's embedding Ψ_i do not appear anywhere in the paper. The crystal is never defined in the paper ("We will not recall the definition").
- BKK (the main crystal reference used) is also tensor-rule-free in the parts I read. Remark 2.6(i) and Cor 2.8 express tensor-product **multiplicities** through ε_i(b) ≤ ⟨h_i,μ⟩, ε_i*(b) ≤ ⟨h_i,λ⟩ (Berenstein–Zelevinsky "good basis" style). That is multiplicity data, not the signature rule.

One structural observation (sketched, computed in A_2). For the special θ with {α : θ(α) > 0} = {α_i}, the polite bijection Φ_θ : B → B_{>0} × B_{≤0} = {ζ_i^n} × B_{≤0} from Def 3.1(iii) is the (mirror) **Kashiwara embedding** B(∞) ↪ B_i ⊗ B(∞), b ↦ (n = ε_i*(b), b·e_i^{(n)}).
Under it, the left crystal operators obey the signature rule. Computed in A_2 (chambers C4, C6 below):
- ε_j(b) = max(ε_j(b_−), ε_j(b_+) − ⟨h_j, wt b_−⟩) for 34/34 elements of height ≤ 5, with j = 1 in C4 and j = 2 in C6 **[corrected: was "ε_1" for both; the .out shows the 34/34 "reversed" match for j=1 in C4 and j=2 in C6]**;
- ẽ_j changes exactly one factor in 22/22 cases for each j.

That the signature rule holds here is textbook for the dual canonical basis: Kashiwara, Duke **71** (1993) 839–858, "The crystal base and Littelmann's refined Demazure character formula", Thm 2.2.1 **[corrected: Duke 69 (1993) is "Global crystal bases of quantum groups"; volume checked via Crossref + TW [Kas93] bibliography; theorem number still [mem]]**. It transfers to every biperfect basis (hence every polite basis) through BKK Thm 2.4 [P]. So inside Baumann's framework the **rank-one** tensor rule is present implicitly, but it is imported, not derived.

## 4. Test case: Sym versus the type-A tableau crystal (d)

Script `sym_vs_crystal.py`. Tableau crystal: LS bracket rule on row reading words (bottom row first), with every output checked to be a valid SSYT. Sym side: sympy bialternant Schur polynomials and leading-term peeling. All results below are **computed**.

**A1. Levi coproduct.** Take Sym(x_1..x_n) → Sym(x_1..x_k) ⊗ Sym(x_{k+1}..x_n), with s_λ ↦ Σ d^λ_{μν} s_μ ⊗ s_ν, for n = 2, 3, |λ| ≤ 4, all k: 28 cases. In all 28, d^λ_{μν} equals the number of highest-weight elements of B(λ) for the colours I∖{k} of weight (μ,ν).
So Δ is exactly the **Levi branching multiplicities** of the crystal. Two consequences:
- The colour sets of the two factors are disjoint ({1..k−1} and {k+1..n−1}). The "tensor product" B(μ) ⊗ B(ν) for gl_k × gl_{n−k} has no signature rule at all, because each colour acts on one factor only. **The variable-split coproduct never sees the signature rule.**
- The crystal edges of the removed colour k, which join different Levi components (e.g. 9 f_k-edges for λ = (3,1), n = 3, k = 1 or 2), are invisible to Δ. Δ only sees the Levi components.

What *is* determined (sketched, a short argument): suppose the vertex set carries the Levi crystal and the whole is known to be ≅ B(λ). Then the f_k-edges are determined up to ∏ (d^λ_{μν})! choices. Reason: extensions correspond to Levi isomorphisms, and Aut B(λ) = 1.
In all 28 cases ∏ d! = 1. These Levi branchings are multiplicity-free, which always holds for the one-variable split (Pieri / Gelfand–Tsetlin). But this "determination" rests on the uniqueness of B(λ), not on reading anything off Δ. And the vertex set itself is not an object of Sym: Sym = K_0, and vertices of the same weight are indistinguishable. For λ = (3,1), n = 3, three weight spaces have dimension ≥ 2.

**A2. Product.** For all 18 pairs μ, ν with |μ|, |ν| ≤ 2 and n = 2, 3, the number of highest-weight elements of B(μ) ⊗ B(ν) (signature rule on concatenated words) of weight λ equals c^λ_{μν} from s_μ s_ν. So the signature rule is *consistent with* the Hopf data.

**B. Does the Hopf data single out the signature rule? No (computed).** Take sl_2 and the vertex set B(m) × B(n). Call a structure "local" if each f-edge changes exactly one factor by one step. Count the normal sl_2-crystal structures that are local; any such structure automatically has the right character s_m s_n.

| m × n | local normal structures | signature conventions (Kashiwara, anti-Kashiwara) |
|---|---|---|
| 1×1, 1×2, 1×3, 2×1, 3×1 | 2 | both present, distinct |
| 2×2, 2×3, 3×2 | **6** | both present |

Hand check, B(1) ⊗ B(1), vertices 1⊗1, 1⊗2, 2⊗1, 2⊗2:
- s_1² = s_2 + s_{11} forces one 3-string and one singleton.
- The weight-0 space {1⊗2, 2⊗1} is 2-dimensional, and the coproduct/character cannot say which vector lies on the string.
- Kashiwara's rule (f acts on the first factor iff φ(b_1) > ε(b_2)) gives 1⊗1 → 2⊗1 → 2⊗2. The anti-Kashiwara rule gives 1⊗1 → 1⊗2 → 2⊗2.

So the multiplicity data in Δ (equivalently, in the product) cannot determine the rule.

**C. What extra structure suffices in rank one (computed).** Restrict to *uniform* rules: f(b⊗b′) = f b ⊗ b′ iff C(φ(b), ε(b′)), else b ⊗ f b′, with C an arbitrary boolean function on {0..3}².
Of the 65536 such rules, normality on all B(m) ⊗ B(n) with m, n ≤ 3 already leaves exactly 2. One is Kashiwara's C = [φ(b) > ε(b′)]. The other differs only in the cell (φ, ε) = (0, 3), where both choices give f = 0, so it is the same operator. Both survive strict associativity on triples.
So within this rule family, *normality plus uniformity* forces the signature rule. Whether associativity is needed in larger families, e.g. C(φ, ε, φ′, ε′), was not tested (2^{256} rules; it needs a smarter search).

**Honest summary for (d).** The coproduct of Sym sees **only multiplicities**: LR coefficients and Levi branching numbers. Sym is K_0 of the representation category, and its Schur basis indexes connected components (irreducibles), not crystal vertices. It cannot see crystal edges, and it does not single out the signature rule among local rules (B: 6 versus 2).
To see edges you need a Hopf algebra whose **basis elements are the crystal vertices**: O(N) = U(n)^* with a perfect, biperfect or polite basis at q = 1, or U_q(n)^* (the quantum shuffle algebra) with its dual canonical basis.
In O(N) the left action e_i· = (id ⊗ ⟨e_i,·⟩)∘Δ is a component of Δ, and §2 shows the polite axioms produce ε_i, ẽ_i from it. That is exactly where "coproduct ⇒ crystal operators" lives.

**A_2 computation (`a2_polite_factorisation.py`).** O = C[x,y,z], with Δz = z⊗1 + x⊗y + 1⊗z, e_1 = ∂_x, e_2 = ∂_y + x∂_z (right: ∂_x + y∂_z, ∂_y). B is the unique biperfect basis {x^a z^b w^c} ∪ {y^a z^b w^c}, w = xy − z (BKK Ex 2.3 [P]); it is polite by Prop 3.12 + uniqueness. 34 elements of height ≤ 5, one θ per Weyl chamber.
- Def 3.1(iii) is verified: the ≤′-extremal component of Δ(b) is a pure tensor b_+ ⊗ b_− of basis elements in every case.
- **Left-locality** (ẽ_j changes exactly one of b_+, b_−):

| chamber | {α > 0 : θ(α) > 0} | left j=1 | left j=2 | right* j=1 | right* j=2 |
|---|---|---|---|---|---|
| C1, C2 | all | 22/22 | 22/22 | 22/22 | 22/22 |
| C3 | {α_1, α_1+α_2} | **14/22** | 22/22 | 22/22 | 22/22 |
| C4 | {α_1} | 22/22 | 22/22 | 22/22 | **14/22** |
| C5 | {α_2, α_1+α_2} | 22/22 | **14/22** | 22/22 | 22/22 |
| C6 | {α_2} | 22/22 | 22/22 | **14/22** | 22/22 |
| C7, C8 | ∅ | 22/22 | 22/22 | 22/22 | 22/22 |

- A non-local example in C3: b = w with Φ(w) = (w, 1), ẽ_1 w = y, and Φ(y) = (1, y). Both factors change.
- The left pattern matches, in all 6 chambers (8 sample θ; C1/C2 and C7/C8 lie in the same chamber) **[corrected: was "all 8 chambers"]**, the criterion **LOC(θ, j): no positive root α with θ(α) ≤ 0 < θ(α + α_j)**, i.e. ad e_j preserves n_{≤0}(θ). The right* pattern matches the mirror criterion.
- Mechanism (sketched): B_{≤0} = B ∩ O^{N>0} is a left U(n)-submodule and hence ẽ_j-stable. B_{>0} = B ∩ ^{N≤0}O is stable under left e_j iff e_j normalizes n_{≤0}.

## 5. Prior art: is "crystal from coproduct" textbook? (e)

1. **Single-crystal operators from a coproduct component: textbook.**
   - Kashiwara, Duke 63 (1991) §3 **[mem]**: defines e_i′ on U_q^−(g) by [e_i, P] = (t_i e_i″P − t_i^{−1} e_i′P)/(q_i − q_i^{−1}), decomposes P = Σ f_i^{(n)} P_n with e_i′P_n = 0, and sets ẽ_i, f̃_i on (L(∞), B(∞)) by shifting n.
   - Lusztig, *Introduction to quantum groups* §1.2.13 **[mem]**: defines r_i, _ir by r(x) = r_i(x) ⊗ θ_i + … + θ_i ⊗ _ir(x), where r is the twisted (braided) coproduct of f. Kashiwara's e_i′ is _ir or r_i up to normalization (Lusztig 3.1.6 **[mem]**).
   - So "ẽ_i from e_i′ = a component of Δ" **is** the textbook answer to "where do crystal operators come from", in the quantum, braided-coproduct setting.
   - The classical O(N) version is the one-line lemma of §1 together with Berenstein–Kazhdan perfect bases [BeKa 2007, Def 5.30 / Thm 5.37, as cited by BKK [P]]. Baumann's Prop 3.12 [P] packages it as a consequence of coproduct compatibility.
2. **The tensor-product (signature) rule from a coproduct: also textbook.** But the coproduct is Δ of the **whole** U_q(g) acting on M ⊗ N, not the coproduct of O(N) or of Sym.
   - Kashiwara, CMP 133 (1990) "Crystalizing the q-analogue…", Thm 1 **[mem]**: if (L_j, B_j) are crystal bases of M_j, then (L_1 ⊗ L_2, B_1 ⊗ B_2) is a crystal base of M_1 ⊗ M_2 with the tensor rule. The proof is the U_q(sl_2) computation with Δ(e_i) = e_i ⊗ t_i^{−1} + 1 ⊗ e_i at q → 0.
   - Textbook treatments **[mem]**: Hong–Kang, *Introduction to Quantum Groups and Crystal Bases*, Thm 4.4.1 (or nearby); Bump–Schilling ch. 2 takes it as a definition.
   - The rank-one embedding B(∞) ↪ B(∞) ⊗ B_i: Kashiwara, Duke **71** (1993) Thm 2.2.1 **[corrected: volume 71, not 69; theorem number still [mem]]**.
   - Berenstein–Kazhdan geometric/unipotent crystals [BeKa, Contemp. Math 433; P: cited as [8]] give a tropicalization route, product of unipotent bicrystals ↦ tensor product. That is another place where a group-multiplication (= coproduct of O) structure produces tensor products of crystals **[mem, hunch-level recall]**.
3. **What is then genuinely open, with this paper as the base:**
   - Baumann's Δ′ factorizations Φ_θ exist for **every** θ, and recover Kashiwara's Ψ_i only when {θ > 0} is a single simple root (or its mirror).
   - Neither paper (nor BKK §2–3.1) asks how the crystal operators behave under Φ_θ for general θ.
   - Possible prior owner: Kamnitzer, Adv. Math 215 (2007) [cited as [13]] describes f̃_i on MV polytopes in terms of which vertices/Lusztig data move. Combined with Pol(b_+) ↔ the "part of Pol(b) up to the θ-vertex" (via Lemma 2.1 / Prop 2.9), Kamnitzer's theorem may already decide the locality question below. **This is the first thing to check** before attempting a proof.

**Verdict on the background question:** the textbook (Kashiwara 1990/1991, Lusztig ch. 1/3 **[mem]**) explains "coproduct → crystal" at both levels:
- e_i′ is a component of the (braided) coproduct of U_q^−, giving ẽ_i;
- Δ of U_q(g) on tensor products at q = 0 gives the signature rule.

The Sym coproduct is the decategorified shadow of these, and provably cannot see the rule (§4 B). The open residue is not "why"; it is the behaviour of crystal operators under the coproduct-defined factorizations at **general** stability parameters, which the polite-basis framework makes available.

---

## 6. ONE question (for a bounded PROVE session)

**Q60 (A_3 polite-factorization locality).** Let G = SL_4 and O = O(N). Let B be the polite basis of O; in type A_3 it is the unique biperfect basis (BKK §2.3), which by Ex 3.2 and Prop 3.12 equals the dual canonical basis at q = 1.

Let θ ∈ Q_R^* be generic: θ(α) ≠ 0 for all roots and θ injective on Q. For b ∈ B let μ_θ(b) be the weight ν ∈ L(b) maximizing θ (the ≤′-top of L_θ(b)). Let Φ_θ(b) = (b_+, b_−) ∈ (B ∩ O_[>0]) × (B ∩ O_[≤0]) be defined by: the component of Δ(b) in bidegree (μ_θ(b), wt b − μ_θ(b)) equals b_+ ⊗ b_− (Def 3.1(iii), Prop 2.9).

Call Φ_θ *ẽ_j-local* if for every b ∈ B with ẽ_j b ≠ 0, exactly one of the two coordinates of Φ_θ(ẽ_j b) differs from the corresponding coordinate of Φ_θ(b), where ẽ_j is the left Kashiwara operator of the perfect basis B.

**Is it true that for every generic θ and every j ∈ I,
Φ_θ is ẽ_j-local ⟺ there is no positive root α with θ(α) ≤ 0 < θ(α + α_j)?**

Evidence: computed true in all 6 chambers of A_2 (8 sample θ, two chambers sampled twice) **[corrected: was "all 8 chambers"]**, for both j (§4). The mechanism for "⇐" on the B_{≤0} side is sketched. A_3 is the first rank where {θ > 0} can be a biconvex set other than the inversion set of a single simple reflection or of its complement, and there the "⇒" direction can fail or hold. A negative answer is as informative as a positive one.

**Why this is not already textbook.**
1. The textbook tensor rules (Kashiwara CMP 1990 Thm 1; Duke 1993 Thm 2.2.1 **[mem]**) concern Δ(U_q(g)) on module tensor products and the rank-one embeddings Ψ_i. In Baumann's language these are only the θ with {θ > 0} = {α_i} (computed: chambers C4/C6 reproduce Ψ_i with the Kashiwara formula, 34/34).
2. Baumann 2401.10340 [P] defines Φ_θ for all θ but never discusses crystal operators under it ("We will not recall the definition [of crystal]"; no tensor rule anywhere in the text).
3. BKK 1905.08460 §2–3.1 [P] treats only ε_i, ε_i* and Saito reflections (Thm 3.2), again rank one.
4. **Live risk of ownership:** Kamnitzer 2007 (MV polytope crystal: f̃_i moves only certain vertices) and Tingley–Webster 2016 (KLR polytopes; θ-cuts = cuspidal/semicuspidal decompositions, where "f_i acts on one factor" statements exist for convex orders, McNamara/Kleshchev–Ram style) **[mem]**. Before attempting a proof, the PROVE session should read Kamnitzer 2007 §3 and Tingley–Webster §2 for a statement of the form "f̃_i changes only the θ-negative part of the Lusztig datum / of the semicuspidal factorization". If one exists, Q60 is owned and the answer to the brief is: textbook = Kashiwara 1990 + Kamnitzer 2007.

---

## Audit + owner check (2026-11-02)

Independent audit, same day. Grades: **computed** (by script), **sketched** (argument written but not checked line by line), **cited-with-locator** (read first-hand in the PDF unless marked otherwise). In-place fixes above are tagged **[corrected]**.

### A1. Where "8 chambers" came from (computed)

- In A_2, θ ∈ R^I with I = {1,2}. The 3 hyperplanes θ(α) = 0 cut the plane into **6** regions, the Weyl chambers. Nothing else is enumerated: the script has no extra hyperplanes and uses no other parameter space.
- `a2_polite_factorisation.py` hard-codes **8 sample θ**, named C1–C8. C1 = (1, .37) and C2 = (.41, 1) both lie in the dominant chamber {θ > 0} = Δ_+. C7 and C8 both lie in the antidominant chamber {θ > 0} = ∅.
- The addendum prints the same {θ>0} set for C1/C2 and for C7/C8, and every count agrees within each pair. So these are duplicate samples. The likely cause is the "(+,+)" sign labels: they come from the signs of (θ(α_1), θ(α_2)), and the sign-quadrant split was only refined by the sign of θ(α_1+α_2) in the mixed quadrants.
- Is it legitimate to reduce to the chamber? Yes. μ_θ(b) = argmax_{L(b)} θ is a vertex of Pol(b) = conv L(b). The normal fan of an MV/GGMS polytope is coarser than the Weyl fan. So Φ_θ depends only on the chamber of θ. The duplicate counts are consistent with this.
- **Correct count: 6 chambers, of which 2 are trivial.** In the dominant chamber Φ_θ(b) = (b, 1) ("on b+: 22/22"). In the antidominant chamber Φ_θ(b) = (1, b). In both, locality is automatic.
- Real evidence: **4 nontrivial chambers × 2 colours = 8 (θ, j) tests.** The criterion predicts NON-local in 2 of them, C3/j=1 and C5/j=2, and the script finds exactly those 2 non-local (14/22 local, 8/22 not).
- **Negative control: present and genuine.** The script really breaks locality. Example: C3, b = w = xy−z, Φ(w) = (w, 1), ẽ_1 w = y, Φ(y) = (1, y). Both coordinates change. I re-read the printout and confirm 8/22 failures in exactly C3/j=1 and C5/j=2.
- The right* table shows the mirror pattern, with failures at C4/j=2 and C6/j=1. **Caveat:** the "mirror criterion" for right* is not coded in the script. I checked it by hand ("no α > 0 with θ(α) > 0 ≥ θ(α+α_j)"), and it matches all 4 nontrivial chambers.
- Genericity of the sample θ: no ties occur for argmax over weights of height ≤ 5, because 0.37, 0.41, 0.43, 0.39, 0.45, 0.35, 0.33 and 0.31 satisfy no relation p·t = q with |p|, |q| ≤ 5. The polite pure-tensor check passed for all 8 samples, so no ties went unnoticed.

### A2. Is "ẽ_j-local" precise, and does the script test it? (computed)

- **Q60's definition is the weak form:** exactly one coordinate of Φ_θ changes. The addendum tests exactly this, as `(tp == bp) != (tm == bm)`.
- The main `run()` tests a **stronger** form: Φ_θ(ẽ_j b) ∈ {(ẽ_j b_+, b_−), (b_+, ẽ_j b_−)}. Its counts are the same (22/22 or 14/22) in every chamber, so in A_2 the weak and strong forms agree.
- **Recommendation:** state Q60 in the strong form. The weak form is not obviously closed: ẽ_j b_+ need not lie in B_{>0} when α_j has θ(α_j) < 0. The strong form also makes the A_3 check sharper.
- Precision issue, harmless here: §6 defines μ_θ(b) as the θ-maximiser of L(b), with weights taken in Q_+. Baumann grades O by Q_−. The script's convention gives Φ_θ(b) = (b, 1) in the dominant chamber, which fits O_[>0] = O there. So the convention is internally consistent, but the write-up should say "degrees in Q_+".
- Minor fix in §3, applied in place: the 34/34 signature-rule match is ε_1 in C4 and **ε_2** in C6.
- No bug in the computation itself.

### A3. A cleaner form of the criterion (computed + sketched)

Write S = {α > 0 : θ(α) > 0} and T = Δ_+ ∖ S. Both are biconvex. Then

**LOC(θ, j) ⟺ α_j ∈ T, or T ∪ {α_j} is biconvex.**

- *Sketch.* If α_j ∈ T and α ∈ T, then α + α_j ∈ T because T is closed, so the criterion cannot fail. If α_j ∈ S, the set T ∪ {α_j} is automatically co-closed, since a simple root is never a sum of two positive roots. It is closed exactly when there is no α ∈ T with α + α_j ∈ S, which is the criterion.
- In Weyl-group language, take T = N(w^{-1}). Then LOC ⟺ α_j ∈ N(w^{-1}), or w^{-1}α_j is a simple root.
- Script `seedQ1/q60_criterion_combinatorics.py` (+ `.out`) checks LOC ⟺ ALT over all chambers:
  - A_2: 12/12. Case split neg / adj / far = 6 / 4 / 2.
  - A_3: 72/72 over all 24 chambers. Case split = 36 / 18 / 18.
  - Here "neg" = (α_j ∈ T), "adj" = (α_j ∈ S and T ∪ {α_j} biconvex), "far" = the predicted non-local cases.

### B. Owner check (cited-with-locator, all PDFs first-hand in `reading/pdfs/`)

Titles and authors below are copied from the arxiv.org abs pages.

**Two ID corrections to the brief:**
1. arXiv **math/0501063** is *not* Kamnitzer. It is Dibaei–Yassemi, "Attached primes of the top local cohomology modules with respect to an ideal (II)".
2. Kamnitzer's "Mirkovic-Vilonen cycles and polytopes" (Kamnitzer, Joel) is **math/0501365**, published as Ann. of Math. 171 (2010) 245–294 (Crossref). Saved as `math_0501365.pdf`.

The other three IDs check out:
- math/0505398: "The crystal structure on the set of Mirkovic-Vilonen polytopes", Kamnitzer, Joel. Published Adv. Math. 215 (2007) 66–93 (Crossref).
- 1210.6921: "Mirkovic-Vilonen polytopes and Khovanov-Lauda-Rouquier algebras", Tingley, Peter; Webster, Ben.
- 1110.3661: "Affine Mirković-Vilonen polytopes", Baumann, Pierre; Kamnitzer, Joel; Tingley, Peter.

**What each paper states:**

- **Kamnitzer math/0505398.**
  - Prop 3.4, attributed to [L2] and [BZ3, Prop 3.6], verbatim: "Let j ∈ I and let i be a reduced word with i_1 = j. … (i) f_j·b has i-Lusztig datum (n_1+1, n_2, …, n_m). (ii) If n_1 = 0, then e_j·b = 0. Otherwise, e_j·b has i-Lusztig datum (n_1−1, n_2, …, n_m)."
  - Thm 3.5: "f_j·P is the unique MV polytope whose set of vertices μ′_• satisfies μ′_e = μ_e − α_j^∨ and μ′_w = μ_w if s_j w < w." The proof does this one vertex μ_w at a time, using a reduced word with i_1 = j that passes through w.
  - Nothing is said about a θ-cut factorization or about the case s_j w > w.
- **Tingley–Webster 1210.6921.**
  - Prop 3.13: "To apply the operator f̃_i to P̃ ∈ P_KLR, choose a convex order with α_i lowest, and read the path … L_1, …, L_h corresponding to increasing roots … If L_h = L_i^k for some k ≥ 1, then f̃_i P̃ = P̃_{A(L_1,…,L_{h−1},L_i^{k−1})}. If L_h ≇ L_i^k, then f̃_i P̃ = 0."
  - Prop 3.4 (Res_F L is simple for each face F) and Prop 3.9 (L = A(L_1, …, L_h) along the path P^≻) give the cut factorization.
  - Cor 2.26 and Cor 3.14 cover Saito reflections, with α_i *greatest*.
  - §3.5 (proof of Thm C): "f̃_i acts in a simple way on the Lusztig datum for any convex order ≻ with α_i minimal."
  - Prop 3.22 / Def 3.21 define face crystals: operators of g_c, where c is a charge, for the simple roots of the face root system.
  - There is **no** statement for α_i in a non-extremal position.
- **BKT 1110.3661.**
  - Thm 4.4: the torsion-pair bijection Ξ : T × F → B. This is the categorical analogue of Φ_θ.
  - Prop 4.5 proves the three-factor compatibility.
  - §1.4 Remark (iii): "if α_i is the smallest element of the order ≼, then φ_i(b) is the α_i-coordinate of Ω_≼(b), and the operators ẽ_i and f̃_i act by incrementing or decrementing this coordinate."
  - §5.5: for w = s_i, Ξ becomes (Z, n) ↦ ẽ_i^n Z, the rank-one Kashiwara-embedding shape.
  - Prop 5.24 gives Lusztig data via Saito reflections.
  - No locality statement for general (θ, j).
- **Kamnitzer math/0501365.** MV cycles and GGMS/BZ data. Its only mention of crystals points to [K] = math/0505398. Nothing relevant.

**Verdict: Q60 is (ii) for one half and (iii) for the rest.**

**(ii), the "neg" half: α_j ∈ T ⇒ strong locality, with ẽ_j acting on b_− only.** This is an easy corollary.
- Locators: Kamnitzer 0505398 Prop 3.4 (= BZ Prop 3.6) together with TW Props 3.4 and 3.9 and Baumann Ex 3.2. Alternatively, TW Prop 3.13 / BKT Rem 1.4(iii).
- Three-line deduction (sketched):
  1. Since T is an inversion set and α_j ∈ T, some convex order lists T first (T is an initial segment) and starts with α_j, because s_j is a descent.
  2. Via Baumann Ex 3.2 (dual canonical basis = KLR simples, through VV+TW) and TW Props 3.4 and 3.9, Φ_θ(b) = (b_+, b_−) is the grouping of the semicuspidal decomposition, or of the Lusztig datum, into its S-part and T-part.
  3. By Prop 3.4 / TW 3.13, ẽ_j changes only the α_j-coordinate, which lies in the T-part. So b_+ is unchanged and b_− ↦ ẽ_j b_−.
- This matches the data: "on b−: 22/22" in every A_2 neg case (C3/j=2, C4/j=2, C5/j=1, C6/j=1, C7, C8).
- **Orientation caveat:** the matching of "T-part ↔ b_−" with TW's "lowest = last factor" was not traced through the conventions. The A_2 data are what fix it.

**(iii), not stated in any of the four papers:**
- the "adj" half: α_j ∈ S and T ∪ {α_j} biconvex ⇒ local. In A_2 this behaves like the signature rule (C4/j=1: 14 cases act on b_+ and 8 on b_−). This is the Kashiwara-embedding shape, twisted by the cut;
- the "⇒" half: "far" ⇒ some b is non-local.

**Not owned by anyone.** No paper discusses ẽ_j for α_j in a non-extremal position relative to a cut. TW and BKT always choose an order with α_i extremal.

### C. Textbook locators

All four locators below are **unverified first-hand**. Project Euclid and Springer return HTML (Incapsula / paywall), not PDFs.

- **Kashiwara, Duke 63 (1991) §3** (e_i′, ẽ_i): bibliographic data confirmed via Crossref. It is "On crystal bases of the Q-analogue of universal enveloping algebras", Duke 63 (1991) 465–516. The section number is unverified.
- **Kashiwara, CMP 133 (1990) Thm 1** (tensor rule): confirmed via Crossref as "Crystalizing the q-analogue of universal enveloping algebras", CMP 133 (1990) 249–260. The theorem number is unverified.
- **Kashiwara, "Duke 69 (1993) Thm 2.2.1"** (B(∞) ↪ B(∞) ⊗ B_i): **wrong volume, corrected in place.**
  - Duke 69 (1993) 455–485 is "Global crystal bases of quantum groups".
  - The embedding paper is "The crystal base and Littelmann's refined Demazure character formula", **Duke 71 (1993) 839–858**. This is confirmed by Crossref and by the bibliographies of TW ([Kas93], cited for the involution at "[Kas93, 2.1.1]") and Kamnitzer 0505398 ([Kas]).
  - "Thm 2.2.1" remains [mem]. It is plausible, since §2.1 = involution per TW.
- Secondary first-hand support for the embedding: TW Prop 1.4's proof writes the map b ↦ (f̃_i^*)^{φ_i^*(b)} b ⊗ ẽ_i^{φ_i^*(b)} b_i, citing [KS97, Prop 3.2.3].

### D. Recommended statement for the PROVE session

**Q60′ (strong form).** Let G be simply laced of finite type, B the dual canonical basis at q = 1 (polite by Baumann Ex 3.2), and θ generic. Set S = {α > 0 : θ(α) > 0} and T = Δ_+ ∖ S. For j ∈ I:

- (a) *[(ii), owned in substance]* If α_j ∈ T, then Φ_θ(ẽ_j b) = (b_+, ẽ_j b_−) for all b.
- (b) *[open]* If α_j ∈ S and T ∪ {α_j} is biconvex, then Φ_θ(ẽ_j b) ∈ {(ẽ_j b_+, b_−), (b_+, ẽ_j b_−)} for all b. The choice between the two is decided by the Kashiwara signature rule applied to the pair (b_−-part restricted along α_j, α_j-string of b_+). Formulating this exactly is part of the task.
- (c) *[open]* If α_j ∈ S and T ∪ {α_j} is not biconvex, then some b has both coordinates changed.

**Suggested route.**
- Write (a) up from Kamnitzer Prop 3.4 + TW Props 3.4 and 3.9. This takes about one page and should first fix the orientation convention.
- For (b), use the three-factor decomposition T | α_j | S′ with S′ = S ∖ {α_j} (Baumann Prop 2.11 / BKT Prop 4.5). This reduces (b) to a two-piece statement about B_{T∪α_j} ≅ B_T × N. That looks like Kashiwara's embedding conjugated by Saito reflections (TW Cor 2.26 / BKT Prop 5.24), and is the likely proof.
- For (c), A_2 already shows the mechanism: b = w in C3. A single explicit witness per "far" pair, together with a reduction to a rank-2 Levi, may suffice.
- **Whether to compute A_3 first: yes, but keep it small.** Nothing in this audit tests (b) or (c) beyond A_2: there are 18 "adj" and 18 "far" (θ, j) pairs in A_3 that have never been tested.
  - The cheapest test is in Lusztig-data coordinates, not via polynomials in O(N_4) with 6 variables. Use Lusztig's piecewise-linear transition maps for A_3 (braid moves of type 121↔212 and commutations).
  - Concretely: ẽ_j = conjugate n_1 ↦ n_1 − 1 to a word adapted to the cut, and compare the S-part and T-part before and after. Do this for heights ≤ 6 or so.
  - That is about one hour of code and settles (b)/(c) empirically before any proof is attempted. It relies on (a)'s identification of Φ_θ with the Lusztig-datum split, which is the sketched TW 3.4/3.9 step.
