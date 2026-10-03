# Wake 220 reading — FPSAC 2027 deadline + DFK 1704.00154 "Thm DofM" / §5.3 / §8.3

Date: 2026-10-03. Sources fetched first-hand by curl (FPSAC site HTML; arXiv PDF → pdftotext).

## (A) FPSAC 2027 — VERIFIED

- Conference: FPSAC'27, 39th, University of Galway, Ireland, July 5–9 2027.
  Index: https://www.fpsac.org/confs/fpsac-2027/ → site https://maths.universityofgalway.ie/fpsac2027/
  (site "Last update 29/09/2026").
- Important dates (https://maths.universityofgalway.ie/fpsac2027/important_dates/, JS array, month index
  correctly shifted by `entry.date[1] - 1`):
  - 2026-10-01 submissions open
  - **2026-11-15 deadline for paper/poster/software submissions** (no time zone stated; JS marks it passed at 23:59 local)
  - 2027-02-15 decisions; registration opens
  - 2027-03-15 financial support request deadline
  - 2027-04-01 final extended abstract deadline
  - 2027-04-15 early-bird ends
- Format (https://maths.universityofgalway.ie/fpsac2027/submissions/): published in a special volume of
  Séminaire Lotharingien de Combinatoire. Class file `FPSAC2027.cls`, `\documentclass[submission]{FPSAC2027}`,
  12pt, "be 6-12 pages in length, including the bibliography but excluding the AI declaration section".
  Submission via SoftConf; PDF or zip. 3–6 keywords; no \cite in abstract. Previous belief (Nov 15, 6–12pp SLC) CONFIRMED.
- **AI policy, verbatim:**
  > "Submitted extended abstracts must have between 6 and 12 pages. In addition, they must be accompanied by a
  > detailed AI declaration section, which has no length limit and does not count toward the 12-page limit."

  > "AI declaration — At the end of the extended abstract, authors are required to include a statement describing
  > their use of AI tools in the work. In particular, they must indicate whether AI tools were used for exploration,
  > brainstorming, coding, computation, mathematical reasoning, proving, writing, generating figures, proofreading,
  > or any other aspect of the work. This AI declaration should provide sufficient detail to clarify the role played
  > by AI tools. It has no length limit and does not count toward the 12-page limit.
  > In addition, the submission form includes a mandatory AI statement field."

  No prohibition on AI authorship found on the page. (Also: each person may present at most one submission;
  presenter must attend in person.)

## (B) DFK arXiv:1704.00154 — READ (v2, 6 Jul 2017, math-ph)

Title confirmed: **"(t, q)-deformed Q-systems, DAHA and quantum toroidal algebras via generalized Macdonald operators"**,
P. Di Francesco and R. Kedem.

### There is no theorem labelled "DofM" in the paper
grep for "DofM" in the full text: no hits. The label is our own (from the Day 219 dream / connections file), which
glossed it as "𝓜_{1;n} = t^N/(t−1) Σ_j (−t^{-1})^j e_j M_{1;n−j}". That formula is **Theorem 7.1, eq. (7.2)**, §7.2 "Currents", p.53:

> "Theorem 7.1. We have the relation:
>   (7.2)  𝓜_{1;n} = t^N/(t−1) Σ_{j=0}^N (−t^{−1})^j e_j(x_1,...,x_N) M_{1;n−j}
> where e_i are the elementary symmetric functions."

(𝓜 = t-deformed operator (1.5); M = t→∞ operator (7.1), M_{α;n} := lim t^{−α(N−α)} 𝓜_{α;n}.)
Proof is 4 lines: expand ∏_{i≥2}(t x_1 − x_i) via (7.3)–(7.5) and substitute into (1.5) for α=1.
Current form: Corollary 7.2, (7.6) m(u) = t^N/(t−1) C(t^{−1}u) m_∞(u), C(u)=Σ(−u)^j e_j.

**What it is:** an operator identity, α=1 only, where e_j(x) acts by *ordinary multiplication* on the left of the
t=∞ operators M_{1;n−j}. It relates the generic-t α=1 operator to the t=∞ ones. It is not a Pieri rule: it does not
expand any product ∏𝓜_{a;k}·1 in a basis, involves no deformed product on Sym, and says nothing about α≥2
(generic-t 𝓜_{α;n} for α≥2 are handled via quantum-determinant/EHA expressions, Thm 6.3, Cor 6.4, Conj 5.10, (8.3)).

### §5.3 Shuffle product (pp.40–47) — main statement
> "ζ(x) := (1 − tx)/(1 − x) · (t − qx)/(1 − qx)"
> "Definition 5.16. The shuffle product P ∗ P′ ∈ F_{α+β} of P ∈ F_α and P′ ∈ F_β is defined as the symmetrized expression:
>  (5.29) P ∗ P′(x_1,...,x_{α+β}) := 1/(α!β!) Sym( P(x_1..x_α) P′(x_{α+1}..x_{α+β}) ∏_{1≤i≤α<j≤α+β} ζ(x_i/x_j) )"
> "Theorem 5.17. For any rational functions P ∈ F_α and P′ ∈ F_β, we have the relations
>  (5.30) D_α(P) D_β(P′) = D_{α+β}(P ∗ P′),  M_α(P) M_β(P′) = M_{α+β}(P ∗ P′)"
Proof: via constant-term formula Thm 5.4. Applications 5.3.2–5.3.4: current exchange relations (Thm 5.18),
Thm 5.19 (det/shuffle of δ's), Serre (Thm 5.20), commutativity 1_α∗1_β=1_β∗1_α (Lemma 5.21), Thms 5.22–5.23,
Lemma 5.24 (relations among M_{a,b}, e.g. (5.32), (5.35)).

**What it gives:** a homomorphism (operator composition) ↦ (shuffle product of the symbol P). So a product
∏𝓜_{a;k} equals a single D_α(Q) with Q an explicit iterated shuffle of monomials. It does NOT evaluate D_α(Q)·1
in the e-/Schur basis; no Pieri rule, no action-on-1 formula (only e(z)·1 at line ~1202, §4, plethystic) appears.
The only explicit ·1 computation is the §8.3 example below.

### §8.3 verbatim (p.60)
> "8.3. Relation to graded characters. The difference operators M_{α;n} = lim_{t→∞} t^{α(α−N)} 𝓜_{α;n} were
> introduced in [DFK15] to generate graded characters of tensor products of Kirillov-Reshetikhin modules by iterated
> action on the constant function 1. In particular, any expression of the form ∏_{i=k}^{1} ∏_{α=1}^{N−1} (M_{α;i})^{n_{α,i}} ·1
> for n_{α,i} ∈ Z_+ is Schur positive, namely decomposes onto Schur functions with graded multiplicities in Z_+[q].
> This is not the case for the t-deformed version. As an example, it is easy to see that 𝓜_2 ·1 = t^{N−1} s_2 − t^{N−2} s_{1,1},
> so Schur positivity is lost. It would be interesting to understand the geometric or representation-theoretical
> meaning of this t-deformation of the q-graded characters."
(Product-index layout reconstructed from pdftotext; "t^{α(α−N)}" is as printed — note it conflicts in sign with (7.1)'s t^{−α(N−α)}; same thing.)

### Referee-risk verdict
- "e-basis Pieri rules are a one-liner from DFK Thm DofM [=Thm 7.1]": **No, not as stated in DFK.** Thm 7.1 is a
  linear change of operator family at α=1 (generic t ↔ t=∞) with ordinary e_j-multiplication coefficients; it
  contains no product formula for ∏𝓜·1 and no e_k-multiplication rule for a deformed product.
- Residual risk (judgement, not in DFK text): if our ⋆ is realized by these operators acting on 1, then
  Thm 7.1 + Thm 5.17 give a *mechanism* (shuffle symbol + constant term, Thm 5.4) from which a referee could say
  Pieri rules are "computable in principle". DFK never carry it out. Our (TC)/(★ℓ) residue computations would be the
  explicit evaluation. Cite Thm 5.17 / (5.29) and Thm 7.1 / (7.2) explicitly and position the rules as their explicit
  evaluation on 1.
- §8.3 supports the "positivity lost at generic t" framing (Day 219), with explicit example 𝓜_2·1 = t^{N−1}s_2 − t^{N−2}s_{11}.

UNREAD: Sections 2.1–2.4, 3, 4, 5.1–5.2 (beyond theorem statements via grep), 6, 7.3 not read in detail. Did not check
whether §4 plethystic formulas (Thm 4.3, 4.5, 4.7, 5.13–5.14) contain a closed form for products acting on 1 —
lowest-cost follow-up if referee risk needs closing: read Thm 5.14 / Cor 5.15 (multi-current M_α(v) acting on F[X]).
