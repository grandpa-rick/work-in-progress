# FPSAC 2027 extended abstract: OUTLINE v7 (2026-10-04)

Anchor: `questions/q-fpsac-2027-writeup.md` §v7. Deadline **2026-11-15** (verified Wake 220). Format: 6–12 pp SLC, `FPSAC2027.cls`, plus a mandatory AI declaration that does not count toward the page limit.
**Robin submits (protocol). Authorship and the wording of the AI declaration are Robin's decisions. They are flagged here and not decided.**
Grades below come from the **registry** (`proofs/registry/hikita-star-*.json`), not from prose.

## Title options
1. *Positivity lost, two valuations survive: Hikita's ⋆-product at the ends of the s-line*
2. *Valuations of Hikita's (s,t)-deformed product of symmetric functions*
3. *An order filtration on Hikita's ⋆-product and the dominance support of e^⋆_λ*

## Abstract draft (~150 words, no \cite allowed)
Hikita introduced a two-parameter deformation ⋆ of the product of symmetric functions. Through Di Francesco–Kedem's quantum Q-system operators, ⋆-products of elementary functions are t-deformed graded characters. Di Francesco and Kedem observed that these lose Schur positivity at generic t. We show that two valuations survive, one at each end of the s-line. Write e^⋆_λ = Σ c_{λμ}(s,t) e_μ. At s=0, c_{λμ} is supported exactly on the dominance up-set of λ and has s-adic valuation n(μ). At s=1, the (1−s)-adic valuation of c_{λμ} is ℓ(λ)−κ(λ,μ), where κ is the maximal number of blocks into which λ and μ split compatibly. The first-order term is a symmetric biderivation (a carré du champ). The leading coefficients on coarsenings are sums over merge histories. Their weights are principal specializations, that is, q-dimensions. Both valuations come from a single degree count. We also give explicit e-basis Pieri rules. *(~150 words.)*

## Sections (12 pp budget)
Paper numbering differs from the session numbering because of a name collision: 217e "Thm A/B" (edges) vs Day 220 "Thm A/B" (s=1). Paper numbering: **T1** = Day 220 Thm 1, **T2** = Day 220 Thm A, **T3** = Day 220 Thm C, **T4** = Day 220 Thm B, **T5** = Day 220 Thm W.

| § | pp | Claims | Registry grade | Proof source (`proofs/`) | Cut to fit |
|---|---|---|---|---|---|
| 1 Intro: ⋆, DFK §8.3 | 1.5 | frame; "positivity lost" (cite DFK 1704 §8.3); dead positivity sweep (one sentence) | `star-structure-constants-no-positive-natural-basis` = dead-end (correct as a negative result) | `scripts/day217/positivity/SUMMARY.txt`, day219 sweep | drop the λ=1^n residual-positivity observation (computed only) or keep it to one line |
| 2 Subset formula + order filtration | 1.5 | (A_k)+(K_k) subset formula; Taylor piece p of E_k has order p | `subset-formula-Ak-Kk` proved; `thmA-premise-subset-formula-taylor` proved | 2026-09-26-day207b…PROVED.md; 2026-10-03-day220-s1-carre-du-champ.md §0 | state the Harrison-rigidity remark in one line, no proof |
| 3 s=0: DS | 2 | dominance up-set support, val_s = n(μ), integrality; s=t=0 → dominance zeta | `full-upset-support-exact-valuation` proved (hostile pass in-session) | 2026-09-30-day214-DS-all-lengths-PROVED.md | Lemmas 1.1/1.2 as a sketch; DS-from-(N) (n≤4 computed) becomes a footnote or is cut |
| 4 s=1: block law | 3 | T1 biderivation/carré du champ; T2 lower bound v ≥ ℓ−κ; T3 equality; T4 coarsenings exact + merge histories; T5 W_k(J) = (−1)^p[n]_t/[k]_t·p_J(1,t,…,t^{k−1}) | `s1-first-order-star-is-biderivation`, `s1-block-valuation-law-day220`, `thmC-block-law-exact`, `thmB-coarsening-exact-valuation`, `thmW-merge-weight-closed-form`: all proved | day220-s1-carre-du-champ.md §§1–5b | full proofs of T1 and T2 only; T3/T4 in sketch form; T5 statement plus one-paragraph proof; worked example (1,1,1)→(3) |
| 5 Pieri rules | 2 | e_k⋆e_r (207b), (TC), (★ℓ), presented as explicit evaluations of DFK 1704 Thm 5.17/(5.29) shuffle on 1 | `ek-star-er-pieri-all-k`, `two-column-gf-rule`, `ell-column-rule` proved | day207b, day209, day212 PROVED files | state the rules only; residue proofs → "full version"; (TC) is subsumed by (★ℓ), so keep only the ★ℓ statement |
| 6 Edges + remarks | 1 | t=0 edge = **DFK15 Cor 5.18 (theirs)**; t=∞ H′ = DFK 1908/15 (theirs); s=0 Thm H **demoted to a remark**; (N) as a remark (folklore-implicit, DFK 1704) | `theorem-B-t0-edge`, `theorem-H-prime…`, `theorem-H-s0-star-is-HL-pieri`, `square-theorem-all-edges-day217e` proved | day217e, day215, day216b/c | no proofs in this section; drop the s=∞ edge (217e Thm A) unless there is space |
| refs + AI decl. | 1 + ∞ | — | — | — | — |

## Novelty status per theorem
| Claim | Status | Locator / owner |
|---|---|---|
| Subset formula | ours as a formula; method cf. Hikita | Hikita 2503.23597 |
| DS (s=0) | **audited, clean** | `reading/2026-10-01-DS-novelty-audit.md`; method cf. Macdonald VI, KN integrality |
| T1–T5 at generic t | **UNAUDITED** (registry note) | `questions/q-block-law-novelty-inverse-HL-monomial.md`; generic t judged "less likely classical" |
| T3 at t=0 | **UNAUDITED, highest risk** | (1−q)-valuation of inverse HL P→m matrix; Macdonald III, Kirillov math/9912094, Wheeler–Zinn-Justin |
| Pieri 207b/(TC)/(★ℓ) | 207b audited, clean on the result, method partly known; ★ℓ vs Chen 2504.17508/Saito still unread | `reading/2026-09-26-novelty-ek-star-er.md`; referee risk: DFK 1704 shuffle + Thm DofM |
| t=0 edge | **NOT OURS**: DFK 1505.01657 Cor 5.18 | must be cited as theirs |
| H′ (t=∞) | NOT OURS in substance | DFK 1908.00806 Thm KNAN / DFK15 |
| (N) | folklore-implicit | DFK 1704.00154, BGHT; Cherednik book ch. 3 unread |
| Thm H (remark) | clean on Hikita; Chen–Lu–Ruan dispute open | `reading/2026-10-01-theorem-H-novelty-audit.md`, Browse 159/160 |
Also cite: Kirillov–Noumi (q-alg/9605004/5), Macdonald III (2.15), (3.2).

## Risks
1. **Novelty gate at t=0.** The t=0 block law may already be classical as a statement about the inverse HL P→m matrix. If it is found, cite it and present T3 as the generic-t deformation, carried by T2 and T4/T5.
2. **DFK15 normalization (PENDING).** The (N)-free claim needs a first-hand match of (5.15)/(5.25)/(5.27) with E_k|_{t=0}. Until then the registry premise `thmC-premise-217e-thmB-t0-edge` is proved *from (N)*, which imports the Cherednik C1–C3 locators that Clio flagged as UNVERIFIED. If the match fails, T3 falls back to "proved on coarsenings (T4) + computed elsewhere".
3. **Chen–Lu–Ruan 2601.13497 Cor 2.10** only matters if Thm H stays as a theorem. v7 demotes it to a remark, so the risk is LOW.
4. **"Computable from DFK"** is a fair referee objection to §5. Answer: the rules are explicit and nobody has computed them. Read Thm DofM before submitting.
5. **T5 has never been cold re-read** (single in-session derivation). It needs a sober re-read before it goes into the abstract.

## Timeline (2026-10-04 → 2026-11-15)
- **Oct 4–8:** run the t=0 novelty gate (risk 1). Do the first-hand DFK15 normalization match (risk 2). Cold re-read of T5. Email Robin this outline and ask about authorship and the AI declaration.
- **Oct 9–12:** scope freeze. LaTeX skeleton in `FPSAC2027.cls`; write §§2–4.
- **Oct 13–20:** §§1, 3, 5, 6; bibliography; first full draft.
- **Oct 21:** send to Clio for review; send the draft to Robin.
- **Oct 22–Nov 5:** revisions; page-cut pass; draft the AI declaration for Robin to edit.
- **Nov 6–10:** final polish; Robin approves.
- **by Nov 13:** Robin submits via SoftConf (keep a 2-day buffer).

## Grade discrepancies noticed (registry vs prose)
- The drunk summary in `day220-s1-carre-du-champ.md` says T5/Thm W is "computed only". §5b and the registry say **proved**, but the registry adds "NOT yet cold re-read". This outline uses *proved, single derivation*.
- The v7 prose says the package is "(N)-free". In the registry, T3's premise is still 217e Thm B, which is **from (N)**. The DFK15 swap is only a *proposal* in the registry notes.
- Check counts for T5 differ between sources: 68/68 (n≤7) vs 112/112 (n≤8).
- The registry note on `theorem-H-prime` still says "novelty UNAUDITED vs Hikita", but Browse 156 recorded it as clean. The note is stale. It does not matter here, because H′ is cited as DFK's.
