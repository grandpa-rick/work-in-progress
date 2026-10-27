# FILL-PLAN: FPSAC 2027 cut, BR-free route (2026-10-27)

WIP plan only. Not committed. Supersedes the "Outline" and "Dependency" parts of `FPSAC-CUT-PLAN.md`,
which still use **pre-Prop-6.16 numbering** (e.g. "Thm 6.5" there is `thm:A` = 6.6 now). This plan uses
labels plus the numbers printed in the current draft.

**Paths.**
- The cut lives in two identical copies:
  - `work-in-progress/papers/lg-og-fpsac-cut/` is git-tracked (cf8868e, 5a2fa86, d118e27). This file is there.
  - `writing/lg-og-fpsac-cut/` is an untracked twin.
  - `projects/papers/lg-og-fpsac-cut/` does not exist.
- Source: `writing/lg-og-negative-degree/main.tex`, 1950 lines, 28 pp. It is byte-identical to
  `work-in-progress/papers/lg-og-negative-degree/main.tex` (HEAD 53675c2).

## 0. Compile status (scratch copy /tmp/cutbuild, latexmk -pdf)

- The build is clean: exit 0, 11 pp. Class `FPSAC2027.cls` is present. It is the Pilaud 2026-09-03 version and
  identical to `work-in-progress/fpsac2027/`. Biblatex runs on the bibtex backend.
- Page layout:
  - The body ends on p. 9.
  - References run from p. 9 to p. 10.
  - The AI declaration starts on p. 11 after `\clearpage`.
  - So 10 pp count toward the limit.
- Warnings:
  - 2 overfull hboxes, at cut lines 375 (18.4pt) and 445 (1.8pt). Both are inside TODO text.
  - 1 overfull vbox (3.3pt).
  - 9 hyperref token warnings.
  - No errors and no undefined references.
- There are 45 visible TODOs. The abstract contains no `\cite`, which is compliant.
- Of the 20 refs.bib entries, 9 are cited so far.

## 1. Section table (target ≤ 12 pp incl. refs; AI declaration uncounted)

| Cut § | Supplied by (main.tex label = printed no., line) | What must change for the BR-free route | Est. pp |
|---|---|---|---|
| Title + abstract | — | **Abstract**: "We give two proofs, via LS paths and via a closed formula for the K-theoretic Pieri coefficients" → one proof (LS + Ehrhart–Macdonald). Add "a K-theoretic reading of the parity rule recovers the Pieri coefficients of the divisor" as a remark-level clause, or drop it. Keep BR by name only, with no \cite (FPSAC rule). | 0.6 |
| 1 Intro | thm:main 1.1 (l.113), ex:intro21 1.2 (l.133), thm:og 1.3 (l.146), thm:window 1.4 (l.158), thm:HZ 1.5 (l.173), cor:OC 1.6 (l.186); ssec:why (l.195–216); ssec:lit (l.240–259) | **Do NOT copy draft l.226–237** (Method). It says the second proof's "recursion extends further … proves Thm HZ". Since PROVE 258 that is stale for the cut: HZ comes from the LS route. **Literature novelty sentence (l.257–259)**: drop "Theorem Q30 as a closed formula" from the novelty claims (see Risk 2). Keep the Dehy footnote hedge. | 2.0 |
| 2 Prelim | lem:ideals 2.1 (l.285), lem:chain 2.2 (l.301), lem:Npos 2.3 (l.329, full), ssec:reciprocity (l.348), ssec:schubert (l.361), ssec:inputs (l.462–500), lem:P 2.5 (l.502), lem:anchor 2.6 (l.519) | **Cited-inputs list.** Main line: (G1) and (LS) only. Move (G0), (BR) and the dictionary lem:dictionary 2.4 (l.398) out of §2 into the §4 remark; they are needed only for Thm 4.8. **Add lem:P 2.5 as a one-sentence statement**, because the BR-free (★) needs polynomiality (H_κ polynomial; also Snapper/Kleiman, Kleiman66 in draft refs). **Add lem:R 4.3** (l.828, 3-line statement + 3-line proof). Pure combinatorics; needed for the explicit form of V(θ), the chain-drop remark and the Thm Z remark. | 1.3 |
| 3 LS proof | lem:bonds 3.1 (l.557), lem:lsparity 3.2 (l.597), prop:model 3.3 (l.615), ssec:proof-main (l.649–665) | No change in content. This proof is already (G1)+(LS) only. Define F_κ in §2 (open TODO). Locator hedge for Lit94 §5.2 (l.494 footnote): keep it, compressed. Cut ex:ls-small first if over. | 1.8 |
| 4 Recursion (★) [**retitle**: "The recursion and the Pieri coefficients"] | prop:star 4.12 (l.1052), lem:peel 4.13 (l.1069), **BR-free proof of prop:star** (l.1101–1111), the non-recovery paragraph (l.1113–1120); thm:Q30 4.8 (l.997), rem:primes 4.9 (l.1013), ex:pieri32 4.10 (l.1021) | **Replace** the cut's current proof skeleton: delete the BR quote TODO, the `prop:BRp1` sketch (9-line budget) and the "Second proof" (13 lines). **Insert**: (i) lem:peel with its full proof (counting, 6 lines). (ii) prop:star **LS-only**: thm:main(b) gives s_j=N_κ(j) for all j≥1 and all κ (κ=∅ included). The peel is then (★) at every integer, both sides are polynomial (lem:P), so (★) holds identically (4 lines). Source: proofs/2026-10-27-star-without-BR.md; registry `q488-oc-star-BR-free`. (iii) **Drop the OG clause** of prop:star. It uses BR Lemma 4.5 via lem:og-corners 4.14, and nothing in the cut needs OG (★). (iv) **Remark (Pieri coefficients)**, ~8 lines: lem:hyperplane 4.1 statement in one line. Thm 4.8 c_ρ=(−1)^{|θ|−1}V(θ) is stated **with \cite{BR12}** [Cor 5.4 + Thm 5.5 at p=1 via the dictionary; p=1 case of BR's Cor 5.6]. Then ex:pieri32, rem:primes in 2 lines, and the warning that (★) alone does not determine c_ρ: {H_ρ} is dependent at κ=(3,1). Last, a one-line pointer: "Thm 4.8 also follows without (BR) by working T-equivariantly: the T-characters are independent though the H_ρ are not (full version)." Source: proofs/2026-10-27-weighted-peel.md Thm E. | 1.5 |
| 5 Window | thm:window proof from counts (l.1177–1182); lem:chaindrop 5.1 (l.1190), lem:A 5.2 (l.1207), ex at l.1273 | The proof from counts needs thm:main(b) on LG and thm:og on OG; keep it (3 lines). **Remark (inductive proof)**: on LG phrase it with V(θ) and lem:R, as the draft now does (l.1197–1201). **Strike the parenthetical "(equivalently c_{ρ_{t*}}≠0, by Thm Q30)" from lem:A(ii)** if lem:A is quoted. The OG inductive proof uses BR Lemma 4.5: say "LG" only, or cite. | 0.5 |
| 6 Half-integers | rem:deps 6.1 (l.1297), lem:W 6.2, lem:C 6.3, prop:bridge 6.4 (l.1347), thm:A 6.6 (l.1370), thm:C 6.7 (l.1385), thm:Z 6.8 (l.1409), lem:layers 6.9 (l.1450), lem:eps 6.10 (l.1466), cor:row 6.11 (l.1481), eq:fibre (l.1505), lem:shift 6.12 (l.1517), lem:onerow 6.13 (l.1536), prop:T 6.17 (l.1622) | **Fix the stale header comment** in the cut ("RISK: depends on (star) hence on Prop 4.6"). Replace it with a 2-line rem:deps: §6 uses (★) from §4 (LS-only) plus Demazure computations, so no (BR) and no (G0). **Thm Z remark**: use the V(θ)-support version (draft l.1425–1440 proof). Do NOT copy l.1434–1436 ("runs with eq:hyperplane … rim condition of BR Thm 5.5") except as an optional clause with \cite. Content is otherwise unchanged. The HZ sketch is 10 lines; prop:T is a statement-only remark. | 1.7 |
| 7 Example + questions | tab:lg3vals Table 2 (l.1866), q:cominuscule 8.1 (l.1895), Canonical modules (l.1904), Higher Pieri (l.1929) | The "higher Pieri c^ν_{λ,p}" question references BR; keep \cite{BR12}. Checks sentence: Demazure m≤5, Clio census m≤6. Do **not** cite rem:degrees (`q488-degree-identity` = computed). | 0.8 |
| References | — | ≈20 entries (below). | 1.2 |
| **Total** | | | **11.4** |

Slack is about 0.6 pp. Cut in this order if over: ex:ls-small, prop:T remark, 2nd/3rd questions, the Thm 4.8 example
(keep the statement).

## 2. Every BR / sub-"proved" dependency touched by the cut

**A. Current skeleton text that rests on BR. Each item must change or carry \cite{BR12} explicitly.**

1. **Abstract.** It claims "two proofs … via a closed formula for the K-theoretic Pieri coefficients". The second proof
   rests on prop:BRp1. That is the p=1 case of BR Cor 5.6, which BR state but do not prove. Rick refereed it against
   BR Thm 5.5.
2. **§1 literature TODO.** The novelty list includes Thm Q30.
3. **§2 cited inputs.** It lists (G0) and (BR), and the dictionary is needed only for Q30.
4. **§4.**
   - The BR quote TODO with Cor 5.4 and Thm 5.5.
   - prop:BRp1.
   - thm:Q30.
   - ex:pieri32.
   - The prop:star OG clause, which uses BR L4.5.
   - The "Second proof … and proof of Theorem og" TODO, which uses BR Lemma 4.5.
5. **§5 window remark.** On OG it needs BR L4.5 for chain drop and c_{ρ1}.
6. **§6.**
   - The header comment says it rests on Prop 4.6. That is stale.
   - The Thm Z remark says "support of the Pieri rule". That is BR unless phrased with V(θ).
7. **AI-declaration TODO.** It lists prop:BRp1 among the items not line-read.
8. **thm:og (OG).** Its BR-free proof is only "classical" (rem:og-classical; LS86, Seshadri, GL96, Hibi87 and
   Proctor84 were **not read first-hand**). The draft's own K-proof uses BR L4.5, and the OG LS route is
   **computed only** (rem:og-ls, 72 cases; type-B bonds lemma not written). The OG half of thm:window inherits this.
   **Decision needed:** state thm:og as classical with citations and say so, or keep the 5-line K-proof with
   \cite[Lemma 4.5]{BR12}.

**B. Registry nodes the planned cut uses (`proofs/registry/q488-converse.json`).**

| Node | Used for | Grade | Note |
|---|---|---|---|
| q488-oc-ls-bijection | lem:bonds, lem:lsparity, prop:model | proved | |
| q488-oc-SR_LG-via-LS-paths | thm:main | proved | Lit94 §5.2 locator second-hand; also lists a "sources: 2609.03195" (CCFL) |
| q488-oc-ls-parity-model | prop:model | proved | |
| q488-oc-even-reciprocity | thm:main(c)/(b) reading | proved | |
| q488-oc-integer-strict-reciprocity | thm:main(b) | proved | **approach text is only "127/127 at m=7"** (computed wording on a proved node; proof is §3) |
| q488-oc-anchor-H-minus-1 | lem:anchor | proved | |
| q488-oc-star-BR-free | prop:star (LS) | proved | children listed are Q45 nodes (not needed by it; registry structure oddity) |
| q488-oc-Zint-by-induction / q488-oc-strip-lemmaA / q488-oc-Ind-chain-drop | §5 remark | proved | approach text mentions BR (pre-258 wording); LG version rerouted via lem:R |
| q488-oc-prop-H-bridge, q488-oc-thmA-half, q488-oc-thmC-M3 | §6 | proved | |
| q488-oc-HZ-half-integer-law | thm:HZ | proved | **approach text still says "I+C=A_oA_e … from Q30", sources 1005.2605** (stale: pre-PROVE-258) |
| q488-oc-thmB-OC | cor:OC | proved | |
| q488-oc-Z-by-induction | thm:Z remark | proved | approach uses BR; BR-free only in the V(θ) reroute |
| q488-oc-half-step-T-ribbons (+T-mobius, T-XP-identities) | prop:T remark | proved | child **q488-oc-T-rim-guess = computed** (discovery only, not used in the proof) |
| q488-oc-Q30-KLG-strict-fillings, q488-oc-BR-dictionary, q488-oc-support-rule-BR | Thm 4.8 remark (with \cite{BR12}) | proved | BR-dependent by design |
| q488-oc-OG-corner-inclusion-exclusion, q488-oc-OG-strict-count, q488-oc-OG-vanishB-BR, q488-oc-ZOG-by-induction, q488-oc-OG-chain-drop | thm:og, OG window | proved | all BR L4.5 |
| q488-Q45-Q30-BR-free (+ Q45-equivariant-peel-lift, Q45-strict-character) | one-line pointer only | proved | **weighted-peel.md §7 says "proved modulo (KT)"** (Thomason locator not written); owner risk Lenart–Shimozono 1203.3237 not checked. Registry grade is stronger than the file. Keep it a pointer only. |

**Below "proved". None may be used as support in the cut.**
- q488-degree-identity (computed): rem:degrees. Excluded.
- q488-oc-NZ-computed-m7 (computed): evidence only. At most it can appear in the checks sentence.
- q488-oc-caseB-boundary-telescoping (computed), q488-oc-strip-lower-bound (computed), q488-oc-Z-perk (computed),
  q488-oc-OG-height-one (computed), q488-oc-NZ-beyond-reduction (computed). Not used.
- q488-oc-T-rim-guess (computed). Not used.
- q488-oc-NZ-edge-sign-coherence and q488-oc-root-multiset-of-top-cell (dead-end). Not used.
- The OG LS twin is computed, recorded inside the SR_LG-via-LS approach text. Not used.

## 3. Bibliography needed (from draft refs.bib; cut refs.bib already has 20)

- **Keep:**
  - BR12, for Thm 4.8, rem:primes and the higher-Pieri question.
  - Lit94, LitICM, for (LS).
  - Dem74, BK05, for (G1).
  - Mac71, BeckRobins, for Ehrhart–Macdonald.
  - Stanley70, EC1, for Stanley reciprocity.
  - Sta86, for the order polytope.
  - DeC79, Chi00, CCFL, Stem97, Dehy98, Mat95, for the literature.
- **Keep if thm:og is presented as classical:** Hibi87, GL96, LS86, Proctor84, and **add Seshadri** (draft has it; the
  cut's refs.bib lacks it).
- **Add:** Kleiman66, for polynomiality of χ, if lem:P is replaced by the Snapper/Kleiman citation.
- **Optional:** Lenart–Shimozono 1203.3237 (not in either refs.bib). Add it only if the T-equivariant pointer names
  K-Chevalley ownership. It has not been read, so it gets no locator.
- **Drop:** Jos85, And85, Ram85, Kum87, Mat88 (the (G1) attribution list; BK05 suffices), Pra91 (rem:degrees excluded).
- **Count:** about 20–22 entries, about 1.2 pp.

## 4. Top risks

1. **OG is the one result that is still not BR-free and not first-hand.** thm:og and the OG half of thm:window rest on
   BR L4.5 or on classical refs nobody read; the LS route is computed only. Decide: cite as classical, or keep the
   BR K-proof.
2. **Thm 4.8 novelty and grade.**
   - prop:BRp1 is the p=1 case of BR's own (unproved) Cor 5.6.
   - Equivariant K-Chevalley (Lenart–Shimozono, Lenart–Postnikov, Pittie–Ram) owns the values in principle.
   - Claiming Q30 as new in the intro is exposed. Present it as "a closed form, compare BR Cor 5.6", with \cite.
   - The T-equivariant BR-free proof is "modulo (KT)" in its own file, though the registry says proved.
3. **Stale dependency text can leak into the cut.**
   - Draft §1.3 (l.234) says HZ comes from the K-theoretic recursion.
   - The registry HZ node says "from Q30".
   - The cut's §4 and §6 comments say "rests on Prop 4.6".
   - The AI declaration must list per result: LS-only (thm:main, (★), §§5–6 on LG), BR-cited (Thm 4.8) and
     classical-cited (thm:og).
   - Separately, the page slack is about 0.6 pp, and Robin's venue decision is still pending, with 19 days left.

## Decisions (Wake 260, 2026-10-27)

- **OG:** the cut keeps the K-theory proof for the OG statements and cites Buch–Ravikumar Lemma 4.5 (and the classical
  refs) explicitly. No claim of a BR-free OG proof. The LG main line is BR-free via (G1)+(LS)+peel+polynomiality
  (PROVE 258, `proofs/2026-10-27-star-without-BR.md`).
- **Thm 4.8 → Remark.** Text: "the values c_ρ are the p=1 case of [BR12, Cor 5.6]; equivariant K-Chevalley
  (Lenart–Shimozono) determines them in principle; a BR-free T-equivariant proof will appear in the long version."
  No novelty claim.
- **AI declaration** sorts results into three lists: proved via LS only / cited from BR / cited classical.
- **Purge stale text:** "two proofs", "from Q30", "rests on Prop 4.6".
- **Budget:** ≤ 12 pp including the bibliography; slack 0.6 pp. Nothing is submitted without Robin's explicit yes.
