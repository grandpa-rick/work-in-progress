# FPSAC 2027 cut of the LG/OG draft: WIP plan (2026-10-26)

**Status: plan only. Not for submission.** Submission, authorship, presenter and venue are Robin's
decision. His yes/no on "(b) arXiv now, (a) FPSAC cut" has been asked three times and is still
unanswered. Deadline: **2026-11-15**.
Skeleton: `/home/agent/projects/writing/lg-og-fpsac-cut/main.tex`. It currently compiles to 11 pp,
and the AI declaration starts on p. 11 (see "Skeleton" below).

## Format rule (checked first-hand)

- Class: `FPSAC2027.cls` (Pilaud 2026-09-03), copied from
  `/home/agent/projects/work-in-progress/fpsac2027/FPSAC2027.cls`. It is 12pt letter with 1in
  margins and `includehead`, uses mathpazo with linespread 1.05, and loads biblatex with
  backend=bibtex.
- Page rule, from curl of `maths.universityofgalway.ie/fpsac2027/submissions/` on 2026-10-26:
  submissions must be **6–12 pages, including the bibliography, excluding the AI declaration**.
  - The AI declaration is mandatory and has no length limit. It sits at the end, in the
    designated section.
  - No changes to the preamble, font, margins or spacing. Packages like savetrees are banned.
  - pdflatex, then bibtex.
- Aside: the old Hikita FPSAC draft (`work-in-progress/fpsac2027/fpsac2027-draft.pdf`, 14 pp) has
  references running onto p. 13. Under the "including bibliography" rule it is about a page over.
  Earlier notes said "12 pp + uncapped AI declaration" and never addressed references.

## How the page estimates were made

| Draft | Class | pdftotext chars/page | TeX bytes/page |
|---|---|---|---|
| this draft | amsart 11pt, 1.1in margins | ≈3,500 | ≈4.4 KB |
| old Hikita draft | FPSAC2027.cls | ≈2,450 | ≈3.1 KB (math-heavy) |

- So **1 draft page ≈ 1.43 FPSAC pages**.
- The draft body (111.6 KB of TeX, 25 pp before the references) would run to about 36 FPSAC pp.
- The 12 pp cap leaves about 10.2 pp of body after title/abstract (0.6 pp) and references
  (about 1.2 pp for 18 refs). The cut must therefore keep **about 28% of the text**.

## The story (one sentence)

> On the Lagrangian Grassmannian, Stanley's order-polynomial reciprocity survives with one change:
> the lattice is $\tfrac12\mathbb Z$ off the diagonal. So $(-1)^{|\kappa|}H_\kappa(-j)$ counts
> strict maps of the shifted diagram into $[2j-1]$ that are even on the diagonal.

The factor 2 is the bond of short-root LS steps. The parity rule is Buch–Ravikumar's "no primed
diagonal". Everything else (integer window, half-integer law) is read off this count or its
recursion.

## Dependency chain (from the draft, numbering as printed)

```
(G1)+(LS) ── L3.1 bonds ─ L3.2 LS-parity ─ P3.3 model ─ Ehrhart–Macdonald ─► Thm 1.1 (SR_LG)   [proof 1, no BR]
(G0)+(G1)+(BR Cor5.4/Thm5.5) ─ L2.4 dictionary ─ L4.2 reflection ─ P4.6 (BR p=1) ─ Thm 4.8 (Q30)
      └─ L4.1 hyperplane ─ P4.12 (★) ─ L4.13 peel ─► Thm 1.1(b) (proof 2);  BR L4.5 ─ L4.14 ─► Thm 1.3 (SR_OG)
Thm 1.1(b)/1.3 + L2.3 ─► Thm 1.4 (window, from counts)        [Prop-4.6-free via proof 1]
(★)+L2.6 anchor + L5.1 chain drop (support only) + L5.2 ─► Thm 1.4 (inductive; Clio uid457)
L2.5 polynomial ext. ─ L6.1/6.2 ─ P6.3 bridge ─ Thm 6.5 (H(-1/2)=2^{-κ1}), Thm 6.6 (-3/2)
Thm 6.6 + L5.1 (support only) ─► Thm 6.7 (Z: vanishing half; Clio uid451)
(★)=A_oA_e (L6.8) + Thm 6.5 ─ L6.9 ─ C6.10 row formula ─ L6.11 shift ─ L6.12 one row ─► Thm 1.5 (HZ) ─► Cor 1.6 ("Thm B")
L6.8 + L6.17 Möbius + L6.18 X,P ─► Prop 6.16 (half-step T)
```

Key structural fact: the **headline Theorem 1.1 and the window Theorem 1.4 can be presented
without Prop 4.6**, by the LS route plus the counting proof. Prop 4.6 is needed for:
- Theorem 4.8 and the second proof of Theorem 1.1;
- all half-integer results except the vanishing half, which Theorem Z covers using support only;
- Prop 6.16.

## Outline with page budget (FPSAC pages; total ≤ 12 including references)

| § | Content | Draft source (TeX KB → FPSAC pp if uncut) | Budget |
|---|---|---|---|
| – | Title, abstract, keywords | – | 0.6 |
| 1 | Intro: motivation (quadric vs P³), Thm 1.1, Ex 1.2, Thm 1.3 (classical, one line), Thm 1.4, Thm 1.5, Cor 1.6; "why the factor two" (10 lines); literature (8 lines) | §1: 12.8 KB → 4.1 | 2.0 |
| 2 | Shifted diagrams, counts $N$, L2.3 (full), reciprocity, Schubert varieties + cited inputs (G0, G1, LS, BR); dictionary L2.4 in one sentence; L2.6 anchor | §2: 16.3 KB → 5.3 | 1.2 |
| 3 | LS proof: L3.1 (sketch), L3.2 (full), P3.3 (sketch), proof of Thm 1.1 (full) | §3: 10.4 KB → 3.4 | 1.8 |
| 4 | L4.1 (sketch), BR quote (8 lines), P4.6 (sketch + case table), Thm 4.8, Ex 4.10, Rem 4.9, P4.12, L4.13, second proof + OG proof | §4: 26.9 KB → 8.7 | 1.8 |
| 5 | Window: proof from counts (full, 3 lines) + remark on inductive proof | §5: 6.0 KB → 1.9 | 0.5 |
| 6 | P6.3 bridge (sketch), Thm 6.5 (full), L6.8, C6.10, HZ proof sketch, Cor 1.6; remarks: Thm Z, Prop 6.16 (statement only) | §6: 28.8 KB → 9.3 | 1.7 |
| 7 | Table 2 (values on LG(3,6)), two readings, checks in one sentence; 3 open questions | §7–8: 7.9 KB → 2.6 | 0.8 |
| – | References (~18 of the 27) | – | 1.2 |
| | **Total** | about 36 if uncut | **11.6** |

That leaves about 0.4 pp slack. Cut first if over, in this order:
1. Example 3.4.
2. The Prop 6.16 remark.
3. The second and third open questions.
4. Cite the Theorem Z remark in one line instead.

## Full proofs / sketches / cited to the full version

- **Full proof (compressed):**
  - Lemma 2.3, Lemma 2.6 (anchor), Lemma 3.2.
  - Proof of Thm 1.1 from Prop 3.3 plus Ehrhart–Macdonald.
  - Lemma 4.13 (peel), Thm 4.8 (from P4.6), proof of Thm 1.4 from the counts.
  - Thm 6.5, Lemma 6.8, Cor 6.10, Cor 1.6.
- **Sketch:**
  - Lemma 3.1 (bonds), Prop 3.3 (bijection), Lemma 4.1 (hyperplane).
  - **Prop 4.6** (induction over BR Thm 5.5: state the p=1 table, one line per case).
  - Second proof of Thm 1.1 and K-proof of Thm 1.3.
  - Prop 6.3 (bridge), Thm 1.5 (fibre + shift + one-row).
- **Statement or remark only, cited to the full version:**
  - Lemma 2.4 (dictionary), Lemma 2.5 (polynomial extension), Lemmas 4.2–4.5 (reflection,
    two-letter counts).
  - Lemma 4.14 (OG corners), Lemmas 5.1–5.2 and the inductive window proof.
  - Lemmas 6.1–6.2, Thm 6.6, Thm 6.7 (Z), Lemmas 6.9, 6.11, 6.12.
  - Remark 6.15 (cohomological formula), **Prop 6.16 + Lemmas 6.17–6.18**, Table 1.
- **Dropped:** §1.3 Method, §1.5 Outline, Remarks 3.5–3.6, Ex 4.11 (dictionary test),
  Remark 4.15 (OG classical; becomes one sentence), §7.1, most of §7.3, Questions 8.3–8.4.

**Examples and figures kept:**
- Ex 1.2 (LG(2,4) quadric);
- Ex 4.10 (κ=(3,2) Pieri expansion, reused for the window: $s_3=1$, $s_4=13$);
- Table 2 (values at negative integers and half-integers on LG(3,6));
- the example κ=(3), a=2 after Cor 6.10.

There are no figures in the draft. If 0.3 pp is freed, add a TikZ shifted diagram of (3,2) with
diagonal boxes marked and one strict even-diagonal filling. A combinatorialist reads this faster
than the definition.

## Risks

1. **Outside verification is thin where the cut leans hardest.** Clio has line-read and endorsed
   only Thm Z (uid451) and the window Thm 1.4, both types (uid457). Nobody else has line-read
   any of the following:
   - **Thm 1.1 itself**, by either proof;
   - **Prop 4.6** (refereed only by Rick against BR Thm 5.5, plus Clio's numerical check on
     1243 shapes);
   - **Thm 1.5 (HZ)**;
   - **Cor 1.6 (= "Thm B")**;
   - **Prop 6.16**, which was proved 2026-10-25, the day before this plan.

   Mitigation:
   - Lead with the LS proof, which does not need Prop 4.6.
   - Prove the window from the counts.
   - Back the vanishing half of HZ with the endorsed Thm Z.
   - Prop 6.16 is a remark only.
   - The nonvanishing and sign parts of HZ and Cor 1.6 still hinge on Prop 4.6. The AI
     declaration must say so per result.
2. **Novelty and locator exposure on the headline.**
   - Dehy 1998 (via Mathieu §4) interprets Demazure dimensions as Ehrhart polynomials. It has
     not been read first-hand, and it may anticipate Thm 1.1(a)/(c).
   - The CCFL lecture-hall work (2609.03195) is the nearest live competitor.
   - The Littelmann Invent. 116 §5.2 locator is quoted second-hand; only the ICM survey was read.
   - In 12 pages the hedges must survive compression. Keep the Dehy footnote in the literature
     paragraph.
3. **Process and page-limit risk.**
   - Robin's venue decision is still pending. FPSAC needs him as submitter and presenter, with
     about 3 weeks left.
   - The cut must stay at or below 12 pp *including references*. The skeleton projects 11.6–11.8 pp
     once the TODO paragraphs are written at their budgeted lengths, so slack is under half a
     page.
   - If the Hikita abstract is also submitted, it is currently over the limit (see the aside
     above).

## Skeleton

`/home/agent/projects/writing/lg-og-fpsac-cut/main.tex`, with `FPSAC2027.cls` and a trimmed
`refs.bib` (20 entries) beside it.
- Statements of kept results are copied verbatim with the same labels and macros.
- Proof placeholders are red `\budget{n}{...}` boxes that reserve the estimated number of lines.
- Build: pdflatex, bibtex, pdflatex ×2. **11 pages**: body plus references end on p. 10, and
  the AI declaration starts on p. 11.
- Prose TODOs are written shorter than their budgets, about 57 lines (≈1.4 pp) short. Only 11 of
  the planned ~18 references are cited yet (≈0.4 pp). The projected final length is about
  11.8 pp, inside the cap.
