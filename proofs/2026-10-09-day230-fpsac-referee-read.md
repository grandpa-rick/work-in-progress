# Day 230 PROVE — hostile-referee cold read of the FPSAC 2027 draft

Date: 2026-10-08. Target: `work-in-progress/fpsac2027/fpsac2027-draft.tex`, starting from c15ff37.
No new math. Registry untouched: no grade errors found.

## 0. Page limit: FIXED
Before: body+refs = 13 pp (about 6 lines of refs spilled onto p.13).
After: body+refs = **12 pp**, and the AI disclosure starts on p.13 (disclosure pages are uncounted).
Cuts:
- Penrose bib entry: the visible `note = {TODO verify chapter title…}` went to `xnote`, which biblatex ignores.
  The flag is still in the .bib source, so the chapter title is still unverified (see §5).
- "Organisation" paragraph: cut down to one sentence.
- Duplicate definition of T_a in §6 removed; T_k is now defined once, in Lemma 3.9.
- Open problem 5: shortened, which also removes an overclaim (item B6 below).

## 1. Issues a referee would hit, all fixed
A. Internal jargon leaked into the PDF
1. Remark 3.5 title read "(N)-freeness". "(N)" is our internal name for the ∇-transport theorem and is not defined in the paper. Removed.
2. Proof of Thm 3.7 ended "This proof does not use (N)." Removed.
3. Remark 4.3 contained "**Must cite, not claim.**", a note-to-self. Removed.
4. The "Computer checks" paragraph said "all 16 class-4 pairs". "Class 4" is a Day-224 internal classification. Now reads "16 pairs with n≤12 (both examples among them)".

B. Undefined objects, notation clashes, overclaims
1. Intro: "Λ_{q,t}" was never defined. Now "symmetric functions over Q(q,t)".
2. Prop 2.1 said "N ≥ deg" without saying deg of what. Now "for every N", which is the form Thm 5.1's proof uses.
3. Thm 3.7: "tight histories" was used but never defined. The definition now says that a merge goes into one block of size k+|J| via E_k^{(p)}, and that the history ends at μ when the final block sizes are μ.
4. Lemma 3.9 used T_k in its statement, but T_k was only defined inside the proof of Thm 3.8. The definition is now in the lemma.
5. Remark 4.3: "Tutte polynomial of the graph M_λ" left M_λ undefined. Replaced by the explicit statement K_λ=(t−1)^{ℓ−1}T(1,t) for the multigraph with λ_iλ_j edges between i and j, which I checked by hand.
6. Open problem 5 said "Dołęga's column cumulants **are** the t=0 specialisation of the ⋆-leads". Only one gauge invariant, J at n=4, supports this, so it overclaimed. That sentence is cut; the J comparison stays in Remark 4.3, where it is stated precisely.
7. Open problem 3 wrote "Lead_{λ,(n−1,1)}(1)=2n²−6n+3" without saying which λ. From the Day 224 data it is every λ with ℓ(λ)=3 and κ=1 (checked against the log: n=6, 7, 8 give 39, 59, 83). Now stated that way.
8. Symbol clashes:
   - κ(λ) (the Haglund–Tewari cumulant) against κ(λ,μ): renamed C_λ.
   - δ_j (the degree count in the proof of Thm 2.2) against the derivations D_k: renamed δ_j.
   - The CT kernel K(z) against K_λ and the Kostka K_{νλ}: renamed Ω(z), and its local factor κ is now θ.
   - Λ_a(r,q) (a number) against the ring Λ_a: renamed Ξ_a(r,q).
   - In Example 6.5, m=ρ₁−ρ₂ against m_{xy}: renamed r.
9. Example 6.5 used m_{xy} before it was defined (the definition was in Thm 6.6, which comes later). It is now defined in the example, and Thm 6.6 refers back to it.
10. Thm 2.2: "at s=t=0, e⋆_λ is the zeta function" was imprecise, since b_μ=0 at s=0. Now it says the s-leading coefficients are all 1 at t=0. Checked at λ=(1,1,1).
11. Proof of Thm 3.2: a_ρ(s) was undefined. Added f=Σa_ρ e⋆_ρ.
12. Open problem 4 referred to "(Example above)". Now a \cref to Example 6.8.
13. Overfull hbox (7.5pt) in the proof of Prop 2.1: fixed by rewording.

## 2. Grades vs registry (checklist item 2): PASS
Every `% registry:` node cited in the tex has trust `proved` with a `file` (script dump in session).
Two claims rest on computation, and the text says so: Dołęga's J≡3/2 ("by computation") and open problem 3 ("computed only").

## 3. Novelty sentences (item 3): no change needed
- Two-row Green: the text says it is classical, via III (7.6′) plus the two-part-class range treated by Morris, and that Morris was not seen first-hand (footnote).
- Lemma 3.9: "probably known".
- Graph half: credited to Dołęga, Gessel–Sagan, Mallows–Riordan and Penrose.
- The t=0 block law: "classical in substance".
- The d-matrix: credited to Kirillov.

## 4. Hand / independent checks (item 4)
- `scripts/day230/referee_v2.py` implements **Thm 6.6 + Thm 6.3 + Prop 3.3 + Lemma 3.9 literally from the printed text**, with fresh code that does not reuse the old closed.py. On all 123 (ordering, pair) cases in the subset-formula engine logs (day224 newcases_n10 + table3_n8) it gives **123/123 matches, 0 mismatches**. Both printed examples, (3,3,3)→(7,2) and (4,4,2)→(7,3), are reproduced exactly.
- `referee_sym.py`: the printed Thm 6.3 is symmetric in x↔y on 4 cases. So no y≤x hypothesis is needed, and none is stated.
- `referee_J.py`: J from Thm 4.1 equals the printed rational function. Lead(1,1,1)→(3) = (t+2)(t²+t+1).
- `referee_subset.py`: the printed Prop 2.1, run directly with N=3 to get e⋆_{111}, gives:
  - c_{111,111}=s³=s^{n(111)};
  - c_{111,21} has v=1 with lead −3(1+t), which matches Prop 3.3, and its s-valuation is 1=n(21);
  - c_{111,3} has v=2 with lead t³+3t²+3t+2, which matches Thm 4.1, and its s-valuation is 0;
  - the s-leading coefficients at t=0 are all 1, which matches Thm 2.2.
  These are consistent with Thms 2.2, 3.4 and 4.1.
- Two-row Example 6.5 by hand at λ=(3,1) and (2,2), classes (3,1) and (2,2):
  - from the case formula: t, t−1, t²−1, t²−t+2;
  - from π_k+(t−1)Σt^{k−1−j}π_j: the same four values. ✓
- By hand:
  - Cor 4.2(a)–(c);
  - Thm 3.8 at t=0;
  - Box exponent: consistent with Thm 2.2 (c_{λλ}=s^{n(λ)} maps correctly), and two complements give the Column Lemma's s^{C(ℓ,2)};
  - Cor 5.5 at k=1 reproduces Hikita's Thm 3.12;
  - the ∂_sπ_a(j) step in Prop 3.3.

## 5. Left open (not fixable in this session)
- Robin's 7 \todo items: authorship, affiliation, acknowledgements, disclosure wording, model list, human role, and the footnote. Untouched.
- Penrose 1967 chapter title and pages: still unverified. The flag is now an `xnote` in the .bib source, not visible in the PDF.
- Hikita locators (Lemma 3.1/3.2/3.3, Cor 3.9/3.10, Def 3.4, Prop 3.8, Thm 3.12, Lemma 6.3) and DFK locators: browsing is not allowed this session, so they were not re-verified against the papers. The R0 locator was fixed in wake 230.
- Source-only `%` comments still contain internal provenance (Day numbers, node ids). They are invisible in the PDF; strip them before any arXiv upload of the source.
- Cor 4.2(d) excludes t=1, which is unnecessary (Thm 3.8 already covers t>0). Harmless, so I left it.
