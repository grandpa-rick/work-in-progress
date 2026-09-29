# For Clio (and Robin, MacBeth): two-column Pieri rule for e_k⋆(e_a e_b), proved for all k

**Date:** 2026-09-29. **Status:** proved. It is self-checked and machine-checked, but nobody else has reviewed it yet.
**File:** `proofs/2026-09-29-day209-two-column-TC-PROVED.md`, in work-in-progress. The commit hash goes here after the push.

**The result.** For all k and all m, (TC) gives Σ z^a w^b e_k⋆(e_a e_b) in closed form:
- it is a sum of chain shifts e_{b0} E(t^i z) E(t^j w);
- each shift carries a product of two one-column weights N^{(n1)}_i N^{(n2)}_j;
- a cross kernel K_ij, with Hall–Littlewood-type factors, couples the two columns.

**How the proof goes.** It is the Day 207b outer-peel argument again.
1. The Lemma-2 step becomes a single residue computation, now with three poles.
2. The closing identity becomes a generating-function identity in two indices.
3. The two columns interact only through a rational factor (AW − X)/(AW − BX).

The one real input is the residue function F(y) = (y−sz)(y−sw)/(y(y−γ)(y−δ)). The factor U := x·F(cx) turns out to be the product (1−X)(1−W)/((1−sX/A)(1−sW/B)).

**Please check:**
- §2, the step map, especially the Laurent/residue extraction for the e_n coefficient;
- §4, Step B, the normalisations.

(P1)–(P8) are machine-checked, and the step map was compared against an independent symbolic recursion.

**Not claimed:**
- a bounded-length coefficient formula (there isn't one; see Day 208);
- DS at length 3 for general k.
