# Wake 257 Q41 — does the OG hyperplane identity separate (3) <-> (2,1)?   Grade: computed

**Verdict: OG hyperplane identity separates the swap: NO.**

Scripts: `q41.py M TRUEMAX` (outputs `out_m3..6.txt`), `isocheck.py` (`out_isocheck.txt`).

## Coefficients c_rho of [O_{X_kappa}]·O^1 on OG(m,2m+1), four independent ways
| m | #kappa | BR eq.(hsumB) == KOG tableaux | == draft Lemma og-corners | == TRUE K-class (G/B Demazure solve, BR-free) |
|---|---|---|---|---|
| 3 | 7  | 7/7  | 7/7  | 7/7 |
| 4 | 15 | 15/15 | 15/15 | 15/15 |
| 5 | 31 | 31/31 | 31/31 | 31/31 |
| 6 | 63 | 63/63 | 63/63 | not run |

Identity H(j)-H(j-1)=Σ c_rho H_rho(j) (H = Omega(P,j+1)) holds on the half-integer grid j∈[-8,12] for all kappa, m=3..6.

## (a) The expansions, m=3,4,5,6 (they are the same at every m)
- kappa=(3):   [O_X(3)]·O^1 = 1·[O_X(2)]   (theta={(0,2)})
- kappa=(2,1): [O_X(21)]·O^1 = 1·[O_X(2)]  (theta={(1,1)})
These are **literally the same expansion**: same rho=(2), same coefficient. The unique poset isomorphism P_(3)→P_(2,1)
(0 of 1 isos respect the diagonal marking) sends D(2) to D(2), so the transported expansion is also equal.
H_(3)=H_(2,1)=C(j+3,3) (Omega model = type-B Demazure chi at j=0..3: 1,4,10,20).

## (b) Swap test
Literal swap (H_(3) with (2,1)'s expansion and vice versa) and iso-relabelled swap: both hold identically, m=3..6.

## (c) All unmarked-isomorphic pairs, m<=6
The only pair is (3)~(2,1) (backtracking iso search, m=3..6; brute-force over all permutations confirms for |kappa|<=7).
H equal: yes. Expansions differ: no.
Structural reason (consistent with the computation): by Lemma og-corners the OG coefficients are (-1)^{|A|-1} over
nonempty antichains A of maximal elements, rho = P_kappa \ A — an unmarked-poset invariant, so no OG Pieri identity
at p=1 can see the marking.

## Negative controls (checker can fail)
- NC1: one sign flipped in (2,1)'s expansion → identity FAILS (all m).
- NC2: LG coefficients (BR Cor 5.4): (3) ↦ 2·(2) − (1), (2,1) ↦ (2). Transport comparator reports DIFFERENT
  (so the comparator detects separation where it exists — LG separates), and plugging the LG coefficients into the OG identity FAILS.
- NC3: a cyclic-shift bijection of the staircase is rejected by the iso test.
