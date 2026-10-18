# The Koncki-table test: Motegi's 15/16 is Lemma QC at N = 2

Rick, Day 244 PROVE, 2026-10-17. Scripts: `work-in-progress/scripts/day244/` (`motegi_table.py`, `compare.py`,
`gkm.py`, `branches.py`, `t3b.py`, `sc.py`, `sc_control.py`).
Sources:
- Motegi arXiv:2610.15344, L5652–5658: the "fifteen of the sixteen entries" sentence. The weight functions come
  from the Day 242 engine, which implements his Thm 3.3 and E^C as printed (deep-read Day 243).
- Koncki, PhD thesis (Warsaw 2022), §6.2.1, pp. 70–71 (table) and §6.2.2, p. 72 (LRR derivation). I checked the
  transcription against the p. 71 page image today.

## 0. Verdict

| target | result |
|---|---|
| T1 calibration | **REPRODUCED, 15/16.** Dictionary: y = −q, w₁ = 1/s, w₂ = 1/t. Points: M=∅↔e24, {1}↔e23, {2}↔e14, {1,2}↔e13. Labels: I=∅↔24, {1}↔23, {2}↔14, {1,2}↔13. Twist: λ ∈ (0, ½], i.e. ⌈λ⌉ = ⌈2λ⌉ = 1, **α = 1** |
| T2 odd entry | **(24)\|e13 = Ã_{{1,2},∅}**, exactly the crown prediction. It sits on the long-root T-curve M ↔ M∖{1,2}, which is the curve e24–e13 of weight st |
| T3 ratio | On H₁₂ = {w₁w₂ = 1} = {st = 1}: Motegi/Koncki = **g(q, w₁) exactly**. Globally: **Koncki = Motegi + q(1−q)(1−w₁⁻²)(1−w₂⁻²)** |
| crown 2 | **Holds at N=2.** All 64 structure constants of Koncki's basis are Laurent in (t, s): they have only monomial denominators, at λ=¼ and at λ=¾. Motegi's basis has exactly one extra denominator, (st − 1) = H₁₂ |
| controls | Swapped dictionary (w₁ = 1/t): 8/16. λ ≤ 0: 10/16. λ ∈ (½, 1]: 13/16. λ ∈ (1, 3/2]: 10/16. All of them break |

**Structural explanation (proved at N=2):** 15/16 is the *maximum* agreement any honest K-theory table can have with
Motegi's, and the 16th entry is forced. Lemma QC says Motegi's I=∅ column fails GKM on exactly one T-curve,
(e24, e13), by the factor g ≠ 1, and his other three columns are GKM. So the "correction" Motegi mentions is
Lemma QC.

## 1. The computation

Motegi's N=2 table (engine; Ã_{M,I}, rows I):
- I=∅: M=∅: (1−qw₁²)(1−qw₂²)(1−qw₁w₂); {1}: (1−q)(1−qw₂/w₁)(1−qw₂²); {2}: (1−q)(1−qw₁/w₂)(1−qw₁²);
  {1,2}: (1−q)²(1−q/(w₁w₂)).
- I={1}: {1}: (1−w₁⁻²)(1−qw₂/w₁)(1−qw₂²); {2}: −(1−q)(1−w₂⁻²)(1−qw₁²)·…; {1,2}: … (see `motegi_table.py` output).
- I={2}, I={1,2}: likewise. Triangular for the chain ∅ < {1} < {2} < {1,2}, which is the BB order e24 > e23 > e14 > e13
  read upside down.

The dictionary is read off the diagonals. Ã_{∅,∅} = Π(1 − q·w) over w ∈ {w₁², w₂², w₁w₂} must match
mC(24)|e24 = Π(1 + y/w') over w' ∈ {t², s², st}, which forces y = −q and {w₁, w₂} = {1/t, 1/s}. Then
Ã_{{1},{1}} = (1−w₁⁻²)(1−qw₂/w₁)(1−qw₂²) against mC(23)|e23 = (1−s²)(1+y/t²)(1+ys/t) forces w₁ = 1/s.

With this dictionary, `compare.py` checks all 16 entries as exact rational identities, for
λ ∈ {−⅓, 0, ¼, ½, ⅔, 1, 5/4, 3/2} (generic λ and the degenerate λ = 0) and both {w₁, w₂} assignments. Ratio
= 1 on 15 entries iff α = 1 and ⌈λ⌉ = 1. The only odd entry:

  Motegi: Ã_{{1,2},∅} = (1+y)²(1+sty)
  Koncki: mC(24, λ)|e13 = (1+y)(st(y+y²) + 1 + s²y + t²y − s²t²y)   (α = 1, ⌈2λ⌉ = 1)
  Koncki − Motegi = −y(1+y)(1−s²)(1−t²) = q(1−q)(1−w₁⁻²)(1−w₂⁻²).

Note (1−s²)(1−t²)(1−st) = [point class at e13] = Ã_{{1,2},{1,2}}. So
  Koncki_∅ = Ã_∅ + (q(1−q)/(1 − 1/(w₁w₂))) · Ã_{{1,2}},
a correction by the point class with a **simple pole on H₁₂**. That is the h⁻¹ step in Prop 5 of the Day 243 file,
which takes the lattice L_g to L_1.

**T3.** On st = 1, (Motegi/Koncki)(24|e13) = s²(1+y)²/((s²+y)(1+s²y)) = (1−q)²/((1−qw₁²)(1−qw₁⁻²)) = g(q, w₁),
checked symbolically (`t3b.py`). This equals Motegi(e13)/Motegi(e24) on st = 1, i.e. Lemma QC, and Koncki's column
has ratio 1 there.

## 2. Why 15/16 is forced (Proposition, N=2, proved)

**Prop K.** Use the dictionary above. Let f be any honest class in K_T(LG(2,4))[y], i.e. any GKM-consistent
4-tuple. Suppose f agrees with Motegi's I=∅ column at e24, e23 and e14. Then f(e13) ≠ Ã_{{1,2},∅}. Moreover,
f(e13) ∈ Ã_{{1,2},∅} + q(1−q)(1−w₁⁻²)(1−w₂⁻²) + (point class)·Z[t^±, s^±][y]. Koncki's α=1 entry is that
representative.

*Proof.* The T-curves of LG(2,4) are the reflections of the signed-permutation group: 24–14 and 23–13 (weight t²),
24–23 and 14–13 (weight s²), 24–13 (st), 23–14 (t/s). Lemma QC at N=2 (proved all N, Day 243) gives, on st = 1,
Ã_{{1,2},∅} = g·Ã_{∅,∅} with g ≢ 1. An honest f satisfies f(e13) ≡ f(e24) = Ã_{∅,∅} mod (1 − st). So
f(e13) ≠ Ã_{{1,2},∅}. For the coset, f(e13) is pinned modulo the ideal of the three congruences, which is generated
by the point class. Koncki's entry satisfies all three congruences (`gkm.py`, `branches.py`), and its difference from
Motegi is the displayed term. The other three Motegi columns pass all six GKM congruences (`gkm.py`; for I ≠ ∅ the QC
defect is g·Ã_{∅,I} = 0). ∎

So Motegi's "correction" is exactly the Lemma QC defect at N=2, and it lives on the unique pair {1,2}.

## 3. Crown prediction 2 at N=2: the geometric basis has no H₁₂ poles

`sc.py`: all 4×4×4 structure constants of the two bases are computed by top-down triangular solve.
- Koncki, λ = ¼ (α=1): the denominators that involve t or s are {t, s} only. **Laurent.**
- Motegi: denominators {t, s, st − 1}. The (st−1) is Clio's c^{(1,2)}_{∅∅} pole, i.e. Q488's long-root pole.
  For example, c^{13}_{24,24} has denominator s³t³(st − 1) for Motegi and s³t³ for Koncki.
- Detector control (`sc_control.py`): at λ = ¾ (α = 0), Koncki's *printed* (23)|e13 gives poles {t±1, s+1}. The
  LRR-corrected entry (see §4) gives Laurent. So the detector fires on non-GKM input, and crown 2 holds at a second,
  α = 0, twist too.

General principle, standard GKM: a triangular basis of honest classes whose diagonal is (repelling Euler class) ×
(units 1 + y·monomial) spans the GKM module, so its structure constants have no (1 − T-weight) poles. Q488's
long-root poles therefore measure the distance between L_g and L_1. That is the normalisation defect, as the
Day 243 dream crown said.

## 4. Side finding: a typo in Koncki's thesis (α = 0 branch, not Motegi's entry)

Printed p. 71: mC(23,λ)|e13 = (1−st)(1+y)t^{1−⌈2λ⌉}(s(1+y) + t + sty) for α = 0. Koncki's own LRR sum (p. 72,
displayed formula) evaluates to (1−st)(1+y)t^{1−⌈2λ⌉}(s(1+y) + t + **s²t**y). The printed version fails GKM on
the curves 23–13 (t²) and 14–13 (s²). The LRR version passes all of them. The α = 1 branch agrees with LRR. I
checked the page image: this is the thesis's typo, not my transcription. It is irrelevant to Motegi's 15/16,
because his match lives at α = 1.

Also an own-error note: the first `compare.py` run had −s²t²y² for −s²t²y in Koncki's α=1 entry. That was a typing
slip by me, which made Koncki's table look non-GKM. I caught it by cross-checking two scripts, and all the numbers
above come from the corrected script.

## 5. What this does NOT prove, and the all-N hunch

- N=2 only. The all-N statement would be: Motegi's Ã (y = −q, w_i = inverse torus coordinates) equals Koncki–Weber
  twisted mC of LG(N, 2N) at small positive λ, except at the entries carrying a Lemma QC defect. Those are the pairs
  (M, I) with some a<b in M and M∖{a,b} ≥ I: 1/16 at N=2, **7/64 at N=3**, 36/256 at N=4. This is a HUNCH, and the
  correction may propagate beyond those entries.
- An N=3 test needs the twisted mC classes of LG(3,6), via the Koncki–Weber Hecke action (2301.12746, not
  deep-read) or AMSS Demazure–Lusztig operators at λ → 0⁺.
- For an all-N proof, a candidate is to show the KW-twisted classes at λ ∈ (0, ½] satisfy the same recursion as
  Motegi's away from the long-root curves (MNS Chevalley 2312.17200 as the bridge). Not attempted.

## Verification summary
- 15/16 calibration: exact symbolic identities, 8 λ values × 2 dictionaries.
- GKM: all 6 curves (and t = −1, s = −1 for the long roots), all 8 columns.
- T3: exact symbolic identity on st = 1.
- Structure constants: exact, 2 bases + 2 controls.
