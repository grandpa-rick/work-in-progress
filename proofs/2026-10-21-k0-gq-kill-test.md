# K0 kill test: Ψ_{κ,0}(M) = Ikeda–Naruse GQ at a fixed point (Wake 249)

**Verdict: PASS, with a dictionary correction.** The match holds with **λ = κ** (reversed to decreasing order), **not λ = κ−1** as the
Day 248 connection file says. Grade: **computed**.

## Statement tested
For κ_1<…<κ_r, M ⊆ [N], |M| = n:
  Ψ_{κ,0}(M) = GQ_λ(x_1,…,x_n | b) at β = −1, with λ = (κ_r,…,κ_1), x_m = 1 − 1/w_m (m ∈ M), b_k = 1 − w_k.
In multiplicative form X = 1 − x, B = 1 − b, this reads X_m = 1/w_m and B_k = w_k.

## Sources (the printed formula)
- IN 1112.5223 Def 2.1 (2.2), fetched and pdftotexted: [[x|b]]^k := (x⊕x)(x⊕b_1)⋯(x⊕b_{k−1}), and
  GQ_λ = 1/(n−r)! Σ_{w∈S_n} w[ [[x|b]]^λ Π_{i≤r} Π_{j=i+1}^{n} (x_i⊕x_j)/(x_i⊖x_j) ].
- Motegi 2610.15344 Def 2.3 (2.17) (local text, L460–557) agrees: Π_{k=1}^{λ_i−1}. Motegi Prop 2.5 has no shift either (i_j = λ_{r+1−j}).
- Why the file's shift is wrong: the factorial factor of GQ has λ_i − 1 factors, i.e. k < λ_i, and our ρ_κ = Π_{k<κ}. So λ_i = κ_j.
  The file wrote "Π_{k<λ_i} ↔ ρ_{κ_j} with λ_j = κ_j − 1", which double-counts the −1.
- Orientation: IN's head order (λ decreasing) ↔ the engine's head order reversed. The pair denominator 1 − x_{lower}/x_{higher} then agrees.

## Implementation
`scripts/wake249/k0_gq_check.py`. The GQ side is written from the printed (2.2) and sums over all of S_n with the 1/(n−r)!.
The engine side is `work-in-progress/scripts/prove248/psi.py` (the T⁰ definition). The log is `scripts/wake249/k0_gq_check.txt`.

## Results (main dictionary, λ = rev κ)
| κ | M | method | result |
|---|---|---|---|
| (1,2) | (1,2,3,4) | sympy cancel(diff) | EQUAL |
| (2,3) | (2,3,4,5) | symbolic | EQUAL |
| (1,3) | (1,2,3,4) | symbolic | EQUAL |
| (2,4) | (2,3,4,5) | symbolic (generic: non-consecutive, no label 1) | EQUAL |
| (2,5) | (1,3,5,6) | 3 random rational points, exact | EQUAL |
| (3,4) | (2,4,5,7) | 3 pts | EQUAL |
| (2) | (1,2,4,5) | 3 pts (r=1) | EQUAL |
| (1,3,4) | (1,2,3,4,5) | 3 pts (r=3) | EQUAL |

The match is exact, with no factor.

## Controls (each must break the claim)
- **λ = κ − 1** (the connection file's dictionary): DIFF in all 5 applicable cases, with a non-constant ratio.
- **X = w** (x_m = 1 − w_m): DIFF in all 8 cases.
- **B = 1/w** (b_k = 1 − 1/w_k): DIFF in 7 of 8. It agrees at κ = (1,2), which is a degenerate point because Ψ_{[r],0} has no w_1 (t3.py).
  It is separated at every other κ.
- **Twist θ ≠ 0** (κ = (2,4), θ = (1,0) and (1,1)) against untwisted GQ: DIFF under every dictionary. This is expected, since twisted Ψ is not plain GQ.
  The moral claim twisted = GQ·[O(θ)] remains untested (that is K1).

## Consequences for the connection file
- Fix the dictionary to **Ψ_{κ,0}(M) = GQ_κ(X_M | B)** with no shift. The P1 prediction there says "κ=[r] ↔ zero part". Under the corrected
  dictionary κ = [r] ↔ λ = (r,…,1), so P1 involves no zero part. That prediction (item 4) and the Lemma Δ "staircase with zero part"
  remark (item 3) need to be restated before K1.
