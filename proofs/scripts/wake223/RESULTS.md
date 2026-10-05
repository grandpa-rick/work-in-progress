# Wake 223 — 2-point test for B, full B(e_a,e_b), and G′ at n = 6 (grade: everything COMPUTED)

B is as in Day 220 Thm 1: B(f,g) = ∂_s(f⋆g)|_{s=1}, and B(e_a,e_b) = M_ab = D_a(e_b) with D_a = Σ_{|A|=a} c_A X_A Δ_A.
L(a,b) := (1−t^{a+b})(t^{ab}−1)/((1−t^a)(1−t^b)) = −[a+b]_t[ab]_t/([a]_t[b]_t).

## Task A
Scripts and logs:
- `bfull.py`: D_a(e_b) with t symbolic, N = a+b variables, full e-expansion. Logs `bfull_n6.log`, `bfull_n7.log`; data `bfull_n7.pkl`.
- `crosscheck_star.py`: independent check against the full s-dependent ⋆ engine (day220/biderivation.py), ∂_s at s=1,
  t = 3/5, n ≤ 6. **9/9 OK.** Log `crosscheck_star_engine.log`.
- `checks.py`: log `checks.log`.
- `second_order.py`: log `second_order_n6.log`.

1. **[e_{a+b}] B(e_a,e_b) = L(a,b): 12/12** (all a ≥ b ≥ 1, a+b ≤ 7). The q-integer form −[a+b][ab]/([a][b]) agrees **12/12**, and the two forms are identical as rational functions. (1,1) gives −(1+t).
2. **B spreads. Full support formula (computed, 12/12, a+b ≤ 7), for a ≥ b:**
   B(e_a,e_b) = b·e_a e_b + Σ_{j=1}^{b} L(a−b+j, j) · e_{a+j} e_{b−j}.
   The support is exactly {e_{a+j}e_{b−j} : 0 ≤ j ≤ b}, which is always a product of at most two e's. It is not just e_{a+b}.
   Example: B(e3,e3) = 3e3e3 + L(1,1)e4e2 + L(2,2)e5e1 + L(3,3)e6.
3. **Principal specialization:** (1−t^n)/∏(1−t^{λ_i}) = ∏ps(p_{λ_i})/ps(p_n) with ps(p_k) = 1/(1−t^k). Checked 65/65 for n ≤ 8, but this is a tautology.
4. **ℓ = 3.** Theorem G was re-checked on 16/16 cases (n ≤ 8): Lead = pref·(w_ab w_ac + w_ab w_bc + w_ac w_bc + w_ab w_ac w_bc) with w_ij = t^{λ_iλ_j} − 1.
   The (s−1)^2 coefficient of e_n in E_aE_bE_c(1) splits into two pieces, both computed directly (7/7, n ≤ 6):
   - Iterated B: [e_n] D_aD_b(e_c) = L(b,c)L(a,b+c) = pref·w_bc(w_ab + w_ac + w_ab w_ac). These are the connected graphs that contain the edge bc.
   - Second-order piece: [e_n] E_a^{(2)}(e_b e_c) = pref·w_ab·w_ac. This is the path centred at a.
   **So B alone does NOT generate the cumulant. The missing graph comes from the order-2 operator E^{(2)}.**

## Task B
`gprime_check.py` uses wake222 `analyze.leads`. Logs `analyze_n6.log` (the stock analyze.py run on n = 6) and `gprime_check.log`. The wake222 `leads_all.pkl` was restored after the run, and the n = 6 copy is saved as `leads_n6.pkl`.
- (i) v = ℓ − κ: **43/43 at n = 6** (1^6 not in the pickle) and **78/78 for n ≤ 6**. No mismatch. Every coarsening matches the history leads ("hist=OK" in every coarsening line).
- (ii) The 9 non-coarsenings at n = 6:
  - (4,2)→(5,1): −[4]
  - (3,3)→(5,1): −(1+t²)²
  - (3,3)→(4,2): −(1+t)
  - (3,2,1)→(4,1,1): −[3]
  - (2,2,2)→(3,2,1): −3(1+t)
  - (2,2,2)→(3,3): (t+2)[3]
  - (2,2,2)→(4,1,1): (t+2)[3]
  - (2,2,2)→(5,1): 2t⁷+3t⁶+6t⁵+6t⁴+9t³+6t²+3t+4 (irreducible over ℚ; value 39 at t = 1)
  - (2,2,1,1)→(3,1,1,1): −(1+t)
- (iii) The fit (t^{λ_recv·k}−1)/(1−t^k) **FAILS**:
  - (3,3)→(4,2): predicted −[3], actual −(1+t).
  - (3,3)→(5,1): predicted −[3]_{t²}, actual −(1+t²)².
  - It passes on the other ℓ−κ=1 cases, with (2,2,2)→(3,2,1) passing only when counted with multiplicity 3.
  - Replacement rule: a single move in which part b gives j cells to part a ≥ b has weight L(a−b+j, j), i.e. the receiver's excess plus the cells moved, paired with the cells moved. This follows from item 2. **The first-order B-rule, [(s−1)^1] e⋆_λ = Σ_{i<j} M_{λ_iλ_j} e_{λ∖{i,j}}, matches 78/78 pairs for n ≤ 6:** the lead agrees when v = 1, and the coefficient vanishes when v ≥ 2.
  - The ℓ−κ = 2 non-coarsenings were not fitted.
- (iv) (2,2,2)→(5,1): v = 2, **Lead(0) = 4. PASS.**

Not computed: B in full for a+b > 7; second-order splitting for n > 6; any G′ fit for the ℓ−κ = 2 non-coarsenings.
