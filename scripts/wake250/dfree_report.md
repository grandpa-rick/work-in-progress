# Wake 250: is there a Motegi-free route to (D)? (dfree1/2/3.py, log dfree_out.txt)

Grade: **computed** throughout. Exact Fractions at the points w=(5,7,11,13,17) and random rationals (seeds 3, 7).

## Independent classes
`kschub.py` computes [O_{X^v}]|_u on Sp_{2n}/P (P omits α_n) by the Graham–Willems subword formula
Σ_{J: Hecke product = v} (−1)^{|J|−ℓ(v)} Π_{j∈J}(1 − e^{β_j}). It uses no Ψ and no F.
- Validation: the support is Bruhat order; χ_T(O)=1; and χ_T(O_{X_u}·O_{X^v}) = [v≤u] for all u, v (n ≤ 5). This holds only for one sign choice in Atiyah–Bott.
- Schubert classes are X_u = w0·X^{w0u}, where w0 = −1 acts by M ↦ M^c and t ↦ 1/t.

## 1. Dictionary
Positions are reversed: κ ↔ −1 at positions p = n+1−k, and e^{ε_p} = w_{n+1−p}. This is Theorem D's matching.
- **Ψ_{ν,0}(M; 1/w) = [O_{X^ν}]|_M exactly, with no normalisation.** Holds for n = 2, 3, 4, 5 and all (ν, M).
- Hence Ψ_{κ^c,0}(L^c; w) = [O_{X_κ}]|_L exactly.
- Negative controls fail. Ψ(·; w) gives 6/16, 27/64, 110/256 and 430/1024 mismatches, and only 2/2^n of the ν are even proportional. Identity positions fail the support test.
- This resolves the Day 249c `ideal.py` mismatch: the wrong base/orientation was used there.

## 2. Chain
- The χ_T-dual basis ξ^v of {O_{X_u}} equals Σ_{x≥v} μ(v,x) O_{X^x}, where μ is the Möbius function of dominance (= (−1)^{Δℓ} or 0).
- **w_M² F_κ(M) = cc(κ) ξ^κ(M)**: 0 mismatches out of 16, 64, 256 and 1024.
- χ_T(O_{X_κ}·w²F_λ) = cc(κ)δ_{κλ}: 0 mismatches.
- The explicit Motegi ℓ_κ(e_L) equals O_{X_κ}|_L w_L²/(cc(κ) eul_L): 0 mismatches for n ≤ 5. It splits into two parts:
  - (i) the dictionary above (with w0);
  - (ii) the one-point identity cc(L)/(F_L(L)Ψ_{L^c,0}(L^c)) = w_L²/eul_L, which compares single-term products.
- F is defined explicitly in eng.py, not through (D): F_κ = w_M^{2r} Ψ_{κ,(κ_j−2j)}.
- **Recursion:** cc(λ)ξ^λ, computed geometrically, satisfies Lemmas A_m(a,b,c) and B(a,b) coefficient for coefficient (8, 24 and 64 instances out of 8, 24 and 64). It agrees with g_[n] at the top. These recursions run downward and determine every class from the top one.

## Where the proof stands
Caution: the proof of Theorem D in §1bis establishes Λ̃_κ(F_λ) = δ without Motegi. Its last step, "Λ̃_κ = ℓ_κ by uniqueness of the dual basis", uses (D) unless ℓ_κ means the abstract dual functional.

To prove (D) for the *explicit* ℓ_κ without Motegi, three inputs are needed beyond §1bis:
- (a) **AB/Demazure:** a^κ_L = [O_{X_κ}]|_L / eul_L. This is standard (Kumar; Graham–Kumar).
- (b) **Dictionary:** Ψ_{ν,0}(·; 1/w) = [O_{X^ν}]. This follows from Ikeda–Naruse's theorem plus K0. Alternatively it could be proved in-house as a θ=0 analogue of Lemmas A/B.
- (c) **The one-point identity (ii).** This is elementary.

The lemma "F = twisted ideal-sheaf basis" is not a missing piece. It already follows from §1bis together with Poincaré duality.
