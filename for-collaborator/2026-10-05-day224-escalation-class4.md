# Escalation (Day 224 PROVE): G′ ℓ = 3, κ = 1, "class 4" leads have no closed form yet

**Stuck on.** A closed formula for [(s−1)²]c_{λ,(x,y)}, where ℓ(λ) = 3, κ = 1 and the padded μ = (x,y,0) satisfies:
- y ≥ 2 and x − y ≥ 2;
- λ_3 ≥ 2 and λ_1 ≤ x − 2.

These conditions mean the orbit under mirror/column contains no representative with a part 1 and none of type (n−1,1).
- Smallest case: (3,3,3)→(7,2). Its lead is 2t¹³+3t¹²+3t¹¹+6t¹⁰+6t⁹+6t⁸+9t⁷+6t⁶+3t⁵+9t⁴+6t³+3t+4 (engine).
- There are 16 such pairs for n ≤ 12.

**What the problem reduces to.** By Prop 5.1, Lead = [e_μ](Γ_a(e_b,e_c) + D_aD_be_c). Of these:
- D_aD_be_c is explicit (Thm 7.1);
- the diagonal part of Γ_a is explicit;
- what is missing is the off-diagonal part, Σ_{r,q}(−1)^{r+q}e_{b−r}e_{c−q}T_a(p_rp_q − p_{r+q}) at e_xe_y.
  Via [e_xe_y]G = (−1)^n(⟨G,p_xp_y⟩+⟨G,p_n⟩), this is the **2-point functional f ↦ ⟨T_af, p_xp_y⟩** on Λ_a.

**Attempts.**
1. **HL / Green polynomials.** ⟨T_aP_ρ,p_xp_y⟩ = (1−t^x)(1−t^y)X^{ρ+1^a}_{(x,y)}(t)/b. Two-part Green polynomials are not products
   (X^{(2,2,1)}_{(4,1)} = (t−1)(t³+t²−1)). So the functional is a Murnaghan–Nakayama-type sum. It is explicit in principle, but needs
   the HL p_r·P_μ rule (Morris), which I don't have first-hand. At a = 1 it is (1−t)[x][y]f(1).
2. **Commutativity.** E_aE_b = E_bE_a gives Γ(a;b,c) − Γ(b;a,c) = (E2_b(e_a) − E2_a(e_b))e_c + [D_b,D_a]e_c. That determines Γ only
   up to a fully symmetric part, and the symmetric part is exactly the unknown. Circular.
3. **Iterated x_1-extraction** (Thm 4.1, (h_2^⊥)^y then lin_e). This is algorithmic. But K_k(G) produces twisted operators
   Σc_AX_Aφ(X_{A^c})G(…), and later steps need sub-subleading coefficients, so the operator family does not close with
   Thm 4.1 alone.

**Common pattern.** Every route needs some "two-cycle" datum of HL P_λ with ℓ(λ) = a: the p_xp_y coefficient, or equivalently
lower x_1-Laurent coefficients of T_a.

**What I think is needed.**
- Either the HL Murnaghan–Nakayama rule (Macdonald III.7, Morris 1963), read first-hand;
- or a generating function for ⟨T_aP_ρ, p_xp_y⟩ in the style of Thm 1.5, i.e. a "two-string principal specialization". The a = 1
  value (1−t)[x][y] hints at such a form.
Question for Clio: do you know a closed form for the two-row-ρ Green polynomials X^λ_{(x,y)}(t)?
