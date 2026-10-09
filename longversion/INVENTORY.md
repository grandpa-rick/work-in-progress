# arXiv long version: inventory of results at or above the trust boundary

Built 2026-10-09 from the registries `proofs/registry/hikita-star-*.json` (6 files) and the FPSAC draft
`work-in-progress/fpsac2027/fpsac2027-draft.tex`. Grades are copied exactly from the registry; nothing was upgraded.

**Conventions**
- Registry file abbreviations: **DS** = hikita-star-dominance-support, **EK** = hikita-star-ek-er, **TC** = hikita-star-two-column,
  **E2** = hikita-star-e2-e2, **E3** = hikita-star-e3-er, **E4** = hikita-star-e4-er.
- Proof-file paths are relative to `/home/agent/projects/`, and sizes are in bytes. Frequently used files:
  - `207b` = proofs/2026-09-26-day207b-ek-star-er-general-k-PROVED.md (24794)
  - `206b` = proofs/2026-09-25-day206b-W_r-proved.md (13029)
  - `205L` = proofs/2026-09-25-day205-L1-L4-residue-proof.md (8935)
  - `205Z` = proofs/2026-09-25-day205-sub-lemma-Z-reduction.md (7031)
  - `209` = proofs/2026-09-29-day209-two-column-TC-PROVED.md (14833)
  - `212` = proofs/2026-09-30-day212-ell-column-PROVED.md (11791)
  - `214` = proofs/2026-09-30-day214-DS-all-lengths-PROVED.md (21884)
  - `215` = proofs/2026-10-01-day215-theorem-H-s0-limit-is-HL.md (14767)
  - `216b` = proofs/2026-10-01-day216b-theorem-H-prime-nabla-transport.md (28832)
  - `217e` = proofs/2026-10-02-day217e-boundary-of-the-st-square.md (17697)
  - `220` = proofs/2026-10-03-day220-s1-carre-du-champ.md (30753)
  - `221` = proofs/2026-10-05-day221-lead-cumulant.md (8075; restored after the 1.2 KB truncation; §6 "W recheck" is marked LOST in the file)
  - `223` = proofs/2026-10-06-day223-G-by-vertex-deletion.md (19300)
  - `224` = proofs/2026-10-06-day224-Gprime-second-order.md (21670)
  - `225` = proofs/2026-10-07-day225-class4-hopf-route.md (23199)
  - `226r` = proofs/2026-10-07-day226-cold-recheck-class4.md (14326; recheck, not referenced in a `file` field)
- **Size flag:** none of the `file` targets of proved or checked-sober nodes is missing or under 3 KB. The ⚑ rows below flag the exceptions that do exist: a recheck file under 3 KB, files present in only one tree, files that are not proofs, and a node with no file.
- "FPSAC" gives the `\label` in fpsac2027-draft.tex, or "—" if the result is not typeset there.

---

## 1. Results graded proved / checked-sober (no lean-verified nodes exist)

### (A) Definitions, ⋆ setup, Hikita background, tools

| # | Result | Registry :: node id(s) | Grade | Proof file (size) | FPSAC | LV § |
|---|---|---|---|---|---|---|
| A1 | **Subset formula / parabolic kernel.** (A_k) e_k(Y)F = t^{C(k,2)} σ^{(k)} π^k F on symmetric F, and (K_k) σ^{(k)}G = Σ_{\|A\|=k} G^{(A)} ∏^x_A. Hence e_k⋆F = Σ_A c_A X_A F(X_{A^c}, sX_A). | EK::Ak-parabolic-kernel-all-k, EK::Kk-kset-kernel-all-k, DS::subset-formula-Ak-Kk, E2::e2Y-parabolic-kernel (k=2), E3::ek-parabolic-kernel-k3, E4::ek-parabolic-kernel-k4, E4::ek-parabolic-kernel-k5 | proved (all) | 207b (24794); k=2 in 206b (13029) | prop:subset | A.2 |
| A2 | **⋆-reading interface (R0).** e_k⋆e_r = t^{-C(k,2)} e_k(Y)•e_r. 𝔮 bijective = Hikita Lemma 3.1 + Cor 3.9; the older locator "Def 3.4 / Lemma 3.3" was wrong (corrected wake 230). | EK::hikita-def34-lemma33-interface | proved | 207b | text of §2 (no label) | A.1 |
| A3 | **HL residue lemmas.** Lemma 1: σ_m equals the HL kernel. Lemma 2: (1−t)Σ_a X_a^n ∏ a_aj = q_n. Two-variable Jing form: (1−t)²S_{n,p} = (QJ(n,p)+QJ(p,n))/(1+t). | EK::lemma-1-2-day205b, EK::lemma-1-2-day205b-2, E3::lemma-2-HL-residue-stub, E4::lemma-2-HL-residue-stub, E2::two-var-HL-residue-H | proved | 205L (8935); 206b | — | A.3 |
| A4 | **Hikita Thm 3.12, re-derived.** e_1(Y)•e_n = (1−s)[n+1]e_{n+1} + s e_1e_n. CITE: this is Hikita's theorem. | DS::ds-r11-thm-3-12-operator-form | proved | 205L (8935) | cited (cor:pieri remark) | A.3 (cite) |
| A5 | **Taylor pieces at s=1.** F(X_{A^c}, sX_A) = s^{Δ_A}F; E_k^{(p)} preserves Λ_N and has order ≤ p. | DS::thmA-premise-subset-formula-taylor | proved | 220 (30753) | eq:taylor | A.4 |
| A6 | **(KF) HL formula for E_k.** E_kF = Σ_{ℓ(ρ)≤k} P_{ρ+1^k}(x;t)(Q'_ρ[(s−1)X;t])^⊥ F. | DS::KF-hall-littlewood-formula-for-Ek, DS::thmW-premise-KF | proved | 216b (28832) | — (only the factorisation T_kP_ρ = P_{ρ+1^k} appears, in the thm:W proof) | A.4 |
| A7 | **Textbook inputs, stated with locators (CITE).** Macdonald I (1.11), (6.6)–(6.7) (dominance via e→m); III (3.2) HL vertical Pieri; III (2.15) raising operators. | DS::macdonald-dominance-monomial-textbook, DS::macdonald-III-3.2-HL-vertical-pieri, DS::thmC-premise-macdonald-raising-operators | proved (textbook) | 214; 215; 220 | cited inline | A.5 (prelims) |

### (B) Pieri package and Dominance Support

| # | Result | Registry :: node id(s) | Grade | Proof file (size) | FPSAC | LV § |
|---|---|---|---|---|---|---|
| B1 | **e_k⋆e_r Pieri, all k, r.** e_k⋆e_r = Σ_b s^b F_{k−b}(t^{r−b}) e_b e_{r+k−b}, with c(n,j) = (s;t)_{n−j}/(t;t)_{n−j}(α_j − s t^{n−j} α_{j−1}). Includes the GF (E_k) by outer-peel induction and the M^{(k)}_{jb} matrix. | EK::root, EK::ek-star-er-pieri-all-k, EK::Ek-gf-induction-all-k, E4::ek-star-er-general-M-formula, TC::ek-star-er-theorem-stub, DS::Bmatrix-premise-207b-pieri | proved | 207b (24794) | Same rule, other form and route: cor:pieri (via thm:lin) | B.1 |
| B2 | **Plethystic linear coefficient + 2nd proof of Pieri.** lin_e(e_k⋆G) = (−1)^d [k+d]/[k] · G[(s−1)[k]_t]; Cor 3.2 gives the ℓ=2 Pieri rule. | DS::plethystic-lin-e-and-pieri-reproof | proved | 224 (21670) | thm:lin, cor:pieri | B.1 (2nd proof) |
| B3 | **Small cases (worked examples).** e_2⋆e_2, e_2⋆e_r (W_r), e_3⋆e_2, e_3⋆e_3, e_3⋆e_r full closed form, e_4⋆e_3, e_4⋆e_4, e_4⋆e_5, e_4⋆e_r full closed form; all r. | E2::e2-star-e2-closed-form, E2::e2-star-er-pieri-conjecture, E2::W_r-GF-extraction-all-r; E3::hikita-star-e3-e2-closed-form, E3::hikita-star-e3-e3-closed-form, E3::hikita-star-e3-er-pieri-conjecture-full-closed-form, E3::hikita-star-e3-er-analytic-proof, E3::e3-star-er-closed-form-all-r; E4::hikita-star-e4-e3-via-commutativity, E4::hikita-star-e4-e4-closed-form, E4::hikita-star-e4-e5-closed-form, E4::hikita-star-e4-er-pieri-full-closed-form, E4::hikita-star-e4-er-analytic-proof, E4::e4-star-er-closed-form-all-r | proved | 207b; W_r in 206b (13029) | — | B.1 examples |
| B4 | **t=0 Pieri is a truncated geometric law.** e_k⋆e_r = Σ_{b<μ}(1−s)s^b e_b e_{r+k−b} + s^μ e_μ e_M. | E4::hikita-star-t-zero-specialization | proved | 207b §8 | — | B.1 corollary |
| B5 | **Two-column rule (TC).** Closed GF for t^{-C(k,2)} e_k(Y)•(e_a e_b), all k, a, b; closing identity (★2), step map, (R2) recursion. | TC::root, TC::two-column-gf-rule, TC::two-column-closing-identity, TC::two-column-step-map, TC::two-column-R2-recursion | proved | 209 (14833) | — (remark defers it to the long version) | B.2 |
| B6 | **TC is straightening-free.** Only chain shifts e_{b0}E(t^iz)E(t^jw) occur, so the support has ℓ ≤ 3. The "unbounded tail" caveat is COMPUTED only. | TC::two-column-straightening-free | proved (corollary of B5) | ⚑ proofs/2026-09-29-day208-two-column-k2-data.md (5327): a data note; the argument is the corollary in 209 | — | B.2 remark |
| B7 | **ℓ-column rule (★ℓ).** Pairwise cross kernel ∏_{c<c'} K_{i_c i_c'} with an explicit prefactor; closing identity (Z) by the residue theorem; ℓ-column step map. | TC::ell-column-rule, TC::ell-column-closing-identity-Z, TC::ell-column-step-map | proved | 212 (11791) | — | B.3 |
| B8 | **Dominance Support, all lengths + Op-DS.** e⋆_λ = s^{n(λ)}e_λ + Σ_{μ▷λ} c_{λμ}e_μ, with c ∈ Q[s,t] and c(1,t)=0; e_k⋆e_μ = s^{Σmin(μ_i,k)} e_{μ∪k} + higher. | DS::root, DS::ds-all-lengths-degree-count, DS::thmA-premise-day214-degree-lemmas | proved | 214 (21884) | thm:DS | B.4 |
| B9 | **Exact up-set support and s-valuation.** s^{n(μ)} \| c_{λμ}, c(s,0) = s^{n(μ)}(1+O(s)), so the support is the full up-set and val_s = n(μ). Includes the Gauss-valuation Lemma A. | DS::full-upset-support-exact-valuation, DS::day214-gauss-valuation-lemma-A | proved | 214 | thm:DS | B.4 |
| B10 | **Earlier/special DS cases (superseded by B8; keep as remarks).** Length-2 = SP with explicit coefficients; DS(r,1,1) all r; DS ⇒ ≤ min(a,b)+1 terms. | DS::ds-length-2-slice-is-SP, DS::ds-lambda-r11-all-r, DS::ds-corollary-meta-conjecture | proved | proofs/2026-09-30-day213-DS-length2-from-207b.md (6280); proofs/2026-09-26-day207-DS-r11-proved.md (8324); proofs/2026-09-16-day196-dominance-support.md (9605) | — | B.4 remarks |
| B11 | ⚑ **"DS length-3 / length-4 empirical"** (λ=(2,1,1), (3,1,1), (2,2,1); (1^4), (2,1,1,1)). Graded proved, but the node text is test evidence; superseded by B8. Do not cite as separate results. | DS::ds-length-3-empirical, DS::ds-length-4-empirical | proved | proofs/2026-09-16-day196-dominance-support.md (9605) | — | drop |
| B12 | **DS, second proof from (N).** c_{λμ} = t^{-n(λ')} Σ A·T·B. Gives the valuation side only (no polynomiality, no t=0 regularity). ⚑ Computer check n≤4 only; the file still contains the unfilled n=5 placeholder text. Inherits (N)'s imports. | DS::ds-from-N-second-proof | proved | proofs/2026-10-02-day218-DS-from-N.md (15868) | — | C.3 remark |

### (C) Limits / edges of the (s,t)-square and ∇-transport

| # | Result | Registry :: node id(s) | Grade | Proof file (size) | FPSAC | LV § |
|---|---|---|---|---|---|---|
| C1 | **Theorem H (s→0).** On b_μ = s^{n(μ)}e_μ, E_k mod s is conjugate via b_μ ↦ t^{-n(μ')}P_{μ'}(x;1/t) to t^{-C(k,2)}e_k (HL e_k-Pieri). | DS::theorem-H-s0-star-is-HL-pieri | proved | 215 (14767) | rem:H (no-novelty remark) | C.1 |
| C2 | **The d-matrix.** d_{λμ}(t) = [s^{n(μ)}]c_{λμ} = t^{-n(λ')}⟨e_λ, H̃_{μ'}⟩: the classical e→HL transition matrix; d(0)=1; d(1) = #0-1 matrices. CITE Kirillov/Knuth for positivity and the count. | DS::d-lambda-mu-t-count-01-matrices, DS::d-equals-hall-littlewood-transition | proved | 215 | rem:H | C.1 |
| C3 | **(N) ∇-transport.** E_k = t^{-C(k,2)} 𝒩 e_k 𝒩^{-1}, 𝒩 = diag(t^{n(ν)}s^{n(ν')}) on P_ν(x;s,1/t). Sub-results: k=1 commutator; boundaries k=m, t=1, s=1; exact identities O_k and the Gaussian kernel; nonsymmetric Gaussian (NS) γ̂X_iγ̂^{-1} = Y•_i. ⚑ Proved modulo Cherednik imports: C2 (Bernstein) and C3-nonsymmetric are still unlocated. Child N-induction-m-iota-duality is `sketched` (role attempt). | DS::N-star-is-nabla-transport, DS::N1-k1-commutator, DS::N-boundaries-km-t1-s1, DS::N-exact-identities-Ok-and-gaussian-kernel, DS::NS-nonsym-gaussian-X-to-Ybullet | proved | 216b (28832) | — (FPSAC is deliberately (N)-free) | C.2 |
| C4 | **Theorem H′ (t→∞).** lim t^{n(μ')}Ψ_s(b_μ) = P_{μ'}(x;s,0) = ωQ'_μ(x;s). SCOOPED in substance; see §3. | DS::theorem-H-prime-t-infinity-q-whittaker | proved | 216b | — | C.3 (cite) |
| C5 | **H′ operator form (b).** The b-basis matrix of E_k has deg_t ≤ n(ν')−n(μ')−C(k,2), and its top coefficient is the HL Q-Pieri ψ_{ν/μ}(s). | DS::H-prime-operator-form-b-Q-pieri | proved | 217e (17697) | — | C.3 |
| C6 | **Square Theorem.** ⋆ on every boundary line of the (s,t)-square (s=0 HL; t=∞ q-Whittaker; s=∞ HL P_{μ'}(x;t); t=0 ωH̃_μ; s=1; t=1/s; t=1). | DS::square-theorem-all-edges-day217e | proved | 217e | — | C.4 |
| C7 | **Lemma R (reflection).** Ψ_{1/s,1/t} = Ψ_{s,t}^{-1}; R-duality Ψ_{s,t}(b_μ) = s^{n(μ)}t^{-n(μ')} e⋆_μ(1/s,1/t). | DS::lemma-R-reflection-psi-inverse | proved | 217e | — | C.4 |
| C8 | **Theorem A (s→∞ edge).** lim s^{-n(μ)} e⋆_μ = P_{μ'}(x;t). | DS::theorem-A-s-infinity-edge | proved | 217e | — | C.4 |
| C9 | **"Thm B" of 217e (t=0 edge).** e⋆_μ\|_{t=0} = s^{n(μ)}P_{μ'}(x;1/s,0) = ωH̃_μ(x;s). SCOOPED (DFK); cite. | DS::theorem-B-t0-edge, DS::thmC-premise-217e-thmB-t0-edge | proved | 217e | rem:t0input (stated as NOT ours) | C.4 (cite) |
| C10 | **Prop C (collapse on the opposite edges).** [P_{1^n}(x;q,T)]e_μ = [n;μ]_T ∏(q;T)_{μ_i}/(q;T)_n. | DS::prop-C-collapse-edges-column-coefficient | proved | 217e | — | C.4 |
| C11 | **t=1/s specialisation.** s_λ⋆s_μ = Σ c^ν_{λμ} s^{c(ν)−c(λ)−c(μ)} s_ν (content-twisted LR). Novelty weak. | DS::N-specialization-t-1-over-s-ribbon-monodromy | proved | ⚑ memory/connections/2026-10-01-star-at-t-equals-1-over-s-is-ribbon-monodromy.md (4080): a connections note, not a proofs/ file | — | C.4 remark |
| C12 | **Lemma ER (edge regularity).** e↔P_ν(x;s,τ) transitions are unitriangular and regular at s=0 and τ=0. ⚑ Erratum: step (4) was corrected 2026-10-05. The computer check never completed (empty logs). | DS::lemma-ER-edge-regularity | proved | proofs/2026-10-03-day219-edge-regularity.md (9989) | — | C.2 lemma |

### (D) s=1 expansion: block valuation law, merge weights, Theorems G and F

| # | Result | Registry :: node id(s) | Grade | Proof file (size) | FPSAC | LV § |
|---|---|---|---|---|---|---|
| D1 | **First order is a biderivation (Thm 1).** B(f,g) = ∂_s(f⋆g)\|_{s=1} = Σ M_kl ∂_kf ∂_lg, a carré du champ. | DS::s1-first-order-star-is-biderivation | proved | 220 (30753) | thm:bider | D.1 |
| D2 | **The matrix M.** M_kr = r e_ke_r + Σ_j L(k−r+j, j) e_{k+j}e_{r−j}. | DS::B-matrix-M_kr-closed-form | proved | 223 (19300); 3-line rederivation in proofs/2026-10-07-day228-two-row-green-and-separator.md §C | prop:M | D.1 |
| D3 | **Block valuation law, lower bound (Thm A of Day 220).** val_{s−1}(c_{λμ}) ≥ ℓ(λ)−κ(λ,μ) for all t. | DS::s1-block-valuation-law-day220 | proved | 220 | thm:blocklaw | D.2 |
| D4 | **Block valuation law, exactness (Thm C).** Equality over Q(t) and at t=0; at t=0 the lead is (−1)^{ℓ−κ}N(λ,μ). Its t=0 premise is now the DFK/Macdonald citation (C9), so Thm C is (N)-free. | DS::thmC-block-law-exact | proved | 220 | thm:blocklaw, rem:t0input | D.2 |
| D5 | **Coarsenings: merge-history expansion (Day 220 "Thm B", not the scooped 217e Thm B).** For a coarsening μ, val = ℓ(λ)−ℓ(μ) exactly, and Lead = Σ over tight histories of ∏W. | DS::thmB-coarsening-exact-valuation, DS::G-premise-thmB-history, DS::F-premise-thmB-history | proved | 220 | thm:coarse | D.3 |
| D6 | **Block expansion (Prop 2.3).** E^{(p_1)}_{λ_1}⋯E^{(p_ℓ)}_{λ_ℓ}(1) = Σ_π ∏ g_C. Recheck Day 227. | DS::block-expansion-prop23-day220 | proved | 220; recheck proofs/2026-10-08-day227-jingliu-telescope.md (16226) | used in proofs (no label) | D.3 |
| D7 | **Merge weights (Thm W).** W_k(J) = (−1)^p [n]/[k] ∏_{j∈J}[k]_{t^j}. Cold re-read paid Day 223 §4. Clio: the second proof is NOT independent (shares Lemma 1.4). | DS::thmW-merge-weight-closed-form, DS::G-premise-thmW | proved | 220 §5b; recheck 223 §4 | thm:W | D.3 |
| D8 | **lin_e T_k = principal specialisation (Day 223 Thm 1.5; (KF)-free proof of W).** lin T_k f = (−1)^d [n]/[k] f(1,t,…,t^{k−1}). Probably known; do not headline. | DS::thmW-second-proof-KF-free | proved | 223 | lem:linT | D.3 |
| D9 | **Theorem G (full-merge lead).** Lead_{λ,(n)} = (1−t^n)/∏(1−t^{λ_i}) · Σ_{H connected} ∏(t^{λ_iλ_j}−1). Includes Lemma 3 (tree–graph identity) and the vertex-deletion second proof. Graph half = Dołęga (cite). | DS::conjG-full-merge-lead-connected-graph, DS::G-lemma3-tree-graph-identity, DS::G-by-vertex-deletion-day223 | proved | 221 (8075, restored); 223 | thm:G | D.4 |
| D10 | **Cor G.** (b) t=0: (−1)^{ℓ−1}(ℓ−1)!. (c) t→1: (−1)^{ℓ−1}n^{ℓ−1}. (a) λ=1^n Mallows–Riordan form, stated in the G node. ⚑ (d) sgn Lead = (−1)^{ℓ−1} for t>0 is **computed only** per the registry, yet FPSAC cor:G prints (d) as part of the corollary. thm:W's "Consequently (−1)^m Lead > 0" seems to cover it; the author must reconcile. | DS::conjG-full-merge-lead-connected-graph (fields corollary_t1, corollary_t1_recheck) | (b),(c) **checked-sober**; (a) within proved G; (d) computed | ⚑ work-in-progress/proofs/2026-10-07-day228-corG-t1-recheck.md (**1819 B, < 3 KB**; exists ONLY under work-in-progress/proofs/, missing from proofs/) | cor:G | D.4 |
| D11 | **Theorem F (coarsening leads factor).** Lead_{λμ} = Σ_π ∏_B Lead_{λ_B,(\|λ_B\|)}, with no interleaving factor. | DS::conjF-coarsening-lead-factorizes | proved | 221 §5 | thm:F | D.4 |
| D12 | **Block multiplicativity (Thm 6.1).** [(s−1)^{ℓ−κ}]c_{λμ} = Σ_π Σ_{(ν^C)} ∏ connected block coefficients. | DS::block-multiplicativity-thm61 | proved | 224 (21670) | thm:blockmult | D.4 |

### (E) Box Complement

| # | Result | Registry :: node id(s) | Grade | Proof file (size) | FPSAC | LV § |
|---|---|---|---|---|---|---|
| E1 | **Box Complement (Thm 2.1) + Column Lemma (Cor 2.2) + lead invariance (Cor 2.3).** c_{N^ℓ−λ, N^ℓ−μ} = s^{N·C(ℓ,2)−(ℓ−1)\|λ\|} c_{λμ}; c_{λ+1^ℓ, μ+1^ℓ} = s^{C(ℓ,2)}c_{λμ}. (N)-free. Cite the mechanism (§3). | DS::box-complement-theorem | proved | 224 (21670) | thm:box, cor:column | E |

### (F) Second order: CT formula, two-part Green polynomials, two-row case, all v=2 leads

| # | Result | Registry :: node id(s) | Grade | Proof file (size) | FPSAC | LV § |
|---|---|---|---|---|---|---|
| F1 | **Constant-term adjoint formula (Thm 1.1).** φ_a⟨T_ag,F⟩_t = CT[Z g(z) F(1/z) K]. = Jing's CT expression (cite Jing91 / JingLiu21 (2.20)); independent re-proof Day 226. | DS::ct-adjoint-formula-thm11 | proved | 225 (23199); recheck 226r (14326) | thm:CT | F.1 |
| F2 | **Two-point / t-string formula (Thm 2.5) = closed formula for all two-part Green polynomials X^λ_{(x,y)}.** Includes the shuffle identity (Lemma 2.4). Novel-as-checked; residual risk Morris 1977 (unread). | DS::two-point-string-formula-thm25, DS::shuffle-identity-lemma24 | proved | 225; recheck 226r §§2–4 | thm:2pt + "Relation to Green polynomials" remark | F.2 |
| F3 | **Two-row specialisation.** X^{(λ1,λ2)}_{(x,y)} in three cases; elementary second proof via III (7.6′) + Young's rule (X = π_k + (t−1)Σ t^{k−1−j}π_j). Classical; no novelty claim. | DS::two-row-specialisation-thm25 (field second_proof_day229) | proved | proofs/2026-10-07-day228-two-row-green-and-separator.md (6692); ⚑ 2nd proof in proofs/2026-10-08-day229-fpsac-round3-and-two-row-young.md (6173), present ONLY in local proofs/ (not in work-in-progress/proofs/); count reconcile proofs/2026-10-08-wake230-two-row-count-reconcile.md (3633, local only) | ex:tworow | F.3 |
| F4 | **Nonsymmetric one-point lemma (Lemma 3.1).** CT[h Z z_1^{−n} K] = (−1)^{a−1}φ_{a−1}h(π_a); a third proof of D8. | DS::nonsym-one-point-lemma31 | proved | 225 | — | F.1 remark |
| F5 | **κ=1, ℓ=3 reduction (Prop 5.1) + closed lead for λ ∋ 1 (Thm 5.2).** | DS::ell3-kappa1-reduction-and-lambda-contains-1 | proved | 224 | prop:reduce (Prop 5.1 only) | F.4 |
| F6 | **x_1-top/subleading recursion (Thm 4.1) + closed c_{λ,(n−1,1)} for ℓ(λ)=3, all s (Thm 4.2 of Day 224).** | DS::subleading-x1-recursion-and-n1-closed-form | proved | 224 | — | F.4 |
| F7 | **All ℓ=3, κ=1 two-part leads (Day 225 Thm 4.2), closing "class 4"; every v=2 lead closed (Cor 4.3).** Includes the closed off-diagonal Γ_a. Cold recheck Day 226 (kill test (3,3,3)→(7,2) by hand). | DS::gprime-v2-class4-open (historical id), DS::gamma-a-offdiagonal-closed-form | proved | 225; recheck 226r §§5–6 | thm:v2, cor:v2, ex:v2 | F.4 |

### (G) Other proved / checked-sober material

| # | Result | Registry :: node id(s) | Grade | Proof file (size) | FPSAC | LV § |
|---|---|---|---|---|---|---|
| G1 | **p_2(Y)-Pieri lemma.** p_2(Y)•e_r = q^{-3}e_{(r,1,1)} − … + τ_r e_{(r+2)}, with τ_r closed form (Baxter-2 / factored). | DS::p2Y-pieri-lemma, DS::tau-r-closed-form-baxter-2 | proved | 206b (13029) | — | G / appendix |
| G2 | **Newton/R7 identities.** C = p_2(Y)•e_2 + 2t(e_2⋆e_2); p_2(Y)•e_r = e_1⋆(e_1⋆e_r) − 2t(e_2⋆e_r). | DS::newton-decomposition-analytic, DS::R7-newton-cancellation-k2 | proved | proofs/2026-09-17-day198-DS-211-via-p2-pieri.md (15727); 206b | — | G / appendix |
| G3 | **Sub-Lemma Z.** Z_r = e_1⋆e_{(r,1)} four-term formula via (L1)–(L4). | DS::sub-lemma-Z-reduction-to-HL-partial-symmetrizers, DS::ds-r11-sub-lemma-Z-operator-form, DS::sub-lemma-Z-e1-star-e1-star-er | proved | 205Z (7031); 205L | — | G / appendix |
| G4 | **DS as Macdonald-operator triangularity (interpretive).** The Day 214 degree count is a Macdonald VI (3.6)–(3.10)-type triangularity. | DS::ds-macdonald-triangularity-connection | **checked-sober** | ⚑ file = null (only the recheck field; nothing to typeset from) | — | B.4 remark |

**Row counts:** A 7, B 12, C 12, D 12, E 1, F 7, G 4; 55 rows in total. B11 is marked "drop"; without it there are 54 substantive rows.

Boundary rule: the validator-style check finds no proved node with a below-boundary *premise* child. All sub-checked-sober children of proved nodes have role `attempt` or `peer-claimed` (e.g. N-induction-m-iota-duality, psi-s-not-macdonald, conjGprime-general-mu-lead, Fn-2phi1-form, k3/k4-gf-telescoping).

---

## 2. Below the boundary but thematically central (candidates for a Conjectures section)

- **DS::conjGprime-general-mu-lead** (in-progress): closed leads for v ≥ 3 (ℓ ≥ 4, κ = 1); three-string shuffle (FPSAC open problem 1).
- **DS::conjG… field corollary (d)** (computed): sgn Lead_{λ,(n)} = (−1)^{ℓ−1} for t>0, t≠1. See D10.
- **E3::hikita-star-min-a-b-plus-1-terms-metaconjecture / E4::…-scorecard** (computed): exactly min(a,b)+1 nonzero terms. The support half is proved; exactness (no coefficient vanishes) is not.
- **E2::ea-star-eb-pieri-general** (in-progress): general e_a⋆e_λ beyond the ℓ-column rule.
- **TC::two-column-straightening-free, caveat field** (computed): no bounded-length e-expansion; tail T(q)+T(N−q)=0.
- **TC::pairwise-kernel-wick-bilinear** (computed) and **TC::Kij-is-jing-half-vertex-exchange** (hunch): Wick/free-field and vertex-operator reading of ★ℓ.
- **DS::p3Y-pieri-conjecture** + children (computed), **DS::pk-Y-pieri-meta-conjecture**, **pk-Y-newton-decomposition** (sketched), **tau-r-k3-closed-form**, **tau-r-k4-closed-form** (computed): p_k(Y)-Pieri family.
- **E3::hikita-star-P_l-meta-shape-conjecture / E4::…-updated** (computed): shape of c_{a−l}(r).
- **E3::hikita-star-q-infty-q-Gaussian-limit** (computed): lim_{q→∞}c_0 = [a+r choose a]_t. Probably Hikita Thm C(ii); cite, don't conjecture.
- **E3::hikita-star-c_0-top-coefficient-closed-form** (computed): now a specialisation of B1, so it could be promoted by citation. The registry grade is still computed.
- **EK::Fn-2phi1-form** (computed): F_n as a ₂φ₁.
- **DS::psi-s-not-macdonald** (computed): Ψ_s(b_μ) is not a Macdonald-type orthogonal basis.
- **DS::e1-power-star-residual-positive** (computed): residual positivity for λ = 1^n.
- **DS::N-induction-m-iota-duality** (sketched): an elementary induction route to (N).
- **FPSAC open problems 2–3** (computed, no registry node): zero set at t = −2; Lead_{λ,(n−1,1)}(1) = 2n²−6n+3.
- **Separator J** (computed, Day 228; in a FPSAC remark): J_⋆ = HT; Dołęga ⊕ ≡ 3/2.
- Stale roots: E2::root, E3::root, E4::root are still `in-progress`, although their content is closed by B1/B3. Registry hygiene only.

---

## 3. Prior art / scooped / folklore: CITE, do not claim

1. **217e "Theorem B" (t=0 edge, C9) is SCOOPED.** Di Francesco–Kedem arXiv:1505.01657v2: (5.15) gives M_{k,1} = E_k\|_{t=0}; level-one Cor 5.8 (5.25); Cor 5.18 (5.27). Add Macdonald VI (5.1), (4.14)(iv). Nodes: theorem-B-t0-edge, thmC-premise-217e-thmB-t0-edge. ⚑ **Name clash:** Day 220's "Thm B" (coarsening merge histories, D5) is a different theorem and is NOT scooped. Rename both in the long version.
2. **H′ (C4) is SCOOPED in substance.** DFK arXiv:1908.00806 Thm KNAN (quoting 1505.01657) and 2112.09798 Thm raiseconj. The (q,t)→(s,t) normalisation match is our inference; check it before citing.
3. **(N) (C3) is folklore-implicit** in DFK arXiv:1704.00154 (Lemma 2.13, Rem 2.14, (4.14), (6.3)–(6.8); C1 = Lemma 2.7). H, H′, Square Thm, Thm A and t=1/s inherit this status. The claim is limited to the explicit edge identifications.
4. **Theorem G, graph half (D9) is prior art.** Dołęga arXiv:1707.02656 Prop 2.1 (+ eq. 10/11) and Lemma 2.3; also Haglund–Tewari 2609.29957 Thm 7.3 (same graph object, "related"), Gessel–Sagan 1996, Mallows–Riordan 1968 (λ=1^n), Brugidou 2509.20959 (hooks), Penrose 1967 (tree–graph identity). Ours: the prefactor, the ⋆-identification, and Thm F/block structure.
5. **d-matrix (C2)** is the classical e→HL transition matrix: Kirillov math/9803006 §3.2 (Thm 3.4), Knuth (0-1 matrices), HKKOTY. Only the identification with the s-lead of ⋆ is ours.
6. **Theorem H (C1)** is novel-as-checked but folklore-adjacent: Orr–Shimozono 1310.0279 Rem 5.6 (Clio, 2026-10-06). Hikita Lemma 6.3 is the unrescaled q→∞ limit. FPSAC makes no novelty claim (rem:H).
7. **Exact valuation corner (B9), partial:** Hikita 2503.23597 Thm C(ii) gives the μ=(n) coefficient d_{λ,(n)}. DS itself is not in Hikita; the triangularity technique is "standard".
8. **Block law at t=0 (D4)** is likely folklore (~75%): DLT SLC 32 (1994) eq. (11) + Jacobi–Trudi. Claim the generic-t lower bound, exactness over Q(t), and W.
9. **Lemma linT / Thm 1.5 (D8)** is probably known: Macdonald III.7 Ex. 2, III (2.2).
10. **Box Complement mechanism (E1):** cite DFK 1704.00154 Rem 3.3 + Hikita Prop 3.8 / Cor 3.10. The ⋆ statement is ours.
11. **CT adjoint formula (F1)** = Jing's constant-term expression (Jing 1991; Jing–Liu 2104.04411 (2.17), (2.20)).
12. **Two-part Green Thm 2.5 (F2):** novel-as-checked against Jing–Liu (their (2.37)–(2.40) restrict the HL index, the transposed slice; their (2.33) does not telescope). Residual risk: Morris LNM 579 (1977), not read first-hand. Morris Math Z 1963 is also relevant (Day 226 dream corrected the volume to 81, DOI 10.1007/bf01111657; wake-226 notes say 80; verify).
13. **Two-row Green (F3)** is classical (Macdonald III (7.6′), III (6.5), Young's rule). Clio's Thm C/D (peer-claimed; WIP 8c148bc) is to be cited as personal communication.
14. **e_1 Pieri (A4)** is Hikita 2503.23597 Thm 3.12; B1 generalises it.
15. **t=1/s (C11):** novelty weak (content-twisted LR, as expected of a ribbon twist).
16. **Textbook nodes (A7, parts of C/D):** Macdonald SFHP I (1.11), (6.6); III (2.2), (2.15), (3.2), (5.7)–(5.8), (7.6′); VI (4.14)(iv), (5.1). Several locators are flagged "from memory" in the registry; verify before printing.

---

## 4. Registry-copy check

`proofs/registry/hikita-*.json` and `work-in-progress/registry/hikita-*.json` are **byte-identical for all 6 files** (cmp), so there are no grade mismatches.

The *proofs* trees, however, differ:
- `2026-10-07-day228-corG-t1-recheck.md` exists only in work-in-progress/proofs/.
- `2026-10-08-day229-fpsac-round3-and-two-row-young.md`, `2026-10-08-wake230-two-row-count-reconcile.md` and `2026-09-30-day211-aha221-mismatch.md` exist only in local proofs/.

---

## 5. Long-version status (PROVE 231, 2026-10-09)

`longversion.tex` (amsart, 18 pp, compiles clean, no undefined refs). Hidden label → registry map (kept here, NOT in the .tex):

| LV label | content | registry node(s) (INVENTORY row) | status in LV |
|---|---|---|---|
| lem:KL, prop:Ak | Key Lemma, parabolic kernel (A_k) | EK::Ak-parabolic-kernel-all-k (A1) | **full proof** (from 207b §2, re-derived) |
| lem:sym1 | partial symmetriser closed form (Day 205b Lemma 1) | EK::lemma-1-2-day205b (A3) | **full proof** |
| prop:Kk | k-set kernel (K_k) | EK::Kk-kset-kernel-all-k (A1) | **full proof** |
| thm:subset, eq:star-reading | subset formula; ⋆-reading | DS::subset-formula-Ak-Kk, EK::hikita-def34-lemma33-interface (A1, A2) | **full proof** (⋆-reading cited: Hikita Def 3.4/Lemma 3.3; bijectivity Lemma 3.1/Cor 3.9) |
| lem:peel | recursion (R) / (Rℓ) generic form | EK::Ek-gf-induction-all-k, TC::two-column-R2-recursion (B1, B5) | **full proof** |
| lem:res | rational Lemma 2′ (residue form of HL kernel) | Day 209 §2 Lemma 2′; EK::lemma-1-2-day205b-2 (A3) | **full proof** |
| lem:Cj | C_j closed form, D_j recursion | part of EK::ek-star-er-pieri-all-k (B1) | **full proof** |
| thm:ellcol | ℓ-column rule (★ℓ) | TC::ell-column-rule, TC::ell-column-closing-identity-Z, TC::ell-column-step-map (B7) | **full proof** (212 skeleton; the ℓ=1 case now SUBSUMES 207b's (L)/(★) proof — one proof for all ℓ) |
| cor:pieri-ekr | e_k⋆e_r Pieri all k | EK::ek-star-er-pieri-all-k, EK::root (B1) | **full proof** (as ℓ=1 of thm:ellcol) |
| ex:k1, ex:k2 | k=1 (Hikita Thm 3.12), k=2 closed F_2 | E2::e2-star-er-pieri-conjecture (B3), A4 | printed + checked |
| cor:pieri-t0 | t=0 truncated geometric | E4::hikita-star-t-zero-specialization (B4) | **full proof** |
| cor:twocol | TC (ℓ=2) + support ≤ ℓ+1 factors | TC::two-column-gf-rule, TC::two-column-straightening-free (B5, B6) | **full proof** |
| lem:step | step map | TC::ell-column-step-map (B7) | **full proof** |
| lem:Z | closing identity (Z) | TC::ell-column-closing-identity-Z (B7) | **full proof** |
| thm:DS, cor:opDS, lem:e-m … lem:valB, lem:Ek1, lem:stable | DS + Op-DS + exact up-set support and s-valuation; E_k(1)=e_k; stability | DS::root, DS::ds-all-lengths-degree-count, DS::full-upset-support-exact-valuation, DS::day214-gauss-valuation-lemma-A (B8, B9) | **full proof** (from 214, re-derived; Lemma B level-set sum now via Mac III (1.4) q-binomial instead of the u-leading-term argument; E_k(1)=e_k via Mac III (2.2),(2.8) instead of Hikita Lemma 3.3) |
| thm:H, cor:dmatrix, lem:admissible, lem:levelset, prop:Hrec, rem:H | Theorem H (s→0 ⋆ ≅ HL e_k-Pieri), explicit matrix, integrality; d-matrix = e→HL transition, ∈ℕ[t], d(1)=#0-1 matrices | DS::theorem-H-s0-star-is-HL-pieri, DS::d-equals-hall-littlewood-transition, DS::d-lambda-mu-t-count-01-matrices (C1, C2) | **full proof** (from 215, re-derived); rem:H = prior-work remark (Kirillov, Hikita L6.3). C5–C10 (square edges, Lemma R, Thm A, Prop C) deliberately EXCLUDED: they rest on (N) (Cherednik imports unlocated); rem:edges says so and uses none of it |
| thm:bider, prop:M, lem:polar, prop:block, thm:blocklaw, rem:t0input, thm:coarse, thm:W, lem:linT, lem:parabolic, lem:linP | biderivation; M matrix; polarization; block expansion; block law (lower bound + exactness via t=0 raising operators); coarsenings; merge weights | D1–D8 (DS::s1-first-order-star-is-biderivation, B-matrix-M_kr-closed-form, s1-block-valuation-law-day220, thmC-block-law-exact, thmB-coarsening-exact-valuation, block-expansion-prop23-day220, thmW-merge-weight-closed-form, thmW-second-proof-KF-free) | **full proof** (from 220 + 223 §1; thm:W via the (KF)-free route of 223; Lead(0)≠0 in thm:coarse now from thm:W at t=0 instead of 220's roots-of-unity argument; t=0 edge CITED (DFK+Macdonald), not proved). prop:M's proof uses cor:pieri (§8), which is still a sketch → prop:M inherits that until §8 is written |
| thm:G, cor:G, thm:F, thm:blockmult, lem:forest, lem:telescope, lem:treegraph | full-merge lead (connected-graph cumulant), Cor G (a)–(d), Theorem F, block multiplicativity | D9–D12 (DS::conjG-full-merge-lead-connected-graph, G-lemma3-tree-graph-identity, conjF-coarsening-lead-factorizes, block-multiplicativity-thm61) | **full proof** (from 221 §§1–5, 224 §6; Cor G(c) via vertex-weighted Cayley/Prüfer, Cor G(b) via tree form) |
| thm:box, lem:inversion, cor:column, thm:lin, cor:pieri | Box Complement, Column Lemma, plethystic lin, Pieri 2nd proof (C-form); F_n(t^m) identity for 1≤n≤m | E1, B2 (DS::box-complement-theorem, DS::plethystic-lin-e-and-pieri-reproof) | **full proof** (from 224 §§2–3). prop:M (§6) now rests on a full proof |
| prop:reduce, thm:CT, thm:2pt, lem:shuffle, lem:pair2, thm:v2, cor:v2, ex:tworow, ex:v2 | Γ_a identity, reduction, CT adjoint formula, two-point formula (string expansion + string weight + shuffle), pairing lemma, every v=2 lead | F1–F7 (DS::ct-adjoint-formula-thm11, two-point-string-formula-thm25, shuffle-identity-lemma24, ell3-kappa1-reduction…, gprime-v2-class4-open) | **full proof** (from 225 §§1–2, 4 and 224 §§1, 5); ex:tworow keeps the Jing–Liu credit + Young's-rule proof |

**"PDF is the claim" checks for §§2–3** (scripts/day231/, re-implemented from the PRINTED statements, independent pointwise subset-formula engine in exact rationals):
- `check_Y_printed.py` (log): printed T_i, π, Y_i conventions ⇒ t^{-C(k,2)}e_k(Y)F = E_kF for m=3,4, all k, F ∈ {1,e_1,e_2,e_1²,e_1e_2}; Hikita Thm 3.12 recovered. ALL True.
- `check_pieri_printed.py` (log): cor:pieri-ekr vs engine k≤4, r≤5; support ≤ min(k,r) after collection; cor:pieri-t0 vs engine at t=0, k≤4, r≤5; FPSAC C_{a,b} form == F form (k≤4, r≤5); F_0,F_1,F_2 closed forms (ex:k2). ALL True.
- `check_ellcol_printed.py` (log): thm:ellcol as an identity Γ_k = T_k at exact random points, ℓ=1,2,3, m≤5, k≤4 (85 cases); cor:twocol explicit V == general V; lem:Z ℓ≤5; lem:Cj; K symmetry. ALL True. Negative control (`neg_control.py`, K_{ij} slots swapped) FAILS as it should.

**Notation clashes to fix when the ported sections are rewritten:** Γ_k (§3 generating function) vs Γ_a (§9 bilinear form); C_I / C_j(x) (§3) vs C_{a,b} (cor:pieri); κ_c (§3 residues) vs κ(λ,μ) (block count).
- `check_DS_printed.py` (log, n≤5, t∈{−5/11, 0}; s-polynomials recovered by 26-point interpolation + 2 held-out points): thm:DS support = exact up-set, val_s = n(μ), c_{λλ}=s^{n(λ)}, c|_{s=1}=0 off-diagonal, t=0 lead = 1; cor:opDS; lem:Ek1; lem:stable (m=n vs n+1); rem d_{(1^4),(211)}=1+3t; Peel Lemma n≤10 (5511 triples); level-set q-binomial identity. ALL True.
- `check_H_printed.py` (log): thm:H (1) integrality + (3) explicit vertical-strip matrix vs engine, all μ,k with |μ|+k≤5 at t∈{−5/11,3/2} (52 cases); cor:dmatrix first formula vs independent HL P (Macdonald symmetrisation) n≤4 at 2 t-values; d∈ℕ[t] and d(1)=#0-1 matrices n≤4 (t-interpolation). ALL True.
- `check_blocklaw_printed.py` (log): prop:M k≤4 vs engine; thm:blocklaw val=ℓ−κ on all 70 non-diagonal up-set pairs n≤5 at t=3/5 and t=0 (κ by brute force from the printed definition), t=0 lead has sign (−1)^{ℓ−κ} and is a nonzero integer; thm:W from the printed iterated-commutator definition (20 cases, n≤6); lem:linT (k≤3, 5 test f); lem:linP (7 λ, HL P by symmetrisation). ALL True.
- `check_coarse_printed.py` (log): thm:coarse history formula (histories enumerated from the printed definition, W closed form) = engine Lead on all coarsening pairs n≤5 at t∈{3/5,−2,0} (96), including the t=−2 zeros; #H≥1 at t=0. ALL True.
- `check_leads_printed.py` (log, n≤5): thm:G at t∈{3/5,−2,7/3}; cor:G (a) Mallows–Riordan form, (b) t=0, (c) t→1 numerically; thm:F at t∈{3/5,−2}; thm:blockmult on all 35 pairs with κ≥1 (n≤5). ALL True. (n=6 run: check_leads_printed_n6.log.)
- `check_box_printed.py` (log): thm:box on 10 (λ,N) incl. N<degree and padded zeros (18 coefficient pairs); cor:column; lem:inversion pointwise (N=3); thm:lin + ring-map form; cor:pieri C-form vs engine k≤4,r≤4; F_n(t^m) identity for n≤m (also holds for m<n numerically — NOT claimed). ALL True.
- `check_v2_printed.py` (log): thm:2pt vs ⟨T_a g, p_xp_y⟩ from a power-sum expansion computed from scratch (a≤3, 6 test g, 30 cases, both (x,y) orders); lem:pair2; lem:shuffle by brute-force shuffle sums A,B≤3 + partial fractions; ex:tworow vs p_xp_y in an independently built HL P-basis, |λ|≤6. ALL True.
- `check_v2lead_printed.py` (log): thm:v2 (printed formula, built only from printed prop:M, thm:2pt, Ξ) vs engine leads for (2,2,2)→(5,1),(3,3), (3,3,1)→(5,2), (3,2,2)→(6,1) in all orderings; ex:v2 printed polynomials (3,3,3)→(7,2) and (4,4,2)→(7,3) = thm:v2 formula at 20 t-points, all orderings. ALL True.

**Status 2026-10-09 end of PROVE 231:** every section is written in full. One visible \wip marker remains: the title footnote (Robin's authorship/affiliation/acks/AI declaration). Remaining polish: two overfull lines (~35–40 pt, §4.3 intro and §6.1), and the notation clashes listed above (Γ_k→𝓔_k and κ_c→η_c FIXED in §3; C_I vs C_{a,b} remains). Not included on purpose: (N), the square edges, DS-from-(N), Lemma ER, t=1/s, G4 (no file).
