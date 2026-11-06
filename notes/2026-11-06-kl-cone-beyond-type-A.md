# The KL cone beyond type A (Wake 279, 2026-11-06)

**Status: computed / sketched — not checked.** Nothing here has had a cold re-derivation.

## Question (Q68)
Goertzen–Williamson, *Kazhdan–Lusztig Basis and Optimization*, arXiv:2604.18894, show that for S_n Specht modules
of shape hook (Thm 4.11), two-column (Thm 4.14) and (n−2,2) (Thm 4.17) the KL basis at v=1 spans the maximal
{1+s}-invariant cone (computed n ≤ 8; fails non-simplicially at (4,4,1), Ex 4.21). Does this persist for left-cell
modules of other finite Coxeter groups?

**Convention.** GW's W-graph convention (Ex 4.7): (1+s)e_x = 0 if s∈D(x), else 2e_x + Σ_{y: s∈D(y)} μ(x,y) e_y.
In the C′ convention ((1+s)c_w = 2c_w for sw<w) maximality already fails in A3 — that is the dual (w0-twisted,
minimal-cone) question. Everything below uses GW's convention. *(computed)*

## Counterexample in B3 *(computed, sympy exact; hand-checked once)*
Irreducible 3-dim left cell {x=s3s1, y=s2s3s1, z=s1s2s3s1}: D(x)={1,3}, D(y)={2}, D(z)={1};
W-graph edges x–y, y–z; μ(x,z)=0. Put v = c_z − c_x. Then
(1+s1)v = 0, (1+s2)v = 2v, (1+s3)v = 2c_z,
so cone{c_x, c_y, v} is invariant, simplicial and strictly contains C_KL. The module is irreducible, so the maximal
cone is unique (Mejstrik–Protasov, via GW Thm 4.5) and strictly contains C_KL. **So the answer is no in type B3.**

## Dominated-twin lemma *(sketched)*
If x≠z in a left cell Γ with D(z) ⊆ D(x) and μ(y,z) ≥ μ(y,x) for every y∈Γ with D(y) ⊄ D(x), then v = c_z − c_x
satisfies (1+s)v ∈ cone(C_KL ∪ {v}) for every s, and cone{c_w (w≠z), v} is simplicial, invariant, and strictly
contains C_KL.
*Sketch.* s∈D(z) ⇒ s∈D(x): both are killed. s∈D(x)∖D(z): (1+s)v = (1+s)c_z ∈ C_KL. s∉D(x):
(1+s)v = 2v + Σ_{y: s∈D(y)} (μ(y,z)−μ(y,x)) c_y, and every such y has D(y)⊄D(x), so the coefficients are ≥ 0.
Replacing c_z by c_z−c_x is a unitriangular basis change, so the cone stays simplicial.

## Scan *(computed; test (A) is depth-limited, so the failures from test (B) and the twin lemma are what count)*
| type | result |
|---|---|
| A3 | all cells pass (agrees with GW) |
| D4 | irreducible cells (dim 3, 4) pass; reducible (dim 10, 14) fail |
| I2(m), m=4,5,6,8 | middle (reducible) cells fail |
| B3 | all non-degenerate cells fail, incl. 6 irreducible 3-dim |
| H3 | all fail, incl. irreducible 5-dim |

Dominated twins appear in all B3 3- and 4-dim cells, every dihedral middle cell and 15/22 H3 cells, and never in
A3, A4, D4.

## Remarks / open
- **Hunch:** simply-laced groups have no dominated twins (twins seem to come from m≥4 strings x, s_j x, s_i s_j x),
  which leaves the simply-laced case of Q68 open.
- Reducible cells (I2(4) middle, B3 4-dim) have two distinct maximal pointed invariant cones (±v on a 1-dim
  submodule), so "the maximal cone" is ill-posed there. *(computed)*
- H3 cells of type 5_2 fail with no dominated twin, so there is a second obstruction we don't understand yet.
- Owner check owed: Mejstrik–Protasov arXiv:2505.01229 not yet read.
