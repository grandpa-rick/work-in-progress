# Day 166 — Slice identities for $G := L_A F_P/u_3$; slice of $L_A F_1$; honest stall

**Date:** 2026-09-04 (Day 166). **Status: HONEST STALL WITH ONE NEW PROVED LEMMA.**

Target (from PROVE.md): prove any ONE of the three-way equivalence:
$\Sigma_0$'s closed form, $R^{(-1)}$'s closed form, or Theorem B (Day 162).

**Result 1 (new, PROVED).** The 2-variable operator $L_A = (\theta+u_1+1)(\theta+u_2+1) - \partial_T$
satisfies $L_A F_0 = 0$ (Day 158 Prop A). Extended to the 3-variable $F_P$ (treating $u_3$ as a
parameter), $L_A F_P$ is **exactly divisible by $u_3$**. Define $G := L_A F_P/u_3 \in
\mathbb Q[u_1,u_2,u_3][[T]]$. Then:
$$G|_{u_1=0}(u_2, u_3) \;=\; -(\theta + u_2 + 1)\, F_0(u_2, u_3),\tag{S1}$$
$$G|_{u_2=0}(u_1, u_3) \;=\; -(\theta + u_1 + 1)\, F_0(u_1, u_3),\tag{S2}$$
where $F_0(x, y) := \sum_{k \ge 0} (T^k/k!) A_k(x) A_k(y)$, $A_k(x) = (x+1)_k$ (Day 158).

**Consequence.** $L_A F_1 = G|_{u_3=0}$ (from Taylor expansion of $L_A F_P = u_3 G$ in $u_3$).
Setting $u_2 = 0$ in $G|_{u_3=0}$ vs setting $u_3 = 0$ in (S2):
$$L_A F_1(u_1, 0) = G(u_1, 0, 0) = -(\theta + u_1 + 1)\, F_0(u_1, 0) = -\sum_{k\ge 0} T^k (u_1+1)_{k+1}.\tag{S3}$$

By $u_1 \leftrightarrow u_2$ symmetry:
$$L_A F_1(0, u_2) = -\sum_k T^k (u_2+1)_{k+1}.\tag{S4}$$

**Numerical certification of (S1)-(S4) to $n\le 5$**: script
`scratch/day166/verify_slice_identity.py`. All 4 identities hold exactly.

**PROOF of (S1)-(S4).** Let $L^{(i,j)}_A := (\theta + u_i + 1)(\theta + u_j + 1) - \partial_T$
(so $L_A = L^{(1,2)}_A$; $L^{(2,3)}_A F_0(u_2, u_3) = 0$ by Day 158 Prop A applied to the
$(u_2, u_3)$ variables).

Restricting the operator identity $L^{(1,2)}_A F_P = u_3 G$ to $u_1 = 0$:
- $L^{(1,2)}_A|_{u_1=0} = (\theta + 1)(\theta + u_2 + 1) - \partial_T$
- $F_P|_{u_1 = 0}(u_2, u_3) = F_P(0, u_2, u_3) = F_0(u_2, u_3)$ by $S_3$-symmetry of $F_P$.
- Hence $L^{(1,2)}_A F_P|_{u_1=0} = [L^{(1,2)}_A|_{u_1=0}] F_0(u_2, u_3)$.

Now $L^{(1,2)}_A|_{u_1=0} - L^{(2,3)}_A = [(\theta + 1)(\theta+u_2+1)] - [(\theta+u_2+1)(\theta+u_3+1)]
= (\theta+u_2+1)[(\theta+1) - (\theta+u_3+1)] = -u_3(\theta+u_2+1)$.

Since $L^{(2,3)}_A F_0(u_2, u_3) = 0$ (Day 158 Prop A), we get:
$$L^{(1,2)}_A|_{u_1=0} F_0(u_2, u_3) = -u_3(\theta+u_2+1) F_0(u_2, u_3).$$

Combined with $L_A F_P = u_3 G$: $u_3 \cdot G|_{u_1=0} = -u_3(\theta+u_2+1) F_0(u_2, u_3)$, giving
(S1) after cancellation. (S2)-(S4) follow by symmetry and by evaluation at further slices. $\square$

**Result 2 (NEW STALL POINT).** From (S3), we know $L_A F_1$ on the $u_2 = 0$ codimension-1
slice. By $u_1 \leftrightarrow u_2$ symmetry, we also know $L_A F_1$ on $u_1 = 0$. But $L_A F_1$
is a 2-variable symmetric polynomial (in $u_1, u_2$), and its two coordinate slices are
NOT enough to determine it uniquely.

Specifically, write $L_A F_1(u_1, u_2) = -\sum_k T^k Q_k(u_1, u_2)$ with $Q_k$ symmetric. Then
$Q_k(x, 0) = Q_k(0, x) = (x+1)_{k+1}$ (from (S3)/(S4)). Numerical data ($n \le 7$):

| $k$ | $Q_k$ in $(E_1, E_2)$-basis |
|---|---|
| 0 | $E_1 + 1$ |
| 1 | $E_1^2 + E_1 E_2 + 3 E_1 + 5 E_2 + 2$ |
| 2 | $E_1^3 + \tfrac{3}{2} E_1^2 E_2 + \tfrac{1}{2} E_1 E_2^2 + 6 E_1^2 + 15 E_1 E_2 + \tfrac{9}{2} E_2^2 + 11 E_1 + \tfrac{45}{2} E_2 + 6$ |
| 3 | $E_1^4 + \tfrac{11}{6} E_1^3 E_2 + E_1^2 E_2^2 + \tfrac{1}{6} E_1 E_2^3 + \ldots$ (rational coefs) |
| ... | (grows in complexity, non-integer rationals throughout) |

No clean pattern emerges. In particular:
- $Q_k \neq (E_1+k+1) Q_{k-1}$ (checked, differences nonzero).
- $Q_k \neq A_{k+1}(u_1) A_{k+1}(u_2)/(k+1)!$ (fails at $k=0$).
- $Q_k \neq A_{k+1}(u_1) + A_{k+1}(u_2) - (k+1)!$ (fails at $k=1$).
- Search over degree-bounded ansatze $Q_k = P(E_1, E_2, T)$ with $P$ of $E$-degree $\le n$
  fails (agent report Day 166, Task 1-2).

**Precise stall.** To prove Theorem B (equivalently $\Sigma_0$'s closed form, equivalently
$R^{(-1)}$'s closed form via the three-way equivalence), we need to determine either:
(a) $Q_k(u_1, u_2)$ in closed form, or
(b) equivalently, $\ell^{\rm top}_{-1}(L_A F_1/F_0) = \Sigma_0$ in closed form.

The natural continuation is to compute $G$ FULLY (as a 3-variable polynomial): from (S1)-(S2)
and additional structure ($u_1 \leftrightarrow u_2$ symmetry of $G$, plus perhaps constraints
from the $S_3$-orbit $L^{(1,3)}_A F_P = u_2 G^{(2)}$, etc.), it may be possible to reconstruct
$G$ and hence $G|_{u_3=0} = L_A F_1$.

## Recommended queue for Day 167+

1. **Reconstruction of $G(u_1, u_2, u_3)$.** Use (S1), (S2), and the $S_3$-orbit identities to
   build up $G$ as a polynomial in $u_1, u_2, u_3$. Concretely: expand $G$ in the basis
   $\{u_1^a u_2^b u_3^c : a, b, c \ge 0\}$, use (S1) to determine $[u_1^0]$-coefficients,
   (S2) to determine $[u_2^0]$-coefficients. Then use additional structure (see step 2).

2. **Additional structural constraints.** The 3-variable $L_A F_P$ is *not* the sum of
   $S_3$-orbit $L^{(i,j)}_A F_P$'s (since $L^{(i,j)}_A F_P = u_k G^{(k)}$ with different $u_k$
   factors). But by $S_3$-symmetry of $F_P$: $L^{(i,j)}_A F_P/u_k$ = permuted-version of $G$.
   Combining three orbits gives 6 slice identities (2 per orbit); check if together they
   uniquely determine $G$.

3. **Alternative: Route (i) prefactor** (Day 163 §3 stall). Compute
   $\partial_{u_3}[\varrho \log(\ell^{\rm top}_3(\mathcal R)/V(u))]|_{u_3=0}$ in closed form.
   Since $\ell^{\rm top}_3(\mathcal R) = V(\nu)$ (Day 148 §5, audited Day 149 §7), this reduces
   to computing $\partial_{u_3} \log[V(\nu)/V(u)]|_{u_3=0}$. Using the $\nu$-system:
   $\nu_i(q + T\nu_i) = u_i$, and the derivatives $\dot\nu_3 = 1/q$,
   $\dot\nu_i = -\nu_i \dot q/R_i$ for $i=1,2$ (Day 166 derivation, §3 below).

4. **Prove $\Sigma_0$ P-recurrence directly.** Day 165 found (via nullspace):
   $(k+3) P_{k+3} = (12k+26)P_{k+2} - (48k+64)P_{k+1} + (64k+32) P_k$
   for the $A_P(W)$ coefficient sequence. If proved structurally (from $\Sigma_0$'s
   definition), the closed form follows by verifying initial conditions.

## §3. Derivation of $\dot\nu_i := \partial_{u_3}\nu_i|_{u_3=0}$

From $\nu_i(q + T \nu_i) = u_i$, differentiate w.r.t. $u_3$:
$\dot\nu_i(q + T\nu_i) + \nu_i(\dot q + T \dot\nu_i) = \delta_{i,3}$

With $\dot q = -T \partial_{u_3} e_1(\nu) = -T(\dot\nu_1 + \dot\nu_2 + \dot\nu_3)$.

**For $i = 3$:** at $u_3 = 0$, $\nu_3 = 0$; equation becomes $\dot\nu_3 \cdot q = 1$, so
$\dot\nu_3 = 1/q$.

**For $i = 1, 2$:** using $u_i/\nu_i = q + T\nu_i$ (from constraint):
$\dot\nu_i R_i = -\nu_i \dot q$ where $R_i := q + 2T\nu_i$;
so $\dot\nu_i = -\nu_i \dot q/R_i$.

**Solving for $\dot q$:** $\dot q = -T[-\nu_1\dot q/R_1 - \nu_2\dot q/R_2 + 1/q]$, i.e.,
$\dot q(1 - T\nu_1/R_1 \cdot (-1)/(-1) - \ldots)$... let me redo:
$\dot q + T\dot q(\nu_1/R_1 + \nu_2/R_2) = -T/q$

Using $T\nu_i/R_i = \tfrac{1}{2} - q/(2R_i)$ (Day 152 (5.1)) and $R_1 + R_2 = 2$,
$R_1 R_2 = 1 - T^2(E_1^2 - 4E_2) = 2 - q^2 - 2E_1T$ (Day 161):

$T(\nu_1/R_1 + \nu_2/R_2) = 1 - q(1/R_1+1/R_2)/2 = 1 - q/R_1R_2$

$\dot q [2 R_1 R_2 - q]/R_1 R_2 = -T/q$

$$\dot q = -\frac{T R_1 R_2}{q(2 R_1 R_2 - q)}. \tag{S5}$$

**Consequence.** From (S5) and $\dot\nu_i$ formulas, we have full first-order in $u_3$
expansion of the $\nu$-system. This is the input to Route (i) continuation.

## §4. What (S3)/(S4) actually give us

Formula (S3): $L_A F_1(x, 0) = -\sum_k T^k (x+1)_{k+1}$. Sum in closed form:

$\sum_k T^k (x+1)_{k+1} = (x+1) \sum_k T^k (x+2)_k$

The inner sum is a divergent (formal) series; recognizing as ${}_1F_0$: $\sum_k (x+2)_k T^k/k!
= (1-T)^{-(x+2)}$, so $\sum_k (x+2)_k T^k = ?$ — this is a Borel-type object, not a rational
function. In particular, $L_A F_1(x, 0)$ has coefficients $-(x+1)_{k+1}$ which grow like
$(k+1)!$; not rational-in-$T$.

Hence $L_A F_1(x, 0)$ is a formal series but not a nice algebraic function. This is
consistent with our failure to find a compact algebraic closed form for $Q_k(u_1, u_2)$.

**Structural obstacle.** $L_A F_1$'s slices grow factorially in $T$; any closed-form
representation must be as a formal series, not a rational function of $Y, q$ etc. This
suggests that reconstructing $L_A F_1$ requires working with $F_0$ (which also grows
factorially) rather than with the algebraic $Y, q$ variables directly.

**Alternative strategy.** Rather than compute $L_A F_1$ in closed form, work with the
COMBINATION $L_A F_1 / F_0$, which has wt-decomposition into $q'$ (top, PROVED),
$\Sigma_0$ (sub-top), etc. The wt-decomposition converts the factorial growth into
algebraic-in-$(Y, q)$ pieces. This is Day 164's original route; the sub-top layer is
still open.

## §5. Registry updates

- `LA-FP-slice-identity-u1-eq-0` (NEW): **`proved`** (this file §PROOF).
  Role: `premise` toward reconstructing $L_A F_1$.
- `LA-F1-slice-u2-eq-0` (NEW): **`proved`** (via §PROOF + $S_2$ symmetry).
  Formula: $L_A F_1(x, 0) = -\sum_k T^k (x+1)_{k+1}$.
- `sigma-0-closed-form` (Day 165): **stays `checked-sober` at N=24, 15 specs.**
- `bar-D-closed-form-E3-zero` (Day 162 Theorem B): **stays `checked-sober` at $n \le 14$.**
- `R-minus-one-closed-form` (Day 162): **stays `checked-sober` at $n \le 14$.**
- `narayana-layer-d1-E3-zero` (C.5): **stays `computed`.**

**No upgrade of Missing Lemma (R) this session.**

## §6. Scripts

Location: `/home/agent/projects/scratch/day166/`.

| Script | Purpose | Verdict |
|---|---|---|
| `verify_slice_identity.py` | Verify (S1), (S2) numerically | ✓ $n \le 5$ (n=6 needs FP at N=7) |
| `explore_Q_n.py` | Search for closed form of $Q_k$ | Failed: ansatze don't fit |

## §7. Discipline scorecard

- **[[feedback_check_convention_before_compute]]**: Applied. Used `FP_coeffs` from
  `scratch/day152/lib.py` (true object). Verified $F_0$ closed form matches library.
- **[[feedback_verify_scripts_implement_what_they_claim]]**: Applied. Verified $L_A F_P$
  is $u_3$-divisible (up to boundary $n=N$).
- **Rule 11 firing #8**: Unfolding the operator identity $L_A F_P = u_3 G$ into slice
  identities (S1)-(S4). Same style as Day 164's Riccati split, one level higher.
- **Honest stall**: The two slice identities are not enough to determine $L_A F_1$
  uniquely. Reported precisely.
- **Feedback for future**: When facing a 3-variable object with $S_3$-symmetry, extract
  all 6 slice identities (3 orbits × 2 slices each) before attempting reconstruction.

## §8. One-line handoff

Day 166 proved a new operator identity $L_A F_P = u_3 G$ with clean slice formulas
$G|_{u_1=0} = -(\theta+u_2+1) F_0(u_2, u_3)$ and $G|_{u_2=0} = -(\theta+u_1+1) F_0(u_1, u_3)$;
extracted $L_A F_1$'s single-variable slices in closed form; failed to reconstruct
$L_A F_1$'s full 2-variable form from these slices; Missing Lemma (R) unchanged.
