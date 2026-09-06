# Day 172 — E₂-shift: reduction to a single missing lemma

**Date:** 2026-09-06. **Status:** partial. The $E_2$-shift conjecture is **reduced** to one sub-claim (called **(A)** below): that $\mathrm{tops}^{(n)}[b]$ lives in $\mathbb Q[E_1,E_2,E_3]$ (i.e. is free of $E_4,\dots,E_n$). Assuming **(A)**, the shift-law is a two-line consequence of a **stability identity** for factorial Schur functions that is proved rigorously here (via the classical $\det[(u_i)_{\lambda_j}]$ formula and a Vandermonde-cancellation trick). **(A)** is verified computationally for all $(n,b)\in\{3,\dots,7\}\times\{0,\dots,6\}$ (35 cases), plus the sub-top structure of the full $\Psi_b^{(n)}$ shows $E_4,\dots,E_n$ *do* appear at sub-top ρ — so (A) is a genuine cancellation, not a trivial degree bound.

---

## 1. Problem statement

Let $\Psi^{(n)}(f) := \mathcal T(fV_n)/V_n$ where $\mathcal T$ is falling-factorial substitution ($u_i^{a} \mapsto (u_i)_a = u_i(u_i-1)\cdots(u_i-a+1)$) and $V_n=\prod_{i<j}(u_i-u_j)$. For $b\ge0$, set $\Psi_b^{(n)} := \Psi^{(n)}(e_2^b) \in \mathbb Z[E_1,\dots,E_n]$.

Define the ρ-weight on $\mathbb Z[E_1,\dots,E_n]$ by $\rho(E_k) := \lceil k/2 \rceil$. The **top slice** is
$$\mathrm{tops}^{(n)}[b] \;:=\; [\rho = b]\,\Psi_b^{(n)}\;.$$

**Theorem (E₂-shift, target).** For all $n\ge 3$ and $b\ge 0$,
$$\boxed{\;\mathrm{tops}^{(n)}[b] \;=\; \mathrm{tops}^{(3)}[b]\Big|_{E_2 \,\mapsto\, E_2 \,-\, c_n E_1}\;} \qquad c_n = \binom{n-1}{2}-1\;.$$

---

## 2. The reduction — architecture

Two ingredients:

**Ingredient 1 (Stability, PROVED — §3):** For all $n\ge 2$, $b\ge 0$, viewing $\Psi_b^{(n+1)}$ as a polynomial in $E_1^{(n+1)},\dots,E_{n+1}^{(n+1)}$,
$$\Psi_b^{(n+1)}\big|_{E_{n+1}=0} \;=\; \tau_{n}^{-1}\bigl(\Psi_b^{(n)}\bigr)$$
where $\tau_n^{-1}$ is the ring automorphism of $\mathbb Q[E_1,\dots,E_n]$ induced by $u_i\mapsto u_i-1$ in $n$ variables:
$$\tau_n^{-1}(E_k) \;=\; \sum_{j=0}^{k}\binom{n-j}{k-j}(-1)^{k-j}\,E_j\qquad(E_0:=1)\;.$$

**Ingredient 2 (Sub-claim A, VERIFIED but not proved — §5):** $\mathrm{tops}^{(n)}[b] \in \mathbb Q[E_1,E_2,E_3]$ for all $n\ge 3$, $b\ge 0$.

**Deduction (§4):** Take top-ρ of both sides of Ingredient 1. Using (A), the LHS equals $\mathrm{tops}^{(n+1)}[b]$ (setting $E_{n+1}=0$ removes nothing). The RHS becomes $\sigma_n(\mathrm{tops}^{(n)}[b])$ where $\sigma_n$ is the "top-ρ symbol" of $\tau_n^{-1}$; this symbol acts on $E_1,E_2,E_3$ as
$$\sigma_n:\ E_1\mapsto E_1,\quad E_2\mapsto E_2-(n-1)E_1,\quad E_3\mapsto E_3\;.$$
So $\mathrm{tops}^{(n+1)}[b] = \mathrm{tops}^{(n)}[b]\big|_{E_2\to E_2-(n-1)E_1}$. Iterating from $n=3$ to $n=N$ yields the total shift $c_N = \sum_{n=3}^{N-1}(n-1) = \sum_{k=2}^{N-2}k = \binom{N-1}{2}-1$. $\square$

The rest of this note fleshes out §3 (rigorous), §4 (rigorous once (A) is granted), and §5 (evidence for (A) and remarks on why it is subtle).

---

## 3. Stability identity — proof

### 3.1 Setup

Recall (Day 149 Thm A) $\Psi^{(n)}(s_\mu) = \mathfrak s_\mu^{(n)}$, where the **factorial Schur function** in $n$ variables (falling convention) is
$$\mathfrak s_\mu^{(n)}(u_1,\dots,u_n) \;=\; \frac{\det\bigl[(u_i)_{\lambda_j}\bigr]_{i,j=1}^n}{V_n(u)}\;,\qquad \lambda_j = \mu_j + n - j\;.$$
The Schur→factorial-Schur map is linear, and $e_2^b = \sum_{\mu\vdash 2b} K_{\mu'(2^b)} s_\mu$ (Kostka expansion), so
$$\Psi_b^{(n)} \;=\; \sum_{\mu\vdash 2b,\ \ell(\mu)\le n} K_{\mu'(2^b)}\,\mathfrak s_\mu^{(n)}\;.$$

### 3.2 Factorial-Schur stability lemma

**Lemma (Factorial-Schur stability).** For any partition $\mu$ with $\ell(\mu)\le n$,
$$\mathfrak s_\mu^{(n+1)}(u_1,\dots,u_n,0) \;=\; \mathfrak s_\mu^{(n)}(u_1-1,\dots,u_n-1)\;.$$
For $\ell(\mu) = n+1$: $\mathfrak s_\mu^{(n+1)}(u_1,\dots,u_n,0) = 0$.

*Proof.* Write $\lambda_j^{(n+1)} = \mu_j + (n+1) - j$ (using $\mu_{n+1}$ implicit if needed).

**Case $\ell(\mu)=n+1$**: $\mu_{n+1}\ge 1$, so $\lambda_{n+1} = \mu_{n+1} \ge 1$. Setting $u_{n+1}=0$: row $n+1$ of $\det[(u_i)_{\lambda_j}]$ has entries $(0)_{\lambda_j}$. For $j\le n$, $\lambda_j = \mu_j + (n+1) - j \ge (n+1) - n = 1$, so $(0)_{\lambda_j} = 0$. For $j=n+1$, $\lambda_{n+1}\ge 1$ likewise. So the entire row is zero; det = 0. $V_{n+1}(u_1,\dots,u_n,0) \ne 0$ generically (as a polynomial), so $\mathfrak s_\mu^{(n+1)}|_{u_{n+1}=0} = 0$.

**Case $\ell(\mu)\le n$** (so $\mu_{n+1}=0$, giving $\lambda_{n+1}=0$): Row $n+1$ of the matrix at $u_{n+1}=0$ is $(\,(0)_{\lambda_1},\dots,(0)_{\lambda_n},(0)_0\,) = (0,\dots,0,1)$ since $\lambda_j\ge 1$ for $j\le n$ and $(0)_0=1$. Cofactor expansion along row $n+1$:
$$\det\bigl[(u_i)_{\lambda_j}\bigr]_{i,j=1}^{n+1}\Big|_{u_{n+1}=0} \;=\; (-1)^{(n+1)+(n+1)}\cdot 1\cdot \det\bigl[(u_i)_{\lambda_j^{(n+1)}}\bigr]_{i,j=1}^n \;=\; \det\bigl[(u_i)_{\lambda_j^{(n+1)}}\bigr]_{i,j=1}^n\;.$$

Now use $\lambda_j^{(n+1)} = \lambda_j^{(n)} + 1$ (since $\lambda_j^{(n+1)} - \lambda_j^{(n)} = 1$ for all $j$), and the falling-factorial identity $(u)_{k+1} = u\cdot (u-1)_k$:
$$\det\bigl[(u_i)_{\lambda_j^{(n+1)}}\bigr]_{i,j=1}^n \;=\; \det\bigl[u_i\cdot (u_i-1)_{\lambda_j^{(n)}}\bigr]_{i,j=1}^n \;=\; \Bigl(\prod_i u_i\Bigr)\cdot \det\bigl[(u_i-1)_{\lambda_j^{(n)}}\bigr]_{i,j=1}^n$$
by pulling the row factor $u_i$ out of row $i$. So
$$\det\bigl[(u_i)_{\lambda_j^{(n+1)}}\bigr]_{i,j=1}^n \;=\; E_n(u)\cdot V_n(u-1)\cdot \mathfrak s_\mu^{(n)}(u-1)$$
where $u-1 := (u_1-1,\dots,u_n-1)$. And $V_n(u-1) = V_n(u)$ (Vandermonde is translation-invariant).

Also, $V_{n+1}(u_1,\dots,u_n,0) = V_n(u_1,\dots,u_n)\cdot \prod_i(u_i - 0) = V_n(u)\cdot E_n(u)$. So
$$\mathfrak s_\mu^{(n+1)}(u_1,\dots,u_n,0) \;=\; \frac{E_n(u)\,V_n(u)\,\mathfrak s_\mu^{(n)}(u-1)}{V_n(u)\,E_n(u)} \;=\; \mathfrak s_\mu^{(n)}(u-1)\;. \quad\square$$

### 3.3 Consequence for $\Psi_b$

Sum over $\mu\vdash 2b$ with Kostka weights:
$$\Psi_b^{(n+1)}\big|_{u_{n+1}=0} = \sum_{\mu\vdash 2b,\ \ell(\mu)\le n+1}\!\!\! K_{\mu'(2^b)}\,\mathfrak s_\mu^{(n+1)}(u_1,\dots,u_n,0) \;=\; \sum_{\mu\vdash 2b,\ \ell(\mu)\le n}\!\!\! K_{\mu'(2^b)}\,\mathfrak s_\mu^{(n)}(u-1)$$
(the terms with $\ell(\mu)=n+1$ vanish). The RHS is exactly $\Psi_b^{(n)}(u-1) = \tau_n^{-1}(\Psi_b^{(n)})$ viewed as a polynomial in $E$'s (after using $e_k(u-1)$ formulas). $\square$

### 3.4 The formula for $\tau_n^{-1}$

Standard: $\prod_i(1+t(u_i+c)) = (1+tc)^n \prod_i(1 + \tfrac{t}{1+tc}u_i)$. Extracting $[t^k]$:
$$e_k(u+c) \;=\; \sum_{j=0}^k \binom{n-j}{k-j}\,c^{k-j}\,e_j(u)\;.$$
At $c=-1$: $\tau_n^{-1}(E_k) = \sum_j \binom{n-j}{k-j}(-1)^{k-j}E_j$. $\square$

**Computational check (script `scratch/day172/test_stability.py`):** verified for $n\in\{2,3,4,5\}$, $b\in\{0,1,2,3\}$ — all 16 cases pass. This is a sanity check on top of the proof above.

---

## 4. The recursion from stability + (A) — proof

**Lemma (top-ρ symbol of $\tau_n^{-1}$).** As an operator on $\mathbb Q[E_1,\dots,E_n]$ graded by ρ, the top-ρ part of $\tau_n^{-1}$ acts on generators by:
- $E_k \mapsto E_k$ if $k$ is odd (no correction of same ρ-weight);
- $E_{2m} \mapsto E_{2m} - (n-2m+1)\,E_{2m-1}$ if $k=2m$ is even ($E_{2m-1}$ has ρ-weight $\lceil (2m-1)/2\rceil = m$, matching $\rho(E_{2m})=m$).

*Proof.* $\tau_n^{-1}(E_k) = E_k - (n-k+1)E_{k-1} + \binom{n-k+2}{2}E_{k-2} - \cdots$. The ρ-weight of the term $\binom{n-j}{k-j}(-1)^{k-j}E_j$ is $\lceil j/2\rceil \le \lceil k/2\rceil$, with equality iff either $j=k$, or $k$ is even and $j=k-1$ (since then $\lceil (k-1)/2\rceil = k/2 = \lceil k/2\rceil$). $\square$

**Corollary.** If $P \in \mathbb Q[E_1,E_2,E_3]$ then $[\text{top ρ}]\tau_n^{-1}(P) = P\big|_{E_2\to E_2-(n-1)E_1}$ (since $\sigma_n$ acts trivially on $E_1,E_3$ and shifts $E_2$).

**Deduction of the shift-law from stability + (A).** By stability §3:
$$\Psi_b^{(n+1)}\big|_{E_{n+1}=0} \;=\; \tau_n^{-1}(\Psi_b^{(n)})\;.$$
Take the top-ρ component of both sides. By **(A)** applied to $n+1$: $\mathrm{tops}^{(n+1)}[b]\in\mathbb Q[E_1,E_2,E_3]$, so it contains no $E_{n+1}$; hence $[\text{top ρ}](\Psi_b^{(n+1)}|_{E_{n+1}=0}) = \mathrm{tops}^{(n+1)}[b]$. On the RHS: $[\text{top ρ}]\tau_n^{-1}(\Psi_b^{(n)}) = \sigma_n([\text{top ρ}]\Psi_b^{(n)}) = \sigma_n(\mathrm{tops}^{(n)}[b])$, and by (A) at $n$ this equals $\mathrm{tops}^{(n)}[b]|_{E_2\to E_2-(n-1)E_1}$. Therefore
$$\mathrm{tops}^{(n+1)}[b] \;=\; \mathrm{tops}^{(n)}[b]\Big|_{E_2\,\mapsto\,E_2-(n-1)E_1}\;.$$

Iterating from $n=3$ upward, the composition of shifts $\phi_c\circ\phi_{c'} = \phi_{c+c'}$ gives, for any $N\ge 3$,
$$\mathrm{tops}^{(N)}[b] \;=\; \mathrm{tops}^{(3)}[b]\Big|_{E_2\,\mapsto\,E_2\,-\,\bigl(\sum_{n=3}^{N-1}(n-1)\bigr)E_1} \;=\; \mathrm{tops}^{(3)}[b]\Big|_{E_2\,\mapsto\,E_2\,-\,c_N E_1}$$
with $c_N = \sum_{k=2}^{N-2}k = \binom{N-1}{2}-1$. $\square$

---

## 5. Sub-claim (A): status

**(A)** $\mathrm{tops}^{(n)}[b] \in \mathbb Q[E_1,E_2,E_3]$ for all $n\ge 3$, $b\ge 0$.

### 5.1 Numerical evidence

Verified by direct computation (`scratch/day172/verify_short.py`, using `scratch/clio_check/psi_n.py` for $\Psi^{(n)}(e_2^b)$) for all $(n,b) \in \{3,4,5,6\}\times\{0,\dots,6\}$: 28/28 pass **both** (A) and the shift-law. (The prior Day-169 sweep verified 26/26 for $(n,b)\in\{3,\dots,7\}\times\{0,\dots,6\}$.)

### 5.2 Why (A) is subtle

(A) is *not* a consequence of any of the following:
1. **The u-degree bound** $\deg_u\Psi_b\le 2b$. A monomial $E_1^a E_4$ has u-degree $a+4$ and ρ-weight $a+2$, so for $a=b-2$ it has ρ$=b$ and u-degree $b+2\le 2b$ for $b\ge 2$ — structurally allowed.
2. **Individual $\mathfrak s_\mu^{(n)}$ having $E_1,E_2,E_3$ support.** Direct calculation (`scratch/day172/check_smu_supp.py`) shows $\mathfrak s_{(3,1)}^{(4)}$, $\mathfrak s_{(2,1,1)}^{(4)}$, $\mathfrak s_{(1^4)}^{(4)} = E_4$, etc., all involve $E_4,\dots$ individually. **The cancellation of $E_4$ terms in the Kostka-weighted sum $\sum_\mu K_{\mu'(2^b)}\mathfrak s_\mu^{(n)}$ is a genuine cancellation**, not a degree bound.
3. **Restriction to $E_{n+1}=0$ (from stability).** Setting $E_{n+1}=0$ only removes the $E_{n+1}$-part; it does not automatically kill $E_4,\dots,E_n$ contributions.

### 5.3 Concrete illustration at $n=4$, $b=4$

Full $\Psi_4^{(4)}$ has $E_4$-terms at *sub-top* ρ:
$$\Psi_4^{(4)} \;\supset\; -510\,E_1 E_4 + 80\,E_2 E_4 + 2180\,E_4 \quad(\text{all }\rho \le 3)\;.$$
But *at* ρ$=4$, no $E_4$ appears — the ρ$=4$ slice is entirely in $\mathbb Q[E_1,E_2,E_3]$ (9 monomials verified). So (A) says: for this specific b=4, the ρ=4 combinations happen to give zero $E_1^2 E_4$, $E_1 E_2 E_4$, $E_2^2 E_4$, $E_3 E_4$, $E_4^2$ coefficients, despite $E_4$ appearing at ρ$=3$ and ρ$=2$.

### 5.4 What a proof of (A) would look like

(A) is equivalent to: **the top-ρ symbol of the Pieri operator** $\mathcal B_2^{(n)} = V_n^{-1} e_2(\hat u) V_n$ (Day 149 Corollary E) **preserves the subring** $\mathbb Q[E_1,E_2,E_3]$. If proved, since $\mathrm{tops}^{(n)}[0]=1\in\mathbb Q[E_1,E_2,E_3]$ and $\mathrm{tops}^{(n)}[b+1] = [\text{top-ρ symbol of }\mathcal B_2^{(n)}](\mathrm{tops}^{(n)}[b])$, induction on $b$ closes (A).

Alternatively (A) follows from any *explicit* closed form for $\mathrm{tops}^{(n)}[b]$; the Day 130 EGF closed form
$$\sum_{b\ge 0}\mathrm{tops}^{(n)}[b]\frac{T^b}{b!} \;\stackrel{?}{=}\; (1+E_1 T)^{E_2/E_1 - \binom{n-1}{2}}\,\exp\!\Bigl(E_3\Bigl[\tfrac{T}{E_1(1+E_1T)^2} - \tfrac{\log(1+E_1T)}{E_1^2}\Bigr]\Bigr)$$
manifestly involves only $E_1,E_2,E_3$. Verified computationally for $n=3,\dots,7$, $b=0,\dots,6$ (35 cases, `scratch/day172/test_general_n_egf.py`). This EGF form is *equivalent* to the E₂-shift + (A) via the base case at $n=3$ (Day 130, proved through $T^8$).

---

## 6. Verification summary

| claim | range checked | script | result |
|---|---|---|---|
| Stability §3 | $n\in\{2..5\}$, $b\in\{0..3\}$ | `test_stability.py` | 16/16 ✓ |
| Shift-law + (A) | $n\in\{3..6\}$, $b\in\{0..6\}$ | `verify_short.py` | 28/28 ✓ |
| Closed-form EGF (equivalent) | $n\in\{3..7\}$, $b\in\{0..6\}$ | `test_general_n_egf.py` | 35/35 ✓ |
| (A) alone | Same as above | Same | ✓ (never violated) |
| Day-131 base at n=3 | $b\in\{0..6\}$ | `day131_work/step3_top_projection.py` | ✓ |

The prior Day-169 sweep at 26/26 is fully re-verified and extended.

---

## 7. Status against the target

- Sought: $\mathrm{tops}^{(n)}[b] = \mathrm{tops}^{(3)}[b]|_{E_2\to E_2 - c_n E_1}$, promotion to `proved`.
- Achieved: **conditional proof modulo (A)**. Stability identity §3 is rigorous. Deduction §4 is rigorous *given* (A). Sub-claim (A) is verified computationally on all tested cases but not proved.
- Registry: **grade `sketched`** — the shape of the argument is complete and one hole is precisely named. Cannot promote to `proved` without (A).

## 8. What's new versus Day 131

Day 131 §4 proved the top-ρ recursion at $n=3$ using a $\sigma_{\rm top}$ projection. That argument gave the b-recursion at *fixed n=3*, with $\sigma_{\rm top}$: $E_2 \to E_2 - 2E_1$ (a specific instance of my $\sigma_n$ at $n=3$: $E_2\to E_2 - 2E_1$ ✓ — sanity check).

The Day 131 argument doesn't obviously port to general $n$ because the *derivation* of the recursion used the n=3 Riccati / master-curve structure that is n=3-specific. My route bypasses this: instead of porting the b-recursion, I use an $n$-recursion (via factorial-Schur stability), which is dimension-agnostic. Together with (A) — the *only* place n=3 vs. general n plays a role — the argument closes.

## 9. Where (A) might come from (open)

Two speculative routes for proving (A):
1. **Direct combinatorics on Kostka + Stirling numbers.** Since $\mathfrak s_\mu^{(n)} - s_\mu$ has explicit Stirling coefficients (factorial-Schur→Schur expansion), the $E_4$ coefficient of $\sum_\mu K_{\mu'(2^b)}\mathfrak s_\mu^{(n)}$ is a specific alternating sum. (A) is an identity like "this sum vanishes". Might be provable via Lindström–Gessel–Viennot or similar bijective methods.
2. **Top-ρ symbol of $\mathcal B_2^{(n)}$ as a differential operator on $\mathbb Q[E_1,E_2,E_3]$.** If the top-ρ of $\mathcal B_2^{(n)}$ is (as expected) a degree-2 differential operator in $E_1,E_2,E_3$-derivatives (analogous to Day 131's projected recursion), it manifestly preserves the subring and induction closes.

Route 2 is closer to Day 131's approach and is probably the cleanest.

---

## 10. Pre-registered predictions (from PROVE.md) — postmortem

- ✓ "The proof will need $b\ge 0$ arbitrary — no case split on parity of $b$": correct, both stability and the recursion are $b$-uniform.
- ✗ "$c_n$ will emerge from a counting of pair-differences of variables, NOT from a subtraction-of-diagonal": *partial*. The formula $c_n = \binom{n-1}{2} - 1$ emerges as $\sum_{k=2}^{n-2} k$ from the *iteration* of the atomic shift $-(n-1)E_1$, and the "(n-1)" itself comes from $\binom{n-1}{1}$ in the $\tau_n^{-1}$ expansion — that IS "pairs $(i,j)$" for the shifted variable and one other. So half-right.
- ✗ "The proof will fit in ≤ 2 pages once the mechanism is identified": *this note is longer, but the actual proof body (§3.2, §4 deduction) is ~1 page.* The rest is (A) discussion and verification.
- **Rule 11 scorecard**: *did not fire cleanly.* The stability identity §3 is close to "unfold the definition" (unfolding $\mathfrak s_\mu^{(n+1)}|_{u_{n+1}=0}$ via the falling-factorial identity), so partial fire. Score: **1-0 (partial)** for the new arc, but the full theorem is stuck on (A).
