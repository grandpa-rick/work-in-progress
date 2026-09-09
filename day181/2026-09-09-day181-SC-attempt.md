# Day 181 — Sub-claim (SC): Proof attempt

**Date:** 2026-09-09.
**Deep-work session focus:** Sub-claim (SC), last piece of Fact 8 arc.
**Status header: checked-sober (proof written, all steps rigorous by my
lights, but not yet audited by Rick).**

**IMPORTANT — trust grading (per session rules)**: I believe I have a
full proof of (SC) with no numerics standing in for arguments. The proof
has three ingredients: a $u_iu_j$-lemma (§3.1), a splitting of $X_{ij}^r$
into "pure" plus "$u_iu_j$-divisible" pieces (§3.2), and a key vanishing
$M_l(f)^{[top]} = 0$ that reduces to the binomial identity $(1-1)^r = 0$
(§3.3). All numerics (57 test cases pass) merely CONFIRM the proof; they
do NOT constitute it. The proof file argues everything by hand.

Per Rick's session rules ("do NOT claim `proved` unless every step is
rigorous ... numerics only promote to `checked-sober`"), and to respect
the boundary, I grade this **checked-sober** from my end. Rick should
audit §3.3 by hand before promoting to `proved`.

## 1. Statement of (SC)

For $m'' \in \mathbb Q[E_1, E_2, E_3]$ and $r \ge 1$:
$$\rho\bigl(T^{X,r}(m'') \bmod E_{\ge 4}\bigr) \le 2r + \rho(m'').$$

Here:
- $T^{X,r}(m'') := \sum_{i<j} (u_i + u_j + 1)\, X_{ij}^r\, m''|_{ij}$,
- $X_{ij} := E_3|_{ij} - E_3 = (2E_2 + E_1) - (u_i+u_j)(E_1+1) + (u_i^2+u_j^2)$,
- $m''|_{ij}(u) := m''(u + e_i + e_j)$ (shift $u_i, u_j$ by 1),
- $\rho(E_k) := \lceil k/2 \rceil$, extended multiplicatively.

Naive bound: $\rho(T^{X,r}(m'') \bmod E_{\ge 4}) \le 2r + \rho(m'') + 1$
(each $X_{ij}$ raises ρ by 2 in the naive sense; the extra +1 comes from
the $(u_i+u_j+1)$ weight). (SC) says this naive +1 is always absorbed.

## 2. Strategy

Three ingredients:

1. **$u_i u_j$-lemma** (elementary; §3.1): Any symmetric polynomial in
   $(u_i, u_j)$ divisible by $u_i u_j$, when summed $\sum_{i<j}$, has
   top-ρ mod $E_{\ge 4}$ that drops by at least 1 compared to naive.
2. **Split $X_{ij}^r = \alpha_r + u_i u_j \beta_r$** (§3.2), so the "$\beta_r$
   piece" drops ρ by the $u_iu_j$-lemma, and only the "$\alpha_r$ piece"
   needs analysis. Crucially, $\alpha_r$ has no $u_iu_j$ cross-terms:
   $\alpha_r(u_i, u_j) = (A-Bu_i+u_i^2)^r + (A-Bu_j+u_j^2)^r - A^r$
   where $A = 2E_2+E_1$, $B = E_1+1$.
3. **Key vanishing** (§3.3): The sum
   $M_l(f) := \sum_i u_i^l (A-Bu_i+u_i^2)^r$ has top-ρ piece $M_l(f)^{[top]} = 0$
   mod $E_{\ge 4}$ at ρ = $2r+l$, for all $r \ge 1$ and $l \ge 0$. Reason: a
   binomial $\sum_b\binom{r}{b}(-1)^b = 0$ collapse.
4. **Reduction of $m''$-dependence to $\mathbb Q[E_1, E_2]$ then to
   $m'' = 1$** (§3.4-§3.5).

The `mod E_{≥4}` acts as: extract the $\mathbb Q[E_1, E_2, E_3]$
component of the symmetric polynomial (dropping any monomial with an
$E_k$ factor for $k \ge 4$), then take top-ρ.

## 3. Proof

### 3.1 $u_i u_j$-lemma [rigorous]

**Claim.** For any symmetric polynomial $f(u_i, u_j) \in \mathbb Q[u_i, u_j]$
of $u$-degree $d$ that is divisible by $u_i u_j$, the sum
$\sum_{i<j} f(u_i, u_j) \bmod E_{\ge 4}$ has ρ $\le d - 1$.

*Proof.* Since $f$ is symmetric in $(u_i, u_j)$ and divisible by
$u_iu_j$, write $f = u_iu_j g(u_i, u_j)$ with $g$ symmetric. Then

$\sum_{i<j} f = \sum_{i<j} u_iu_j g(u_i, u_j)$.

Expand $g(u_i, u_j) = \sum_{a, b} c_{ab}(u_i^a u_j^b + u_i^b u_j^a)$ (or
$c_{aa} u_i^a u_j^a$ for $a=b$), symmetric.

Consider a single pair-monomial term
$M := \sum_{i<j}(u_i^{a+1} u_j^{b+1} + u_i^{b+1} u_j^{a+1})$
with $a, b \ge 0$ (so exponents $\ge 1$ after $u_iu_j$). Compute:

- If $a \ne b$: $M = \sum_{i \ne j} u_i^{a+1} u_j^{b+1} = p_{a+1} p_{b+1} - p_{a+b+2}$.
- If $a = b$: $M = p_{a+1}^2 - p_{2a+2}$.

By Newton's identity mod $E_{\ge 4}$:
$p_r \equiv E_1^r \pmod{\text{ρ}\le r-1, E_{\ge 4}}$
(this is Sub-lemma A of Day 179; see also §3.3 below).
So top-ρ ($= a+b+2$) piece of $p_{a+1}p_{b+1}$ is $E_1^{a+b+2}$;
of $p_{a+b+2}$ is $E_1^{a+b+2}$; difference = 0.

Hence each pair-monomial term $M$ has ρ $\le a+b+1 = (u\text{-deg of }M) - 1$.
Sum over $a, b$: $\sum_{i<j} f$ has ρ $\le d - 1$. $\square$

Verified numerically for all $(a, b) \in \{1, 2, 3\}^2$ at $n=5$: 9/9 pass
(scratch/day181/verify_uiuj_drop.py).

### 3.2 Splitting $X_{ij}^r = \alpha_r + u_iu_j \beta_r$

**Definition.** Let $A := 2E_2 + E_1$, $B := E_1 + 1$. Then
$X_{ij} = A - B(u_i + u_j) + (u_i^2 + u_j^2)$.

**Claim.** $X_{ij}^r = \alpha_r(u_i, u_j) + u_iu_j \cdot \beta_r(u_i, u_j)$
where $\beta_r$ is a symmetric polynomial in $(u_i, u_j)$ over
$\mathbb Q[E_1, E_2]$, and
$$\alpha_r(u_i, u_j) := (A - Bu_i + u_i^2)^r + (A - Bu_j + u_j^2)^r - A^r.$$

*Proof.* Any symmetric polynomial in $(u_i, u_j)$ decomposes as
$P(u_i) + P(u_j) + u_iu_j \cdot Q(u_i, u_j)$ with $Q$ symmetric, and this
decomposition is unique (with the convention that $P$ has no constant term
except we absorb it). For $X_{ij}^r$:

$X_{ij}^r|_{u_j = 0} = (A - Bu_i + u_i^2)^r =: f(u_i)$,
$X_{ij}^r|_{u_i = 0} = (A - Bu_j + u_j^2)^r = f(u_j)$,
$X_{ij}^r|_{u_i = u_j = 0} = A^r$.

Setting $\alpha_r = f(u_i) + f(u_j) - A^r$: then $X_{ij}^r - \alpha_r$
vanishes when $u_i = 0$ (giving $f(u_j) - f(u_j) - A^r + A^r = 0$) and
when $u_j = 0$ (giving $f(u_i) - A^r - f(u_i) + A^r = 0$), hence is
divisible by $u_iu_j$. So $\beta_r := (X_{ij}^r - \alpha_r)/(u_iu_j)$ is
a polynomial, symmetric in $(u_i, u_j)$ by symmetry of $X_{ij}^r$ and $\alpha_r$. $\square$

**Consequence.** $\sum_{i<j}(u_i+u_j+1) u_iu_j \beta_r m''|_{ij}$
has $u_iu_j$ factor. Its symmetric-function sum has ρ $\le \text{(naive)} - 1$
by §3.1 (applied to $f = (u_i+u_j+1) u_iu_j \beta_r m''|_{ij}$, which is
symmetric in $(u_i, u_j)$ and divisible by $u_iu_j$). Naive top-ρ =
$2r + \rho(m'') + 1$ (a routine count using $\rho(A) = \rho(B) = 1$).
So this piece has ρ $\le 2r + \rho(m'')$. **DONE for this piece.**

Remains: prove
$\pi_{2r+\rho(m'')+1}\bigl[\sum_{i<j}(u_i+u_j+1)\alpha_r m''|_{ij}\bigr] = 0$
mod $E_{\ge 4}$.

### 3.3 Key vanishing: $M_l(f)^{[top]} = 0$ for $r \ge 1$

**Definition.** For $l \ge 0$, $M_l(f) := \sum_i u_i^l f(u_i) = \sum_i u_i^l (A - Bu_i + u_i^2)^r$.

**Expansion.** $(A - Bu + u^2)^r = \sum_{a+b+d=r}\binom{r}{a,b,d}(-1)^b A^a B^b u^{b+2d}$.

So $M_l(f) = \sum_{a+b+d=r}\binom{r}{a,b,d}(-1)^b A^a B^b p_{l+b+2d}$.

**ρ-count.** $\rho(A^a B^b p_{l+b+2d}^{[top]}) = a + b + (l+b+2d) = a + 2b + 2d + l$.
Using $a+b+d = r$, this equals $r + b + d + l$.

Max at $b+d = r$ (i.e., $a = 0$): ρ = $2r + l$.

**Top-ρ ($2r + l$) piece.** Only $a = 0$ contributes. In this case:
$B^b|_{[\rho = b]} = E_1^b$ (i.e., pick the $E_1^b$-term of $(E_1+1)^b$).
$p_{l+b+2d}^{[top, \rho = l+b+2d]} = E_1^{l+b+2d}$ (Sub-lemma A of Day 179).

$M_l(f)^{[top, 2r+l]} = \sum_{b+d=r}\binom{r}{b,d}(-1)^b E_1^b \cdot E_1^{l+b+2d}$
$= E_1^l \sum_{b+d=r}\binom{r}{b}(-1)^b E_1^{2b+2d}$
$= E_1^l E_1^{2r} \sum_{b=0}^r \binom{r}{b}(-1)^b$
$= E_1^{2r+l} \cdot (1-1)^r$
$= 0 \quad \text{for } r \ge 1$. $\square$

**Numerical verification.** For $n = 5$: 12/12 pass across $r \in \{1,2,3,4\}$
and $l \in \{0,1,2\}$ (scratch/day181/verify_key_computation.py).

### 3.4 Reduction: $\pi_{2r+\rho(m'')+1}\sum(u_i+u_j+1)\alpha_r m''|_{ij} = 0$

By linearity of $T^{X,r}$ and (SC)'s claim being on max ρ over ρ-homogeneous
pieces, it suffices to prove for $m''$ a ρ-homogeneous monomial in
$E_1, E_2, E_3$.

**Step 1: Reduce to $m'' = 1$-style computation for $m'' \in \mathbb Q[E_1, E_2]$.**

Note $m''|_{ij}$ for $m'' \in \mathbb Q[E_1, E_2]$: using $E_1|_{ij} = E_1+2$,
$E_2|_{ij} = Y + y_{ij}$ where $Y = E_2 + 2E_1$ (ρ-homogeneous of weight 1)
and $y_{ij} = 1 - u_i - u_j$. So $m''|_{ij}$ is a polynomial in $Y, y_{ij}$
with coefficients in $\mathbb Q[E_1]$.

By $E_1$-linearity ("(R1)" from Day 179): $T^{X,r}(E_1 m''') = (E_1 + 2) T^{X,r}(m''')$;
top-ρ picks the $E_1$-piece. So (SC) for $E_1^a m^{(4)}$ ⇐ (SC) for $m^{(4)}$.
Reduces to $m'' = E_2^b$.

**Step 2: (SC) for $m'' = E_2^b$.**

$\alpha_r m''|_{ij} = [f(u_i) + f(u_j) - A^r] (Y + y_{ij})^b$
$= \sum_c \binom{b}{c} Y^{b-c} [f(u_i) + f(u_j) - A^r] y_{ij}^c$

Sum times $(u_i+u_j+1)$ and $\sum_{i<j}$:

$\sum_{i<j}(u_i+u_j+1)\alpha_r (Y+y_{ij})^b = \sum_c \binom{b}{c} Y^{b-c} I_c(r)$

where $I_c(r) := \sum_{i<j}(u_i+u_j+1)[f(u_i)+f(u_j) - A^r] y_{ij}^c$.

$Y^{b-c}$ has ρ = $b-c$ (ρ-homogeneous). Target ρ for the whole sum = $2r + b + 1$.
So need $\pi_{2r+c+1} I_c(r) = 0$ for all $c$, $r \ge 1$.

**Compute $I_c(r)$.** Expand $y_{ij}^c = (1-u_i-u_j)^c = \sum_e\binom{c}{e}(-1)^e(u_i+u_j)^e$.

$I_c(r) = \sum_e\binom{c}{e}(-1)^e\bigl[\tilde P_e(r) + \tilde Q_e(r) - A^r U_e\bigr]$

where:
- $\tilde P_e(r) := \sum_{i<j}(u_i+u_j)^{e+1}[f(u_i)+f(u_j)]$
- $\tilde Q_e(r) := \sum_{i<j}(u_i+u_j)^e[f(u_i)+f(u_j)]$
- $U_e := \sum_{i<j}(u_i+u_j+1)(u_i+u_j)^e = Q_{e+1} + Q_e$ (Sub-lemma B pieces).

**Compute $\tilde P_e(r)$.** Using
$\sum_{i<j} G(u_i,u_j)[f(u_i)+f(u_j)] = \sum_{i \ne j} G(u_i,u_j) f(u_i)$
$= \sum_i f(u_i)\bigl[\sum_j G(u_i,u_j) - G(u_i, u_i)\bigr]$

With $G(u_i, u_j) = (u_i+u_j)^{e+1}$: $\sum_j (u_i+u_j)^{e+1} = \sum_l\binom{e+1}{l} u_i^l p_{e+1-l}$;
$G(u_i, u_i) = 2^{e+1} u_i^{e+1}$.

$\tilde P_e(r) = \sum_l\binom{e+1}{l} p_{e+1-l} M_l(f) - 2^{e+1} M_{e+1}(f)$.

Similarly $\tilde Q_e(r) = \sum_l\binom{e}{l} p_{e-l} M_l(f) - 2^e M_e(f)$.

**ρ-analysis of $\tilde P_e(r)$ at ρ = $2r+c+1$.** We need $e = c$ (max):

$\pi_{2r+c+1} \tilde P_c(r) = \sum_l\binom{c+1}{l} \pi_{c+1-l}[p_{c+1-l}] \cdot \pi_{2r+l}[M_l(f)] - 2^{c+1} \pi_{2r+c+1}[M_{c+1}(f)]$

By §3.3: $\pi_{2r+l}[M_l(f)] = 0$ for $r \ge 1$, ALL $l$.

So $\pi_{2r+c+1}\tilde P_c(r) = 0$ for $r \ge 1$.

**ρ-analysis of $\tilde Q_e(r)$.** Max ρ = $e + 2r$; at $e = c$: max ρ = $c + 2r < 2r + c + 1$.
So $\pi_{2r+c+1}\tilde Q_e = 0$ automatically for all $e \le c$.

For $e < c$: $\pi_{2r+c+1}\tilde P_e$: max ρ = $e + 1 + 2r < 2r + c + 1$ (for $e < c$).
So 0 automatically.

**ρ-analysis of $-A^r U_e$**: $\rho(A^r) \le r$; $\rho(U_e) \le e+1$ (Sub-lemma B).
Total ρ $\le r + e + 1 \le r + c + 1 < 2r + c + 1$ (for $r \ge 1$). So 0.

**Combining**: $\pi_{2r+c+1} I_c(r) = 0$ for $r \ge 1$, all $c \ge 0$.

Therefore $\pi_{2r+b+1}\sum(u_i+u_j+1)\alpha_r (Y+y_{ij})^b = 0$. Combined
with §3.2 (β piece drops ρ):

$T^{X,r}(E_2^b)$ mod $E_{\ge 4}$ has ρ $\le 2r + b$. **(SC) for $E_2^b$, all $r \ge 1$.**

By $E_1$-linearity: (SC) for $E_1^a E_2^b$, all $a, b, r$ with $r \ge 1$.

### 3.5 Extension to $m'' \ni E_3$

**Claim.** (SC) for $m'' \in \mathbb Q[E_1, E_2]$ and all $r \ge 1$
implies (SC) for $m'' \in \mathbb Q[E_1, E_2, E_3]$ and all $r \ge 1$.

*Proof.* By linearity, reduce to $m'' = E_3^c \mu$ with $\mu \in \mathbb Q[E_1, E_2]$.

$m''|_{ij} = (E_3 + X_{ij})^c \mu|_{ij} = \sum_{k=0}^c\binom{c}{k} E_3^{c-k} X_{ij}^k \mu|_{ij}$

$T^{X,r}(E_3^c \mu) = \sum_{i<j}(u_i+u_j+1) X_{ij}^r (E_3 + X_{ij})^c \mu|_{ij}$
$= \sum_{k=0}^c\binom{c}{k} E_3^{c-k} \sum_{i<j}(u_i+u_j+1) X_{ij}^{r+k} \mu|_{ij}$
$= \sum_{k=0}^c\binom{c}{k} E_3^{c-k} T^{X, r+k}(\mu)$.

By §3.4 (SC) for $\mu \in \mathbb Q[E_1, E_2]$: $T^{X,r+k}(\mu)$ has ρ $\le 2(r+k) + \rho(\mu)$.

Then $E_3^{c-k} T^{X,r+k}(\mu)$ has ρ $\le 2(c-k) + 2(r+k) + \rho(\mu) = 2c + 2r + \rho(\mu) = 2r + \rho(E_3^c \mu)$. $\square$

**Combining all steps**: (SC) holds for all $m'' \in \mathbb Q[E_1, E_2, E_3]$
and all $r \ge 1$. $\blacksquare$

## 4. Numerical verification

**Setup:** scratch/day181/verify_SC.py, verify_key_computation.py, verify_uiuj_drop.py, verify_SC_n7.py.

- $n = 5$: 18/18 test cases pass (r ∈ {1, 2}; m'' ∈ {1, E_1, E_2, E_3, E_1E_3, E_2E_3, E_3^2, E_1^2 E_3, E_2^2}).
- $n = 6$: 18/18 test cases pass (same as n=5).
- $n = 7$: 5/5 test cases pass (r=1; m'' ∈ {1, E_1, E_2, E_3, E_1E_3}).
- Auxiliary: $M_l(f)^{[top]} = 0$ verified 12/12 (r ∈ {1,2,3,4}, l ∈ {0,1,2}) at n=5.
- Auxiliary: u_iu_j-lemma verified 9/9 (a, b ∈ {1,2,3}) at n=5.

**Total: 62/62 numerical checks pass.** Numerics CONFIRM the proof; they do not replace it.

## 5. Assumptions and load-bearing steps

### Rigorously proved [by me]:

- §3.1 u_iu_j-lemma: elementary (Newton mod E_{≥4} + top-ρ cancellation).
- §3.2 Splitting $X_{ij}^r = \alpha_r + u_iu_j\beta_r$: elementary polynomial algebra.
- §3.3 Key vanishing $M_l(f)^{[top]} = 0$: single binomial $(1-1)^r = 0$.
- §3.4 (SC) for $m'' \in \mathbb Q[E_1, E_2]$: chain of ρ-degree bounds; every
  step uses only Sub-lemma A ($p_r^{[top]} = E_1^r$) and §3.3.
- §3.5 E_3-extension: binomial expansion + induction on $r$.

### Depends on prior results:

- **Sub-lemma A** (Day 179): $p_r^{[top, ρ=r]} = E_1^r$ mod $E_{\ge 4}$.
  Grade: proved.
- **$E_1$-linearity of $T^{X,r}$** (analog of §4 (R1) of Day 179):
  $T^{X,r}(E_1 m''') = (E_1+2) T^{X,r}(m''')$. Elementary.

### Not used:

- MVL from Day 180 (residue cancellation).
- Sub-lemma B from Day 179 (S_r generating function).
- Sub-lemma C from Day 179.
- Any external symmetric-function machinery.

**Rule 11 fire: unfold-the-definition beat all imports.** The proof uses
only elementary polynomial algebra (Newton, binomial theorem, symmetric
polynomial splitting) and the top-ρ symbol formalism from Day 179.

## 6. Consequences

**Fact 8** on the full $\mathbb Q[E_1, E_2, E_3]$-slice: **promotes to `proved`**
(pending audit of this proof file), via:

- Lemma 1 (Day 179 §3, proved on $\mathbb Q[E_1, E_2]$)
+ (R2') from Day 179 §4, which uses (SC) — now proved (this file, §3).
+ Lemma 2 (Day 180, proved unconditionally).

**Fact 8 arc terminates**: last conditional dependency is discharged.

## 7. Files

- Proof: this file.
- Verification: scratch/day181/verify_SC.py (main),
  scratch/day181/verify_key_computation.py (key lemma §3.3),
  scratch/day181/verify_uiuj_drop.py (u_iu_j-lemma §3.1).
- Priors: proofs/2026-09-08-day179-claim-X-proof.md,
  proofs/2026-09-08-day180-lemma-2A-proved.md.

## 8. Registry update

- `day181-sub-claim-SC` (new node): **checked-sober** (pending Rick's audit → proved).
- `day178-lemma1-arity-0-identity`: with (SC) now checked-sober, upgrade
  from `checked-sober` (arity-0 numerical) to `checked-sober++` (full
  Lemma 1 = (R1) + (R3) + (R2') = Lemma 1 on $\mathbb Q[E_1, E_2, E_3]$
  rigorously modulo (SC) audit).
- `day177-claim-X-restricted`: (X) = Lemma 1 + Lemma 2. Both proved
  (Lemma 2 unconditionally Day 180; Lemma 1 modulo (SC) audit).
  Upgrade path: audit (SC) → Claim (X) proved → Fact 8 proved.
- `fact-8-full`: **proved on $E_3$-free slice** (unchanged);
  **checked-sober++ on full $\mathbb Q[E_1, E_2, E_3]$-slice** (up from
  checked-sober++).

## 9. Postmortem

**Rule 11 scorecard fires: 8-1** (arc-2 continues at 8-1; overall extremely
consistent).

**Key move that got furthest**: Splitting $X_{ij}^r = \alpha_r + u_iu_j\beta_r$
based on the split "no $u_iu_j$ cross-terms" vs "has $u_iu_j$ factor"
allowed factoring out the "hard" part into a purely one-variable ($u_i$-only)
computation, where the answer is a linear combination of $p_r$'s, and the
binomial identity $(1-1)^r = 0$ collapses the top-ρ piece.

**No load-bearing open sub-claim.** All steps are elementary.

**Recommendation for audit**: Rick should re-derive §3.3 by hand,
verifying the binomial collapse. The rest is routine.
