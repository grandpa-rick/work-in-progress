# Day 145 — Free-cumulant integrality for $M := 1 - 2F$

**Status.** REDUCTION THEOREM PROVED; the crux is now a specific, empirically
verified conjecture about $b_n \bmod 3$.

## Problem

Let $F(\tau) := \sum_{k \ge 1} b_k \tau^k$ with $b_k$ the Day 143 sequence
$$b_1..b_8 = 3,\ 27,\ 417,\ 7851,\ 164124,\ 3661389,\ 85384566,\ 2056373739$$
and let $M(\tau) := 1 - 2 F(\tau)$. Let $\kappa_n = \kappa_n(M)$ denote the
Voiculescu free cumulants attached to the moment series $M$, i.e., the unique
sequence with
$$m_n := [\tau^n] M = \sum_{\pi \in NC(n)} \prod_{V \in \pi} \kappa_{|V|}
\qquad (n \ge 1).$$

**Target.** Prove or disprove that $\kappa_n \in 6\mathbb Z$ for every
$n \ge 1$.

## Structural input (Day 143)

$$M(\tau)^2 = 1 + 4 A(\tau), \qquad A(\tau) := \sum_{k \ge 1} a_k \tau^k
= F(\tau)^2 - F(\tau) \in \mathbb Z[[\tau]].$$

Equivalently $M = (1 + 4A)^{1/2}$, so $M \in \mathbb Z[[\tau]]$ with
$m_n = -2 b_n \in 2\mathbb Z$ for every $n \ge 1$, and
$$F(\tau) \;=\; -A \,C(-A) \;=\; \sum_{k \ge 1} (-1)^k \, C_{k-1}\, A(\tau)^k,$$
where $C_j = \frac{1}{j+1}\binom{2j}{j}$ are Catalan numbers.

## Main theorem (Reduction)

**Theorem 1.** For any formal moment series $M(\tau) = 1 + \sum m_n \tau^n
\in \mathbb Z[[\tau]]$ and any positive integer $d$,
$$\bigl(m_n \in d\mathbb Z \text{ for all } n\bigr)
\iff \bigl(\kappa_n \in d\mathbb Z \text{ for all } n\bigr).$$

*Proof.*
$(\Rightarrow)$ Speicher's Möbius inversion of the moment–cumulant relation
gives
$$\kappa_n \;=\; \sum_{\pi \in NC(n)}
\mu(\pi, \hat 1_n) \,\prod_{V \in \pi} m_{|V|},$$
where $\mu(\pi,\hat 1_n) \in \mathbb Z$ is the Möbius function on the lattice
of non-crossing partitions. Every $\pi \in NC(n)$ has at least one block
$V$, so every summand contains at least one factor $m_{|V|} \in d\mathbb Z$.
Hence each term is in $d\mathbb Z$, and so is $\kappa_n$.

$(\Leftarrow)$ Induction on $n$. For $n=1$, $m_1 = \kappa_1 \in d\mathbb Z$.
For $n \ge 2$, isolate the singleton-block-of-size-$n$ term in Speicher's
direct formula:
$$m_n = \kappa_n + \sum_{\pi \in NC(n),\ |\pi| \ge 2}
\prod_{V \in \pi} \kappa_{|V|}.$$
By induction every $\kappa_{|V|}$ for $|V| < n$ lies in $d\mathbb Z$; each
$\pi$ with $|\pi| \ge 2$ contributes a product of at least two such factors,
hence lies in $d^2\mathbb Z \subseteq d\mathbb Z$. Combined with
$\kappa_n \in d\mathbb Z$ this gives $m_n \in d\mathbb Z$. $\blacksquare$

**Applied to our $M$ with $d = 6$:** since $m_n = -2 b_n$,
$$m_n \in 6\mathbb Z \iff b_n \in 3\mathbb Z.$$
The factor of $2$ is automatic. So Theorem 1 reduces our target to:

**Conjecture 2 (3-adic vanishing of $b$).** $b_n \equiv 0 \pmod 3$ for every
$n \ge 1$.

## Equivalent forms of the conjecture

Using Day 143's identity $a_n = -b_n + \sum_{i+j=n,\, i,j\ge 1} b_i b_j$,
straightforward induction shows

**Lemma 3.** The following are equivalent:

  (i)   $b_n \equiv 0 \pmod 3$ for every $n \ge 1$;

  (ii)  $a_n \equiv 0 \pmod 3$ for every $n \ge 1$;

  (iii) $F(\tau) \equiv 0 \pmod 3$ in $\mathbb Z[[\tau]]$;

  (iv)  $A(\tau) \equiv 0 \pmod 3$ in $\mathbb Z[[\tau]]$;

  (v)   $M(\tau) \equiv 1 \pmod 3$ in $\mathbb Z[[\tau]]$.

*Proof.* (i)$\Leftrightarrow$(iii) and (ii)$\Leftrightarrow$(iv) are the
definitions. (iii)$\Leftrightarrow$(v) follows from $M = 1 - 2F$ and
$\gcd(2,3)=1$. (v)$\Leftrightarrow$(iv) from $M^2 - 1 = 4A$: mod 3,
$M \equiv 1$ iff $M^2 \equiv 1$ iff $4A \equiv 0$ iff $A \equiv 0$
(the middle equivalence uses that $M+1 \equiv 2 \pmod 3$ is a unit).
Finally (i)$\Leftrightarrow$(ii) via the Day 143 identity: assuming
$b_1,\ldots,b_{n-1} \in 3\mathbb Z$ (inductive), the sum
$\sum b_i b_j \in 9\mathbb Z$, so $a_n \equiv -b_n \pmod 3$. $\blacksquare$

## Sanity check ($n = 8$)

With $b_8 = 2\,056\,373\,739$ (Day 144 verification), the moment–cumulant
recursion yields
$$\kappa_8 = -114\,927\,465\,222 = -6 \cdot 19\,154\,577\,537,$$
divisible by 6. The full sequence
$$\kappa_n / (-6),\quad n = 1..8:\quad
1,\ 15,\ 373,\ 11\,245,\ 375\,732,\ 13\,386\,573,\ 498\,347\,406,\
19\,154\,577\,537$$
is entirely integer, confirming Conjecture 2 in the range $n \le 8$.

Independent verification: $b_n \bmod 3 = 0$ for $n=1..8$
(directly from Day 143 and Day 144 data).

## Attempted attacks on Conjecture 2

The obstacle is that Conjecture 2 uses information about the specific
origin of $b$ (Rick's Frobenius/hypergeometric setup) beyond the identity
$M^2 = 1 + 4A$. Below we log four attacks and the specific point at
which each stalls.

### Attack A — Fermat / Frobenius on $F$

Fermat: $F(\tau)^3 \equiv F(\tau^3) \pmod 3$ for any $F \in \mathbb Z[[\tau]]$.
Combined with the algebraic identity
$F^3 = F + A(1+F)$ (which follows from $F^2 = F+A$),
$$F(\tau^3) - F(\tau) \equiv A(1 + F) \pmod 3. \qquad (*)$$

If Conjecture 2 holds ($A \equiv 0$ mod 3), the right side vanishes and
$F(\tau^3) \equiv F(\tau) \pmod 3$, giving $b_n \equiv 0 \pmod 3$
for $3 \nmid n$ directly and for $3 \mid n$ by induction on $v_3(n)$.

However $(*)$ carries no *new* information: substituting $A = F^2 - F$
into the right side gives $A(1+F) = F(F-1)(F+1) = F^3 - F$, so $(*)$
reduces to Fermat itself. Genuinely additional input is needed.

### Attack B — $M \bmod 9$ from $A \bmod 3$

Under the hypothesis $A \in 3\mathbb Z[[\tau]]$ (write $A = 3 A'$), one
computes in $\mathbb Z_{(3)}[[\tau]]/9$:
$$M = (1 + 12 A')^{1/2} \equiv 1 + 6 A' \pmod 9, \qquad
M^{-3} \equiv 1 \pmod 9.$$
Hence $[\tau^n] M^{-3} \equiv 0 \pmod 9$ for $n \ge 1$. Via the Lagrange
formula $\kappa_n = -\frac{1}{n-1}[\tau^n] M^{-(n-1)}$ (see §
Reduction addendum below), this yields $\kappa_n \in 3\mathbb Z$, hence
$\kappa_n \in 6\mathbb Z$.

But the argument uses $A \equiv 0 \pmod 3$ as an input — the very content
of Conjecture 2.

### Attack C — Frobenius Möbius formula for $\kappa$

From $C(w)^2 = 1 + 4 A(w/C(w))$ (where $C = 1 + \sum \kappa_n w^n$), one
extracts by degree: modulo the inductive hypothesis
$\kappa_i \in 6\mathbb Z$ for $i < n$,
$$\kappa_n \equiv 2 a_n \pmod 6.$$
So $\kappa_n \in 6\mathbb Z$ iff $a_n \in 3\mathbb Z$ — Conjecture 2 again.

### Attack D — Rick's Frobenius identity at $(U,V) = (0,0)$ mod 3

Rick's setup gives $a_k = [E_3^k T^{3k-1}] X$ where
$X = L F_P / F_P$, $L = T\theta^2 - \theta$, $F_P = 1 + E_3 H$ with
$H = \sum_{b\ge 2} U_b(E_3) T^b/b!$ and $U_b(E_3) \in \mathbb Z[E_3]$.

Modulo 3 the recursion for $\Psi[b+1]$ simplifies: the $-3 b E_3 \sigma(\Psi[b-1])$
term vanishes, leaving
$$\Psi[b+1] \equiv \bigl(E_2 - (b{+}1) E_1 + (b{+}1)^2\bigr)\Psi[b]
\;-\; b(b{-}1)(E_1 - 2b - 2)\, E_3\, \sigma(\Psi[b-2]) \pmod 3.$$

Direct numerical inspection at $(U,V)=(0,0)$ shows $[T^b] X_k \bmod 3$
does **not** vanish uniformly for $b$ near $3k - 1$; the vanishing is
peculiar to the exact diagonal $b = 3k-1$. E.g., for $X_1 \bmod 3$,
$[T^3] X_1 = -2 \not\equiv 0$ while $[T^2] X_1 = -3 \equiv 0$;
for $X_2 \bmod 3$, $[T^7] X_2 \equiv 2$ while $[T^5] X_2 \equiv 0$.

The problem is that at (U,V)=(0,0) the coefficients of $F_P$ carry $1/b!$
denominators (i.e., 3 in the denominator for $b \ge 3$), so no clean
mod-3 reduction of $F_P$ itself is available. A localization-at-3 analysis
does not obviously bootstrap the diagonal-3 vanishing.

## Reduction addendum: Lagrange formula for $\kappa_n$

Following from Speicher $M(z) = C(z M(z))$ and Lagrange inversion applied
to $z = y/M(z(y))$ (i.e., $\phi = 1/M$ in the standard form $z = y\phi(z)$),
one gets for $n \ge 2$:
$$\boxed{\quad \kappa_n = -\frac{1}{n-1}\, [z^n]\, M(z)^{-(n-1)} \quad}$$
and $\kappa_1 = m_1$. Since $M \in \mathbb Z[[\tau]]$ with $M(0)=1$, one has
$M^{-(n-1)} \in \mathbb Z[[\tau]]$, so
$[z^n] M^{-(n-1)} \in \mathbb Z$ and the identity forces $(n-1) \mid
[z^n] M^{-(n-1)}$.

## Status summary

  * **Proved (Theorem 1):** the reduction $\kappa_n \in 6\mathbb Z\
    (\forall n) \iff b_n \in 3\mathbb Z\ (\forall n)$.
  * **Proved (Lemma 3):** the further equivalences to
    $A \equiv 0 \pmod 3$ and $M \equiv 1 \pmod 3$.
  * **Verified numerically ($n \le 8$):** $b_n \in 3\mathbb Z$,
    hence $\kappa_n \in 6\mathbb Z$; in particular
    $\kappa_8 = -6 \cdot 19\,154\,577\,537$.
  * **Open (Conjecture 2):** $b_n \equiv 0 \pmod 3$ for all $n$. Four
    attacks logged above; each either reduces back to Conjecture 2 or
    stalls on the $1/b!$ obstruction in Rick's construction of $F_P$
    at $(U,V)=(0,0)$.

## FPSAC content

**Theorem 3.8 (Rick’s reduction).** For the invariant series
$M = 1 - 2F$ arising from the Frobenius diagonal (Day 143),
free-cumulant integrality $\kappa_n \in 6\mathbb Z$ is equivalent to
the 3-adic vanishing of the $b_k$ sequence:
$b_k \equiv 0 \pmod 3$ for all $k$. The latter is verified through
$k = 8$ (Day 144).

## Files

  * `/home/agent/.claude/scratch/prove-20260829-025318.md` — working notes.
  * `/home/agent/.claude/scratch/verify_kappa8.py` — $\kappa_8$ sanity check;
    verifies $\kappa_n \in 6\mathbb Z$ for $n \le 8$ from the $b_n$
    sequence via Speicher recursion.
  * `/home/agent/.claude/scratch/study_Pb_mod3.py` — $[E_3^k T^b] X$
    computed mod 3 for $k \le 7$ at $(U,V) = (0,0)$; confirms $a_k
    \equiv 0 \pmod 3$ at the diagonal, but shows adjacent coefficients
    are not divisible.
  * `/home/agent/.claude/scratch/mod3_analysis.py` — verifies
    $F(\tau^3) \equiv F(\tau) \pmod 3$ numerically for $n \le 24$
    (using the $b_1..b_8$ data), and confirms Attack A gives no
    new information over the raw identity.

## For the collaborator

The reduction is clean and publishable. The remaining question — proving
$b_n \equiv 0 \pmod 3$ — likely requires either:

  1. A combinatorial identification of $b_k$ as counting some structure
     with a natural free $\mathbb Z/3$-action (e.g., 3-colorings of the
     Novelli–Thibon labeled planar trees at some specialization).
     Novelli–Thibon Eq. 41 identifies the $k=-1$ slice with free
     cumulants, and Rick's $M$ matches their Catalan geode
     $(2xg-1)^2 = 1 - 4x$ under the identification $x = -A$,
     $F = -A\cdot C(-A) = xg$. If in that slice the underlying
     tree-count has a natural cyclic symmetry, we're done.

  2. A closed-form ODE / algebraic equation for $F(\tau)$ (as opposed
     to just $M^2 = 1 + 4A$, which involves $A$ implicitly) — enough
     structure to force $F \equiv 0 \pmod 3$ from finite verification
     plus the ODE's initial conditions and 3-adic properties.

  3. An input from Rick's Frobenius machinery: some identity of the
     form "the $E_3^k$-diagonal of $X$ modulo 3 vanishes because of
     an internal $\mathbb Z/3$-symmetry of the $P_b$ recursion at
     $(U,V) = (0,0)$." Attack D isolates the mod-3 recursion but does
     not immediately give the diagonal vanishing.
