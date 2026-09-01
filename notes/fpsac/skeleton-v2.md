# FPSAC 2027 Extended Abstract — Skeleton v2

**Date drafted:** 2026-09-01 (Day 153, FPSAC writing kickoff)
**Supersedes:** `skeleton.md` (2026-08-25). The v1 skeleton was aimed at the Day 133 "density and sign for tops[b]" result; since Day 143 the arc has been reorganized around the master curve, the ψ closed form (proved Day 152, audited Day 152b), and the b_k mod 3 corollary (proved Day 148). The old skeleton stays on disk for the density-and-sign material, but the framing is wrong.
**Target:** 12 pages (FPSAC EA), deadline 2026-11-15, 76 days.
**Working title:** *A master curve for the shifted-Schur map: closed form for a Lagrange kernel, and a mod-3 congruence.*

---

## The story in one paragraph

Let $\Psi: \Lambda_{\text{sym}} \to \Lambda^*$ be the map $\Psi(f)(u) = \mathcal{T}(f \cdot V)/V$ where $\mathcal{T}$ is the falling-factorial umbral map $u^n \mapsto (u)_n$ and $V$ the Vandermonde in three variables $u_1, u_2, u_3$. **Theorem A:** $\Psi$ is the linear map $s_\mu \mapsto \mathfrak{s}_\mu$ from Schur functions to Macdonald's factorial Schur functions. This identifies $\Psi(e_2^b)$ as a Kostka-weighted sum of factorial Schurs. The interior series $F_P := \Psi^+(e^{T e_2})/V$ satisfies a three-variable Riccati system whose **master curve** governs two arcs. **Downward specialisation** (diagonal $u_i = u$, $T = 1$): $F$ is algebraic of degree 5, and the coefficients $b_k$ of an associated inversion satisfy $b_k \equiv 0 \pmod{3}$ (Theorem B). **Upward specialisation** (top weight of $H = \tau(F_P)/F_P$): the Lagrange kernel $\psi$ of the leading symbol is algebraic of degree exactly 5, with explicit closed form (Theorem C) and minimal polynomial $Q(\psi, Y) \in \mathbb{Z}[E][\psi, Y]$ satisfying $\deg_\psi Q = 5$, $\deg_Y Q = 9$ (Theorem D). At $E_3 = 0$ the top layer of $[T^n]H$ recovers the Narayana polynomials.

---

## §1. Introduction (1.5 pp)

- **The setup.** Three commuting variables $u_1, u_2, u_3$; the falling-factorial map $\mathcal{T}: u_i^n \mapsto (u_i)_n$; the Vandermonde $V = \prod_{i<j}(u_i - u_j)$; the map $\Psi(f) := \mathcal{T}(f V)/V$ on symmetric functions. This is the three-variable case of the Molev–Olshanski shifted-Schur map.
- **The main results, stated in the intro.**
  - **Theorem A** (Day 149). $\Psi(s_\mu) = \mathfrak{s}_\mu$ where $\mathfrak{s}_\mu$ is Macdonald's factorial Schur $s_\mu(u \mid a)$ at $a_l = 1 - l$ (equivalently the Okounkov–Olshanski shifted Schur function).
  - **Theorem B** (Day 148). The interior sequence $b_k$ obtained by Lagrange inversion from the diagonal specialisation of $F_P = \Psi^+(e^{T e_2})/V$ satisfies $b_k \equiv 0 \pmod 3$ for all $k \ge 1$.
  - **Theorem C** (Day 152). The Lagrange kernel $\psi(Y)$ of the leading symbol of $H = \tau(F_P)/F_P$ satisfies $Y = T \psi(Y)$, and has the closed form
    $$\psi = \frac{4q(q+2)}{(q+1)^2(2q+1 - 2E_1 T) + \Delta_2 T^2},$$
    where $q$ satisfies the master quintic (see Theorem D).
  - **Theorem D** (Day 152). $\psi$ is algebraic of degree exactly 5 over $\mathbb{Q}(E_1, E_2, E_3)(Y)$, with minimal polynomial $Q(\psi, Y) = \sum_{j=0}^{9} Y^j C_j(\psi, E)$ monic in $\psi$, $E$-homogeneous of degree $j$ in $C_j$, and irreducible.
- **Positioning.** Frame against Jing–Rozhkovskaya vertex-operator normal ordering, Molev–Sagan factorial Schur theory, and Speicher's free cumulant machinery. The novelty is scalar-level, in three variables, and closed-form. None of the three parent frameworks contains an explicit $\psi$.
- **Open problem to advertise at the end:** Conjecture P — $E$-positivity of $[T^n]H$ with minimum coefficient $n+1$ — verified to $n \le 16$, unresolved. The paper explains why the algebraic equation for $\psi$ does not deliver positivity by itself: $Q$ has leading term $-16 E_3^3 Y^9$. The correct handle is filtration-by-defect.
- **Backing files:** `proofs/2026-08-30-day149-H2-PROVED.md`, `proofs/2026-08-30-day148-bk-mod3-SOLVED.md`, `proofs/2026-08-31-day152-psi-closed-form-PROVED.md`, `proofs/2026-08-31-day152b-audit-of-psi-closed-form.md`.

## §2. Object dictionary and conventions (1 pp — but only if needed; may be merged into §1)

- The three normalisation knobs (Day 150b): falling vs rising factorial ($\mathcal T$ vs $\mathcal T^+$), shifted vs plain variables ($u_i$ vs $x_i = u_i + n - i$), ordinary vs shifted Vandermonde. Eight frames; we work in one.
- Fix: $\mathcal T$ takes $u_i^n \to (u_i)_n$ (falling); $x_i = u_i$ (plain); $V = \prod_{i<j}(u_i - u_j)$ (ordinary). Everything else is a translate.
- **Backing file:** `notes/object-dictionary.md` (in WIP repo).
- **NOTE:** if page count is tight, this section drops into a single paragraph and a footnote pointing at the WIP repo dictionary. FPSAC readers know what a factorial Schur function is; we don't need to re-derive it.

## §3. Ψ is Schur → factorial Schur, and its consequences (2 pp)

- **Statement of Theorem A** with the one-line proof: $s_\mu V = \sum_w \operatorname{sgn}(w) \prod_i u_i^{\lambda_{w(i)}}$, and $\mathcal T$ acts monomial-wise, giving $\det[(u_i)_{\lambda_j}]/V = \mathfrak s_\mu$ by Macdonald.
- **Corollary A1:** $\Psi^+(e_2^b) = \sum_{\mu \vdash 2b, \ell(\mu) \le 3} K_{\mu', (2^b)} \mathfrak s_\mu$ (Kostka numbers).
- **Corollary A2:** $\Psi^+(e_3 \cdot f) = E_3 \cdot \tau(\Psi^+(f))$, i.e. $\tau$ acts through $\Psi^+$ as multiplication by $e_3$. Rewriting the whole $\Psi$-recursion this way turns it into the $e_2$-Pieri rule on factorial Schurs, with all structure constants in $\{0, 1\}$.
- **The master curve.** Introduce the Riccati system on $F_P = \Psi^+(e^{T e_2})/V$: this is the three-variable Horn hypergeometric structure, spelled out. The system has a shared master curve which specialises two ways.
- **Proof sketch of Theorem A**: ~4 lines, done. The identifying observation (that $\Psi$ is a known thing) is arguably the entire point of §3.
- **Backing:** Day 149 §5, the Kostka table for $b \le 6$ from `core.py`.

## §4. Downward specialisation: $b_k \equiv 0 \pmod 3$ (2.5 pp)

- **The specialisation.** Set $u_1 = u_2 = u_3 = u$, $T = 1$; only the symmetric part of $F_P$ survives; the diagonal series $F(\vartheta)$ satisfies an explicit algebraic equation
  $$F(F-1)^3 (4F - 3) = \vartheta (2F - 3)^2$$
  of degree 5 in $F$ over $\mathbb{Q}(\vartheta)$ (Day 148 Thm 4.1).
- **The substitution $F = 3G$** kills the leading 3's, giving a monic Lagrange inversion whose kernel $\Phi_A(G) = (2G-1)^2/((3G-1)^3(4G-1))$ has coefficients in $\mathbb{Z}[G]$. Then Lagrange inversion delivers integer coefficients $b_k/3 \in \mathbb{Z}$ directly.
- **The two lemmas used:** (H2) $\deg_{E_3}[T^n]H \le \lfloor n/3 \rfloor$ (Day 149 Thm 2) and a normalisation identity from Day 145 (free-cumulant reduction, Speicher Möbius formula) — both unconditional now.
- **Historical remark (footnote).** This closes a five-session mod-3 arc that started with a Rubine-inspired hypergeometric hunt (dead), a Dwork Frobenius attempt (tautological), and a p-adic split-point analysis (Conjecture C, false). The correct proof needed only a closed form for $F_P$, which nobody looked for.
- **Backing:** `proofs/2026-08-30-day148-bk-mod3-SOLVED.md`.

## §5. Upward specialisation: the Lagrange kernel ψ (3.5 pp)

- **The setup.** $H := \tau(F_P)/F_P$ where $\tau: u_i \mapsto u_i + 1$. Weight grading on $\mathbb{Q}(u_1,u_2,u_3)$ with $\operatorname{wt}(u_i) = 1$; let $\ell_k^{\text{top}}$ be the projection onto weight $k$. Day 149 Thm 1 gives $\operatorname{wt}(\log F_P) \le 1$ (six-line induction on Horn $t$-degree). Set $Y := \int_0^T \ell_0^{\text{top}}(H) dT$ and define $\psi$ by $Y = T \psi(Y)$.
- **Theorem C: closed form.** In the auxiliary variable $q$ satisfying the master quintic:
  $$\psi = \frac{4q(q+2)}{(q+1)^2 (2q + 1 - 2E_1 T) + \Delta_2 T^2}, \qquad \Delta_2 = E_1^2 - 4 E_2.$$
- **The master quintic (Theorem C.2).** Define $R_i := q + 2 T \nu_i$ where $\nu_i$ are the roots of the diagonal Riccati $\nu_i (1 - T(e_1(\nu) - \nu_i)) = u_i$; then $R_i^2 = q^2 + 4 T u_i$, and
  $$\prod_i (R_i^2 - q^2) = 64 T^3 E_3 \Longleftrightarrow (2q^2 - e_2(R) q + e_3(R))(2q^3 + 2q^2 + e_2(R) q + e_3(R)) = 64 T^3 E_3,$$
  with $e_2(R) = -q^2 + 2q + 2 - 2 E_1 T$ and $e_3(R)$ **rational** in $q, T, E$ because the constraint $e_2(u) = E_2$ is linear in it. This is the master quintic in two lines — no radicals, no resultants.
- **Theorem D: minimal polynomial.** Eliminating $q$ against the $\psi$-relation gives $\operatorname{Res}_q = -2048 \psi^9 (3\psi^2 + 2 E_1 Y \psi - \Delta_2 Y^2)^2 (-Q)$ with $Q(\psi, Y)$ monic in $\psi$, $E$-homogeneous coefficient-wise, and irreducible. Irreducibility has a one-line certificate: $Q$ monic in $\psi$, so any factorisation has monic factors, so degree is preserved by specialisation; at $(E_1, E_2, E_3, Y) = (0, 0, 1, 1)$, $Q$ becomes $\psi^5 - \psi^4 + 22 \psi^3 - 89 \psi^2 + 152 \psi - 32$, irreducible mod 5 by distinct-degree factorisation.
- **Two auxiliary results (P1) and (P2).** State each with the six-line proof.
  - (P1) $\log \ell_0^{\text{top}}(H) = \partial \Xi$ where $\Xi = \ell_1^{\text{top}}(\log F_P)$ and $\partial = \sum \partial_{u_i}$.
  - (P2) $\theta \Xi = (P - E_1)/2$ where $\theta = T \, d/dT$ and $P = e_1(\nu)$, via applying $\ell_1^{\text{top}}$ to the multivariate Riccati system.
- **Narayana at $E_3 = 0$.** With $E_3 \to 0$ the top layer of $[T^n]H$ becomes $(n+1) \cdot W_n(E_1, E_2)$ where $W_n$ is the Narayana polynomial; on the further specialisation $E_2 = 0$ this recovers Catalan.
- **Backing:** `proofs/2026-08-31-day152-psi-closed-form-PROVED.md`, `proofs/2026-08-31-day152b-audit-of-psi-closed-form.md`.

## §6. Conjecture P and what remains (1 pp)

- **Statement of Conjecture P (Day 149 §6.2).** $[T^n]H \in \mathbb{Z}_{\ge 0}[E_1, E_2, E_3]$ for all $n$, with minimum coefficient exactly $n+1$ attained at $E_1^n$. Verified to $n \le 16$, zero negatives.
- **Why $Q$ does not deliver it.** $Q$ has leading term $-16 E_3^3 Y^9$, i.e. $\psi$ is not manifestly $E$-positive from its algebraic equation. The positivity has to come from a filtration, not from a single closed form.
- **The filtration.** Stratify $[T^n]H$ by defect $d = n - (a + 2b + 3c)$ over $E_1^a E_2^b E_3^c$. Layer $d = 0$ is the leading symbol, positive by Narayana. Prove Conjecture P layer by layer downward.
- **Two related open problems.** (i) The slice $\psi|_{E_1 = E_2 = 0}$ gives an unidentified sequence $1, 2, 5, 34, 334, 3958, 52599, 755256, \dots$ — not in OEIS, not P-recursive of order $\le 4$, algebraic of degree 5 in $W = E_3 Y^3$. Find a combinatorial model. (ii) $H_1$: $H \in \mathbb{Z}[E][[T]]$. Implied by Conjecture P; verified independently to $T^{16}$.

---

## Length budget

| Section          | Pages | Notes                                                                 |
|------------------|-------|-----------------------------------------------------------------------|
| §1 Intro         | 1.5   | Includes statement of all four theorems.                              |
| §2 Dictionary    | 0.5   | Compressed. Full dictionary is a WIP repo file, cited.                |
| §3 Ψ = Schur→FS  | 2.0   | Includes the master curve introduction.                               |
| §4 b_k mod 3     | 2.5   | Historical remark in a footnote.                                      |
| §5 ψ closed form | 3.5   | The technical centrepiece.                                            |
| §6 Open probs    | 1.0   | Conjecture P + slice + H_1.                                           |
| References       | 1.0   |                                                                       |
| **Total**        | **12**|                                                                       |

## Writing plan (76 days)

- **Weeks 1–2 (Sep 1–14):** §1 draft (statement-level, no proofs), §2 dictionary compression. Get the narrative arc right.
- **Weeks 3–4 (Sep 15–28):** §3 draft with Theorem A proof. This is the shortest and cleanest section; write it first as a warm-up.
- **Weeks 5–6 (Sep 29–Oct 12):** §5 draft. The technical centrepiece. Refactor the Day 152 proof for readability. The two-line master quintic derivation is the aesthetic peak of the paper.
- **Weeks 7–8 (Oct 13–26):** §4 draft. Compress the b_k mod 3 proof to fit; cite Day 148 for full details.
- **Weeks 9 (Oct 27–Nov 2):** §6, references, citation check.
- **Weeks 10 (Nov 3–9):** Send to Clio for review. Six days round-trip if she is fast.
- **Week 11 (Nov 10–15):** Final edits, submit.

## What to send to Clio when

- **End of week 4:** §1 + §3 draft. She has factorial-Schur expertise; the identification Theorem A is her natural reviewing zone.
- **End of week 6:** §5 draft. She has already reviewed the Day 133 material; ψ is new to her.
- **End of week 8:** full draft. This is the honest review request.

## Things that would improve the paper if they landed in time

- Conjecture P proved at layer d=1 (down one level from the top-Narayana result).
- Identification of the slice sequence 1, 2, 5, 34, 334, ...
- A combinatorial interpretation of Q(ψ, Y) — the fact that $\psi^{9}$ appears with sign shift suggests something Lagrange-flavoured.

**None of these are required for the paper.** The paper stands on Theorems A, B, C, D.

## What NOT to include

- **Kerov character polynomials.** Killed Day 151. If T1 revives with a stable-frame lift, we can add a §7; not before.
- **Dwork / p-adic Frobenius.** Tautological. Dead end from Day 147, unrelated to the current arc.
- **The Rubine template.** Never applied.
- **Two-Lagrange-kernels shortcut.** Died Day 150b.
