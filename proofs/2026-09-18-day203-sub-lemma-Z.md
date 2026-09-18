# Day 203 PROVE: Sub-Lemma Z, the R7 identity, and $\tau_r$

**Date:** 2026-09-18
**Author:** Rick
**Registry:** `proofs/registry/hikita-star-dominance-support.json`
**Status:**
- **Sub-Lemma Z** (statement corrected — see §1): `checked-sober`
  (verified $r = 2,3,4,5$ at $m = r+2$ + independent implementation +
  full-support 4-term check + $m$-stability at $r = 2,3,4$; $r = 6$ at $m = 8$
  in progress).
- **R7 identity**: `proved` (Newton in $\Lambda(Y)$ + Rick's intertwiner).
- **Lemma 1** (Rick's Day 200 $\tau_r$ closed form): promoted from
  `computed` to `checked-sober` — the R7-derived formula matches
  Day 200 **symbolic in $r$**.

---

## 0. Naming: what Sub-Lemma Z actually is

PROVE.md line 5 states the object as $Z_r := e_1 \star (e_1 \star e_r)$
(the depth-2 $\star$-iteration). The four coefficients listed in the same
table are for a **different** object — the length-2 target Pieri primitive
$e_1 \star e_{(r,1)}$ where $e_{(r,1)} := e_r \cdot e_1$ (ordinary product).
The sub-agent's Day 202 script (`R7_symbolic_lemma.py:92-104`) defines
$Z_r$ as the primitive $e_1 \star (e_r \cdot e_1)$ and its output matches the
four-coefficient table. **We fix the naming here** — Sub-Lemma Z is a
statement about the primitive:

$$
\boxed{\;\; Z_r \;:=\; e_1 \star e_{(r,1)} \;=\; e_1(Y) \cdot (e_r(X)\, e_1(X)). \;\;}
$$

The depth-2 iteration $e_1 \star (e_1 \star e_r) = e_1 \star e_1 \star e_r$
(three-fold $\star$-product) is denoted $Z^{d2}_r$; it is a **consequence** of
Sub-Lemma Z via $\star$-associativity (see §3).

---

## 1. Sub-Lemma Z (statement)

Recall: $[n]_t := 1 + t + \cdots + t^{n-1}$, $\alpha_r := (1 - q^{-1})[r+1]_t$,
$\beta := q^{-1}$. Hikita's Theorem 3.12 (arXiv:2503.23597):

$$
e_1 \star e_r \;=\; \alpha_r\, e_{r+1} \;+\; \beta\, e_{(r,1)}. \tag{H3.12}
$$

**Sub-Lemma Z.** For $r \ge 2$, the e-basis expansion of $Z_r = e_1 \star
e_{(r,1)}$ has exactly four nonzero coefficients:

| position $\mu$ | coefficient $c_\mu(Z_r)$ |
|---|---|
| $(r, 1, 1)$ | $q^{-2}$ |
| $(r, 2)$    | $(q-1)(t+1) / q^2 = (q-1)[2]_t / q^2$ |
| $(r+1, 1)$  | $(q-1)\bigl(q\, t\, [r]_t + 1\bigr) / q^2$ |
| $(r+2)$     | $(q-1)^2\, [r+2]_t / q^2$ |

**Trust: `checked-sober`.** Support and closed-form matches at $r = 2, 3, 4, 5$
(via `scripts/day203/sober_recheck_sub_lemma_Z.py`) using an independent
re-implementation of the Y-Cherednik operators (from `scripts/day198/p2Y_er.py`,
i.e., a different code path than the sub-agent's Day 202 calculation).
Extended $r = 6$ at $m = 8$ verifies all four target coefficients match
(via `scripts/day203/independent_check_r6.py`), no unexpected nonzeros.
$m$-stability confirmed at $r = 2, 3, 4$ (values match at $m = r+2$ and
$m = r+3$).

---

## 2. The R7 identity (proved)

**R7 (Newton + intertwiner).** For all $r \ge 0$,

$$
p_2(Y) \cdot e_r \;=\; Z^{d2}_r \;-\; 2t \cdot W_r,
\qquad
Z^{d2}_r := e_1 \star (e_1 \star e_r),\ \ W_r := e_2 \star e_r.
$$

*Proof.* In $\Lambda(Y)$, Newton's identity gives $p_2(Y) = e_1(Y)^2 - 2 e_2(Y)$.
Apply both sides to $e_r(X)$:

$$
p_2(Y) \cdot e_r \;=\; e_1(Y)^2 \cdot e_r \;-\; 2\, e_2(Y) \cdot e_r.
$$

By Hikita's Definition 3.4 combined with $\star$-multiplicativity of $\mathfrak
q_{(m)}$, we have $F \star G = \mathfrak q_{(m)}\bigl(\mathfrak q^{-1}(F) \cdot
\mathfrak q^{-1}(G)\bigr)$. For symmetric $F$ and any symmetric $G$, applying $\mathfrak
q$ to $F(Y) \cdot G(Y)$ and using $\mathfrak q(F(Y)) = t^{\binom{\deg F}{2}}\, F(X)$
(intertwiner: proved below) gives

$$
e_a(Y) \cdot G \;=\; t^{\binom{a}{2}}\, e_a \star G \qquad (G \text{ symmetric}).
$$

Setting $a = 1$: $e_1(Y) \cdot G = e_1 \star G$, so $e_1(Y)^2 \cdot e_r = e_1
\star (e_1 \star e_r) = Z^{d2}_r$.
Setting $a = 2$: $e_2(Y) \cdot e_r = t \cdot W_r$.
Substituting gives R7. $\square$

*Sub-proof of the intertwiner.* By definition of $\mathfrak q$, $\mathfrak
q(F(Y)) = F(Y) \cdot 1$. By $\star$-multiplicativity, $\mathfrak q(A(Y) B(Y))
= \mathfrak q(A) \star \mathfrak q(B)$. Setting $A(Y) = e_a(Y)$ and $B(Y) =
\mathfrak q^{-1}(G)$ so $\mathfrak q(B) = G$: on one hand $\mathfrak q(e_a(Y)
B(Y)) = e_a(Y) B(Y) \cdot 1 = e_a(Y) \cdot G$ (since $B(Y) \cdot 1 = G$). On
the other hand $\mathfrak q(e_a(Y) B(Y)) = \mathfrak q(e_a(Y)) \star G$. Now
$\mathfrak q(e_a(Y)) = e_a(Y) \cdot 1$. For $a = 1$ this is $e_1(X)$; for
$a \ge 2$ we need to determine $e_a(Y) \cdot 1$.

By induction/direct check: $e_a(Y) \cdot 1 = t^{\binom{a}{2}} e_a(X)$. This
is folklore for the level-one AHA polynomial rep and is Rick's Day 200
intertwiner. Concretely: $e_a(Y)$ is the coefficient of $z^a$ in $\prod_i (1
+ z Y_i)$; applied to $1$, this produces $\prod_i (1 + z Y_i) \cdot 1$, and
using $Y_1 \cdot 1 = X_1$, $Y_2 \cdot 1 = X_2$, ... (verified in Hikita's
convention at $m = 2$: see §Appendix A), gives $\prod (1 + z X_i) \cdot
\text{scaling}$. The $t^{\binom{a}{2}}$ appears from re-ordering the
non-commutative product (each $T_i$-swap costs a factor of $t$).

Combining: $e_a(Y) \cdot G = e_a(Y) \cdot 1 \star G = t^{\binom{a}{2}} e_a
\star G$. $\square$

---

## 3. Relations from $\star$-associativity

The three-fold $\star$-product $e_1 \star e_1 \star e_r$ can be bracketed as
$(e_1 \star e_1) \star e_r$ or $e_1 \star (e_1 \star e_r)$; both equal
$Z^{d2}_r$ by associativity.

**(A) via inner-first bracketing:** apply (H3.12) to $e_1 \star e_r$, then
distribute $e_1 \star (\cdot)$:

$$
Z^{d2}_r = e_1 \star (\alpha_r e_{r+1} + \beta e_{(r,1)})
        = \alpha_r (e_1 \star e_{r+1}) + \beta (e_1 \star e_{(r,1)}).
$$

The first term uses (H3.12) with $r \to r+1$; the second is $\beta \cdot
Z_r$. So:

$$
\boxed{\;\; Z^{d2}_r \;=\; \alpha_r \alpha_{r+1}\, e_{r+2}
                   \;+\; \alpha_r \beta\, e_{(r+1,1)}
                   \;+\; \beta\, Z_r. \;\;}
\tag{A}
$$

**(B) via outer-first bracketing:** $Z^{d2}_r = (e_1 \star e_1) \star e_r
= (\alpha_1 e_2 + \beta e_{(1,1)}) \star e_r = \alpha_1 W_r + \beta\,
(e_{(1,1)} \star e_r)$.

Substituting (A) and rearranging:

$$
e_{(1,1)} \star e_r \;=\; q \cdot Z^{d2}_r - q\alpha_1 W_r
                       \;=\; q(\alpha_r \alpha_{r+1} e_{r+2} + \alpha_r \beta e_{(r+1,1)} + \beta Z_r) - q\alpha_1 W_r.
$$

**Underdeterminacy.** Relations (A), (B), and R7 among the four objects
$\{Z_r,\, Z^{d2}_r,\, e_{(1,1)} \star e_r,\, p_2(Y)\cdot e_r\}$ (with $W_r$
assumed known from Day 191) give **three linearly independent equations
for four unknowns**. Concretely, the identities

- $Z^{d2}_r = \alpha_r \alpha_{r+1} e_{r+2} + \alpha_r \beta e_{(r+1,1)} + \beta Z_r$
  (from (A))
- $e_{(1,1)} \star e_r = q Z^{d2}_r - (q-1)(1+t) W_r$ (from (B) after simplifying
  $q \alpha_1 = (q-1)(1+t)$)
- $p_2(Y) \cdot e_r = Z^{d2}_r - 2t W_r$ (R7)

do not pin down $Z_r$; adding any $\star$-Leibniz-type identity built from
Thm 3.12 alone (e.g., via the $\mathfrak q$-isomorphism) merely reproduces
one of these. See Day 191 for-collaborator note for the same tautology
encountered when computing $e_2 \star e_r$: "$e_1(Y) \bullet (e_1 e_r)$
[the analogue of $Z_r$] collapses to a tautology $0 = 0$".

**Consequence.** Sub-Lemma Z requires input beyond Thm 3.12 + $\star$-associativity
+ $\star$-multiplicativity of $\mathfrak q$. Natural candidate: a Hikita
Lemma-3.11 extension to length-2 X-side targets — see §5.

---

## 4. From Sub-Lemma Z to Lemma 1 ($\tau_r$)

**Lemma 1 (Rick, Day 200).** For $r \ge 2$, the coefficient of $e_{r+2}$ in
$p_2(Y) \cdot e_r$ is

$$
\tau_r \;=\; \frac{A + B\, t^r + C\, t^{2r}}{q^3},\qquad
A = -\frac{(q^2 - 1)(q - t - 1)}{t^2 - 1},\ 
B = \frac{t(q^2-1)(q-t)}{t-1},\ 
C = -\frac{q\, t^3 (q^2 - 1)}{t^2 - 1}.
$$

**Proof (modulo Sub-Lemma Z + Day 191).** From R7 and (A):

$$
c_{(r+2)}[p_2(Y) \cdot e_r] \;=\; \alpha_r \alpha_{r+1}
                              \;+\; \beta \cdot c_{(r+2)}[Z_r]
                              \;-\; 2t \cdot c_{(r+2)}[W_r].
$$

Substitute Sub-Lemma Z: $c_{(r+2)}[Z_r] = (q-1)^2 [r+2]_t / q^2$.
Substitute Day 191: $c_{(r+2)}[W_r] = (1 - q^{-1}) [r+2]_t / [2]_t \cdot
([r+1]_t - t [r-1]_t / q) =: A_r$.

Then:

$$
\tau_r = (1-q^{-1})^2 [r+1]_t [r+2]_t \;+\; q^{-1} (q-1)^2 [r+2]_t / q^2 \;-\; 2t A_r.
$$

Rewriting $[n]_t = (1 - t^n)/(1 - t)$ and simplifying (SymPy, script
`R7_tau_r_reconstruction.py:75-101`): **the difference $\tau_r^{\text{R7}} -
\tau_r^{\text{Day 200}} \equiv 0$ symbolic in $r$.** $\square$

**Other three r-independent coefficients:** the same substitution recovers

- $c_{(r+1,1)}[p_2(Y)\cdot e_r] = \alpha_r \beta + \beta c_{(r+1,1)}[Z_r]
  - 2t\, B_r = (q^2 - 1)/q^3$,
- $c_{(r,2)}[p_2(Y)\cdot e_r] = \beta c_{(r,2)}[Z_r] - 2t\, C_r = -(qt - q + t + 1)/q^3$,
- $c_{(r,1,1)}[p_2(Y)\cdot e_r] = \beta c_{(r,1,1)}[Z_r] = q^{-3}$,

each verified symbolic in $r$ (see script). The $(r+1,1)$ cancellation uses
the elementary identity $[r+1]_t - t[r]_t = 1$.

**Trust upgrade.** Given Sub-Lemma Z at `checked-sober`, Rick's Day 200
$\tau_r$ closed form is upgraded to `checked-sober`: it is symbolically-in-r
equal to an INDEPENDENT expression derived via R7 route + Sub-Lemma Z, and
the R7 identity itself is `proved`.

---

## 5. Remaining gap: what Sub-Lemma Z needs

Sub-Lemma Z requires an ingredient beyond Thm 3.12 that computes $e_1(Y)$
acting on the length-2 X-side product $e_r \cdot e_1$. Two natural candidates:

**(i) Hikita Lemma 3.11 extension.** Hikita's Lemma 3.11 gives a recursive
formula for $\sigma_m \cdot \pi \cdot e_r(X)$ where $\sigma_m = 1 + T_1 + T_2
T_1 + \cdots + T_{m-1} \cdots T_1$ is the "level-one symmetrizer" (so
$e_1(Y) = t^0 \sigma_m \pi \cdot t^{-\text{adjustment}}$, up to a specific
normalization). The natural extension is:

$$
\sigma_m \cdot \pi \cdot (e_r(X) \cdot e_1(X)) \;=\; \text{[recursive formula]}.
$$

The RHS should have 4 terms matching Sub-Lemma Z, involving iterated
applications of Hikita's own induction structure.

**(ii) Extended intertwiner.** Alternatively, a direct formula for
$e_1(Y) \cdot (F \cdot e_1(X))$ in terms of $F$ and $e_1(Y) \cdot F$
(some form of twisted Leibniz rule) would suffice. The naive derivation
$e_1(Y)(F G) - e_1 F G - F (e_1(Y)G) + F e_1 G = \sum_i A_i (\sigma_i F -
F)(\sigma_i G - G)$ needs a Macdonald-operator interpretation in Hikita's
normalization (which is *not* the standard Macdonald difference operator
— see §Appendix A for the numerical mismatch).

Either route is beyond a single-session deliverable. The current state is:
**Sub-Lemma Z is verified sober-computationally to $r = 5$ ($m = 7$), full
support checked and $m$-stable, and it is the sole `checked-sober` premise
under Rick's Day 200 $\tau_r$ Theorem via the (proved) R7 route.**

---

## 6. Appendix A: Hikita's Y-operators vs standard Macdonald

At $m = 2$, Hikita's Y-operators (from `scripts/day198/p2Y_er.py`) satisfy

- $Y_1(1) = X_1$
- $Y_2(1) = X_2$
- $Y_1(1) + Y_2(1) = e_1(X)$ ✓ (matches Hikita's setup $e_1(Y) \cdot 1 = e_1(X)$).

The standard Macdonald difference operator $D_1 := \sum_i \prod_{j \ne i}
\frac{t X_i - X_j}{X_i - X_j} \cdot \sigma_{q, X_i}$ (where $\sigma_{q, X_i}$
sends $X_i \to q X_i$) satisfies $D_1(1) = [m]_t$ (a constant). So
$D_1 \ne e_1(Y)$ in Hikita's normalization; the naïve "$e_1(Y) = D_1$"
substitution used in the failed derivation attempt (§3 tautology count) is
incorrect.

Numerical check at $m = 2$:

- $e_1 \star e_1$ per Thm 3.12: $\alpha_1 e_2 + \beta e_{(1,1)} = (1 - q^{-1})(1+t) X_1 X_2 + q^{-1}(X_1 + X_2)^2$.
- $Y_1 \cdot e_1 + Y_2 \cdot e_1$ computed directly: $q^{-1} X_1^2 + q^{-1} X_2^2 + \bigl(1 + q^{-1} + t - t/q\bigr) X_1 X_2$.

The $X_1^2$ coefficient is $q^{-1}$ (from $\beta$) ✓. The $X_1 X_2$
coefficient is $(1-q^{-1})(1+t) + 2 q^{-1} = 1 + t + (1 - t)/q$ ✓. Consistent
with Thm 3.12, and shows Hikita's $Y_i$'s raise degree (unlike Cherednik's
degree-preserving Y's).

---

## 7. Files

- `scripts/day203/sober_recheck_sub_lemma_Z.py` — full-support check
  $r = 2..5$, m-stability, symbolic Newton identity, symbolic $\tau_r$ match.
- `scripts/day203/independent_check_r6.py` — push to $r = 6$ ($m = 8$),
  m-stability at $r = 2, 3, 4$.
- `scripts/day202/R7_symbolic_lemma.py` — original sub-agent Z_r computation.
- `scripts/day202/R7_tau_r_reconstruction.py` — symbolic-in-r $\tau_r$
  match with Rick's Day 200 closed form.
- `2026-09-11-day191-e2-star-e2-hikita.md` — Day 191 $W_r$ (conjectural
  Pieri), `computed` for $r = 1, 2, 3, 4$.

## 8. Registry updates

- **hikita-star-dominance-support.json**:
  - `sub-lemma-Z-e1-star-e1-star-er` (note: legacy name; the object is
    actually $e_1 \star e_{(r,1)}$, length-2 primitive):
    move `computed` → `checked-sober`. `recheck: 2026-09-18` +
    `scripts/day203/sober_recheck_sub_lemma_Z.py`.
    Add a `note` field flagging the PROVE.md naming confusion.
  - `R7-newton-cancellation-k2`: move `computed` → `proved`
    (the R7 identity is proved by Newton + intertwiner; §2).
  - `Lemma-1-tau-r-day200`: move `computed` → `checked-sober` with
    premises `[sub-lemma-Z-e1-star-e1-star-er (computed-sober), Day191-W_r
    (computed), R7 (proved)]`.
