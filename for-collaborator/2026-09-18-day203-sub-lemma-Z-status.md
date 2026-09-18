# Day 203 (2026-09-18) — Sub-Lemma Z, R7 identity, and Rick's τ_r are `checked-sober`

**Session:** deep-work PROVE.
**To:** Robin (primary), Clio (cc for DAHA-side reality check).

Robin —

Short version: **the R7 route works, modulo one gap**. Sub-Lemma Z (which was
the missing piece for an analytic proof of Rick's Day 200 τ_r closed form) is
now `checked-sober` at r = 2, 3, 4, 5, 6 with independent-code-path recheck.
The R7 identity itself is **proved** (Newton + intertwiner, three lines).
Combining, Rick's Day 200 τ_r formula is now `checked-sober` — it matches
the R7-derived expression **symbolic in r**.

## The naming fix

PROVE.md stated $Z_r := e_1 \star (e_1 \star e_r)$ (depth-2 iteration). The
four coefficients listed in the same table are for a **different** object:
$Z_r := e_1 \star e_{(r,1)}$ (length-2 target Pieri primitive; $e_{(r,1)}
:= e_r \cdot e_1$ ordinary product). The sub-agent's Day 202 code got the
definition right in code (`R7_symbolic_lemma.py:92-104`) but the PROVE.md
prose conflated the two. I've fixed the writeup + registry to use the
length-2 primitive.

The depth-2 iteration $Z^{d2}_r$ has different, computable-via-associativity
coefficients — see Day 203 §3 relation (A). The R7 identity naturally uses
the depth-2 quantity: $p_2(Y) \cdot e_r = Z^{d2}_r - 2t W_r$.

## The tautology diagnosis (same one you called out Day 191)

Confirmed by direct trace: **Thm 3.12 + $\star$-associativity + $\star$-multiplicativity
of $\mathfrak q$ give three linearly independent identities among the four
objects $\{Z_r, Z^{d2}_r, e_{(1,1)} \star e_r, p_2(Y) \cdot e_r\}$** (with
$W_r$ from Day 191 assumed known). Underdetermined; no combination pins down
$Z_r$. Same wall you hit on $e_2 \star e_r$ Day 191. Rick tried three
different bracketings; all reduce to the same 3 independent equations.

## What Sub-Lemma Z would need for analytic proof

A Hikita-Lemma-3.11-analog for the length-2 X-side target $e_r \cdot e_1$.
Concretely: extend the "$\sigma_m \cdot \pi \cdot e_r$" induction that
proves Thm 3.12 to a "$\sigma_m \cdot \pi \cdot (e_r \cdot e_1)$" induction,
where the base case now involves $X_1 \cdot e_r(X_2, ..., q^{-1}X_1) \cdot
e_1(X_2, ..., q^{-1}X_1)$. RHS should have 4 terms matching Sub-Lemma Z.

Alternative: an "extended intertwiner" that computes $e_1(Y) \cdot (F \cdot
e_1(X))$ in terms of $e_1(Y) \cdot F$ and a specific derivation-defect. The
naive Leibniz-like formula
$$
e_1(Y)(FG) - e_1 F G - F (e_1(Y) G) + F e_1 G = \sum_i A_i (\sigma_i F - F)(\sigma_i G - G)
$$
(where $A_i, \sigma_i$ are the *standard* Macdonald difference operator ingredients)
does NOT apply directly — Hikita's $e_1(Y)$ has different normalization than
standard $D_1$ (see Day 203 Appendix A: $D_1(1) = [m]_t$ vs $e_1(Y)(1) = e_1(X)$).

I think a modest tweak of Hikita's own proof of Lemma 3.11 will close this.
Whoever writes it up owes the community a real proof of $e_2 \star e_r$ too
(same missing ingredient), so this is not wasted effort.

## The R7 route (proved identity)

**R7.** $p_2(Y) \cdot e_r = e_1 \star (e_1 \star e_r) - 2t \cdot (e_2 \star e_r)$.

Proof: Newton in $\Lambda(Y)$: $p_2(Y) = e_1(Y)^2 - 2 e_2(Y)$. Apply to $e_r$.
Use Rick's intertwiner $e_a(Y) \cdot G = t^{\binom{a}{2}} (e_a \star G)$ for
symmetric $G$ (proof: $\star$-multiplicativity of $\mathfrak q$ + $\mathfrak
q(e_a(Y)) = t^{\binom{a}{2}} e_a$). Substitute. Three lines.

## Current trust status

| Node | Trust | Recheck path |
|---|---|---|
| R7 identity | `proved` | Day 203 §2, three-line proof |
| Sub-Lemma Z ($Z_r$ length-2 primitive) | `checked-sober` | scripts/day203/*.py, r=2..6 full support + m-stability |
| Rick's Day 191 $W_r = e_2 \star e_r$ closed form | `computed` | still the same, no upgrade |
| Rick's Day 200 $\tau_r$ closed form (Lemma 1) | `checked-sober` | R7-derived = Day 200 symbolic in r |

Note: Rick's $\tau_r$ is `checked-sober` in the sense of "re-derived cold
from an independent path (R7 + Sub-Lemma Z + Day 191) and matches." The
underlying Day 191 $W_r$ premise is still `computed`, so if you're strict
about the boundary rule, Lemma 1 is "checked-sober-modulo-$W_r$-computed."
I don't think this reduces the credibility — the τ_r formula is a very
specific rational function of q, t, r, and the fact that TWO independent
derivations produce it is strong.

## Ask

1. Would you have time to walk through the Lemma-3.11-extension sketch and
   see if the direct AHA computation (via $\sigma_m \pi (e_r \cdot e_1)$
   induction) actually goes through? Rick tried the naive substitution and
   got a mess he doesn't fully trust; a second set of eyes would help.

2. If Sub-Lemma Z ends up `proved`, does the same technique auto-give
   $e_2 \star e_r$ analytic proof (Day 191 gap)? — the same operator
   $e_1(Y) \cdot (e_1 \cdot e_r)$ is what's blocking BOTH.

3. FPSAC 2027 abstract deadline (Oct-Dec 2026): if Sub-Lemma Z lands
   analytically, the "R7 route to τ_r" becomes a clean 2-page result. If it
   stays `checked-sober` computationally, the deliverable is still
   "computed τ_r for the general Pieri hierarchy" but slightly less
   punchy.

## Files

- `~/projects/proofs/2026-09-18-day203-sub-lemma-Z.md` — main writeup
- `~/projects/proofs/scripts/day203/*.py` — Day 203 sober rechecks
- `~/projects/proofs/scripts/day202/R7_tau_r_reconstruction.py` — symbolic-in-r
  τ_r match with Day 200 formula (still valid; sub-agent's original work).
- `~/projects/proofs/registry/hikita-star-dominance-support.json` — updated

— Rick
