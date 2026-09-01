# work-in-progress

Rick's stuff, mid-flight.

**Nothing in here is claimed to be correct.** Not the theorems, not the numbers, not the
prose. This repo exists so the work is *visible early*, not so it's *right*. If you read
something here and act on it without checking it yourself, that's on you. A push to this
repo is a push, not a claim. Some of it will turn out to be garbage — that's the whole
point of having a place where garbage is allowed.

Protocol §3 governs what lands here. §3.2 wants source *and* compiled PDF; right now
everything is Markdown with LaTeX inline and there is no `.tex`, so there is nothing to
compile yet. When a paper exists, the PDF ships with it.

## What's actually in here

```
proofs/    the five live write-ups, newest last-modified first
registry/  trust grades — the only place that says how much I believe anything
notes/     object-dictionary.md, aka §1 of the FPSAC paper
```

### `proofs/`

One arc, seven days of it — b_k mod 3 (Days 143–148), then the (H2) cleanup and $\Psi$
identification (Day 149), then $\psi$ closed form (Day 152) and its audit (Day 152b).
The entry point was a sequence $b_k = 3, 27, 417, 7851, 164124, 3661389, 85384566,
2056373739, \dots$ that fell out of an interior series $F_P$ attached to a $\Psi$-recursion
on three variables. Not in OEIS. Not P-recursive at low order. Five sessions to find out
what it was; two more to move up the chain to the leading symbol of $H$.

- **`2026-08-28-day143-invariant-quadratic-identity.md`** — the entry point.
  $a_k = -b_k + \sum_{i+j=k} b_i b_j$, equivalently $(1-2F)^2 = 1+4A$. So $1+4A$ is a
  perfect square in $\mathbb{Q}[[\tau]]$. Verified $k \le 7$, predicted $b_7 = 85384566$
  and the prediction held.
- **`2026-08-29-day145-free-cumulant-integrality.md`** — reduction theorem, proved: for any
  integral moment series, $d \mid m_n$ for all $n$ iff $d \mid \kappa_n$ for all $n$
  (Speicher Möbius inversion, both directions). Applied to $M = 1-2F$ this turns
  "free cumulants land in $6\mathbb{Z}$" into "$b_n \equiv 0 \bmod 3$". The reduction is
  a theorem. What it reduces *to* was, at that point, a conjecture.
- **`2026-08-29-day146-bk-mod3-master-equation.md`** — the master equation
  $L F_P = E_3 T^2[-3 + T(E_1+6+2\theta)]\tau(F_P)$. Reduces $3 \mid b_k$ to integrality of
  $\mathcal{H} = \tau(F_P)/F_P$'s diagonal. **Read the corrections block at the top.** Day 147
  audited this document and killed things in it: a claimed equivalence in §6.3 and §11 was
  false and is corrected in place; §9's Dwork verification was run with $\varsigma = \mathrm{id}$,
  which is not a Frobenius lift, and the conclusion drawn from it is **retracted**. Retracted
  text is marked, not deleted. That's deliberate. You should be able to see what I got wrong.
- **`2026-08-30-day148-bk-mod3-SOLVED.md`** — the theorem. $F$ is algebraic of degree 5:
  $F(F-1)^3(4F-3) = \vartheta(2F-3)^2$. Substituting $F = 3G$ the nines cancel and you get a
  Lagrange inversion with an integral kernel, so $b_k/3 \in \mathbb{Z}$ outright.
  **$b_k \equiv 0 \bmod 3$ for all $k \ge 1$.** The four-session arc collapsed the moment a
  closed form for $F_P$ turned up, and nobody had looked for one. Lesson noted.
- **`2026-08-30-day149-H2-PROVED.md`** — cleanup and identification. (H2) proved twice, the
  second proof stronger ($\deg_u [T^n] H \le n$), so Day 146's Theorem 2 and Day 148's
  corollary are now **unconditional**. Plus the thing I should have known forty days ago:
  $\Psi$ is the map $s_\mu \mapsto \mathfrak{s}_\mu$ from Schur to **factorial Schur**
  functions, and $\tau$ is multiplication by $e_3$ through it.
- **`2026-08-31-day152-psi-closed-form-PROVED.md`** — $\psi$, the Lagrange kernel of the
  leading symbol, is algebraic of degree exactly $5$ over $\mathbb{Q}(E)(Y)$, and its minimal
  polynomial $Q(\psi, Y)$ is now a **theorem**. The two Day-149 §4 statements that Day 151
  had taken on faith — $\log \ell_0^{\rm top}(H) = \partial \Xi$ and $\theta \Xi = (P-E_1)/2$ —
  are both proved in six-line arguments each from the weight grading. Better closed form
  drops out ($\psi = 4q(q+2) / [(q+1)^2(2q+1-2E_1 T) + \Delta_2 T^2]$, regular at $E_3 = 0$,
  no $T^3$ in the denominator); the master quintic falls out of the identity $\prod_i (R_i^2 - q^2)
  = 64 T^3 E_3$ in two lines instead of a resultant chain; irreducibility of $Q$ needs no
  polynomial-factorisation black box. Rule 11 scorecard for the day: unfolding $3$–$0$,
  importing $0$–$9$. No papers opened.
- **`2026-08-31-day152b-audit-of-psi-closed-form.md`** — adversarial audit of the Day 152
  proof, second session same day. **No error found.** Every step re-derived by hand from
  the raw definition $F_P = \mathcal{T}^+(e^{T e_2} V)/V$, including both imports (Day 148
  Lemma 2.1 / Theorem 2.2 / (H) / (R), and Day 149 Theorem 1 (a),(b)). Frame check passed —
  Day 152 §5 and Day 149 §4 are the same $\psi$, same normalisation. New symbolic-in-$E$
  pipeline (strictly stronger than Day 152's `check5.py`, which was numeric only) verifies
  $Q(\psi, Y) = 0$ through $Y^{12}$ and reproduces the pre-registered $[Y^9]\psi$. Theorem D
  now depends on no polynomial-factorisation algorithm at all.

### What is still open

1. **(H1)**: $H = \tau(F_P)/F_P \in \mathbb{Z}[E_1,E_2,E_3][[T]]$. Verified to $T^{16}$.
   Not proved. It's now a completely explicit question about Kostka numbers and factorial
   Schur functions — no $p$-adics anywhere — and the entire difficulty is one $b!$.
2. **Conjecture P**: $[T^n]H$ is coefficientwise non-negative, minimum $n+1$. Verified to
   $n \le 16$, zero negatives. Strictly stronger than (H1) and implies it. The mechanism has
   to be cancellation *between* factorial Schur functions, since the individual ones aren't
   positive. Find the basis where it's obvious. The $\psi$ closed form does **not** deliver
   $E$-positivity by itself: $Q$ has leading term $-16 E_3^3 Y^9$, not manifestly positive.
3. **The slice $1,2,5,34,334,3958,\dots$**: this is $[Y^n]\psi|_{E_1=E_2=0} / E_3^{\lceil ? \rceil}$
   (see the Day 152 proof for the exact statement); it satisfies an explicit quintic in
   $W = E_3 Y^3$, is **not** in OEIS, is **not** P-recursive of order $\le 4$, and is not
   $\Phi_{\rm A}$. Wanted: identify it, or prove it can't be identified with a standard
   family. The pre-registered Catalan prediction from Day 150b died at coefficient 4.

Dead ends worth knowing about, so nobody wastes a week: Conjecture C ($F_P$ $\ell$-integral at
separable $E$) is **FALSE**, 38108 violations out of 14835 tested pairs. The $J$-fraction of
$\sum H_n T^n$ is neither integral nor positive. The Dwork route is tautological.

### `registry/`

**This is where the actual trust grades live. Believe the registry, not the prose.**

Grades are ordered: `hunch < sketched < computed < checked-sober < proved < lean-verified`,
plus `dead-end` (abandoned, must carry a reason), `in-progress`, and `unclassified` (on disk,
never re-checked). `hunch` means the 2am itch with no evidence. `computed` means the machine
agreed on every case I tried — evidence, not proof. `checked-sober` means I re-derived it cold
and it survived, and it *must* carry a `recheck` field with the date and path or the validator
errors. `proved` *must* point at the written proof file. `lean-verified` names a sorry-free
declaration.

The boundary rule: a node may claim `checked-sober` or above only if every non-dead-end child
is at least `checked-sober`. **A proof standing on a hunch is a hunch with paperwork.**
Dead ends carry their own evidence level too — a `reason` is a claim, and a wrong one silently
prunes a live branch, which is the most expensive error in a search.

`registry/README.md` has the full format. Four JSON trees currently:
`beta-prime-mod8.json` (116 nodes) and `strict-axis-closed-form.json` (10 nodes) are from
**earlier arcs**, both still `in-progress`. `bk-mod3.json` (24 nodes) is the b_k arc,
status `proved` since Day 148. `conjecture-P.json` (29 nodes) is the live one covering the
$\psi$/master-curve/Conjecture-P arc; status `in-progress`. As of Day 152 the sub-node
`psi-closed-form-degree5` is `proved`, with an independent audit at `checked-sober`;
`kerov-character-polynomial-bridge` and `psi-catalan-prediction` are `dead-end`.

### `notes/object-dictionary.md`

Every symbol, its ring, its normalisation. Written because the same failure mode fired about
ten times — "two names, one object" or "one name, two objects" — and every single firing was
a disagreement about one of exactly three binary normalisation knobs. Eight frames; the project
walked between them for forty days without a map. This is §1 of the FPSAC paper, so it isn't
overhead. If you're reading the proofs and two formulas look the same "up to normalisation,"
come here and name which knob differs before you write anything down.
