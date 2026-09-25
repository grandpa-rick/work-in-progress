#!/usr/bin/env python3
"""
Referee check for MacBeth, "Tool-calling protocols as directed containers" (wip 6a66fe3).

Question: do comonad distributive laws  lam : Store_S Store_T => Store_T Store_S
exist for the bare store comonads (Delta S, Delta T of the note, Def. 5), with
|S| = |T| = 2?  The note (Sec. 5) says "without a valid lambda there is no
composite protocol at all" and that at the smallest re-entrant model a weld exists
iff eps = 0.  If distributive laws between the bare stores exist unconditionally,
then eps must encode extra data (prescribed mutual actions) that the note never
specifies.

Natural transformations between container extensions = container morphisms, so
lam is determined by: for each source shape (s, g : S->T) a target shape
(t, h : T->S) and a backward position map rho : T x S -> S x T.  All four
distributive-law axioms are checked on generic (symbolic) elements, which suffices
by naturality.
"""
import itertools

S = (0, 1); T = (0, 1)

def eps(v):            # Store counit: (a, f) -> f[a]
    a, f = v; return f[a]
def delta(v):          # Store comultiplication
    a, f = v; return (a, tuple((a2, f) for a2 in range(len(f))))
def fmap(phi, v):
    a, f = v; return (a, tuple(phi(x) for x in f))

src_shapes = [(s, g) for s in S for g in itertools.product(T, repeat=len(S))]
tgt_shapes = [(t, h) for t in T for h in itertools.product(S, repeat=len(T))]
tpos = [(t2, s2) for t2 in T for s2 in S]
spos = [(s2, t2) for s2 in S for t2 in T]

def apply(lam, v):
    """lam: dict src_shape -> (tgt_shape, rho dict tpos->spos).  v in Store_S Store_T Y."""
    s, F = v
    g = tuple(F[s2][0] for s2 in S)
    (t, h), rho = lam[(s, g)]
    val = lambda sp: F[sp[0]][1][sp[1]]
    Psi = tuple((h[t2], tuple(val(rho[(t2, s2)]) for s2 in S)) for t2 in T)
    return (t, Psi)

def generic(shape):
    s, g = shape
    return (s, tuple((g[s2], tuple(('x', s2, t2) for t2 in T)) for s2 in S))

def ax_counits(lam, sh):
    v = generic(sh); out = apply(lam, v)
    a = fmap(eps, out) == eps(v)                       # H eps^G . lam = eps^G H
    b = eps(out) == fmap(eps, v)                       # eps^H G . lam = G eps^H
    return a and b

def ax_comult(lam, sh):
    v = generic(sh)
    # (c)  H delta^G . lam = lam_G . G lam . delta^G_H
    lhs = fmap(delta, apply(lam, v))
    rhs = apply(lam, fmap(lambda w: apply(lam, w), delta(v)))
    if lhs != rhs: return False
    # (d)  delta^H_G . lam = H lam . lam_H . G delta^H
    lhs = delta(apply(lam, v))
    rhs = fmap(lambda w: apply(lam, w), apply(lam, fmap(delta, v)))
    return lhs == rhs

# per-shape candidates satisfying the counit axioms
cands = {}
for sh in src_shapes:
    cs = []
    for tsh in tgt_shapes:
        for img in itertools.product(spos, repeat=len(tpos)):
            rho = dict(zip(tpos, img))
            lam = {sh: (tsh, rho)}
            try:
                if ax_counits(lam, sh): cs.append((tsh, rho))
            except KeyError:
                pass
    cands[sh] = cs
print("counit-compatible candidates per source shape:", [len(cands[sh]) for sh in src_shapes])

laws = []
for choice in itertools.product(*[cands[sh] for sh in src_shapes]):
    lam = dict(zip(src_shapes, choice))
    if all(ax_comult(lam, sh) for sh in src_shapes):
        laws.append(lam)
print("distributive laws Store_2 Store_2 => Store_2 Store_2:", len(laws))
for lam in laws:
    print("  shape map:", {sh: lam[sh][0] for sh in src_shapes})

# Closed form of one law (T = Z/2 viewed as a group):
#   lam(s, F) = (g(s), t' |-> (s, s' |-> f_{s'}(g(s') + t' - g(s))))
# i.e. shape (s,g) |-> (g(s), const_s), rho(t', s') = (s', g(s') + t' - g(s)).
closed = {}
for sh in src_shapes:
    s, g = sh
    closed[sh] = ((g[s], tuple(s for _ in T)),
                  {(t2, s2): (s2, (g[s2] + t2 - g[s]) % len(T)) for t2 in T for s2 in S})
print("closed-form 'relative displacement' law satisfies all four axioms:",
      all(ax_counits(closed, sh) and ax_comult(closed, sh) for sh in src_shapes),
      "| is among the laws found:", closed in laws)
