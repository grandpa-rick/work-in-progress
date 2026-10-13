"""Day 238 (2026-10-13): checks for MacBeth UID 378 (nerve A / uniqueness B, e4d1785 / stamp fed2448)
and UID 379 (crown Stage 2, e4d1785; p.1 still prints 4e7264d).  Rick.

  Part 0  re-run the 2026-10-11 r2 scripts (reviews/scripts/2026-10-11-{crown,nerve}-r2/) and diff each
          fresh log against the committed log.  The printed theorem statements those scripts encode did
          not change (statement diff in Part 3), so identical logs = the r2 verdicts still apply.
  Part 1  crown Rem. 20 / Sec. 5 intro (l.421-423, l.497-498): "the functors extending Psi with no condition
          are exactly the ones that act trivially on automorphisms".  Counterexample: a functor extending
          Psi on objects (skeleton <=5 nodes, both variances) that acts NON-trivially (and non-canonically)
          on Aut(3-leaf star) through the sign character.  Negative controls that must FAIL:
          (i) canonical action on the star (no fixed branch for the constants), (ii) sign action with the
          constant parked on a moved branch.
  Part 2  crown l.483-486: "read backward, R_f is a left-total injective partial map br(T') -> br(T)".
          As a partial map b' |-> b: count morphisms where it is non-total and where it is non-injective.
  Part 3  statement drift: sentence-normalised diff of pdftotext of old vs new PDFs, restricted to
          Theorem/Lemma/Proposition/Corollary/Definition/Example blocks.
"""
import itertools, os, re, subprocess, sys, difflib

WIP = "/home/agent/projects/work-in-progress"
PEER = "/home/agent/projects/peers/macbeth/proofs"

# ------------------------------------------------------------------ Part 0
print("== Part 0: re-run 2026-10-11 r2 scripts, compare with committed logs")
for d, s in [("2026-10-11-crown-r2", "crown_r2_checks"), ("2026-10-11-nerve-r2", "check_N1_and_levels"),
             ("2026-10-11-nerve-r2", "check_magma_and_splitting"), ("2026-10-11-nerve-r2", "check_noteB_rank1_lift")]:
    p = os.path.join(WIP, "reviews/scripts", d)
    out = subprocess.run([sys.executable, s + ".py"], cwd=p, capture_output=True, text=True, timeout=1200).stdout
    old = open(os.path.join(p, s + ".log")).read()
    same = out.strip() == old.strip()
    print(f"   {s}: fresh log {'IDENTICAL' if same else 'DIFFERS'} to committed log ({len(out.splitlines())} lines)")
    if not same:
        for l in difflib.unified_diff(old.splitlines(), out.splitlines(), lineterm="", n=0):
            print("      ", l)

# ------------------------------------------------------------------ shared tree machinery (as in crown_r2_checks)
def trees(n):
    if n == 1:
        yield (None,); return
    for par in itertools.product(*[range(i) for i in range(1, n)]):
        yield (None,) + par

def anc(par, v):
    out = [v]
    while par[v] is not None:
        v = par[v]; out.append(v)
    return out

def leaves(par):
    kids = {p for p in par if p is not None}
    return [v for v in range(len(par)) if v not in kids]

def branches(par):
    return {l: frozenset(anc(par, l)) for l in leaves(par)}

def embeddings(p, q):           # parent-preserving, root->root, injective
    n, m = len(p), len(q)
    for img in itertools.permutations(range(m), n):
        if img[0] != 0: continue
        if all(q[img[v]] == img[p[v]] for v in range(1, n)):
            yield img

def compose(jf, jg):
    return tuple(jg[x] for x in jf)

def canon(par):
    n = len(par); kids = {v: [w for w in range(n) if par[w] == v] for v in range(n)}
    def c(v): return "(" + "".join(sorted(c(w) for w in kids[v])) + ")"
    return c(0)

skel = {}
for n in range(1, 6):
    for t in trees(n):
        skel.setdefault(canon(t), t)
S = list(skel.values())
BR = [branches(t) for t in S]
nS = len(S)
H = {(a, b): list(embeddings(S[a], S[b])) for a in range(nS) for b in range(nS)}
STAR = S.index(skel[canon((None, 0, 0, 0))])
star_leaves = sorted(BR[STAR])                     # three leaves of the 3-leaf star
L1, L2, L3 = star_leaves

def sign(perm):                                     # perm: dict leaf->leaf
    xs = list(perm); s = 1
    for i in range(len(xs)):
        for k in range(i + 1, len(xs)):
            if (xs[i] - xs[k]) * (perm[xs[i]] - perm[xs[k]]) < 0: s = -s
    return s

def make_F(mode, park):
    """mode: 'sign' (odd automorphisms of the star swap branches L1,L2, fix L3), 'canonical' (star acts by
    its own branch permutation), all other trees act trivially. Non-isos: constant at c_T, positions -> root.
    park: the branch of the star used as constant target."""
    cst = {a: min(BR[a]) for a in range(nS)}; cst[STAR] = park
    def aut_phi(a, j):
        if a != STAR: return {x: x for x in BR[a]}
        perm = {l: j[l] for l in BR[a]}
        if mode == "canonical": return perm
        return {L1: L2, L2: L1, L3: L3} if sign(perm) < 0 else {x: x for x in BR[a]}
    def F(a, b, j, variance):
        if a == b:
            phi = aut_phi(a, j)
            if variance == "contra":                  # contravariant: act by the inverse
                phi = {v: k for k, v in phi.items()}
            # positions: chain_{phi x} -> chain_x, root->root, leaf->leaf
            if a != STAR:
                rho = {x: {v: v for v in BR[a][x]} for x in BR[a]}
            else:
                rho = {x: {0: 0, phi[x]: x} for x in BR[a]}
            return phi, rho
        if variance == "co":
            c = cst[b]; return {x: c for x in BR[a]}, {x: {v: 0 for v in BR[b][c]} for x in BR[a]}
        c = cst[a]; return {x: c for x in BR[b]}, {x: {v: 0 for v in BR[a][c]} for x in BR[b]}
    return F

def ccomp(m1, m2):
    (p1, r1), (p2, r2) = m1, m2
    return ({s: p2[p1[s]] for s in p1}, {s: {v: r1[s][r2[p1[s]][v]] for v in r2[p1[s]]} for s in p1})

def test(F, variance):
    bad = pairs = 0
    for a, b, c in itertools.product(range(nS), repeat=3):
        for jf in H[(a, b)]:
            for jg in H[(b, c)]:
                pairs += 1
                lhs = F(a, c, compose(jf, jg), variance)
                rhs = ccomp(F(a, b, jf, variance), F(b, c, jg, variance)) if variance == "co" \
                    else ccomp(F(b, c, jg, variance), F(a, b, jf, variance))
                if lhs != rhs: bad += 1
    ident_bad = sum(F(a, a, tuple(range(len(S[a]))), variance)[0] != {x: x for x in BR[a]} for a in range(nS))
    return pairs, bad, ident_bad

print("\n== Part 1: crown Rem. 20 'escapes are exactly the ones acting trivially on automorphisms'")
print(f"   skeleton <=5 nodes: {nS} objects; 3-leaf star has {len(H[(STAR, STAR)])} automorphisms")
for variance in ("co", "contra"):
    F = make_F("sign", L3)
    pairs, bad, ib = test(F, variance)
    nontriv = sum(F(STAR, STAR, j, variance)[0] != {x: x for x in BR[STAR]} for j in H[(STAR, STAR)])
    noncanon = sum(F(STAR, STAR, j, variance)[0] != ({l: j[l] for l in BR[STAR]} if variance == "co"
                   else {j[l]: l for l in BR[STAR]}) for j in H[(STAR, STAR)])
    print(f"   [{variance}] sign action, constant parked on fixed branch: {pairs} composable pairs, "
          f"composition failures {bad}, identity failures {ib}; automorphisms acting non-trivially {nontriv}/6, "
          f"non-canonically {noncanon}/6")
    for label, mode, park in [("NEG (i) canonical star action", "canonical", L3),
                              ("NEG (ii) sign action, constant on moved branch", "sign", L1)]:
        p2, b2, i2 = test(make_F(mode, park), variance)
        print(f"   [{variance}] {label}: composition failures {b2} (must be > 0)")

print("\n== Part 2: crown l.483-486 'read backward, R_f is a left-total injective partial map br(T')->br(T)'")
OBJ = [t for n in range(1, 5) for t in trees(n)]
tot = nontotal = noninj = 0; ex = None
for a, b in itertools.product(range(len(OBJ)), repeat=2):
    Bs, Bt = branches(OBJ[a]), branches(OBJ[b])
    for j in embeddings(OBJ[a], OBJ[b]):
        tot += 1
        back = {}
        for y in Bt:
            xs = [x for x in Bs if {j[v] for v in Bs[x]} <= Bt[y]]
            assert len(xs) <= 1
            if xs: back[y] = xs[0]
        if len(back) < len(Bt): nontotal += 1
        if len(set(back.values())) < len(back):
            noninj += 1
            ex = ex or (OBJ[a], OBJ[b], j, back)
print(f"   {tot} T-morphisms (labelled trees <=4 nodes): backward partial map non-total in {nontotal}, "
      f"non-injective in {noninj}; it is always single-valued (assert passed) and surjective onto br(T)")
print(f"   first non-injective example: T={ex[0]}, T'={ex[1]}, j={ex[2]}, b'->b: {ex[3]}")

print("\n== Part 3: statement drift (Thm/Lemma/Prop/Cor/Def/Example blocks), old -> new PDF")
def txt(pdf):
    return subprocess.run(["pdftotext", "-layout", pdf, "-"], capture_output=True, text=True).stdout
def blocks(t):
    """all occurrences 'Theorem N (' / 'Theorem N.' etc. -> list of normalised windows up to the next proof/env."""
    t = re.sub(r"\s+", " ", t)
    out = {}
    for h in re.finditer(r"(Theorem|Lemma|Proposition|Corollary) (\d+)(?= \([^()]{1,90}(\([^()]*\))?[^()]{0,40}\)\.)", t):
        key = h.group(0)            # statement heads only: 'Theorem 13 (Title).'  (not in-prose references)
        body = t[h.start(): h.start() + 700]
        m = re.search(r" (Proof\.|Reason\.|Computation\.|Remark \d+ \(|Theorem \d+ \(|Lemma \d+ \(|Proposition \d+ \(|Definition \d+\.|Example \d+ \()", body[len(key):])
        body = body[: len(key) + m.start()] if m else body
        out.setdefault(key, []).append(re.sub(r"[^A-Za-z]", "", body))
    return out
pairs = [("Note A (nerve)", "2026-10-11-composition-is-a-nerve-level-two-datum-6089268.pdf",
          "2026-10-13-uid378-composition-is-a-nerve-level-two-datum.pdf"),
         ("Note B (uniqueness)", "2026-10-11-container-derivative-uniqueness-6089268.pdf",
          "2026-10-13-uid378-container-derivative-uniqueness.pdf"),
         ("NEG control: Note B 75603ad -> 6089268 (Prop 12 WAS restated; must be flagged)",
          "2026-10-09-container-derivative-uniqueness-75603ad.pdf",
          "2026-10-11-container-derivative-uniqueness-6089268.pdf")]
for name, o, n in pairs:
    bo, bn = blocks(txt(os.path.join(PEER, o))), blocks(txt(os.path.join(PEER, n)))
    keys = sorted(bo)
    missing = [k for k in keys if k not in bn]
    worst = {}
    for k in keys:
        if k not in bn: continue
        # every old occurrence must reappear verbatim (letters only) among the new occurrences
        worst[k] = min(max(difflib.SequenceMatcher(None, x, y).ratio() for y in bn[k]) for x in bo[k])
    changed = {k: round(v, 3) for k, v in worst.items() if v < 1.0}
    print(f"   {name}: {len(keys)} Thm/Lemma/Prop/Cor labels; missing in new {missing}; "
          f"labels whose text is not reproduced verbatim (best similarity): {changed}")
# Note B Theorem 1 residual (0.999): the 700-char window runs past the statement into page furniture.
# Cut at the statement's end (next paragraph 'The question is Cockett') and compare exactly:
def thm1(pdf):
    t = re.sub(r"\s+", " ", txt(os.path.join(PEER, pdf))); i = t.index("Theorem 1 (Uniqueness)")
    b = t[i:i + 900]; return re.sub(r"[^A-Za-z]", "", b[: b.index("The question is")])
print("   Note B Theorem 1 cut at statement end: identical =",
      thm1("2026-10-11-container-derivative-uniqueness-6089268.pdf") == thm1("2026-10-13-uid378-container-derivative-uniqueness.pdf"))
