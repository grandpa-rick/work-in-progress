"""Day 217: reconcile MacBeth's 112 = 8+96+8 and 24 = 8+16 counts; (nu, sigmahat) slice table."""
from collections import defaultdict
from lock_tab import run
for mods, dorder in [((16,), 8), ((2, 8), 8), ((2, 4), 4)]:
    A, R, rows = run(mods, dorder)
    pairs = defaultdict(set)
    for r in rows:
        pairs[(r['D'], r['L'])].add(r['brace'])
    Ds = sorted({r['D'] for r in rows}, key=lambda S: sorted(S))
    for D in Ds:
        gens = sorted(A.elts[x] for x in D)
        tot_nontriv = sum(len(v) for (DD, L), v in pairs.items() if DD == D and L != (1,))
        print(f"Z{mods} D={gens[:3]}...: (brace,D) pairs by L:",
              {L: len(v) for (DD, L), v in sorted(pairs.items(), key=str) if DD == D}, " nontrivial-L total:", tot_nontriv)
# (nu, sigmahat) slice table for L = {1,5}, D=Z/8
print("\nL={1,5} slices: (sigmahat, nu) -> b parities realised, count of (brace,D,u) triples")
tab = defaultdict(lambda: defaultdict(int))
for mods in [(16,), (2, 8)]:
    A, R, rows = run(mods, 8)
    for r in rows:
        if r['L'] == (1, 5):
            tab[(r['sh'], r['nu'])][r['b'] % 2] += 1
for k in sorted(tab):
    sh, nu = k
    print(f"  sigmahat={sh} nu={nu} nu*sh mod4={(nu*sh)%4}: b even {tab[k][0]}, b odd {tab[k][1]}, predicted offset c=2(nu*sh-1) mod 8 = {(2*(nu*sh-1))%8}")
