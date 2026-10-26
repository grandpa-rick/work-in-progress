"""(a) Prop 5.1 TFAE of mirror-LA (6df0043) and (d) the corrected M1/M2 table. Rick, Wake 258.
(1) Norm_K associative; (1') (RA#); (2) (mLA) for all c, x in ImK.ImK; (3) K(gh) theta_K gh, theta_K as printed in Sec 5.
theta_K is also recomputed by an independent refinement algorithm (largest congruence inside ~K)."""
import sys, time
from itertools import product
from common import *

def p51_report(t, n, e, K):
    th = theta_explicit(t, n, K)
    return dict(c1=norm_assoc(t, n, K), c1r=RAsharp(t, n, K), c2=mLA_restricted(t, n, K),
                c3=p51_three(t, n, K, th), theta_eq=(th == largest_congruence_in(t, n, K)),
                mLA=mLA_full(t, n, K), closed=closed(t, n, K))

print("=== (d) M1/M2 table: M={e,a,b}, aa=b, ab=ba=a, bb=b; K: e,a,b -> e,a,a ===")
e, a, b = 0, 1, 2
T3 = [[0, 1, 2], [1, 2, 1], [2, 1, 2]]
K3 = [0, 1, 1]
print(" monoid (assoc, unit e):", is_assoc(T3, 3) and all(T3[0][x] == x == T3[x][0] for x in range(3)))
print(" a^3 = a:", T3[T3[a][a]][a] == a, "| b = a^2:", T3[a][a] == b, "| b is identity of {a,b} (M = C2 + adjoined unit):", T3[b][a] == a and T3[b][b] == b)
print(" (K1),(U),(LA):", is_normalizer(T3, 3, 0, K3))
print(" ~K classes {e},{a,b} two-sided congruence:", theta_explicit(T3, 3, K3) == [[K3[x] == K3[y] for y in range(3)] for x in range(3)])
print(" Prop 5.1 report:", p51_report(T3, 3, 0, K3))
OLD = [[0, 1, 2], [1, 2, 1], [2, 1, 1]]   # pre-456 printed: b.b = a
bad = sum(OLD[OLD[x][y]][z] != OLD[x][OLD[y][z]] for x in range(3) for y in range(3) for z in range(3))
print(" negative control: pre-456 table (b.b=a) non-associative triples:", bad)

print("\n=== (a) paper examples ===")
SEP1 = [[0, 1, 2], [1, 1, 2], [2, 1, 2]]; K1 = [0, 0, 2]
print(" Thm 3.1 Sep.1 (closed, not mLA):", is_assoc(SEP1, 3), is_normalizer(SEP1, 3, 0, K1), p51_report(SEP1, 3, 0, K1))
# TR Sec 6 handed model M={e,z,a,b}: z zero, aa=z, ab=a, bb=b, ba=z (TR's own pending Lean note); K=id on e,a,b, K(z)=b
E, Z, A, B = 0, 1, 2, 3
TR = [[0, 1, 2, 3], [1, 1, 1, 1], [2, 1, 1, 2], [3, 1, 1, 3]]; KTR = [0, 3, 2, 3]
print(" TR handed model (ba=z):", is_assoc(TR, 4), is_normalizer(TR, 4, 0, KTR), p51_report(TR, 4, 0, KTR))
TRc = [r[:] for r in TR]; TRc[B][A] = A
print(" TR control ba=a: monoid?", is_assoc(TRc, 4))
# group examples (S3 witness H=<(12)>, T=A3; C4 H={0,2}, T={0,1})
S3 = sym(3); e3 = ident(S3); inv3 = inverses(S3, e3)
subs = subgroups(S3); A3 = [H for H in subs if len(H) == 3][0]
for H in [H for H in subs if len(H) == 2]:
    for T in transversals(S3, H, e3):
        if set(T) == set(A3):
            K = rep_map(S3, H, T); print(" S3 witness H order 2, T=A3:", p51_report(S3, 6, e3, K))
    break
C4 = cyclic(4); K = rep_map(C4, frozenset({0, 2}), (0, 1)); print(" C4 H={0,2} T={0,1}:", p51_report(C4, 4, 0, K))

print("\n=== (a) exhaustive: all monoids with unit 0, order <= N, all normalizers ===")
N = int(sys.argv[1]) if len(sys.argv) > 1 else 4
for n in range(1, N + 1):
    t0 = time.time(); mons = monoids_unit0(n)
    st = dict(monoids=len(mons), normalizers=0, assoc=0, mism_1_1r=0, mism_1_2=0, mism_1_3=0, theta_alg_mism=0,
              between=0, ctrl_noLA_K=0, ctrl_noLA_mism_1_3=0, ctrl_simK_mism=0)
    for t in mons:
        for Kt in product(range(n), repeat=n - 1):
            K = [0] + list(Kt)
            if not is_normalizer(t, n, 0, K, require_LA=False): continue
            r = p51_report(t, n, 0, K)
            if not is_normalizer(t, n, 0, K):        # control: drop the (LA) hypothesis of Prop 5.1
                st['ctrl_noLA_K'] += 1; st['ctrl_noLA_mism_1_3'] += (r['c1'] != r['c3']); continue
            st['normalizers'] += 1; st['assoc'] += r['c1']
            st['mism_1_1r'] += r['c1'] != r['c1r']; st['mism_1_2'] += r['c1'] != r['c2']; st['mism_1_3'] += r['c1'] != r['c3']
            st['theta_alg_mism'] += not r['theta_eq']
            st['between'] += r['c1'] and not r['mLA'] and not r['closed']
            st['ctrl_simK_mism'] += (not r['c1'])   # control: (3) with ~K in place of theta_K is always true
    print(f" order {n}: {st}  [{time.time()-t0:.1f}s]", flush=True)
