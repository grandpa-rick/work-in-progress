"""Day 242 Q473 checks. Run: python3 run_checks.py  (output -> run_checks.out)"""
import random
from gdmw_modp import (P, inv, Params, G, subsets, rank_modp, preceq, contained,
                       dominance_nocard)

par = Params(seed=242)
R = random.Random(4730)


def rnd():
    return R.randrange(2, P - 1)


print("== C1 stability (Clio Lemma 12): G(x_1..x_{L-1},1) == G(x_1..x_{L-1}) ==")
ok = tot = 0
for L in range(1, 6):
    for k in subsets(6, L):
        xs = [rnd() for _ in range(L - 1)]
        tot += 1
        ok += G(par, k, xs + [1]) == G(par, k, xs)
print(f"  {ok}/{tot}")
print("== C1' stability at x=-1 (h_inf(-1)=0, zeta(-1)=1): ==")
ok = tot = 0
for L in range(1, 5):
    for k in subsets(6, L):
        xs = [rnd() for _ in range(L - 1)]
        tot += 1
        ok += G(par, k, xs + [P - 1]) == G(par, k, xs)
print(f"  {ok}/{tot}")


def Phi(par, k, kp, shift=0):
    """Phi_k(G_kp) via stability = G_kp(xi_1..xi_d), xi_i = 1/y_{k_i}."""
    xi = [inv(par.y[m]) for m in k]
    return G(par, kp, xi, shift)


def Dk(par, k):
    q = par.q
    xi = [inv(par.y[m]) for m in k]
    d = len(k)
    v = 1
    for i in range(d):
        v = v * par.hinf(xi[i]) % P
        for j in range(i + 1, d):
            v = v * (1 - xi[i] * xi[j]) % P * inv(1 - q * xi[i] * xi[j]) % P
            v = v * (xi[j] - q * xi[i]) % P * inv(xi[j] - xi[i]) % P
        for j in range(1, k[i]):
            r = par.y[j] * inv(par.y[k[i]]) % P
            v = v * (1 - r) % P * inv(1 - q * r) % P
    return v


print("== C2 evaluation triangularity, kappa,kappa' subset [6], |.|<=4 (incl. empty row) ==")
S = subsets(6, 4)
mis = {"preceq": 0, "contained": 0, "dominance_nocard": 0}
zeros = nz = diag_ok = diag_formula = 0
tab = {}
for k in S:
    for kp in S:
        v = Phi(par, k, kp)
        tab[(k, kp)] = v
        if v == 0:
            zeros += 1
        else:
            nz += 1
        for name, rel in (("preceq", preceq), ("contained", contained),
                          ("dominance_nocard", dominance_nocard)):
            mis[name] += (v != 0) != rel(kp, k)
    diag_ok += tab[(k, k)] != 0
    diag_formula += tab[(k, k)] == Dk(par, k)
print(f"  cells {len(S)**2}: zero {zeros}, nonzero {nz}")
print(f"  mismatches vs Clio order preceq: {mis['preceq']}")
print(f"  NEG CONTROL wrong order 'containment': mismatches {mis['contained']}")
print(f"  NEG CONTROL wrong order 'entrywise, no cardinality': mismatches {mis['dominance_nocard']}")
print(f"  diagonal nonzero {diag_ok}/{len(S)}, equals closed form (8): {diag_formula}/{len(S)}")

print("== C2' direct (no stability lemma) for |kappa| = L = 4: Phi = plain evaluation ==")
S4 = [k for k in S if len(k) == 4]
bad = 0
for k in S4:
    for kp in S:
        bad += (tab[(k, kp)] != 0) != preceq(kp, k)
print(f"  rows {len(S4)}, mismatches {bad}")

print("== NEG CONTROL planted phi error (prod to j<=m): diagonal must die ==")
dz = sum(Phi(par, k, k, shift=1) == 0 for k in S if k)
print(f"  diagonal zero in {dz}/{len(S)-1} nonempty rows (expected all)")

print("== Count audit of Clio's L=3 tables (kappa,kappa' subset [4], |.|<=3) ==")
S3 = subsets(4, 3)
full = sum(preceq(kp, k) for k in S3 for kp in S3)
ev_nz = sum(preceq(kp, k) for k in S3 if k for kp in S3)
res_nz = sum(preceq(k, kp) for k in S3 if k for kp in S3)
print(f"  full 15x15 #pairs kp<=k = {full} (= #pairs k<=kp by transposition)")
print(f"  rows kappa nonempty (14x15=210): eval nonzero predicted {ev_nz} (Clio 109),"
      f" residue nonzero predicted {res_nz} (Clio 95)")

print("== C3 rank at generic point (two seeds), kappa subset [6], |kappa|<=L ==")
for seed in (1, 2):
    pr = Params(seed=1000 + seed)
    for L in (2, 3, 4, 5):
        fam = subsets(6, L)
        N = len(fam) + 6
        pts = [[rnd() for _ in range(L)] for _ in range(N)]
        M = [[G(pr, k, x) for x in pts] for k in fam]
        rk = rank_modp(M)
        line = f"  seed {seed} L={L}: #functions {len(fam)}, #points {N}, rank {rk}"
        if seed == 1 and L in (3, 4):
            # NEG CONTROL 1: replace one row by a combination of others -> rank must drop by 1
            M2 = [row[:] for row in M]
            i0, i1, i2 = len(fam) - 1, 3, 7
            M2[i0] = [(2 * u - 5 * v) % P for u, v in zip(M[i1], M[i2])]
            line += f" | planted-dependency rank {rank_modp(M2)}"
            # NEG CONTROL 2: degenerate y_1 = 0 -> phi_1 == 0 -> G_kappa == 0 for 1 in kappa
            pz = Params(seed=1000 + seed, ys={1: 0})
            Mz = [[G(pz, k, x) for x in pts] for k in fam]
            n1 = sum(1 in k for k in fam)
            line += f" | y1=0 rank {rank_modp(Mz)} (predicted {len(fam)-n1})"
        print(line)
