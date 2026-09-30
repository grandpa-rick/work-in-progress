# Lemma (greedy peel): kappa ⊴ rho (same size), rho_k > rho_{k+1} (so rho-1^k is a partition).
# Then sort(kappa - 1 on its k largest parts) ⊴ rho - 1^k  (requires kappa to have >= k nonzero parts, automatic).
from ds_pointeval import parts, dominates
bad = 0; tot = 0
for n in range(1, 13):
    P = list(parts(n))
    for rho in P:
        r = list(rho) + [0]
        for k in range(1, len(rho)+1):
            if not r[k-1] > r[k]: continue
            rt = tuple(x for x in (r[i]-1 if i < k else r[i] for i in range(len(r))) if x > 0)
            for kap in P:
                if not dominates(rho, kap): continue
                tot += 1
                if len(kap) < k: bad += 1; print('short', rho, k, kap); continue
                kt = tuple(sorted((x for x in [kap[i]-1 if i < k else kap[i] for i in range(len(kap))] if x > 0), reverse=True))
                if not dominates(rt, kt): bad += 1; print('FAIL', rho, k, kap, kt, rt)
print('checked', tot, 'bad', bad)
