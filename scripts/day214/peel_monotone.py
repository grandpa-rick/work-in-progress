from ds_pointeval import parts, dominates
def peel(kap, k):
    return tuple(sorted((x for x in [kap[i]-1 if i < k else kap[i] for i in range(len(kap))] if x > 0), reverse=True))
bad = tot = 0
for n in range(1, 12):
    P = list(parts(n))
    for lam in P:
        for kap in P:
            if not dominates(lam, kap): continue
            for k in range(1, len(lam)+1):
                tot += 1
                if not dominates(peel(lam, k), peel(kap, k)): bad += 1; print('FAIL', lam, kap, k) if bad < 5 else None
print('checked', tot, 'bad', bad)
