"""Re-run of the Sec. 3(b) j-sum negative control (Clio UID 380 F3), 2026-10-13.

Same control as section3b.py (drop the (t^-j - t^j) G_{1,y-1-j} sum), but reported on
BOTH row sets: |lambda|<=8 (what the 10-11 note printed) and |lambda|<=10 (the 155-row sweep),
stratified by whether the j-sum is non-empty (y > lambda_2) or empty (y <= lambda_2).
"""
import sympy as sp
from section3b import t, two_row_partitions, ground_truth, G1_coeffs

def dropped(lam, xx, yy):
    G = G1_coeffs(lam); g = lambda k: G.get(k, 0)
    Pr = sum(c * t**k for k, c in G.items())
    brace = (g(yy - 1) / (1 - t)**2) / t - Pr / (1 - t**2) * (1 + t**-yy)
    b = (1 - t)**2 if lam[0] > lam[1] else (1 - t) * (1 - t**2)
    return sp.expand(sp.cancel(b * brace))

for NMAX in (8, 10):
    rows = brk = ne = ne_brk = e_brk = 0
    for n in range(2, NMAX + 1):
        for xx in range(1, n):
            yy = n - xx
            gt = ground_truth(n, xx, yy)
            for lam in two_row_partitions(n):
                rows += 1
                nonempty = yy > lam[1]
                b = dropped(lam, xx, yy) != gt[lam]
                brk += b; ne += nonempty
                ne_brk += (b and nonempty); e_brk += (b and not nonempty)
    print(f'|lambda|<={NMAX}: rows {rows}; control breaks {brk}; j-sum non-empty {ne} '
          f'(broken {ne_brk}); empty {rows-ne} (broken {e_brk})')
