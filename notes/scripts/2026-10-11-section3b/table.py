"""Emit the column-B table rows (|lambda|<=8, y<=x) as LaTeX, recomputed from ground truth."""
import sympy as sp
from section3b import ground_truth, two_row_partitions, S3b, t
rows = []
for n in range(2, 9):
    for yy in range(1, n // 2 + 1):
        xx = n - yy
        gt = ground_truth(n, xx, yy)
        for lam in two_row_partitions(n):
            X = gt[lam]
            assert S3b(lam, xx, yy)[0] == X
            rows.append(f"$({lam[0]},{lam[1]})$ & $({xx},{yy})$ & ${sp.latex(sp.Poly(X, t).as_expr())}$ \\\\")
half = (len(rows) + 1) // 2
open('table_left.tex', 'w').write('\n'.join(rows[:half]) + '\n')
open('table_right.tex', 'w').write('\n'.join(rows[half:]) + '\n')
print(len(rows), half)
