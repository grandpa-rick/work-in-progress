from lead2 import *
def Xgreen(lam, rho, t):  # coefficient of P_lam in p_rho
    n = sum(lam); co = S(n, t).to_hl({tuple(rho): Fr(1)}); return co.get(lam, 0)
for n in range(3, 7):
    for lam in parts(n):
        row = []
        for x in range(n-1, 0, -1):
            y = n - x
            if x < y: break
            p = interp(lambda t: Xgreen(lam, (x, y), t), 2*n*n)
            row.append(f'({x},{y}):{p}')
        print(lam, ' X_n=', interp(lambda t: Xgreen(lam, (n,), t), 2*n*n), ' | ', '  '.join(row), flush=True)
