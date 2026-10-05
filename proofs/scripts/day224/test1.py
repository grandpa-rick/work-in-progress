from ops import *
import sympy
t = Fr(3)
ok = 0; tot = 0
for k in range(1, 6):
    for r in range(1, 6):
        if k+r > 8: continue
        A = toe(Dop(k, r, t), t); B = {kk: v for kk, v in Mthm(k, r, t).items() if v}
        tot += 1; ok += (A == B)
        if A != B: print('D mismatch', k, r, A, B)
print('D vs Thm7.1', ok, tot)
# E2 on e_b e_c vs log (2,2,2) at t=3
ts = sympy.symbols('t')
log = {(6,): (ts**2 + 1)**2*(ts**2 - ts + 1)*(ts**2 + ts + 1), (5, 1): ts**6 + 2*ts**5 + ts**4 + 4*ts**3 + 2*ts**2 + ts + 2, (4, 2): -(ts**2 + 1)*(ts**2 - 2*ts + 5), (4, 1, 1): ts**2 + ts + 1, (3, 3): ts**2 + ts + 1, (3, 2, 1): -5*ts - 6, (2, 2, 2): 6}
X = toe(E2prod(2, 2, 2, t), t)
print(X); print({k: v.subs(ts, 3) for k, v in log.items()})
log2 = {(6,): (ts + 1)*(ts**2 - ts + 1)**2*(ts**2 + ts + 1)**2, (5, 1): ts*(2*ts**4 + 2*ts**3 + 3*ts**2 + 2*ts + 3), (4, 2): -(ts + 1)*(ts**2 + 2), (4, 1, 1): -2*(ts**2 + ts + 1), (3, 2, 1): 3}
print(toe(E2prod(3,2,1,t),t)); print({k: v.subs(ts, 3) for k, v in log2.items()})
