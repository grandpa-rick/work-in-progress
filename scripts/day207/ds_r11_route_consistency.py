"""
Day 207: symbolic-in-r consistency of the two routes to DS(r,1,1).
Direct route:  C_r = (1-s)[r+1] e_1(Y).e_{r+1} + s e_1(Y).(e_1 e_r)   (Thm 3.12 op form + Sub-Lemma Z)
Day 198 route: C_r = p_2(Y).e_r + 2t W_r                             (Newton + Lemma 1 + W_r, Day 206b)
[n]_t written via u := t^r, so the check is valid for every r >= 2.
"""
import sympy as sp
q, t, u = sp.symbols('q t u')  # u = t^r
s = 1/q
def br(k):  # [r+k]_t
    return (1 - u*t**k)/(1 - t)
two = 1 + t
# direct route
e1Y_er1 = {'r+2': (1-s)*br(2), 'r+1,1': s}
Z = {'r+2': (1-s)**2*br(2), 'r+1,1': (1-s)*(t*br(0)+s), 'r,2': s*(1-s)*two, 'r,1,1': s**2}
direct = {}
for k, v in e1Y_er1.items(): direct[k] = direct.get(k, 0) + (1-s)*br(1)*v
for k, v in Z.items(): direct[k] = direct.get(k, 0) + s*v
closed = {'r+2': (1-s)**2*br(2)*(br(1)+s), 'r+1,1': s*(1-s)*(br(1)+t*br(0)+s),
          'r,2': s**2*(1-s)*two, 'r,1,1': s**3}
# Day 198 route
tau = -(q**2-1)*br(2)*(q*u*t - q + t + 1)/(q**3*two)
p2 = {'r,1,1': 1/q**3, 'r,2': -(q*t-q+t+1)/q**3, 'r+1,1': (q**2-1)/q**3, 'r+2': tau}
W = {'r,2': s**2, 'r+1,1': s*(1-s)*br(0), 'r+2': (1-s)*br(2)/two*(br(1) - s*(br(0)-1))}
d198 = {k: p2.get(k, 0) + 2*t*W.get(k, 0) for k in closed}
ok = True
for k in closed:
    a = sp.simplify(direct[k] - closed[k]); b = sp.simplify(d198[k] - closed[k])
    print(k, 'direct-closed:', a, ' day198-closed:', b)
    ok &= (a == 0 and b == 0)
print('ALL OK' if ok else 'FAIL')
