"""Day 207b: the proved general Pieri formula, specialized. Coefficient of e_{r+n} e_b (n=k-b) is s^b F_n(t^{r-b}).
(1) k=1: Thm 3.12. (2) k=2: W_r (Day 206b). (3) k=3: Day 193 c_3..c_0. (4) k=4: Day 195 c_4..c_0. All symbolic in u=t^r, s=1/q.
(5) r=0: sum_b s^b F_{k-b}(t^{-b}) e_b e_{k-b} = e_k (Hikita Lemma 3.3 normalization), k<=7."""
import sympy as sp
s, t, u, q, x = sp.symbols('s t u q x')
def poch(a, n):
    r = sp.Integer(1)
    for i in range(n): r *= (1 - a*t**i)
    return r
def alpha(j):
    if j < 0: return sp.Integer(0)
    r = sp.Integer(1)
    for i in range(1, j+1): r *= (s - t**i)
    return r/poch(t, j)
def c(n, j):
    if j < 0 or n < j: return sp.Integer(0)
    return poch(s, n-j)/poch(t, n-j)*(alpha(j) - s*t**(n-j)*alpha(j-1))
def F(n, w): return sum(c(n, j)*w**j for j in range(n+1))
def coef(k, b): return s**b*F(k-b, u*t**(-b))   # coefficient of e_{r+k-b} e_b, u = t^r
ok = True
def chk(tag, a, b):
    global ok
    d = sp.cancel(sp.together((a - b).subs(s, 1/q)))
    ok &= (d == 0); print(tag, 'OK' if d == 0 else f'FAIL {d}')
br = lambda k: (1 - u*t**k)/(1 - t)
bc = lambda n: sum(t**i for i in range(n))
# k=1: e_1*e_r = (1-s)[r+1] e_{r+1} + s e_1 e_r
chk('k=1 b=0', coef(1, 0), (1 - 1/q)*br(1)); chk('k=1 b=1', coef(1, 1), 1/q)
# k=2: W_r
chk('k=2 b=2', coef(2, 2), 1/q**2); chk('k=2 b=1', coef(2, 1), (1/q)*(1 - 1/q)*br(0))
chk('k=2 b=0', coef(2, 0), (1 - 1/q)*br(2)/bc(2)*(br(1) - (1/q)*(br(0) - 1)))
# k=3: Day 193
chk('k=3 b=3', coef(3, 3), 1/q**3); chk('k=3 b=2', coef(3, 2), (q-1)/q**3*br(-1))
chk('k=3 b=1', coef(3, 1), (q-1)/q**3*br(1)/bc(2)*(q*br(0) - t*br(-2)))
chk('k=3 b=0', coef(3, 0), (q-1)/q**3*br(3)/(bc(2)*bc(3))*(br(1)*br(2)*q**2 - t*bc(2)*br(-1)*br(1)*q + t**3*br(-2)*br(-1)))
# k=4: Day 195
pre = (q-1)/q**4
P4 = br(1)*br(2)*br(3)*q**3 - t*bc(3)*br(-1)*br(1)*br(2)*q**2 + t**3*bc(3)*br(-2)*br(-1)*br(1)*q - t**6*br(-3)*br(-2)*br(-1)
for b, tgt in [(4, 1/q**4), (3, pre*br(-2)), (2, pre*br(0)/bc(2)*(q*br(-1) - t*br(-3))),
               (1, pre*br(2)/(bc(2)*bc(3))*(br(1)*br(0)*q**2 - t*bc(2)*br(-2)*br(0)*q + t**3*br(-3)*br(-2))),
               (0, pre*br(4)/(bc(2)*bc(3)*bc(4))*P4)]:
    chk(f'k=4 b={b}', coef(4, b), tgt)
# r=0 normalization
for k in range(1, 8):
    # e_b e_{k-b} and e_{k-b} e_b are the same monomial: pair b with k-b
    good = sp.cancel(F(k, 1) + s**k - 1) == 0 and all(sp.cancel(s**b*F(k-b, t**(-b)) + (s**(k-b)*F(b, t**(-(k-b))) if 2*b != k else 0)) == 0 for b in range(1, k))
    ok &= good; print(f'r=0 normalization k={k}', 'OK' if good else 'FAIL')
print('ALL OK' if ok else 'SOME FAIL')
