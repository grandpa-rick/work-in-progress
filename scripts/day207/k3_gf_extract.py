"""Extract e-basis coefficients of t^{-3} e_3(Y).e_r from the telescoped GF, symbolic in u=t^r,
and compare with the Day 193 closed form (s = 1/q)."""
import sympy as sp, pickle, sys
sys.path.insert(0, '/home/agent/projects/scripts/day207')
s, t, z, v, u, q = sp.symbols('s t z v u q')
I1, Qsym = pickle.load(open('/home/agent/projects/scripts/day207/k3_gf_I1.pkl', 'rb'))
Qa, Qb, Qc = Qsym[sp.factor(-z)], Qsym[sp.factor(-t*z)], Qsym[sp.factor(-t**2*z)]
P = sp.Poly(I1, Qa, Qb, Qc)
chain = {(0,0,0): 0, (1,0,0): 1, (1,1,0): 2, (1,1,1): 3}
C = {}
for mon, coef in P.terms():
    assert mon in chain, ('NON-TELESCOPING monomial', mon)
    C[chain[mon]] = sp.expand(coef)
print('telescoping: all Q-monomials are chain prefixes ->', sorted(C))
# GF:  [3]!(1-t)^3 sum_r z^r t^{-3}e3(Y)e_r = sum_j C_j(z) E(t^j z)
e = sp.symbols('e0:8'); E1, E2, E3 = e[1], e[2], e[3]
# q_n in e's: Q(v) = E(-t v)/E(-v)   (3 generic e's enough up to q3)
Ev = lambda c: 1 + c*E1*v + c**2*E2*v**2 + c**3*E3*v**3
ser = sp.series(Ev(-t)/Ev(-1), v, 0, 4).removeO()
qsub = {sp.Symbol(f'q{n}'): sp.expand(ser.coeff(v, n)) for n in (1, 2, 3)}
# symbolic e_{r+a}: symbols R_a ; [z^r] z^{-a}... : term coef*z^b*E(t^j z) -> coef * t^{j(r-b)} e_{r-b}
R = {a: sp.Symbol(f'R{a}') for a in range(-3, 4)}   # R_a = e_{r+a}
tot = 0
for j, Cj in C.items():
    Cj = sp.expand(Cj.subs(qsub))
    for b in range(-3, 4):
        cb = Cj.coeff(z, b) if b != 0 else Cj.subs(z, 0) if not Cj.has(1/z) else sp.expand(Cj*z**3).coeff(z, 3)
        cb = sp.expand(Cj*z**3).coeff(z, b+3)
        if cb == 0: continue
        tot += cb * u**j * t**(-j*b) * R[-b]
fac = (1-t)**3*(1+t)*(1+t+t**2)
tot = sp.expand(tot)
# collect by e-monomials
res = sp.Poly(tot, R[0], R[1], R[2], R[3], R[-1], R[-2], R[-3], E1, E2, E3)
coeffs = {}
for mon, c in res.terms():
    cc = sp.factor(sp.cancel(c / fac))
    coeffs[mon] = cc
    names = ['e_r','e_{r+1}','e_{r+2}','e_{r+3}','e_{r-1}','e_{r-2}','e_{r-3}','e1','e2','e3']
    print(' * '.join(f'{n}^{k}' if k > 1 else n for n, k in zip(names, mon) if k), ':', cc)
# Day 193 closed form, [n] = (1-t^n)/(1-t) with t^r = u
def br(k):  # [r+k]
    return (1 - u*t**k)/(1-t)
def bc(n): return sum(t**i for i in range(n))
c3 = 1/q**3
c2 = (q-1)/q**3*br(-1)
c1 = (q-1)/q**3*br(1)/bc(2)*(q*br(0) - t*br(-2))
c0 = (q-1)/q**3*br(3)/(bc(2)*bc(3))*(br(1)*br(2)*q**2 - t*bc(2)*br(-1)*br(1)*q + t**3*br(-2)*br(-1))
target = {(1,0,0,0,0,0,0,0,0,1): c3, (0,1,0,0,0,0,0,0,1,0): c2, (0,0,1,0,0,0,0,1,0,0): c1, (0,0,0,1,0,0,0,0,0,0): c0}
allmons = set(coeffs) | set(target)
ok = True
for mon in allmons:
    d = sp.simplify(coeffs.get(mon, 0).subs(s, 1/q) - target.get(mon, 0))
    if d != 0: ok = False; print('MISMATCH', mon, sp.factor(d))
print('match Day 193 (symbolic in u=t^r, generic r):', ok)
pickle.dump((C, coeffs), open('/home/agent/projects/scripts/day207/k3_gf_coeffs.pkl', 'wb'))
