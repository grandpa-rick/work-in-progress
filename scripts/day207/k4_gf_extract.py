"""Day 207 k=4: extract M_{jb} and e-basis coefficients of t^{-6} e_4(Y).e_r from the telescoped GF,
symbolic in u = t^r, and compare with Day 195 e_4*e_r closed form (s = 1/q).
GF: (1-t)(1-t^2)(1-t^3)(1-t^4) sum_r z^r t^{-6}e_4(Y).e_r = sum_j C_j(z) E(t^j z),  C_j = coefficient of chain monomial j."""
import sympy as sp, pickle, sys
K = int(sys.argv[1]) if len(sys.argv) > 1 else 4
s, t, z, v, u, q = sp.symbols('s t z v u q')
I1, Qsym = pickle.load(open(f'/home/agent/projects/scripts/day207/k{K}_gf_I1_gen.pkl', 'rb'))
Qs = [Qsym[sp.factor(-t**i * z)] for i in range(K)]
assert set(Qsym.values()) == set(Qs), Qsym
P = sp.Poly(I1, *Qs)
chain = {tuple([1]*j + [0]*(K-j)): j for j in range(K+1)}
C = {}
for mon, coef in P.terms():
    assert mon in chain, ('NON-TELESCOPING monomial', mon)
    C[chain[mon]] = sp.expand(coef)
print('telescoping: Q-monomials are exactly chain prefixes ->', sorted(C))
e = sp.symbols('e0:10')
Ev = lambda c: 1 + sum(c**i * e[i] * v**i for i in range(1, K+1))
ser = sp.series(Ev(-t)/Ev(-1), v, 0, K+1).removeO()
qsub = {sp.Symbol(f'q{n}'): sp.expand(ser.coeff(v, n)) for n in range(1, K+1)}
# M_{jb}: C_j(z) = sum_b M_{jb} e_b z^{b-K}  (check linearity in e_b)
M = {}
linear = True
for j in range(K+1):
    Cj = sp.expand(C[j].subs(qsub) * z**K)
    assert not Cj.has(1/z)
    for b in range(0, 2*K+1):
        cb = sp.expand(Cj.coeff(z, b))
        if cb == 0: continue
        if b > K: linear = False; print('  z-power beyond K', j, b); continue
        Mjb = sp.cancel(cb / (e[b] if b > 0 else 1))
        if any(Mjb.has(e[i]) for i in range(1, K+1)): linear = False; print('  NONLINEAR', j, b, sp.factor(cb))
        M[(j, b)] = sp.factor(Mjb)
print('C_j linear in single e_b (M_{jb} scalars):', linear)
print('support {(j,b)}:', sorted(M))
for key in sorted(M): print(f'M_{key[0]}{key[1]} =', M[key])
# e-basis extraction, symbolic r: term M e_b z^{b-K} E(t^j z) -> [z^r]: M e_b t^{j(r+K-b)} e_{r+K-b}
R = {a: sp.Symbol(f'R{a}') for a in range(0, K+1)}   # R_a = e_{r+a}
tot = 0
for (j, b), Mjb in M.items():
    tot += Mjb * (e[b] if b > 0 else 1) * u**j * t**(j*(K-b)) * R[K-b]
fac = sp.Mul(*[(1 - t**i) for i in range(1, K+1)])
coeffs = {}
for a in range(K+1):
    for bb in range(K+1):
        c = sp.factor(sp.cancel(sp.expand(tot).coeff(R[a]).coeff(e[bb]) if bb > 0 else sp.expand(tot).coeff(R[a]).subs({e[i]: 0 for i in range(1, K+1)})) / fac)
        if c != 0: coeffs[(a, bb)] = c; print(f'coef of e_(r+{a}) e_{bb}:', c)
pickle.dump((C, M, coeffs), open(f'/home/agent/projects/scripts/day207/k{K}_gf_coeffs.pkl', 'wb'))
if K == 4:
    br = lambda k: (1 - u*t**k)/(1 - t)
    bc = lambda n: sum(t**i for i in range(n))
    pre = (q-1)/q**4
    c4 = 1/q**4
    c3 = pre*br(-2)
    c2 = pre*br(0)/bc(2)*(q*br(-1) - t*br(-3))
    c1 = pre*br(2)/(bc(2)*bc(3))*(br(1)*br(0)*q**2 - t*bc(2)*br(-2)*br(0)*q + t**3*br(-3)*br(-2))
    P4 = br(1)*br(2)*br(3)*q**3 - t*bc(3)*br(-1)*br(1)*br(2)*q**2 + t**3*bc(3)*br(-2)*br(-1)*br(1)*q - t**6*br(-3)*br(-2)*br(-1)
    c0 = pre*br(4)/(bc(2)*bc(3)*bc(4))*P4
    target = {(0, 4): c4, (1, 3): c3, (2, 2): c2, (3, 1): c1, (4, 0): c0}
    ok = True
    for key in set(coeffs) | set(target):
        d = sp.cancel(coeffs.get(key, 0).subs(s, 1/q) - target.get(key, 0))
        print(f'  {key}: diff = {d}')
        if d != 0: ok = False
    print('P4 match Day 195 (symbolic in u=t^r):', ok)
