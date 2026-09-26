"""Day 207 k=4 (general k): iterated-Lemma-2 GF engine, generalizing k3_gf_engine.py.
  (1-t)^k [k]_t! sum_r z^r t^{-C(k,2)} e_k(Y).e_r = E(z) * I_1,
  I_1 = Omega_{x1} Omega_{x2} ... Omega_{xk} [ prod_i phi(x_i) ],  phi(u) = u(1+s u z)/(1+u z),
  Omega_{x_i}: ambient P(v) = Q(v) prod_{l<i} h(x_l, v),  h(u,v) = (1-uv)/(1-tuv); Omega[u^n] = [v^n]P, n>=1.
Engine asserts simple poles and f(0)=0 at every stage (Omega from k3_gf_engine, unchanged).
Usage: python3 k4_gf_engine.py K"""
import sympy as sp, sys, pickle, time
sys.path.insert(0, '/home/agent/projects/scripts/day207')
import k3_gf_engine as E3
from k3_gf_engine import Omega, phi, h, Qsym, s, t, z, v

def run(k):
    xs = sp.symbols(f'x1:{k+1}')
    I = sp.Integer(1)
    for i in range(k, 0, -1):
        t0 = time.time()
        amb = sp.Integer(1)
        for l in range(i - 1): amb *= h(xs[l], v)
        I = Omega(sp.expand(phi(xs[i - 1]) * I), xs[i - 1], amb)
        print(f'stage x{i} done: {len(sp.Add.make_args(I))} terms [{time.time()-t0:.1f}s]'); sys.stdout.flush()
    return I, dict(Qsym)

if __name__ == '__main__':
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    I1, Q = run(k)
    print('Q symbols:', Q)
    pickle.dump((I1, Q), open(f'/home/agent/projects/scripts/day207/k{k}_gf_I1_gen.pkl', 'wb'))
    P = sp.Poly(I1, *Q.values())
    for mon, coef in P.terms():
        print(mon, sp.factor(coef))
