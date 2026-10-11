#!/usr/bin/env python3
"""Independent check of Clio Q410 Prop 17: Korff's Psi (arXiv:1110.6356, (5.13), Def 5.5)
equals Macdonald's psi_{Lambda/mu} (SFHP III (5.8'), from memory) times one wrap factor.

Implemented in KORFF'S OWN COORDINATES (not via Clio's dictionary):
  (n,k) = Korff's (n,k); A+_{k,n} = partitions with exactly k parts, each in [1,n].
  lambda/d/mu with d in {0,1}; conjugate loop lambda'[0]_{i+l n} = lambda'_i - l k  ((5.8)ff),
  lambda'[d] = lambda'[0] + d.  theta_i = lambda'[d]_i - mu'[0]_i.
  Horizontal strip <=> theta_i in {0,1} for all i (periodic, so i=1..n suffices).
  J = {1<=i<=n : theta_i = 0, theta_{i+1} = 1}, theta_{n+1} from the loop (not assumed = d).
  Psi = prod_{i in J} (1 - t^{m_i(mu)}),  m_i(mu) = #parts of mu equal to i.
  Phi = prod_{i in I} (1 - t^{m_i(lambda)}), I: theta_i=1, theta_{i+1}=0.
Macdonald side (ordinary partitions, row-based interlacing test, independent of theta code):
  Lambda = lambda with d extra rows of length n.  Lambda/mu horizontal strip iff
  Lambda_1 >= mu_1 >= Lambda_2 >= mu_2 >= ...;  psi_{Lambda/mu} = prod_{j in J'} (1-t^{m_j(mu)}),
  J' = {j>=1 : Lambda'_j - mu'_j = 0, Lambda'_{j+1} - mu'_{j+1} = 1}.
"""
import itertools, sys
import sympy as sp

t = sp.symbols('t')


def parts_exact(k, n):
    # partitions with exactly k parts in [1,n], weakly decreasing
    for c in itertools.combinations_with_replacement(range(n, 0, -1), k):
        yield tuple(c)


def conj(p, length):
    return [sum(1 for x in p if x >= j) for j in range(1, length + 1)]


def mult(p, i):
    return sum(1 for x in p if x == i)


def loopconj(pc, n, k, i):
    # pc = conjugate list (index 1..n stored 0..n-1); lambda'[0]_{i + l n} = lambda'_i - l k
    l, r = divmod(i - 1, n)
    return pc[r] - l * k


def korff(lam, mu, d, n, k):
    lc, mc = conj(lam, n), conj(mu, n)
    th = {i: loopconj(lc, n, k, i) + d - loopconj(mc, n, k, i) for i in range(1, n + 2)}
    if any(th[i] not in (0, 1) for i in range(1, n + 1)):
        return None
    J = [i for i in range(1, n + 1) if th[i] == 0 and th[i + 1] == 1]
    I = [i for i in range(1, n + 1) if th[i] == 1 and th[i + 1] == 0]
    Psi = sp.prod([1 - t**mult(mu, i) for i in J])
    Phi = sp.prod([1 - t**mult(lam, i) for i in I])
    return th, J, I, sp.expand(Psi), sp.expand(Phi)


def macdonald_psi(Lam, mu):
    L = list(Lam); m = list(mu) + [0] * (len(Lam) - len(mu))
    if len(m) > len(L):
        return None
    # interlacing L1>=m1>=L2>=m2>=...
    for i in range(len(L)):
        if not (L[i] >= m[i]):
            return None
        if i + 1 < len(L) and not (m[i] >= L[i + 1]):
            return None
    W = max(L) + 2
    Lc, mc = conj(L, W), conj(m, W)
    th = [Lc[j] - mc[j] for j in range(W)]  # th[j-1] = theta'_j
    assert all(x in (0, 1) for x in th)
    Jp = [j for j in range(1, W) if th[j - 1] == 0 and th[j] == 1]
    return sp.expand(sp.prod([1 - t**mult(mu, j) for j in Jp]))


def bfac(p, n):
    out = 1
    for i in range(1, n + 1):
        for r in range(1, mult(p, i) + 1):
            out *= (1 - t**r)
    return sp.expand(out)


def main():
    ranges = [(n, k) for n in range(2, 7) for k in range(1, 5)]
    tot = 0; claim_fail = 0; gate_fail = 0; d1 = 0
    ctrl = {'C0 no wrap': 0, 'C1 wrap iff d=1': 0, 'C2 exp m_1(mu)': 0,
            'C3 exp m_n(lambda)': 0, 'C4 wrap at 1->0 (th_n=1,th_1=0)': 0}
    for (n, k) in ranges:
        A = list(parts_exact(k, n))
        for mu in A:
            for lam in A:
                for d in (0, 1):
                    r = korff(lam, mu, d, n, k)
                    if r is None:
                        continue
                    th, J, I, Psi, Phi = r
                    tot += 1; d1 += d
                    assert th[n + 1] == th[1] == d  # Korff's theta_{n+1}=theta_1, and = d
                    # gate: Korff (5.16) b_lam Psi = b_mu Phi
                    if sp.expand(bfac(lam, n) * Psi - bfac(mu, n) * Phi) != 0:
                        gate_fail += 1
                    Lam = tuple([n] * d + list(lam))
                    psi = macdonald_psi(Lam, mu)
                    assert psi is not None, (n, k, lam, mu, d)
                    wrap_on = (th[n] == 0 and th[1] == 1)
                    W = (1 - t**mult(mu, n)) if wrap_on else 1
                    if sp.expand(psi * W - Psi) != 0:
                        claim_fail += 1
                    preds = {
                        'C0 no wrap': psi,
                        'C1 wrap iff d=1': psi * ((1 - t**mult(mu, n)) if d == 1 else 1),
                        'C2 exp m_1(mu)': psi * ((1 - t**mult(mu, 1)) if wrap_on else 1),
                        'C3 exp m_n(lambda)': psi * ((1 - t**mult(lam, n)) if wrap_on else 1),
                        'C4 wrap at 1->0 (th_n=1,th_1=0)': psi * ((1 - t**mult(mu, n)) if (th[n] == 1 and th[1] == 0) else 1),
                    }
                    for key, p in preds.items():
                        if sp.expand(p - Psi) != 0:
                            ctrl[key] += 1
    print(f"(n_K,k_K) in n=2..6, k=1..4: cylindric horizontal strips = {tot} (d=1: {d1})")
    print(f"Korff (5.16) gate failures: {gate_fail}")
    print(f"CLAIM Psi = psi_Mac * wrap failures: {claim_fail}/{tot}")
    for key, v in ctrl.items():
        print(f"  control {key}: fails {v}/{tot}")


if __name__ == '__main__':
    main()
