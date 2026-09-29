"""Explicit witnesses + sanity checks.
(1) enumerator count cross-check: Z2xZ4 should give 28 labelled braces (Rick's 09-25 w4 script).
(2) MacBeth's |G|=16 witnesses G0, G1: validity, L, [nu], sigma, beta, Omega.
(3) print explicit braces for the obstruction examples.
"""
from enum_braces import AbGroup, enumerate_braces

G = AbGroup((2, 4))
br, A = enumerate_braces(G)
print("labelled braces on Z2xZ4:", len(br))
G8 = AbGroup((8,))
print("labelled braces on Z8:", len(enumerate_braces(G8)[0]))
G16 = AbGroup((16,)); G28 = AbGroup((2, 8))
print("labelled braces on Z16:", len(enumerate_braces(G16)[0]), " on Z2xZ8:", len(enumerate_braces(G28)[0]))


def check_brace(elts, add, lam):
    """lam(a) returns a function G->G. verify lam_a additive bijective and lam_{a o b} = lam_a lam_b."""
    for a in elts:
        La = lam(a)
        assert len({La(x) for x in elts}) == len(elts)
        for x in elts:
            for y in elts:
                assert La(add(x, y)) == add(La(x), La(y))
    for a in elts:
        for b in elts:
            c = add(a, lam(a)(b))
            for x in elts:
                assert lam(c)(x) == lam(a)(lam(b)(x))
    return True


def report(name, elts, add, lam, D, neg, dmult):
    circ = lambda a, b: add(a, lam(a)(b))
    inv = {a: next(b for b in elts if circ(a, b) == elts[0]) for a in elts}
    check_brace(elts, add, lam)
    assert all(lam(a)(d) in D for a in elts for d in D)
    assert all(circ(a, b) == circ(b, a) for a in D for b in D)
    out = [a for a in elts if a not in D]
    Lset = sorted({dmult(lam(d)) for d in D})
    nus = sorted({dmult(lam(k)) for k in out})
    sig = sorted({dmult(lambda d, k=k: circ(inv[k], circ(d, k))) for k in out})
    twoD = {add(d, d) for d in D}
    beta0 = any(add(k, k) in twoD for k in out)
    om0 = any(add(k, k) == elts[0] and circ(k, k) == elts[0] for k in out)
    print(f"{name}: valid brace, D ideal, L={Lset}, [nu]={nus}, sigma(mult)={sig}, "
          f"beta={'0' if beta0 else '!=0'}, Omega={'0' if om0 else '!=0'}")


# MacBeth G0, G1 on Z2 x Z8, D = {(0,x)}
E = [(i, x) for i in range(2) for x in range(8)]
add = lambda a, b: ((a[0] + b[0]) % 2, (a[1] + b[1]) % 8)
U = lambda a: (pow(5, a[1] % 2, 8) * pow(3, a[0], 8)) % 8
D = [(0, x) for x in range(8)]
def dm(f):
    return f((0, 1))[1]  # multiplier on D = <(0,1)>, returns image of generator (for sigma may be non-mult)
for tflag, name in [(0, "MacBeth G0"), (1, "MacBeth G1")]:
    lam = lambda a, tf=tflag: (lambda y: (y[0], (U(a) * y[1] + y[0] * 4 * a[0] * tf) % 8))
    report(name, E, add, lam, D, None, dm)

# Obstruction example O8 (trivial kernel, RY regime): Z2xZ4, D = 0xZ4, nu = 1, sigma = -1.
# lambda_(i,x)(j,y) = (j, y + 2 j x)   [lambda_d(k) = k - 2d on D, identity on D]
E4 = [(i, x) for i in range(2) for x in range(4)]
add4 = lambda a, b: ((a[0] + b[0]) % 2, (a[1] + b[1]) % 4)
lamO8 = lambda a: (lambda y: (y[0], (y[1] + 2 * y[0] * a[1]) % 4))
report("O8 witness (beta=0 realises nu=1,sigma=-1)", E4, add4, lamO8, [(0, x) for x in range(4)], None, dm)

# Order 16, sigma = id, L = {1,5}:
# (Z16)  lambda_a = mult by u(a);  search the enumerator output for sigma=id, nu in {3,7}
from cross_eq_c_H2 import analyse

# W16a: Z/16, lambda_a = mult by (1+2a) mod 16; D = 2Z/16 ~ Z/8 (gen 2)
E16 = list(range(16)); add16 = lambda a, b: (a + b) % 16
lam16 = lambda a: (lambda y: ((1 + 2 * a) * y) % 16)
dm16 = lambda f: (f(2) // 2) % 8
report("W16a Z/16, lambda_a=1+2a", E16, add16, lam16, [2 * x for x in range(8)], None, dm16)
# W16b: Z2xZ8, lambda_(i,x)(j,y) = (j, 5^x y); D = 0xZ8
lamb = lambda a: (lambda y: (y[0], (pow(5, a[1], 8) * y[1]) % 8))
report("W16b Z2xZ8, lambda=(j,5^x y)", E, add, lamb, D, None, dm)
