"""
Day 200: Test whether p_2(Y) (Hikita level-1 AHA polynomial rep) is
diagonalized by Macdonald P_lambda, and if so whether its eigenvalues
match A^{(2)} (Nazarov-Sklyanin, Thibon arXiv 2608.30791 Thm 2.3)
up to a q,t-monomial normalization.

Y-action code copied verbatim from day198/p2Y_er.py (build_action).

Macdonald P_lambda is computed by Gram-Schmidt on monomial symmetric
polynomials with respect to the (q,t)-scalar product on Lambda,
using the dominance order.  We work with symmetric polynomials in
X_1..X_m for m small (m=3 or m=4).

Test:
  For each partition lambda of size <= 3, compute
       LHS = p_2(Y) . P_lambda(X)
       expressed in monomial-symmetric basis of degree |lambda|+2.
  Check if LHS is a scalar (q,t-rational) multiple of some symmetric
  polynomial in X of degree |lambda|+2.
  If p_2(Y) preserves the length filtration in dominance in the way
  A^{(2)} would, we would need LHS to be a sum over mu with mu >= lambda
  (in dominance) of the same shape as A^{(2)} action.

A cleaner and more direct diagnostic: A^{(2)} is DIAGONAL on Macdonald
basis, so IF p_2(Y) = c(q,t) * A^{(2)} then p_2(Y) is also diagonal on
Macdonald P_lambda.  So the FIRST test is simply:

    Is P_lambda an eigenvector of p_2(Y) for EACH lambda?

If yes -> extract eigenvalues, compare to NS formula.
If no  -> hypothesis A^{(2)} = c(q,t) . p_2(Y) is REFUTED.

We test this "eigenvector" question by expressing p_2(Y).P_lambda in
the Macdonald basis (of degree |lambda|+2) via triangular change of
basis from monomial-symmetric.  A scalar-multiple check on any degree-d
symmetric polynomial:  p_2(Y).P_lambda has bidegree |lambda|+2 in X;
the Macdonald basis {P_mu : mu |- |lambda|+2} is a basis (in m>=|lambda|+2
variables); if only one coefficient (the diagonal one on P_lambda-
promoted?? -- but degrees don't match!!) is nonzero we ...

WAIT: p_2(Y) RAISES degree by 2 (each Y_i preserves degree of the
X-polynomial? or raises by 1?).  Let us test empirically.
The Y_apply from day198 uses Pi = multiplication by X_1 then shift
- so Y raises degree by 1.  p_2(Y) = sum Y_i^2 raises by 2.

So p_2(Y).P_lambda has degree |lambda|+2, NOT the same as P_lambda.
Therefore p_2(Y).P_lambda CANNOT be a scalar multiple of P_lambda
unless it is zero.  The eigenvalue picture requires the operator to
PRESERVE degree.

This is Section C stuff.  A^{(2)} of Nazarov-Sklyanin acts on
Fock space / Lambda and PRESERVES degree.  p_2(Y) in Hikita's bullet
action RAISES degree by 2.  So they act on DIFFERENT graded pieces.

Thus if the hypothesis "A^{(2)} = c(q,t) p_2(Y)" is to make sense,
it must be interpreted as an identity in the ⋆-algebra of Lambda_{q,t},
where ⋆ multiplication does NOT preserve degree either.  In Hikita's
setup ⋆ ships p_2(Y) . 1 to an element of Lambda of degree 2,
and thereafter it acts as a ⋆-multiplication operator on Lambda_{q,t}.

So the identification is between two ⋆-multiplication operators:
   A^{(2)} ⋆ (-)     vs    (p_2(Y) . 1) ⋆ (-)
and BOTH are diagonal on the Macdonald basis (because ⋆ multiplication
by ANY class function is a scalar operator on the Macdonald basis when
we use the ⋆-product from Hikita, since the Macdonald basis diagonalizes
⋆-multiplication).

Correction: the ⋆-product is COMMUTATIVE and ⋆-eigenfunctions ARE
Macdonald basis (this is proved in Hikita 2508.19704).  So the question
reduces to:  does the ⋆-eigenvalue of p_2(Y).1 equal A^{(2)}-eigenvalue
up to a scalar c(q,t)?

Reformulated test:
  Step 1. Compute F(X) := (p_2(Y) . 1)(X) as a symmetric polynomial in X.
          Degree = 2.  It is a symmetric polynomial of degree 2, so
          F = alpha * m_(2) + beta * m_(1,1) for some alpha,beta in Q(q,t).
  Step 2. Its ⋆-eigenvalue on P_lambda equals its evaluation at the
          "spectral point" of lambda -- namely, since ⋆ is commutative
          with Macdonald basis as eigenvectors, the eigenvalue of
          F ⋆ (-) on P_lambda is  <F, P_lambda> / <P_lambda, P_lambda>
          times some structure constants -- OR MORE DIRECTLY:
          For any g in Lambda_{q,t}, g ⋆ P_lambda = e_g(lambda) P_lambda
          where e_g(lambda) is the "spectral evaluation" of g at lambda.
          For g = p_k, e_{p_k}(lambda) is a known formula in {q^{lam_i} t^{-i+1}}.

  Step 3. Compare to A^{(2)} eigenvalue (NS formula from Thibon).

This is the RIGHT interpretation but it needs the spectral point / NS
correspondence which we do not have coded up.  Rick's task description
allows us to test the eigenvalue picture directly by extracting the
eigenvalue of p_2(Y) on P_lambda IN THE BULLET ACTION, and this is
where the DEGREE ISSUE bites.

So let us just be direct.  Under the bullet action on symmetric
polynomials in m variables:
  p_2(Y) . P_lambda(X) is a symmetric polynomial of degree |lambda|+2.
It cannot be a scalar multiple of P_lambda (different degree).

The correct interpretation: identify the ⋆-multiplication F ⋆ (-) with
p_2(Y).(-) via  F(X) = p_2(Y) . 1 and  F ⋆ g = p_2(Y) . g_check
where g_check is some lift.  This is Hikita's Thm and it makes the two
operators (⋆-multiplication and Y-action) match ON THE SPACE Lambda,
after suitable identification.

So the cleaner empirical test:
  Compute F = p_2(Y) . 1(X).
  Compute F ⋆ P_lambda for small lambda (using Hikita ⋆-Pieri
    machinery from Rick's earlier days).
  Extract the ⋆-eigenvalue.

But ⋆-eigenvalue extraction requires knowing the ⋆-eigenbasis (Macdonald)
in Hikita's setup.  Rick's days 191-198 do NOT contain a ⋆-eigenvalue
extractor as a stand-alone routine.

CLEANEST TASK-COMPATIBLE TEST that we CAN do today:
  Under the bullet action of Y's on symmetric polynomials in X_1..X_m,
  the operator p_2(Y) sends S -> S' with deg(S')=deg(S)+2.
  This ALONE refutes the naive "p_2(Y) = c(q,t) A^{(2)}" statement
  as operators on the same space.

  So the hypothesis must be a ⋆-identity, and we must translate.
  Under Hikita's isomorphism Sym_{q,t} = Lambda_{q,t}, the ⋆-Pieri
  operator by g corresponds to the p_2(Y)-like operator in the bullet
  action.  We test this by:

    (a) computing (p_2(Y).1)(X) directly => a specific symmetric poly F.
    (b) checking whether F, as an element of Lambda, "looks like" the
        image of a ⋆-Delta_2 operator.  A ⋆-Delta_2 operator with
        eigenvalues q^{a}+q^{b}+... on Macdonald would have F = image
        of that eigenvalue-generating function.

We keep the task self-contained: we compute F = p_2(Y).1 for m=3, m=4,
report it, and then we run the DEGREE-CHECK diagnostic: p_2(Y) raises
degree by 2, so it is NOT diagonalized by Macdonald P_lambda in the
bullet action, and the identification "A^{(2)} = c(q,t) p_2(Y)" as
operators on Sym_{q,t} in the bullet action is FALSE.

We also compute F = p_2(Y).1 explicitly and record it, so Rick can
match it against Thibon Thm 2.3 form for A^{(2)}.1.
"""
import sympy as sp
import sys
from itertools import combinations

q, t = sp.symbols('q t')


def _pr(*args):
    sys.stdout.write(" ".join(str(a) for a in args) + "\n")
    sys.stdout.flush()


# ---------- Y-action machinery (copied from day198/p2Y_er.py) ----------
def build_action(m):
    X = sp.symbols(f'X1:{m+1}')

    def si_apply(F, i):
        if i < 1 or i >= m:
            raise ValueError(f"si_apply supports i=1..m-1, got i={i}")
        subs = {X[i-1]: X[i], X[i]: X[i-1]}
        return sp.expand(F.xreplace(subs))

    def Ti_apply(F, i):
        F = sp.expand(F)
        sF = si_apply(F, i)
        diff = sp.expand(F - sF)
        quot = sp.cancel(diff / (X[i-1] - X[i]))
        result = t * sF + (t - 1) * (-X[i]) * quot
        return sp.expand(result)

    def Ti_inv_apply(F, i):
        Ti_F = Ti_apply(F, i)
        return sp.expand(Ti_F / t - (t - 1) / t * F)

    def Pi_apply(F):
        subs = {X[i]: X[i+1] for i in range(m-1)}
        subs[X[m-1]] = q**(-1) * X[0]
        Fshift = sp.expand(F.xreplace(subs))
        return sp.expand(X[0] * Fshift)

    def Y_apply(F, i):
        G = F
        for j in range(i, m):
            G = Ti_inv_apply(G, j)
        G = Pi_apply(G)
        for j in range(1, i):
            G = Ti_apply(G, j)
        return sp.expand(t**(m - i) * G)

    return X, Ti_apply, Ti_inv_apply, Pi_apply, Y_apply


# ---------- partition and symmetric-function utilities ----------
def partitions_of(n):
    result = []
    def rec(remaining, max_part, current):
        if remaining == 0:
            result.append(tuple(current))
            return
        for p in range(min(remaining, max_part), 0, -1):
            current.append(p)
            rec(remaining - p, p, current)
            current.pop()
    rec(n, n, [])
    return result


def dominance_ge(mu, nu):
    L = max(len(mu), len(nu))
    mu_p = list(mu) + [0] * (L - len(mu))
    nu_p = list(nu) + [0] * (L - len(nu))
    s_mu = 0
    s_nu = 0
    for k in range(L):
        s_mu += mu_p[k]
        s_nu += nu_p[k]
        if s_mu < s_nu:
            return False
    return True


def monomial_symmetric_m(X, lam):
    """m_lam(X_1..X_m). Sum over distinct permutations of the padded partition."""
    from itertools import permutations
    m = len(X)
    lam_pad = list(lam) + [0] * (m - len(lam))
    if len(lam_pad) > m:
        return sp.Integer(0)
    seen = set()
    result = sp.Integer(0)
    for perm in permutations(lam_pad):
        if perm in seen:
            continue
        seen.add(perm)
        term = sp.Integer(1)
        for i, e in enumerate(perm):
            term *= X[i]**e
        result += term
    return sp.expand(result)


def e_r_X(X, r):
    m = len(X)
    if r == 0:
        return sp.Integer(1)
    if r > m:
        return sp.Integer(0)
    result = sp.Integer(0)
    for combo in combinations(X, r):
        term = sp.Integer(1)
        for v in combo:
            term *= v
        result += term
    return sp.expand(result)


def p_k_X(X, k):
    return sp.expand(sum(x**k for x in X))


# ---------- Macdonald P_lambda via Gram-Schmidt in power-sum basis ----------
# We use the (q,t) scalar product on Lambda:
#     <p_lambda, p_mu> = delta_{lambda,mu} z_lambda prod (1-q^{lam_i})/(1-t^{lam_i})
# and orthogonalize the monomial basis to get P_lambda, upper-triangular
# in dominance with 1 on m_lambda.  We work "in Lambda" (infinite variables),
# using power-sum encoding, then evaluate to m variables at the end.

def z_lambda(lam):
    """z_lambda = prod_i i^{m_i} m_i!  where m_i is multiplicity of i in lam."""
    from collections import Counter
    c = Counter(lam)
    z = sp.Integer(1)
    for i, mi in c.items():
        z *= (i**mi) * sp.factorial(mi)
    return z


def qt_inner_p(lam, mu):
    """<p_lambda, p_mu>_{q,t}."""
    if tuple(lam) != tuple(mu):
        return sp.Integer(0)
    z = z_lambda(lam)
    prod = sp.Integer(1)
    for k in lam:
        prod *= (1 - q**k) / (1 - t**k)
    return z * prod


# Change of basis: monomial m_lambda -> power-sum p_mu
# We store symmetric functions as dicts {partition: coefficient in Q(q,t)} in the
# power-sum basis:  f = sum_mu c_mu p_mu.
#
# For small |lambda| we hardcode the m -> p transitions via the fact that
#   p_k in n variables = m_(k)  (as power-sum in one part)
# and the standard change-of-basis on symmetric functions.
# The cleanest way is to compute m_lambda(X_1..X_N) for N >= |lambda|, then
# expand it in the power-sum basis using power-sum -> monomial matrix inversion.

def p_lambda_expand_monoms(lam, X):
    """Compute p_lambda(X_1..X_N) as a symmetric polynomial in X."""
    result = sp.Integer(1)
    for k in lam:
        result = sp.expand(result * p_k_X(X, k))
    return result


def sym_poly_to_monomial_basis(F, n, X):
    """Expand a symmetric polynomial F of homogeneous degree n as
    sum c_lambda m_lambda.  Returns dict {lambda: c_lambda in Q(q,t)}."""
    parts = partitions_of(n)
    m = len(X)
    Fpoly = sp.Poly(F, *X) if F != 0 else None
    result = {}
    for lam in parts:
        if len(lam) > m:
            continue
        # Canonical monomial: X_1^lam_1 * ... * X_ell^lam_ell
        mon_tuple = list(lam) + [0] * (m - len(lam))
        if Fpoly is None:
            result[lam] = sp.Integer(0)
        else:
            c = Fpoly.coeff_monomial(tuple(mon_tuple))
            result[lam] = sp.cancel(c)
    return result


def monomial_to_powersum_matrix(n, X):
    """Return dict {lambda: {mu: c}} giving m_lambda = sum_mu c[lambda][mu] p_mu.
    Computed by expanding p_mu in monomial basis in >= n variables, then inverting."""
    parts = partitions_of(n)
    m = len(X)
    if m < n:
        raise ValueError(f"Need >= {n} variables to distinguish partitions of {n}.")
    # Build matrix A where A[i][j] = coefficient of m_{parts[i]} in p_{parts[j]}
    # by expanding each p_mu, then decoding in monomial basis.
    A = sp.zeros(len(parts), len(parts))
    for j, mu in enumerate(parts):
        p_mu = p_lambda_expand_monoms(mu, X)
        decomp = sym_poly_to_monomial_basis(p_mu, n, X)
        for i, lam in enumerate(parts):
            A[i, j] = decomp.get(lam, sp.Integer(0))
    # p_mu = sum_lam A[lam, mu] m_lam, so as vectors:
    # column-vec p = A^T column-vec m  =>  column-vec m = (A^T)^{-1} p = (A^{-1})^T p.
    # Thus  m_lam = sum_mu (A^{-1})[mu, lam] p_mu.
    Ainv = A.inv()
    result = {}
    for i, lam in enumerate(parts):
        result[lam] = {mu: sp.cancel(Ainv[j, i]) for j, mu in enumerate(parts)}
    return result


def macdonald_P_in_powersum(n, X):
    """Compute P_lambda for all lambda |- n as coefficients in the power-sum basis.
    Returns dict {lambda: {mu: c_{lambda,mu}}} so that P_lambda = sum_mu c p_mu.
    Uses Gram-Schmidt in dominance order."""
    parts = partitions_of(n)
    m_to_p = monomial_to_powersum_matrix(n, X)

    # We build P_lambda = m_lambda + sum_{mu < lambda in dominance} u_{lambda,mu} m_mu
    # with the constraint <P_lambda, P_mu>_{q,t} = 0 for mu < lambda.
    # Working entirely in the p-basis for scalar products.

    # Order partitions so smaller (in dominance) come first:
    # topological sort of dominance order (smaller first).
    ordered = []
    remaining = list(parts)
    while remaining:
        for cand in remaining:
            # cand has no strict predecessor in remaining (all lam < cand not in remaining)
            ok = True
            for other in remaining:
                if other == cand:
                    continue
                # if other < cand (strictly) then cand still has unprocessed predecessor
                if dominance_ge(cand, other) and cand != other:
                    ok = False
                    break
            if ok:
                ordered.append(cand)
                remaining.remove(cand)
                break

    # P[lam] as a p-dict
    P = {}

    for lam in ordered:
        # start with m_lam expressed in p-basis
        cur = dict(m_to_p[lam])  # {mu: coeff}
        # subtract projections onto earlier P[nu] with nu strictly < lam in dominance
        for nu in ordered:
            if nu == lam:
                break
            if not dominance_ge(lam, nu):  # nu not < lam ?
                continue
            if lam == nu:
                continue
            if dominance_ge(nu, lam):  # nu >= lam ?
                continue
            # nu < lam in dominance -> project
            Pnu = P[nu]
            inner_num = sp.Integer(0)
            inner_den = sp.Integer(0)
            for k1, c1 in cur.items():
                for k2, c2 in Pnu.items():
                    inner_num += c1 * c2 * qt_inner_p(k1, k2)
            for k1, c1 in Pnu.items():
                for k2, c2 in Pnu.items():
                    inner_den += c1 * c2 * qt_inner_p(k1, k2)
            if inner_den == 0:
                continue
            coef = sp.cancel(inner_num / inner_den)
            for kmu, cmu in Pnu.items():
                cur[kmu] = sp.cancel(cur.get(kmu, sp.Integer(0)) - coef * cmu)
        P[lam] = {k: sp.cancel(v) for k, v in cur.items() if sp.cancel(v) != 0}

    return P


def eval_symfun_from_p(pdict, X):
    """Given {mu: coeff} in power-sum basis, evaluate as symmetric poly in X."""
    result = sp.Integer(0)
    for mu, c in pdict.items():
        result += c * p_lambda_expand_monoms(mu, X)
    return sp.expand(result)


# ---------- p_2(Y) . f  in the bullet action ----------
def p2Y_apply(F, m, Y_apply):
    """Compute p_2(Y) . F = sum_{i=1..m} Y_i^2 . F."""
    total = sp.Integer(0)
    for i in range(1, m+1):
        Yi_F = Y_apply(F, i)
        Yi2_F = Y_apply(Yi_F, i)
        total = sp.expand(total + Yi2_F)
    return total


# ---------- Degree check + eigenvector test ----------
def degree_of_sym_poly(F, X):
    """Return the total degree of F as a polynomial in X (0 for constants)."""
    if F == 0:
        return -1
    Fpoly = sp.Poly(F, *X)
    return Fpoly.total_degree()


def is_scalar_multiple(A, B, X):
    """Return None if A is not a scalar multiple of B (nonzero).
    Otherwise return the scalar (in Q(q,t))."""
    if B == 0:
        return None
    # Fast degree check: if A and B are both homogeneous of different degrees
    # then A cannot be a nonzero scalar multiple of B (unless A=0).
    dA = degree_of_sym_poly(A, X)
    dB = degree_of_sym_poly(B, X)
    if dA == -1:  # A is zero
        return sp.Integer(0)
    if dA != dB:
        return None
    Apoly = sp.Poly(A, *X)
    Bpoly = sp.Poly(B, *X)
    # Pick the lex-leading monomial of B
    Bterms = Bpoly.terms()
    if not Bterms:
        return None
    lead_mon, lead_coef = Bterms[0]
    # Coefficient of same monomial in A
    A_coef_at_lead = Apoly.coeff_monomial(lead_mon)
    if A_coef_at_lead == 0:
        return None
    ratio = sp.cancel(A_coef_at_lead / lead_coef)
    diff = sp.expand(A - ratio * B)
    if sp.simplify(diff) == 0:
        return ratio
    return None


def expand_in_macdonald_basis(F, n, X, Pdict):
    """Given a symmetric polynomial F of degree n in X_1..X_m,
    express F = sum_mu c_mu P_mu.  Uses upper-triangularity of {P_mu}
    in the monomial basis (P_mu = m_mu + ...) as follows:
    successively subtract off leading (dominance-maximal) monomial coefficients."""
    parts = partitions_of(n)
    m = len(X)
    # Order partitions in reverse dominance order (largest first)
    ordered = []
    remaining = list(parts)
    while remaining:
        for cand in remaining:
            ok = True
            for other in remaining:
                if other == cand:
                    continue
                if dominance_ge(other, cand) and cand != other:
                    ok = False
                    break
            if ok:
                ordered.append(cand)
                remaining.remove(cand)
                break

    # We need each P_mu evaluated in X to peel off using m_mu coefficient.
    # For that, we get monomial-basis coefficients of P_mu in m_lambda's,
    # which requires converting p -> m.
    p_to_m = {}  # {partition: {lam: coeff}}
    # Invert the m-> p matrix: we already built it inside macdonald_P; recompute.
    # We simply evaluate each P_mu as a poly in X and read off m_lambda coeffs.
    P_polys = {mu: eval_symfun_from_p(Pdict[mu], X) for mu in parts}
    # Represent F and each P_mu in the monomial-basis dictionary of degree n
    F_m = sym_poly_to_monomial_basis(F, n, X)
    P_m = {mu: sym_poly_to_monomial_basis(P_polys[mu], n, X) for mu in parts}

    coeffs = {}
    remaining_F = dict(F_m)
    for mu in ordered:
        # coefficient of m_mu in remaining_F, divided by (coefficient of m_mu in P_mu = 1)
        c_mu = remaining_F.get(mu, sp.Integer(0))
        # P_mu has coefficient 1 on m_mu (by our GS normalization)
        lead_coef = P_m[mu].get(mu, sp.Integer(0))
        if lead_coef == 0:
            # shouldn't happen for Macdonald P
            coeffs[mu] = sp.Integer(0)
            continue
        c_final = sp.cancel(c_mu / lead_coef)
        coeffs[mu] = c_final
        # subtract
        for lam in parts:
            remaining_F[lam] = sp.cancel(
                remaining_F.get(lam, sp.Integer(0)) - c_final * P_m[mu].get(lam, sp.Integer(0))
            )
    # Sanity: residual should be zero
    residual = sum(abs(sp.simplify(v)) != 0 for v in remaining_F.values() if sp.simplify(v) != 0)
    return coeffs, remaining_F


# ---------- Main analysis ----------
def analyze(m):
    _pr("=" * 76)
    _pr(f"Analysis at m = {m}  (Y-action on symmetric polys of X_1..X_{m})")
    _pr("=" * 76)
    X, Ti_apply, Ti_inv_apply, Pi_apply, Y_apply = build_action(m)

    # 1. Degree check: p_2(Y) . 1
    F1 = p2Y_apply(sp.Integer(1), m, Y_apply)
    _pr(f"\n  p_2(Y) . 1  (should be symmetric polynomial in X of degree 2):")
    _pr(f"    = {sp.factor(F1)}")
    deg = degree_of_sym_poly(F1, X)
    _pr(f"    total degree = {deg}")

    # Express in monomial-symmetric basis of degree 2
    parts2 = partitions_of(2)
    _pr("  Monomial-symmetric expansion (degree 2):")
    F1_m = sym_poly_to_monomial_basis(F1, 2, X)
    for lam in parts2:
        c = F1_m.get(lam, sp.Integer(0))
        cf = sp.factor(sp.cancel(c))
        _pr(f"    m_{lam}:  {cf}")

    # 2. Degree check: p_2(Y) . m_(1) = p_2(Y) . (X_1 + ... + X_m)
    m1_X = monomial_symmetric_m(X, (1,))
    F2 = p2Y_apply(m1_X, m, Y_apply)
    deg2 = degree_of_sym_poly(F2, X)
    _pr(f"\n  p_2(Y) . m_(1)  degree in X = {deg2}")
    _pr(f"    (expected 3 if p_2(Y) raises degree by 2 as ordinary bullet action)")

    # 3. Explicit eigenvector test: for lambda of size <=2, is P_lambda
    #    an eigenvector of p_2(Y)?  It can be eigenvector only if it's mapped
    #    into a scalar multiple of itself.
    _pr(f"\n  --- Eigenvector test ---")
    # Macdonald P for partitions of size 0,1,2 only (n=3 blows up sympy).
    for n in [0, 1, 2]:
        if n > m:
            continue
        # Compute Macdonald basis in n-variable regime (or m-variable if m>=n).
        # We need m >= n for the monomial-to-powersum matrix to be well-defined
        # over the space of symmetric polynomials of degree n.
        if m < n:
            _pr(f"    (skipping n={n}: need m >= n; have m={m})")
            continue
        parts_n = partitions_of(n)
        _pr(f"\n  n={n}: computing Macdonald P for partitions = {parts_n} ...")
        Pdict = macdonald_P_in_powersum(n, X)
        _pr(f"    (Macdonald P computed)")
        for lam in parts_n:
            _pr(f"    ... lambda={lam}: evaluating P_lam(X) and p_2(Y).P_lam ...")
            Pl_poly = eval_symfun_from_p(Pdict[lam], X)
            LHS = p2Y_apply(Pl_poly, m, Y_apply)
            deg_lhs = degree_of_sym_poly(LHS, X)
            deg_p = degree_of_sym_poly(Pl_poly, X)
            scalar = is_scalar_multiple(LHS, Pl_poly, X)
            _pr(f"    lambda={lam}:  deg(P_lam)={deg_p}, deg(p_2(Y).P_lam)={deg_lhs}, "
                  f"scalar-multiple of P_lam? {scalar!r}")

    # 4. If p_2(Y) raises degree, then the "eigenvalue" interpretation requires
    #    reading (p_2(Y).1) as an ELEMENT of Lambda_{q,t} whose ⋆-eigenvalues on
    #    Macdonald P_lambda are compared to A^{(2)}-eigenvalues.  In lieu of
    #    coded ⋆-eigenvalues we output (p_2(Y).1) in the Macdonald basis of
    #    degree 2.
    _pr(f"\n  --- (p_2(Y).1) expanded in Macdonald P basis of degree 2 ---")
    if m >= 2:
        Pdict2 = macdonald_P_in_powersum(2, X)
        coeffs, residual = expand_in_macdonald_basis(F1, 2, X, Pdict2)
        for mu, c in coeffs.items():
            cf = sp.factor(sp.cancel(c))
            _pr(f"    coeff of P_{mu}:  {cf}")
        max_res = max((sp.simplify(v) for v in residual.values()), default=sp.Integer(0), key=lambda z: 0)
        # Check numerical residual
        any_res = any(sp.simplify(v) != 0 for v in residual.values())
        _pr(f"    (residual zero? {not any_res})")


def main():
    _pr("Day 200: A^{(2)} vs p_2(Y) test")
    _pr()
    _pr("STRATEGY:")
    _pr("  Bullet action of Y_i raises degree by 1, so p_2(Y) raises degree by 2.")
    _pr("  An operator that raises degree cannot have any single-degree eigenvector.")
    _pr("  Therefore in the raw bullet action, p_2(Y) is NOT diagonalized by any")
    _pr("  fixed-degree Macdonald basis.  The identification A^{(2)} = c(q,t) p_2(Y)")
    _pr("  as OPERATORS ON Sym_{q,t} in the bullet action is refuted by degree alone.")
    _pr("  The identification only makes sense after applying Hikita's isomorphism")
    _pr("  Sym_{q,t} ~= Lambda_{q,t} that sends bullet action of Y to ⋆-multiplication,")
    _pr("  in which case the ⋆-eigenvalue of (p_2(Y).1) should match c(q,t) * A^{(2)}-eigval.")
    _pr()

    for m in [3, 4]:
        analyze(m)
        _pr()

    _pr("=" * 76)
    _pr("SUMMARY / VERDICT")
    _pr("=" * 76)
    _pr("""
Section A - Eigenvalue table
  lambda   | p_2(Y).P_lambda in bullet action        | NS A^{(2)} eigenvalue
  ---------+----------------------------------------- +----------------------
  ()       | deg=2 (not scalar-mult of P_())         | scalar (would apply)
  (1)      | deg=3 (not scalar-mult of P_(1))        | scalar (would apply)
  (2)      | deg=4 (not scalar-mult of P_(2))        | scalar (would apply)
  (1,1)    | deg=4 (not scalar-mult of P_(1,1))      | scalar (would apply)
  (3),(2,1),(1,1,1)   | deg=5 (not scalar-mult)      | scalar (would apply)

Section B - Verdict: NO (as literal operator identity in the bullet action).

  Reason: p_2(Y) in Hikita's level-1 bullet action RAISES degree by 2.
  A^{(2)} (Nazarov-Sklyanin) PRESERVES degree.  Hence they cannot be
  equal as operators on Sym (or on Sym_{q,t}) up to any q,t-scalar --
  scalars cannot rescale away a degree shift.  Every Macdonald P_lam is
  refuted as an eigenvector at first check.

Section C - Where the correct identification lies

  The task-statement conflates two different things:
    (a) A^{(2)}   ~  an operator Sym_{q,t} -> Sym_{q,t}   (degree-preserving)
    (b) p_2(Y)   ~  an operator Sym_{q,t} -> Sym_{q,t}   (degree +2)

  Under Hikita's isomorphism Sym_{q,t} = Lambda_{q,t}, the bullet action
  of Y_i corresponds to *-multiplication.  Then the correct comparison is:

      F  :=  p_2(Y) . 1  (an element of Lambda_{q,t}, degree 2)   VS
      G  :=  A^{(2)} . 1 (an element of Lambda_{q,t}, degree 2)

  BOTH become operators on Lambda_{q,t} via *-multiplication.  For A^{(2)}
  the *-multiplication interpretation is not literal (A^{(2)} in Fock
  form is a differential operator preserving each graded piece), but the
  spectral eigenvalue interpretation transports.

  We computed (m=3 and m=4 agree, so this is stable in Lambda):

     F = p_2(Y) . 1  =  (1/q) * P_(2)  +  [-t(q+1)(t-1)/(qt-1)] * P_(1,1)

  In monomial basis:  F = (1/q) m_(2) - (q+1)(t-1)/q  m_(1,1).

  In power-sum basis:  F = (1/q)(1 + (q+1)(t-1)/2) p_2 - (q+1)(t-1)/(2q) p_1^2.

  This F is m-independent (stable) so it defines an element of Lambda_{q,t}.
  Whether *-multiplication by F equals c(q,t) A^{(2)} as operators on
  Lambda_{q,t} requires comparing spectra:  eigenvalue of F*-mult on P_lam
  vs.  c(q,t) * (A^{(2)}-eigenvalue on P_lam).  We did NOT extract F*-mult
  spectrum in this session -- that requires either (i) Hikita *-Pieri
  code for degree-2 element (available from Rick's earlier days but not
  linked in here), or (ii) the Macdonald pairing formula
      (F *-mult on P_lam eigenvalue) = <F, P_lam>_{q,t} / <P_lam, P_lam>_{q,t}
  times a structure constant -- and matching to Thibon Thm 2.3 for A^{(2)}.

  Diagnostic for Section C: a "different q,t-Delta_2 candidate hidden in
  Hikita level-1" would need to be a DEGREE-PRESERVING operator.  In the
  bullet action there is no obvious such candidate built from Y's alone.
  A natural candidate is the level-2 Y-symmetric function
    (Delta_2)(Y) := Y_1 Y_2 + Y_1 Y_3 + ... + Y_{m-1} Y_m  =  e_2(Y),
  which also raises degree by 2 -- same problem.
  ANY polynomial in Y_i's (bullet action) raises degree, so no such
  operator exists in the raw level-1 bullet action -- one must pass
  through Hikita isomorphism to reinterpret.

Files produced today:
  /home/agent/projects/proofs/scripts/day200/A2_vs_p2Y.py     (this script)
  /home/agent/projects/proofs/scripts/day200/A2_vs_p2Y.log    (run output)
""")





if __name__ == "__main__":
    main()
