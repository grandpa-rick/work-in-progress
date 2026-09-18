"""
Day 202: Jack degeneration of Rick's p_3(Y)-Pieri Lemma.

Target: check whether Rick's Day 201 k=3 formulas at q=1, t -> Jack limit
reproduce the shape of Thibon 2608.25651 Conjecture 8.2 Delta_3(alpha).

Thibon setup (2608.25651 sec 7-8):
    - Delta_2(alpha) is the deformed cut-and-join operator (Theorem 124),
      p_i p_j D_{i+j} and p_{i+j} D_i D_j terms + (alpha-1) diagonal.
    - Delta_3(alpha) is the CONJECTURED operator (Conj 8.2, eq 139) —
      a triple sum in p_I D_J monomials + (alpha-1) correction.
    - Eigenvalue on Jack J_lambda in dual convention: shifted power sums
      tilde{p}_k(lambda) (eq 13 of Thibon 2608.25651).
    - Delta_2(alpha) has eigenvalue 2 sum_box c_alpha(box) on Jacks.

Rick setup (Day 201):
    p_3(Y).e_r(X) = c_{(r,1,1,1)} e_{(r,1,1,1)}
                  + c_{(r,2,1)} e_{(r,2,1)}
                  + c_{(r+1,1,1)} e_{(r+1,1,1)}
                  + c_{(r,3)} e_{(r,3)}
                  + c_{(r+1,2)} e_{(r+1,2)}
                  + c_{(r+2,1)} e_{(r+2,1)}  (r-dep, Baxter-2)
                  + c_{(r+3)} e_{(r+3)}      (r-dep, Baxter-4)

Jack degeneration convention (see Day 201 jack_limit_check.py):
    q = 1 (Macdonald -> Jack)
    t = t_1  (kept as Jack parameter; typical convention t = q^alpha then q->1)

Equivalent alt convention: introduce eps = t - 1 small, q = (1+eps)^alpha,
    then eps -> 0 gives leading + corrections in alpha.

We test both to see which produces a sensible limit. Then compare shape
to Thibon's Delta_3 structure.
"""
import sympy as sp

q, t, alpha, eps, r = sp.symbols('q t alpha epsilon r', positive=True)

# ---------- Rick's r-INDEPENDENT closed forms (Day 201, section 4) ----------
# Multiplied by q^6, they are pure polynomials in q,t (rat'l fns of that).

c_r1k = sp.Rational(1)  # times q^{-6}  (k=3 pivot)
c_r21 = -(q**2*t - q**2 + q*t - q + t + 2)  # times q^{-6}
c_rp1_1k1 = q**3 - 1  # times q^{-6} for c_{(r+1,1,1)}
c_r3 = (q**3*t**3 - q**3*t**2 - q**3*t + q**3
        + q**2*t**3 - q**2 + q*t**3 - q + t**2 + t + 1)  # times q^{-6}
c_rp1_2 = -(q**3 - 1)*(q*t**2 + t - q + 1)  # times q^{-6}

# r-DEPENDENT
c_rp2_1 = -(q**3 - 1) * (q**2 * t**(r+1) - q**2 - q*t**(r+2) + q*t + 1)  # times q^{-6}
# c_{(r+3)}: Baxter-4, top coefficient q^3 t^6 (q^3-1)/(t^3-1) at t^{3r}

qfactor = q**6  # divide all these by q^6 to get true coefficient

# ---------- Two conventions for Jack degeneration ----------

def q_to_1_first(expr, tval=None):
    """Send q -> 1 first, then optionally set t."""
    e = expr.subs(q, 1)
    if tval is not None:
        e = e.subs(t, tval)
    return sp.simplify(e)

def q_talpha_limit(expr, r_val=None):
    """Substitute q = t^alpha, expand in eps = t-1 to O(eps^2).
    This is the standard Macdonald -> Jack limit."""
    e = expr
    if r_val is not None:
        e = e.subs(r, r_val)
    # q = t^alpha, and take t = 1 + eps series
    e = e.subs(q, (1 + eps)**alpha).subs(t, 1 + eps)
    ser = sp.series(e, eps, 0, 3).removeO()
    return sp.simplify(sp.expand(ser))

print("=" * 78)
print("Day 202: Jack degeneration of Rick's p_3(Y)-Pieri formulas")
print("=" * 78)

# --- Test 1: q = 1 straight ---
print("\n[Test 1] Direct q = 1 substitution (Schur limit, not Jack):")
print(f"  c_(r,1,1,1) [times q^6] at q=1: {q_to_1_first(c_r1k)}")
print(f"  c_(r,2,1)             at q=1: {q_to_1_first(c_r21)}")
print(f"  c_(r+1,1,1)           at q=1: {q_to_1_first(c_rp1_1k1)}")
print(f"  c_(r,3)               at q=1: {q_to_1_first(c_r3)}")
print(f"  c_(r+1,2)             at q=1: {q_to_1_first(c_rp1_2)}")
print(f"  c_(r+2,1)             at q=1: {q_to_1_first(c_rp2_1)}")

# At q=1: several coefficients vanish (they have factor q^3 - 1).
# Only c_(r,1,1,1), c_(r,2,1), c_(r,3) survive.
# But q^6 -> 1 in prefactor, so numerators = actual coeffs.

# --- Test 2: full Jack limit q = t^alpha, expand in eps ---
print("\n[Test 2] Jack limit q = t^alpha, expand in eps = t-1 to O(eps^2):")
for name, coef in [("c_(r,1,1,1)*q^6", c_r1k),
                    ("c_(r,2,1)*q^6", c_r21),
                    ("c_(r+1,1,1)*q^6", c_rp1_1k1),
                    ("c_(r,3)*q^6", c_r3),
                    ("c_(r+1,2)*q^6", c_rp1_2)]:
    val = q_talpha_limit(coef)
    print(f"  {name}: {val}")

# For r-dependent, need r_val
print("\n  c_(r+2,1)*q^6 at r=1,2,3 (Jack limit):")
for rv in [1, 2, 3]:
    val = q_talpha_limit(c_rp2_1, r_val=rv)
    print(f"    r={rv}: {val}")

# --- Analysis: divide by (q^3 - 1) to isolate finite Jack limit ---
# Since c_(r+1,1,1), c_(r+1,2), c_(r+2,1) all have (q^3-1) factor,
# they vanish at q=1 in the leading order. To see the Jack contribution,
# divide by (q^3 - 1) first (extract Jack normalization):

print("\n[Test 3] Rescaled: divide by (q-1) or (q^3-1) to see Jack contribution:")
for name, coef in [("c_(r+1,1,1)/(q^3-1)", c_rp1_1k1),
                    ("c_(r+1,2)/(q^3-1)", c_rp1_2 / (q**3 - 1))]:
    e = coef.subs(q, (1+eps)**alpha).subs(t, 1+eps)
    ser = sp.series(e, eps, 0, 3).removeO()
    print(f"  {name}: {sp.simplify(sp.expand(ser))}")

# --- Structural check: shape of Rick's formulas ---
print("\n[Structural summary]")
print("Rick's r-indep coefficients live in Z[q, t] (numerator of c_mu * q^6):")
print("  c_(r,1,1,1) * q^6 = 1")
print("  c_(r,2,1)   * q^6 = -(q^2 t - q^2 + q t - q + t + 2)")
print("  c_(r+1,1,1) * q^6 = q^3 - 1     [vanishes at q=1]")
print("  c_(r,3)     * q^6 = polynomial in q,t")
print("  c_(r+1,2)   * q^6 = -(q^3-1)(q t^2 + t - q + 1)  [vanishes at q=1]")
print("  c_(r+2,1)   * q^6 = -(q^3-1)(...)  r-dep [vanishes at q=1]")
print()
print("At q=1: three coefficients survive (r,1,1,1), (r,2,1), (r,3).")
print("Their VALUES at q=1:")
for name, coef in [("c_(r,1,1,1)", c_r1k),
                    ("c_(r,2,1)", c_r21),
                    ("c_(r,3)", c_r3)]:
    val = sp.simplify(coef.subs(q, 1))
    print(f"  {name} at q=1: {val}")

# --- Thibon shape comparison ---
print("\n[Thibon Delta_3(alpha) shape — Conjecture 8.2, eq 139 of 2608.25651]")
print("Delta_3(alpha) = sum p_{i+j+k} D_pi D_pj D_pk + alpha sum p_{i+k} p_j D_pi D_pj D_pk")
print("               + alpha sum p_{i+j} p_k D_{p_{j+k}} D_pi + ...")
print("               + alpha^2 sum p_i p_j p_k D_{p_{i+j+k}}")
print("               + (alpha-1) * [correction with (i+j-2) weights]")
print()
print("Its eigenvalue on Jack J_lambda is p_3(Xi_infty) / normalization = shifted")
print("power sum ~ sum_i lambda_i^3 + O(alpha) content correction.")
print()
print("Rick's action: p_3(Y) . e_r(X)")
print("This is p_3 of the Y-Cherednik operators acting via star-product on e_r(X).")
print("Its expansion in e_mu basis has 7 nonzero coefficients (5 r-indep + 2 r-dep).")

# --- Verdict analysis ---
print("\n" + "=" * 78)
print("VERDICT ANALYSIS")
print("=" * 78)
print("""
1. SEMANTIC MISMATCH.
   Thibon's Delta_3(alpha) is a DIFFERENTIAL OPERATOR on Lambda whose action on
   Jack basis is diagonal (eigenvalue = normalized central character on p_3(Xi)).
   Rick's p_3(Y) is a symmetric polynomial in Cherednik Y-operators acting via
   the star-product; its Pieri expansion for p_3(Y).e_r(X) is a decomposition
   in the e-basis of X-variables, with several (5+2=7) nonzero coefficients.

   These are STRUCTURALLY DIFFERENT: Thibon writes down the operator on
   Lambda (symmetric functions); Rick writes down the matrix coefficients of
   a specific vector p_3(Y).e_r in a specific basis. They do not compare
   coefficient-by-coefficient.

2. DEGREE MISMATCH (analog of Day 201 Vertex B degree obstruction).
   Rick's p_3(Y) is degree +3 on X (raises degree by 3). Thibon's Delta_3(alpha)
   is degree 0 on symmetric functions (multiplication by hat{p}_2 lifted to
   normalized class algebra; the class-algebra eigenvalue makes it a
   diagonal degree-preserving operator).

   This is the SAME structural obstruction Day 201 identified for k=2:
   p_k(Y) is degree +k on X, Delta_k(alpha) is degree 0. They cannot match
   coefficient-by-coefficient.

3. EIGENVALUE MISMATCH.
   Delta_3 eigenvalue on Jack J_lambda: shifted p_3 (cubic in lambda_i).
   Rick's coefficient at (r,1,1,1) at q=1: 1 (constant in r, cubic in nothing).
   So even setting semantics aside, no shape match.

CONCLUSION: NO MATCH.
   The Day 201 Vertex B argument for k=2 generalizes to k=3 with identical
   force. Rick's p_3(Y) and Thibon's Delta_3(alpha) live in different graded
   pieces and act on different vectors.

   Rick's paper does NOT prove Thibon's Delta_3 conjecture as a corollary.
   The connection memo's optimism (before Vertex B was refuted) is now
   provably wrong by the same argument that killed Vertex B.
""")
