"""Clio 416 sec 5: Rick's three Theorem-Psi sec-6 'hunches' as identities. With near-miss controls."""
import sympy as sp
v, q, x = sp.symbols('v q x')
psi = lambda v, q: (1 - q*v**2)/(1 - q)
print('bar: psi(v;1/q) - v^2 psi(1/v;q) =', sp.simplify(psi(v, 1/q) - v**2*psi(1/v, q)))
print('  control (1-q v^3)/(1-q):', sp.simplify(((1-v**3/q)/(1-1/q)) - v**2*((1-q/v**3)/(1-q))) != 0)
print('q->oo limit:', sp.limit(psi(v, q), q, sp.oo))
chi = lambda t: t**2
print('chi=v^2: chi(x)chi(1/x)=', sp.simplify(chi(x)*chi(1/x)), ' chi(1)=', chi(1), ' chi(-1)=', chi(-1))
print('  control chi=v: chi(x)chi(1/x)=', sp.simplify(x*(1/x)), ' chi(-1)=', -1, '(rejected)')
print('bar(psi) = chi*psi with chi=v^2 evaluated at 1/v? psi(v;1/q)/psi(1/v;q) =', sp.simplify(psi(v,1/q)/psi(1/v,q)))
# Which chi does bar actually produce?  Full bar (q->1/q, v->1/v) vs q-only bar.
full = sp.simplify(psi(1/v, 1/q)/psi(v, q)); qonly = sp.factor(sp.simplify(psi(v, 1/q)/psi(v, q)))
print('full bar  psi(1/v;1/q)/psi(v;q) =', full)
print('q-only bar psi(v;1/q)/psi(v;q) =', qonly, '; chi(x)chi(1/x)=', sp.simplify(qonly*qonly.subs(v, 1/v)), ' chi(+-1)=', qonly.subs(v,1), qonly.subs(v,-1))
