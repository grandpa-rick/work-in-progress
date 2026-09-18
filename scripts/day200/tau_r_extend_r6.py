"""
Day 200: try to extend Day 198 SymPy data to r=6 (m=8).
If it runs in reasonable time, we get an independent verification of the r=6
prediction from the Day 200 closed form.
"""
import sys
sys.path.insert(0, '/home/agent/projects/proofs/scripts/day198')
from p2Y_er_extended import run_case
import sympy as sp

q, t = sp.symbols('q t')

# Predicted r=6 tau from Day 200:
tau_r6_predicted = -(q**2 - 1) * (t**2 + 1) * (t**4 + 1) * (q*t**7 - q + t + 1) / q**3

print("Day 200 predicted tau_6(q,t):")
print(f"  {sp.factor(tau_r6_predicted)}")
print()
print("Running run_case(r=6, m=8) ...")
factored = run_case(r=6, m=8)
# The top-row coefficient is factored[(8,)]
top = factored.get((8,), sp.Integer(0))
print()
print(f"[EXTRACTED] coefficient of e_(8) = tau_6 (SymPy) = {top}")
diff = sp.simplify(sp.expand(top - tau_r6_predicted))
print(f"[CHECK] Day 200 closed form matches Day 198-style direct SymPy?  diff = {diff}   {'OK' if diff == 0 else 'FAIL'}")
