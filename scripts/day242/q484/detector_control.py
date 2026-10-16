import sympy as sp
from engine import is_laurent, q
w = sp.symbols('w1:3')
tests = [(1/(w[0]*w[1]-1), False), ((q+w[1])/w[0]**3, True), ((w[0]*w[1]-1)/((w[0]*w[1]-1)*w[0]), True),
         ((w[0]+1)/(q-w[0]*w[1]), False), (1/(1-q), True), (sp.Integer(0), True)]
print('Laurent-detector control:', sum(is_laurent(e, w) == t for e, t in tests), '/', len(tests))
