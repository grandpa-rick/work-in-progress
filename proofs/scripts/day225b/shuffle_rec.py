# symbolic check (X=t^A, Y=t^B independent) that closed Sh satisfies the last-letter recursion. grade: computed (identity in Q(t,w,X,Y))
import sympy as sp
t, w, X, Y = sp.symbols('t w X Y')
def s(X, Y, pref):  # Sh/qbin with t^{-AB} replaced by pref
    return pref*(1-w)*(1-w*Y/X)/((1-w/X)*(1-Y*w))
lhs = s(X, Y, 1)
rhs = (1-w*t/X)/(1-w*Y*t/X)*(1-X)/(1-X*Y)*s(X/t, Y, Y) + (1/X)*(1-w*Y/t)/(1-w*Y/(t*X))*(1-Y)/(1-X*Y)*s(X, Y/t, X)
print('recursion identity holds:', sp.simplify(lhs - rhs) == 0)
