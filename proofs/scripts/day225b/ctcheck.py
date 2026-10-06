# Day 225 PROVE: test adjoint formula
#  <T_a g, p_x p_y> (Hall) = (1-t^x)(1-t^y)/phi_a(t) * CT_{|z1|<..<|za|}[ g(z) Z p_x(1/z) p_y(1/z) prod_{i<j}(z_j-z_i)/(z_j-t z_i) ]
import sys; sys.path.insert(0,'../day224')
import numpy as np, itertools
from fractions import Fraction as Fr
from twopoint_helpers import *
def ct(a, nu, x, y, t, M=40):
    rs = [2.0**i for i in range(a)]
    grids = [r*np.exp(2j*np.pi*np.arange(M)/M) for r in rs]
    Z = np.meshgrid(*grids, indexing='ij')
    g = np.ones_like(Z[0])
    for r in nu: g = g*sum(z**r for z in Z)
    F = g*np.prod(Z, axis=0)*sum(z**(-x) for z in Z)*sum(z**(-y) for z in Z)
    K = np.ones_like(Z[0])
    for i in range(a):
        for j in range(i+1, a): K = K*(Z[j]-Z[i])/(Z[j]-t*Z[i])
    return (F*K).mean()
def phia(a, t):
    p = 1
    for i in range(1, a+1): p *= (1-t**i)
    return p
if __name__ == '__main__':
    for t in (0.3, 0.45):
        for a, nu, x, y in [(1,(2,),2,1),(2,(),1,1),(2,(1,),2,1),(2,(2,),2,2),(2,(1,1),2,2),(2,(2,),3,1),(2,(3,),3,2),(2,(2,1),4,1),(3,(2,),3,2),(3,(2,1),4,2),(3,(3,2),5,3),(2,(2,2),3,3)]:
            v = ct(a, nu, x, y, t)*(1-t**x)*(1-t**y)/phia(a, t)
            e = float(phi_engine(a, nu, x, y, Fr(t).limit_denominator(1000)))
            print(t, a, nu, x, y, round(v.real, 9), round(e, 9), abs(v-e) < 1e-7)
