import sys; sys.path.insert(0,'../day224')
from ops import *
def hall(f, g): return sum(v*g[k]*zee(k) for k, v in f.items() if k in g)
def phi_engine(a, nu, x, y, t):
    f = {tuple(nu): Fr(1)}
    return hall(T(a, f, t), {tuple(sorted((x, y), reverse=True)): Fr(1)})
