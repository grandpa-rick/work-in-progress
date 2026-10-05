from lead2 import *
def hall(f, g):  # Hall scalar product in p-basis
    return sum(v*g[k]*zee(k) for k, v in f.items() if k in g)
def phi(a, nu, x, y, t):
    f = {tuple(nu): Fr(1)} if nu else {(): Fr(1)}
    return hall(T(a, f, t), {tuple(sorted((x, y), reverse=True)): Fr(1)})
for a in (1, 2, 3):
    for d in range(0, 4):
        for nu in parts(d):
            n = a + d
            for x in range(1, n):
                y = n - x
                if x < y: continue
                p = interp(lambda t: phi(a, nu, x, y, t) * ( (1-t**x)*(1-t**y) ), 2*n*n+10)
                print(f'a={a} nu={nu} (x,y)=({x},{y}):  phi*(1-t^x)(1-t^y) =', p, flush=True)
