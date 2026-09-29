import sys, pickle, time
from ek_two_col import *
k = int(sys.argv[1]); N = int(sys.argv[2])
data = {}
for a in range(0, N+1):
    for b in range(a, N+1):
        r, dt = compute(k, a, b)
        data[(a, b)] = r
        print(f'k={k} (a,b)=({a},{b}) m={a+b+k} {dt:.1f}s support={sorted(r, reverse=True)}', flush=True)
        pickle.dump(data, open(f'data_k{k}_N{N}.pkl', 'wb'))
