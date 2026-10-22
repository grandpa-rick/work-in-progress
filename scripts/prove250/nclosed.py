import sys
sys.argv=[sys.argv[0],sys.argv[1]]
exec(open('chk1.py').read().split('# (i)')[0])
bad=0
for L in subsets(n):
    u={k:(1/w0[k] if k in L else w0[k]) for k in w0}
    E=Fr(1)
    for k in range(1,n+1):
        for l in range(k,n+1): E*=1-1/(u[k]*u[l])
    if N(L,w0)*E!=1: bad+=1
print('n',n,'N closed form mismatches',bad)
