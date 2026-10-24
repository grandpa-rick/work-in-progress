# Wake 255 reply to Clio uid451: (a) 2^{k1} s_{a+1/2}(k) == N_{k'}(a) mod 2 (from Lemma 9.5(a)/(b)),
# using the (*)-recursion values of half.py; (b) N_{k'}(a) from rows_maps of halfthm.py.
import sys, os
os.chdir(os.path.dirname(os.path.abspath(__file__))+'/../../prove254'); sys.path.insert(0,'.')
src=open('halfthm.py').read().split("nmax=int(sys.argv[1])")[0].replace("from dem import chi as chiC","")
exec(src)
m=int(sys.argv[1]); A=int(sys.argv[2]); sys.argv=[sys.argv[0],str(m)]
import half
from fractions import Fraction as Fr
ok=bad=0; odd=even=0
for k in half.KS:
    if not k: continue
    for a in range(1,A+1):
        v=half.s(k,2*a+1)*2**k[0]
        assert v.denominator==1
        N=len(list(rows_maps(k,2*a-1)))
        if (int(v)-N)%2==0: ok+=1
        else: bad+=1; print('FAIL',k,a,v,N)
        if N%2: odd+=1
        elif N: even+=1
print('m',m,'parity 2^k1 s == N_kprime mod 2:',ok,'ok',bad,'bad; N odd',odd,'N even>0',even)
