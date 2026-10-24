# HZ (sum_g Phi, Lemma 9.5 closed forms) vs (*)-recursion values (half.py), larger m. Checks the algebra of Sec. 9.
import sys
src=open('halfthm.py').read().split("nmax=int(sys.argv[1])")[0].replace("from dem import chi as chiC","")
exec(src)
m=int(sys.argv[1]); A=int(sys.argv[2])
sys.argv=[sys.argv[0],str(m)]
import half
ok=bad=0
for k in half.KS:
    if not k: continue
    n1=k[0]
    for a in range(1,A+1):
        M=2*a-1
        tot=sum((Phi_pos(n1,M,u) if M>=n1-1 else Fr((-1)**n1,2**n1)) for u in rows_maps(k,M))
        if tot==half.s(k,2*a+1): ok+=1
        else: bad+=1; print('FAIL',k,a)
print('m',m,'HZ closed form vs (*):',ok,'ok',bad,'bad')
