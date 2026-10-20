"""Cross-check: exact flint residue vs prove247/mp.py mod-p residue at mp's random point, plus P sign change witness."""
import sys,random; sys.path.insert(0,'../prove247')
import mp, eng
from fractions import Fraction
N=4; eng.setup(N)
for key,c in [((1,2,(),(),(1,3,4)),5),((1,2,(),(3,),(1,2,4)),0),((3,4,(1,),(2,),(1,2,3,4)),6),((1,2,(),(),(1,2,4)),5),((1,2,(),(3,),(1,2,4)),1)]:
    a,b,l,m,k=key
    s=mp.Sys(N,1,a,b); vm=mp.A(s,a,b,l,m,k,c)
    RR=random.Random(1); vals=[RR.randrange(2,mp.P-1) for _ in range(N)]
    r=eng.A(eng.Sys(N,a,b),a,b,l,m,k,c)
    def evp(poly):
        tot=0
        for e,co in poly.to_dict().items():
            co=Fraction(str(co)); t=co.numerator*mp.inv(co.denominator)
            for i,ei in enumerate(e): t=t*pow(vals[i],ei,mp.P)
            tot+=t
        return tot%mp.P
    ve=evp(r.n)*mp.inv(evp(r.d))%mp.P
    print(key,c,'modp',vm,'exact',ve,'MATCH' if vm==ve else 'MISMATCH')
