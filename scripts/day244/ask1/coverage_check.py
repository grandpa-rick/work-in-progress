# Reproduces Clio 10-17 reply 4.4 coverage list (100,91.7,81.9,72.8,64.8,58.0,52.3) under the reading:
# domain = {(M,I,a<b): M>=I (Motegi order), {a,b} cap M = empty}  (= her nonvanH hypotheses, reply l.363)
# stratum = max I <= b, I = basis label (column) of A_{M,I}, max(empty)=0.
import itertools
def subsets(N): return [c for r in range(N+1) for c in itertools.combinations(range(1,N+1),r)]
def geq(M,I):
    s,r=len(M),len(I)
    return s>=r and all(M[s-r+j]>=I[j] for j in range(r))
for N in range(2,9):
    S=subsets(N); tot=ins=0; tot2=ins2=0
    for a,b in itertools.combinations(range(1,N+1),2):
        for M in S:
            if a in M or b in M: continue
            for I in S:
                if not geq(M,I): continue
                tot+=1; ins+= (max(I) if I else 0)<=b
    print(N, round(100*ins/tot,1))
