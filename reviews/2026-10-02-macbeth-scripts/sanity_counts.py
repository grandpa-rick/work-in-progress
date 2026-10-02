from braces import *
import time
for mods in [(8,),(2,4),(16,),(2,8)]:
    A=AbGroup(mods); t=time.time()
    R=regular_subgroups(A)
    bad=sum(not check_brace(A,l,brace_ops(A,l)[0]) for l in R)
    print(mods,len(A.automorphisms()),len(R),'bad',bad,round(time.time()-t,1))
