import sys, itertools
sys.path.insert(0,'../prove249d'); from dem import Di, word
def strict(m): return [tuple(sorted(c,reverse=True)) for L in range(m+1) for c in itertools.combinations(range(1,m+1),L)]
def n(k,j,m):
    f={tuple([-j]*m):1}
    for i in reversed(word(set(k),m)): f=Di(f,i,m)
    s=(-1)**sum(k)
    return {e:s*c for e,c in f.items()}
if __name__=="__main__":
    m=int(sys.argv[1])
    for k in strict(m):
        for j in range(1,4):
            d=n(k,j,m)
            sh={tuple(e[t]+j for t in range(m)):c for e,c in d.items()}
            print(k,j,sum(d.values()),sorted(sh.items()))
