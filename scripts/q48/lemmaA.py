"""Q48 checks. (A) down(z_rho, lam) == z_lam  iff  lam/rho antichain (rho ⊆ lam).
(P) P^down_{z_rho,z_kappa} == { (kappa\\O ; kappa) : V-map o-sets O on kappa/rho }, with weights, all kappa,rho.
(C) negative control: untwisted dictionary breaks (P)."""
import sys
src=open('../q47/ls12_cor34_check.py').read().split('allok = True')[0]
sys.argv=['x',sys.argv[1]]+sys.argv[2:]
exec(src)
cp=lambda k: tuple(-x for x in v_of(k)) if not CTRL else v_of(k)
zz={k:minrep[cp(k)] for k in strict}
part=lambda c: [k for k in strict if cp(k)==c][0]
def sdown(w,c):
  try: return down(w,c)
  except ValueError: return None
def sdownp(w,p):
  for c in p[0]:
    w=sdown(w,c)
    if w is None: return None
  return w
def anti(S): return not any(x!=y and x[0]<=y[0] and x[1]<=y[1] for x in S for y in S)
badA=0
for r in strict:
  for l in strict:
    if not D(r)<=D(l): continue
    if (sdown(zz[r],cp(l))==zz[l]) != anti(D(l)-D(r)): badA+=1
print('m',m,'CTRL',CTRL,'Lemma A failures:',badA)
badP=0; tot=0
for k in strict:
  for r in strict:
    if not D(r)<=D(k): continue
    ps=set()
    for p in paths:
      if cleq(p[0][0],coset[zz[r]]) and sdownp(zz[r],p)==zz[k]:
        ps.add(tuple(part(c) for c in p[0]))
    th=D(k)-D(r); pred=set()
    for o in itertools.product([0,1],repeat=len(th)):
      O={x for x,b in zip(sorted(th),o) if b}; A=th-O
      if any(x[0]==x[1] for x in O): continue
      if not(anti(O) and anti(A)): continue
      if any(x!=y and x[0]<=y[0] and x[1]<=y[1] for x in O for y in A): continue  # A must be an order ideal of theta
      lam=D(k)-O
      lp=[l for l in strict if D(l)==lam][0]
      pred.add((k,) if not O else (lp,k))
    tot+=len(pred)
    if ps!=pred: badP+=1
print('  Prop P failures:',badP,'of pairs; total V-maps',tot)
