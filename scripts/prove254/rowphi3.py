from rowphi import *
import itertools
def classify(cond,name,Moffs=(0,1,2)):
  for Moff in Moffs:
    neg=0; tot=0; zero=0; ex=[]
    for n in range(1,6):
      M=n+Moff
      for k in range(0,n):
        for u in itertools.product(range(1,M+2),repeat=k):
          if not cond(u): continue
          p=Phi(n,M,u); tot+=1
          if p<0:
              neg+=1
              if len(ex)<3: ex.append((n,M,u,str(p)))
          if p==0: zero+=1
    print(name,'M-n=',Moff,'tested',tot,'neg',neg,'zero',zero,ex)
inc=lambda u: all(u[i]<u[i+1] for i in range(len(u)-1))
winc=lambda u: all(u[i]<=u[i+1] for i in range(len(u)-1))
classify(lambda u: inc(u) and (not u or u[0]>=2),'strict, u1>=2')
classify(lambda u: inc(u) and (not u or u[0]>=2) and (not u or u[0]%2==0),'strict, u1 even')
classify(lambda u: winc(u) and all(u[b]>=b+2 for b in range(len(u))),'weak, u_b>=b+1 (1-idx)')
classify(lambda u: all(u[b]>=b+2 for b in range(len(u))),'any, u_b>=b+1')
classify(lambda u: all(u[b]>=b+1 for b in range(len(u))),'any, u_b>=b (1-idx: u_b>=b+1)')
classify(lambda u: winc(u) and all(u[b]>=b+1 for b in range(len(u))),'weak, u_b>=b+1')
classify(lambda u: winc(u) and (not u or u[0]>=2),'weak, u1>=2')
