import sympy as sp
t,s,y=sp.symbols('t s y')
def gk(f,nb):  # nb: dict curve->neighbor value
  out=[]
  for wt,val in nb.items():
    sbs={'t2':[{t:1},{t:-1}],'s2':[{s:1},{s:-1}],'st':[{t:1/s}]}[wt]
    out.append(all(sp.simplify((f-val).subs(sb))==0 for sb in sbs))
  return out
B24={0:lambda c2:(1+y)*(s*t)**(1-c2)*(s*t*(1+y)-y+s**2*y+t**2*y+s**2*t**2*y**2),
     1:lambda c2:(1+y)*(s*t)**(1-c2)*(s*t*(y+y**2)+1+s**2*y+t**2*y-s**2*t**2*y)}
B23={0:lambda c2:(1-s*t)*(1+y)*t**(1-c2)*(s*(1+y)+t+s*t*y),
     1:lambda c2:(1-s*t)*(1+y)*t**(1-c2)*(1+s*t*(1+y)+s**2*y)}
for c in [0,1,2]:
  for c2 in [2*c-1,2*c]:
    al=2*c-c2
    n24={'st':(1+y/t**2)*(1+y/s**2)*(1+y/(s*t)),'t2':(1+y/t**2)*(1+y*s/t)*(1+y)*s**(2-2*c),'s2':(1+y/s**2)*(1+y*t/s)*(1+y)*t**(2-2*c)}
    n23={'st':0,'t2':(1-s**2)*(1+y/t**2)*(1+y*s/t),'s2':(1-t**2)*(1+y/s**2)*(1+y)*(s/t)**(c2-1)}
    print('ceil',c,'ceil2',c2,'alpha',al,' 24 branch0',gk(B24[0](c2),n24),'branch1',gk(B24[1](c2),n24),
          '| 23 branch0',gk(B23[0](c2),n23),'branch1',gk(B23[1](c2),n23))
