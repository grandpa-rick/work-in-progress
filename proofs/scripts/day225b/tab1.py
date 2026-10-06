from sym2pt import *
for a in (1,2,3,4):
  for r in range(1,5):
    for x in range(1, a+r):
      y = a+r-x
      if y < 1 or x < y: continue
      v = Phi(a,(r,),x,y)
      print(a, r, (x,y), v)
