import sys; sys.path.insert(0,'../day242/q484')
from engine import *
S=subsets(2)
for I in S:
  for M in S:
    v,w=W('C',2,I,M)
    print('I',I,'M',M,':',v)
