import sys; sys.path.insert(0,"../prove247b"); from scert5 import *
for x in [((),(1,2),(1,5,6)),((1,),(1,2),(1,5,6)),((1,2),(2,),(1,5,6)),((1,2),(1,2),(1,5,6))]:
    print(x,start('C',6,*x,None,None),PSTAT)
