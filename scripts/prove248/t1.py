from psi import *
# (1)-identity check
print(Psi((1,),(0,),(1,2,3)))
for M in [(2,3),(1,2,3),(2,3,4),(1,2,3,4)]:
    print(M, Psi((1,2),(0,0),M))
