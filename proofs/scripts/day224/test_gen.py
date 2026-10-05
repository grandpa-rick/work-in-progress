from second_gen import *
t = Fr(3)
for lam in [(2,2,2), (3,2,1)]:
    A = second_general(lam, t)[2]; B = second(*lam, t)
    print(lam, A == B)
