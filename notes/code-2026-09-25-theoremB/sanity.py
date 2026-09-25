from check_theoremB import *
n,k=5,2; st=setup(n,k); A=[gen_a(n,st,i) for i in range(n)]
H=[h_op(n,st,A,r) for r in range(n)]; P=[None]+[p_heap(n,st,A,e) for e in range(1,n)]
S0=frozenset({0,1})  # empty partition: beads at 0..k-1
print('h2 on empty:', H[2].M[S0]); print('h3 on empty:', H[3].M[S0])
print('p2 on empty:', P[2].M[S0]); print('p4 on (3,3)=beads{3,4}:', P[4].M[frozenset({3,4})])
print('nnz P:', [sum(len(P[e].M[S]) for S in st) for e in range(1,n)])
D=Op(st)
for e in range(1,n): D=D+P[e].after(H[n-e])
print('D ==  wrong-sign?', D==Op.identity(st,{1:-(-1)**(k-1)*(n-k)}), ' D==right?', D==Op.identity(st,{1:(-1)**(k-1)*(n-k)}))
# rigidity: t generic breaks commutation? use t=2 numerically
