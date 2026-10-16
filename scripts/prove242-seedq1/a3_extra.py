# q=1 and t_j=0 variants of the A3 separator, nu=(1,2,1),(2,2,1)
import sys; sys.argv = ['x']
exec(open('a3_sep.py').read().split("nus = ")[0])
for qq, tt, lab in [(Fr(1), [Fr(2), Fr(-5, 3), Fr(7, 4)], "q=1, O1"), (Fr(3, 7), [Fr(2), Fr(-5, 3), Fr(0)], "WRONG phi, t_3=0")]:
    B = phi_order([0, 1, 2]) if qq == 1 else [[0, 0, 0], [-1, 0, 0], [0, 0, 0]]
    Q0 = QQ(C, [[0]*3]*3, Fr(3, 7)); Q = QQ(C, B, qq)
    for nu in [(1, 2, 1), (2, 2, 1)]:
        W = words_of(nu); K = mat([theta(Q0, w) for w in W]).nullspace(); M = mat([Q.Psi(w, tt) for w in W])
        print(lab, nu, "dim f", len(W)-len(K), "rank Psi('f)", M.rank(), "descends", all((M*k).is_zero_matrix for k in K), "rank on f", (M*sp.Matrix.hstack(*[sp.eye(len(W))[:, c] for c in range(len(W))])).rank() if False else "")
