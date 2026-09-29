"""Day 209: closing identity (★2) for (TC), reduced to a rational identity.
A=t^i, B=t^j symbolic; X,W GF variables with rho=W/X=z/w.
All C's divided by D_i D_j X^i W^j G(t^2X) G(t^2W)."""
import sympy as sp
s,t,X,W,A,B=sp.symbols('s t X W A B')
rho=W/X
def Chat(Y,a):   # C_i(Y)/(D_i Y^i G(t^2 Y)) with a=t^i, for Y=X (unshifted arg)
    return (1-s*t*Y)*(1-s*Y/a)/((1-t*Y)*(1-Y))
def Chat_t(Y,a): # C_i(tY)/(D_i Y^i G(t^2Y))
    return a*(1-s*t*Y/a)/(1-t*Y)
# C_{i-1}(tY)/(D_i Y^i G(t^2 Y)) = D_{i-1}/D_i * (a/t) Y^{-1} (1-s t^2 Y/a)/(1-tY)
def Cm1_t(Y,a):
    return (1-a)/(t*s-a) * (a/t)/Y * (1-s*t**2*Y/a)/(1-t*Y)
def K(a,b):  # symbolic K would need products; use ratios only
    pass
kg=lambda a,b:(a-s)*(a*rho-s)/(a*(a*rho-b))
kd=lambda a,b:(b-s*rho)*(b-s)/(b*(b-a*rho))
# K_{i-1,j}/K_{ij}, K_{i,j-1}/K_{ij}
rKi = B*(A/t*rho-B)*(A/t*rho-1/t)/((A/t*rho-s)*(A/t*rho-B/t))
# K_{i,j-1}/K_ij: by symmetry z<->w (rho->1/rho, swap roles). derive directly below
# K_ij = prod_{p<i}(t^p rho - s)/(t^p rho - B) * prod_{r<j}(s rho - t^r)/(A rho - t^r)
# K_{i,j-1}/K_ij = prod_{p<i}(t^p rho - B)/(t^p rho - B/t) * (A rho - t^{j-1})/(s rho - t^{j-1})
# prod_{p<i}(t^p rho - B)/(t^p rho - B/t) = prod t^{-1}... telescoping: (t^p rho - B) = t (t^{p-1} rho - B/t)
#   => prod_{p<i} t (t^{p-1}rho - B/t)/(t^p rho - B/t) = t^i (rho/t - B/t)/(t^{i-1} rho - B/t)
rKj = A*(rho/t-B/t)/(A/t*rho-B/t) * (A*rho-B/t)/(s*rho-B/t)
LHS = Chat(X,A)*Chat(W,B) - Chat_t(X,A)*Chat_t(W,B)
RHS = Chat(X,A)*Chat(W,B)*(kg(A,B)*B*X/(s*(1-s*X/A)) + kd(A,B)*A*W/(s*(1-s*W/B))) \
    - rKi*s/B**2*Cm1_t(X,A)*Chat_t(W,B)*kg(A/t,B)*t*B*X/(s*(1-s*t**2*X/A)) \
    - rKj*s/A**2*Chat_t(X,A)*Cm1_t(W,B)*kd(A,B/t)*t*A*W/(s*(1-s*t**2*W/B))
print(sp.factor(sp.together(LHS-RHS)))
