def Psi1(s,t,S,w):
    """Psi_{(s),t}(S) = sum_{m in S, m>=s} (1-x^-2) prod_{k<s}(1-w_k/x) prod_{u in S\\m} (1-1/(w_u x))/(1-w_u/x) * x^t, x=w_m"""
    tot=0
    for m in S:
        if m<s: continue
        x=w[m]; e=(1-x**-2)*x**t
        for k in range(1,s): e*=(1-w[k]/x)
        for u in S:
            if u!=m: e*=(1-1/(w[u]*x))/(1-w[u]/x)
        tot+=e
    return tot
