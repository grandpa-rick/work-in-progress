#!/usr/bin/env python3
"""Points strictly inside Korff's stated domain: n_K=2k>2 (k>=2) and #variables n <= k_K = w = 2l+2.
alpha=(1^n): every strip non-constant, nu(P)=n; report ord_{t=1} c_alpha vs n."""
from warnaar_calpha import calpha, ord1
for k in (2,3):
    for l in range(0,4):
        w=2*l+2
        for n in range(2, w+1):
            c=calpha(k,l,[1]*n); o=ord1(c)
            print(f"k={k} l={l} w={w} n={n}: ord={o} required>={n}  {'FAIL(identity impossible)' if o is not None and o<n else ('vacuous c=0' if o is None else 'not informative')}")
