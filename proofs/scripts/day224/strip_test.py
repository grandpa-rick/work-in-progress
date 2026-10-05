import re, sympy
t = sympy.symbols('t')
L = {}
for f in ['../wake222/gprime/analyze_n2-5.log', '../wake223/analyze_n6.log']:
    for line in open(f):
        m = re.match(r'(\([\d, ]*\)) -> (\([\d, ]*\)):.* v=(\d+) .*Lead=(.*)$', line.strip())
        if m: L[(eval(m.group(1)), eval(m.group(2)))] = (int(m.group(3)), sympy.sympify(m.group(4)))
print(len(L))
ok = bad = miss = 0
for (lam, mu), (v, ld) in sorted(L.items()):
    l = len(lam)
    if len(mu) != l or lam == mu: continue
    lt = tuple(x-1 for x in lam if x > 1); mt = tuple(x-1 for x in mu if x > 1)
    if (lt, mt) not in L: miss += 1; print('missing', lam, mu, '->', lt, mt, v, sympy.factor(ld)); continue
    v2, ld2 = L[(lt, mt)]
    good = sympy.expand(ld-ld2) == 0 and v == v2
    ok += good; bad += not good
    print(lam, mu, '->', lt, mt, 'OK' if good else f'FAIL {sympy.factor(ld)} vs {sympy.factor(ld2)} v={v},{v2}')
print('ok', ok, 'bad', bad, 'missing', miss)
