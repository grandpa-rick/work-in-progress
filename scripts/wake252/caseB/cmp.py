# compare fields of two gate outputs (dicts parsed, order-insensitive)
import re, ast, sys
def norm(f):
    d = {}
    for L in open(f):
        if L.startswith('('):
            k = re.match(r'(\(.*?\))', L).group(1)
            g = lambda key: ast.literal_eval(re.search(key + r" (\[.*?\])", L).group(1))
            K = ast.literal_eval(re.search(r"'?K'? (\{.*?\})", L).group(1))
            nv = ast.literal_eval(re.search(r"nonvanishing terms (\{.*?\})", L).group(1))
            d[k] = (g('zH'), g('zb'), g('need'), K, nv)
    return d
a = norm(sys.argv[1]); b = norm(sys.argv[2]); common = [k for k in b if k in a]
print(f'identical (zH,zb,need,K,nonvanishing) on common kappa: {sum(a[k] == b[k] for k in common)}/{len(common)}')
for k in common:
    if a[k] != b[k]: print('DIFF', k, a[k], b[k])
