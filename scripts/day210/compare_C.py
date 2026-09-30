import pickle, sys
from fit_Pell import exact_compare
res = pickle.load(open(sys.argv[1], 'rb')); ell = int(sys.argv[2])
ok = all(exact_compare(ell, k, res[k]) for k in sorted(res) if k > 0)
print('ALL OK' if ok else 'FAIL')
