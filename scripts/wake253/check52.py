# independent Demazure (type C) check of the sign-exception value H_{52}(-5/2) on LG(5,10)
import sys
from fractions import Fraction as Fr
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/prove249d'); from dem import chi
for k in [(5, 2), (4, 2), (3, 2)]:
    vals = [chi(set(k), 5, j)[0] for j in range(sum(k) + 1)]; n = len(vals)
    def H(x):
        s = Fr(0)
        for i in range(n):
            t = Fr(vals[i])
            for l in range(n):
                if l != i: t *= Fr(x - l, i - l)
            s += t
        return s
    print(k, 'Demazure H(-5/2) =', H(Fr(-5, 2)), ' H(-3/2) =', H(Fr(-3, 2)))
