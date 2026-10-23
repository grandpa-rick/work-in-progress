# negative control: the LG (Buch-Ravikumar Cor 5.4) Pieri rule must FAIL against the true OG K-class.
import sys, random
sys.argv = ['og.py', sys.argv[1]]
import og
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake252/caseB'); import brpieri
n = int(sys.argv[1]); random.seed(1); LAMS = []
for _ in range(6 * 2 ** n):
    par = random.randint(0, 1); LAMS.append(tuple(2 * random.randint(-og.RR, og.RR) + par for _ in range(n)))
cache = {}; KS = [k for k in og.strict_parts(n) if k]
print(f'm={n}: LG rule == true OG class: {sum(og.true_kclass(k, n, LAMS, cache) == brpieri.pieri_dim(k, n) for k in KS)}/{len(KS)} (should be small)')
