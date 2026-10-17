import sys
from num4 import run
for seed in (11,23):
  for a,b in [(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]:
    out=run(4,a,b,seed=seed)
    print(seed,(a,b),out[((),(3,),(1,2,4))],flush=True)
