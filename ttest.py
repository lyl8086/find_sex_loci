#/bin/python3
import sys
import numpy as np 
from scipy.stats import ttest_ind

d = sys.argv[1]
a = int(sys.argv[2])
b = int(sys.argv[3])
if d == '-':
    d = sys.stdin
# scaffold1	257	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0
# the python slice include the start and not the end, i.e. [,)
gp1 = 2+a
gp2 = 2+a+b
for line in d:
    p    = line.split("\t")
    a1   = list(map(int,p[2:gp1]))
    a2   = list(map(int,p[gp1:gp2]))
    if np.sum(a1)<=60 and np.sum(a2)<=60:
        continue
    b1,b2,b3,b4 = np.mean(a1),np.std(a1,ddof = 1), np.mean(a2),np.std(a2,ddof = 1) 
    t, P = ttest_ind(a1, a2, equal_var=False)
    if P > 0.05:
        continue
    print(p[0],p[1],P,b1,b2,b3,b4,sep="\t")
    