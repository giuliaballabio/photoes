
import numpy as np
from math import e

## original grid logarithmically spaced
r = np.logspace(np.log(0.0300785), np.log(9.97397), num=1113, base=e , dtype=None)
r = np.around(r, decimals=7)

with open('grid.txt', 'w+') as f:
    for item in r:
        f.write("%s\n" % item)

## extended grid to 30 Rg, logarithmically spaced
r_new = np.logspace(np.log(0.0300785), np.log(30.0), num=1324, base=e , dtype=None)
r_new = np.around(r_new, decimals=7)

with open('grid_r_30.dat', 'w+') as f2:
    for item2 in r_new:
        f2.write("%s\n" % item2)

## Formula to derive the new number of cells
dtheta = (np.pi)/600 + 1.
nr = np.log(10./0.03)/np.log(dtheta)
