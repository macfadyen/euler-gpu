import sys
import math
import numpy as np
import matplotlib.pyplot as plt

x, rho, px, py, e = np.loadtxt(sys.argv[1]).T
ni = len(x)

plt.plot(x,e)
plt.show()
