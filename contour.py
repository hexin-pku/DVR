import numpy as np
import pandas as pd
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt

x,p = np.mgrid[-10:10:201j,-10:10:201j] # generate with the boundary
x = x[1:-1,1:-1]
p = p[1:-1,1:-1]
rho = pd.read_csv('rho.dat',header=None, sep='\s+').values

side = 80
x = x[side:-side,side:-side]
p = p[side:-side,side:-side]
rho = rho[side:-side,side:-side]

fig, ax = plt.subplots()
CS = ax.contour(x, p, rho)
ax.clabel(CS, inline=1, fontsize=10)
ax.set_title('Wigner density function')
ax.set_xlabel('X')
ax.set_ylabel('P')
savename = str(input('give a save name: ')) + '.png'
plt.savefig(savename)
plt.show()
