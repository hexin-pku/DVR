import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
from matplotlib import cm
import pandas as pd


x,p = np.mgrid[-10:10:201j,-10:10:201j] # generate with the boundary
x = x[1:-1,1:-1]
p = p[1:-1,1:-1]

rho = pd.read_csv('rho.dat',header=None, sep='\s+').values

poly = np.polynomial.laguerre.Laguerre((1))
rhoexact = np.zeros((199,199))

for i in range(len(x)):
    for j in range(len(p)):
        xi = x[i,0]**2 + p[0,j]**2
        # (-1)^n (1/\[Pi]) E^-(p^2 + q^2) LaguerreL[n, 2 (p^2 + q^2)]
        rhoexact[i,j] = 1/np.pi * poly(2*xi) * np.exp(-xi)# L_n(x)
        
#print(len(rhoexact[0]))

ax = plt.subplot(111, projection='3d')
ax.set_title('denisty');
#ax.plot_surface(x,p,rho,rstride=2, cstride=1, cmap=cm.jet)
ax.plot_surface(x,p,rho-rhoexact,rstride=2, cstride=1, cmap=cm.jet)
#ax.plot_surface(x,p,-rhoexact,rstride=2, cstride=1, cmap=cm.jet)

ax.set_xlabel('x')
ax.set_ylabel('p')
ax.set_zlabel('D')
plt.show()
