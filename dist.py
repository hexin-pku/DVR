import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
from matplotlib import cm
import pandas as pd

Nx = 200; xmin=-10; xmax=10; dx=(xmax-xmin)/Nx
Np = 200; pmin=-10; pmax=10; dp=(pmax-pmin)/Np

Vs = pd.read_csv('eigv.dat',header=None,sep='\s+').values.T
Es = pd.read_csv('eige.dat',header=None,sep='\s+').values

print(Es)

rho = np.zeros((Nx-1,Np-1))
beta = 8
cutoff = 10

for i in range(1,Nx):
    x = xmin + (xmax-xmin)*i/Nx
    for j in range(1,Np-1):
        p = pmin + (pmax-pmin) * j / Np
        tmp = 0
        for k in range(cutoff):
            for d in range(min(i,Nx-i)):
                tmp += Vs[k,i+d]*Vs[k,i-d]*np.cos(2*d*dx*p)*2*dx
            rho[i,j] += tmp * np.exp(-beta*Es[k,0])

exit(-1)

x,p = np.mgrid[-10:10:199j,-10:10:199j] #100j

ax = plt.subplot(111, projection='3d')
ax.set_title('denisty');
ax.plot_surface(x,p,rho,rstride=2, cstride=1, cmap=cm.jet)

ax.set_xlabel('x')
ax.set_ylabel('p')
ax.set_zlabel('D')
plt.show()
