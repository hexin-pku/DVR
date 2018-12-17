import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
from matplotlib import cm
import pandas as pd


x,p = np.mgrid[-10:10:199j,-10:10:199j] #100j为设置曲面平滑度

a = pd.read_csv('dist.dat', header=None, sep='\s+')
a = a.values


ax = plt.subplot(111, projection='3d')
ax.set_title('denisty');
ax.plot_surface(x,p,a,rstride=2, cstride=1, cmap=cm.jet)

ax.set_xlabel('x')
ax.set_ylabel('p')
ax.set_zlabel('D')
plt.show()
