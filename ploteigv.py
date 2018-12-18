#!/usr/bin/env python
# coding=utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

cut = int(input('give a cut number (from 0): '))

a=pd.read_csv('eigv.dat',header=None,sep= '\s+').values
e=pd.read_csv('eige.dat',header=None,sep= '\s+').values
x = np.linspace(-10,10,len(e)+2)
x = x[1:-1]

for i in range(cut):
    ps = a[i] + e[i]
    plt.plot(x,ps,'-',label='state %d-th'%i)

plt.plot(x,0.5*x**2,'k--',label='PES')
plt.ylim((-0.5,e[cut]+1))
    
plt.legend(loc=1)

savename = str(input('give a save name: ')) + '.png'
plt.savefig(savename)
plt.show()

