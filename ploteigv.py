#!/usr/bin/env python
# coding=utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

a=pd.read_csv('eigv.dat',header=None,sep= '\s+').values
a=a.T
for i in range(len(a)):
    lenrange = len(a[i])
    x = np.linspace(-10,10,lenrange)
    ps = a[i]
    plt.plot(x,ps,'-',label='%d'%i)

plt.show()

