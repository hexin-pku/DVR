#!/usr/bin/env python
# coding=utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

a=pd.read_csv('eigv.dat',header=None,sep= '\s+').values
for i in range(5):
    lenrange = len(a[i])
    x = np.linspace(-10,10,lenrange)
    ps = a[i]
    plt.plot(x,ps,'-',label='state %d-th'%i)
plt.legend(loc=1)
plt.show()

