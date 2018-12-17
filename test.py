import numpy as np

a=np.linspace(-10,10,2001)
e=np.exp(-a**2)*np.cos(8*a)
s=np.exp(-a**2)*np.sin(10*a)
print(np.mean(e),np.mean(s))
