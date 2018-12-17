import numpy as np
import pandas as pd

def VHO(x):
    return 0.5*x*x

# the DVR solvor
xmin=-10;xmax=10;Nx=200; dx=(xmax-xmin)/Nx
pmin=-10;pmax=10;Np=200; dp=(pmax-pmin)/Np

coeff = 0.5/(xmax-xmin)**2 * 0.5 * np.pi**2

H = np.zeros((Nx-1,Nx-1))
R = np.zeros((Nx-1,Np-1))

for i in range(Nx-1):
    print(np.pi)
    print((2*Nx**2+1)/3,coeff)
    H[i,i] = coeff * ( (2*Nx**2+1)/3 - 1/np.sin(np.pi*(i+1)/Nx)**2 )
    print(coeff * ( (2*Nx**2+1)/3), coeff*( 1/np.sin(np.pi*(i+1)/Nx)**2 ))
    exit()
    for j in range(i+1,Nx-1):
        H[i,j] = coeff * (-1)**(i-j) * ( 1/np.sin(np.pi*(i-j)/(2*Nx))**2 
                     - 1/np.sin(np.pi*(i+j+2)/(2*Nx))**2  )
        H[j,i] = H[i,j]
    print(H[i,i])
    exit()
    H[i,i] += VHO(xmin+(i+1)*dx)
outH = pd.DataFrame(H)
outH.to_csv('H_py.dat',header=None, index=None, sep=' ')

Es, Vs = np.linalg.eigh(H)
Vs = Vs/np.sqrt(dx) # normalization
print('check the eigenvalues: ', Es[0:6])
print('check the vector normalization : ',np.sum(Vs[:,0]**2)*dx)
outeige = pd.DataFrame(Es)
outeigv = pd.DataFrame(Vs)
outeige.to_csv('eige_py.dat', header=None, index=None, sep=' ')
outeigv.to_csv('eigv_py.dat', header=None, index=None, sep=' ')  

# the wigner calculator
cutoff = 2
beta = 8
for i in range(Nx-1):
    mrange = min(i,Nx-2-i)
    for j in range(Np-1):
        p = pmin + (j+1) * dp
        for k in range(cutoff):
            tmp = 0
            for m in range(-mrange, mrange+1):
                tmp += Vs[i-m,k]*Vs[i+m,k]*np.cos(2*dx*m*p) * 2*dx
            R[i,j] += tmp * np.exp(-beta*Es[k])
Z = 0
for k in range(cutoff):
    Z += np.exp(-beta*Es[k])
R = R/(2*np.pi*Z)

print('check distribution normalization', np.sum(R)*dx*dp)  

out = pd.DataFrame(R)
out.to_csv('rho_py.dat', header=None, index=None, sep=' ')  

            
