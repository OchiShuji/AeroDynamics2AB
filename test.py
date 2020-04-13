import numpy as np 
from matplotlib import pyplot as plt 


n = 256
c = 1
a = 0.6
vinft = 1
gamma = 5
r = np.linspace(0.01,1.0,n)
theta = np.linspace(0.0,2*np.pi,n)
x,y,z,zeta,xi,eta,f,psi,phi = [],[],[],[],[],[],[],[],[]
for i in range(0,len(r)):
    for j in range(0,len(theta)):
        x.append(r[i]*np.cos(theta[j]))
        y.append(r[i]*np.sin(theta[j]))
for i in range(0,len(x)):
    z.append(complex(x[i],y[i]))
for i in range(0,len(z)):
    zeta.append(z[i]+(c**2)/z[i])
    xi.append(zeta[i].real)
    eta.append(zeta[i].imag)
for i in range(0,len(zeta)):
    f.append(complex(vinft*(zeta[i]+a**2/zeta[i]),gamma*np.log(zeta[i])/2/np.pi))
    psi.append(f[i].imag)
    phi.append(f[i].real)
X = np.array(x).reshape(256,256)
Y = np.array(y).reshape(256,256)
Z = np.array(z).reshape(256,256)
Psi = np.array(psi).reshape(256,256)
Phi = np.array(phi).reshape(256,256)
Xi = np.array(xi).reshape(256,256)
Eta = np.array(eta).reshape(256,256)

#plt.contour(Xi,Eta,Psi,50, cmap="winter", linewidths=0.5)
#plt.contour(Xi,Eta,Psi,1, cmap="winter", linewidths=2)
plt.show()
plt.plot(xi,eta)
plt.show()
plt.plot(Phi)
plt.show()