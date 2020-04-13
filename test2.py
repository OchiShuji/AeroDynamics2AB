import numpy as np 
from matplotlib import pyplot as plt 

n = 256
r = np.linspace(0.01,3,n)
theta = np.linspace(0,2*np.pi,n)
alpha = 5 * np.pi /180
beta = 20 * np.pi /180
c = 0.5
vinft = 1
a = 0.6
gamma = 5
x = []
y = []
zeta = []
z = []

for i in range(0,n):
    for j in range(0,n):
        zeta.append(complex(r[i]*np.cos(theta[j]),r[i]*np.sin(theta[j])))
    
for i in range(0,len(zeta)):
    z.append(zeta[i]+a**2/zeta[i])
    x.append(z[i].real)
    y.append(z[i].imag)

X = np.array(x)
Y = np.array(y)
XX = X.reshape(256,256)
YY = Y.reshape(256,256)

plt.plot(XX,YY)
plt.show()