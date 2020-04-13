import numpy as np 
from matplotlib import pyplot as plt 

n = 256
m = 2
r = np.linspace(0.1,10,n)
theta = np.linspace(0,2*np.pi,n)
x,y,psi= [],[],[]
for i in range(0,len(r)):
    for j in range(0,len(theta)):
        x.append(r[i]*np.cos(theta[j]))
        y.append(r[i]*np.sin(theta[j]))
        psi.append(r[i]**m*np.sin(m*theta[j]))
X = np.array(x).reshape(256,256)
Y = np.array(y).reshape(256,256)
psipsi = np.array(psi).reshape(256,256)


plt.axes([0.025, 0.025, 0.95, 0.95])

plt.contour(X, Y, psipsi, 75, cmap="winter", linewidths=0.5)
plt.contour(X, Y, psipsi, 1,linewidths=2)

#plt.clabel(C, inline=1, fontsize=10)

plt.xlim(X.min() * 1.3, X.max() * 1.3)
plt.xticks(())
plt.yticks(())
plt.show()
print(X,Y)
