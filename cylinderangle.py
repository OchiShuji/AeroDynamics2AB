import numpy as np 
from matplotlib import pyplot as plt 


'''円筒+循環+迎角'''
alpha = 5 * np.pi /180
beta = 0 * np.pi /180
c = 0.5
vinft = 1
a = 1
gamma = 2

n = 256
x = np.linspace(-3, 3, n)
y = np.linspace(-3, 3, n)
X,Y = np.meshgrid(x, y)
z,Z,f,psi = [],[],[],[]
print(X)
print(Y)

for j in range(0,len(x)):
    for i in range(0,len(y)):
        Z.append(complex(x[i],y[j]))

for i in range(0,len(Z)):
    z.append((Z[i] - a*complex(np.cos(np.pi-beta),np.sin(np.pi-beta)) - c)/complex(np.cos(alpha),np.sin(alpha)))

for i in range(0,len(z)):
    f.append(complex(vinft*(z[i]+a**2/z[i]),gamma*np.log(z[i])/2/np.pi))
    psi.append(f[i].imag)
psi = np.array(psi)
psipsi = psi.reshape(256,256)

#print(psipsi)
#print(len(psi))

plt.axes([0.025, 0.025, 0.95, 0.95])

plt.contour(X, Y, psipsi, 750, cmap="winter", linewidths=0.5)
plt.contour(X, Y, psipsi, 1,linewidths=2)

plt.xlim(X.min() * 1.3, X.max() * 1.3)
plt.xticks(())
plt.yticks(())
ax = plt.gca()  
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
plt.show()


