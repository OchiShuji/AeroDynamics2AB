import numpy as np 
from matplotlib import pyplot as plt 


'''翼型'''
alpha = 5 * np.pi /180
beta = 20 * np.pi /180
c = 0.5
vinft = 1
a = 0.6
gamma = 5

n = 256
x = np.linspace(-3, 3, n)
y = np.linspace(-3, 3, n)
X,Y = np.meshgrid(x, y)
z,Z,f,psi = [],[],[],[]
zeta=[]

for j in range(0,len(x)):
    for i in range(0,len(y)):
        zeta.append(complex(x[i],y[j]))


for i in range(0,len(zeta)):
    delta = zeta[i]*zeta[i]-4*c**2
    red = delta.real
    imd = delta.imag
    absd = np.sqrt(red*red+imd*imd)
    Z.append(zeta[i]*0.5+np.sqrt((red+absd)/2)*0.5+0.5*complex(0.0,imd/np.sqrt((red+absd)/2)))


for i in range(0,len(Z)):
    z.append((Z[i] - a*complex(np.cos(np.pi-beta),np.sin(np.pi-beta)) - c)/complex(np.cos(alpha),np.sin(alpha)))
   

for i in range(0,len(z)):
    f.append(complex(vinft*(z[i]+a**2/z[i]),gamma*np.log(z[i])/2/np.pi))
    psi.append(f[i].imag)
psi = np.array(psi)
psipsi = psi.reshape(256,256)
XX = np.array(X).reshape(256,256)
YY = np.array(Y).reshape(256,256)

#print(psipsi)
#print(len(psi))

plt.axes([0.025, 0.025, 0.95, 0.95])

plt.contour(XX, YY, psipsi, 50, cmap="winter", linewidths=0.5)
plt.contour(XX, YY, psipsi, 1,linewidths=2)

#plt.xlim(X.min() * 1.3, X.max() * 1.3)
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