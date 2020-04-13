import numpy as np 
from matplotlib import pyplot as plt 

'''円筒+循環'''
vinft = 1
a = 1
gamma = 2

def psi(x,y):
    return vinft*(y-a**2*y/(x**2+y**2))+gamma*np.log(np.sqrt(x**2+y**2)/a)/2/np.pi

n = 256
x = np.linspace(-3, 3, n)
y = np.linspace(-3, 3, n)
X,Y = np.meshgrid(x, y)
print(Y)

plt.axes([0.025, 0.025, 0.95, 0.95])

plt.contour(X, Y, psi(X, Y), 500, cmap="winter", linewidths=0.5)
plt.contour(X, Y, psi(X, Y), 1,linewidths=2)

#plt.clabel(C, inline=1, fontsize=10)

plt.xlim(X.min() * 1.3, X.max() * 1.3)
plt.xticks(())
plt.yticks(())
plt.show()


