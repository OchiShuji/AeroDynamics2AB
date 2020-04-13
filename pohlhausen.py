import numpy as np 
from matplotlib import pyplot as plt

'''Blasius解の導出'''
dlt = 0.1
eta = np.arange(0,10.0,dlt)
N = eta.shape[0]
X = np.zeros((N,3))
f2_init = []
a = 0.3
b = 0.4
uinft = 1.0

def blasius(f):
    df = f[1]
    d2f = f[2]
    d3f = -0.5*f[0]*f[2]
    deriv = np.array([df,d2f,d3f])
    return deriv

j = 0
while abs(X[N-1,1]-1) > 0.00000000001:
    f2_init.append((a+b)/2)
    X[0,:] = [0.0,0.0,f2_init[j]]
    for i in range(1,N):
        k1 = blasius(X[i-1,:])
        k2 = blasius(X[i-1,:]+0.5*dlt*k1)
        k3 = blasius(X[i-1,:]+0.5*dlt*k2)
        k4 = blasius(X[i-1,:]+dlt*k3)
        X[i,:] = X[i-1,:]+dlt/6 * (k1+2*k2+2*k3+k4)
    if X[N-1,1] > 1.0:
        b = (a + b)*0.5
    else:
        a = (a + b)*0.5
    j = j + 1
    
print(f2_init)
print(X[N-1,1])

u = np.array(uinft*X[:,1])
bl_thickness = 0.995*np.ones(N)

delta = []
for i in range(0,len(eta)):
    delta.append(abs(bl_thickness[i]-u[i]))
delta_min = min(delta)
eta_bl = eta[delta.index(delta_min)]
print(eta_bl)

'''Pohlhausen解の導出'''
pr = 1
dtheta_init = []
Y = np.zeros((N,2))
c = 0
d = 1

def pohlhausen(i,theta):
    dtheta = theta[1]
    d2theta = -0.5*pr*X[i,0]*theta[1]
    deriv2 = np.array([dtheta,d2theta])
    return deriv2

j = 0
while abs(Y[N-1,1]-1) > 0.00000000001:
    dtheta_init.append((c+d)/2)
    Y[0,:] = [0.0,j]
    for i in range(1,N):
        k1 = pohlhausen(i,Y[i-1,:])
        k2 = pohlhausen(i,Y[i-1,:]+0.5*dlt*k1)
        k3 = pohlhausen(i,Y[i-1,:]+0.5*dlt*k2)
        k4 = pohlhausen(i,Y[i-1,:]+dlt*k3)
        Y[i,:] = Y[i-1,:]+dlt/6 * (k1+2*k2+2*k3+k4)
    if Y[N-1,1] > 1.0:
        d = (c + d)*0.5
    else:
        c = (c + d)*0.5
    j = j + 1
    
print(dtheta_init)
print(Y[N-1,1])


plt.plot(eta,u)
plt.plot(eta,bl_thickness,color="black",linewidth=0.5,linestyle="--")
plt.plot(eta_bl*np.ones(N),np.arange(0,1.0,0.01),color="black",linewidth=0.5,linestyle="--")
ax = plt.gca() 
ax.spines["right"].set_color("none")
ax.spines["top"].set_color("none")
ax.xaxis.set_ticks_position('bottom')
ax.spines["bottom"].set_position(("data",0))
ax.yaxis.set_ticks_position("left")
ax.spines["left"].set_position(("data",0))
plt.xticks([0,2,4,5.3,6,8,10],["0","2","4","5.3","6","8","10"])
plt.xlabel(r"$\eta=\frac{y}{x}\sqrt{Re_x}$")
plt.ylabel(r"$\frac{u}{U_{\infty}}=\frac{df}{d\eta}$")
plt.yticks([0.0,0.2,0.4,0.6,0.8,0.995],["0.0","0.2","0.4","0.6","0.8","0.995"])
plt.show()