import numpy as np 
from matplotlib import pyplot as plt 

r = np.linspace(0.25,1.5,25)
theta = np.linspace(0,2*np.pi,25)
x,y,u_x,u_y = [],[],[],[]
for i in range(0,len(r)):
    for j in range(0,len(r)):
        x.append(r[i]*np.cos(theta[j]))
        y.append(r[i]*np.sin(theta[j]))
        if r[i] <= 1:
            u_x.append(1/r[i]*(-np.sin(theta[j])))
            u_y.append(1/r[i]*(np.cos(theta[j])))
        else:
            u_x.append(r[i]*(-np.sin(theta[j])))
            u_y.append(r[i]*(np.cos(theta[j])))

plt.quiver(x,y,u_x,u_y,linewidths=.25)
plt.show()
