# -*- coding: utf-8 -*-
"""
Created on Sat Aug  4 01:07:27 2018

(5)用のやつです
あまり有効なやり方が思いつかなかったのでごり押しでやりました
少し重いですがその時は一つ二つ消せばいいと思います

@author:雨宮潤
"""

import numpy as np
import matplotlib.pyplot as plt

a=0.6
v0=1.0
rad=np.pi/180.0
alfa=10*rad
beta=18*rad
c025=0.5

n=181
m=181
dr=0.015
dtheta=2.0*rad
x=np.zeros((n,n))
y=np.zeros((n,n))
X=np.zeros((n,n))
Y=np.zeros((n,n))
Xze=np.zeros((n,n))
Yze=np.zeros((n,n))
stream=np.zeros((n,n))
poten=np.zeros((n,n))
#--------------------------------

gama=4.0*np.pi*v0*a*np.sin(alfa+beta) #クッタ条件
Zc=a*np.exp(1j*(np.pi-beta))+c025



def potential(s,t):
    z=s+1j*t
    f=v0*(z+a**2/z)+(1j*gama/(2*np.pi))*np.log(z)
    return f

def circle_fc(s,t):
    Z=s+1j*t
    z=Z*np.exp(1j*alfa)+Zc
    return z


def Zeta(s,t):
    z=s+1j*t
    Z=z*np.exp(1j*alfa)+Zc
    Zeta=Z+c025**2/Z
    return Zeta

def pres(s,t):
    z=s+1j*t
    Z=z*np.exp(1j*alfa)+Zc
    f=np.exp(-1j*alfa)*v0*(1-a**2/z**2+1j*gama/(v0*2*np.pi*z))/(1-c025**2/Z**2)
    Cp=1-np.abs(f)**2/v0**2
    return Cp


for i in range(n):
    for k in range(m):
        #z平面上で同心円状に計算点を配置
        theta0=0.0
        theta0=-beta
        rrr=a+dr*i
        theta=theta0+dtheta*k
        x[i,k]=rrr*np.cos(theta)
        y[i,k]=rrr*np.sin(theta)
        poten[i,k]=np.real(potential(x[i,k],y[i,k]))
        stream[i,k]=np.imag(potential(x[i,k],y[i,k]))
        Z=np.exp(1j*alfa)*(x[i,k]+1j*y[i,k])+Zc
        Xze[i,k]=np.real(Z+c025**2.0/Z)
        Yze[i,k]=np.imag(Z+c025**2.0/Z)


Cp=pres(x,y)


Cmc5=np.zeros(n)
x0=np.linspace(-1.00,0,181)

for l in range(n-1):
    dxw=Xze[0,l+1]-Xze[0,l]
    dyw=Yze[0,l+1]-Yze[0,l]
    dnx=dyw
    dny=-dxw
    cpm=(Cp[0,l+1]+Cp[0,l])/2.0
    fx=-cpm*dnx
    fy=-cpm*dny
    for h in range(n):
        Cmc5[h]=Cmc5[h]+((x0[h]-Xze[0,l])*fy+Yze[0,h]*fx)/(0.5*1*1.0**2*16*c025**2)



#----------------------------------
#----------------------------------


alfa=0*rad
gama=4.0*np.pi*v0*a*np.sin(alfa+beta)

for i in range(n):
    for k in range(m):
        theta0=0.0
        theta0=-beta
        rrr=a+dr*i
        theta=theta0+dtheta*k
        x[i,k]=rrr*np.cos(theta)
        y[i,k]=rrr*np.sin(theta)
        poten[i,k]=np.real(potential(x[i,k],y[i,k]))
        stream[i,k]=np.imag(potential(x[i,k],y[i,k]))
        Z=np.exp(1j*alfa)*(x[i,k]+1j*y[i,k])+Zc
        Xze[i,k]=np.real(Z+c025**2.0/Z)
        Yze[i,k]=np.imag(Z+c025**2.0/Z)


Cp=pres(x,y)


Cmc0=np.zeros(n)
x0=np.linspace(-1.00,0,181)

for l in range(n-1):
    dxw=Xze[0,l+1]-Xze[0,l]
    dyw=Yze[0,l+1]-Yze[0,l]
    dnx=dyw
    dny=-dxw
    cpm=(Cp[0,l+1]+Cp[0,l])/2.0
    fx=-cpm*dnx
    fy=-cpm*dny
    for h in range(n):
        Cmc0[h]=Cmc0[h]+((x0[h]-Xze[0,l])*fy+Yze[0,h]*fx)/(0.5*1*1.0**2*16*c025**2)

#-----------------------------------------------------------
#----------------------------------------------------------


alfa=3*rad
gama=4.0*np.pi*v0*a*np.sin(alfa+beta)

for i in range(n):
    for k in range(m):
        theta0=0.0
        theta0=-beta
        rrr=a+dr*i
        theta=theta0+dtheta*k
        x[i,k]=rrr*np.cos(theta)
        y[i,k]=rrr*np.sin(theta)
        poten[i,k]=np.real(potential(x[i,k],y[i,k]))
        stream[i,k]=np.imag(potential(x[i,k],y[i,k]))
        Z=np.exp(1j*alfa)*(x[i,k]+1j*y[i,k])+Zc
        Xze[i,k]=np.real(Z+c025**2.0/Z)
        Yze[i,k]=np.imag(Z+c025**2.0/Z)


Cp=pres(x,y)


Cmc3=np.zeros(n)
x0=np.linspace(-1.00,0,181)

for l in range(n-1):
    dxw=Xze[0,l+1]-Xze[0,l]
    dyw=Yze[0,l+1]-Yze[0,l]
    dnx=dyw
    dny=-dxw
    cpm=(Cp[0,l+1]+Cp[0,l])/2.0
    fx=-cpm*dnx
    fy=-cpm*dny
    for h in range(n):
        Cmc3[h]=Cmc3[h]+((x0[h]-Xze[0,l])*fy+Yze[0,h]*fx)/(0.5*1*1.0**2*16*c025**2)

#------------------------------------------------
#------------------------------------------------
        

alfa=-3*rad
gama=4.0*np.pi*v0*a*np.sin(alfa+beta)

for i in range(n):
    for k in range(m):
        theta0=0.0
        theta0=-beta
        rrr=a+dr*i
        theta=theta0+dtheta*k
        x[i,k]=rrr*np.cos(theta)
        y[i,k]=rrr*np.sin(theta)
        poten[i,k]=np.real(potential(x[i,k],y[i,k]))
        stream[i,k]=np.imag(potential(x[i,k],y[i,k]))
        Z=np.exp(1j*alfa)*(x[i,k]+1j*y[i,k])+Zc
        Xze[i,k]=np.real(Z+c025**2.0/Z)
        Yze[i,k]=np.imag(Z+c025**2.0/Z)


Cp=pres(x,y)


Cmc_3=np.zeros(n)
x0=np.linspace(-1.00,0,181)

for l in range(n-1):
    dxw=Xze[0,l+1]-Xze[0,l]
    dyw=Yze[0,l+1]-Yze[0,l]
    dnx=dyw
    dny=-dxw
    cpm=(Cp[0,l+1]+Cp[0,l])/2.0
    fx=-cpm*dnx
    fy=-cpm*dny
    for h in range(n):
        Cmc_3[h]=Cmc_3[h]+((x0[h]-Xze[0,l])*fy+Yze[0,h]*fx)/(0.5*1*1.0**2*16*c025**2)

#----------------------------------------------------
#----------------------------------------------------
        

alfa=-5*rad
gama=4.0*np.pi*v0*a*np.sin(alfa+beta)

for i in range(n):
    for k in range(m):
        theta0=0.0
        theta0=-beta
        rrr=a+dr*i
        theta=theta0+dtheta*k
        x[i,k]=rrr*np.cos(theta)
        y[i,k]=rrr*np.sin(theta)
        poten[i,k]=np.real(potential(x[i,k],y[i,k]))
        stream[i,k]=np.imag(potential(x[i,k],y[i,k]))
        Z=np.exp(1j*alfa)*(x[i,k]+1j*y[i,k])+Zc
        Xze[i,k]=np.real(Z+c025**2.0/Z)
        Yze[i,k]=np.imag(Z+c025**2.0/Z)


Cp=pres(x,y)


Cmc_5=np.zeros(n)
x0=np.linspace(-1.00,0,181)

for l in range(n-1):
    dxw=Xze[0,l+1]-Xze[0,l]
    dyw=Yze[0,l+1]-Yze[0,l]
    dnx=dyw
    dny=-dxw
    cpm=(Cp[0,l+1]+Cp[0,l])/2.0
    fx=-cpm*dnx
    fy=-cpm*dny
    for h in range(n):
        Cmc_5[h]=Cmc_5[h]+((x0[h]-Xze[0,l])*fy+Yze[0,h]*fx)/(0.5*1*1.0**2*16*c025**2)



    
plt.plot(x0,Cmc5,color="black",linewidth=0.75)
plt.plot(x0,Cmc3,color="black",linewidth=0.75)
plt.plot(x0,Cmc0,color="black",linewidth=0.75)
plt.plot(x0,Cmc_3,color="black",linewidth=0.75)
plt.plot(x0,Cmc_5,color="black",linewidth=0.75)
plt.xlim([-0.55,-0.4])
plt.ylim([-1.6,-0.5])
plt.xlabel("x")
plt.ylabel("Cm")
plt.show()