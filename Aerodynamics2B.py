# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 15:30:46 2018

課題二のやつです
翼まわりの流線の図を作成するものです
プリントのやつをpython用に翻訳しただけですが適当に使ってください
何かあったら雨宮までお願いします

@author:雨宮潤
"""

import numpy as np
import matplotlib.pyplot as plt

#初期条件です。ここをいじって色々な翼について調べていきます
a=0.6
gama=5.0
v0=1.0
rad=np.pi/180.0
alfa=10*rad
beta=18*rad
c025=0.5
#-----------------------------------

#計算範囲及び精度-----------------
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
cgama=1j*gama/2.0*np.pi

#複素ポテンシャルを求める関数です
def potential(s,t):
    z=s+1j*t
    f=v0*(z+a**2/z)+(1j*gama/(2*np.pi))*np.log(z)
    return f


def circle_fc(s,t):
    Z=s+1j*t
    z=Z*np.exp(1j*alfa)+Zc
    return z

#ζ座標に変換するやつです
def Zeta(s,t):
    z=s+1j*t
    Z=z*np.exp(1j*alfa)+Zc
    Zeta=Z+c025**2/Z
    return Zeta

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

#翼表面部分はi=0が対応してるらしいです
Wing=Zeta(x[0,:],y[0,:])
#-------------------------------

#ｚ座標系からζ座標系の圧力分布に変換する関数です
def pres(s,t):
    z=s+1j*t
    Z=z*np.exp(1j*alfa)+Zc
    f=np.exp(-1j*alfa)*v0*(1-a**2/z**2+1j*gama/(v0*2*np.pi*z))/(1-c025**2/Z**2)
    Cp=1-np.abs(f)**2/v0**2
    return Cp

Cp=pres(x,y)

#翼に働く空気力を計算します
cxp=0
cyp=0
Z=circle_fc(x,y)
X=np.real(Z)
Y=np.imag(Z)

Af=0.0
Bf=0.0

#翼表面の圧力分布を積分します
for l in range(n-1):
    dxw=Xze[0,l+1]-Xze[0,l]
    dyw=Yze[0,l+1]-Yze[0,l]
    dnx=dyw
    dny=-dxw
    cpm=(Cp[0,l+1]+Cp[0,l])/2.0
    fx=-cpm*dnx
    fy=-cpm*dny
    Af=Af+Xze[0,l]*fy-Yze[0,l]*fx
    Bf=Bf+fy
    cxp=cxp-cpm*dnx
    cyp=cyp-cpm*dny
    
#抵抗係数、揚力係数を求めます
cxp=cxp/(4.0*c025)
cyp=cyp/(4.0*c025)
cdp=cxp*np.cos(alfa)+cyp*np.sin(alfa)
clp=cyp*np.cos(alfa)-cxp*np.sin(alfa)

print("Cd=",end="")
print(cdp)
print("Cl=",end="")
print(clp)

#風圧中心を求めます
xcp=Af/Bf
print("xcp=",end="")
print(xcp)
print(1+xcp)
print((1+xcp)/2)

sub_x = np.linspace(-3,3,256)
sub_y = np.linspace(-3.3,0.75,256)
sub_X,sub_Y = np.meshgrid(sub_x,sub_y)
def bar(x,y):
        return  y
#以下にあるコードは(1)～(3)までのやつなので、適宜　"　を消して使ってください

plt.axes([0.04, 0.04, 0.75, 0.9])
C = plt.contour(Xze,Yze,Cp,30,cmap="jet",linewidths=0.5)
plt.xticks([],[])
plt.yticks([],[])
plt.plot(np.real(Wing),np.imag(Wing),color='k')
plt.axes([0.9, 0.04, 0.05, 0.9])
plt.xticks([],[])
plt.yticks([-3.3,-3.0,-2.7,-2.4,-2.1,-1.8,-1.5,-1.2,-0.9,-0.6,-0.3,0.0,0.3,0.6],["-3.3","-3.0","-2.7","-2.4","-2.1","-1.8","-1.5","-1.2","-0.9","-0.6","-0.3","0.0","0.3","0.6"],fontsize=7)
plt.contour(sub_X,sub_Y,bar(sub_X,sub_Y),30,cmap="jet",linewidths=0.5)
plt.ylabel("Cp")
plt.show()

plt.axes([0.2, 0.025, 0.6, 0.95])
plt.plot(np.real(Wing),np.imag(Wing),color='k')
plt.plot(Xze[0],Cp[0],color='b')
plt.yticks([-3,-2,-1,0,1,2,3],["-3","-2","-1","0","1","2","3"])
plt.xlim(min(Xze[0])*2,max(Xze[0])*2)
plt.xticks([],[])
plt.ylabel("Cp")
plt.show()


plt.contour(Xze,Yze,stream,75,cmap="winter",linewidths=0.5)
plt.plot(np.real(Wing),np.imag(Wing),color='k')
ax = plt.gca()  
ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
ax.yaxis.set_major_locator(plt.MultipleLocator(0.5))
ax.grid(which='major', axis='x', linewidth=0.75, linestyle='--',color="black")
ax.grid(which='major', axis='y', linewidth=0.75, linestyle='--',color="black")
plt.show()