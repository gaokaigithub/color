import numpy as np
from xyz import xyz
from tem import tem
import math

#xyz实现颜色配色，参数也全部改为列表
#求出uk,vk,t
def uvkt(l,Al):
    x = xyz(l,Al)[0]
    y = xyz(l,Al)[1]
    uk = 4*x/(-2*x+12*y+3)
    vk = 6*y/(-2*x+12*y+3)
    t = tem(x,y)
    return [uk,vk,t]


#fm函数
def fm(a,b,c,t):
    m = 10**4/t
    f = a+b*m+c*m**2
    return f
#求出ur,vr,Uri,Vri,Wri
def rdata(l,Al):
    x = xyz(l,Al)[0]
    y = xyz(l,Al)[1]
    t = tem(x,y)
    if t<= 5000:
        a,b,c = np.loadtxt("fm.csv",delimiter=",",usecols=(0,1,2),unpack=True,skiprows=1)
    elif t<=10000:
        a,b,c = np.loadtxt("fm.csv",delimiter=",",usecols=(3,4,5),unpack=True,skiprows=1)
    else:
        a,b,c = np.loadtxt("fm.csv",delimiter=",",usecols=(6,7,8),unpack=True,skiprows=1)
    ur = fm(a[0],b[0],c[0],t)
    vr = fm(a[1],b[1],c[1],t)
    Uri = []
    Vri = []
    Wri = []
    for ui in range(2,10):
        Uri.append(fm(a[ui],b[ui],c[ui],t))
    for vi in range(16,24):
        Vri.append(fm(a[vi],b[vi],c[vi],t))
    for wi in range(30,38):
        Wri.append(fm(a[wi],b[wi],c[wi],t))

    return [ur,vr,Uri,Vri,Wri]


#求uki,vki函数，此处用到了功率分布函数，我需要xyz函数导出来
def ki(l,Al):
    a1,a2,a3,a4,a5,a6,a7,a8 = np.loadtxt("β.csv",delimiter=",",
                                         usecols=(1,2,3,4,5,6,7,8),unpack=True)
    xd,yd,zd = np.loadtxt("xyz.csv",delimiter=",",usecols=(1,2,3),unpack=True)
    aall = [a1,a2,a3,a4,a5,a6,a7,a8]
    p = xyz(l,Al)[-1]
    k = xyz(l,Al)[-2]
    uki = []
    vki = []
    Yki = []
    for i in range(8):
        g = aall[i]
        Xk = 0
        Yk = 0
        Zk = 0
        for t in range(len(p)):
            d = k*g[t]*p[t]*xd[t]*5
            e = k*g[t]*p[t]*yd[t]*5
            f = k*g[t]*p[t]*zd[t]*5
            Xk = Xk+d
            Yk = Yk+e
            Zk = Zk+f
        uk = 4*Xk/(Xk+15*Yk+3*Zk)
        vk = 6*Yk/(Xk+15*Yk+3*Zk)
        uki.append(uk)
        vki.append(vk)
        Yki.append(Yk)
    return [uki,vki,Yki]

#计算c,d
def cd(u,v):
    c = (4-u-10*v)/v
    d = (1.708*v+0.404-1.481*u)/v
    return c,d
#计算ukip,vkip,很明显这里有问题，int和list
def kip(l,Al):
    ur = rdata(l,Al)[0]
    vr = rdata(l,Al)[1]
    uk = uvkt(l,Al)[0]
    vk = uvkt(l,Al)[1]
    uki = ki(l,Al)[0]
    vki = ki(l,Al)[1]
    cr,dr = cd(ur,vr)
    ck,dk = cd(uk,vk)
    ukips = []
    vkips = []
    for i in range(len(uki)):
        cki,dki = cd(uki[i],vki[i])
        cf = cr*cki/ck
        df = dr*dki/dk
        ukip = (10.872+0.404*cf-4*df)/(16.518+1.481*cf-df)
        vkip = 5.520/(16.518+1.481*cf-df)
        ukips.append(ukip)
        vkips.append(vkip)
    return ukips,vkips
#求Wki,Uki,Vki
def WUV(l,Al):
    Yki = ki(l,Al)[-1]
    ukips,vkips = kip(l,Al)
    ur = rdata(l,Al)[0]
    vr = rdata(l,Al)[1]
    Wki = []
    Uki = []
    Vki = []
    for i in range(len(Yki)):
        wki = 25*(Yki[i])**(1/3)-17
        uki = 13*wki*(ukips[i]-ur)
        vki = 13*wki*(vkips[i]-vr)
        Wki.append(wki)
        Uki.append(uki)
        Vki.append(vki)
    return [Uki,Vki,Wki]

def ra(l,Al):
    rd = rdata(l,Al)
    kd = WUV(l,Al)
    Uri = rd[2]
    Vri = rd[3]
    Wri = rd[4]
    Uki = kd[0]
    Vki = kd[1]
    Wki = kd[2]
    ris = []
    for i in range(len(Uri)):
        u = (Uri[i]-Uki[i])**2
        v = (Vri[i]-Vki[i])**2
        w = (Wri[i]-Wki[i])**2
        e = math.sqrt(u+v+w)
        ri = 100-4.6*e
        ris.append(ri)
    Ra = sum(ris)/8
    x = xyz(l,Al)[0]
    y = xyz(l,Al)[1]
    t = tem(x, y)
    return Ra,t,x,y




