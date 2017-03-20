import numpy as np
from xyz import xyz
from ra import ra
from ler import ler
import math
import time

# 求解一元二次方程
def xs(a,b,c):
    d = b**2-4*a*c
    if d>=0:
        k1 = (-b - math.sqrt(d))/(2*a)
        k2 = (-b + math.sqrt(d))/(2*a)
        return k1,k2
    else:
        return None

# 根据功率分布混合光
def twomix(p1,p2,A):
    xd, yd, zd = np.loadtxt("xyz.csv", delimiter=",", usecols=(1, 2, 3), unpack=True)
    X1 = 0
    Y1 = 0
    Z1 = 0
    p = []
    for j in range(len(p1)):
        ap = p1[j]+A*p2[j]
        p.append(ap)
    for i in range(len(p)):
        d = p[i]*xd[i]*5
        e = p[i]*yd[i]*5
        f = p[i]*zd[i]*5
        X1 = X1+d
        Y1 = Y1+e
        Z1 = Z1+f
    k = 100/Y1
    X = k*X1
    Y = 100
    Z = k*Z1
    x = X/(X+Y+Z)
    y = Y/(X+Y+Z)
    return x,y

# 求三基色混合目标色坐标的功率分布，l是三基色参数，g目标色坐标，r精度
def twoopt(l,g,r):
    xb = xyz(l[0:2],[1])[0]
    yb = xyz(l[0:2],[1])[1]
    xg = xyz(l[2:4],[1])[0]
    yg = xyz(l[2:4],[1])[1]
    xr = xyz(l[4:6],[1])[0]
    yr = xyz(l[4:6],[1])[1]
    c1 = g[0]
    c2 = g[1]
    a = (c1-xb)**2 - r**2
    b = 2*(yb-c2)*(c1-xb)
    c = (yb-c2)**2 - r**2
    if xs(a,b,c):
        k1,k2 = xs(a,b,c)
    else:
        return None
    k3 = (yg-yr)/(xg-xr)
    x1 = (yr-yb+k1*xb-k3*xr)/(k1-k3)
    x2 = (yr-yb+k2*xb-k3*xr)/(k2-k3)
    if x1<x2:
        xp = x1
        xq = x2
    else:
        xp = x2
        xq = x1
    # 以上，求出xp,xq，下面求ok1即r的比例参数
    l2s = []
    l1 = 0
    ok1 = 0
    for i in range(10):
        x0 = xyz(l[2:6],[1,0.1*2**i])[0]
        if xp < x0 < xq:
            ok1 = 0.1*2**i
        else:
            if x0 < xp:
                l1 = 0.1 * 2 ** i
            if x0 > xq:
                l2s.append(0.1 * 2 ** i)
    if ok1 == 0:
        l2 = l2s[0]
        while True:
            j = (l1+l2)/2
            x0 = xyz(l[2:6],[1,j])[0]
            if xp < x0 < xq:
                ok1 = j
                break
            else:
                if x0<xp:
                    l1 = j
                if x0>xq:
                    l2 = j
    # 以上，求出g,r的比例参数，下面求b的比例参数
    xp0 = xyz(l[2:6],[1,ok1])[0]
    yp0 = xyz(l[2:6],[1,ok1])[1]
    p1 = xyz(l[2:6],[1,ok1])[-1]
    p2 = xyz(l[0:2],[1])[-1]
    kb = (yb-yp0)/(xb-xp0)
    a1 = 1+kb**2
    b1 = 2*kb*(yb-c2-kb*xb)-2*c1
    d1 = (yb-c2-kb*xb)**2-r**2+c1**2
    if xs(a1,b1,d1):
        x3,x4 = xs(a1,b1,d1)
    else:
        return None
    if x3>x4:
        xbr = x4
        xbs = x3
    else:
        xbr = x3
        xbs = x4
    # 以上求出xbr和xbs，下面求b的比例参数
    l4s = []
    l3 = 0
    ok2 = 0
    for i in range(10):
        x5 = twomix(p1,p2,0.1*2**i)[0]
        if xbr < x5 < xbs:
            ok2 = 0.1 * 2 ** i
        else:
            if x5 < xbr:
                l3 = 0.1 * 2 ** i
            if x5 > xbs:
                l4s.append(0.1 * 2 ** i)
    if ok2 == 0:
        l4 = l4s[0]
        while True:
            j = (l3 + l4) / 2
            x5 = twomix(p1,p2,j)[0]
            if xbr < x5 < xbs:
                ok2 = j
                break
            else:
                if x5 < xbr:
                    l3 = j
                if x5 > xbs:
                    l4 = j
    # 以上求出b的比例参数
    return [ok2,1,ok1]
if __name__ == '__main__':
    start = time.clock()
    al = twoopt([460,25,540,36,615,13],[0.28,0.35],0.00001)
    end = time.clock()
    r = ra([460,25,540,36,615,13],al)
    p = xyz([460,25,540,36,615,13],al)[-1]
    kl = ler(p)
    print(r)
    print(al)
    print(kl)
    print(end-start)










