import numpy as np
from cur import curs
#求出色坐标以及三刺激值、色温,K值做1处理,使用1931的cie值
def xyz(l,Al):
    xd,yd,zd = np.loadtxt("xyz.csv",delimiter=",",usecols=(1,2,3),unpack=True)
    X1 = 0
    Y1 = 0
    Z1 = 0
    li = int(len(l) / 2)
    p = []
    for i in range(li):
        a, b = l[2*i:2*i + 2]
        A = Al[i]
        pp = curs(a, b,A=A)
        if len(p) == 0:p = pp
        else:
            for j in range(len(pp)):
                p[j] = p[j]+pp[j]

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
    return [x,y,X,Y,Z,k,p]
