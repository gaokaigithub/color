import numpy as np


#高斯函数,A是功率比例
def cur(a,b,c,A):
    d = (a-b)**2/c**2
    s = A*np.exp(-4*np.log(2)*d)
    return s

#得到高斯函数模拟数据
def curs(e,f,A ):
    g = np.arange(380,781,5)
    h = []
    for i in g:
        h.append(cur(i,e,f,A))
    return h


