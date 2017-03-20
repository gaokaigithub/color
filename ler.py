import numpy as np
def ler(p):
    yd = np.loadtxt("xyz.csv", delimiter=",", usecols=(2,), unpack=True)
    s1 = 0
    s2 = 0
    for i in range(len(yd)):
        s1 = s1+p[i]
        s2 = s2+yd[i]*p[i]
    kl = 683*(s2/s1)
    return kl



