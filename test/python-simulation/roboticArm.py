from random import uniform, gauss
from math import sqrt, pi, sin, cos

def truncGaussian( mean, sigma, a, b):
    r = gauss(mean, sigma)
    while (r < a or r > b):
        r = gauss(mean, sigma)
    return r

def simulateOnce(nSteps):
    angles = [10, 60, 110, 160, 140, 100, 60, 20, 10, 0]
    x = truncGaussian(0.0, 0.025, -0.1, 0.1)
    y = truncGaussian(0.0, 0.01, -0.5, 0.5)
    for j in range(nSteps):
        for ang in angles:
            d  = 1 + uniform(-0.02, 0.02)
            t = (1 + truncGaussian(0.0, 0.01, -0.05, 0.05))*ang * pi/180
            x = x + d * cos(t)
            y = y + d * sin(t)
    return x

def simulateAll(nSamples):
    ex = 0.0
    n1 = 0
    for j in range(nSamples):
        x = simulateOnce(100)
        ex = ex + x
        if x >= 272.0:
            n1 += 1
    print('E(x) = %f ' %(ex/nSamples))
    print('P(x >= 272.0) = %d/%d' %(n1, nSamples))

if __name__ == '__main__':
    simulateAll(100000)

