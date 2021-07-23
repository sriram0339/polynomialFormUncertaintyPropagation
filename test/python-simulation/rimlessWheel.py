from random import gauss, uniform
from math import pi, sin, cos
import numpy as np
from scipy.io import savemat

def randn():
    return gauss(0,1)

def truncGaussian( mean, sigma, a, b):
    r = gauss(mean, sigma)
    while (r < a or r > b):
        r = gauss(mean, sigma)
    return r

nSims = 100000;
nSteps = 5000;

def runSimulations(nSims, nSteps):
    ex = 0
    n1 = 0
    n2 = 0
    n3 = 0
    n4 = 0
    meanParam = 4.0
    simVals = np.zeros((nSims,1))
    for i in range(nSims):
        x = 0.2 * uniform(0,1) - 0.1
        theta = pi/6
        g = 10.0
        for j in range(nSteps):
            w = pi/180 * truncGaussian(meanParam, 1.5, meanParam-9, meanParam+9)
            beta1 = theta/2 + w
            beta2 = theta/2 - w
            x = cos(theta)**2 * ( x + 2*g*(1- cos(beta1))) - 2*g * (1-cos(beta2))
        ex = ex + x
        if (x <= 0):
            n1 += 1
        simVals[i,0] = x
    print('E(x) = %f' % (ex/nSims))
    print('P(x <= 0) = %d/%d' %(n1, nSims))
    savemat('rimlesswheel-model.mat', {'x':simVals})

if __name__ == '__main__':
    runSimulations(10000,2000)

