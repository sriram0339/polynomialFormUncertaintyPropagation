from random import uniform
import numpy as np
from scipy.io import savemat


def simulateOnce(nSteps):
    x1 = uniform(19.0, 21.1)
    x2 = uniform(19.0, 21.1)
    x3 = uniform(19.0, 21.1)
    x4 = uniform(19.0, 21.1)
    x5 = uniform(19.0, 21.1)
    x6 = uniform(19.0, 21.1)
    for j in range(nSteps):
        x1n =  x1 + (0.5 - 9.963e-6*x1*x5 - 1.925e-5*x1)*36
        x2n =  x2 + (1.5e-3 +  1.5e-2*(x1**2/(547600 + x1**2)) - 8e-4*x2)*36 
        x3n =  x3 + (8e-4*x2 - 1.444e-4*x3)*36
        x4n =  x4 + (1.66e-2*x3 - 9e-4*x4)*36
        x5n =  x5 + (9e-4*x4 - 1.66e-7*x4**2 - 9.963e-6*x5*x6)*36
        x6n =  x6 + (0.5 - 3.209e-5*x6 - 9.963e-6*x5*x6)*36
        (x1, x2, x3, x4, x5, x6) = (x1n , x2n , x3n, x4n, x5n, x6n)
    return (x1, x6)


def simulateAll(nSamples):
    e1 = 0.0
    e6 = 0.0
    nq1 = 0
    nq2 = 0

    for j in range(nSamples):
        (x1, x6) = simulateOnce(16)
        e1 = e1 + x1
        e6 = e6 + x6
        if x1 >= 285:
            nq1 += 1
        if x6 >= 285:
            nq2 += 1
    print('E(x1) = %f' % (e1/nSamples))
    print('E(x6) = %f' % (e6/nSamples))
    print('P(x1 >= 285) = %d/%d' % (nq1, nSamples))
    print('P(x6 >= 285) = %d/%d' % (nq2, nSamples))

if __name__ == '__main__':
    simulateAll(100000)


