from random import uniform, gauss
from math import sqrt, pi

def truncGaussian( mean, sigma, a, b):
    r = gauss(mean, sigma)
    while (r < a or r > b):
        r = gauss(mean, sigma)
    return r

def simulateOnce(nSteps):
    x = uniform(-0.1, 0.1)
    theta = uniform(-0.1, 0.1)
    dxdt = 0.0
    dthetadt = 0.0
    delta = 0.1
    tau = 0.01 * sqrt(0.1)
    for j in range(nSteps):
        w1 =  tau * truncGaussian(0.0, 1.0, -10.0, 10.0)
        w2 = tau*truncGaussian(0.0, 1.0, -10.0, 10.0)
        w3 = tau * truncGaussian(0.0, 1.0, -10.0, 10.0)
        w4 = tau * truncGaussian(0.0, 1.0, -10.0, 10.0)
        feedback = 10.0 * x - 289.83 * theta + 19.53 * dxdt - 63.25*dthetadt
        ddx = 0.1 * feedback + 0.98 * theta
        tmp0 = theta *(-0.75 * theta - 0.1 * feedback)
        #tmp0 = tmp0 * theta
        # thetaSq = theta ** 2
        tmp1 = -0.05 * dthetadt**2 * theta**2
        ddx = tmp0 + tmp1
        ddtheta = 0.2 * feedback + 21.56* theta
        tmp2 = (-5.75*theta - 0.12 * feedback)* theta**2
        tmp3 = -0.1 * dthetadt**2
        ddtheta += tmp2 + tmp3
        x = x + delta * dxdt + w1
        theta = theta + delta * dthetadt + w2
        dxdt = dxdt + delta * ddx + w3
        dthetadt = dthetadt + delta * ddtheta + w4
    return (x, theta)
 
def simulateAll(nSamples):
    ex = 0.0
    etheta = 0.0
    qx = 0
    qtheta = 0
    for j in range(nSamples):
        (x, th) = simulateOnce(8)
        ex = ex + x
        etheta = etheta + th
        if x >= 2.0:
            qx += 1
        if th >= pi/6:
            qtheta += 1
    print('E(x) = %f' %(ex/nSamples))
    print('E(theta) = %f' % (etheta/nSamples))
    print('P(x>= 2.0) = %d/%d' %(qx, nSamples))
    print('P(theta>= pi/6) = %d/%d' %(qtheta, nSamples))


if __name__ == '__main__':
    simulateAll(100000)





