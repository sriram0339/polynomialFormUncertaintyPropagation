from random import uniform

def simulateOnce(nSteps):
    x1  = uniform(1.0, 1.2)
    y1 = uniform(0.85, 0.95)
    x2 = uniform(1.3, 1.4)
    y2 = uniform(2.1, 2.3)
    x3 = uniform(-0.1, 0.1)
    y3 = uniform(0.3, 0.5)

    for i in range(nSteps):
        x1n = x1  + 0.05 * y1 + uniform(-0.01,0.01)
        y1n = y1 + 0.05 * (0.5 * (1.0 - x1*x1) * y1 - x1 + 0.05 * x2 )
        x2n = x2 + 0.05 * y2 + uniform(-0.01,0.01)
        y2n = y2 + 0.05 * (0.33 * (1.0 - x2*x2) * y2 - x2 + 0.05 * x3 )
        x3n = x3 + 0.05 * y3 + uniform(-0.01,0.01)
        y3n = y3 + 0.05 * (0.45 * (1.0 - x3*x3) * y3 - x3 )
        x1 = x1n
        x2 = x2n
        x3 = x3n
        y1 = y1n
        y2 = y2n
        y3 = y3n
    return y3

def simulateAll(nSamples):
    ey3 = 0
    maxy3 = -1000000.0
    miny3 = 1000000.0
    nq1 = 0
    nq2 = 0
    for i in range(nSamples):
        y3 = simulateOnce(15)
        ey3 = ey3 + y3
        miny3, maxy3 = min(miny3, y3), max(maxy3, y3)
        if y3 >= 0.6:
            nq1 = nq1 + 1
        if y3 <= 0.2:
            nq2 = nq2 + 1
    print('E(y3) = %f' % (ey3/nSamples))
    print('rng(y3) = [%f, %f]' % (miny3, maxy3))
    print('P(y3 >= 0.6) = %d/ %d' % (nq1, nSamples))
    print('P(y3 <= 0.2) = %d/ %d' % (nq2, nSamples))

if __name__ == '__main__':
    simulateAll(100000)
