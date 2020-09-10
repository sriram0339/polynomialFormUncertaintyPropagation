from random import uniform, gauss

def truncGaussian( mean, sigma, a, b):
    r = gauss(mean, sigma)
    while (r < a or r > b):
        r = gauss(mean, sigma)
    return r

def simulateOnce(nSteps):
    x1 = uniform (1.0, 1.2)
    x2 = uniform (0.85, 1.05)
    x3 = uniform (1.3, 1.5)
    x4 = uniform (2.25, 2.55)
    x5 = uniform (0.4, 0.7)
    x6 = uniform (-0.2, 0.2)
    x7 = uniform (0.3, 0.55)
    for j in range(nSteps):
        x1n = x1 + 0.14 * x3 - 0.09 * x1 
        x2n = x2 + 0.25 * x5 - 0.15 * x2 
        x3n = x3 + 0.06 * x7 - 0.08 * x2 * x3 
        x4n = x4 + 0.2 - 0.13 * x3 * x4 
        x5n = x5 + 0.07 * x1 - 0.1 * x4 * x5 
        x6n = x6 + 0.03 * x1 - 0.31 * x6 
        x7n = x7 + 0.18 * x6 - 0.15 * x2* x7
        (x1, x2, x3, x4, x5, x6, x7) = (x1n, x2n, x3n, x4n, x5n, x6n, x7n)
    return (x1, x2)

def simulateAll(nSamples):
    e1 = 0.0
    e2 = 0.0
    nq1 = 0
    nq2 = 0

    for j in range(nSamples):
        (x1, x2) = simulateOnce(25)
        e1 = e1 + x1
        e2 = e2 + x2
        if x1 <= 0.7:
            nq1 += 1
        if x2 >= 0.95:
            nq2 += 1
    print('E(x1) = %f' % (e1/nSamples))
    print('E(x2) = %f' % (e2/nSamples))
    print('P(x1 <= 0.7) = %d/%d' % (nq1, nSamples))
    print('P(x2 >= 0.95) = %d/%d' % (nq2, nSamples))

if __name__ == '__main__':
    simulateAll(100000)




