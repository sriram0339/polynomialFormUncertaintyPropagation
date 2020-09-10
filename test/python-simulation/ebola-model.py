from random import uniform, gauss

def truncGaussian( mean, sigma, a, b):
    r = gauss(mean, sigma)
    while (r < a or r > b):
        r = gauss(mean, sigma)
    return r

def simulateOnce(nSteps):
    s = truncGaussian(0.7,0.02,0.6,0.8)
    e = uniform(0.2, 0.4)
    i = uniform(0.0, 0.04)
    r = uniform(0.0, 0.04)
    c = uniform(0.0, 0.04)
    for j in range(nSteps):
        sn = s - (s * 0.35 * i) * 0.5
        en = e + ( (s * 0.35 * i) - (0.28)*e) * 0.5
        ine = i + (0.28 * e - 0.29 * i) * 0.5
        rn = r + ( 0.29 * i) * 0.5
        cn = c + 0.28 * e * 0.5
        s = sn
        i = ine
        r = rn
        e = en
        c = cn
    return (i, e)

def simulateAll(nSamples):
    ei = 0.0
    ee = 0.0
    nq1 = 0
    nq2 = 0

    for j in range(nSamples):
        (i, e) = simulateOnce(25)
        ei = ei + i
        ee = ee + e
        if i >= 0.13:
            nq1 += 1
        if i <= 0.05:
            nq2 += 1
    print('E(i) = %f' % (ei/nSamples))
    print('E(e) = %f' % (ee/nSamples))
    print('P(i >= 0.13) = %d/%d' % (nq1, nSamples))
    print('P(i <= 0.05) = %d/%d' % (nq2, nSamples))

if __name__ == '__main__':
    simulateAll(100000)
