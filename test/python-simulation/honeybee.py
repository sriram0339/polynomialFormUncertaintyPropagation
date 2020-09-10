from random import uniform, gauss

def truncGaussian( mean, sigma, a, b):
    r = gauss(mean, sigma)
    while (r < a or r > b):
        r = gauss(mean, sigma)
    return r

def simulateOnce(nSteps):
    x = truncGaussian(475, 5, 450, 500)
    y1 = uniform(350, 400)
    y2 = uniform(100, 150) 
    z1 = truncGaussian(35, 1.5, 20, 50)
    z2 = truncGaussian(35, 1.5, 20, 50)
    for j in range(nSteps):
        xn = x + 0.1 *(-0.001 * x * y1 - 0.001 * x * y2)
        y1n = y1 + 0.1 *( 0.001 * x * y1 - 0.3 * y1 + 0.5 * 0.001 * y1 * z1 + 0.7 * 0.001 * y1 * z2 )
        y2n = y2 + 0.1 * (0.001 * x * y2 - 0.3 * y2 + 0.5 * 0.001 * y2 * z2 + 0.7 * 0.001 * y2 * z1)
        z1n = z1 + 0.1 * (0.3 * y1 - 0.5 * 0.001 * y1 * z1 - 0.7 * 0.001 * y2 * z1)
        z2n = z2 + 0.1 * (0.3 * y2 - 0.5 * 0.001 * y2 * z2 - 0.7 * 0.001 * y1 * z2)
        (x,y1, y2, z1, z2) = (xn, y1n, y2n, z1n, z2n)
    return (z1, z2)



def simulateAll(nSamples):
    e1 = 0.0
    e2 = 0.0
    nq1 = 0
    nq2 = 0

    for j in range(nSamples):
        (z1, z2) = simulateOnce(25)
        e1 = e1 + z1
        e2 = e2 + z2
        if z1 >= 265:
            nq1 += 1
        if z2 <= 60:
            nq2 += 1
    print('E(z1) = %f' % (e1/nSamples))
    print('E(z2) = %f' % (e2/nSamples))
    print('P(z1 >= 265) = %d/%d' % (nq1, nSamples))
    print('P(z2 <= 60) = %d/%d' % (nq2, nSamples))

if __name__ == '__main__':
    simulateAll(1000000)

