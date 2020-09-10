from random import uniform

def simulateOnce(nSteps):
    u1 =  uniform(-0.5, -0.3)
    u2  = uniform (0.4, 0.5)
    u3 = uniform (-0.2, 0)
    u4  = uniform (-0.2, 0)
    u5 = uniform( 0.55, 0.75)
    u6  = uniform (0.1, 0.3)
    u7  = uniform(0.55, 0.75)
    u8  = uniform(-0.19, 0.19)
    u9  = uniform(-0.6, -0.4)
    u10  = uniform(-0.19, 0.19)
    for i in range(nSteps):
        u1n = u1 + 0.1 *( 0.1 * (u2 - u1)  - u1 * (u1 - 1) * (u1 - 0.6))
        u2n = u2 + 0.1 *( 0.1 * (u1 + u3 - 2 * u2)  - u2 * (u2 - 1) * (u2 - 0.6))
        u3n = u3 + 0.1 *( 0.1 * (u2 + u4 - 2 * u3)  - u3 * (u3 - 1) * (u3 - 0.6))
        u4n = u4 + 0.1 *( 0.1 * (u5 + u3 - 2 * u4)  - u4 * (u4 - 1) * (u4 - 0.6))
        u5n = u5 + 0.1 *( 0.1 * (u6 + u4 - 2 * u5)  - u5 * (u5 - 1) * (u5 - 0.6))
        u6n = u6 + 0.1 *( 0.1 * (u7 + u5 - 2* u6)  - u6 * (u6 - 1) * (u6 - 0.6))
        u7n = u7 + 0.1 *( 0.1 * (u8 + u6 - 2* u7)  - u7 * (u7 - 1) * (u7 - 0.6))
        u8n = u8 + 0.1 *( 0.1 * (u9 + u7 - 2* u8)  - u8 * (u8 - 1) * (u8 - 0.6)) 
        u9n = u9 + 0.1 *( 0.1 * (u10 + u8 - 2 *u9)  - u9 * (u9 - 1) * (u9 - 0.6))
        u10n = u10 + 0.1 *( 0.1 * (u9 - u10)  - u10 * (u10 - 1) * (u10 - 0.6))
        u1 = u1n
        u2 = u2n
        u3 = u3n
        u4 = u4n
        u5 = u5n
        u6 = u6n
        u7 = u7n
        u8 = u8n
        u9 = u9n 
        u10 = u10n
    return (u1, u8)

def runSamples(nSamples):
    eu1= 0.0
    eu8 = 0.0
    u1_max = -1000000000.0
    u1_min = 1000000000.0
    num_u1_event = 0
    num_u8_event = 0
    for j in range(nSamples):
        (u1,u8) = simulateOnce(15)
        eu1 = eu1 + u1
        eu8 = eu8 + u8
        u1_max, u1_min = max(u1, u1_max), min(u1, u1_min)
        if u1 >= 0.0:
            num_u1_event = num_u1_event + 1
        if u8 <= 0:
            num_u8_event = num_u8_event + 1
    print(' E(u1) = %f ' % (eu1/nSamples))
    print(' rng(u1) = [%f, %f]' %(u1_min, u1_max))
    print(' E(u8) = %f ' % (eu8/nSamples))
    print(' P(u1 >= 0) = %d/%d' % (num_u1_event, nSamples))
    print(' P(u8 <= 0) = %d/%d' % (num_u8_event, nSamples))

if __name__ == '__main__':
    runSamples(100000)

