import math
import scipy
from scipy import special
import cmath

pi = 3.1415926

# complex number for exp(i*pi/4)
a = complex(1 / math.sqrt(2), 1 / math.sqrt(2))
# complex number for exp(i*3*pi/4)
b = complex(-1 / math.sqrt(2), 1 / math.sqrt(2))
# the complex number
i = complex(0, 1)

# Dimensionality of space
N = 3
# lengths of the rectangular region
d = [1]*N
# order of approximation
n = 1000


def S(x):
    output = scipy.special.fresnel(x)
    return output[0]


def C(x):
    output = scipy.special.fresnel(x)
    return output[1]


def sqrt(x):
    if x >= 0:
        return math.sqrt(x)
    else:
        return i * math.sqrt(-x)


def cos(x):
    return math.cos(x)


def sin(x):
    return math.sin(x)


def diagonal(d):
    sum = 0
    for x in d:
        sum += x ** 2
    return math.sqrt(sum)


def area(d):
    prod = 1
    for x in d:
        prod *= x
    return prod

# diagonal across the rectangular region
D = diagonal(d)
# area of the rectangular region
AD = area(d)

def convertBinary(N, n):
    binary = [0] * N
    x = n
    if n == 0:
        return binary
    else:
        i = 0
        while (x > 0):
            binary[i] = x % 2
            x = math.floor(x / 2)
            i += 1
        return binary


def size(d):
    i = 0
    for x in range(N):
        if d[x] == 1:
            i += 1
    return i

def product(N, n, d, bin):
    prod = 1
    for k in range(N):
        if (bin[k] == 0):
            prod *= sqrt(2) * C(d[k] * sqrt(2 * n) / (D)) - D * sin(pi * n * (d[k] ** 2) / (D ** 2)) / (pi * sqrt(n))
        else:
            prod *= sqrt(2) * S(d[k] * sqrt(2 * n) / (D)) - 2 * D * (sin(pi * n * (d[k] ** 2) / ((D ** 2) * 2)) ** 2) / (pi * sqrt(n))
    return prod

def sum(N, n, d, q):
    sum = 0
    N_set = []
    for x in range(1, N + 1):
        N_set.append(x)

    binaryPset = [[] * N] * (2 ** N)
    for x in range(2 ** N):
        binaryPset[x] = convertBinary(N, x)

    if (q == 0):
        for x in range(2 ** N):
            if (size(binaryPset[x]) % 2 == 0):
                sum += ((-1) ** (size(binaryPset[x]) / 2)) * product(N, n, d, binaryPset[x])
        return sum
    else:
        for x in range(2 ** N):
            if (size(binaryPset[x]) % 2 == 1):
                sum += ((-1) ** ((size(binaryPset[x]) - 1) / 2)) * product(N, n, d, binaryPset[x])
        return sum


def fourierSeries(N, n, d):
    approx = 0
    for x in range(1, n):
        phi_1 = (-1) * S(sqrt(2 * x)) / ((sqrt(x) ** 3) * sqrt(2 * D)) + sin(pi * x * D) / x
        phi_2 = (-1) * C(sqrt(2 * x)) / ((sqrt(x) ** 3) * sqrt(2 * D)) + cos(pi * x * D) / x
        approx += phi_1 * sum(N, n, d, 0) - phi_2 * sum(N, n, d, 1)
    return approx


result = (D**N)*fourierSeries(N, n, d)/(pi*(AD**2))+D
print(result)
