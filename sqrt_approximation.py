import matplotlib.pyplot as plt
import scipy
from scipy import special
import numpy as np
import cmath
import math

pi = 3.1415926

# complex number for exp(i*pi/4)
a = complex(1 / math.sqrt(2), 1 / math.sqrt(2))
# complex number for exp(i*3*pi/4)
b = complex(-1 / math.sqrt(2), 1 / math.sqrt(2))
# the complex number
i = complex(0, 1)


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


D = sqrt(3)
n = 10


def fourierSeries(x):
    sum = 0
    for i in range(1, n):
        phi_1 = (-1) * S(sqrt(2 * i)) / ((sqrt(i) ** 3) * sqrt(2 * D)) + sin(pi * i * D) / i
        phi_2 = (-1) * C(sqrt(2 * i)) / ((sqrt(i) ** 3) * sqrt(2 * D)) + cos(pi * i * D) / i
        sum += phi_1 * cos(i * pi * x / (D ** 2)) - phi_2 * sin(i * pi * x / (D ** 2))
    return 2 * sum + D


t1 = np.arange(0, D, .01)
out = np.vectorize(fourierSeries)
plt.plot(t1,out(t1))
plt.show()

