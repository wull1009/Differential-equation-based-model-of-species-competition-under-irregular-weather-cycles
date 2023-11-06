import numpy
from scipy.integrate import odeint
from matplotlib import pyplot
import math


def f(t, u):
    x1, x2, x3, x4, x5, x6, y = u
    r1 = 0.22 * 0.26
    r2 = 0.75 * 0.13
    r3 = 0.67 * 0.05
    r4 = 0.35 * 0.34
    r5 = 0.6 * 0.62
    r6 = 0.56 * 0.56
    k = 1
    dx1dt = -r1 * (1 - x1 / 500) * x1 + (
            -1.784875 * 2 * x1 * x2 / 3 - 1.96301 * 2 * x1 * x3 / 3 + 1.674715 * x1 * x4 / 3 + 0.307152 * x1 * x5 + 0.645408 * x1 * x6) * 0.0002825 + 0.01 * r1 * x1 / (1 + math.exp(-k * (y - 10)))
    dx2dt = -r2 * (1 - x2 / 500) * x2 + (
            + 1.784875 * x1 * x2 / 3 + 0.260596 * x2 * x3 + 3.08484 * x2 * x4 / 3 + 0.952032 * x2 * x5 + 0.932085 * x2 * x6) * 0.0005825 + 0.01 * r2 * x2 / (1 + math.exp(-k * (y - 10)))
    dx3dt = -r3 * (1 - x3 / 500) * x3 + (
            + 1.96301 * x1 * x3 / 3 + 0.260596 * x2 * x3 + 2.740332 * x3 * x4 / 3 + 2.372211 * x3 * x5 / 3 + 2.015143 * x3 * x6 / 3) * 0.000825 + 0.01 * r3 * x3 / (1 + math.exp(-k * (y - 10)))
    dx4dt = r4 * (1 - x4 / 500) * x4 + (
            - 1.674715 * x1 * x4 * 2 / 3 - 3.08484 * x2 * x4 * 2 / 3 - 2.740332 * x3 * x4 * 2 / 3 + x5 * x4 * 0.377994 + x6 * x4 * 0.137895) * 0.00002825 + 0.01 * r4 * x4 / (1 + math.exp(-k * (y - 10)))
    dx5dt = r5 * (1 - x5 / 300) * x5 + (
            0.307152 * x1 * x5 + 0.952032 * x2 * x5 - 2.372211 * x3 * x5 * 2 / 3 + x4 * x5 * 0.377994 + x6 * x5 * 0.217119) * 0.00000825 + 0.01 * r5 * x5 / (1 + math.exp(-k * (y - 10)))
    dx6dt = r6 * (1 - x6 / 300) * x6 + (
            + 0.645408 * x1 * x6 + 0.932085 * x2 * x6 - 2.015143 * x3 * x6 * 2 / 3 + x6 * x4 * 0.137895 + x6 * x5 * 0.217119) * 0.00000825 - 0.01 * r6 * x6 / (1 + math.exp(-k * (y - 10)))
    dydt = 0.004306 * 5 * (t ** 4) - 0.1814 * 4 * (t ** 3) + 2.77 * 3 * (t ** 2) - 18.59 * 2 * t + 52.1
    dudt = [dx1dt, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt, dydt]
    return dudt


def main():
    year = 15
    u0 = [10, 58, 22, 42, 14, 12, 10]
    t = numpy.linspace(1, 15, 15)
    result = odeint(f, u0, t, tfirst=True)
    print(result)
    pyplot.plot(range(0, year), result[:, 0], 'r', label="ACWR")
    pyplot.plot(range(0, year), result[:, 1], 'g', label='FEMI')
    pyplot.plot(range(0, year), result[:, 2], 'b', label="HOMU")
    pyplot.plot(range(0, year), result[:, 3], 'c', label="SACO")
    pyplot.plot(range(0, year), result[:, 4], 'm', label="PLER")
    pyplot.plot(range(0, year), result[:, 5], 'y', label="URLI")
    pyplot.xlabel("year")
    pyplot.ylabel("number of each species")
    pyplot.legend()
    pyplot.grid()
    pyplot.show()


if __name__ == "__main__":
    main()
