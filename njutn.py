import numpy


def f(x):
    return x ** 3 + 0.1 - numpy.sin(x)


def df(x):
    return 3 * x**2 - numpy.cos(x)


def ddf(x):
    return 6 * x + numpy.sin(x)


def newton(a, b, precision):

    current_x = b
    if ddf(a) * f(a) > 0:
        current_x = a
    last_x = 10 * current_x
    while(numpy.abs(last_x - current_x) > precision):
        last_x = current_x
        current_x = last_x - f(last_x) / df(last_x)

    return current_x


print(newton(0.67, 1, 10**(-4)))
