"""Exact rational point arithmetic for basis changes, without rank algorithms."""
from fractions import Fraction as Q


def negate(a, p):
    return None if p is None else (p[0], -p[1]-a[0]*p[0]-a[2])


def add(a, p, q):
    if p is None:
        return q
    if q is None:
        return p
    x, y = p
    xx, yy = q
    if x == xx:
        if y+yy+a[0]*x+a[2] == 0:
            return None
        slope = (3*x*x+2*a[1]*x+a[3]-a[0]*y)/(2*y+a[0]*x+a[2])
    else:
        slope = (yy-y)/(xx-x)
    intercept = y-slope*x
    xr = slope*slope+a[0]*slope-a[1]-x-xx
    return xr, -(slope+a[0])*xr-intercept-a[2]


def multiply(a, p, n):
    if n < 0:
        return multiply(a, negate(a, p), -n)
    result = None
    while n:
        if n & 1:
            result = add(a, result, p)
        p = add(a, p, p)
        n >>= 1
    return result


def combine(a, points, coefficients):
    if len(points) != len(coefficients):
        raise ValueError("Wrong number of coefficients")
    result = None
    for p, n in zip(points, coefficients):
        if int(n) != Q(n):
            raise ValueError("Basis coefficients must be integers")
        result = add(a, result, multiply(a, p, int(n)))
    return result


def on_curve(a, p):
    if p is None:
        return True
    x, y = p
    return y*y+a[0]*x*y+a[2]*y == x**3+a[1]*x*x+a[3]*x+a[4]
