"""Exact line projection of E through a rational point; no descent or database."""
from fractions import Fraction as Q
import math
from point_arithmetic import on_curve


def invariants(a):
    a1,a2,a3,a4,a6=map(Q,a)
    return a1*a1+4*a2, a1*a3+2*a4, a3*a3+4*a6


def quartic(a, anchor):
    a=list(map(Q,a)); x0,y0=map(Q,anchor)
    if not on_curve(a,(x0,y0)): raise ValueError('Anchor off curve')
    b2,b4,_=invariants(a);v0=2*y0+a[0]*x0+a[2]
    return [b2*b2-8*b2*x0-48*x0*x0-32*b4,32*v0,-2*(12*x0+b2),Q(0),Q(1)]


def to_quartic(a, anchor, point):
    a=list(map(Q,a)); x0,y0=map(Q,anchor);x,y=map(Q,point)
    if not on_curve(a,(x,y)): raise ValueError('Point off curve')
    if x==x0: raise ValueError('Vertical/tangent exceptional point; handled separately')
    b2,_,_=invariants(a);v0=2*y0+a[0]*x0+a[2];v=2*y+a[0]*x+a[2]
    t=(v-v0)/(x-x0);z=8*x-t*t+b2+4*x0
    f=quartic(a,anchor)
    if z*z!=sum(c*t**i for i,c in enumerate(f)): raise ArithmeticError('Quartic identity')
    return t,z


def from_quartic(a, anchor, point):
    a=list(map(Q,a));x0,y0=map(Q,anchor);t,z=map(Q,point)
    if z*z!=sum(c*t**i for i,c in enumerate(quartic(a,anchor))): raise ValueError('Point off quartic')
    b2,_,_=invariants(a);v0=2*y0+a[0]*x0+a[2]
    x=(t*t-b2-4*x0+z)/8;y=(v0+t*(x-x0)-a[0]*x-a[2])/2
    if not on_curve(a,(x,y)): raise ArithmeticError('Inverse identity')
    return x,y


def small_points(a, numerator_bound=1024, denominator_bound=8):
    """Exact search x=n/d^2, y=m/d^3 with necessary square-residue filters.

    Increasing both bounds exhausts rational points on an integral model.
    The filters cannot discard a square. This finite implementation has no
    claim that the chosen bounds suffice for any particular rank.
    """
    a=list(map(int,a));b2,b4,b6=map(int,invariants(a))
    filters=[(p,{x*x%p for x in range(p)}) for p in (16,9,5,7,11)]
    for d in range(1,denominator_bound+1):
        for k in range(numerator_bound+1):
            for n in ((0,) if k==0 else (k,-k)):
                if math.gcd(n,d)!=1: continue
                value=4*n**3+b2*n*n*d*d+2*b4*n*d**4+b6*d**6
                if value<0 or any(value%p not in residues for p,residues in filters): continue
                z=math.isqrt(value)
                if z*z!=value: continue
                ynum=-z-a[0]*n*d-a[2]*d**3
                if ynum%2: continue
                point=Q(n,d*d),Q(ynum,2*d**3)
                if not on_curve(list(map(Q,a)),point): raise ArithmeticError('Denominator square identity')
                yield point,{'n':n,'d':d}
