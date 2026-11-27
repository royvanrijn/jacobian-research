"""Directed sharp scalar comparison with the undamped angular variance gain.

The variance gain is (6 rho d / (25 pi)) Var_{mu_1}, d=1-rho^2.
Its two changes are a diagonal increment and a positive update to mode one.
"""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import factorial
from pathlib import Path
import sys

from flint import arb, arb_mat, ctx
import numpy as np


folder = Path(__file__).resolve().parent
if not (folder/"check_three_cell_global_scalar_through_half.py").is_file():
    folder = folder.parent/"research/plurality-stability"
sys.path.insert(0,str(folder))
from check_three_cell_global_scalar_through_half import radial_coefficients, upper_float
from check_three_cell_sharp_positive_arc_envelope import sharp_moments


def enclose(q):
    return arb(q.numerator)/q.denominator


class Quadrature:
    """Exact monomial integration of fixed piecewise-linear basis functions."""
    def __init__(self, panels, radius, order):
        self.radius=arb(radius)
        self.width=self.radius/panels
        self.nodes=[j*self.width for j in range(panels+1)]
        self.order=order
        self.endpoint_weights=arb_mat(panels+1,order+1)
        self.panel_weights=arb_mat(panels,order+1)
        for i,(left,right) in enumerate(zip(self.nodes,self.nodes[1:])):
            left_power,right_power=left,right
            mass=right-left
            for n in range(order+1):
                left_power*=left
                right_power*=right
                moment=(right_power-left_power)/(n+2)
                left_weight=(right*mass-moment)/self.width
                right_weight=(moment-left*mass)/self.width
                self.endpoint_weights[i,n]+=left_weight
                self.endpoint_weights[i+1,n]+=right_weight
                self.panel_weights[i,n]=mass
                mass=moment

    def radial_moments(self,d,kappa,beta,row):
        orders=(self.order,self.order//2,0)
        exponents=(beta,arb(1),1-kappa)
        assert all(exponent>0 for exponent in exponents)
        node_values=[[] for _ in exponents]
        panel_errors=[[] for _ in exponents]
        for j,x in enumerate(self.nodes):
            polynomial, first, second = row[-1],arb(0),arb(0)
            for coefficient in reversed(row[:-1]):
                second=second*x+2*first
                first=first*x+polynomial
                polynomial=polynomial*x+coefficient
            assert polynomial>0
            values=[8*d*(-exponent*x**2).exp()/polynomial for exponent in exponents]
            if j:
                log_slope=first/left_polynomial
                curvature=second/left_polynomial
                for errors,exponent,left_value in zip(panel_errors,exponents,left_values):
                    second_bound=left_value*((2*exponent*x+log_slope)**2
                                             +2*exponent+curvature+log_slope**2)
                    error=(second_bound*self.width**2/8).upper()
                    errors.append((-error).union(error))
            for output,value in zip(node_values,values):
                output.append(value)
            left_polynomial,left_values=polynomial,values

        radius=self.radius
        leading=max(range(len(row)),key=lambda j:upper_float(row[j]*radius**j))
        result=[]
        for values,errors,order,exponent in zip(node_values,panel_errors,orders,exponents):
            values_row=arb_mat(1,len(values),values)
            errors_row=arb_mat(1,len(errors),errors)
            enclosure=values_row*self.endpoint_weights+errors_row*self.panel_weights
            moments=[]
            for n in range(order+1):
                power=n-leading-1
                decay=exponent-arb(max(0,power))/(2*radius**2)
                assert decay>0
                tail=4*d*radius**power*(-exponent*radius**2).exp()/(row[leading]*decay)
                moments.append(enclosure[0,n]+arb(0).union(tail))
            result.append(moments)
        tt,uu,(ww,)=result
        assert ww>0
        return tt,uu,ww


def gram_bounds(modes,tt,uu,ww):
    size=len(modes)
    roots=[mode.sqrt() for mode in modes]
    gram=arb_mat(size,size)
    hs2=arb(0)
    for i in range(size):
        for j in range(i+1):
            entry=roots[i]*roots[j]*(tt[i+j]-uu[i]*uu[j]/ww)
            gram[i,j]=gram[j,i]=entry
            hs2+=(1 if i==j else 2)*entry**2

    midpoint=np.array([[float(gram[i,j].mid()) for j in range(size)] for i in range(size)])
    _,vectors=np.linalg.eigh(midpoint)
    trial=arb_mat(size,1,[arb(round(float(x)*2**44))/2**44 for x in vectors[:,-1]])
    image=gram*trial
    norm=(trial.transpose()*trial)[0,0]
    rayleigh=(trial.transpose()*image)[0,0]/norm
    residual=image-trial*rayleigh
    residual2=sum((residual[i,0]**2 for i in range(size)),arb(0))/norm
    bound_second=(hs2-rayleigh**2).upper().sqrt()
    assert rayleigh>bound_second
    return rayleigh,rayleigh+residual2/(rayleigh-bound_second)


def bound(lower,upper,angular,row_order,quadrature):
    rho=enclose(lower).union(enclose(upper))
    d=1-rho**2
    kappa=rho**2/d
    beta=1+kappa
    rate=2*rho/d
    row,modes=radial_coefficients(d,beta,rate,*angular,row_order)
    plain_row=row[0]
    tilt=2/(arb.pi().sqrt()*d*d.sqrt())
    gain=6*rho*d/(25*arb.pi())
    row[0]+=2*d*gain*tilt
    modes[1]+=gain*tilt**2
    rayleigh,finite=gram_bounds(modes,*quadrature.radial_moments(d,kappa,beta,row))

    n=len(modes)
    small_ratio=4*rate/(n+1)
    eta=1+arb(49)/25*kappa
    outer_ratio=rate/eta
    assert 0<small_ratio<1 and 0<outer_ratio<1
    small=16*(4*rate)**n/(plain_row*factorial(n)*(1-small_ratio))
    large=4096*outer_ratio**n/(55*d*eta*(1-outer_ratio))
    return rayleigh,finite,small+large


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--points",type=Q,nargs="+")
    p.add_argument("--start",type=Q,default=Q(29,50))
    p.add_argument("--stop",type=Q,default=Q(31,50))
    p.add_argument("--width",type=Q,default=Q(1,1000))
    p.add_argument("--panels",type=int,default=2048)
    p.add_argument("--modes",type=int,default=144)
    p.add_argument("--rows",type=int,default=72)
    p.add_argument("--radius",type=int,default=12)
    p.add_argument("--cap",type=Q,default=Q(1999,1000))
    args=p.parse_args()
    ctx.prec=110
    assert (arb.pi()/16).cos()**2>arb(24)/25
    assert (2*arb.pi()/9).cos()**2<arb(3)/5
    assert (arb(18)/25).exp()>2
    start=min(args.points) if args.points else args.start
    assert start>=Q(29,50)
    assert enclose(start)**2/(1-enclose(start)**2)>arb(1)/2
    assert 0<args.rows<=args.modes
    angular=sharp_moments(args.modes)
    quadrature=Quadrature(args.panels,args.radius,2*args.modes)
    inputs=(Path(__file__),folder/"check_three_cell_global_scalar_through_half.py",
            folder/"check_three_cell_sharp_positive_arc_envelope.py")
    settings=dict(precision_bits=ctx.prec,mode_order=args.modes,row_order=args.rows,
                  panels=args.panels,radial_limit=args.radius,cap=str(args.cap),
                  source_sha256={p.name:sha256(p.read_bytes()).hexdigest() for p in inputs})

    if args.points:
        for point in args.points:
            rayleigh,finite,tail=bound(point,point+args.width,angular,args.rows,quadrature)
            print(json.dumps(dict(start=str(point),width=str(args.width),rayleigh_lo=float(rayleigh.lower()),
                                  finite_upper=upper_float(finite),tail_upper=upper_float(tail),
                                  total_upper=upper_float(finite+tail))),flush=True)
        return
    count=(args.stop-args.start)/args.width
    assert count.denominator==1 and count>0
    slabs=[]
    def certify(left,right,depth=0):
        _,finite,tail=bound(left,right,angular,args.rows,quadrature)
        if finite+tail<enclose(args.cap):
            slabs.append(dict(rho=[str(left),str(right)],finite_upper=upper_float(finite),
                              tail_upper=upper_float(tail),total_upper=upper_float(finite+tail)))
            return
        assert depth<8,(str(left),str(right),upper_float(finite+tail))
        midpoint=(left+right)/2
        certify(left,midpoint,depth+1)
        certify(midpoint,right,depth+1)

    for index in range(count.numerator):
        certify(args.start+index*args.width,args.start+(index+1)*args.width)
    print(json.dumps(dict(settings,status="pass",rho_bounds=[str(args.start),str(args.stop)],
                          total_upper=max(s["total_upper"] for s in slabs),slabs=slabs),indent=2))


if __name__=="__main__":
    main()
