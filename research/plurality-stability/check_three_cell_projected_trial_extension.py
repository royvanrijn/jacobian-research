"""Directed projected-trial comparison through correlation 2589/4000.

The proof and replay command are in THREE_CELL_PROJECTED_TRIAL_EXTENSION.md.
The finite operator, anchored trial, and degree-one variance credit are
reused from the excess-purity checker.
"""

import argparse
from bisect import bisect_left
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from flint import arb, ctx
from math import isfinite

from check_three_cell_excess_purity_extension import (
    Quadrature, block, enclose, geometry, polynomial, sharp_moments, upper_float,
)


def median_geometry(lo,hi,row,coefficients,center,norm,kappa):
    d_min=1-hi**2
    d_max=1-lo**2
    beta=1+kappa
    gamma1=arb.pi().sqrt()/2
    gamma3=3*gamma1/2
    norm1=gamma1/(beta*beta.sqrt())
    densities=[8*arb.pi()/(lo*d_min*(arb(11)/2+arb(24)/25)),None]
    qcoef=row[2]+2*(row[0]*row[4]).sqrt()+2*(row[1]*row[3]).sqrt()
    densities[1]=4/(d_min*d_min.sqrt()*gamma3*qcoef)
    limit=arb(4)
    width=limit/2400
    mass3,second3,mean3=arb(0),arb(0),arb(0)
    tail_nodes=[]
    for i in range(2400):
        x=(i*width).union((i+1)*width)
        square=x**2
        q=polynomial(row,x)
        numerator=(-kappa*square).exp()*polynomial(coefficients,x)-center
        mass1=2/norm1*square*(-beta*square).exp()*width
        ratio=2*square/(3*d_min)
        z=(4*d_min/norm.sqrt())*(polynomial(coefficients,x)-center*(kappa*square).exp())/q
        upper=float("inf") if i==0 else upper_float(abs(z)/(i*width))
        tail_nodes.append((upper,mass1.upper(),(mass1*ratio).upper()))
        mass3+=mass1*ratio
        second3+=ratio*32*d_min**2/(norm*norm1)*(-(1-kappa)*square).exp()*(numerator/q)**2*width
        mean3+=ratio*8*d_min/(norm.sqrt()*norm1)*x*(-square).exp()*numerator/q*width
    variance3=(d_min/d_max*(mass3*second3-mean3**2)).lower()
    tail_nodes.sort(key=lambda v:v[0])
    y=beta.sqrt()*limit
    tail1=y.erfc()+2*y*(-y**2).exp()/arb.pi().sqrt()
    tail3=tail1+4*y**3*(-y**2).exp()/(3*arb.pi().sqrt())
    totals=[(tail1.upper(),tail3.upper())]
    for _,w1,w3 in reversed(tail_nodes):
        last1,last3=totals[-1]
        totals.append(((last1+w1).upper(),(last3+w3).upper()))
    clip_prefix=[arb(0)]
    for value,w1,_ in tail_nodes:
        if isfinite(value):
            cap_moment=(clip_prefix[-1]+arb(value)**2*w1).upper()
        else:
            cap_moment=arb("inf")
        clip_prefix.append(cap_moment)
    return densities,variance3,[v[0] for v in tail_nodes],list(reversed(totals)),clip_prefix


def clip_components(amplitude,pole,slope,radius,origin):
    if amplitude==0:
        return arb(0),arb(0)
    decay=slope+(arb(2)/3).sqrt()/amplitude
    endpoint=pole/decay
    if endpoint<radius:
        first=2*origin*pole**3/(3*decay)
        second=origin*pole**3*(3*decay-slope)/(3*decay**2)
    else:
        x=min(endpoint.lower(),radius.lower())
        first=2*origin*(pole**2*x-pole*decay*x**2+decay**2*x**3/3)
        second=2*origin*(pole**2*x-pole*(decay+slope)*x**2/2+decay*slope*x**3/3)
    return max(arb(0),first.lower()),max(arb(0),second.lower())


def amplitude_upper(r,pole,slope,radius,origin,projected=False):
    circ=(arb(2)/3).sqrt()
    denom=2*origin*pole**3-3*slope*r**2
    cutoff=3*r**2/(2*origin*pole**2)
    if not projected and denom>0 and cutoff<=radius:
        return (3*circ*r**2/denom).upper()

    def comparison(amplitude):
        if projected:
            first,second=clip_components(amplitude,pole,slope,radius,origin)
            return first+second**2>r**2
        decay=slope+circ/amplitude
        z=pole-decay*radius
        positive=max(z.upper(),arb(0))
        return 2*origin*(pole**3-positive**3)/(3*decay)>r**2

    upper=arb(1)
    if not comparison(upper):
        return None
    lower=arb(0)
    for _ in range(18):
        mid=(lower+upper)/2
        if comparison(mid):
            upper=mid
        else:
            lower=mid
    return upper


def median_lower(amplitude,r,lo,hi,variance1,median,error_residual=None):
    densities,variance3,nodes,totals,_=median
    variances=[variance1,variance3]
    coeff1lo=hi*(1-hi**2)/(2*arb.pi())
    coeff1hi=lo*(1-lo**2)/(2*arb.pi())
    base3lo=8*lo**3*(1-lo**2)*9*arb.pi()/16
    base3hi=8*hi**3*(1-hi**2)*9*arb.pi()/16
    coeff3=arb(7)/1000
    residual3=arb(7)/(54*arb.pi()**2)-coeff3
    weightslo=[coeff1lo,coeff3*base3lo]
    weightshi=[coeff1hi,coeff3*base3hi]
    var_lo=[]
    for B,V in zip(densities,variances):
        spread=V.sqrt()-(B-V).sqrt()*r
        var_lo.append(spread**2 if spread>0 else arb(0))
    purity3=arb(2)/3-amplitude**2*densities[1]*(1+r**2)
    total=weightslo[0]*var_lo[0]/50+weightslo[1]*var_lo[1]/2
    if purity3>0:
        total+=residual3*base3lo*purity3*var_lo[1]
    cost=None
    error=r if error_residual is None else error_residual
    for u in [arb(1)/100,arb(1)/50,arb(1)/30,arb(1)/25,arb(1)/20,arb(1)/15,arb(2)/25,arb(1)/10,arb(3)/25,arb(7)/50]:
        threshold=arb(6).sqrt()*(arb(1)/6-u)/amplitude
        position=bisect_left(nodes,-upper_float(-threshold))
        charges=[]
        for B,omega,w in zip(densities,weightshi,totals[position]):
            m=B*amplitude*error**2/(6*u)+(B*w/2).sqrt()*error
            m=min(m.upper(),(B/2).sqrt().upper()*error)
            charges.append(5*omega*m**2)
        charge=sum(charges,arb(0)).upper()
        cost=charge if cost is None else min(cost,charge)
    return (total-cost).lower()


def quartic_excess(raw_upper,coefficient,spread,circumradius,chi,fourth):
    purity=arb(2)/3-raw_upper
    if fourth:
        mean_bound=(chi*raw_upper/(1+chi)).sqrt().upper()
        quartic_factor=1-3*raw_upper
        if mean_bound<circumradius/2:
            quartic_factor=arb(4)/3-3*raw_upper-2*circumradius*mean_bound+2*mean_bound**2
        purity=max(purity.lower(),quartic_factor.lower())
    extra=purity-arb(12)/25
    return (coefficient*extra*spread**2).lower() if extra>0 and spread>0 else arb(0)


def residual_gate(e,residual,delta,lo,hi,geometry,steps=128,median=None):
    assert delta>0
    discriminant=residual**2+delta*e
    if discriminant<0:
        return True, upper_float(e+residual**2/delta), 0.0
    maximum=((residual+discriminant.sqrt())/delta).upper()
    pole,slope,radius,origin,variance,raw=geometry
    assert pole>0 and variance>0 and raw>0
    d_min=1-hi**2
    density=(8*arb.pi()/(lo*d_min*(arb(11)/2+arb(24)/25))).upper()
    coefficient=hi*d_min/(2*arb.pi())
    circumradius=arb(2).sqrt()/arb(3).sqrt()
    chi=(arb.pi()*d_min*d_min.sqrt()/(2*(1-hi**2/d_min).sqrt())-1).upper()
    maximum_gap=-1e5
    worst=0.0

    for number in range(steps):
        a=number*maximum/steps
        b=(number+1)*maximum/steps
        turning=residual/delta
        if turning<a:
            scalar=e+2*residual*a-delta*a**2
        elif turning>b:
            scalar=e+2*residual*b-delta*b**2
        else:
            scalar=e+residual**2/delta
        if scalar<0:
            continue
        amplitude=amplitude_upper(b,pole,slope,radius,origin,projected=True)
        if amplitude is None:
            return False, None, float(b)
        if raw<density/(1+b**2):
            second=(raw.sqrt()+b*(density-raw).sqrt())**2
        else:
            second=density*(1+b**2)
        spread=variance.sqrt()-(density-variance).sqrt()*b
        excess=quartic_excess((amplitude**2*second).upper(),coefficient,spread,circumradius,chi,True)
        best=excess.lower()
        if median is not None and not scalar<best:
            best=max(best,median_lower(amplitude,b,lo,hi,variance,median))
        if not scalar<best:
            lower_by_amplitude=[]
            for k in range(16):
                small=k*amplitude/16
                large=(k+1)*amplitude/16
                first,_=clip_components(small,pole,slope,radius,origin)
                if first>b**2:
                    continue
                eff=(b**2-first).upper().sqrt()
                raw_proj=(large*(raw.sqrt()+density.sqrt()*eff))**2
                if median is not None:
                    _,_,nodes,totals,prefix=median
                    index=bisect_left(nodes,float(circumradius/large))
                    clip_raw=(large**2*prefix[index]+arb(2)/3*totals[index][0]).upper()
                    raw_proj=min(raw_proj.upper(),((clip_raw.sqrt()+density.sqrt()*eff*large)**2).upper())
                raw_upper=min(raw_proj.upper(),(large**2*second).upper())
                extra=quartic_excess(raw_upper,coefficient,spread,circumradius,chi,True)
                lower=max(best,extra.lower())
                if median is not None and not scalar<lower:
                    lower=max(lower,median_lower(large,b,lo,hi,variance,median,error_residual=eff))
                lower_by_amplitude.append(lower)
            if lower_by_amplitude:
                best=max(best,min(lower_by_amplitude))
        gap=upper_float(scalar-best)
        if gap>maximum_gap:
            maximum_gap=gap
            worst=float(b)
    return maximum_gap<0,maximum_gap,worst


def certificate(point=None):
    ctx.prec=120
    assert (arb.pi()/16).cos()**2>arb(24)/25
    assert (2*arb.pi()/9).cos()**2<arb(3)/5
    assert (arb(18)/25).exp()>2
    modes,rows,panels,radius=168,54,4096,14
    plain=sharp_moments(modes)
    amended=(list(plain[0]),plain[1])
    amended[0][1]+=24/(25*arb.pi()**2)
    quadrature=Quadrature(panels,radius,2*modes)
    start,stop=Fraction(513,800),Fraction(2589,4000)
    width=Fraction(1,4000)
    if point is not None:
        assert start<=point<=stop
        bounds=[(point,point)]
    else:
        bounds=[(start+i*width,start+(i+1)*width) for i in range(int((stop-start)/width))]

    observations=[]
    for left,right in bounds:
        lo,hi=enclose(left),enclose(right)
        e,residual,delta,*trial=block(lo,hi,plain,amended,rows,quadrature)
        shape=geometry(lo,hi,*trial)
        median=median_geometry(lo,hi,*trial)
        passed,gap,location=residual_gate(e,residual,delta,lo,hi,shape,median=median)
        assert passed,(str(left),str(right),gap,location)
        observations.append(dict(rho=[str(left),str(right)],e_upper=upper_float(e),
                                 residual_upper=upper_float(residual),delta_lower=-upper_float(-delta),
                                 gap_upper=gap,residual_at_max=location))

    folder=Path(__file__).resolve().parent
    sources=(Path(__file__),folder/"check_three_cell_excess_purity_extension.py",
             folder/"check_three_cell_sharp_undamped_extension.py",
             folder/"check_three_cell_global_scalar_through_half.py",
             folder/"check_three_cell_sharp_positive_arc_envelope.py")
    return dict(status="pass",rho_bounds=[str(bounds[0][0]),str(bounds[-1][1])],
                precision_bits=ctx.prec,mode_order=modes,row_order=rows,radial_limit=radius,
                radial_panels=panels,correlation_slabs=len(bounds),slab_width=str(width if point is None else 0),
                max_gap_upper=max(q["gap_upper"] for q in observations),
                source_sha256={p.name:sha256(p.read_bytes()).hexdigest() for p in sources},
                slabs=observations)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--point",type=Fraction)
    args=parser.parse_args()
    print(json.dumps(certificate(args.point),indent=2))
