"""Directed excess-purity trial comparison through correlation 513/800.

The proof and replay command are in THREE_CELL_EXCESS_PURITY_EXTENSION.md.
Floating-point diagonalization chooses dyadic feature coefficients only.
"""
import argparse
from fractions import Fraction
from math import factorial
import json
from pathlib import Path
from hashlib import sha256

from flint import arb, arb_mat, ctx
import numpy as np

from check_three_cell_global_scalar_through_half import radial_coefficients, upper_float
from check_three_cell_sharp_positive_arc_envelope import sharp_moments
from check_three_cell_sharp_undamped_extension import Quadrature, enclose


def polynomial(coefficients, x):
    value = arb(0)
    for coefficient in reversed(coefficients):
        value = value*x+coefficient
    return value


def tail(lo, hi, angular, size):
    rho = lo.union(hi)
    d = 1-rho**2
    kappa = rho**2/d
    rate = 2*rho/d
    row = radial_coefficients(d, 1+kappa, rate, *angular, 1)[0][0]
    short_ratio=4*rate/(size+1)
    eta=1+arb(49)/25*kappa
    long_ratio=rate/eta
    assert 0<short_ratio<1 and 0<long_ratio<1
    short=16*(4*rate)**size/(row*factorial(size)*(1-short_ratio))
    long=4096*long_ratio**size/(55*d*eta*(1-long_ratio))
    return (short+long).upper()


def imbalance_moment(d, kappa, row, radius):
    cutoff=arb(radius)
    width=cutoff/2048
    prefactor=8*d*width
    eta=1-kappa
    total=arb(0)
    for i in range(2048):
        x=(i*width).union((i+1)*width)
        total += (prefactor*x**4*(-eta*x**2).exp()/polynomial(row,x)).upper()
    decay=eta-arb(3)/(2*cutoff**2)
    assert decay>0
    return (total+4*d*cutoff**3*(-eta*cutoff**2).exp()/(row[0]*decay)).upper()


def block(lo, hi, plain, amended, row_order, quadrature):
    d=1-hi**2
    kappa=hi**2/d
    row,modes=radial_coefficients(d,1+kappa,2*hi/d,*amended,row_order)
    row=[b*(lo/hi)**n for n,b in enumerate(row,1)]
    tt,uu,ww=quadrature.radial_moments(d,kappa,1+kappa,row)
    size=len(modes)
    roots=[k.sqrt() for k in modes]
    gram=arb_mat(size,size)
    hs2=arb(0)
    for i in range(size):
        for j in range(i+1):
            entry=roots[i]*roots[j]*(tt[i+j]-uu[i]*uu[j]/ww)
            gram[i,j]=gram[j,i]=entry
            hs2+=(1 if i==j else 2)*entry**2
    midpoint=np.array([[float(gram[i,j].mid()) for j in range(size)] for i in range(size)])
    _,vectors=np.linalg.eigh(midpoint)
    r=vectors[:,-1]*np.sign(vectors[0,-1])
    rational=[arb(round(float(v)*2**46))/2**46 for v in r]
    direction=arb_mat(size,1,rational)
    image=gram*direction
    norm=(direction.transpose()*image)[0,0]
    assert norm>0
    rayleigh=(image.transpose()*image)[0,0]/norm
    residual_vector=image-direction*rayleigh.mid()
    residual=((residual_vector.transpose()*gram*residual_vector)[0,0]/norm).upper().sqrt().upper()
    complement=arb(0)
    for i in range(size):
        for j in range(i+1):
            complement+=(1 if i==j else 2)*(gram[i,j]-image[i,0]*image[j,0]/norm)**2
    complement=complement.upper().sqrt()
    flux_vector=arb_mat(size,1,[k*u/ww.sqrt() for k,u in zip(roots,uu)])
    contact=(flux_vector.transpose()*flux_vector)[0,0]
    flux=(flux_vector.transpose()*gram*flux_vector)[0,0].upper().sqrt()
    leakage=abs((hi**2-lo**2)/d*(imbalance_moment(d,kappa,row,quadrature.radius)/ww).sqrt()).upper()
    charge=(2*flux*leakage+contact*leakage**2).upper()
    truncation=tail(lo,hi,plain,size)
    e=(rayleigh-2+truncation+charge).upper()
    delta=(2-truncation-charge-complement).lower()
    coefficients=[q*k for q,k in zip(rational,roots)]
    center=sum((a*u for a,u in zip(coefficients,uu)),arb(0))/ww
    return e,residual,delta,row,coefficients,center,norm,kappa


def geometry(lo,hi,row,coefficients,center,norm,kappa):
    d_min=1-hi**2
    d_max=1-lo**2
    radius=arb(1)/4
    slope=arb(13)/10
    scale=4*d_min/norm.sqrt()
    width=radius/600
    pole=arb(2)
    for i in range(600):
        x=(i*width).union((i+1)*width)
        value=scale*(polynomial(coefficients,x)-center*(kappa*x**2).exp())/polynomial(row,x)
        pole=min(pole,(value+slope*x).lower())
    origin=(row[0]*(-radius**2/d_min).exp()/(4*(d_min*d_max).sqrt())).lower()
    beta=1+kappa
    normalization=arb.pi().sqrt()/(2*beta*beta.sqrt())
    width=arb(12)/6000
    tilt_mass,low_second,raw,mean=arb(0),arb(0),arb(0),arb(0)
    for i in range(6000):
        x=(i*width).union((i+1)*width)
        square=x**2
        den=polynomial(row,x)
        numerator=(-kappa*square).exp()*polynomial(coefficients,x)-center
        square_term=32*d_min**2/(norm*normalization)*(-(1-kappa)*square).exp()*(numerator/den)**2*width
        raw+=square_term
        if i<2000:
            tilt_mass+=2/normalization*square*(-beta*square).exp()*width
            low_second+=square_term
            mean+=8*d_min/(norm.sqrt()*normalization)*x*(-square).exp()*numerator/den*width

    # A tail for the tilted raw second moment. The largest radial monomial
    # exponent bounds the whole absolute polynomial past the cutoff.
    end=arb(12)
    order=len(coefficients)-1
    assert beta>arb(2*order-1)/(2*end**2)
    last_coeff=polynomial([abs(a) for a in coefficients],end)
    den=polynomial(row,end)
    square_tail=last_coeff**2*(-beta*end**2).exp()/(2*end*(beta-arb(2*order-1)/(2*end**2)))
    constant_tail=center**2*(-(1-kappa)*end**2).exp()/(2*(1-kappa)*end)
    raw+=arb(0).union(64*d_min**2*(square_tail+constant_tail)/(norm*normalization*den**2))
    variance=(d_min/d_max*(tilt_mass*low_second-mean**2)).lower()
    return pole,slope,radius,origin,variance,raw.upper()


def residual_gate(e,residual,delta,lo,hi,geometry,steps=128):
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
    mean_ratio=(arb.pi()*d_min**2/(2*(1-2*hi**2).sqrt())-1).upper()
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
        denom=2*origin*pole**3-3*slope*b**2
        cutoff=3*b**2/(2*origin*pole**2)
        if denom<=0 or cutoff>radius:
            return False, None, float(b)
        amplitude=3*circumradius*b**2/denom
        if raw<density/(1+b**2):
            second=(raw.sqrt()+b*(density-raw).sqrt())**2
        else:
            second=density*(1+b**2)
        raw_upper=(amplitude**2*second).upper()
        purity=arb(2)/3-raw_upper
        mean_bound=(mean_ratio*raw_upper/(1+mean_ratio)).sqrt().upper()
        quartic_factor=1-3*raw_upper
        if mean_bound<circumradius/2:
            quartic_factor=arb(4)/3-3*raw_upper-2*circumradius*mean_bound+2*mean_bound**2
        quartic_factor=max(purity.lower(),quartic_factor.lower())
        spread=variance.sqrt()-(density-variance).sqrt()*b
        extra=quartic_factor-arb(12)/25
        # The absorbed variance is already charged in the spectral block.
        excess=coefficient*extra*spread**2 if extra>0 and spread>0 else arb(0)
        gap=upper_float(scalar-excess)
        if gap>maximum_gap:
            maximum_gap=gap
            worst=float(b)
    return maximum_gap<0,maximum_gap,worst


def certificate(point=None):
    ctx.prec=120
    assert (arb.pi()/16).cos()**2>arb(24)/25
    assert (2*arb.pi()/9).cos()**2<arb(3)/5
    assert (arb(18)/25).exp()>2
    modes, rows, panels, radius = 168, 54, 4096, 14
    plain=sharp_moments(modes)
    amended=(list(plain[0]),plain[1])
    amended[0][1]+=24/(25*arb.pi()**2)
    quadrature=Quadrature(panels,radius,2*modes)
    start,stop=Fraction(3,5),Fraction(513,800)
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
        passed,gap,location=residual_gate(e,residual,delta,lo,hi,shape)
        assert passed,(str(left),str(right),gap,location)
        observations.append(dict(rho=[str(left),str(right)],e_upper=upper_float(e),
                            residual_upper=upper_float(residual),delta_lower=-upper_float(-delta),
                            pole_lower=-upper_float(-shape[0]),origin_lower=-upper_float(-shape[3]),
                            variance_lower=-upper_float(-shape[4]),second_upper=upper_float(shape[5]),
                            gap_upper=gap))
    folder=Path(__file__).resolve().parent
    sources=(Path(__file__),folder/"check_three_cell_sharp_undamped_extension.py",
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
