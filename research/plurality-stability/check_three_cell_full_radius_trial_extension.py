"""Directed full-radius clipped comparison through correlation 2611/4000.

See THREE_CELL_FULL_RADIUS_TRIAL_EXTENSION.md. This checker reuses the
anchored spectral block and the degree-one variance credit.
"""
import argparse
from bisect import bisect_left
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from flint import arb, ctx

from check_three_cell_excess_purity_extension import (
    Quadrature, block, enclose, geometry, polynomial, sharp_moments, upper_float,
)
from check_three_cell_projected_trial_extension import (
    median_geometry, median_lower, quartic_excess,
)


def clip_geometry(lo,hi,row,coefficients,center,norm,kappa):
    d_min=1-hi**2
    beta=1+kappa
    width=arb(4)/2400
    cells=[]
    for i in range(2400):
        right=((i+1)*width).upper()
        x=(i*width).union((i+1)*width)
        square=x**2
        q=polynomial(row,x)
        xphi=(4*d_min/norm.sqrt())*(polynomial(coefficients,x)-center*(kappa*square).exp())/q
        xi=max(arb(0),abs(xphi).lower())
        if xi>0:
            weight=(q*(-beta*square).exp()*width/(2*d_min)).lower()
            cells.append((-upper_float(-xi/right),weight*xi**2,
                          weight*xi*right,weight*right**2))
    cells.sort(key=lambda v:v[0])
    suffix=[(arb(0),arb(0),arb(0))]
    for _,c0,c1,c2 in reversed(cells):
        a0,a1,a2=suffix[-1]
        suffix.append((a0+c0,a1+c1,a2+c2))
    scale=(arb(2)/3*(1-lo**2)/d_min).sqrt().upper()
    return [v[0] for v in cells],list(reversed(suffix)),scale


def central_geometry(lo,hi,row,coefficients,center,norm,kappa):
    d_min=1-hi**2
    d_max=1-lo**2
    rate_lo=(1+lo)/d_min
    rate_hi=(1+hi)/d_min
    pi=arb.pi()
    pref_lo=4*rate_lo*rate_lo.sqrt()/pi.sqrt()
    pref_hi=4*rate_hi*rate_hi.sqrt()/pi.sqrt()
    ratio=8*(1+hi)*(1+hi).sqrt()/(pi*d_min).sqrt()
    width=arb(4)/1600
    mass,second,mean=arb(0),arb(0),arb(0)
    data=[]
    for i in range(1600):
        a=i*width
        x=a.union((i+1)*width)
        square=x**2
        q=polynomial(row,x)
        xphi=4*d_min/(norm.sqrt()*q)*(polynomial(coefficients,x)-center*(kappa*square).exp())
        density=(pref_lo*(-rate_hi*square).exp()).union(pref_hi*(-rate_lo*square).exp())
        mass+=density*square*width
        second+=density*xphi**2*width
        mean+=density*x*xphi*width
        if i:
            local_mass=(density*square*width).lower()
            density_bound=(ratio*(-lo*a**2/d_min).exp()/polynomial(row,a)).upper()
            data.append((upper_float(abs(xphi)/a),local_mass,density_bound))
    variance=(d_min/d_max*(mass*second-mean**2)).lower()
    density=(ratio/row[0]).upper()
    data.sort(key=lambda record:record[0])
    masses=[arb(0)]
    bounds=[arb(0)]
    for _,w,b in data:
        masses.append((masses[-1]+w).lower())
        bounds.append(max(bounds[-1],b))
    return variance,density,[record[0] for record in data],masses,bounds


def central_lower(amplitude,axis,orthogonal,r,hi,central):
    if central is None:
        return arb(0)
    variance,density,levels,masses,bounds=central
    spread=variance.sqrt()-(density-variance).sqrt()*r
    if not spread>0:
        return arb(0)
    central_radius=1/(12*(arb(2)/3).sqrt())
    coefficient=arb(3)/100*hi*arb.pi()*(1-hi)/(2*(1+hi)**2)
    best=arb(0)
    for k in range(1,10):
        threshold=central_radius*k/(10*amplitude*(1+axis))
        position=bisect_left(levels,-upper_float(-threshold))
        margin=central_radius*(10-k)/10
        good_mass=masses[position]-bounds[position]*(amplitude*orthogonal/margin)**2
        best=max(best,good_mass.lower())
    return (coefficient*best*spread**2).lower()


def clip_components(amplitude,pole,slope,radius,origin,table=None):
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
    first,second=first.lower(),second.lower()
    if table is not None:
        nodes,values,scale=table
        v=scale/amplitude
        c0,c1,c2=values[bisect_left(nodes,upper_float(v))]
        first=max(first,(c0-2*v*c1+v**2*c2).lower())
        second=max(second,(c0-v*c1).lower())
    return max(arb(0),first),max(arb(0),second)


def amplitude_upper(r,pole,slope,radius,origin,projected=False,table=None):
    circ=(arb(2)/3).sqrt()
    denom=2*origin*pole**3-3*slope*r**2
    cutoff=3*r**2/(2*origin*pole**2)
    if not projected and denom>0 and cutoff<=radius:
        return (3*circ*r**2/denom).upper()

    def comparison(amplitude):
        if projected:
            first,second=clip_components(amplitude,pole,slope,radius,origin,table)
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


def residual_gate(e,residual,delta,lo,hi,geometry,median,central,table):
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

    for number in range(128):
        a=number*maximum/128
        b=(number+1)*maximum/128
        turning=residual/delta
        if turning<a:
            scalar=e+2*residual*a-delta*a**2
        elif turning>b:
            scalar=e+2*residual*b-delta*b**2
        else:
            scalar=e+residual**2/delta
        if scalar<0:
            continue
        amplitude=amplitude_upper(b,pole,slope,radius,origin,projected=True,table=table)
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
                first,second_clip=clip_components(small,pole,slope,radius,origin,table)
                if first>b**2:
                    continue
                eff=(b**2-first).upper().sqrt()
                orthogonal=eff
                axis=arb(0)
                available=b**2-first-second_clip**2
                if available<0:
                    continue
                orthogonal=available.upper().sqrt()
                axis=eff
                raw_proj=(large*(raw.sqrt()+density.sqrt()*eff))**2
                if median is not None:
                    _,_,nodes,totals,prefix=median
                    index=bisect_left(nodes,float(circumradius/large))
                    clip_raw=(large**2*prefix[index]+arb(2)/3*totals[index][0]).upper()
                    raw_proj=min(raw_proj.upper(),((clip_raw.sqrt()+density.sqrt()*eff*large)**2).upper())
                    raw_proj=min(raw_proj.upper(),((clip_raw.sqrt()+large*(axis*raw.sqrt()+density.sqrt()*orthogonal))**2).upper())
                raw_upper=min(raw_proj.upper(),(large**2*second).upper())
                extra=quartic_excess(raw_upper,coefficient,spread,circumradius,chi,True)
                lower=max(best,extra.lower())
                if median is not None and not scalar<lower:
                    lower=max(lower,median_lower(large,b,lo,hi,variance,median,error_residual=eff))
                    # The amplitude in this call occurs in the guard and the quadratic
                    # cost only. Both remain upper bounds with the axial enlargement.
                    lower=max(lower,median_lower(large*(1+axis),b,lo,hi,variance,median,error_residual=orthogonal))
                if central is not None and not scalar<lower:
                    lower=max(best,(lower+central_lower(large,axis,orthogonal,r=b,hi=hi,central=central)).lower())
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
    start,stop=Fraction(2589,4000),Fraction(2611,4000)
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
        central=central_geometry(lo,hi,*trial)
        table=clip_geometry(lo,hi,*trial)
        passed,gap,location=residual_gate(e,residual,delta,lo,hi,shape,median=median,central=central,table=table)
        assert passed,(str(left),str(right),gap,location)
        observations.append(dict(rho=[str(left),str(right)],e_upper=upper_float(e),
                                 residual_upper=upper_float(residual),delta_lower=-upper_float(-delta),
                                 gap_upper=gap,residual_at_max=location))

    folder=Path(__file__).resolve().parent
    sources=(Path(__file__),folder/"check_three_cell_projected_trial_extension.py",
             folder/"check_three_cell_excess_purity_extension.py",
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
