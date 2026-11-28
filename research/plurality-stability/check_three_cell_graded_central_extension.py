"""Directed graded central-reserve comparison through 2661/4000.

See THREE_CELL_GRADED_CENTRAL_EXTENSION.md. This checker reuses the
signed-axis projection and the anchored scalar operator.
"""
import argparse
from bisect import bisect_left
from fractions import Fraction
from functools import cache
from hashlib import sha256
import json
from pathlib import Path

from flint import arb
import check_three_cell_full_radius_trial_extension as base
from check_three_cell_signed_axis_extension import (
    cap_amplitude, clip_bound, directions, raw_clip, signed_tables,
)


def richer_central_geometry(lo,hi,row,coefficients,center,norm,kappa):
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
        q=base.polynomial(row,x)
        xphi=4*d_min/(norm.sqrt()*q)*(base.polynomial(coefficients,x)-center*(kappa*square).exp())
        density=(pref_lo*(-rate_hi*square).exp()).union(pref_hi*(-rate_lo*square).exp())
        mass+=density*square*width
        second+=density*xphi**2*width
        mean+=density*x*xphi*width
        if i:
            local_mass=(density*square*width).lower()
            density_bound=(ratio*(-lo*a**2/d_min).exp()/base.polynomial(row,a)).upper()
            raw_lower=(density*xphi**2*width*d_min/d_max).lower()
            data.append((base.upper_float(abs(xphi)/a),local_mass,density_bound,max(arb(0),raw_lower)))
    variance=(d_min/d_max*(mass*second-mean**2)).lower()
    density=(ratio/row[0]).upper()
    data.sort(key=lambda record:record[0])
    masses=[arb(0)]
    bounds=[arb(0)]
    seconds=[arb(0)]
    for _,w,b,s in data:
        masses.append((masses[-1]+w).lower())
        bounds.append(max(bounds[-1],b))
        seconds.append((seconds[-1]+s).lower())
    return variance,density,[record[0] for record in data],masses,bounds,seconds


@cache
def central_coefficients():
    pi=arb.pi()
    f=lambda s:(1-(pi*s).cos())/pi**2
    slope=arb(3).sqrt()/(2*pi)
    origin=arb(2)/3
    result=[]
    def pos_coeff(h,other):
        rn=2*(f(origin-h)-f(origin)+slope*h)/h**2
        rp=2*(f(origin+h/2)-f(origin)-slope*h/2)/(h/2)**2
        assert rn>other
        return (-3/(2*pi**2)-max(((2*rn+rp)/3).upper(),((rn+other)/2).upper())).lower()/4
    thresholds=zip((arb(7)/20,arb(3)/8,arb(2)/5,arb(5)/12),(411,376,342,319),(831,770,706,663))
    for h,one_lower,both_lower in thresholds:
        coeff=pos_coeff(h,-arb(12)/25)
        both=pos_coeff(2*(h-arb(1)/3),-arb(1)/2)
        assert coeff>arb(one_lower)/10000 and both>arb(both_lower)/10000
        result.append((h,coeff,both))
    return result


def graded_central(amplitude,axis,orthogonal,r,hi,central):
    variance,density,levels,masses,bounds,seconds=central
    spread=variance.sqrt()-(density-variance).sqrt()*r
    if not spread>0:
        return arb(0)
    root=(arb(2)/3).sqrt()
    prefactor=hi*arb.pi()*(1-hi)/(2*(1+hi)**2)
    best=arb(0)
    for h,coeff,both in central_coefficients():
        margin=(h-arb(1)/3)/root
        for j in range(1,10):
            stop=margin*j/(10*amplitude*(1+axis))
            index=bisect_left(levels,-base.upper_float(-stop))
            m=masses[index]-bounds[index]*(amplitude*orthogonal/(margin*(10-j)/10))**2
            good=max(arb(0),m.lower())
            trial_good=seconds[index]-bounds[index]*(orthogonal*j/((10-j)*(1+axis)))**2
            extra=arb(0)
            if trial_good>0 and both>coeff:
                good_norm=trial_good.sqrt()-bounds[index].sqrt()*r
                if good_norm>0:
                    extra=2*(both-coeff)*good*good_norm**2
            best=max(best,(coeff*2*good/(2-good)*spread**2+extra).lower())
    return (prefactor*best).lower()


def test(e,R,delta,lo,hi,geom,med,central,table,full_clip,angle_panels=24):
    pole,slope,radius,origin,V,raw=geom
    d_min=1-hi**2
    B=(8*arb.pi()/(lo*d_min*(arb(11)/2+arb(24)/25))).upper()
    coefficient=hi*d_min/(2*arb.pi())
    circ=(arb(2)/3).sqrt()
    chi=(arb.pi()*d_min**2/(2*(1-2*hi**2).sqrt())-1).upper()
    horizon=((R+(R**2+delta*e).sqrt())/delta).upper()
    directions_list=list(directions(angle_panels))
    largest=-100.0
    worst=None
    for k in range(128):
        r0=k*horizon/128
        r1=(k+1)*horizon/128
        scalar=e+2*R*r1-delta*r0**2
        if scalar<0:
            continue
        A_full=base.amplitude_upper(r1,pole,slope,radius,origin,True,full_clip)
        assert A_full is not None
        if raw<B/(1+r1**2):
            raw_factor=(raw.sqrt()+r1*(B-raw).sqrt())**2
        else:
            raw_factor=B*(1+r1**2)
        spread=V.sqrt()-r1*(B-V).sqrt()
        if scalar<base.quartic_excess((A_full**2*raw_factor).upper(),coefficient,spread,circ,chi,True):
            continue
        for params,median in directions_list:
            cap=cap_amplitude(r1,A_full,params,table)
            for j in range(12):
                left=j*cap/12
                upper=(j+1)*cap/12
                one,two=clip_bound(left,params,table)
                c1,c2=base.clip_components(left,pole,slope,radius,origin,full_clip)
                one,two=max(one,c1),max(two,c2)
                avail=r1**2-one-two**2
                if avail<0:
                    continue
                rem=avail.upper().sqrt()
                axial=(r1**2-one).upper().sqrt()
                clipped=raw_clip(upper,params,table,raw)
                upper_raw=min((upper**2*raw_factor).upper(),((clipped.sqrt()+upper*(axial*raw.sqrt()+B.sqrt()*rem))**2).upper())
                excess=base.quartic_excess(upper_raw,coefficient,spread,circ,chi,True)
                if scalar<excess:
                    gap=base.upper_float(scalar-excess)
                else:
                    # Retaining the universal median bound is enough for this probe.
                    nonlinear=base.median_lower(upper*(1+axial),r1,lo,hi,V,med,error_residual=rem)
                    anchor=graded_central(upper,axial,rem,r1,hi,central)
                    gap=base.upper_float(scalar-max(excess,nonlinear)-anchor)
                if gap>largest:
                    largest=gap;worst=[float(r1),float(upper),*[float(v) for v in params]]
    return largest,worst


def certificate(point=None):
    base.ctx.prec=120
    assert (arb.pi()/16).cos()**2>arb(24)/25
    assert (2*arb.pi()/9).cos()**2<arb(3)/5
    assert (arb(18)/25).exp()>2
    modes,rows,panels,radius=168,54,4096,14
    plain=base.sharp_moments(modes)
    amended=(list(plain[0]),plain[1])
    amended[0][1]+=24/(25*arb.pi()**2)
    quadrature=base.Quadrature(panels,radius,2*modes)
    start,stop=Fraction(2641,4000),Fraction(2661,4000)
    width=Fraction(1,4000)
    if point is not None:
        assert start<=point<=stop
        bounds=[(point,point)]
    else:
        bounds=[(start+i*width,start+(i+1)*width) for i in range(int((stop-start)/width))]
    observations=[]
    for left,right in bounds:
        lo,hi=base.enclose(left),base.enclose(right)
        e,residual,delta,*trial=base.block(lo,hi,plain,amended,rows,quadrature)
        shape=base.geometry(lo,hi,*trial)
        median=base.median_geometry(lo,hi,*trial)
        central=richer_central_geometry(lo,hi,*trial)
        full_clip=base.clip_geometry(lo,hi,*trial)
        signed=signed_tables(lo,hi,*trial)
        gap,worst=test(e,residual,delta,lo,hi,shape,median,central,signed,full_clip)
        assert gap<0,(str(left),str(right),gap,worst)
        observations.append(dict(rho=[str(left),str(right)],gap_upper=gap,worst_bin=worst))

    folder=Path(__file__).resolve().parent
    sources=(Path(__file__),folder/"check_three_cell_signed_axis_extension.py",
             folder/"check_three_cell_full_radius_trial_extension.py",
             folder/"check_three_cell_projected_trial_extension.py",
             folder/"check_three_cell_excess_purity_extension.py",
             folder/"check_three_cell_sharp_undamped_extension.py",
             folder/"check_three_cell_global_scalar_through_half.py",
             folder/"check_three_cell_sharp_positive_arc_envelope.py")
    return dict(status="pass",rho_bounds=[str(bounds[0][0]),str(bounds[-1][1])],
                precision_bits=base.ctx.prec,mode_order=modes,row_order=rows,radial_limit=radius,
                radial_panels=panels,correlation_slabs=len(bounds),slab_width=str(width if point is None else 0),
                residual_bins=128,amplitude_bins=12,direction_bins=24,
                max_gap_upper=max(q["gap_upper"] for q in observations),
                source_sha256={p.name:sha256(p.read_bytes()).hexdigest() for p in sources},
                slabs=observations)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--point",type=Fraction)
    args=parser.parse_args()
    print(json.dumps(certificate(args.point),indent=2))
