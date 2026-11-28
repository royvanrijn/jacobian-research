"""Directed signed-axis feasibility comparison through 2641/4000.

The proof and replay command are in THREE_CELL_SIGNED_AXIS_EXTENSION.md.
The anchored spectral block and the original one-sided reserve are
reused from the full-radius trial checker.
"""

import argparse
from bisect import bisect_left
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from flint import arb
import check_three_cell_full_radius_trial_extension as base


def signed_tables(lo,hi,row,coefficients,center,norm,kappa):
    d_u=1-hi**2
    beta=1+kappa
    norm_mu=arb.pi().sqrt()*d_u*d_u.sqrt()/2
    panel=arb(4)/2400
    clip=[[],[]]
    moment=[[],[]]
    for i in range(2400):
        x=(i*panel).union((i+1)*panel)
        upper=((i+1)*panel).upper()
        q=base.polynomial(row,x)
        xf=4*d_u/(norm.sqrt()*q)*(base.polynomial(coefficients,x)-center*(kappa*x**2).exp())
        if xf>0:
            sign=0
        elif xf<0:
            sign=1
        else:
            continue
        xi=abs(xf).lower()
        w=(q*(-beta*x**2).exp()*panel/(2*d_u)).lower()
        mu=(2*x**2*(-beta*x**2).exp()*panel/norm_mu).lower()
        phi=xi/upper
        key=-base.upper_float(-phi)
        clip[sign].append((key,w*xi**2,w*xi*upper,w*upper**2))
        moment[sign].append((key,mu*phi**2,mu))

    def cumulative(table):
        table.sort(key=lambda v:v[0])
        sums=[tuple(arb(0) for _ in table[0][1:])]
        for row in reversed(table):
            sums.append(tuple(x+y for x,y in zip(sums[-1],row[1:])))
        return [v[0] for v in table],list(reversed(sums))

    return [cumulative(t) for t in clip],[cumulative(t) for t in moment],((1-lo**2)/d_u).sqrt().upper()


def directions(panels):
    pi=arb.pi()
    root2=arb(2).sqrt()
    root6=arb(6).sqrt()
    for i in range(panels):
        l=pi*i/(3*panels)
        u=pi*(i+1)/(3*panels)
        a=(l.cos()/root6+l.sin()/root2).lower()
        b=max(arb(0),(3*u.cos()/root6-u.sin()/root2).lower())
        ap=(2*u.cos()/root6).lower()
        bp=max(arb(0),(root2*l.sin()).lower())
        median=max(abs(-l.cos()/root6+l.sin()/root2).upper(),abs(-u.cos()/root6+u.sin()/root2).upper())
        yield (a,b,ap,bp),median


def clip_bound(A,params,table):
    if A==0:
        return arb(0),arb(0)
    clips,_,scale=table
    w=scale/A
    f1,f2=arb(0),arb(0)
    for (a,b),(keys,records) in zip([params[:2],params[2:]],clips):
        for k,c,factor in [(3*a,arb(1)/6,a/2),(b,arb(1)/2,b/2)]:
            if k<=0:
                continue
            i=bisect_left(keys,base.upper_float(w/k))
            x0,x1,x2=records[i]
            f1+=c*(k**2*x0-2*k*w*x1+w**2*x2)
            f2+=factor*(k*x0-w*x1)
    return max(arb(0),f1.lower()),max(arb(0),f2.lower())


def raw_clip(A,params,table,unclipped):
    _,moments,_=table
    reduction=arb(0)
    for (a,b),(keys,records) in zip([params[:2],params[2:]],moments):
        for k,c in [(3*a,arb(1)/6),(b,arb(1)/2)]:
            if k<=0:
                continue
            i=bisect_left(keys,base.upper_float(1/(k*A)))
            x0,x1=records[i]
            reduction+=c*((k*A)**2*x0-x1)
    return (A**2*unclipped-reduction).upper()


def cap_amplitude(r,upper,params,table):
    def comparison(A):
        clip,forced=clip_bound(A,params,table)
        return clip+forced**2>r**2
    if not comparison(upper):
        return upper
    lower=arb(0)
    for _ in range(15):
        mid=(lower+upper)/2
        if comparison(mid):
            upper=mid
        else:
            lower=mid
    return upper


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
                    anchor=base.central_lower(upper,axial,rem,r1,hi,central)
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
    start,stop=Fraction(2611,4000),Fraction(2641,4000)
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
        central=base.central_geometry(lo,hi,*trial)
        full_clip=base.clip_geometry(lo,hi,*trial)
        directed=signed_tables(lo,hi,*trial)
        gap,worst=test(e,residual,delta,lo,hi,shape,median,central,directed,full_clip)
        assert gap<0,(str(left),str(right),gap,worst)
        observations.append(dict(rho=[str(left),str(right)],gap_upper=gap,worst_bin=worst))

    folder=Path(__file__).resolve().parent
    sources=(Path(__file__),folder/"check_three_cell_full_radius_trial_extension.py",
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
