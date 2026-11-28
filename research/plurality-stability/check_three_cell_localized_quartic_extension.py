"""Directed localized-quartic comparison through correlation 2681/4000.

See THREE_CELL_LOCALIZED_QUARTIC_EXTENSION.md. The spectral block,
signed-axis feasibility, and absorbed degree-one variance are unchanged.
"""
import argparse
from bisect import bisect_left
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from flint import arb
import check_three_cell_full_radius_trial_extension as base
from check_three_cell_signed_axis_extension import (
    cap_amplitude, clip_bound, directions, raw_clip, signed_tables,
)
from check_three_cell_graded_central_extension import richer_central_geometry

upper_float=base.upper_float


def tilt_slopes(lo,hi,row,coefficients,center,norm,kappa,variances,densities):
    d_u=1-hi**2
    beta=1+kappa
    pi=arb.pi()
    radius=arb(12)
    grid=3600
    dx=radius/grid
    rates=[beta,beta,((1+lo)*beta).union((1+hi)*beta)]
    degrees=[1,3,1]
    gamma0=pi.sqrt()/2
    gammas=[gamma0,3*gamma0/2,gamma0]
    normalizations=[2*rate**(1+n//2)*rate.sqrt()/gamma for rate,n,gamma in zip(rates,degrees,gammas)]
    cells=[]
    means=[arb(0),arb(0),arb(0)]
    for j in range(grid):
        x=(j*dx).union((j+1)*dx)
        square=x**2
        q=base.polynomial(row,x)
        xphi=4*d_u/(norm.sqrt()*q)*(base.polynomial(coefficients,x)-center*(kappa*square).exp())
        pairs=[]
        for i,(rate,n,c) in enumerate(zip(rates,degrees,normalizations)):
            d=c*x**(n-1)*(-rate*square).exp()*dx
            weight=2*d_u*c*x**(n-1)*(-(rate-beta)*square).exp()/q
            means[i]+=d*x*xphi
            pairs.append((d,weight))
        cells.append((x,xphi,pairs))

    upper_poly=base.polynomial([abs(a) for a in coefficients],radius)
    order=len(coefficients)-1
    denominator=base.polynomial(row,radius)
    values=[]
    for i,(rate,n,c) in enumerate(zip(rates,degrees,normalizations)):
        pow1=2*order+n-1
        pow2=n-1
        decay1=rate-arb(max(0,pow1-1))/(2*radius**2)
        decay2=rate-2*kappa-arb(max(0,pow2-1))/(2*radius**2)
        decaym=rate-arb(n)/(2*radius**2)
        assert decay1>0 and decay2>0 and decaym>0
        tail_second=16*d_u**2*c/(norm*denominator**2)*radius**(n-2)*(upper_poly**2*(-rate*radius**2).exp()/decay1+center**2*(-(rate-2*kappa)*radius**2).exp()/decay2)
        tail_mass=c*radius**n*(-rate*radius**2).exp()/(2*decaym)
        tail_mean=(tail_second*tail_mass).sqrt().upper()
        means[i]+=(-tail_mean).union(tail_mean)
        values.append(2*densities[i]*(tail_second+abs(means[i])**2*tail_mass))
    for x,xphi,pairs in cells:
        for i,(d,weight) in enumerate(pairs):
            values[i]+=d*weight*(xphi-means[i]*x)**2

    slopes=[]
    for value,V,B in zip(values,variances,densities):
        direct=max(arb(0),(value-V**2).upper()).sqrt().upper()/V.sqrt().lower()
        slopes.append(min((B-V).sqrt().upper(),direct.upper()))
    return slopes


def local_quartic_tables(lo,hi,row,coefficients,center,norm,kappa):
    d_u=1-hi**2
    factor=d_u/(1-lo**2)
    beta=1+kappa
    normalization=arb.pi().sqrt()*d_u*d_u.sqrt()/2
    step=arb(4)/2400
    data=[]
    for i in range(1,2400):
        a=i*step
        x=a.union((i+1)*step)
        q=base.polynomial(row,x)
        xphi=4*d_u/(norm.sqrt()*q)*(base.polynomial(coefficients,x)-center*(kappa*x**2).exp())
        density=2/normalization*(-beta*x**2).exp()*step
        data.append((base.upper_float(abs(xphi)/a),(density*x**2).lower(),
                     (factor*density*xphi**2).lower(),density*x*xphi,
                     (8/(arb.pi().sqrt()*d_u.sqrt()*base.polynomial(row,a))).upper()))
    data.sort(key=lambda x:x[0])
    m=[arb(0)];s=[arb(0)];u=[arb(0)];Bs=[arb(0)]
    for _,m0,s0,u0,b0 in data:
        m.append((m[-1]+m0).lower())
        s.append((s[-1]+s0).lower())
        u.append(u[-1]+u0)
        Bs.append(max(Bs[-1],b0))
    return [d[0] for d in data],m,s,u,Bs


def local_quartic(A,axis,rem,r,raw,variance,coefficient,table):
    tags,masses,seconds,means,ratios=table
    rstar=(arb(2)/3).sqrt()
    best=arb(0)
    for tau in [arb(1)/20,arb(3)/40,arb(1)/10,arb(1)/8,arb(3)/20]:
        one=1-(rstar+tau)**2/2
        both=1-2*tau**2
        for j in (2,3,4,5,6):
            q=arb(j)/10
            i=bisect_left(tags,-base.upper_float(-q*tau/(A*(1+axis))))
            B=ratios[i]
            bad=B*(A*rem/(tau*(1-q)))**2
            good=(masses[i]-bad).lower()
            if not 2*good*one>arb(12)/25:
                continue
            base_gain=(2*good*one-arb(12)/25)*variance
            # Restrict the trial variance to selected good cells. Its conditional
            # mean is bounded by the whole raw moment and the surviving mass.
            conditional=(raw/good).sqrt()
            trial_loss=B*rem**2/(1-q)**2*(q/(1+axis)+A*conditional/tau)**2
            s=seconds[i]-abs(means[i])**2/masses[i]-trial_loss
            extra=arb(0)
            if s>0:
                spread=s.sqrt()-B.sqrt()*r
                if spread>0:
                    extra=2*good*(both-one)*spread**2
            best=max(best,(coefficient*(base_gain+extra)).lower())
    return best


def median_lower(amplitude,r,lo,hi,variance1,median,error_residual=None,slopes=None):
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
    for i,(B,V) in enumerate(zip(densities,variances)):
        deviation=(B-V).sqrt().upper()
        if slopes is not None:
            deviation=min(deviation,slopes[i])
        spread=V.sqrt()-deviation*r
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


def graded_central(amplitude,axis,orthogonal,r,hi,central,slope0=None):
    variance,density,levels,masses,bounds,seconds=central
    spread=variance.sqrt()-((density-variance).sqrt() if slope0 is None else slope0)*r
    if not spread>0:
        return arb(0)
    pi=arb.pi()
    root=(arb(2)/3).sqrt()
    f=lambda s:(1-(pi*s).cos())/pi**2
    slope=arb(3).sqrt()/(2*pi)
    origin=arb(2)/3
    prefactor=hi*pi*(1-hi)/(2*(1+hi)**2)
    best=arb(0)
    def pos_coeff(h,other):
        rn=2*(f(origin-h)-f(origin)+slope*h)/h**2
        rp=2*(f(origin+h/2)-f(origin)-slope*h/2)/(h/2)**2
        return (-3/(2*pi**2)-max(((2*rn+rp)/3).upper(),((rn+other)/2).upper())).lower()/4
    for h in [arb(7)/20,arb(3)/8,arb(2)/5,arb(5)/12]:
        coeff=pos_coeff(h,-arb(12)/25)
        both=pos_coeff(2*(h-arb(1)/3),-arb(1)/2)
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


def test(e,R,delta,lo,hi,geom,med,central,table,full_clip,local,slopes,angle_panels=24):
    c,C,X,origin_weight,V,raw=geom
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
        A_full=base.amplitude_upper(r1,c,C,X,origin_weight,True,full_clip)
        assert A_full is not None
        if raw<B/(1+r1**2):
            raw_factor=(raw.sqrt()+r1*(B-raw).sqrt())**2
        else:
            raw_factor=B*(1+r1**2)
        spread=V.sqrt()-r1*((B-V).sqrt() if slopes is None else slopes[0])
        if scalar<base.quartic_excess((A_full**2*raw_factor).upper(),coefficient,spread,circ,chi,True):
            continue
        for params,median in directions_list:
            cap=cap_amplitude(r1,A_full,params,table)
            for amplitude_index in range(12):
                left=amplitude_index*cap/12
                upper=(amplitude_index+1)*cap/12
                one,two=clip_bound(left,params,table)
                c1,c2=base.clip_components(left,c,C,X,origin_weight,full_clip)
                one,two=max(one,c1),max(two,c2)
                avail=r1**2-one-two**2
                if avail<0:
                    continue
                rem=avail.upper().sqrt()
                axial=(r1**2-one).upper().sqrt()
                clipped=raw_clip(upper,params,table,raw)
                upper_raw=min((upper**2*raw_factor).upper(),((clipped.sqrt()+upper*(axial*raw.sqrt()+B.sqrt()*rem))**2).upper())
                excess=base.quartic_excess(upper_raw,coefficient,spread,circ,chi,True)
                if local is not None and not scalar<excess and spread>0:
                    excess=max(excess,local_quartic(upper,axial,rem,r1,raw,spread**2,coefficient,local))
                if scalar<excess:
                    gap=base.upper_float(scalar-excess)
                else:
                    # Retaining the universal median bound is enough for this probe.
                    nonlinear=median_lower(upper*(1+axial),r1,lo,hi,V,med,error_residual=rem,slopes=slopes)
                    anchor=graded_central(upper,axial,rem,r1,hi,central,slope0=None if slopes is None else slopes[2])
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
    start,stop=Fraction(2661,4000),Fraction(2681,4000)
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
        table=signed_tables(lo,hi,*trial)
        local=local_quartic_tables(lo,hi,*trial)
        slopes=tilt_slopes(lo,hi,*trial,[shape[4],median[1],central[0]],[*median[0],central[1]])
        gap,worst=test(e,residual,delta,lo,hi,shape,median,central,table,full_clip,local,slopes)
        assert gap<0,(str(left),str(right),gap,worst)
        observations.append(dict(rho=[str(left),str(right)],gap_upper=gap,worst_bin=worst,
                                 slope_bounds=[upper_float(value) for value in slopes]))

    folder=Path(__file__).resolve().parent
    sources=(Path(__file__),folder/"check_three_cell_graded_central_extension.py",
             folder/"check_three_cell_signed_axis_extension.py",
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
