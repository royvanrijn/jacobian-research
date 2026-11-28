"""Directed signed-axis and shifted clipping bootstrap for the sharp arc bound."""

from bisect import bisect_left
from fractions import Fraction as Q
import argparse
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import inspect
import json
from pathlib import Path

from flint import arb, ctx

import check_q_directmean_bootstrap_spent as spent

enclose=spent.enclose
upper_float=spent.upper_float
poly=spent.poly


def suffix(records):
    records.sort(key=lambda v:v[0])
    out=[tuple(arb(0) for _ in records[0][1:])]
    for rec in reversed(records):
        out.append(tuple(a+b for a,b in zip(rec[1:],out[-1])))
    return [r[0] for r in records],list(reversed(out))


def boxes(count):
    pi=arb.pi()
    sqrt6=arb(6).sqrt()
    sqrt2=arb(2).sqrt()
    for i in range(count):
        left=i*pi/(3*count)
        right=(i+1)*pi/(3*count)
        yield (
            (left.cos()/sqrt6+left.sin()/sqrt2).lower(),
            max(arb(0),(3*right.cos()/sqrt6-right.sin()/sqrt2).lower()),
            (2*right.cos()/sqrt6).lower(),
            max(arb(0),(sqrt2*left.sin()).lower()),
        )


def clip_values(amplitude,box,tables):
    if amplitude==0:
        return arb(0),arb(0)
    threshold=1/amplitude
    energy,normal=arb(0),arb(0)
    for (a,b),(keys,values) in zip((box[:2],box[2:]),tables):
        for k,weight,slope in ((3*a,arb(1)/6,a/2),(b,arb(1)/2,b/2)):
            if not k>0:
                continue
            i=bisect_left(keys,(threshold/k).upper())
            x0,x1,x2=values[i]
            value=(k**2*x0-2*k*threshold*x1+threshold**2*x2).lower()
            component=(k*x0-threshold*x1).lower()
            energy+=weight*max(arb(0),value)
            normal+=slope*max(arb(0),component)
    return energy.lower(),normal.lower()


def raw_direction(amplitude,box,tables,untruncated):
    total=arb(0)
    for (a,b),(keys,values) in zip((box[:2],box[2:]),tables):
        for k,weight in ((3*a,arb(1)/6),(b,arb(1)/2)):
            if not k>0:
                continue
            i=bisect_left(keys,(1/(k*amplitude)).upper())
            moment,mass=values[i]
            total+=weight*max(arb(0),((amplitude*k)**2*moment-mass).lower())
    return (amplitude**2*untruncated-total).upper()


def clip_moment_bound(amplitude,table):
    keys,masses,seconds,tail=table
    index=bisect_left(keys,(amplitude**2).upper())
    return (masses[index]+amplitude**2*seconds[index]+tail).upper()


def local_stage(rho,d,kappa,beta,tilt,row,modes,ww,uu,tt,image,norm,roots,residual,delta,excess,slope,
                vec,clip_panels,raw_panels,angle_bins,use_shift):
    coefficients=[roots[i]*vec[i,0] for i in range(len(modes))]
    center=sum((a*b for a,b in zip(coefficients,uu)),arb(0))/ww
    scale=4*d/norm.sqrt()
    signed=[[],[]]
    moments=[[],[]]
    clip_pts=max(1200,clip_panels)
    dx=arb(4)/clip_pts
    for i in range(clip_pts):
        x=(i*dx).union((i+1)*dx)
        xmax=((i+1)*dx).upper()
        q=poly(row,x)
        assert q>0
        xf=scale*(poly(coefficients,x)-center*(kappa*x**2).exp())/q
        if xf>0:
            sign=0
        elif xf<0:
            sign=1
        else:
            continue
        xi=abs(xf).lower()
        phi=xi/xmax
        key=phi.lower()
        radial=(q*(-beta*x**2).exp()*dx/(2*d)).lower()
        mass=(2*tilt*x**2*(-beta*x**2).exp()*dx).lower()
        signed[sign].append((key,radial*xi**2,radial*xi*xmax,radial*xmax**2))
        moments[sign].append((key,mass*phi**2,mass))
    signed=[suffix(v) for v in signed]
    moments=[suffix(v) for v in moments]

    radius=arb(12)
    leading=max(range(len(row)),key=lambda i:upper_float(row[i]*radius**i))
    assert row[leading]>0
    def monomial_tail(p,eta):
        power=p-2*leading
        decay=eta-arb(max(power-1,0))/(2*radius**2)
        assert decay>0
        return radius**(power-1)*(-eta*radius**2).exp()/(2*decay*row[leading]**2)

    # Whole raw tail bound from an absolute polynomial with fixed degree.
    order=len(coefficients)-1
    magnitude=poly([abs(c) for c in coefficients],radius)
    tail=2*center**2*monomial_tail(0,1-kappa)+2*magnitude**2*monomial_tail(2*order,beta)/radius**(2*order)
    width=radius/raw_panels
    raw=arb(0)
    panels=[]
    for i in range(raw_panels):
        x=(i*width).union((i+1)*width)
        q=poly(row,x)
        assert q>0
        raw_panel=(32*tilt*d*d/norm)*(-(1-kappa)*x*x).exp()*((-kappa*x*x).exp()*poly(coefficients,x)-center)**2/q**2*width
        raw+=raw_panel
        if i<raw_panels//3:
            sat=(arb(4)/3*tilt*x*x*(-beta*x*x).exp()*width).upper()
            upper=raw_panel.upper()
            assert upper>0
            key=(sat/upper).upper()
            panels.append((key,sat,upper))
    M=(raw+32*tilt*d*d/norm*arb(0).union(tail)).upper()
    x0=(raw_panels//3)*width*beta.sqrt()
    clip_tail=(arb(2)/3)*(x0.erfc()+2*x0*(-x0**2).exp()/arb.pi().sqrt())
    panels.sort(key=lambda v:v[0])
    masses=[arb(0)]; seconds=[arb(0)]
    for _,m,s in panels:
        masses.append(masses[-1]+m)
    for _,m,s in reversed(panels):
        seconds.append(seconds[-1]+s)
    moment_table=[p[0] for p in panels],masses,list(reversed(seconds)),clip_tail

    B=(4*d*tilt/row[0]).upper()
    h=slope.upper()
    if not M<B/(1+h*h):
        return None
    raw_coefficient=(M.sqrt()+h*(B-M).sqrt())**2
    raw_mean=tilt*image[1,0]/(roots[1]*norm.sqrt())
    mean_norm=tilt**2*(tt[2]-uu[1]**2/ww)
    mean_coefficient=abs(raw_mean)+h*(mean_norm-raw_mean**2).upper().sqrt()
    chi=arb.pi()*d*d.sqrt()/(2*(1-kappa).sqrt())-1
    assert chi>0
    simplex=(arb(2)/3).sqrt()
    result=arb(4)/3
    raw_profile=arb(0)
    mean_profile=arb(0)
    max_amplitude=Q(0)
    shifts=[Q(1,25),Q(2,25),Q(3,25),Q(4,25),Q(1,5)] if use_shift else []

    for direction in boxes(angle_bins):
        def clip(amplitude):
            one,two=clip_values(amplitude,direction,signed)
            out=(one+two**2).lower()
            if shifts:
                for rational in shifts:
                    shift=enclose(rational)
                    shifted,_=clip_values((1+shift)*amplitude,direction,signed)
                    out=max(out,((1+shift)**2*shifted-shift**2).lower())
            return out

        left,right=Q(0),Q(2,5)
        if not clip(enclose(right))>h*h:
            return None
        while right-left>Q(1,15000):
            mid=(left+right)/2
            if clip(enclose(mid))>h*h:
                right=mid
            else:
                left=mid
        max_amplitude=max(max_amplitude,right)
        lower=Q(0)
        for j in range(1,33):
            upper=right*j/32
            A=enclose(upper)
            spent,aligned=clip_values(enclose(lower),direction,signed)
            available=(h*h-spent).upper()
            free=available-aligned**2
            if free<0:
                break
            available=max(arb(0),available)
            clipped=min(clip_moment_bound(A,moment_table),raw_direction(A,direction,moments,M))
            assert clipped>0
            allowance=B.sqrt()*available.sqrt()
            if aligned**2>available*M/B and free>=0:
                allowance=(M.sqrt()*aligned+(B-M).sqrt()*free.sqrt()).upper()
            candidate=(clipped.sqrt()+A*allowance)**2
            upper_raw=min((A**2*raw_coefficient).upper(),candidate.upper())
            mu_bound=min((A*mean_coefficient).upper(),(chi*upper_raw/(1+chi)).sqrt().upper())
            if not mu_bound<simplex/2:
                return None
            result=min(result,(arb(4)/3-2*simplex*mu_bound+2*mu_bound**2-3*upper_raw).lower())
            raw_profile=max(raw_profile,upper_raw)
            mean_profile=max(mean_profile,mu_bound)
            lower=upper

    return dict(status="improved",gain_lower=result,squared_ratio_upper=upper_float(h*h),
                amplitude_bound=str(max_amplitude),mu_raw_profile_upper=upper_float(raw_profile),
                mu_mean_upper=upper_float(mean_profile))


def build_probe(angle_bins,use_shift):
    prefix=inspect.getsource(spent.gain_bound).split("    radius=(arb(2)/3).sqrt()")[0]
    code=prefix+("    return local_stage(rho,d,kappa,beta,tilt,row,modes,ww,uu,tt,image,norm,roots,"
                 f"residual,delta,excess,slope,vec,clip_panels,raw_panels,{angle_bins},{use_shift!r})\n")
    env=dict(vars(spent),local_stage=local_stage)
    exec(compile(code,"<signed-bootstrap-scratch>","exec"),env)
    return env["gain_bound"]


def initialize():
    global checker, quad, angular
    ctx.prec=110
    checker=build_probe(16,True)
    quad=spent.HermiteQuadrature(512,12,384)
    angular=spent.sharp_moments(192)


def cover(interval,depth=0):
    left,right=interval
    rho=enclose(left).union(enclose(right))
    assert rho*rho/(1-rho*rho)>arb(1)/2 and rho*rho/(1-rho*rho)<1
    gain=Q(12,25)
    stages=[]
    for _ in range(24):
        stage=checker(rho,gain,quad,angular,96,1200,4096,arb(1999)/1000)
        if stage is None:
            break
        if stage["status"]=="negative":
            return [dict(rho=[str(left),str(right)],stages=stages,coefficient=str(gain),
                         total_upper=stage["total_upper"])]
        lower=stage.pop("gain_lower")
        updated=Q(int(float((1000*lower).floor().lower())),1000)
        if updated<=gain:
            break
        assert enclose(updated)<lower
        stages.append(dict(gain=str(gain),updated=str(updated),gain_lower=str(lower),**stage))
        gain=updated
    assert depth<6,(str(left),str(right),stages)
    mid=(left+right)/2
    return cover((left,mid),depth+1)+cover((mid,right),depth+1)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start",type=Q,default=Q(329,500))
    parser.add_argument("--stop",type=Q,default=Q(671,1000))
    parser.add_argument("--width",type=Q,default=Q(1,16000))
    parser.add_argument("--workers",type=int,default=1)
    args=parser.parse_args()
    slabs=(args.stop-args.start)/args.width
    assert slabs.denominator==1 and slabs>0 and args.workers>0
    intervals=[(args.start+i*args.width,args.start+(i+1)*args.width) for i in range(slabs.numerator)]
    if args.workers==1:
        initialize()
        accepted=[leaf for interval in intervals for leaf in cover(interval)]
    else:
        with ProcessPoolExecutor(max_workers=args.workers,initializer=initialize) as pool:
            accepted=[leaf for result in pool.map(cover,intervals) for leaf in result]
    import check_three_cell_sharp_undamped_extension as base
    import check_three_cell_sharp_variance_bootstrap_anchor as anchor
    sources=[Path(__file__),Path(spent.__file__),Path(base.__file__),Path(anchor.__file__),
             base.folder/"sharper_hermite_radial_quadrature.py",
             base.folder/"check_three_cell_global_scalar_through_half.py",
             base.folder/"check_three_cell_sharp_positive_arc_envelope.py"]
    print(json.dumps(dict(status="pass",rho_bounds=[str(args.start),str(args.stop)],
                          slabs=accepted,total_upper=max(x["total_upper"] for x in accepted),
                          directions=16,clip_panels=1200,clip_limit=4,amplitude_bins=32,
                          raw_panels=4096,raw_limit=12,mode_order=192,row_order=96,panels=512,
                          precision_bits=110,initial_coefficient="12/25",cap="1999/1000",
                          preliminary_radial_reserve="1/100000",
                          source_sha256={path.name:sha256(path.read_bytes()).hexdigest() for path in sources}),indent=2))


if __name__=="__main__":
    main()
