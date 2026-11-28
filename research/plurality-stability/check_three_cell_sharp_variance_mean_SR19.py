"""Reserved sharp scalar comparison using the tilted mean and the SR19 bound."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

from flint import arb,ctx
import check_three_cell_sharp_undamped_extension as undamped
import check_three_cell_sharp_variance_bootstrap_anchor as bootstrap
import sharper_hermite_radial_quadrature as integration


def moment_bound(d,kappa,beta,row,tilt,coeff,center,norm,slope):
    """Bound the raw tilted second moment divided by the trial amplitude squared."""
    density=(4*d*tilt/row[0]).upper()
    width=arb(2)/1024
    lost_mass=arb(0)
    for i in range(1024):
        x=(i*width).union((i+1)*width)
        denominator=bootstrap.poly(row,x)
        assert denominator>0
        numerator=bootstrap.poly(coeff,x)-center*(kappa*x**2).exp()
        mass=8*d*(-beta*x**2).exp()*numerator**2/(norm*denominator)
        difference=max(arb(0),(density-4*d*tilt/denominator).lower())
        lost_mass += difference*mass*width
    trial_moment=(density-lost_mass.lower()).upper()
    assert trial_moment>0
    if trial_moment<density/(1+slope**2):
        return (trial_moment.sqrt()+slope*(density-trial_moment).sqrt())**2
    return density*(1+slope**2)


def orthogonal_clip(d,kappa,beta,row,coeff,center,norm,amplitude):
    width=arb(1)/768
    ratio=((arb(2)/3).sqrt()/amplitude).upper()
    radial=arb(0)
    normal=arb(0)
    for i in range(768):
        x=(i*width).union((i+1)*width)
        denominator=bootstrap.poly(row,x)
        assert denominator>0
        numerator=bootstrap.poly(coeff,x)-center*(kappa*x**2).exp()
        value=max(arb(0),abs(4*d*numerator/(denominator*norm.sqrt())).lower())
        gap=max(arb(0),(value-ratio*x).lower())
        weight=(denominator*(-beta*x**2).exp()*width/(2*d)).lower()
        radial+=weight*gap**2
        normal+=weight*value*gap
    return (radial+normal**2).lower()


def certify_clip(rho,angular,rows,quadrature,current,raw_bound,new_gain):
    d,kappa,beta,tilt,row,modes,tail=bootstrap.parts(rho,angular,rows,current)
    tt,uu,ww=quadrature.radial_moments(d,kappa,beta,row)
    coeff,center,norm,_,_,_,slope=bootstrap.trial(modes,tt,uu,ww,tail)
    bound=moment_bound(d,kappa,beta,row,tilt,coeff,center,norm,slope)
    raw=undamped.enclose(raw_bound)
    amplitude=(raw/bound).sqrt().lower()
    assert amplitude>0
    clipped=orthogonal_clip(d,kappa,beta,row,coeff,center,norm,amplitude)
    if not clipped>slope**2:
        return None

    raw_mean=sum((a*(tt[i+1]-uu[i]*uu[1]/ww) for i,a in enumerate(coeff)),arb(0))
    mean_trial=tilt*raw_mean/norm.sqrt()
    mean_norm=tilt**2*(tt[2]-uu[1]**2/ww)
    complementary=(mean_norm-mean_trial**2).upper().sqrt()
    mean_direct=(amplitude*(abs(mean_trial)+slope*complementary)).upper()
    assert kappa<1
    chi=arb.pi()*d*d.sqrt()/(2*(1-kappa).sqrt())-1
    mean_chi=(chi*raw/(1+chi)).sqrt().upper()
    mean_upper=min(mean_direct,mean_chi)
    assert mean_upper<(arb(2)/3).sqrt()/2
    gain=arb(4)/3-2*(arb(2)/3).sqrt()*mean_upper+2*mean_upper**2-3*raw
    if not gain>undamped.enclose(new_gain):
        return None
    return dict(current=str(current),max_second_moment=str(raw_bound),new_gain=str(new_gain),
                moment_gain_lower=str(gain.lower()),mean_upper=undamped.upper_float(mean_upper),
                h_squared_upper=undamped.upper_float(slope**2),
                amplitude_cap_lower=str(amplitude),clipped_squared_lower=str(clipped),
                clipping_margin_lower=str((clipped-slope**2).lower()))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start",type=Q,default=Q(641,1000))
    parser.add_argument("--stop",type=Q,default=Q(129,200))
    parser.add_argument("--width",type=Q,default=Q(1,4000))
    parser.add_argument("--panels",type=int,default=512)
    parser.add_argument("--modes",type=int,default=204)
    parser.add_argument("--rows",type=int,default=104)
    parser.add_argument("--cap",type=Q,default=Q(1999,1000))
    args=parser.parse_args()
    ctx.prec=110
    assert 0<args.cap<Q(199999,100000)
    assert 0<args.rows<=args.modes
    assert args.start>=Q(29,50)
    assert undamped.enclose(args.start)**2/(1-undamped.enclose(args.start)**2)>arb(1)/2
    assert (arb.pi()/16).cos()**2>arb(24)/25
    assert (2*arb.pi()/9).cos()**2<arb(3)/5
    assert (arb(18)/25).exp()>2
    angular=undamped.sharp_moments(args.modes)
    quadrature=integration.HermiteQuadrature(args.panels,12,2*args.modes)
    steps=((Q(12,25),Q(1,6),Q(13,20)), (Q(13,20),Q(7,200),Q(11,10)))
    count=(args.stop-args.start)/args.width
    assert count.denominator==1 and count>0
    slabs=[]
    def certify(left,right,depth=0):
        rho=undamped.enclose(left).union(undamped.enclose(right))
        stages=[]
        for current,raw_bound,new_gain in steps:
            stage=certify_clip(rho,angular,args.rows,quadrature,current,raw_bound,new_gain)
            if stage is None:
                break
            stages.append(stage)
        if len(stages)==len(steps):
            d,kappa,beta,_,row,modes,tail=bootstrap.parts(rho,angular,args.rows,steps[-1][2])
            _,finite=undamped.gram_bounds(modes,*quadrature.radial_moments(d,kappa,beta,row))
            final=finite+tail
            if final<undamped.enclose(args.cap):
                slabs.append(dict(rho=[str(left),str(right)],stages=stages,
                                  finite_upper=undamped.upper_float(finite),
                                  tail_upper=undamped.upper_float(tail),
                                  total_upper=undamped.upper_float(final)))
                return
        assert depth<8,(str(left),str(right),len(stages))
        midpoint=(left+right)/2
        certify(left,midpoint,depth+1)
        certify(midpoint,right,depth+1)

    for index in range(count.numerator):
        certify(args.start+index*args.width,args.start+(index+1)*args.width)
    paths=(Path(__file__),Path(undamped.__file__),Path(bootstrap.__file__),Path(integration.__file__),
           undamped.folder/"check_three_cell_global_scalar_through_half.py",
           undamped.folder/"check_three_cell_sharp_positive_arc_envelope.py")
    print(json.dumps(dict(status="pass",rho_bounds=[str(args.start),str(args.stop)],cap=str(args.cap),
                          steps=[list(map(str,s)) for s in steps],preliminary_radial_reserve="1/100000",
                          mode_order=args.modes,row_order=args.rows,panels=args.panels,radial_limit=12,
                          clip_panels=768,clip_limit=1,tilted_moment_panels=1024,tilted_moment_limit=2,
                          precision_bits=ctx.prec,
                          total_upper=max(s["total_upper"] for s in slabs),slabs=slabs,
                          source_sha256={p.name:sha256(p.read_bytes()).hexdigest() for p in paths}),indent=2))


if __name__=="__main__":
    main()
