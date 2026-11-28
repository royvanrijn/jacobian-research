"""Directed two-stage SR14/SR16 variance boost for the sharp scalar envelope."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

from flint import arb,ctx
from sharper_hermite_radial_quadrature import HermiteQuadrature
import check_three_cell_sharp_undamped_extension as undamped
import check_three_cell_sharp_variance_bootstrap_anchor as bootstrap


def check_stage(rho,angular,rows,quadrature,current,raw_bound,new_gain):
    d,kappa,beta,tilt,row,modes,tail=bootstrap.parts(rho,angular,rows,current)
    coeff,center,norm,_,_,_,slope=bootstrap.trial(
        modes,*quadrature.radial_moments(d,kappa,beta,row),tail)
    moment=undamped.enclose(raw_bound)
    density=4*d*tilt/row[0]
    cap=(moment/(density*(1+slope**2))).sqrt().lower()
    clipped=bootstrap.clip_value(d,kappa,beta,row,coeff,center,norm,cap)
    if clipped<=slope**2:
        return None
    assert kappa<1
    chi=arb.pi()*d*d.sqrt()/(2*(1-kappa).sqrt())-1
    mean_upper=(chi*moment/(1+chi)).sqrt().upper()
    radius=(arb(2)/3).sqrt()
    assert mean_upper<radius/2
    improved=arb(4)/3-2*radius*mean_upper+2*mean_upper**2-3*moment
    assert improved>undamped.enclose(new_gain)
    return dict(current=str(current),max_second_moment=str(raw_bound),new_gain=str(new_gain),
                moment_gain_lower=str(improved.lower()),
                h_squared_upper=undamped.upper_float(slope**2),
                amplitude_cap_lower=str(cap),clipped_squared_lower=str(clipped),
                clipping_margin_lower=str((clipped-slope**2).lower()))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start",type=Q,default=Q(63,100))
    parser.add_argument("--stop",type=Q,default=Q(641,1000))
    parser.add_argument("--width",type=Q,default=Q(1,1000))
    parser.add_argument("--panels",type=int,default=512)
    parser.add_argument("--modes",type=int,default=192)
    parser.add_argument("--rows",type=int,default=96)
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
    quadrature=HermiteQuadrature(args.panels,12,2*args.modes)
    steps=((Q(12,25),Q(47,300),Q(31,50)), (Q(31,50),Q(7,150),Q(1)))
    count=(args.stop-args.start)/args.width
    assert count.denominator==1 and count>0
    slabs=[]
    def certify(left,right,depth=0):
        rho=undamped.enclose(left).union(undamped.enclose(right))
        stages=[]
        for current,raw_bound,new_gain in steps:
            stage=check_stage(rho,angular,args.rows,quadrature,current,raw_bound,new_gain)
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
    import sharper_hermite_radial_quadrature as integration
    paths=(Path(__file__),Path(undamped.__file__),Path(bootstrap.__file__),Path(integration.__file__),
           undamped.folder/"check_three_cell_global_scalar_through_half.py",
           undamped.folder/"check_three_cell_sharp_positive_arc_envelope.py")
    print(json.dumps(dict(status="pass",rho_bounds=[str(args.start),str(args.stop)],cap=str(args.cap),
                          steps=[list(map(str,s)) for s in steps],preliminary_radial_reserve="1/100000",
                          mode_order=args.modes,row_order=args.rows,panels=args.panels,radial_limit=12,
                          clip_panels=768,clip_limit=1,precision_bits=ctx.prec,
                          total_upper=max(s["total_upper"] for s in slabs),slabs=slabs,
                          source_sha256={p.name:sha256(p.read_bytes()).hexdigest() for p in paths}),indent=2))


if __name__=="__main__":
    main()
