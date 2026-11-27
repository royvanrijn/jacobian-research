"""Directed sharp variance bootstrap via the weighted simplex-amplitude constraint."""
import argparse,json
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as Q
from math import factorial

from flint import arb,arb_mat,ctx
import numpy as np
from check_three_cell_sharp_undamped_extension import Quadrature,enclose,gram_bounds,radial_coefficients,sharp_moments,upper_float


def parts(rho,angular,row_order,boost):
    d=1-rho**2
    kappa=rho**2/d
    beta=1+kappa
    rate=2*rho/d
    row,modes=radial_coefficients(d,beta,rate,*angular,row_order)
    plain_row=row[0]
    tilt=2/(arb.pi().sqrt()*d*d.sqrt())
    gain=rho*d*enclose(boost)/(2*arb.pi())
    row[0]+=2*d*gain*tilt
    modes[1]+=gain*tilt**2
    n=len(modes)
    small_ratio=4*rate/(n+1)
    eta=1+arb(49)/25*kappa
    outer_ratio=rate/eta
    assert 0<small_ratio<1 and 0<outer_ratio<1
    small=16*(4*rate)**n/(plain_row*factorial(n)*(1-small_ratio))
    large=4096*outer_ratio**n/(55*d*eta*(1-outer_ratio))
    return d,kappa,beta,tilt,row,modes,small+large


def trial(modes,tt,uu,ww,tail):
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
    eigenvalues,vectors=np.linalg.eigh(midpoint)
    eigen_anchor=arb(round(float(eigenvalues[-1])*2**44))/2**44
    vec=arb_mat(size,1,[arb(round(float(x)*2**44))/2**44 for x in vectors[:,-1]])
    image=gram*vec
    norm=(vec.transpose()*image)[0,0]
    norm_vec=(vec.transpose()*vec)[0,0]
    assert norm>0
    direction_lower=(norm/norm_vec).lower()
    raw_residual=image-vec*eigen_anchor
    raw_residual_squared=(raw_residual.transpose()*raw_residual)[0,0]
    complementary=(hs2-direction_lower**2).upper().sqrt()
    assert direction_lower>complementary
    lambda_upper=(norm/norm_vec+raw_residual_squared/(norm_vec*(direction_lower-complementary))).upper()
    residual=(raw_residual_squared*lambda_upper/norm).upper().sqrt()
    reserve=arb(1)/100000
    delta=(2-tail-reserve-complementary).lower()
    excess=max(arb(0),(lambda_upper+tail+reserve-2).upper())
    assert delta>0
    slope=(residual+(residual**2+delta*excess).sqrt())/delta

    coefficients=[roots[i]*vec[i,0] for i in range(size)]
    center=sum((coef*u for coef,u in zip(coefficients,uu)),arb(0))/ww
    return coefficients,center,norm,excess,residual,delta,slope


def poly(coefficients,x):
    total=arb(0)
    for coefficient in reversed(coefficients):
        total=total*x+coefficient
    return total


def clip_value(d,kappa,beta,row,coefficients,center,norm,cap):
    width=arb(1)/768
    cap_scale=(arb(2)/3).sqrt()/cap
    total=arb(0)
    for i in range(768):
        x=(i*width).union((i+1)*width)
        denominator=poly(row,x)
        assert denominator>0
        normalized=4*d*(poly(coefficients,x)-center*(kappa*x**2).exp())/(norm.sqrt()*denominator)
        excess=max(arb(0),(abs(normalized)-cap_scale*x).lower())
        total+=denominator*(-beta*x**2).exp()*excess**2*width/(2*d)
    return total.lower()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--start",type=Q,default=Q(311,500))
    p.add_argument("--stop",type=Q,default=Q(253,400))
    p.add_argument("--width",type=Q,default=Q(1,2000))
    p.add_argument("--panels",type=int,default=8192)
    p.add_argument("--modes",type=int,default=168)
    p.add_argument("--rows",type=int,default=80)
    p.add_argument("--purity",type=Q,default=Q(311,500))
    p.add_argument("--cap",type=Q,default=Q(19999,10000))
    args=p.parse_args()
    ctx.prec=110
    assert 0<args.purity<Q(2,3)
    assert 0<args.rows<=args.modes
    assert args.start>=Q(29,50)
    assert enclose(args.start)**2/(1-enclose(args.start)**2)>arb(1)/2
    assert (arb.pi()/16).cos()**2>arb(24)/25
    assert (2*arb.pi()/9).cos()**2<arb(3)/5
    assert (arb(18)/25).exp()>2
    angular=sharp_moments(args.modes)
    quadrature=Quadrature(args.panels,12,2*args.modes)
    count=(args.stop-args.start)/args.width
    assert count.denominator==1 and count>0
    slabs=[]
    def certify(left,right,depth=0):
        rho=enclose(left).union(enclose(right))
        d,kappa,beta,tilt,row,modes,tail=parts(rho,angular,args.rows,Q(12,25))
        coeff,center,norm,excess,residual,delta,slope=trial(
            modes,*quadrature.radial_moments(d,kappa,beta,row),tail)
        sup_density=4*d*tilt/row[0]
        cap=((arb(2)/3-enclose(args.purity))/(sup_density*(1+slope**2))).sqrt().lower()
        assert cap>0
        clipped=clip_value(d,kappa,beta,row,coeff,center,norm,cap)
        if clipped>slope**2:
            _,_,_,_,improved_row,improved_modes,tail2=parts(rho,angular,args.rows,args.purity)
            _,improved=gram_bounds(improved_modes,*quadrature.radial_moments(d,kappa,beta,improved_row))
            combined=improved+tail2
            if combined<enclose(args.cap):
                slabs.append(dict(rho=[str(left),str(right)],h_squared_upper=upper_float(slope**2),
                    amplitude_cap_lower=str(cap),clipped_squared_lower=str(clipped),
                    clipping_margin_lower=str((clipped-slope**2).lower()),
                    finite_upper=upper_float(improved),tail_upper=upper_float(tail2),
                    total_upper=upper_float(combined)))
                return
        assert depth<8,(str(left),str(right),upper_float(slope**2),str(clipped))
        midpoint=(left+right)/2
        certify(left,midpoint,depth+1)
        certify(midpoint,right,depth+1)

    for index in range(count.numerator):
        certify(args.start+index*args.width,args.start+(index+1)*args.width)

    import check_three_cell_sharp_undamped_extension as base
    paths=(Path(__file__),Path(base.__file__),base.folder/"check_three_cell_global_scalar_through_half.py",
           base.folder/"check_three_cell_sharp_positive_arc_envelope.py")
    print(json.dumps(dict(status="pass",rho_bounds=[str(args.start),str(args.stop)],cap=str(args.cap),
        purity=str(args.purity),preliminary_constant="12/25",preliminary_radial_reserve="1/100000",mode_order=args.modes,row_order=args.rows,
        panels=args.panels,radial_limit=12,clip_panels=768,clip_limit=1,precision_bits=ctx.prec,
        total_upper=max(s["total_upper"] for s in slabs),slabs=slabs,
        source_sha256={f.name:sha256(f.read_bytes()).hexdigest() for f in paths}),indent=2))


if __name__=="__main__":
    main()
