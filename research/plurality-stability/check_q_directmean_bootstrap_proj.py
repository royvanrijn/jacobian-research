"""Directed scratch probes of SR19 with a projected tilted mean."""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial
import argparse
from hashlib import sha256
import json
import numpy as np
from flint import arb, arb_mat, ctx
from check_three_cell_sharp_undamped_extension import enclose,gram_bounds,radial_coefficients,sharp_moments,upper_float
from sharper_hermite_radial_quadrature import HermiteQuadrature
from check_three_cell_sharp_variance_bootstrap_anchor import parts, poly, clip_value

def gain_bound(rho,boost,quad,angular,rows,clip_panels,raw_panels,cap_norm):
    d,kappa,beta,tilt,row,modes,tail=parts(rho,angular,rows,boost)
    tt,uu,ww=quad.radial_moments(d,kappa,beta,row)
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
    if lambda_upper+tail<cap_norm:
        return dict(status="negative", total_upper=upper_float(lambda_upper+tail))
    # Second spectral bound after subtracting the trial projection, now using
    # coefficient-space covariance of the projected physical-feature columns.
    compressed=gram-image*image.transpose()/norm
    hs2c=sum(((1 if a==b else 2)*compressed[a,b]**2 for a in range(size) for b in range(a+1)), arb(0))
    cmat=np.array([[float(compressed[i,j].mid()) for j in range(size)] for i in range(size)])
    ev, vmat=np.linalg.eigh(cmat)
    vv=arb_mat(size,1,[arb(round(float(x)*2**44))/2**44 for x in vmat[:,-1]])
    vi=compressed*vv
    vnorm=(vv.transpose()*vv)[0,0]
    zray=(vv.transpose()*vi)[0,0]/vnorm
    zcomp=(hs2c-zray**2).upper().sqrt()
    if zray>zcomp:
        zres=vi-vv*zray
        update=(zray+(zres.transpose()*zres)[0,0]/(vnorm*(zray-zcomp))).upper()
        complementary=min(complementary.upper(),update)
    delta=(2-tail-reserve-complementary).lower()
    excess=max(arb(0),(lambda_upper+tail+reserve-2).upper())
    assert delta>0
    slope=(residual+(residual**2+delta*excess).sqrt())/delta
    coefficients=[roots[i]*vec[i,0] for i in range(size)]
    center=sum((coef*u for coef,u in zip(coefficients,uu)),arb(0))/ww
    radius=(arb(2)/3).sqrt()
    width=arb(1)/clip_panels
    # Cache the pole and mass bounds. All neglected x>1 terms are nonnegative.
    lower=[]
    for i in range(clip_panels):
        x=(i*width).union((i+1)*width)
        Qx=poly(row,x)
        assert Qx>0
        xphi=4*d*(poly(coefficients,x)-center*(kappa*x**2).exp())/(norm.sqrt()*Qx)
        weight=Qx*(-beta*x*x).exp()*width/(2*d)
        lower.append(((radius*x).upper(),abs(xphi).lower(),weight.lower()))
    def clip(cap):
        radial=arb(0)
        normal=arb(0)
        for x,phi,weight in lower:
            u=max(arb(0),(phi-x/cap).lower())
            radial+=weight*u*u
            normal+=weight*phi*u
        return (radial+normal**2).lower()
    left,right=Q(1,1000),Q(1,3)
    if clip(enclose(right))<=slope**2:return None
    while right-left>Q(1,15000):
        mid=(left+right)/2
        if clip(enclose(mid))>slope**2:
            right=mid
        else:
            left=mid
    cap=enclose(right)
    clipping=clip(cap)
    assert clipping>slope**2
    # Explicit tail bound for the tilted raw moment.
    radius=arb(12)
    leading=max(range(len(row)),key=lambda j:upper_float(row[j]*radius**j))
    assert row[leading]>0
    def int_monomial_tail(p,eta):
        power=p-2*leading
        decay=eta-arb(max(power-1,0))/(2*radius**2)
        assert decay>0
        return radius**(power-1)*(-eta*radius**2).exp()/(2*decay*row[leading]**2)
    tail2=2*center**2*int_monomial_tail(0,1-kappa)
    for i,ci in enumerate(coefficients):
        for j in range(i+1):
            cj=coefficients[j]
            tail2+=2*(1 if i==j else 2)*abs(ci*cj)*int_monomial_tail(i+j,beta)
    width=arb(12)/raw_panels
    integral=arb(0)
    for i in range(raw_panels):
        x=(i*width).union((i+1)*width)
        Qx=poly(row,x)
        if not Qx>0:
            raise ValueError('bad row')
        integrand=(-kappa*x*x).exp()*poly(coefficients,x)-center
        integral+=((-(1-kappa)*x*x).exp()*integrand**2/Qx**2)*width
    M=(32*tilt*d*d*(integral+arb(0).union(tail2))/norm).upper()
    B=(4*d*tilt/row[0]).upper()
    h=slope.upper()
    if not M<B/(1+h*h):return None
    projected=cap**2*(M.sqrt()+h*(B-M).sqrt())**2
    raw_mean=tilt*image[1,0]/(roots[1]*norm.sqrt())
    centered_mean_norm=tilt**2*(tt[2]-uu[1]**2/ww)
    mu_bound=(cap*(abs(raw_mean)+h*(centered_mean_norm-raw_mean**2).upper().sqrt())).upper()
    chi=arb.pi()*d*d.sqrt()/(2*(1-kappa).sqrt())-1
    alternative=(chi*projected/(1+chi)).sqrt().upper()
    m=min(mu_bound,alternative)
    simplex=(arb(2)/3).sqrt()
    if not m<simplex/2:return None
    result=(arb(4)/3-2*simplex*m+2*m*m-3*projected).lower()
    return dict(status='improved', gain_lower=result, squared_ratio_upper=upper_float(h*h),
                amplitude_bound=str(right), clip_lower=str(clipping),
                clip_slack_lower=str((clipping-h*h).lower()),
                physical_complement_upper=upper_float(complementary),
                mu_raw_trial_upper=upper_float(M), mu_raw_profile_upper=upper_float(projected),
                mu_mean_upper=upper_float(m))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start",type=Q,default=Q(129,200))
    parser.add_argument("--stop",type=Q,default=Q(81,125))
    parser.add_argument("--width",type=Q,default=Q(1,4000))
    parser.add_argument("--panels",type=int,default=512)
    parser.add_argument("--modes",type=int,default=192)
    parser.add_argument("--rows",type=int,default=96)
    parser.add_argument("--clip-panels",type=int,default=768)
    parser.add_argument("--raw-panels",type=int,default=4096)
    parser.add_argument("--cap",type=Q,default=Q(1999,1000))
    args=parser.parse_args()
    ctx.prec=110
    assert 0<args.cap<Q(199999,100000)
    assert 0<args.rows<=args.modes
    assert (arb.pi()/16).cos()**2>arb(24)/25
    assert (2*arb.pi()/9).cos()**2<arb(3)/5
    assert (arb(18)/25).exp()>2
    assert enclose(args.start)**2/(1-enclose(args.start)**2)>arb(1)/2
    assert enclose(args.stop)**2/(1-enclose(args.stop)**2)<1
    quad=HermiteQuadrature(args.panels,12,2*args.modes)
    angular=sharp_moments(args.modes)
    steps=(args.stop-args.start)/args.width
    assert steps.denominator==1 and steps>0
    slabs=[]
    def certify(left,right,depth):
        rho=enclose(left).union(enclose(right))
        boost=Q(12,25)
        records=[]
        for _ in range(5):
            stage=gain_bound(rho,boost,quad,angular,args.rows,args.clip_panels,args.raw_panels,enclose(args.cap))
            if stage is None:break
            if stage["status"]=="negative":
                slabs.append(dict(rho=[str(left),str(right)],stages=records,final_coefficient=str(boost),total_upper=stage["total_upper"]))
                return
            lower=stage.pop("gain_lower")
            updated=Q(int(float((1000*lower).floor().lower())),1000)
            if updated<=boost:break
            assert lower>enclose(updated)
            records.append(dict(current=str(boost),new_gain=str(updated),gain_lower=str(lower),**stage))
            boost=updated
        assert depth<9,(str(left),str(right),records)
        mid=(left+right)/2
        certify(left,mid,depth+1)
        certify(mid,right,depth+1)
    for i in range(steps.numerator):
        l=args.start+i*args.width
        certify(l,l+args.width,0)
    import check_three_cell_sharp_undamped_extension as base
    import check_three_cell_sharp_variance_bootstrap_anchor as parent_checker
    files=(Path(__file__),Path(__file__).with_name("sharper_hermite_radial_quadrature.py"),
           Path(base.__file__),Path(parent_checker.__file__),
           base.folder/"check_three_cell_global_scalar_through_half.py",
           base.folder/"check_three_cell_sharp_positive_arc_envelope.py")
    print(json.dumps(dict(status="pass",rho_bounds=[str(args.start),str(args.stop)],cap=str(args.cap),
        preliminary_radial_reserve="1/100000",mode_order=args.modes,row_order=args.rows,panels=args.panels,
        radial_limit=12,clip_panels=args.clip_panels,clip_limit=1,raw_panels=args.raw_panels,raw_limit=12,
        precision_bits=ctx.prec,total_upper=max(s["total_upper"] for s in slabs),slabs=slabs,
        source_sha256={f.name:sha256(f.read_bytes()).hexdigest() for f in files}),indent=2))

if __name__=="__main__":main()
