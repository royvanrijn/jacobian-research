#!/usr/bin/env python3
"""Exact supporting inequalities for the lower-two weighted saddle proof."""
from __future__ import annotations
import argparse, hashlib, json, math
from fractions import Fraction as F
from pathlib import Path
if not __debug__:
    raise SystemExit('Refusing optimized execution: proof assertions are required.')
HERE=Path(__file__).resolve().parent
RECORDS=[]

def rat(x): return F(x)
def show(x):
    x=F(x); return str(x.numerator) if x.denominator==1 else f'{x.numerator}/{x.denominator}'
def check(label, left, right, relation='<'):
    left,right=F(left),F(right)
    ok=left<right if relation=='<' else left<=right if relation=='<=' else left==right
    if not ok: raise AssertionError((label,show(left),relation,show(right)))
    RECORDS.append({'label':label,'left':show(left),'relation':relation,'right':show(right)})
def exp_bounds(x):
    x=F(x)
    if x<0:
        if -x>100: return F(0),1/exp_bounds(F(100))[0]
        lo,hi=exp_bounds(-x); return 1/hi,1/lo
    term=F(1); total=term
    # After 256 terms the remaining term ratios are <=x/258.
    for k in range(1,257): term*=x/k; total+=term
    nxt=term*x/257; q=x/258
    assert 0<=q<1
    return total,total+nxt/(1-q)
def eu(x): return exp_bounds(x)[1]
def el(x): return exp_bounds(x)[0]
def su(x):
    x=F(x); assert x>=0
    scale=10**12; k=math.isqrt(x.numerator*scale*scale//x.denominator)
    assert F(k*k,scale*scale)<=x<F((k+1)**2,scale*scale)
    return F(k+1,scale)
def sl(x):
    x=F(x); scale=10**12; k=math.isqrt(x.numerator*scale*scale//x.denominator)
    return F(k,scale)
def powshift(poly,l,h):
    n=len(poly)-1; out=[F(0)]*(n+1)
    for k,c in enumerate(poly):
        for j in range(k+1):out[j]+=c*math.comb(k,j)*l**(k-j)*(h-l)**j
    return out

def bernstein(poly,l,h):
    a=powshift(poly,F(l),F(h)); n=len(poly)-1
    return [sum((a[i]*F(math.comb(k,i),math.comb(n,i)) for i in range(k+1)),F(0)) for k in range(n+1)]
def sparsepoly(n,terms):
    p=[F(0)]*(n+1)
    for i,c in terms.items():p[i]=F(c)
    return p

def padd(a,b):
    n=max(len(a),len(b));return [(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(n)]
def pscale(a,c):return [v*c for v in a]
def pmul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,u in enumerate(a):
        for j,v in enumerate(b):out[i+j]+=u*v
    return out
def ppow(a,n):
    out=[F(1)]
    for _ in range(n):out=pmul(out,a)
    return out
def peq(a,b):
    d=padd(a,pscale(b,-1));return all(v==0 for v in d)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',type=Path,default=HERE/'lower_two_weighted_saddle_verification.json'); ap.add_argument('--repo-root',type=Path); args=ap.parse_args()
    repo=args.repo_root.resolve() if args.repo_root else next((p for p in (HERE,*HERE.parents) if (p/'research/MATH_STATUS.json').is_file()),None)
    if repo is None: raise SystemExit('Use --repo-root when running outside the repository.')
    ledger_path=HERE/'lower_two_weighted_saddle_source_ledger.json'
    ledger=json.loads(ledger_path.read_text())
    expected={
        'proof':('local','ORDINARY_LOWER_TWO_STRONG_FACTORIAL.md'),
        'checker':('local',Path(__file__).name),
        'constants':('local','lower_two_weighted_saddle_constants.json'),
        'sparse':('repository','research/factorial/image-rigidity/lower_exponent_initial_and_boundary.tex'),
        'prime_gap':('repository','research/factorial/image-rigidity/LOWER_TWO_PRIME_GAP_SFC.tex'),
        'fourth_gap':('repository','research/factorial/image-rigidity/LOWER_TWO_FOURTH_POWER_GAP_SFC.tex'),
        'adjacent':('repository','research/factorial/image-rigidity/UNIFORM_PADE_PARABOLA_SFC.tex')}
    assert set(ledger['files'])==set(expected),'missing or extra mandatory ledger entry'
    for key,(location,name) in expected.items():
        rec=ledger['files'][key];assert rec['location']==location and rec['path']==name
        source=(HERE if location=='local' else repo)/name
        assert hashlib.sha256(source.read_bytes()).hexdigest()==rec['sha256'],f'stale source: {key}'

    # Stationary branches, curvature/local estimates and monotonic ratio guards.
    check('e lower 163/60',F(163,60),el(1));check('e upper 87/32',eu(1),F(87,32));check('e lower 8/3',F(8,3),el(1))
    check('upper Rouche D8',F(85761906966528,93776847890045),1)
    check('upper Rouche successive ratio',F(1215,2119),1)
    check('lower Rouche D8',F(235,64)*F(1479,2048)**6,1)
    check('lower Rouche successive ratio',F(64,47)*F(87,32)*F(19,81),1)
    check('upper height radial',F(6144,71981),F(7,80))
    check('lower height',F(1,240)+F(1,46),F(1,32))
    q=F(69984,588245); y=q/(1-q)
    check('upper local third derivative',2*F(8,7)**3+64*(y+3*y*y+2*y**3),18)
    check('upper local q ratio',F(3,8)*F(9,8)**3,1)
    check('upper local D2q ratio',F(9,8)**2*F(3,8)*F(9,8)**3,1)
    q=F(5,64)*(F(87,32)*F(9,4)/8)**6;y=q/(1-q)
    check('lower local third derivative',F(4,7)**3*(4+216*(y+3*y*y+2*y**3)),2)
    check('lower local d3q ratio',F(7,6)**3*F(783,128)/18,1)
    check('half-power comparison',2,F(9,8)**6)
    check('upper central q step',F(4,7),1)
    check('upper central first correction ratio',F(32,27)*F(4,7),1)
    check('upper central second correction ratio',F(49,36)*F(4,7)**2,1)
    check('upper left h ratio',F(28,27)*F(6,7),1)
    check('lower right first correction ratio',F(7,6)*F(17,32)*F(5,4)/(1+F(17,128)),1)
    check('lower right second correction ratio',F(7,6)*F(17,32)**2,1)
    check('upper vertical height',(F(27,32)*eu(F(59,400))),1)
    check('upper vertical rational alternative',F(27,32)*F(400,341),1)
    check('upper far q',F(3,10)*64*F(3,10)**6,F(1,50))
    q=F(5,64)*(F(87,32)*F(251,100)/8)**6;y=q/(1-q)
    check('lower central curvature',F(189003,200000),F(39,20)-30*y-36*y*y)
    check('lower left small q',F(5,64)*(F(87,32)*F(151,100)/8)**6,F(1,500))
    check('lower left exponential',F(145,256)*F(501,500)*eu(F(1,2)+F(1,32)+F(1,100)),1)
    check('lower left exp shortcut',eu(F(11,20)),F(7,4))
    # The real-envelope extrema are endpoints of the indicated monotone pieces.
    for i,(l,h) in enumerate([(F(7,2),F(4)),(F(4),F(9,2)),(F(9,2),F(19,4)),(F(19,4),F(5)),(F(5),F(5))],1):
        env=(l/2)**2*eu(2-l)+F(5,4)*(h/8)**8*eu(8-h)
        check(f'lower far envelope interval {i}',env*eu(F(1,32)+F(1,196)+F(1,200)),1)
    # Full polynomial positivity certificates, not interpolation/sampling.
    A=234827336908800000;B=18755369578009
    polys=[('upper central',sparsepoly(14,{14:-5*B*B,12:1248*B*B,8:10*A*B,6:-9216*A*B,2:-5*A*A,0:2208*A*A}),F(7),F(10)),
           ('upper left complex',sparsepoly(15,{15:-101*B*B,14:800*B*B,12:-200*B*B,8:-600*A*B,3:101*A*A,2:-200*A*A}),F(5),F(7))]
    a=6022586125125;b=49539595901075456
    polys.append(('lower right complex',sparsepoly(15,{15:-152*a*a,14:1280*a*a,12:-15*a*a,8:-960*a*b,3:152*b*b,2:-320*b*b}),F(5,2),F(7,2)))
    certs=[]
    for name,p,l,h in polys:
        coeff=bernstein(p,l,h); assert all(v>0 for v in coeff)
        certs.append({'label':name,'interval':[show(l),show(h)],'power_coefficients':[show(v) for v in p],'bernstein_coefficients':[show(v) for v in coeff]})
    assert sum(len(c['bernstein_coefficients']) for c in certs)==47
    # Independently clear the rational Mobius/curvature denominators using
    # polynomial arithmetic, so positivity is linked to the actual estimates.
    x=[F(0),F(1)]; X6=pscale(ppow(x,6),B); AA=[F(A)]
    mm=padd(X6,pscale(AA,-1)); plus=padd(X6,AA)
    uc=padd(pmul(padd([F(1248)],pscale(ppow(x,2),-5)),ppow(mm,2)),
            padd(pscale(pmul(AA,mm),-6720),pscale(ppow(AA,2),-5760)))
    assert peq(uc,polys[0][1])
    # UL: 100*x^3*(1-r^2)*(lower_bound-1/100), r=A/(B*x^6).
    base=padd(pscale(ppow(x,2),800),pscale(ppow(x,3),-101))
    ul=padd(pmul(base,padd(ppow(X6,2),pscale(ppow(AA,2),-1))),
            padd(pscale(pmul(ppow(x,2),pmul(AA,mm)),-600),pscale(ppow(X6,2),-200)))
    assert peq(ul,polys[1][1])
    # LR: -160*x^3*(1-r^2)*(upper_bound+1/20), r=a*x^6/b.
    R=pscale(ppow(x,6),a);BB=[F(b)];minus=padd(BB,pscale(R,-1))
    lr=padd(pmul(padd(pscale(ppow(x,3),152),pscale(ppow(x,2),-320)),padd(ppow(BB,2),pscale(ppow(R,2),-1))),
            padd(pscale(pmul(ppow(x,2),pmul(R,minus)),-960),pscale(ppow(R,2),-15)))
    assert peq(lr,polys[2][1])
    # Local upper decay: exact numerator in z=D-8.
    z=[F(8),F(1)]
    decay=padd(padd(pscale(ppow(z,3),6),pscale(ppow(z,2),-55)),padd(pscale(z,82),[F(-21)]))
    assert peq(decay,[F(187),F(354),F(89),F(6)])
    # Numerical packet budgets: monotonicity reductions are proved in the note.
    check('EIU D<=17 local',147/sl(32000),F(33,40))
    C=F(4000,17)
    check('EIU D<=17 completion',14/sl(C)*eu(-C/28),F(1,1000))
    check('EIU D<=17 central',16/sl(C)*eu(-C/32),F(1,1000))
    check('EIU D>=18 local',147/sl(72000),F(29,40))
    check('EIU D>=18 completion',14/sl(128)*eu(F(-128,28)),F(3,200))
    check('EIU D>=18 central',16/sl(128)*eu(-4),F(3,100))
    check('EIU left',100/sl(32000)*eu(F(-32000,5120)),F(1,500))
    check('EIU vertical',F(1,2)*su(500)*eu(-80),F(1,10000))
    check('EIU far',10/sl(32000)*eu(-25),F(1,10000))
    check('EIU low sum',F(33,40)+F(2,1000)+F(1,500)+F(2,10000),F(5,6))
    check('EIU high sum',F(29,40)+F(3,200)+F(3,100)+F(1,500)+F(2,10000),F(5,6))
    check('EJU local bracket',F(5,3)*F(1000,999)*(F(33,40)+F(5,8)),F(5,2))
    check('EJU sqrtD/n guard',7/sl(F(128*999,1000)),F(5,8))
    C=F(128*999,1000)
    check('EJU completion',14/sl(C)*eu(-C/28),F(1,50))
    check('EJU central',32/sl(128)*eu(-2),F(2,5))
    check('EJU left',200/sl(32000)*eu(F(-32000,10240)),F(1,20))
    check('EJU vertical',F(1,2)*su(500)*eu(-40),F(1,1000))
    check('EJU far',10/sl(32000)*eu(F(-32000,2560)),F(1,1000))
    check('EJU sum',F(5,2)+F(1,50)+F(2,5)+F(1,20)+F(2,1000),3)
    lowerI=[24/sl(4000),48/sl(4000)*eu(F(-4000,768)),160/(3*sl(4000))*eu(F(-12000,2560)),F(3,2)*su(4000)*eu(-40),10/sl(4000)*eu(F(-36000,1280)),F(4000,128)*su(4000)*eu(-20)]
    capsI=[F(38,100),F(5,1000),F(1,100),F(1,1000),F(1,1000),F(1,1000)]
    for i,(v,c) in enumerate(zip(lowerI,capsI),1):check(f'EIL term {i}',v,c)
    check('EIL sum',sum(capsI),F(2,5))
    lowerJ=[240/(7*sl(3999)),48/sl(3999)*eu(F(-3999,768)),320/(3*sl(4000))*eu(F(-12000,5120)),F(3,2)*su(4000)*eu(F(17,8)+F(1,100)-40),10/sl(3999)*eu(F(-9*3999,1280)),F(4000,128)*su(4000)*eu(-20)]
    capsJ=[F(545,1000),F(5,1000),F(165,1000),F(1,1000),F(1,1000),F(1,1000)]
    for i,(v,c) in enumerate(zip(lowerJ,capsJ),1):check(f'EJL term {i}',v,c)
    check('EJL sum',sum(capsJ),F(18,25))
    check('inner height separation',eu(F(7,32)),F(5,4))
    check('outer height separation',eu(F(27,160)),F(6,5))
    check('inner zero-free ratio',54*eu(-400),1)
    check('outer zero-free ratio',4000*eu(-200),1)
    check('saddle amplitude separation',150,F(5,6)*el(6-F(119,160)))
    check('weighted pair ratio',F(120,150),1)
    check('sparse overlap d46',49152,el(F(2116,192)))
    check('sparse overlap derivative D48',0,F(1,4)-F(1,48)-F(1,48**2))
    assert sum(max(4000,128*(d+2))-2 for d in range(6,46))==176816
    assert max(max(4000,128*(d+2))-1 for d in range(6,46))==6015
    assert max((d+2)*(max(4000,128*(d+2))-1) for d in range(6,46))==282705
    # Deliberate false controls: old broad central concavity really fails.
    D=8;x=6;w=20000;q=F(w,x**6)
    curvature=-F(D,x*x)+F(6,x*x)*q*(7+q)/(1+q)**2
    assert curvature==F(2952569,78100002)>0
    bad=polys[0][1][:];bad[0]=-sum(abs(c)*10**(i+1) for i,c in enumerate(bad))
    assert not all(c>0 for c in bernstein(bad,7,10))
    assert not F(5,2)+F(1,50)+F(2,5)+F(1,20)+F(2,1000)<F(29,10)
    # Persist every exact polynomial certificate alongside the receipt.
    constants=HERE/'lower_two_weighted_saddle_constants.json'
    payload={'schema_version':1,'bernstein_certificates':certs,'onset':{'minimum_terminal':4000,'terminal_per_D':128,'minimum_gap':6},'finite_domain':{'minimum_gap':6,'maximum_gap':45,'case_count':176816,'maximum_terminal':6015,'maximum_factorial':282705}}
    assert constants.is_file(),'missing mandatory exact constants input'
    assert json.loads(constants.read_text())==payload,'constants input differs from independently reconstructed exact data'
    proof=HERE/'ORDINARY_LOWER_TWO_STRONG_FACTORIAL.md'
    assert proof.exists()
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    receipt={'schema_version':1,'status':'PASS','scope':'Exact supporting inequalities for the weighted actual consecutive lower-two saddle proof; finite coprimality is a separate computation.','counts':{'exact_comparisons':len(RECORDS),'bernstein_coefficients':47,'full_bernstein_polynomials':3,'cleared_polynomial_identities':4,'negative_controls':3},'files':[{'path':p.name,'sha256':sha(p)} for p in (proof,Path(__file__).resolve(),constants,ledger_path)],'comparisons':RECORDS}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','counts':receipt['counts'],'output':str(args.output)},sort_keys=True))
if __name__=='__main__':main()
