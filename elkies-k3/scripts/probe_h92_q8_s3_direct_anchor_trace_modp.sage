#!/usr/bin/env sage -python
"""
Direct anchored Abel-Jacobi trace for H92 q8 S3 modulo p.

Uses the exact II*_E8_1 branch-point zero.  The degree-52 S3 divisor is mapped
birationally from the q8 quartic directly to the anchored canonical D13 curve,
then summed using L(53 O).  No binary-quartic covariant 2-cover, IV* origin
subtraction, or halving is used.
"""

import argparse
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix, sage_eval


def locate_repo(explicit=None):
    candidates=[]
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd=Path.cwd().resolve()
    candidates += [cwd,*cwd.parents]
    home=Path.home()
    candidates += [
        home/'Documents'/'jacobian-research', home/'jacobian-research',
        home/'src'/'jacobian-research', home/'git'/'jacobian-research',
        home/'projects'/'jacobian-research',
    ]
    seen=set()
    for candidate in candidates:
        try: candidate=candidate.resolve()
        except Exception: continue
        if candidate in seen: continue
        seen.add(candidate)
        if (candidate/'elkies-k3/scripts').is_dir() and (candidate/'artifacts/generated-results').is_dir():
            return candidate
    raise SystemExit('Could not locate jacobian-research; pass --repo PATH')


parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--repo',type=Path)
parser.add_argument('--prime',type=int,default=100003)
parser.add_argument('--tau',type=int,default=2)
parser.add_argument('--output',type=Path)
args=parser.parse_args()

ROOT=locate_repo(args.repo)
GEN=ROOT/'artifacts/generated-results'
LOCAL=ROOT/'artifacts/local/elkies-k3'
CORE=ROOT/'elkies-k3/scripts/elliptic_neighbor_compiler.sage'
Q6=GEN/'elkies-k3-h92-q6-child-jacobian.json'
Q8=LOCAL/'q8-corrected2cover-qq-child.json'
if not Q8.exists():
    Q8=GEN/'elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json'
BRIDGE=LOCAL/'q6-third-to-q8-bridge.json'
ANCHOR=LOCAL/'q8-d13-branch-anchor.json'
EQNS=LOCAL/'q8-equation-ns-divisor.json'
OUTPUT=(args.output.resolve() if args.output else LOCAL/f'q8-s3-direct-anchor-trace-mod-{args.prime}-tau-{args.tau}.json')
for path in (CORE,Q6,Q8,BRIDGE,ANCHOR,EQNS):
    if not path.exists(): raise SystemExit(f'Missing prerequisite: {path}')

scope={}
exec(compile(CORE.read_text(),str(CORE),'exec'),scope)
squarefree_binary_quartic=scope['squarefree_binary_quartic']

q6=json.loads(Q6.read_text())
q8=json.loads(Q8.read_text())
bridge=json.loads(BRIDGE.read_text())
anchor=json.loads(ANCHOR.read_text())
eqns=json.loads(EQNS.read_text())
assert q6['status']=='PASS_EXACT_E8_E6_CHILD_JACOBIAN'
assert q8['status']=='PASS_EXACT_CORRECTED_Q8_D13_CHILD'
assert bridge['status']=='PASS_EXACT_Q6_THIRD_TO_Q8_DEGREE52'
assert anchor['status']=='PASS_EXACT_D13_BRANCH_ANCHOR'
assert eqns['status']=='PASS_EXACT_Q8_EQUATION_NS_DIVISOR'
assert eqns['q8_equation_fibre']['S3_degree']==52
assert eqns['q8_equation_fibre']['root_data']==[13,312,4]

p=ZZ(args.prime)
if not p.is_prime() or p in (2,3): raise ValueError('prime must be odd and !=3')
F=GF(p)
R=PolynomialRing(F,'T'); T=R.gen(); K=R.fraction_field()


def modq(value):
    value=QQ(value); den=ZZ(value.denominator())
    if den%p==0: raise ZeroDivisionError(f'denominator divisible by {p}: {value}')
    return F(ZZ(value.numerator()))/F(den)


def poly_strings(values): return R([modq(v) for v in values])

def rf_from_bridge(entry):
    return K(poly_strings(entry['numerator_coefficients_low_to_high']))/K(poly_strings(entry['denominator_coefficients_low_to_high']))

# q6 curve + exact S3.
A6=poly_strings(q6['minimal_short_weierstrass']['A_coefficients_low_to_high'])
B6=poly_strings(q6['minimal_short_weierstrass']['B_coefficients_low_to_high'])
section3=bridge['third_section_canonical_q6']
x3=rf_from_bridge(section3['x']); y3=rf_from_bridge(section3['y'])
assert y3**2==x3**3+K(A6)*x3+K(B6)
U3_num=poly_strings(bridge['q8_parameter_on_third']['numerator_coefficients_low_to_high'])
U3_den=poly_strings(bridge['q8_parameter_on_third']['denominator_coefficients_low_to_high'])
assert max(U3_num.degree(),U3_den.degree())==52

# Corrected marked section and q8 pencil frame.
mdata=q8['marking']['section']
sx=K(poly_strings(mdata['x_numerator_coefficients_low_to_high']))/K(poly_strings(mdata['x_denominator_coefficients_low_to_high']))
sy=K(poly_strings(mdata['y_numerator_coefficients_low_to_high']))/K(poly_strings(mdata['y_denominator_coefficients_low_to_high']))
assert sy**2==sx**3+K(A6)*sx+K(B6)


def monic_power_root(value,exponent):
    out=R.one()
    for factor,mult in value.factor():
        assert mult%exponent==0
        out*=factor.monic()**(mult//exponent)
    return out.monic()

nx,dx=R(sx.numerator()),R(sx.denominator())
ny,dy=R(sy.numerator()),R(sy.denominator())
h=monic_power_root(dx,2)
assert h==monic_power_root(dy,3) and h.degree()==10
QQTR=PolynomialRing(QQ,'T')
def reduce_factor_string(text):
    src=QQTR(text); return R([modq(v) for v in src.list()]).monic()
ii=reduce_factor_string(next(item for item in q6['finite_fibres'] if item['kodaira']=='II*')['factor'])
iv=reduce_factor_string(next(item for item in q6['finite_fibres'] if item['kodaira']=='IV*')['factor'])
M=(ii**2*iv**2).monic()
normalizer=(ny*dx*(h*dy).inverse_mod(nx)).mod(nx)
p_fun=-sy/sx
rho=(normalizer*nx.inverse_mod(M)).mod(M)
pairs=[]
for entry in q8['rr']['kernel_polynomials']:
    sp=R([modq(v) for v in QQTR(entry['s']).list()])
    tp=R([modq(v) for v in QQTR(entry['t']).list()])
    Bcoef=K(sp)/K(h)
    Acoef=(-K(sp)*p_fun/K(h)-K(sp)*K(normalizer)/K(nx)+K(sp*rho)+K(tp*M))
    pairs.append((Acoef,Bcoef))
(A0,B0),(A1,B1)=pairs

# D13 specialization.
tau=F(args.tau)
A13p=[modq(v) for v in q8['child']['minimal_A_coefficients_low_to_high']]
B13p=[modq(v) for v in q8['child']['minimal_B_coefficients_low_to_high']]
def eval_coeff(values,at): return sum(v*at**i for i,v in enumerate(values))
A13=eval_coeff(A13p,tau); B13=eval_coeff(B13p,tau)
E13=EllipticCurve(F,[0,0,0,A13,B13])
if not E13.discriminant(): raise ArithmeticError('tau gives singular D13 child fibre')

H=R(U3_num-tau*U3_den)
if H.degree()!=52: raise ArithmeticError(f'S3 q8 fibre degree dropped to {H.degree()}')
H=H.monic()
if H.gcd(H.derivative()).degree()!=0: raise ArithmeticError('degree-52 S3 fibre is not etale')
print(f'Q8S3DIRECT|prime={p}|tau={int(tau)}|degree=52|stage=setup|status=PASS',flush=True)

q8_m=-(A1-tau*A0)/(B1-tau*B0)
radicand=q8_m**4-6*sx*q8_m**2-8*sy*q8_m-3*sx**2-4*K(A6)
quartic,square_factor=squarefree_binary_quartic(radicand,R)
if quartic.degree()!=4: raise ArithmeticError(f'q8 quartic degree dropped to {quartic.degree()}')


def reduce_mod_H(value):
    value=K(value); num=R(value.numerator()); den=R(value.denominator())
    if den.gcd(H).degree()!=0: raise ZeroDivisionError('denominator not invertible modulo degree-52 fibre')
    return (num*den.inverse_mod(H))%H

m3=(y3+sy)/(x3-sx)
if reduce_mod_H(q8_m-m3): raise ArithmeticError('q8 chord mismatch on S3 degree-52 fibre')
w3=(2*x3+sx-q8_m**2)/square_factor
wA=reduce_mod_H(w3)
if (wA*wA-quartic)%H: raise ArithmeticError('S3 q8 quartic square-root mismatch')

# Specialize exact branch-anchor data over QQ(U).
QUQ=PolynomialRing(QQ,'U'); UQ=QUQ.gen(); KUQ=QUQ.fraction_field()
def parse_u(text): return KUQ(sage_eval(str(text),locals={'U':UQ}))
def spec_u(text):
    v=parse_u(text); num=QUQ(v.numerator()); den=QUQ(v.denominator())
    n=sum(modq(c)*tau**i for i,c in enumerate(num.list()))
    d0=sum(modq(c)*tau**i for i,c in enumerate(den.list()))
    if not d0: raise ZeroDivisionError('branch-anchor U denominator vanished')
    return n/d0

tii=modq(QQ(anchor['zero']['old_base_T']))
coef=anchor['quartic_to_anchor']['shifted_coefficients']
a=spec_u(coef['a_r4']); b=spec_u(coef['b_r3']); c=spec_u(coef['c_r2']); d=spec_u(coef['d_r1'])
rpoly=T-tii
branch_poly=d*rpoly+c*rpoly**2+b*rpoly**3+a*rpoly**4
if branch_poly.degree()!=4: raise ArithmeticError('anchored branch quartic degree dropped')
scale=quartic[4]/branch_poly[4]
if quartic != scale*branch_poly:
    raise ArithmeticError('squarefree quartic and anchored quartic are not scalar-equivalent')
if not scale.is_square(): raise ArithmeticError(f'quartic scale is nonsquare: {scale}')
scale_root=scale.sqrt()
wbA=wA/scale_root
if (wbA*wbA-branch_poly)%H: raise ArithmeticError('branch W conversion failed')

rA=rpoly%H
if rA.gcd(H).degree()!=0: raise ZeroDivisionError('S3 divisor meets branch zero at chosen tau')
rInv=rA.inverse_mod(H)
Xa=(d*rInv)%H
Ya=(d*wbA*rInv**2)%H

urst=anchor['anchor_to_canonical']['urst']
u=spec_u(urst[0]); rr=spec_u(urst[1]); ss=spec_u(urst[2]); tt=spec_u(urst[3])
if not u: raise ZeroDivisionError('anchor isomorphism u vanished')
xA=((Xa-rr)/(u**2))%H
yA=((Ya-ss*(Xa-rr)-tt)/(u**3))%H
if (yA*yA-xA*xA*xA-A13*xA-B13)%H:
    # The only remaining global ambiguity is W -> -W.  Try the opposite once;
    # if the stored branch-anchor orientation is opposite to square_factor's.
    Ya=(-d*wbA*rInv**2)%H
    yA=((Ya-ss*(Xa-rr)-tt)/(u**3))%H
    if (yA*yA-xA*xA*xA-A13*xA-B13)%H:
        raise ArithmeticError('direct anchored S3 images miss D13 child')
    w_sign=-1
else:
    w_sign=1
print(f'Q8S3DIRECT|quartic_scale={scale}|W_sign={w_sign}|stage=direct_transport|status=PASS',flush=True)

# Sum the 52 D13 points via L(53 O).
one=R.one(); xp=[one]
for _ in range(26): xp.append((xp[-1]*xA)%H)
columns=list(xp)+[(yA*xp[e])%H for e in range(26)]
assert len(columns)==53
Eval=matrix(F,52,53,lambda row,col: columns[col][row])
ker=Eval.right_kernel().basis_matrix()
if ker.nrows()!=1: raise ArithmeticError(f'L(53O) trace kernel dimension {ker.nrows()}')
rel=ker[0]
XR=PolynomialRing(F,'X'); X=XR.gen()
Afun=sum(rel[i]*X**i for i in range(27))
Bfun=sum(rel[27+i]*X**i for i in range(26))
Rint=Afun**2-(X**3+A13*X+B13)*Bfun**2
if Rint.degree()!=53: raise ArithmeticError(f'residual intersection degree {Rint.degree()}, expected 53')
root_sum=-Rint[52]/Rint[53]

def newton_power_sums(poly):
    n=poly.degree(); assert poly[n]==1
    sums=[F(n)]
    for k in range(1,n):
        total=F(k)*poly[n-k]
        for j in range(1,k): total += poly[n-j]*sums[k-j]
        sums.append(-total)
    return sums
ps=newton_power_sums(H)
trace_x=sum(xA[i]*ps[i] for i in range(52))
xQ=root_sum-trace_x
bQ=Bfun(xQ)
if not bQ: raise ArithmeticError('trace residual has B(x_Q)=0')
yQ=-Afun(xQ)/bQ
Qres=E13(xQ,yQ)
AJ=-Qres
if AJ.is_zero():
    point={'zero':True}
    xout=yout='ZERO'
else:
    ax,ay=AJ.xy(); point={'zero':False,'x':int(ax),'y':int(ay)}; xout=str(int(ax)); yout=str(int(ay))
print(f'Q8S3DIRECT_RESULT|prime={p}|tau={int(tau)}|AJ_x={xout}|AJ_y={yout}|status=PASS_DIRECT_ANCHORED_AJ',flush=True)

payload={
    'schema':'elkies-k3.h92-q8-s3-direct-anchor-trace-modp.v1',
    'status':'PASS_DIRECT_ANCHORED_Q8_S3_AJ_TRACE',
    'prime':int(p),'tau':int(tau),'degree':52,'quartic_scale':str(scale),'W_sign':w_sign,
    'AJ':point,
    'method':'II*_E8_1 branch-point birational map; direct L(53O) sum; no covariant 2-cover or halving',
    'inputs':{'equation_ns_status':eqns['status'],'branch_anchor_status':anchor['status']},
}
OUTPUT.parent.mkdir(parents=True,exist_ok=True)
OUTPUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\\n')
print(f'OUTPUT|{OUTPUT}',flush=True)
