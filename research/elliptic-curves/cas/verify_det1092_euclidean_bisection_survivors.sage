#!/usr/bin/env sage-python
"""Independent exact trace, Euclidean-line and rational nonsquare replay.

Manual Q(t) group law, reversed word order, coefficient identities and
integer inequalities; no producer imports, factorization, is_square,
nullspace solver, new prime or point/parameter search. 25-second cap.
"""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,vector
signal.alarm(25)
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_euclidean_bisection_survivors_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
PRE=ART/'det1092_modular_bisection_exclusion_v1'
MOD=ART/'det1092_modular_bisection_direct_v1'
REG=ART/'det1092_rational_bisection_index_v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bind(record):
    for name,h in record['inputs'].items():assert sha(ROOT/name)==h,name
protocol=read(OUT/'protocol.json');summary=read(OUT/'summary.json')
bind(protocol);bind(summary)
assert protocol['masks']==[61,107,111] and protocol['limits']['incidence_tests']==27
assert protocol['limits']['new_primes']==protocol['limits']['new_parameters']==0
parent=read(PARENT);R=PolynomialRing(QQ,'t');K=R.fraction_field()
def dec(r):return K(R(r['numerator']))/R(r['denominator'])
a1,a2,a3,a4,a6=[dec(r) for r in parent['a_invariants']]
b2=a1*a1+4*a2;b4=a1*a3+2*a4;b6=a3*a3+4*a6
A=R(b4/2-b2*b2/48);B=R(b6/4-b2*b4/24+b2**3/864)
basis=[]
for pair in parent['basis_weierstrass_coordinates']:
    x,y=map(dec,pair);xx=x+b2/12;yy=y+(a1*x+a3)/2
    assert yy*yy==xx**3+A*xx+B;basis.append((xx,yy))
def add(P,Q):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    if x==u:
        if y==-v:return None
        assert y==v and y
        slope=(3*x*x+A)/(2*y)
    else:slope=(v-y)/(u-x)
    z=slope*slope-x-u
    return (z,-y+slope*(x-z))
def mul(n,P):
    if n<0:return mul(-n,(P[0],-P[1]))
    Q=None
    while n:
        if n&1:Q=add(Q,P)
        n//=2;P=add(P,P)
    return Q
G=matrix(ZZ,parent['generic_height_gram']);outcomes=[]
for mask in protocol['masks']:
    data=read(OUT/('orbit-%d.json'%mask));tr=read(PRE/('trace-%d.json'%mask))
    assert data['mask']==tr['mask']==mask and data['word']==tr['word']
    assert data['trace_sha256']==sha(PRE/('trace-%d.json'%mask))
    w=vector(ZZ,data['word']);assert w*G*w==10
    P=None
    for n,Q in reversed(list(zip(w,basis))):P=add(P,mul(-int(n),Q))
    assert P is not None
    h,nx,ny=[R(tr[k]) for k in ['pole','nx','ny']]
    assert h.is_monic() and h.degree()==3 and nx.gcd(h).degree()==0
    assert P==(K(nx)/h**2,K(ny)/h**3)
    assert ny*ny==nx**3+A*nx*h**4+B*h**6
    f0,f1,f2=map(R,data['line']);q=R(data['branch_quadratic'])
    assert f2==h and f1.degree()<6 and f0.degree()<=9
    assert f0*h*h+f1*nx+ny==0  # Unique inverse-mod-h^2 line, without inversion.
    assert h**6*q==f1**4-6*nx*f1**2-8*ny*f1-3*nx**2-4*A*h**4
    assert q.degree()==2 and q.gcd(q.derivative()).degree()==0
    x0,x1,y0,y1=map(R,data['elliptic_maps'])
    assert x1==h/2 and y1==-f1/2
    assert 2*h*h*x0==f1*f1-nx
    assert h*y0+f0+f1*x0==0
    # Direct polynomial quotient-ring substitution, not separate copied
    # constant/cross coefficient formulas from the producer.
    RW=PolynomialRing(K,'W');W=RW.gen()
    XX=x0+x1*W;YY=y0+y1*W
    assert (YY*YY-XX**3-A*XX-B)%(W*W-q)==0
    # Descent to the old fixed151 line; no new prime is introduced.
    Fp=GF(151);Rp=PolynomialRing(Fp,'t')
    vec=vector(Fp,[f0[i] for i in range(10)]+[f1[i] for i in range(6)]+[h[i] for i in range(4)])
    old=next(r for r in read(MOD/('orbit-%d.json'%mask))['trials'] if r['prime']==151)
    kv=vector(Fp,old['kernel']);index=next(i for i,v in enumerate(kv) if v)
    scale=vec[index]/kv[index];assert scale and vec==scale*kv
    exposure=read(OUT/('incidence-%d.json'%mask))
    assert exposure['curve_sha256']==sha(OUT/('orbit-%d.json'%mask))
    assert [{'label':c['label'],'parameter':c['parameter']} for c in exposure['cases']]==protocol['roster']
    for case in exposure['cases']:
        tt=QQ(case['parameter']);value=q(tt);cert=case['square_test']
        assert value==QQ(cert['value']) and -16*(4*A(tt)**3+27*B(tt)**2)
        assert value>0  # All27 cases are real-soluble, not sign obstructions.
        n,d=ZZ(cert['numerator']),ZZ(cert['denominator'])
        a,b=ZZ(cert['numerator_floor_sqrt']),ZZ(cert['denominator_floor_sqrt'])
        assert d>0 and n.gcd(d)==1 and value==QQ(n)/d
        assert a>=0 and b>=0 and a*a<=n<(a+1)**2 and b*b<=d<(b+1)**2
        assert a*a!=n or b*b!=d
        assert cert['square'] is False and case['outcome']=='EXACT_NONSPLIT'
        assert 'short_points' not in case
    target=next(c for c in exposure['cases'] if QQ(c['parameter'])==0)
    outcomes.append({'mask':mask,'exact_generic_trace_word':True,'q_degree':2,
        'real_soluble_control_fibres':9,'rationally_nonsplit_control_fibres':9,
        'target302_value':target['square_test']['value']})
assert summary['trials']==27 and summary['nonsplit_count']==27 and not summary['splits']
# Assemble the completed fifteen-orbit panel, keeping prior certificate
# dependencies explicit rather than claiming to repeat all old calculations.
modular=read(MOD/'independent-replay.json');bind(modular)
mod_summary=read(MOD/'summary.json');bind(mod_summary)
old_regressions=read(REG/'controls.json');bind(old_regressions)
closed=[]
for m in mod_summary['new_masks']:
    if m in protocol['masks']:source='new exact rational nonsquare certificate'
    else:
        assert m in modular['new_excluded']
        assert any(r['mask']==m for r in modular['certificates'])
        source='previous independently replayed modular non-incidence'
    closed.append({'mask':m,'source':source,'at302':'EXACT_NONSPLIT'})
for m in mod_summary['regression_masks']:
    row=next(r for r in old_regressions['rows'] if r.get('orbit',r.get('mask'))==m)
    case=next(c for c in row['cases'] if QQ(c['parameter'])==0)
    assert not case['split']
    closed.append({'mask':m,'source':'previous exact conic/control certificate','at302':'EXACT_NONSPLIT'})
assert len(closed)==len({r['mask'] for r in closed})==15
report={'status':'PASS_INDEPENDENT_EUCLIDEAN_CONICS_AND_CLOSED_15_ORBIT_PANEL',
    'classification':'verified exact construction and exclusions; general sixth-order cancellation is proved in the note',
    'cases':outcomes,'new_incidence_decisions':27,'all27_real_soluble':True,
    'closed_panel':closed,'full_atlas_executed':False,
    'boundary':'Three new exact conic equations, no rational target/control split or seed. Prior modular and old-regression proofs are explicit dependencies. No new rank, prime, rational parameter or later exceptional input.',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [OUT/'protocol.json',OUT/'summary.json',
        MOD/'independent-replay.json',REG/'controls.json',Path(__file__)]}}
dest=OUT/'independent-replay.json';encoded=json.dumps(report,indent=2,sort_keys=True)+'\n'
if dest.exists():assert read(dest)==report
else:
    with dest.open('x') as stream:stream.write(encoded)
print(report['status'],'three exact trace words,27 nonsquares,15 closed orbits',flush=True)
