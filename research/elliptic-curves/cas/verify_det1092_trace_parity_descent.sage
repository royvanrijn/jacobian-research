#!/usr/bin/env sage-python
"""Symbolic and pinned-input checks for the trace-parity descent theorem.

The cohomological equivalence and uniform geometric genus are written proofs.
Only generic coordinates are used in algebra; old split-control conclusions
are read as retrospective evaluation. No exception-point reconstruction.
"""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,matrix,vector,EllipticCurve
signal.alarm(25)
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_trace_parity_descent_v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def provenance(row):
    for name,digest in row['inputs'].items():assert sha(ROOT/name)==digest
protocol=read(OUT/'protocol.json');provenance(protocol)
parent=read(ART/'curve302_recovered_mw17_parent_v1.json')
geometry=read(ART/'curve302_parent_geometric_picard19_v1.json')
orbits=read(ART/'curve302_parent_degree2_multisection_lattice_v1.json')
conic=read(ART/'det1092_trace_twist_kummer_v1/construction.json')
conic_replay=read(ART/'det1092_trace_twist_kummer_v1/replay.json');provenance(conic_replay)
norm8=read(ART/'det1092_norm8_seed_cover_v2/generic.json');provenance(norm8)
first=read(ART/'det1092_rr_net_halving_gate_v1.json')
first_replay=read(ART/'det1092_rr_net_halving_gate_replay_v1.json');provenance(first_replay)
assert geometry['status']=='PASS'
assert geometry['geometric_generic_MW_rank']==geometry['arithmetic_generic_MW_rank']==17
assert geometry['full_geometric_basis_is_displayed_rational_basis'] is True
assert conic_replay['status']=='PASS_INDEPENDENT_TRACE_TWIST_KUMMER_OBSTRUCTION'
assert first_replay['status']=='PASS_INDEPENDENT_HALVING_GENUS9_AND_SINGULAR_INCIDENCE'

# Universal chord-norm identity, checked as an exact polynomial identity.
B=PolynomialRing(QQ,names=('a','m','u','v','d','n'))
a,m,u,v,d,n=B.gens();S=PolynomialRing(B,'x');x=S.gen()
cx=m*m-a-2*u;line=m*x+n;g=(x-u)**2-v*v*d
f=line*line+(x-cx)*g
assert f.is_monic() and f.degree()==3 and f[2]==a
assert g*(cx-x)-line*line==-f
assert (u-x)**2-v*v*d==g

# Even-trace construction: set f(r)=d*e^2 and verify the cubic norm.
B=PolynomialRing(QQ,names=('a','b','r','d','e'))
a,b,r,d,e=B.gens();S=PolynomialRing(B,'x');x=S.gen()
c=d*e*e-r**3-a*r*r-b*r;f=x**3+a*x*x+b*x+c
alpha=d*(r-x)
assert matrix(B,3,3,lambda i,j:((alpha*x**j)%f)[i]).det()==d**4*e**2

# Universal halving quartic. No special centre is fitted or enumerated.
B=PolynomialRing(QQ,names=('a','cx','cy'));a,cx,cy=B.gens()
b=cy*cy-cx**3-a*cx
S=PolynomialRing(B,'m');m=S.gen()
H=m**4-6*cx*m*m-8*cy*m-3*cx*cx-4*a
assert H.discriminant()==256*(-16*(4*a**3+27*b*b))
xx=(m*m-cx)/2;yy=m*(xx-cx)-cy
assert (yy*yy-xx**3-a*xx-b)%H==0
assert (2*m*yy-3*xx*xx-a)%H==0

R=PolynomialRing(QQ,'t');K=R.fraction_field()
def dec(v):return K(R(v['numerator']))/R(v['denominator'])
ai=list(map(dec,parent['a_invariants']));E=EllipticCurve(K,ai)
A=-E.c4()/48;BB=-E.c6()/864
disc=E.discriminant();dn=disc.numerator()
assert dn.degree()==24 and dn.gcd(dn.derivative())==1
assert R(A).degree()==8 and R(BB).degree()==12
assert -16*(4*R(A)[8]**3+27*R(BB)[12]**2)!=0
G=matrix(ZZ,parent['generic_height_gram']);assert G.nrows()==17 and G.det()==1092
nh=K(R(norm8['pole_h']))
ncx=K(R(norm8['nx']))/nh**2;ncy=K(R(norm8['ny']))/nh**3
ccx,ccy=map(dec,conic['trace'])
centres=[('historical-first-centre',first['centre_word'],*map(dec,first['centre_short'])),
         ('norm8-orbit20124',norm8['selection']['word'],ncx,ncy),
         ('conic-orbit8044',conic['trace_word'],ccx/4+E.b2()/12,ccy/8)]
rows=[]
for label,word,cx,cy in centres:
    w=vector(ZZ,word);h=w*G*w
    assert h in (8,10) and any(v%2 for v in w)
    assert QQ(h)/4<4 and cy*cy==cx**3+A*cx+BB
    assert cx.denominator().gcd(dn)==cy.denominator().gcd(dn)==1
    rows.append({'label':label,'trace_word':list(map(int,w)),
        'trace_parity':[int(v%2) for v in w],'height':int(h),
        'generic_branch_class_descends':False,'halving_cover_degree':4,'halving_cover_genus':9})
assert [r['height'] for r in rows]==[10,8,10]
# Fresh equation-side norm check for the existing conic, with no specialization.
u,v=map(dec,conic['quadratic_point']['X']);n,m=map(dec,conic['line']);d=dec(conic['d'])
S=PolynomialRing(K,'x');x=S.gen();a,b,c=map(dec,conic['cubic_coefficients'])
f=x**3+a*x*x+b*x+c;g=(x-u)**2-v*v*d;line=m*x+n
assert g*(ccx-x)-line*line==-f
assert f(ccx)==ccy*ccy and v and d

assert orbits['status']=='PASS_COMPLETE_DEGREE2_TRANSLATION_QUOTIENT'
assert orbits['normalization']['total_translation_orbits']==2**17
assert orbits['rational_bisections']['translation_orbits']==40917
hist={int(k):v for k,v in orbits['enumeration']['minimum_norm_histogram'].items()}
assert hist=={0:1,4:1218,6:24875,8:63922,10:40917,12:139}
assert sum(hist.values())==2**17

split=read(ART/'det1092_split_descent_v1/independent-replay.json')
assert split['status']=='PASS_INDEPENDENT_HALVING_CYCLE_REPLAY';provenance(split)
summary=[]
for row in split['cases']:
    summary.append({k:row[k] for k in ['family','label','status','steps']})
assert len(summary)==38
assert sum(r['status']=='NEW_INDEPENDENT_DIRECTION' and r['steps']==0 for r in summary)==17
assert sum(r['status']=='NEW_INDEPENDENT_DIRECTION' and r['steps']==1 for r in summary)==8
assert sum(r['status']=='INHERITED_RATIONAL_SPAN' for r in summary)==13
report={'status':'PASS_TRACE_PARITY_DESCENT_IDENTITIES_AND_INPUTS',
    'classification':'new written deductions with symbolic checks and reused independent certificates',
    'universal_norm_identity':True,'constructive_even_trace_direction':True,
    'universal_halving_discriminant_identity':True,'generic_centres':rows,
    'nonzero_trace_parity_labels':2**17-1,'rational_bisection_orbits_with_obstruction':40917,
    'uniform_geometric_conclusion':'For every nonzero M17 parity, the connected halving cover has degree4 and genus9; inherited M17 remains2-saturated under any nonconstant base change of genus below9.',
    'diagnostic_control_decisions_reused':summary,
    'boundary':'No ordinary descended generic class can restrict to the moving odd-trace branch. A separately fitted class may still specialize to a known branch at one parameter. Split-branch independence and member selection are separate.',
    'no_point_searches':True,'no_new_parameters':True,'formal_verification':False,
    'inputs':protocol['inputs'],'protocol_sha256':sha(OUT/'protocol.json'),
    'checker_sha256':sha(Path(__file__))}
dest=OUT/'replay.json'
if dest.exists():assert read(dest)==report
else:
    with dest.open('x') as stream:json.dump(report,stream,indent=2,sort_keys=True);stream.write('\n')
print(report['status'],[(r['label'],r['height']) for r in rows],flush=True)
