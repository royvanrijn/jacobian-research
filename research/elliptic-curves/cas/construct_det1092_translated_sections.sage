#!/usr/bin/env sage-python
"""Six generic curves: minimal-degree alternate multisections, translated +/-B.

Only equation/generic-section inputs. Each --case has a25-second external
cap and freezes all six cases before computation. No target evaluation.
"""
import argparse,hashlib,json,time
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,vector
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_translated_sections_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
PENCIL=ART/'det1092_norm8_seed_cover_v2/generic.json'
RESTRICTIONS=ART/'det1092_genus1_picard_image_v1/generic-restrictions.json'
OLDPROTO=ART/'det1092_pencil_multiples_v2/protocol.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    p=OUT/name
    if p.exists():assert json.loads(p.read_text())==data
    else:
        with p.open('x') as stream:json.dump(data,stream,indent=2,sort_keys=True);stream.write('\n')
def rec(a):return {'numerator':list(map(str,a.numerator().list())),'denominator':list(map(str,a.denominator().list()))}
parser=argparse.ArgumentParser();parser.add_argument('--case',type=int,required=True);args=parser.parse_args()
parent=json.loads(PARENT.read_text());pencil=json.loads(PENCIL.read_text());restrictions=json.loads(RESTRICTIONS.read_text())
primes=json.loads(OLDPROTO.read_text())['ramification_certificate_primes']
G=matrix(ZZ,parent['generic_height_gram']);w=vector(ZZ,[0]*14+[1,-1,0])
degrees=[int(G[i,i]-(G*w)[i]) for i in range(17)]
assert degrees==[row['intersection_degree'] for row in restrictions['sections']]
minimum=min(d for d in degrees if d>1);indices=[i for i,d in enumerate(degrees) if d==minimum]
cases=[{'section':i,'translation':sign} for i in indices for sign in [-1,1]]
assert minimum==2 and indices==[1,7,9] and len(cases)==6 and 0<=args.case<6
OUT.mkdir(exist_ok=True)
save('protocol.json',{'classification':'prospective generic-only rational-curve construction',
    'selection':'All original displayed generic sections with smallest intersection degree greater than one against the alternate fibre; translate once by each sign of B. No target point or parameter enters.',
    'intersection_degrees':degrees,'selected_degree':minimum,'cases':cases,
    'limits':{'seconds_per_case':25,'curves':6,'translations_per_curve':1,'point_searches':0,
              'exceptional_point_inputs':0,'parameter_evaluations':0},
    'certificate_prime_pool':primes,
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [PARENT,PENCIL,RESTRICTIONS,OLDPROTO,Path(__file__)]}})
case=cases[args.case];R=PolynomialRing(QQ,'u');u=R.gen();L=R.fraction_field()
def decode(row):return L(R(row['numerator']))/R(row['denominator'])
source=restrictions['sections'][case['section']];z=decode(source['z']);Wsource=decode(source['W'])
assert max(z.numerator().degree(),z.denominator().degree())==2
origin=pencil['sections'][0]
def compose(row):return L(R(row['numerator'])(z))/R(row['denominator'])(z)
t0=compose(origin['t_of_z']);s=compose(origin['W_of_z']);S=PolynomialRing(L,'t');t=S.gen()
f=S([R(row)(z) for row in pencil['quartic_t_coefficients_in_z']]);assert Wsource**2==f(u)
q0,q1,q2,q3,q4=f(t+t0).list();assert q0==s*s and s
J=EllipticCurve(L,[0,q2,0,q1*q3-4*s*s*q4,s*s*q3*q3+q1*q1*q4-4*s*s*q2*q4])
bx=q1*q1/(4*s*s)-q2;by=-(q1*bx+2*s*s*q3)/(2*s);B=J([bx,by])
du=u-t0;assert du
x=(2*s*(Wsource+s)+q1*du)/(du*du)
y=((x*x-4*s*s*q4)*du-q1*x-2*s*s*q3)/(2*s);Q=J([x,y])
target=Q+case['translation']*B;assert not target.is_zero();xx,yy=target[:2]
v=(2*s*yy+q1*xx+2*s*s*q3)/(xx*xx-4*s*s*q4)
T=t0+v;W=xx*v*v/(2*s)-s-q1*v/(2*s)
assert W*W==f(T)
degree=max(T.numerator().degree(),T.denominator().degree());assert degree>1
save('map-%02d.json'%args.case,{'status':'EXACT_GENERIC_TRANSLATION_MAP',
    'classification':'generic-only construction; independent rank replay required',
    'case':case,'source_parameter':'u on the original generic section',
    'alternate_parameter':rec(z),'source_W':rec(Wsource),
    't_of_u':rec(T),'W_of_u':rec(W),'degree':int(degree),
    'map_to_original':pencil['inverse_map'],'checker_sha256':sha(Path(__file__))})
print(args.case,case,'exact map degree',degree,flush=True)
# Reuse the unchanged small prime pool only to certify smooth ramification.
E=EllipticCurve(L,[decode(a) for a in parent['a_invariants']]);delta=E.discriminant()
assert delta.denominator().degree()==0;delta=delta.numerator()
N,D=T.numerator(),T.denominator();ram=N.derivative()*D-N*D.derivative()
attempts=[];witness=None
for p in primes:
    Rp=PolynomialRing(GF(p),'u')
    try:nn,dd,rr,ds=map(Rp,[N,D,ram,delta])
    except (ValueError,ZeroDivisionError):attempts.append({'p':p,'gate':'denominator'});continue
    if any(a.degree()!=b.degree() for a,b in [(nn,N),(dd,D),(rr,ram),(ds,delta)]):
        attempts.append({'p':p,'gate':'degree_drop'});continue
    if nn.gcd(dd).degree():attempts.append({'p':p,'gate':'map_bad'});continue
    bad=dd*sum((ds[i]*nn**i*dd**(ds.degree()-i) for i in range(ds.degree()+1)),Rp.zero())
    radical=rr.squarefree_part();common=radical.gcd(bad);lower=radical.degree()-common.degree()
    attempts.append({'p':p,'gate':'checked','smooth_ramification_lower_bound':int(lower)})
    if lower>0:
        witness={'p':p,'radical_mod_p':list(map(int,radical.list())),
                 'bad_gcd_mod_p':list(map(int,common.list())),
                 'smooth_ramification_lower_bound':int(lower)};break
save('rank-%02d.json'%args.case,{'status':'GENERIC_RANK_AT_LEAST18' if witness else 'UNKNOWN',
    'classification':'verified ramification application, pending independent replay',
    'case':case,'degree':int(degree),'attempts':attempts,'witness':witness,
    'proof':'Translation by a generic section is an automorphism of the minimal elliptic K3, so the image curve is smooth rational. A multisection ramified above a smooth original fibre supplies an independent section after base change.',
    'boundary':'No specialization or independence at a particular rational parameter has been tested.'})
print('ramification',None if witness is None else witness['p'],flush=True)
