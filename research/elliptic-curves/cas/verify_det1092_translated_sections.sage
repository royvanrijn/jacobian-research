#!/usr/bin/env sage-python
"""Manual chord-law replay, ramification proof and independent302 root exclusion.

One --case per25-second process. No factorization or elliptic point search.
"""
import argparse,hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,gcd,lcm
from rational_root_lattice_certificate import root_certificate
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_translated_sections_v1';PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
PENCIL=ART/'det1092_norm8_seed_cover_v2/generic.json'
HELPER=Path(__file__).with_name('rational_root_lattice_certificate.py')
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    path=OUT/name
    if path.exists():assert read(path)==data
    else:
        with path.open('x') as stream:json.dump(data,stream,indent=2,sort_keys=True);stream.write('\n')
parser=argparse.ArgumentParser();parser.add_argument('--case',type=int,required=True);arg=parser.parse_args();assert 0<=arg.case<6
protocol=read(OUT/'protocol.json');d=read(PARENT);p=read(PENCIL)
for path,expected in protocol['inputs'].items():assert sha(ROOT/path)==expected
mp=OUT/('map-%02d.json'%arg.case);rp=OUT/('rank-%02d.json'%arg.case);ip=OUT/('302-incidence-%02d.json'%arg.case)
saved=read(mp);rank=read(rp);inc=read(ip);i=saved['case']['section'];sign=saved['case']['translation']
assert saved['case']==protocol['cases'][arg.case]
R=PolynomialRing(QQ,'u');u=R.gen();L=R.fraction_field()
def dec(a):return L(R(a['numerator']))/R(a['denominator'])
aa=[dec(a) for a in d['a_invariants']];a1,a2,a3,a4,a6=aa;b2=a1*a1+4*a2
xold,yold=[dec(a) for a in d['basis_weierstrass_coordinates'][i]]
assert yold*yold+a1*xold*yold+a3*yold==xold**3+a2*xold*xold+a4*xold+a6
X=xold+b2/12;Y=yold+(a1*xold+a3)/2;h=R(p['pole_h']);shift=R(p['shift'])
cx=L(R(p['nx']))/(h*h);cy=L(R(p['ny']))/(h**3);m=(Y+cy)/(X-cx)
z=(h*m+shift)/(h*h);Wsource=(2*X+cx-m*m)/h
assert z==dec(saved['alternate_parameter']) and Wsource==dec(saved['source_W'])
origin=p['sections'][0]
def compose(row):return L(R(row['numerator'])(z))/R(row['denominator'])(z)
t0=compose(origin['t_of_z']);s=compose(origin['W_of_z']);S=PolynomialRing(L,'t');t=S.gen()
f=S([R(row)(z) for row in p['quartic_t_coefficients_in_z']]);assert f(u)==Wsource**2
q0,q1,q2,q3,q4=f(t+t0).list();assert q0==s*s
a4j=q1*q3-4*s*s*q4;a6j=s*s*q3*q3+q1*q1*q4-4*s*s*q2*q4
du=u-t0;x=(2*s*(Wsource+s)+q1*du)/(du*du)
y=((x*x-4*s*s*q4)*du-q1*x-2*s*s*q3)/(2*s)
bx=q1*q1/(4*s*s)-q2;by=-(q1*bx+2*s*s*q3)/(2*s)
assert y*y==x**3+q2*x*x+a4j*x+a6j and by*by==bx**3+q2*bx*bx+a4j*bx+a6j
slope=(sign*by-y)/(bx-x);xn=slope*slope-q2-x-bx;yn=-y+slope*(x-xn)
assert yn*yn==xn**3+q2*xn*xn+a4j*xn+a6j
v=(2*s*yn+q1*xn+2*s*s*q3)/(xn*xn-4*s*s*q4)
T=t0+v;W=xn*v*v/(2*s)-s-q1*v/(2*s)
assert T==dec(saved['t_of_u']) and W==dec(saved['W_of_u']) and W*W==f(T)
# Inverse chord restores Q: a direct birationality check on this source.
back_slope=(-sign*by-yn)/(bx-xn);back_x=back_slope**2-q2-xn-bx
back_y=-yn+back_slope*(xn-back_x);assert back_x==x and back_y==y
N,D=T.numerator(),T.denominator();degree=max(N.degree(),D.degree());assert degree==saved['degree'] and N.gcd(D).degree()==0
delta=EllipticCurve(L,aa).discriminant();assert delta.denominator().degree()==0;delta=delta.numerator()
critical=N.derivative()*D-N*D.derivative();cert=rank['witness'];prime=cert['p']
Rp=PolynomialRing(GF(prime),'u');nn,dd,rr,ds=map(Rp,[N,D,critical,delta])
assert all(a.degree()==b.degree() for a,b in [(nn,N),(dd,D),(rr,critical),(ds,delta)]) and nn.gcd(dd).degree()==0
radical=Rp(cert['radical_mod_p']);assert radical.gcd(radical.derivative()).degree()==0 and rr%radical==0
assert pow(radical,int(rr.degree()),rr)==0
bad=dd*sum((ds[k]*nn**k*dd**(ds.degree()-k) for k in range(ds.degree()+1)),Rp.zero())
common=radical.gcd(bad);assert common.monic()==Rp(cert['bad_gcd_mod_p']).monic()
assert radical.degree()-common.degree()==cert['smooth_ramification_lower_bound']>0
# Exact rational-root exclusion at the original target t=0, independently of QQ factorization.
poly=N*lcm([c.denominator() for c in N.list()]);poly/=gcd([ZZ(c) for c in poly.list()]);poly=R(poly)
if poly.leading_coefficient()<0:poly=-poly
assert poly==R(inc['primitive_polynomial']) and poly.degree()==degree and poly[0]
root_proof=root_certificate(poly,protocol['certificate_prime_pool']);assert not root_proof['rational_roots']
report={'status':'PASS_GENERIC_RANK18_AND_EXACT302_EXCLUSION',
    'classification':'independent generic construction and retrospective incidence replay',
    'case':saved['case'],'degree':int(degree),'ramification_prime':prime,
    'smooth_ramification_lower_bound':cert['smooth_ramification_lower_bound'],
    'chord_and_inverse_identities':True,'root_proof':root_proof,
    'inputs':{str(path.relative_to(ROOT)):sha(path) for path in [PARENT,PENCIL,OUT/'protocol.json',mp,rp,ip,HELPER,Path(__file__)]},
    'boundary':'Rank18 over this rational base change, not on302. No other specialization or residual representative is classified.'}
save('independent-%02d.json'%arg.case,report)
print(arg.case,saved['case'],'degree',degree,'root prime',root_proof['prime'],report['status'],flush=True)
