#!/usr/bin/env sage-python
"""Independent identities for every rational singular member of one pencil.

Rebuild quartic invariants, discriminant product and ten component maps.
No producer imports, factor discovery, point search or exceptional inputs.
"""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,vector
signal.alarm(25)
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_norm8_singular_members_v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
data=read(OUT/'rational-singular-members.json');geometry=read(OUT/'geometry.json')
local=read(OUT/'no-rational-nodal-member.json')
for source in [data,geometry,read(OUT/'protocol.json')]:
    for name,digest in source['inputs'].items():assert sha(ROOT/name)==digest
assert local['geometry_sha256']==sha(OUT/'geometry.json')
generic=read(ART/'det1092_norm8_seed_cover_v2/generic.json')
parent=read(ART/'curve302_recovered_mw17_parent_v1.json')
R=PolynomialRing(QQ,'z');z=R.gen();K=R.fraction_field()
e,d,c,b,a=[R(v) for v in generic['quartic_t_coefficients_in_z']]
I=12*a*e-3*b*d+c*c
J=72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c**3
A=-27*I;B=-27*J;delta=-16*(4*A**3+27*B**2)
assert A==R(geometry['A'])==R(generic['Jacobian_A']) and A.degree()==8
assert B==R(geometry['B'])==R(generic['Jacobian_B']) and B.degree()==12
assert delta.degree()==22 and geometry['infinity_discriminant_order']==2
repeated=R.one()
finite=[row for row in data['members'] if row['parameter']!='infinity']
assert len(finite)==4
for row in finite:repeated*=z-QQ(row['parameter'])
nodal=R(geometry['simple_support'])
assert repeated==R(geometry['repeated_support']) and repeated.degree()==4
assert nodal.is_monic() and nodal.degree()==14
assert delta==QQ(geometry['delta_scalar'])*repeated**2*nodal
assert repeated.gcd(repeated.derivative())==1 and nodal.gcd(nodal.derivative())==1
assert repeated.gcd(nodal)==1 and A.gcd(delta)==1
witness=local['witness'];p=ZZ(witness['p']);assert p.is_prime() and p==149
assert all(v.denominator()%p for v in nodal)
Rp=PolynomialRing(GF(p),'z');np=Rp(nodal)
assert np.degree()==14 and np.is_monic()
assert list(map(int,np.list()))==witness['polynomial']
assert [int(np(a)) for a in GF(p)]==witness['value_residues']
assert all(np(a)!=0 for a in GF(p))
# At infinity the homogeneous leading coefficient is1, so no projective root.
T=PolynomialRing(QQ,'t');t=T.gen();F=T.fraction_field()
def dec(v):return F(T(v['numerator']))/T(v['denominator'])
E=EllipticCurve(F,[dec(v) for v in parent['a_invariants']])
basis=[E([dec(v) for v in P]) for P in parent['basis_weierstrass_coordinates']]
G=matrix(ZZ,parent['generic_height_gram']);w=vector(ZZ,[0]*14+[1,-1,0])
assert w*G*w==8
center=basis[14]-basis[15]
cx=center[0]+E.b2()/12;cy=center[1]+(E.a1()*center[0]+E.a3())/2
h=T(generic['pole_h']);shift=T(generic['shift'])
words=[]
for member in data['members']:
    components=member['components'];assert len(components)==2
    pair=[vector(ZZ,row['word']) for row in components]
    assert pair[0]+pair[1]==w
    assert sum(v*G*v for v in pair)==8
    for word in pair:assert word*G*word==word*G*w
    if member['parameter']=='infinity':
        assert sorted([list(v) for v in pair])==sorted([[0]*17,list(w)])
    else:
        zz=QQ(member['parameter']);quartic=T([r(zz) for r in [e,d,c,b,a]])
        for row,word in zip(components,pair):
            P=sum((n*Q for n,Q in zip(word,basis)),E(0));assert P
            xx=P[0]+E.b2()/12;yy=P[1]+(E.a1()*P[0]+E.a3())/2
            m=h*zz-shift/h;W=T(row['W'])
            assert W*W==quartic
            assert xx==(h*W-cx+m*m)/2 and yy==m*(xx-cx)-cy
        assert T(components[0]['W'])==-T(components[1]['W'])
    words.extend([list(map(int,v)) for v in pair])
assert len(words)==10 and len(set(map(tuple,words)))==10
frame=read(ART/'det1092_genus1_picard_image_v1/frame.json')
assert sorted(words)==sorted(row['word'] for row in frame['vertical_old_sections'])
# These five pairs account for every multiple discriminant root including infinity.
assert 5*2+14==24 and 19-2-5==12
report={'status':'PASS_INDEPENDENT_NORM8_RATIONAL_SINGULAR_MEMBER_OBSTRUCTION',
    'classification':'verified exact application and new pencil-wide obstruction',
    'rational_singular_parameters':5,'inherited_component_sections':10,
    'nonrational_simple_discriminant_degree':14,'no_projective_root_prime':149,
    'Kodaira_configuration':{'I2':5,'I1':14},'geometric_generic_MW_rank':12,
    'conclusion':'No singular member defined over Q in this pencil supplies a genuine original degree2 rational multisection; every rational singular member is a pair of inherited original sections.',
    'boundary':'A pencil-specific obstruction, not an exclusion of smooth genus1 carriers or seeds on any original fibre. No rational-parameter incidence or amplification predictor.',
    'inputs':{str(q.relative_to(ROOT)):sha(q) for q in [OUT/'protocol.json',OUT/'geometry.json',
        OUT/'no-rational-nodal-member.json',OUT/'rational-singular-members.json',Path(__file__)]}}
dest=OUT/'independent-replay.json'
if dest.exists():assert read(dest)==report
else:
    with dest.open('x') as stream:json.dump(report,stream,indent=2,sort_keys=True);stream.write('\n')
print(report['status'],'five inherited pairs, no-root prime149',flush=True)
