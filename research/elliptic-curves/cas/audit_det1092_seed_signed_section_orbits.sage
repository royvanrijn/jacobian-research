#!/usr/bin/env sage-python
"""All34 signed displayed generic-section curves on the fixed seed carrier.

This is retrospective: it uses the already fixed first-seed carrier label,
not later points. Fibrewise translations preserve that label. No translated
orbit is enumerated, and no elliptic point search is run. Cap25seconds.
"""
import hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector,gcd,lcm
from rational_root_lattice_certificate import root_certificate
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_signed_source_orbits_v1';PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
PENCIL=ART/'det1092_norm8_seed_cover_v2/generic.json'
CARRIER=ART/'det1092_norm8_seed_cover_v2/equation-only-cover-01.json'
OLDPROTO=ART/'det1092_pencil_multiples_v2/protocol.json'
HELPER=Path(__file__).with_name('rational_root_lattice_certificate.py')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    path=OUT/name
    if path.exists():assert json.loads(path.read_text())==data
    else:
        with path.open('x') as stream:json.dump(data,stream,indent=2,sort_keys=True);stream.write('\n')
def rec(a):return {'numerator':list(map(str,a.numerator().list())),'denominator':list(map(str,a.denominator().list()))}
d=json.loads(PARENT.read_text());p=json.loads(PENCIL.read_text());carrier=json.loads(CARRIER.read_text())
zstar=QQ(carrier['cover_parameter_z']);primes=json.loads(OLDPROTO.read_text())['ramification_certificate_primes']
OUT.mkdir(exist_ok=True)
save('protocol.json',{'classification':'retrospective source-curve incidence obstruction',
    'sources':'O and both signs of each of the17 displayed generic sections',
    'carrier_label':str(zstar),'limits':{'seconds':25,'signed_curves':34,
        'prime_pool':primes,'orbit_enumeration':0,'point_searches':0,'later_points':0},
    'inputs':{str(path.relative_to(ROOT)):sha(path) for path in [PARENT,PENCIL,CARRIER,OLDPROTO,HELPER,Path(__file__)]},
    'implication':'A signed source curve with no rational point on this carrier cannot supply one after any rational fibrewise translation. Degree-one curves specialize inside the full Picard image.'})
R=PolynomialRing(QQ,'u');u=R.gen();K=R.fraction_field()
def dec(v):return K(R(v['numerator']))/R(v['denominator'])
E=EllipticCurve(K,[dec(a) for a in d['a_invariants']]);basis=[E([dec(a) for a in row]) for row in d['basis_weierstrass_coordinates']]
G=matrix(ZZ,d['generic_height_gram']);w=vector(ZZ,[0]*14+[1,-1,0]);C=basis[14]-basis[15]
cx=C[0]+E.b2()/12;cy=C[1]+(E.a1()*C[0]+E.a3())/2
h=R(p['pole_h']);shift=R(p['shift']);cases=[]
for i,base in enumerate(basis):
    for sign in [-1,1]:
        P=sign*base;xx=P[0]+E.b2()/12;yy=P[1]+(E.a1()*P[0]+E.a3())/2
        slope=(yy+cy)/(xx-cx);zmap=(h*slope+shift)/(h*h)
        degree=int(G[i,i]-sign*(G*w)[i]);assert max(zmap.numerator().degree(),zmap.denominator().degree())==degree
        row={'index':i,'sign':sign,'degree':degree,'z_map':rec(zmap)}
        if degree==0:
            assert zmap!=zstar;row['status']='VERTICAL_AWAY_FROM_CARRIER'
        elif degree==1:row['status']='INHERITED_ALTERNATE_SECTION'
        else:
            f=R(zmap.numerator()-zstar*zmap.denominator());infinity=f.degree()<degree
            f*=lcm([a.denominator() for a in f.list()]);f/=gcd([ZZ(a) for a in f.list()]);f=R(f)
            if f.leading_coefficient()<0:f=-f
            roots=[];g=f
            if not f[0]:
                roots.append('0')
                while not g[0]:g=R(g/u)
            proof=root_certificate(g,primes) if g.degree()>0 else None
            if proof:roots+=proof['rational_roots']
            row.update(primitive_polynomial=list(map(str,f.list())),rational_preimages=roots,
                       infinity_preimage=infinity,root_proof=proof,
                       status='RATIONAL_SOURCE_POINT_REQUIRES_CLASS_TEST' if roots or infinity else 'NO_RATIONAL_SOURCE_POINT')
        save('source-%02d-%s.json'%(i,'minus' if sign<0 else 'plus'),row);cases.append(row)
        print(i,sign,degree,row['status'],flush=True)
complete=all(r['status']!='RATIONAL_SOURCE_POINT_REQUIRES_CLASS_TEST' for r in cases)
save('summary.json',{'status':'PASS_ENTIRE_FIXED_SOURCE_TRANSLATION_ORBITS_EXCLUDED' if complete else 'PARTIAL_CLASS_TEST_NEEDED',
    'classification':'retrospective theorem application, not a prospective source selector',
    'counts':{key:sum(r['status']==key for r in cases) for key in sorted({r['status'] for r in cases})},
    'source_count':34,'zero_section':'vertical and disjoint from the smooth seed carrier',
    'first_seed_outside_full_picard_image_dependency':'det1092_genus1_picard_image_v1/replay.json',
    'boundary':'Only the signed displayed source curves and the group preserving this alternate fibration are covered. Not all original MW sections, alternating fibration words, or all representatives of the same original residual direction are excluded.',
    'checker_sha256':sha(Path(__file__))})
