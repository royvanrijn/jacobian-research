#!/usr/bin/env sage-python
"""Independent covariant/group-law comparison on32 pre-existing controls.

The universal proof is frozen first. This is retrospective evaluation only:
10 old split branches,5 old generic points,9 blinded generic reconstructions,
and8 old null-fibre generic points. No point search or new specialization.
"""
import hashlib,json
from pathlib import Path
from sage.all import QQ,PolynomialRing,EllipticCurve

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_reciprocal_quartic_v1'
OLD=ART/'det1092_split_descent_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
proof=read(OUT/'generic.json')
assert proof['status']=='PASS_UNIVERSAL_COVER_CLASS_AND_GENERIC_RECIPROCAL_IDENTITY'
assert proof['checker_sha256']==sha(ROOT/'elliptic-curves/cas/verify_det1092_reciprocal_quartic.sage')
incidence=read(OLD/'incidence-replay.json');descent=read(OLD/'independent-replay.json')
assert incidence['descent_replay_sha256']==sha(OLD/'independent-replay.json')
assert descent['status']=='PASS_INDEPENDENT_HALVING_CYCLE_REPLAY'
for path,expected in descent['inputs'].items():assert sha(ROOT/path)==expected
protocol=read(OLD/'protocol.json')
conic=read(ART/'det1092_orbit8044_rank18_base_change_v2.json')
inputs={str(p.relative_to(ROOT)):sha(p) for p in
    [OUT/'generic.json',OLD/'protocol.json',OLD/'incidence-replay.json',OLD/'independent-replay.json',Path(__file__)]}
R=PolynomialRing(QQ,'m');m=R.gen();rows=[]


def check(E,Z,P,label,family,status,source=None,core_status=None):
    assert E.a1()==E.a2()==E.a3()==0
    assert not Z.is_zero() and not P.is_zero()
    a,b=E.a4(),E.a6();c,d=Z[:2]
    assert P[0]!=c
    slope=(P[1]+d)/(P[0]-c)
    ordinate=2*P[0]+c-slope*slope
    f=m**4-6*c*m*m-8*d*m-3*c*c-4*a
    assert ordinate*ordinate==f(slope) and ordinate
    # Rebuild Hessian coefficients from the quartic, not the producer's G.
    e,D,C,B,A=f.list()
    g0=B*B/16-A*C/6;g1=B*C/12-A*D/2
    g2=C*C/12-B*D/8-A*e;g3=C*D/12-B*e/2;g4=D*D/16-C*e/6
    g=g0*m**4+g1*m**3+g2*m*m+g3*m+g4
    # Independently use homogeneous derivatives for the sextic covariant.
    f_y=4*f-m*f.derivative();g_y=4*g-m*g.derivative()
    hh=(f.derivative()*g_y-f_y*g.derivative())/8
    image=E([g(slope)/(ordinate*ordinate),hh(slope)/(2*ordinate**3)])
    assert image==Z-2*P
    assert image-Z==2*(-P)  # exact half proves the common Kummer class
    conjugate=E([( -ordinate-c+slope*slope)/2,
                  slope*((-ordinate-c+slope*slope)/2-c)-d])
    assert P+conjugate==Z
    row={'family':family,'label':label,'original_classification':status,
         'covering_class':'delta(inherited trace Z)',
         'exact_equalities':['Phi(P)=Z-2P','Phi(P)-Z=2*(-P)','P+conjugate(P)=Z'],
         'Sha_image':'zero','raw_quartic_Jacobian':'original E',
         'image_coordinates':[str(v) for v in image[:2]]}
    if source is not None:
        row['source']=str(source.relative_to(ROOT));row['source_sha256']=sha(source)
    if core_status is not None:row['core16_classification']=core_status
    rows.append(row)


for family in ['302',*protocol['parameters']]:
    folder=OLD/family;frame=read(folder/'generic-frame.json')
    E=EllipticCurve(QQ,list(map(QQ,frame['curve'])))
    basis=[E(list(map(QQ,p))) for p in frame['basis']]
    Z=basis[14]-basis[15] if family=='302' else sum((int(n)*P for n,P in zip(conic['lift']['trace_word'],basis)),E(0))
    inputs[str((folder/'generic-frame.json').relative_to(ROOT))]=sha(folder/'generic-frame.json')
    for i in range(2):
        source=folder/('branch-%d.json'%i);saved=read(source)
        check(E,Z,E(list(map(QQ,saved['original']))),'branch-%d'%i,family,saved['status'],source)
    check(E,Z,basis[0],'generic-section-0',family,'INHERITED_RATIONAL_SPAN')
    if family=='302':
        for arm in protocol['blind_arms']:
            source=folder/(arm+'-full17.json');saved=read(source)
            core=read(folder/(arm+'-core16.json'))
            assert core['original']==saved['original']
            check(E,Z,E(list(map(QQ,saved['original']))),arm,family,saved['status'],source,core['status'])

# Unchanged null controls use generic functions only, with no exceptional point.
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
ROSTER=ART/'det1092_rr_generic_point_controls_v2/protocol.json'
parent=read(PARENT);cases=read(ROSTER)['cases']
inputs.update({str(p.relative_to(ROOT)):sha(p) for p in [PARENT,ROSTER]})
def at(v,t):return R(v['numerator'])(t)/R(v['denominator'])(t)
for case in cases:
    if case['parameter']=='0':continue
    t=QQ(case['parameter']);old=EllipticCurve(QQ,[at(v,t) for v in parent['a_invariants']])
    E=EllipticCurve(QQ,[-old.c4()/48,-old.c6()/864])
    def section(i):
        xx,yy=[at(v,t) for v in parent['basis_weierstrass_coordinates'][i]]
        return E([xx+old.b2()/12,yy+(old.a1()*xx+old.a3())/2])
    check(E,section(14)-section(15),section(0),'generic-section-0',case['label'],'INHERITED_RATIONAL_SPAN')

assert len(rows)==32
assert sum(r['original_classification']=='NEW_INDEPENDENT_DIRECTION' for r in rows)==8
assert sum(r.get('core16_classification')=='NEW_INDEPENDENT_DIRECTION' for r in rows)==9
result={'classification':'independent retrospective verification, not a selection experiment',
        'status':'PASS_COVARIANT_CLASS_CONTROLS','count':len(rows),
        'independent_over_full17':8,'inherited_over_full17':24,
        'old_blinded_core16_gains':9,'rows':rows,'inputs':inputs,
        'boundary':'All covering classes are inherited and have zero Sha image despite different marked-point independence. No claim that the fixed-z carrier Jacobian equals the original fibre.'}
path=OUT/'control-replay.json'
if path.exists():assert read(path)==result
else:
    with path.open('x') as stream:json.dump(result,stream,indent=2,sort_keys=True);stream.write('\n')
print(result['status'],'controls',len(rows),'new',result['independent_over_full17'],'inherited',result['inherited_over_full17'])
