#!/usr/bin/env sage-python
"""Exact symbolic intersection calculation, plus generic-parent applications.

No new curve construction, specialization, point search, or lattice census.
Uses independent unramified-map replays as dependencies, not numerical heights.
"""
import hashlib,json
from pathlib import Path
from sage.all import QQ, PolynomialRing, matrix

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_pencil_multiples_v2'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

R=PolynomialRing(QQ,['d','a','vnorm','wnorm','wv','g'])
d,a,vnorm,wnorm,wv,g=R.gens()
# Original NS = <O,F> + (-MW), O^2=-2, O.F=1, F^2=0.
# C=d O+a F+phi(w); S_v=O+(vnorm/2)F+phi(v).
C_O=a-2*d
S_O=vnorm/2-2
C_S=-2*d+d*vnorm/2+a-wv
chi=2*d
new_height=2*chi+2*C_O
cross=chi+C_O+d*S_O-C_S
C_square=-2*d*d+2*d*a-wnorm
assert new_height==2*a and cross==wv
K=R.fraction_field()
schur=K(new_height)-K(wnorm)/d
assert schur==2*d+K(C_square)/d
adjunction_schur=2*d+K(2*g-2)/d
assert schur.subs(wnorm=2*d*a-2*d*d-(2*g-2))==adjunction_schur
parent=json.loads((ART/'curve302_recovered_mw17_parent_v1.json').read_text())
G=matrix(QQ,parent['generic_height_gram'])
assert G.nrows()==17 and G.det()==1092 and G.is_positive_definite()
rows=[]
for n in [2,3]:
    replay=OUT/('replay%d.json'%n)
    r=json.loads(replay.read_text())
    assert r['checker_sha256']==sha(ROOT/'elliptic-curves/cas/verify_det1092_pencil_multiples.sage')
    assert r['construction_sha256']==sha(OUT/('multiple%d.json'%n))
    assert r['status']=='PASS_INDEPENDENT_RANK18_RATIONAL_BASE_CHANGE'
    dd=QQ(r['degree']); assert r['smooth_ramification_points']==2*dd-2
    perpendicular=2*dd-2/dd
    rows.append({'n':n,'degree':int(dd),'genus':0,
                 'schur':str(perpendicular),
                 'displayed_rank18_determinant':str(dd**17*G.det()*perpendicular),
                 'old_basis_min_diagonal_scaled':str(dd*min(G.diagonal())),
                 'replay_sha256':sha(replay)})
assert [r['schur'] for r in rows]==['399/10','3363/29']
assert QQ(adjunction_schur.subs(d=2,g=0))==3
assert QQ(adjunction_schur.subs(d=2,g=1))==4
assert QQ(adjunction_schur.subs(d=2,g=2))==5
result={'classification':'new deduction from established intersection-height theory',
        'status':'PASS_MULTISECTION_SCHUR_IDENTITY',
        'hypotheses':['rootless elliptic K3 with full NS=U+(-MW)',
                      'smooth multisection of degree d and genus g',
                      'cover unramified over every singular elliptic fibre'],
        'identity':'Schur=2*d+(2*g-2)/d',
        'applications':rows,
        'degree2_regressions':{'genus0':'3','genus1':'4','genus2':'5'},
        'boundary':'Generic heights in the pulled-back subgroup, not specialized height values or amplification predictions.',
        'checker_sha256':sha(Path(__file__))}
path=OUT/'height-identity.json'
if path.exists():assert json.loads(path.read_text())==result
else:
    with path.open('x') as stream:json.dump(result,stream,indent=2,sort_keys=True);stream.write('\n')
print(json.dumps(result,sort_keys=True))
