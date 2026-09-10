#!/usr/bin/env sage-python
"""Explicit native twist sections and height-eight certificates on14 fixed covers."""
import argparse,hashlib,json,runpy
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,prime_range
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
INPUT=ART/'curve302_class1_fibre_carriers_v1.json'
TRACE=ART/'x1092_class1_realization_trace_v1.json'
OUTPUT=ART/'curve302_class1_twist_sections_v1.json'
WORK=ROOT/'artifacts/local/elliptic-curves/curve302-class1-twist-sections-v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):
    if p.exists():assert read(p)==d,p
    else:
        p.parent.mkdir(parents=True,exist_ok=True)
        with p.open('x')as out:json.dump(d,out,indent=2,sort_keys=True);out.write('\n')
def compute(check=False):
    checker=ROOT/'elliptic-curves/cas/verify_curve302_class1_bridge.sage'
    bindings={str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),INPUT,TRACE,checker]}
    if not check:save(WORK/'protocol.json',{'bindings':bindings,'limits':{'wall_seconds':120,'workers':1,'cases':14,'good_prime_bound':2000,'new_parameters':0,'section_searches':0}})
    data,tr=map(read,[INPUT,TRACE]);R=PolynomialRing(QQ,'t');F=R.fraction_field()
    A,B=R(tr['A']),R(tr['B']);delta=4*A**3+27*B**2
    assert A.degree()==8 and B.degree()==12 and delta.degree()==24
    entries=[];polys=[]
    for row in data['carriers']:
        q=R(row['branch_quartic']);parse=lambda c:F(R(c['numerator']))/R(c['denominator'])
        x0,x1,y0,y1=[parse(row['source_maps'][k])for k in ['x0','x1','y0','y1']]
        r=y0*y0/(x1*x1*q)-2*x0
        y=y0*(x0-r)/(x1*q)-y1
        X=R(q*r);Y=R(q*q*y)
        assert X.degree()==8 and Y.degree()==12
        assert Y*Y==X**3+q*q*A*X+q**3*B
        rec={'id':row['id'],'twist_convention':'Y^2=X^3+q(t)^2*A(t)*X+q(t)^3*B(t)',
             'X':list(map(str,X.list())),'Y':list(map(str,Y.list())),
             'degree_X':8,'degree_Y':12,'section_height':8,
             'geometric_section_height_lower_bound':4,'native_section_primitive':True,
             'geometric_twist_rank_interval':[1,22],'exact_twist_rank':'UNKNOWN'}
        entries.append(rec);polys.append((q,X))
        if not check:save(WORK/(row['id']+'.json'),rec)
    witness=None
    for pp in prime_range(3,2001):
        p=int(pp)
        if any(c.denominator()%p==0 for poly in [delta,*[v for q,X in polys for v in [q,X]]]for c in poly.list()):continue
        D=delta.change_ring(GF(p))
        if D.degree()!=24 or D.gcd(D.derivative()).degree()!=0:continue
        rows=[]
        for q,X in polys:
            Q,V=q.change_ring(GF(p)),X.change_ring(GF(p))
            if Q.degree()!=4 or V.degree()!=8:break
            bits=[int(Q.resultant(Q.derivative())),int(Q.resultant(D)),int(Q.resultant(V))]
            if not all(bits):break
            rows.append(bits)
        if len(rows)==14:
            witness={'prime':p,'delta_squarefree_resultant':int(D.resultant(D.derivative())),
                     'q_squarefree_q_delta_q_X_resultants':rows};break
    assert witness is not None
    # Replay the modular proof by independent Sylvester determinants.
    if check:
        sylvester=runpy.run_path(str(checker))['sylvester_resultant']
        p=witness['prime'];assert ZZ(p).is_prime(proof=True);K=GF(p);D=delta.change_ring(K)
        assert int(sylvester(D,D.derivative()))==witness['delta_squarefree_resultant']!=0
        for (q,X),expected in zip(polys,witness['q_squarefree_q_delta_q_X_resultants']):
            Q,V=q.change_ring(K),X.change_ring(K)
            assert [int(sylvester(Q,Q.derivative())),int(sylvester(Q,D)),int(sylvester(Q,V))]==expected
    result={'schema':'curve302.class1-twist-sections.v1','bindings':bindings,'sections':entries,
            'finite_certificate':witness,'status':'PASS_FOURTEEN_PRIMITIVE_HEIGHT_EIGHT_SECTIONS',
            'geometric_fibres':'24 I1 plus4 I0*; smooth at infinity','chi':4,'h11':40,
            'trivial_lattice_rank':18,'geometric_MW_height_lower_bound':4,
            'boundary':'An explicit primitive rank-one subgroup in each twist, not an exact rank-one theorem. Rank at least2 and a second independent section remain UNKNOWN.'}
    print('PASS14 polynomial twist sections of height8; all nonzero sections have height at least4; prime',witness['prime'],flush=True)
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args()
    result=compute(args.check)
    if args.check:assert read(OUTPUT)==result
    else:save(OUTPUT,result)
