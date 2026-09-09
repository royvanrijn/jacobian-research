#!/usr/bin/env sage-python
"""Independent exact unit-intersection proof and separate sealed-control comparison.

Manual reversed Q(t) additions; no producer imports or rational-point search.
Every action has a25-second cap. Controls are read only by the compare action,
after the generic construction and its independent replay have been sealed.
"""
import argparse,csv,hashlib,itertools,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector,gcd,lcm
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_unit_intersection_obstructions_v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,row):
    p=OUT/name
    if p.exists():assert read(p)==row,'immutable replay changed'
    else:
        with p.open('x') as s:json.dump(row,s,indent=2,sort_keys=True);s.write('\n')
def replay():
    protocol=read(OUT/'protocol.json');summary=read(OUT/'summary.json')
    for name,h in {**protocol['inputs'],**summary['inputs']}.items():assert sha(ROOT/name)==h
    parent=read(ART/'curve302_recovered_mw17_parent_v1.json');G=matrix(ZZ,parent['generic_height_gram'])
    with (ART/'curve302_parent_degree2_multisection_orbits_v1.tsv').open() as stream:
        source=next(r for r in csv.DictReader(stream,delimiter='\t') if int(r['orbit_mask'])==8044)
    assert source['category']=='rational'
    w=vector(ZZ,source['parent_MW17_w'].split());assert w*G*w==10
    selection=read(OUT/'selection.json');assert selection['word']==list(w)==summary['word']
    # Independently generate the same bounded word set by vector addition.
    unit=[vector(ZZ,[int(i==j) for j in range(17)]) for i in range(17)]
    pool=[sgn*u for u in unit for sgn in [-1,1]]
    pool += [a*unit[i]+b*unit[j] for i in range(17) for j in range(i+1,17) for a,b in itertools.product([-1,1],repeat=2)]
    assert len(pool)==578 and len({tuple(v) for v in pool})==578
    assert [r['word'] for r in selection['exposures']]==[list(v) for v in pool]
    expected=[]
    for v,row in zip(pool,selection['exposures']):
        norm,pair=v*G*v,w*G*v;chosen=norm==4 and pair==3
        assert row['norm']==norm and row['pairing']==pair and row['selected']==chosen
        if chosen:expected.append(v)
    assert [list(v) for v in expected]==selection['selected'] and len(expected)==summary['count']<=16
    R=PolynomialRing(QQ,'t');K=R.fraction_field();t=R.gen()
    def dec(v):return K(R(v['numerator']))/R(v['denominator'])
    original=EllipticCurve(K,[dec(v) for v in parent['a_invariants']])
    A,B=R(-original.c4()/48),R(-original.c6()/864)
    basis=[]
    for row in parent['basis_weierstrass_coordinates']:
        x,y=map(dec,row);X,Y=x+original.b2()/12,y+(original.a1()*x+original.a3())/2
        assert Y*Y==X**3+A*X+B;basis.append((X,Y))
    def add(P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P;u,v=Q
        if x==u and y==-v:return None
        slope=(v-y)/(u-x) if x!=u else (3*x*x+A)/(2*y)
        z=slope*slope-x-u
        return z,slope*(x-z)-y
    def word(v):
        total=None
        for n,P in reversed(list(zip(v,basis))):
            for unused in range(abs(n)):total=add(total,(P[0],P[1] if n>0 else -P[1]))
        return total
    conic=read(OUT/'conic.json')
    h,N,V,m,g,b,k,q=[R(conic[key]) for key in ['h','nx','ny','m','g','b','k','q']]
    assert word(-w)==(K(N)/h**2,K(V)/h**3)
    assert h.is_monic() and h.degree()==3 and gcd(h,N)==1
    assert all(f.degree()<=n for f,n in zip([N,V,m,g,b,k,q],[10,15,5,9,4,3,2]))
    assert m*N+V==h*h*g and m*m-N==h*h*b and m*b-2*g==h*h*k
    assert 4*m*k-3*b*b-4*A==h*h*q and q.degree()==2 and gcd(q,q.derivative())==1
    x0,x1,y0,y1=[R(v) for v in conic['maps']]
    assert (x0,x1,y0,y1)==(b/2,h/2,-h*k/2,-m/2)
    assert y0*y0+y1*y1*q==x0**3+3*x0*x1*x1*q+A*x0+B
    assert 2*y0*y1==3*x0*x0*x1+x1**3*q+A*x1
    U,V0=[R(v) for v in conic['bezout']];assert U*h+V0*m==1
    aa,bb,cc,dd=map(QQ,read(ART/'det1092_reduced_parameter_chart_v1/generic-proof.json')['parameter_matrix'])
    verified=[];S=PolynomialRing(ZZ,'s');s=S.gen();D=S(1)
    for index,v in enumerate(expected):
        data=read(OUT/('intersection-%02d.json'%index))
        assert data['index']==index and data['word']==list(v) and data['companion_word']==list(w-v)
        assert v*G*v-w*G*v==data['intersection_number']==1
        X,Y=map(R,word(v));assert all(f.degree()<=n for f,n in [(X,4),(Y,6)])
        assert [R(a) for a in data['section']]==[X,Y]
        lift=R(data['ordinate_lift']);assert lift==U*(2*X-b)-V0*(2*Y+h*k)
        tests=[2*X-b-h*lift,2*Y+h*k+m*lift,lift*lift-q]
        divisor=R(data['intersection_divisor'])
        assert divisor==gcd(tests).monic()
        # This frozen14-case application has no infinity collision. General
        # unit-intersection rationality does not require that restriction.
        assert divisor.degree()==1 and data['original_parameter'] is not None
        tau=QQ(data['original_parameter']);omega=QQ(data['omega'])
        assert divisor(tau)==0 and lift(tau)==omega and omega and omega*omega==q(tau)
        assert all(f(tau)==0 for f in tests)
        assert 4*A(tau)**3+27*B(tau)**2 and data['smooth_fibre']
        assert data['conclusion']=='BOTH_BRANCHES_INHERITED'
        P=[X(tau),Y(tau)];Q=[x0(tau)-x1(tau)*omega,y0(tau)-y1(tau)*omega]
        assert P==list(map(QQ,data['branch'])) and Q==list(map(QQ,data['companion']))
        companion=word(w-v);assert [f(tau) for f in companion]==Q
        assert P[1]**2==P[0]**3+A(tau)*P[0]+B(tau)
        assert Q[1]**2==Q[0]**3+A(tau)*Q[0]+B(tau)
        reduced=QQ(data['reduced_parameter']);assert (aa*reduced+bb)/(cc*reduced+dd)==tau
        D*=reduced.denominator()*s-reduced.numerator()
        verified.append(dict(index=index,reduced_parameter=str(reduced),original_parameter=str(tau),
            word=list(map(int,v)),companion_word=list(map(int,w-v)),intersection_number=1,
            branch_count=2,conclusion='BOTH_BRANCHES_EXACT_GENERIC_WORDS'))
    assert len({v['original_parameter'] for v in verified})==len(verified)==summary['distinct_parameters']
    assert gcd(D,D.derivative())==1 and D.degree()==len(verified)
    polynomial=dict(coefficients=list(map(str,D.list())),degree=int(D.degree()),
        convention='primitive positive-leading integral polynomial in reduced parent parameter s',
        factors=[str(QQ(v['reduced_parameter']).denominator()*s-QQ(v['reduced_parameter']).numerator()) for v in verified],
        conclusion='Every rational zero is smooth and both orbit8044 branches are inherited; outside its zeros UNKNOWN.')
    assert gcd(D.list())==1
    save('collision-polynomial.json',polynomial)
    return dict(status='PASS_INDEPENDENT_UNIT_INTERSECTION_OBSTRUCTIONS',classification='new deduction/application',
        count=len(verified),cases=verified,checker_sha256=sha(Path(__file__)),
        inputs={str(p.relative_to(ROOT)):sha(p) for p in [OUT/'protocol.json',OUT/'selection.json',OUT/'summary.json',OUT/'collision-polynomial.json']},
        scope='Complete for the578-word support<=2 selector, not all norm4 sections or all dependent splits.',
        control_parameters_read=0,point_searches=0)
def compare():
    proof=read(OUT/'independent-replay.json');assert proof['status']=='PASS_INDEPENDENT_UNIT_INTERSECTION_OBSTRUCTIONS'
    for name,h in proof['inputs'].items():assert sha(ROOT/name)==h
    panel=ART/'det1092_rational_bisection_index_v1/controls.json'
    old=ART/'det1092_euclidean_seed_admission_v1/protocol.json'
    roster=[dict(label=r['label'],original_parameter=r['parameter']) for r in read(panel)['rows'][0]['cases']]
    roster += [r for r in read(old)['roster'] if 'reduced_parameter' in r]
    aa,bb,cc,dd=map(QQ,read(ART/'det1092_reduced_parameter_chart_v1/generic-proof.json')['parameter_matrix'])
    predicted={QQ(v['original_parameter']):v for v in proof['cases']}
    tests=[]
    for row in roster:
        if 'original_parameter' in row:tau=QQ(row['original_parameter'])
        else:
            s=QQ(row['reduced_parameter']);tau=(aa*s+bb)/(cc*s+dd)
        hit=predicted.get(tau)
        tests.append(dict(label=row['label'],original_parameter=str(tau),
            status='CERTIFIED_INHERITED_COLLISION' if hit else 'OUTSIDE_THIS_OBSTRUCTION_NO_CONCLUSION',
            witness_index=hit['index'] if hit else None))
    assert len(tests)==11 and len({r['original_parameter'] for r in tests})==11
    result=dict(status='PASS_POST_SEAL_CONTROL_COMPARISON',cases=tests,
        inputs={str(p.relative_to(ROOT)):sha(p) for p in [OUT/'independent-replay.json',panel,old]},
        boundary='Only addresses read after generic predictions were sealed. No comparison output changes the construction. Nonmatching controls are not certified independent.')
    save('control-comparison.json',result)
    print([(r['label'],r['status']) for r in tests],flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['write','check','compare']);args=parser.parse_args()
    signal.alarm(25)
    if args.action=='compare':compare()
    else:
        result=replay()
        if args.action=='write':save('independent-replay.json',result)
        else:assert read(OUT/'independent-replay.json')==result
        print(result['status'],result['count'],flush=True)
