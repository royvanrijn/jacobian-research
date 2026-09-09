#!/usr/bin/env sage-python
"""Independent Gram, RR, odd-divisor descent and conic-map replay.

Does not import the producers or repeat rational kernel/factor discovery.
"""
import csv,hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,vector,gcd,prime_range
signal.alarm(25)
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_rational_bisection_index_v1'
SELECT=OUT/'index-and-selection.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
selection=read(OUT/'index-and-selection.json')
for name,digest in selection['inputs'].items():assert sha(ROOT/name)==digest
parent=read(ART/'curve302_recovered_mw17_parent_v1.json')
G=matrix(ZZ,parent['generic_height_gram']);assert G.det()==1092
assert all(v%2==0 for v in G.diagonal())
rad=vector(ZZ,selection['radical_word'])
assert matrix(GF(2),G).rank()==16 and any(rad)
assert all(v%2==0 for v in G*rad) and rad*G*rad==180
assert (rad*G*rad)%4==0
stored=[]
with (ART/'curve302_parent_degree2_multisection_orbits_v1.tsv').open() as stream:
    for row in csv.DictReader(stream,delimiter='\t'):
        if row['category']=='rational':
            w=list(map(int,row['parent_MW17_w'].split()))
            stored.append((sum(map(abs,w)),int(row['orbit_mask']),w))
assert len(stored)==40917 and min(v[0] for v in stored)==3
chosen=sorted(v for v in stored if v[0]==3)
assert [v[1] for v in chosen]==[8044,47755,103186]
R=PolynomialRing(QQ,'t');F=R.fraction_field();t=R.gen()
def dec(v):return F(R(v['numerator']))/R(v['denominator'])
E=EllipticCurve(F,[dec(v) for v in parent['a_invariants']])
basis=[E([dec(v) for v in P]) for P in parent['basis_weierstrass_coordinates']]
U=PolynomialRing(QQ,'u');K=U.fraction_field()
def decu(v):return K(U(v['numerator']))/U(v['denominator'])
Xring=PolynomialRing(F,'x');x=Xring.gen()
Sring=PolynomialRing(F,'s');s=Sring.gen()
prime_barrier=read(ART/'det1092_prime_division_barrier_v1/replay.json')
assert prime_barrier['status']=='PASS_PRIME_DIVISION_BARRIER_INPUTS_AND_IDENTITIES'
controls=read(OUT/'controls.json')
for name,digest in controls['inputs'].items():assert sha(ROOT/name)==digest
results=[]
for _,orbit,ww in chosen:
    data=read(OUT/('orbit-%d.json'%orbit));w=vector(ZZ,ww)
    assert data['word']==ww and w*G*w==10
    for name,digest in data['inputs'].items():assert sha(ROOT/name)==digest
    trace=-sum((v*P for v,P in zip(w,basis)),E(0))
    assert trace==E([dec(v) for v in data['negative_trace']])
    f0,f1,f2=[R(v) for v in data['line_coefficients']]
    assert f0.degree()<=9 and f1.degree()<=5 and f2.degree()<=3 and f2
    assert gcd([f0,f1,f2]).degree()==0
    assert f0+f1*trace[0]+f2*trace[1]==0
    pole=R(trace[0].denominator()).sqrt()
    nx=R(trace[0]*pole**2);ny=R(trace[1]*pole**3)
    assert [pole.degree(),nx.degree(),ny.degree()]==[3,10,15]
    columns=[term*t**j for bound,term in [(9,pole**3),(5,nx*pole),(3,ny)] for j in range(bound+1)]
    interpolation=matrix(QQ,19,20,lambda i,j:columns[j][i])
    proof_prime=None
    for p in prime_range(3,212):
        if any(v.denominator()%p==0 for v in interpolation.list()):continue
        if matrix(GF(p),interpolation).rank()==19:proof_prime=int(p);break
    assert proof_prime is not None
    a1,a2,a3,a4,a6=E.a_invariants();line=f0+f1*x
    eliminated=line**2-a1*x*line*f2-a3*line*f2-f2*f2*(x**3+a2*x*x+a4*x+a6)
    c,b,a=[dec(v) for v in data['residual_coefficients']]
    assert eliminated==(x-trace[0])*(a*x*x+b*x+c)
    q=R(data['q']);h=R(data['h'])
    assert q.degree()==2 and q.discriminant()!=0 and b*b-4*a*c==h*h*q
    x0,x1,y0,y1=[R(v) for v in data['elliptic_quadratic_maps']]
    assert x1 and x0==-b/(2*a) and x1==h/(2*a)
    assert f0+f1*x0+f2*y0==0 and f1*x1+f2*y1==0
    xx=x0+x1*s;yy=y0+y1*s
    assert (yy*yy+a1*xx*yy+a3*yy-xx**3-a2*xx*xx-a4*xx-a6)%(s*s-q)==0
    index=data['odd_section_index'];degrees=vector(ZZ,G.diagonal())-G*w
    assert degrees[index]==data['intersection_degree'] and degrees[index]%2
    assert index==next(i for i,v in enumerate(degrees) if v%2)
    point=data['rational_point']
    if point['type']=='finite':
        divisor=R(point['starting_divisor']);ordinate=R(point['starting_ordinate'])
        assert divisor.degree()==degrees[index] and divisor.is_monic()
        P=basis[index]
        for v,expected in [(P[0],x0+x1*ordinate),(P[1],y0+y1*ordinate)]:
            assert gcd(R(v.denominator()),divisor)==1
            assert (R(v.numerator())-expected*R(v.denominator()))%divisor==0
        assert (ordinate*ordinate-q)%divisor==0
        for step in data['norm_descent_steps']:
            assert divisor==R(step['divisor']) and ordinate==R(step['ordinate'])
            quotient=R(step['quotient'])
            assert ordinate*ordinate-q==divisor*quotient
            assert 0<quotient.degree()<divisor.degree() and quotient.degree()%2
            divisor=quotient.monic();ordinate%=divisor
        assert divisor.degree()==1
        t0=QQ(point['t']);s0=QQ(point['s'])
        assert divisor(t0)==0 and ordinate(t0)==s0 and s0*s0==q(t0)
    else:
        assert QQ(point['leading_sqrt'])**2==q[2]
    T=decu(data['base_map']);W=decu(data['conic_ordinate'])
    assert W*W==q(T) and max(T.numerator().degree(),T.denominator().degree())==2
    panel=next(r for r in controls['rows'] if r['orbit']==orbit)
    assert len(panel['cases'])==9
    for case in panel['cases']:
        v=q(QQ(case['parameter']));assert str(v)==case['q_value']
        assert not case['split'] and not v.is_square()
        witness=case['local_nonsquare_witness'];assert witness
        p=ZZ(witness['p']);valuation=ZZ(v.valuation(p));unit=v/p**valuation
        residue=ZZ(unit.numerator()*unit.denominator().inverse_mod(p)%p)
        assert int(valuation)==witness['valuation'] and int(residue)==witness['unit_residue']
        assert valuation%2 or residue not in {ZZ(j*j%p) for j in range(int(p))}
    results.append({'orbit':orbit,'RR_rank_lower_bound_prime':proof_prime,
                    'odd_divisor_degree':int(degrees[index]),'norm_descent_steps':len(data['norm_descent_steps']),
                    'rational_base_change_degree':2,'function_field_rank_lower_bound':18,
                    'old_control_splits':0})
report={'status':'PASS_INDEPENDENT_RATIONAL_BISECTION_INDEX_AND_THREE_MAPS',
    'classification':'new conic-index theorem and independently replayed generic constructions',
    'all_40917_bisection_orbits_are_Q_rational':True,
    'uniform_proof':'The mod2 radical has dimension1 and norm0 mod4; norm10 vectors cannot have even Gram column. An odd-degree rational section intersection and a degree2 fibre intersection give a degree1 rational divisor on each smooth genus0 bisection. Genus0 Riemann--Roch gives a rational point.',
    'constructions':results,'boundary':'Only three fixed stored-coordinate representatives were constructed; all27 old-control fibres are nonsplit. This is not an exclusion of the remaining40914 translation orbits, nor a302 seed selector.',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [SELECT,OUT/'controls.json',Path(__file__),
             ART/'det1092_prime_division_barrier_v1/replay.json',*[OUT/('orbit-%d.json'%i) for i in [8044,47755,103186]]]}}
dest=OUT/'independent-replay.json'
if dest.exists():assert read(dest)==report
else:
    with dest.open('x') as stream:json.dump(report,stream,indent=2,sort_keys=True);stream.write('\n')
print(report['status'],[(r['orbit'],r['odd_divisor_degree']) for r in results],flush=True)
