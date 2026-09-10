#!/usr/bin/env sage-python
"""Independent exact replay of class1 point transport and fourteen carriers.

Uses Python rational group arithmetic for public words, inverse maps for
transport, exhaustive finite elliptic groups for rank18, and Sylvester
determinants for the polynomial-resultant certificate.
"""
import argparse,hashlib,json,sys
from fractions import Fraction as F
from pathlib import Path
from sage.all import QQ,ZZ,GF,EllipticCurve,PolynomialRing,matrix,pari
from dataclasses import asdict
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
PROTOCOL=ART/'curve302_class1_point_transport_protocol_v1.json'
INPUT=ART/'curve302_class1_point_transport_v1.json'
CARRIERS=ART/'curve302_class1_fibre_carriers_v1.json'
OUTPUT=ART/'curve302_class1_bridge_verification_v1.json'
WORK=ROOT/'artifacts/local/elliptic-curves/curve302-class1-bridge-replay-v1'
sys.path.insert(0,str(Path(__file__).parent))
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
import icarm_curve302 as curve
from half_lattice_pointed_sieve import linear_combination_python
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from mod2_reduction_independence import Mod2ReductionSignature,combined_mod2_rank
from certify_compact_r17_candidates import short_curve_has_no_rational_2_torsion_modular_certificate
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def evaluate(coeffs,t):
    ans=F(0)
    for c in reversed(coeffs):ans=ans*t+F(c)
    return ans
def dec(c,t):return evaluate(c['numerator'],t)/evaluate(c['denominator'],t)
def sylvester_resultant(f,g):
    m,n=f.degree(),g.degree();K=f.base_ring()
    a=list(reversed(f.list()));b=list(reversed(g.list()))
    rows=[[K(0)]*i+a+[K(0)]*(n-1-i)for i in range(n)]
    rows += [[K(0)]*i+b+[K(0)]*(m-1-i)for i in range(m)]
    return matrix(K,rows).det()
def finite_rank(model,points,primes,torsion):
    # Point identities are checked once with PARI; exhaustive finite groups
    # then evaluate reductions. No repeated large rational-identity checks.
    E=pari.ellinit([pari(QQ(c))for c in model])
    for point in points:assert pari.ellisoncurve(E,pari([QQ(c)for c in point]))
    cache=ReductionCache(MemoryFactStore());signatures=[]
    for p in primes:
        row,table,_=cache.quotient(model,p)
        values=[]
        for point in points:
            reduced=tuple((c.numerator%p)*pow(c.denominator,-1,p)%p for c in point)
            assert reduced in table;values.append(table[reduced])
        signatures.append(Mod2ReductionSignature(p,row['group_order'],row['doubled_order'],row['dimension'],
            tuple(tuple((mask>>i)&1 for mask in values)for i in range(row['dimension']))))
    assert combined_mod2_rank(signatures,len(points))==len(points)
    assert short_curve_has_no_rational_2_torsion_modular_certificate(model,torsion)
    return {'rank_lower_bound':len(points),'no_rational_2_torsion_prime':torsion,'signatures':[asdict(s)for s in signatures]}
def verify(case_index=None,collect=False):
    protocol,data,carriers=map(read,[PROTOCOL,INPUT,CARRIERS])
    for record in [protocol,carriers]:
        for name,digest in record['bindings'].items():assert sha(ROOT/name)==digest,name
    assert data['protocol_sha256']==sha(PROTOCOL) and data['complete']
    source=read(ART/'curve302_recovered_mw17_parent_v1.json')
    prefix='x1092_class1_realization_'
    tr,rr,eq,normalized,compact=[read(ART/(prefix+s+'_v1.json'))for s in ['trace','rr','equation','parent','compact_parent']]
    M=matrix(ZZ,source['basis_embedding_in_public_D'])
    C=matrix(ZZ,[r['public_word']for r in protocol['rows']]).transpose()
    assert M.augment(C).det()==protocol['integral_basis_determinant']==-1
    assert len(data['cases'])==14 and [r['id']for r in data['cases']]==[r['id']for r in protocol['rows']]
    reports=[]
    output_bindings={str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),PROTOCOL,INPUT,CARRIERS]}
    for index,(r,p) in enumerate(zip(data['cases'],protocol['rows'])):
        if case_index is not None and index!=case_index:continue
        if collect:
            previous=read(WORK/('case-%02d.json'%index));assert previous['bindings']==output_bindings
            assert previous['status']=='PASS' and len(previous['point_transport_and_rank_cases'])==1
            assert previous['point_transport_and_rank_cases'][0]['id']==r['id']
            reports.extend(previous['point_transport_and_rank_cases']);continue
        assert r['protocol_sha256']==sha(PROTOCOL) and r['public_word']==p['public_word']
        target=linear_combination_python(curve.short_coefficients(),curve.SHORT_POINTS,p['public_word'])
        sx,sy=target[0]/36,target[1]/216
        assert list(map(F,r['source_short_point']))==[sx,sy]
        x,y=map(F,r['source_point']);assert sx==x+F(5,12) and sy==y+(x+1)/2
        u,s=F(r['u']),F(r['s']);assert dec(compact['u_of_s'],s)==u
        model=tuple(dec(c,s)for c in compact['a_invariants'])
        assert tuple(map(F,r['compact_model']))==model
        generic=tuple(tuple(dec(c,s)for c in pt)for pt in compact['basis_weierstrass_coordinates'])
        assert generic==tuple(tuple(map(F,pt))for pt in r['generic_points'])
        point=tuple(map(F,r['compact_point']))
        # Recover pointed cubic coordinates from the compact chart.
        eb=tuple(dec(c,u)for c in eq['a_invariants'])
        assert eb[0]==eb[2]==0
        gauge=dec(normalized['gauge_from_pointed_cubic'],u)
        den=evaluate(compact['u_of_s']['denominator'],s)
        scale=gauge*den**2/F(compact['scale']);assert scale==F(r['compact_scale'])
        X=point[0]/scale**2-eb[1]/3;Y=point[1]/scale**3
        assert [X,Y]==list(map(F,r['pointed_cubic_point']))
        assert Y*Y==X**3+eb[1]*X*X+eb[3]*X+eb[4]
        R=PolynomialRing(QQ,'t');qc=R([QQ(dec(c,u))for c in eq['quartic_coefficients']])
        t0,q0=[dec(c,u)for c in eq['quartic_zero']]
        shifted=qc(R.gen()+QQ(t0));e,d,c,b,a=[F(str(shifted[i]))for i in range(5)]
        z=(2*q0*Y+d*X+2*e*b)/(X*X-4*e*a);t=z+t0
        v=(X*z*z-d*z-2*e)/(2*q0)
        assert t==0 and [t,v]==list(map(F,r['quartic_point']))
        sf=evaluate([dec(c,u)for c in eq['radical_square_factor']],t)
        ds=evaluate([dec(c,u)for c in eq['radical_denominator_sqrt']],t)
        h,nx,ny=[evaluate(tr[k],t)for k in ['h','Nx','Ny']]
        a0,b0,a1,b1=[evaluate(rr[k],t)for k in ['a0','b0','a1','b1']]
        m=(a1-u*a0)/((u*b0-b1)*h);xp,yp=nx/h**2,ny/h**3
        rx=(m*m-xp+v*sf/ds)/2;ry=m*(rx-xp)-yp
        assert (rx,ry)==(sx,sy)
        assert u==(a1*(rx*h*h-nx)+b1*(ry*h**3+ny))/(a0*(rx*h*h-nx)+b0*(ry*h**3+ny))
        # Full finite groups, independent of the producer's cubic Legendre map.
        primes=[b['prime']for b in r['finite_blocks']]
        assert len(set(primes))==len(primes)and all(3<=p<=2000 for p in primes)
        proof=finite_rank(model,(*generic,point),primes,r['no_rational_2_torsion_prime'])
        assert proof['rank_lower_bound']==18
        assert r['generic_mod2_rank']==17 and r['augmented_mod2_rank']==18
        assert r['conclusion']=='OUTSIDE_GENERIC_RATIONAL_SPAN'
        reports.append({'id':r['id'],'rank_lower_bound':18,'finite_group_proof':proof})
        print('PASS inverse transport and independent finite rank18',r['id'],flush=True)
    if case_index is not None:
        assert len(reports)==1
        return {'schema':'curve302.class1-bridge-case-verification.v1','status':'PASS',
                'bindings':output_bindings,'point_transport_and_rank_cases':reports}
    assert len({r['s']for r in data['cases']})==data['distinct_B_parameters']==14
    assert data['outside_generic_rational_span_count']==14
    # Reconstruct whole A-carriers, including their rational point, trace and
    # degree2 map. The frozen modular certificate proves every genus assertion.
    R=PolynomialRing(QQ,'t');Fq=R.fraction_field();t=R.gen()
    A,B,h,nx,ny=[R(tr[k])for k in ['A','B','h','Nx','Ny']]
    xp,yp=Fq(nx/h**2),Fq(ny/h**3)
    assert yp*yp==xp**3+A*xp+B
    polys=[]
    for record,r in zip(carriers['carriers'],data['cases']):
        assert record['id']==r['id']and record['u']==r['u']
        q=R(record['branch_quartic']);assert q.degree()==4
        assert q==R([QQ(dec(c,F(r['u'])))for c in eq['quartic_coefficients']])
        parse=lambda c:Fq(R(c['numerator']))/R(c['denominator'])
        x0,x1,y0,y1=[parse(record['source_maps'][k])for k in ['x0','x1','y0','y1']]
        assert x1!=0
        m=y1/x1
        assert m*m==2*x0+xp and y0==m*(x0-xp)-yp
        assert y0*y0+y1*y1*q==x0**3+3*x0*x1*x1*q+A*x0+B
        assert 2*y0*y1==3*x0*x0*x1+x1**3*q+A*x1
        v=QQ(record['rational_point_at_A_zero'][1]);assert v and q(0)==v*v
        at0=lambda f:f.numerator()(0)/f.denominator()(0)
        assert [at0(x0)+at0(x1)*v,at0(y0)+at0(y1)*v]==list(map(QQ,r['source_short_point']))
        polys.append(q)
    cert=carriers['finite_certificate'];p=cert['prime'];assert ZZ(p).is_prime(proof=True)and 3<=p<=2000
    for q in polys:assert all(c.denominator()%p for c in q.list())
    qs=[q.change_ring(GF(p))for q in polys];assert all(q.degree()==4 for q in qs)
    for q,value in zip(qs,cert['squarefree_resultants']):
        assert int(sylvester_resultant(q,q.derivative()))==value!=0
    assert len(cert['pairs'])==91
    assert [(r['i'],r['j'])for r in cert['pairs']]==[(i,j)for i in range(14)for j in range(i+1,14)]
    for r in cert['pairs']:
        assert int(sylvester_resultant(qs[r['i']],qs[r['j']]))==r['resultant_mod_p']!=0
        assert r['joint_normalization_genus']==5
    assert carriers['shared_quadratic_squareclasses']==carriers['genus_zero_or_one_joint_pairs']==0
    # Generic automorphism-return theorem applies to this smooth non-special-j
    # fibre: Aut(E302,O)={+1,-1}, and the generic17 basis is the full MW lattice.
    assert curve.weierstrass_invariants()[4]and curve.weierstrass_invariants()[5]
    print('PASS fourteen genus1 carriers;91 genus5 joint normalizations; ordinary fibre automorphism signs only',flush=True)
    return {'schema':'curve302.class1-bridge-verification.v1','status':'PASS',
            'bindings':output_bindings,
            'point_transport_and_rank_cases':reports,'distinct_B_parameters':14,
            'outside_generic_rational_span_count':14,'carrier_genus':1,'joint_genus':5,'pair_count':91,
            'boundary':'Native B-fibres only. Full twist ranks, other shared carriers and seed propagation remain UNKNOWN.'}
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true')
    parser.add_argument('--case',type=int);parser.add_argument('--collect',action='store_true');args=parser.parse_args()
    assert not(args.case is not None and args.collect)
    result=verify(args.case,args.collect)
    target=OUTPUT if args.case is None else WORK/('case-%02d.json'%args.case)
    target.parent.mkdir(parents=True,exist_ok=True)
    if args.check:assert read(target)==json.loads(json.dumps(result))
    else:
        with target.open('x')as stream:json.dump(result,stream,indent=2,sort_keys=True);stream.write('\n')
