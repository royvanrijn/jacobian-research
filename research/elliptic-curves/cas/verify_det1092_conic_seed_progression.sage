#!/usr/bin/env sage-python
"""Independent complete finite groups and uniform congruence specialization.

No constructor import or cubic Kummer backend. Enumerate selected finite
groups and their doubled subgroups using elementary integer formulas.
Verify the global conic map,18 point identities, polynomial congruence
conditions and all nine exact nonsplitting controls.25-second cap.
"""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,gcd,prime_range,prod
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
DIR=ART/'det1092_conic_seed_progression_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def provenance(d):
    for path,digest in d['inputs'].items():assert sha(ROOT/path)==digest,(path,'HASH_MISMATCH')
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True,default=int)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def residue(a,p):
    a=QQ(a);assert a.denominator()%p
    return int(a.numerator()*a.denominator().inverse_mod(p)%p)
def finite_group(p,a,b,c):
    square_roots={}
    for y in range(p):square_roots.setdefault(y*y%p,[]).append(y)
    points=[None]+[(x,y) for x in range(p) for y in square_roots.get((x*x*x+a*x*x+b*x+c)%p,[])]
    point_set=set(points)
    def add(P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P;xx,yy=Q
        if x==xx and (y+yy)%p==0:return None
        if x==xx:
            assert y==yy
            slope=(3*x*x+2*a*x+b)*pow(2*y,-1,p)%p
        else:slope=(yy-y)*pow(xx-x,-1,p)%p
        xxx=(slope*slope-a-x-xx)%p
        yyy=(slope*(x-xxx)-y)%p
        result=(xxx,yyy);assert result in point_set
        return result
    doubled={add(P,P) for P in points}
    labels={P:0 for P in doubled};representatives=[None];dimension=0
    for P in points:
        if P in labels:continue
        new=[]
        for j,Q in enumerate(representatives):
            R=add(P,Q);new.append(R)
            for D in doubled:
                V=add(R,D);value=j|(1<<dimension)
                assert V not in labels
                labels[V]=value
        representatives+=new;dimension+=1
    assert set(labels)==point_set and len(representatives)==2**dimension
    assert len(points)==len(doubled)*2**dimension
    return points,doubled,labels,representatives,dimension
def verify():
    paths=[DIR/name for name in ['protocol.json','input.json','finite-code.json','progression.json']]
    protocol,inp,finite,progression=map(read,paths)
    for d in [protocol,inp,finite,progression]:provenance(d)
    assert protocol['limits']['conic_parameter']==inp['conic_parameter']==0
    assert protocol['limits']['new_specialization_addresses']==1
    assert protocol['limits']['finite_prime_bound']==1009
    parent=read(ART/'curve302_recovered_mw17_parent_v1.json')
    cover=read(ART/'det1092_orbit8044_rank18_base_change_v2.json');provenance(cover)
    roster=read(ART/'det1092_rr_generic_point_controls_v2/protocol.json')
    R=PolynomialRing(QQ,'t');K=R.fraction_field();U=PolynomialRing(QQ,'u');F=U.fraction_field();u=U.gen()
    def dec(d,R,K):return K(R(d['numerator']))/R(d['denominator'])
    q=R(cover['curve_over_Q']['q_coefficients'])
    t0,s0=map(QQ,cover['curve_over_Q']['rational_point'])
    assert q.degree()==2 and q.discriminant() and s0*s0==q(t0)
    T=F(t0)+(q.derivative()(t0)-2*s0*u)/(u*u-q[2])
    S=F(s0)+u*(T-t0)
    assert max(T.numerator().degree(),T.denominator().degree())==2
    assert S*S==q(T) and (S-s0)/(T-t0)==u
    assert T==dec(inp['base_map'],U,F) and S==dec(inp['conic_ordinate'],U,F)
    assert str(T(0))==inp['parameter']==progression['base_parameter_at_u_zero']
    def pull(d):
        v=dec(d,R,K)
        return F(v.numerator()(T))/v.denominator()(T)
    ai=[pull(d) for d in parent['a_invariants']]
    functions=[[pull(d) for d in P] for P in parent['basis_weierstrass_coordinates']]
    c,b,a=[pull(d) for d in cover['lift']['residual_coefficients']]
    f0,f1,f2=[R(d)(T) for d in cover['lift']['line_coefficients']]
    h=R(cover['splitting']['discriminant_square_factor'])(T)
    assert b*b-4*a*c==h*h*q(T)
    x=(-b+h*S)/(2*a);y=-(f0+f1*x)/f2
    assert a*x*x+b*x+c==0 and f0+f1*x+f2*y==0
    functions.append([x,y])
    assert ai==[dec(d,U,F) for d in inp['weierstrass_functions']]
    assert functions==[[dec(d,U,F) for d in P] for P in inp['point_functions']]
    a1,a2,a3,a4,a6=ai
    for xx,yy in functions:
        assert yy*yy+a1*xx*yy+a3*yy==xx**3+a2*xx*xx+a4*xx+a6
    point_polys=[[U(list(map(QQ,d))) for d in P] for P in inp['projective_point_polynomials']]
    pairs=[[U(list(map(QQ,d))) for d in pair] for pair in inp['coefficient_polynomial_pairs']]
    Tpair=[U(list(map(QQ,d))) for d in progression['base_map_polynomial_pair']]
    for polys in point_polys+pairs+[Tpair]:
        assert all(c in ZZ for poly in polys for c in poly)
        assert gcd(polys).degree()==0
    assert len(point_polys)==18 and len(pairs)==5
    for P,(xx,yy) in zip(point_polys,functions):
        assert P[0]==xx*P[2] and P[1]==yy*P[2]
    for pair,v in zip(pairs,ai):assert pair[0]==v*pair[1]
    assert Tpair[0]==T*Tpair[1]
    E=EllipticCurve(QQ,[a(0) for a in ai]);assert E.discriminant()
    assert list(map(str,E.a_invariants()))==inp['equation_at_zero']
    points=[E([f(0) for f in P]) for P in point_polys]
    assert [list(map(str,P)) for P in points]==inp['points_at_zero']
    assert [r['p'] for r in finite['exposures']]==list(map(int,prime_range(3,1010)))
    good={};exposure_counts={}
    for r in finite['exposures']:
        p=r['p'];status=r['status']
        exposure_counts[status]=exposure_counts.get(status,0)+1
        if any(a.denominator()%p==0 for a in E.a_invariants()):
            assert status=='SKIP_COEFFICIENT_DENOMINATOR';continue
        if E.discriminant().valuation(p)!=0:
            assert status=='SKIP_BAD_REDUCTION';continue
        assert status=='PASS_CUBIC_KUMMER_REDUCTION' and r==read(DIR/f'prime-{p}.json')
        aa,bb,cc=[residue(v,p) for v in [E.b2(),8*E.b4(),16*E.b6()]]
        roots=[v for v in range(p) if (v**3+aa*v*v+bb*v+cc)%p==0]
        assert roots==r['rational_roots'] and [cc,bb,aa,1]==r['cubic']
        good[p]=r
    selected=progression['selected_proof_primes'];assert sorted(set(selected))==selected
    assert all(p in good for p in selected)
    assert finite['no_rational_2_torsion_prime'] in selected
    rows=[];group_certificates=[]
    for p,stability in zip(selected,progression['stability']):
        assert stability['p']==p and stability['exponent']==1
        # The actual certificate needs only squarefree congruences. All
        # displayed maps have integral polynomial coefficients, denominators
        # that are units at0, and nonzero projective reductions there.
        for pair in pairs:
            assert residue(pair[1](0),p)!=0
        assert residue(Tpair[1](0),p)!=0
        reductions=[]
        for P in point_polys:
            coords=[residue(f(0),p) for f in P]
            assert any(coords)
            if coords[2]==0:
                assert coords[0]==0 and coords[1]
                reductions.append(None)
            else:
                xx=coords[0]*pow(coords[2],-1,p)%p
                yy=coords[1]*pow(coords[2],-1,p)%p
                xc=4*xx%p
                yc=(8*yy+4*residue(ai[0](0),p)*xx+4*residue(ai[2](0),p))%p
                reductions.append((xc,yc))
        cc,bb,aa,one=good[p]['cubic'];assert one==1
        allpoints,doubled,labels,reps,dim=finite_group(p,aa,bb,cc)
        assert all(P in labels for P in reductions)
        codes=[labels[P] for P in reductions]
        block=[[(c>>j)&1 for c in codes] for j in range(dim)]
        rows.extend(block)
        # Independent group cosets and recorded cubic characters must span
        # the same subspace on the18 ordered points.
        A=matrix(GF(2),len(block),18,[x for r in block for x in r])
        B=matrix(GF(2),len(good[p]['rows']),18,[x for r in good[p]['rows'] for x in r])
        assert A.row_space()==B.row_space()
        if p==finite['no_rational_2_torsion_prime']:
            assert len(allpoints)%2==1 and not good[p]['rational_roots']
        group_certificates.append({'p':p,'group_order':len(allpoints),'double_order':len(doubled),
            'all_points':allpoints,'doubled_points':sorted(doubled,key=lambda P:(P is not None,P)),
            'coset_representatives':reps,'point_reductions':reductions,'column_codes':codes,
            'independent_matrix_rows':block})
    assert len(progression['stability'])==len(selected)
    M=matrix(GF(2),len(rows),18,[x for r in rows for x in r])
    assert M.rank()==18 and M[:,:17].rank()==17
    B=prod(ZZ(p) for p in selected)
    assert str(B)==progression['B'] and B.nbits()==progression['B_bits']==113
    controls=[]
    assert len(progression['controls'])==len(roster['cases'])==9
    for row,case in zip(progression['controls'],roster['cases']):
        assert row['label']==case['label'] and row['parameter']==case['parameter']
        value=q(QQ(case['parameter']));assert str(value)==row['q_value']
        assert value>=0
        n,d=value.numerator(),value.denominator()
        a,b=ZZ(row['numerator_floor_sqrt']),ZZ(row['denominator_floor_sqrt'])
        assert a*a<=n<(a+1)**2 and b*b<=d<(b+1)**2
        split=a*a==n and b*b==d
        assert split==row['split'] and not split
        controls.append({'label':row['label'],'parameter':row['parameter'],'cover_split':False,
                         'fibre_rank_or_other_seed_existence':'NOT_EXCLUDED'})
    result={'status':'PASS_INDEPENDENT_ORACLE_FREE_RANK18_PROGRESSION',
        'classification':'new constructive deduction and verified application',
        'B':str(B),'base_parameter_at_u_zero':str(T(0)),
        'statement':'For every integer n, u=B*n gives a smooth rational fibre at t(u), with the17 inherited sections and the explicit conic lift linearly independent. Thus each fibre has rank>=18 and the conic point lies outside the inherited rational span.',
        'uniformity_proof':'Every coefficient denominator and the base-map denominator is a unit modulo each proof prime at0; every projective point triple is nonzero. Integral polynomial maps at u=B*n have exactly the same reductions. Complete finite-group quotient rank18 and an odd-order reduction remain unchanged for every integer n.',
        'infinitely_many_distinct_parameters':'The rational map t(u) has degree2 and no pole on this progression; each value has at most two preimages.',
        'rank_at_zero':18,'generic_subgroup_rank_at_zero':17,
        'no_rational_2_torsion_prime':finite['no_rational_2_torsion_prime'],
        'finite_group_certificates':group_certificates,'exposure_counts':exposure_counts,
        'controls':controls,'limits':protocol['limits'],
        'boundary':'One extra direction from the existing conic. The cover is nonsplit at302 and all eight controls, so this does not explain the historical302 seed or predict their ranks. No amplification, catalogue comparison or exact-rank upper bound.',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths+[Path(__file__)]}}
    retain(DIR/'replay.json',result)
    print(result['status'],'rank18 for every integer n','B='+str(B),flush=True)
if __name__=='__main__':
    signal.alarm(25);verify()
