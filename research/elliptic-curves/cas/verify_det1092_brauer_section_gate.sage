#!/usr/bin/env sage-python
"""Independent NS intersections, rational meeting points and local-image replay.

Uses only the generic part of the earlier finite-ring helper; never calls
its historical-point evaluator. No Brauer/Selmer/point search.25-second cap.
"""
import hashlib,json,runpy,signal,itertools
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,vector,identity_matrix,lcm
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
DIR=ART/'det1092_brauer_section_gate_v3';LOCAL=ART/'det1092_seed_local_code_v5'
HELPER=ROOT/'elliptic-curves/cas/verify_det1092_seed_local_code_v2.sage'
helper=runpy.run_path(str(HELPER));Algebra=helper['Algebra'];signs_at_roots=helper['signs_at_roots']
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def provenance(d):
    for p,h in d['inputs'].items():assert sha(ROOT/p)==h,(p,'HASH_MISMATCH')
def retain(p,d):
    s=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==s
    else:p.write_text(s)
def verify():
    parentpath=ART/'curve302_recovered_mw17_parent_v1.json'
    geometrypath=ART/'curve302_parent_geometric_picard19_v1.json'
    paths=[parentpath,geometrypath,DIR/'protocol.json',DIR/'frame.json',
      LOCAL/'protocol.json',HELPER,Path(__file__)]
    parent,geometry,protocol,frame,localprotocol=map(read,paths[:5])
    provenance(protocol);provenance(frame);provenance(localprotocol)
    assert frame['status']=='CONNECTED_DEGREE_ONE_SECTION_GRAPH'
    G=matrix(QQ,parent['generic_height_gram']);NS=matrix(QQ,19)
    NS[0,0]=-2;NS[0,1]=NS[1,0]=1
    for i in range(17):
        for j in range(17):NS[i+2,j+2]=-G[i,j]
    O=vector(QQ,[1,0]+[0]*17)
    sections=[O]+[vector(QQ,[1,G[i,i]/2]+[int(i==j) for j in range(17)]) for i in range(17)]
    assert len(frame['pairs'])==153 and len(frame['tree'])==17
    for row,(i,j) in zip(frame['pairs'],itertools.combinations(range(18),2)):
        assert row['i']==i and row['j']==j
        assert sections[i]*NS*sections[j]==row['intersection_degree']
        assert 2*row['intersection_degree']+4==QQ(row['difference_height'])
    visited={0};W=[]
    for row in frame['tree']:
        i,j=row['i'],row['j'];assert i in visited and j not in visited
        assert sections[i]*NS*sections[j]==1;visited.add(j)
        W.append(list(sections[j][2:]-sections[i][2:]))
    W=matrix(ZZ,W);assert len(visited)==18 and W.det()==QQ(frame['norm6_tree_word_determinant'])==1
    R=PolynomialRing(QQ,'t');t=R.gen();K=R.fraction_field()
    def dec(d):return K(R(d['numerator']))/R(d['denominator'])
    ai=list(map(dec,parent['a_invariants']));basis=[list(map(dec,P)) for P in parent['basis_weierstrass_coordinates']]
    assert ai[:3]==[1,1,1]
    def at_infinity(f,weight):
        degree=f.numerator().degree()-f.denominator().degree()
        assert degree<=weight
        return QQ(0) if degree<weight else f.numerator().leading_coefficient()/f.denominator().leading_coefficient()
    def curve_identity(coeff,x,y):
        a1,a2,a3,a4,a6=coeff
        assert y*y+a1*x*y+a3*y==x**3+a2*x*x+a4*x+a6
        assert (2*y+a1*x+a3)!=0 or (a1*y-3*x*x-2*a2*x-a4)!=0
    meetings=[]
    for row in frame['tree']:
        i,j=row['i'],row['j'];record={'from':row['from'],'to':row['to']}
        if i==0:
            x,y=basis[j-1];dx=x.denominator();dy=y.denominator()
            h=dx.gcd(dx.derivative()).monic()
            assert h.degree()==1 and dx.monic()==h*h and dy.monic()==h**3
            tau=-h[0];assert x.numerator()(tau) and y.numerator()(tau)
            record.update({'parameter':str(tau),'surface_chart':'finite-base zero section',
              'point':'O','pole_orders':[2,3]})
        else:
            xi,yi=basis[i-1];xj,yj=basis[j-1]
            g=(xi-xj).numerator().gcd((yi-yj).numerator())
            if g.degree()==1:
                tau=-g[0]/g[1]
                assert all(f.denominator()(tau) for f in [xi,yi,xj,yj])
                x,y=xi(tau),yi(tau);assert x==xj(tau) and y==yj(tau)
                curve_identity([f(tau) for f in ai],x,y)
                record.update({'parameter':str(tau),'surface_chart':'finite-base affine',
                  'point':[str(x),str(y)]})
            else:
                assert g.degree()==0
                x,y=at_infinity(xi,4),at_infinity(yi,6)
                assert x==at_infinity(xj,4) and y==at_infinity(yj,6)
                infai=[at_infinity(f,w) for f,w in zip(ai,[2,4,6,8,12])]
                curve_identity(infai,x,y)
                record.update({'parameter':'infinity','surface_chart':'s=1/t, X=s^4*x,Y=s^6*y',
                  'point':[str(x),str(y)]})
        meetings.append(record)
    retain(DIR/'rational-intersections.json',{'classification':'explicit rational points joining every old section to O',
      'tree':meetings,'intersection_number_each':1,'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths[:4]}})
    reports=[]
    for index in range(9):
        gp=LOCAL/f'case-{index:02d}-generic.json';ip=LOCAL/f'case-{index:02d}-input.json';rp=LOCAL/f'case-{index:02d}-replay.json'
        generic,inp,oldreplay=map(read,[gp,ip,rp]);paths +=[gp,ip,rp]
        for d in [generic,inp,oldreplay]:provenance(d)
        assert generic['case']==inp['case']==localprotocol['cases'][index]
        tau=QQ(inp['case']['parameter']);a1,a2,a3,a4,a6=[f(tau) for f in ai]
        d=tau.denominator();d*=lcm([((16*a4+8)*d**8).denominator(),((64*a6+16)*d**12).denominator()])
        assert str(d)==inp['scale_denominator']
        P=PolynomialRing(QQ,'x');x=P.gen();f=x**3+5*d**4*x*x+(16*a4+8)*d**8*x+(64*a6+16)*d**12
        assert list(map(str,f.list()))==inp['polynomial']==generic['polynomial']
        assert f.discriminant() and all(a in ZZ for a in f)
        coords=[]
        for xp,yp in basis:
            X,Y=4*d**4*xp(tau),d**6*(8*yp(tau)+4*xp(tau)+4)
            assert f(X)==Y*Y and Y;coords.append([str(X),str(Y)])
        assert coords==inp['generic_cubic_points'];betas=[QQ(X)-x for X,Y in coords]
        algebra=Algebra(f,[P(list(map(QQ,b))) for b in generic['local_order_basis']])
        assert [r['place'] for r in generic['local']]==localprotocol['places']
        local=[]
        for row in generic['local'][:-1]:
            p=row['place'];assert ZZ(p).is_prime(proof=True)
            order=algebra.local_order(p);rank=row['generic_rank'];pivots=row['basis_indices'];cc=row['generic_coordinates']
            C=matrix(GF(2),17,rank,[v for word in cc for v in word]).transpose()
            assert len(cc)==17 and len(pivots)==rank and C.rank()==rank
            assert C.matrix_from_columns(pivots)==identity_matrix(GF(2),rank)
            for j,beta in enumerate(betas):
                denominator=P(1)
                for k,bit in zip(pivots,cc[j]):
                    assert bit in [0,1]
                    if bit:denominator=(denominator*betas[k])%f
                cert=row['relation_square_tests'][j];assert cert['expected_square'] is True
                algebra.square_test((beta*denominator.inverse_mod(f))%f,cert,p)
            assert [r['mask'] for r in row['basis_nonsquare_tests']]==list(range(1,1<<rank))
            for cert in row['basis_nonsquare_tests']:
                v=P(1)
                for j,k in enumerate(pivots):
                    if cert['mask']>>j&1:v=(v*betas[k])%f
                assert cert['expected_square'] is False;algebra.square_test(v,cert,p)
            full=order['prime_count']-1+(p==2)
            assert full==rank==row['full_local_point_dimension']
            local.append({'place':p,'full_point_quotient_dimension':int(full),'generic_rank':int(rank)})
        real=generic['local'][-1];signs=[signs_at_roots(f,QQ(X)) for X,Y in coords]
        assert real['root_signs']==signs
        rank=matrix(GF(2),signs).rank();assert rank==int(len(signs[0])==3)==real['generic_rank']
        local.append({'place':'infinity','full_point_quotient_dimension':int(rank),'generic_rank':int(rank)})
        assert oldreplay['all_generic_local_images_full']
        assert [(r['place'],r['generic_rank']) for r in local]==[(r['place'],r['generic_rank']) for r in oldreplay['local']]
        reports.append({'index':index,'parameter':str(tau),'places':local,
          'normalized_surface_Brauer_2primary_evaluations':'ZERO_ON_THE_WHOLE_LOCAL_ELLIPTIC_FIBRE_AT_ALL_LISTED_PLACES'})
    result={'status':'PASS_RATIONAL_SECTION_CONNECTIVITY_AND_BRAUER_2PRIMARY_LOCAL_BLINDNESS',
      'classification':'new exact application of Brauer evaluation to independently verified geometry and local images',
      'section_pairs_checked':153,'rational_tree_intersections':17,'tree_word_determinant':1,
      'all_generic_sections':'EVERY_NORMALIZED_SURFACE_BRAUER_CLASS_RESTRICTS_TO_ZERO',
      'algebraic_Brauer_quotient':'ZERO_BY_FULL_RATIONAL_TORSION_FREE_PICARD19',
      'cases':reports,'case_place_pairs':189,
      'theorem':'Normalize alpha in Br(X) by its constant restriction to O. The connected rational-intersection tree kills its restriction to every generic basis section. Brauer evaluation on the smooth generic elliptic fibre is a group homomorphism, hence kills every generic section. At the189 certified pairs, generic sections surject onto E(Q_v)/2, and therefore modulo every2^n. Any normalized2-primary Brauer evaluation is consequently zero on all local points there.',
      'boundary':'This neither computes Br(X) nor excludes all transcendental Brauer detection. A nonzero evaluation on a rational302 point would require at least two primes outside the declared footprint, both good for E302, by global reciprocity. The five-place compatibility character is not such an unramified surface Brauer evaluation. No new point or Selmer class is constructed.',
      'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths+[DIR/'rational-intersections.json']}}
    retain(DIR/'replay.json',result)
    print(result['status'],'17 rational intersections;189 complete local pairs',flush=True)
if __name__=='__main__':
    signal.alarm(25);verify()
