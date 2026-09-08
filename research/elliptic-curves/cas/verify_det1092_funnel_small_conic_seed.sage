#!/usr/bin/env sage-python
"""Standalone finite-group proof of the smaller conic M18 specialization.

No search constructor, Kummer library or repository rank backend is imported.
The elementary complete-group kernel is retained from the independent uniform
conic checker. All equations, generic sections and conic incidence are checked.
"""
import argparse,hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,gcd,lcm
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def residue(q,p):
    q=QQ(q);assert q.denominator()%p
    return int(q.numerator()*q.denominator().inverse_mod(p)%p)
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
def verify(run,output):
    seed=run/'seeds/conic-small-01/m18.json';packet=read(seed)
    parent_path=ART/'det1092_reduced_parameter_chart_v1/reduced-parent.json'
    original_path=ART/'curve302_recovered_mw17_parent_v1.json'
    chart_path=ART/'det1092_reduced_parameter_chart_v1/generic-proof.json'
    cover_path=ART/'det1092_orbit8044_rank18_base_change_v2.json'
    parent,original,chart,cover=map(read,[parent_path,original_path,chart_path,cover_path])
    assert chart['export_sha256']==sha(parent_path)
    solve=read(run/'conic-solve.json');s=QQ(packet['parameter'])
    R=PolynomialRing(QQ,'t')
    def val(d,t):return R(d['numerator'])(t)/R(d['denominator'])(t)
    ai=[val(d,s)*s.denominator()**k for d,k in zip(parent['a_invariants'],[2,4,6,8,12])]
    assert ai==list(map(QQ,packet['curve'])) and ai[:3]==[0,0,0]
    E=EllipticCurve(QQ,ai);assert E.discriminant()
    points=[tuple(map(QQ,P)) for P in packet['points']]
    assert len(points)==18 and all(E(list(P)) for P in points)
    expected=[(val(x,s)*s.denominator()**4,val(y,s)*s.denominator()**6)
              for x,y in parent['basis_weierstrass_coordinates']]
    assert points[:17]==expected
    a,b,c,d=map(ZZ,chart['parameter_matrix']);tau=(a*s+b)/(c*s+d)
    assert str(tau)==packet['original_parameter']
    q0,q1,q2=map(ZZ,cover['curve_over_Q']['q_coefficients'])
    q=[q0*d*d+q1*b*d+q2*b*b,2*q0*c*d+q1*(a*d+b*c)+2*q2*a*b,
       q0*c*c+q1*a*c+q2*a*a]
    content=gcd(q);assert content.is_square()
    q=[v/content for v in q]
    G=matrix(QQ,[[2*q[2],q[1],0],[q[1],2*q[0],0],[0,0,-2]])
    v=matrix(QQ,3,1,solve['vector'])
    assert G==matrix(QQ,solve['gram']) and (v.transpose()*G*v)[0,0]==0 and v[0,0]/v[1,0]==s
    old=EllipticCurve(QQ,[val(a,tau) for a in original['a_invariants']])
    h=c*s.numerator()+d*s.denominator();w=QQ(chart['weierstrass_u'])
    X,Y=points[17];x=X*w**2/h**4-old.b2()/12
    y=Y*w**3/h**6-(old.a1()*x+old.a3())/2
    assert old([x,y])
    rc,rb,ra=[val(v,tau) for v in cover['lift']['residual_coefficients']]
    f0,f1,f2=[R(v)(tau) for v in cover['lift']['line_coefficients']]
    assert ra*x*x+rb*x+rc==0 and f0+f1*x+f2*y==0
    proof=packet['proof'];primes=sorted(set([int(v['prime']) for v in proof['signatures']]+[proof['no_rational_2_torsion_prime']]))
    rows=[];certificates=[];torsion=False
    for p in primes:
        assert ZZ(p).is_prime(proof=True) and p>2 and E.discriminant()%p
        A,B=residue(ai[3],p),residue(ai[4],p)
        allpoints,doubled,labels,reps,dim=finite_group(p,0,A,B)
        reductions=[]
        for x,y in points:
            z=lcm(x.denominator(),y.denominator());coords=[ZZ(x*z),ZZ(y*z),ZZ(z)]
            common=gcd(coords);coords=[v//common for v in coords];X,Y,Z=[int(v%p) for v in coords]
            assert X or Y or Z
            if not Z:
                assert X==0 and Y;reductions.append(None)
            else:
                reductions.append((X*pow(Z,-1,p)%p,Y*pow(Z,-1,p)%p))
        codes=[labels[P] for P in reductions]
        rows.extend([[(code>>bit)&1 for code in codes] for bit in range(dim)])
        if p==proof['no_rational_2_torsion_prime']:
            assert len(allpoints)%2==1;torsion=True
        certificates.append(dict(prime=p,finite_order=len(allpoints),dimension=dim,
            all_points=allpoints,doubled_points=sorted(doubled,key=lambda P:(-1,-1) if P is None else P),
            representatives=reps,point_reductions=reductions,column_codes=codes))
    M=matrix(GF(2),rows);assert torsion and M.ncols()==18 and M.rank()==18 and M[:,:17].rank()==17
    j=E.j_invariant()
    paths=[seed,parent_path,original_path,chart_path,cover_path,run/'conic-solve.json',Path(__file__)]
    result=dict(status='PASS_STANDALONE_SMALL_CONIC_M18',parameter=str(s),original_parameter=str(tau),
        rank_lower_bound=18,inherited_rank=17,j_numerator_bits=int(abs(j.numerator()).nbits()),
        j_denominator_bits=int(j.denominator().nbits()),proof_primes=primes,finite_group_certificates=certificates,
        matrix_rank=18,inputs={str(p.relative_to(ROOT)):sha(p) for p in paths},
        argument='All18 points lie on the exact specialized parent; first17 equal its generic sections and the extra point lies on the fixed conic multisection. Complete finite group quotients by doubling have18 independent columns; one odd-order good reduction excludes rational2-torsion. Infinite descent proves rational independence. No full rank, conductor, novelty or amplification claim.')
    payload=json.dumps(result,indent=2,sort_keys=True,default=int)+'\n'
    output.parent.mkdir(parents=True,exist_ok=True)
    if output.exists():assert output.read_text()==payload
    else:output.write_text(payload)
    print('PASS_STANDALONE_SMALL_CONIC_M18',str(s),M.rank(),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();verify(a.run.resolve(),a.output.resolve())
