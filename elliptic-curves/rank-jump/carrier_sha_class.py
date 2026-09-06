#!/usr/bin/env python3
"""Identify one fixed locally soluble parameter carrier in rational Kummer space."""
import argparse
from pathlib import Path
import subprocess
import retrospective as r

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'CARRIER_SHA_CLASS_PROTOCOL.json'
SOURCE=r.OUT/'rank_jump_disjoint_soluble_carriers_v1.json'
OUTPUT=r.OUT/'rank_jump_carrier_sha_class_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-carrier-sha-class-v1'


def worker():
    from sage.all import QQ,PolynomialRing,FractionField,Curve,Jacobian,EllipticCurve,pari
    from sage.version import version
    policy=r.read(PROTOCOL);old=r.read(SOURCE)['rows'][0]
    R=PolynomialRing(QQ,'v');v=R.gen();f=R(old['geometry']['quartic_coefficients'])
    S=PolynomialRing(QQ,names=('v','w'));sv,sw=S.gens()
    h=Jacobian(Curve(sw**2-sum(f[i]*sv**i for i in range(5))),morphism=True);E=h.codomain()
    K=FractionField(R);W=PolynomialRing(K,'w');L=W.quotient(W.gen()**2-f,'w');w=L.gen()
    coords=[sum(QQ(c)*L(v)**ij[0]*w**ij[1] for ij,c in p.dict().items()) for p in h.defining_polynomials()]
    X,Y=coords[0]/coords[2],coords[1]/coords[2]
    assert X[1]==0 and Y[0]==0
    N=R(X[0]*f);G=R(Y[1]*f*f)
    assert N.degree()==4 and G.degree()==6
    assert G**2==N**3+E.a4()*N*f*f+E.a6()*f**3
    assert f.gcd(N)==1
    Q=PolynomialRing(QQ,'theta');theta=Q.gen();cubic=theta**3+E.a4()*theta+E.a6()
    A=Q.quotient(cubic,'theta');th=A.gen();V=PolynomialRing(A,'v');vv=V.gen()
    numerator=V(N)-th*V(f);beta=numerator[4];normal=numerator/beta
    b=normal[3]/2;c=(normal[2]-b*b)/2
    assert (vv*vv+b*vv+c)**2==normal
    mapping={'quartic':list(map(str,f.list())),'Jacobian_model':list(map(str,E.a_invariants())),
             'N':list(map(str,N.list())),'G':list(map(str,G.list())),
             'beta':list(map(str,beta.lift().list())),
             'square_polynomial_coefficients': [list(map(str,x.lift().list())) for x in (c,b,A(1))],
             'map_degree':4,'syzygy_verified':True,'cubic_square_identity_verified':True}
    r.write_new(WORK/'map.json',mapping);print('map and cubic square identity complete',flush=True)
    Em=EllipticCurve(QQ,list(map(QQ,old['geometry']['minimal_Jacobian_model'])))
    pari.allocatemem(policy['limits']['pari_stack_bytes'],silent=True)
    pari.setrand(policy['limits']['pari_random_seed'])
    PE=pari.ellinit(Em.a_invariants());known=[[QQ(x) for x in P] for P in old['descent']['points']]
    ans=pari.ellrank(PE,policy['limits']['pari_effort'],known)
    lo,hi,ct=map(int,ans[:3]);assert (lo,hi,ct)==(2,2,2)
    returned=ans[3]
    if len(returned)>=2:returned=pari.ellsaturation(PE,returned,2)
    points=[Em([QQ(x) for x in P]) for P in returned]
    torsion=[P for P in Em.torsion_points() if P]
    assert len(torsion)==1 and 2*torsion[0]==Em(0)
    iso=Em.isomorphism_to(E);allpoints=[iso(P) for P in points+torsion]
    model=list(map(str,E.a_invariants()));strings=[list(map(str,P[:2])) for P in allpoints]
    blocks=[];beta_signature=0;column=0
    for p in r.primes(policy['limits']['finite_character_prime_bound']):
        roots=r.roots_at(model[3],model[4],p)
        if not roots or len(roots)!=3:continue
        vals=[(r.mod(str(N[4]),p)-z*r.mod(str(f[4]),p))%p for z in roots]
        if 0 in vals:continue
        blocks.append((p,roots))
        for b in vals:beta_signature|=int(pow(b,(p-1)//2,p)==p-1)<<column;column+=1
    signatures=[r.point_signature(model,P,blocks) for P in strings]
    dim=r.rank(signatures);joint=r.rank(signatures+[beta_signature])
    complete=dim==3;obstructed=complete and joint>dim
    result={'schema':'rank-jump.carrier-sha-class.v1',
            'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (PROTOCOL,SOURCE,Path(__file__),HERE/'retrospective.py')},
            'software':{'sage':version,'pari':str(pari('version()'))},
            'mapping':mapping,'seeded_effort_one_ellrank':str(ans),
            'minimal_model_points_after_2_saturation':[list(map(str,P[:2])) for P in points],
            'minimal_model_torsion_point':list(map(str,torsion[0][:2])),
            'transported_rational_points_including_torsion':strings,
            'proof_primes':[p for p,roots in blocks], 'rational_Kummer_fingerprints':signatures,
            'beta_fingerprint':beta_signature,'rational_Kummer_image_dimension':dim,
            'rational_Kummer_space_complete':complete,'dimension_after_beta':joint,
            'carrier_nonzero_Sha_class_proved':obstructed,
            'carrier_global_solubility':'NO' if obstructed else 'UNKNOWN',
            'boundary':'Exact rank2 and rational2torsion dimension1 bound the full rational Kummer dimension by3. A fourth independent character image proves the particular cover class is not rational. Equal fingerprints would not prove solubility.'}
    r.write_new(WORK/'result.json',result)
    print('generators',len(points),'rational dimension',dim,'with beta',joint,'Sha obstruction',obstructed,flush=True)


def run():
    WORK.mkdir(parents=True,exist_ok=True);log=WORK/'worker.log'
    if not log.exists():
        with log.open('x') as f:
            try:
                p=subprocess.run(['sage','-python',str(Path(__file__).resolve()),'worker'],cwd=r.ROOT,stdout=f,stderr=f,timeout=60)
                status='COMPLETE' if p.returncode==0 else 'FAILED'
            except subprocess.TimeoutExpired:status='TIMEOUT'
        r.write_new(WORK/'execution.json',{'status':status})
    result=r.read(WORK/'result.json') if (WORK/'result.json').exists() else {'status':'UNKNOWN','execution':r.read(WORK/'execution.json'),'log':log.read_text()}
    r.write_new(OUTPUT,result);print(log.read_text());print(r.read(WORK/'execution.json'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['run','worker']);a=p.parse_args()
    worker() if a.mode=='worker' else run()
