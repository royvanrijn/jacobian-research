#!/usr/bin/env sage -python
"""Bounded, deterministic global2-descent on lexicographic paired base[0,1].

Trust boundary: PARI's complete ellrank/ell2cover implementations. Exact
normalization, local witnesses, Lucas primality certificates, points and cover
maps are also exported for an independent standard-library arithmetic replay.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import resource
import time
from sage.all import (EllipticCurve, PolynomialRing, QQ, ZZ, gcd, lcm,
                       pari, sage_eval, set_random_seed, version)

ROOT=Path(__file__).resolve().parents[2]
PILOT=ROOT/"artifacts/generated-results/hill1152-transfer-20261113-v2/result.json"
JACS=ROOT/"artifacts/generated-results/x1152-pair-jacobians-20261114-v1.json"


def rows(f):return [str(f[i]) for i in range(max(0,int(f.degree()))+1)]


def run():
    resource.setrlimit(resource.RLIMIT_CPU,(30,35))
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    set_random_seed(11520001);pari.setrand(11520001)
    started=time.process_time()
    pilot=json.loads(PILOT.read_text());jacs=json.loads(JACS.read_text())
    pair=list(min(map(tuple,pilot["locally_surviving_disjoint_pairs"])))
    assert pair==[0,1]
    raw=[list(map(QQ,pilot["quadratic_covers"][i]["q"])) for i in pair]
    primitive=[];contents=[]
    for q in raw:
        den=lcm([z.denominator() for z in q]);ints=[ZZ(z*den) for z in q]
        g=gcd(ints);content=QQ(g)/den
        assert content>0 and content.is_square()
        primitive.append([ZZ(z/g) for z in ints]);contents.append(content)
    q,r=primitive
    m=q[2]*r[0]+q[0]*r[2]-QQ(q[1]*r[1])/2
    dq=q[1]**2-4*q[2]*q[0];dr=r[1]**2-4*r[2]*r[0]
    resultant=m*m-QQ(dq*dr)/4
    assert resultant.denominator()==1
    E=EllipticCurve(QQ,[0,4*m,0,dq*dr,0])
    k=2/(contents[0]*contents[1]).sqrt()
    source=next(x for x in jacs["pairs"] if x["pair"]==pair)
    a,b=QQ(source["jacobian_a_invariants"][1]),QQ(source["jacobian_a_invariants"][3])
    assert k*k*a==E.a2() and k**4*b==E.a4()
    assert not QQ(E.a2()**2-4*E.a4()).is_square()
    # Prime certificates use the elementary Lucas full-order criterion.
    primes={}
    def prime_certificate(p):
        p=ZZ(p);key=str(p)
        if key in primes:return
        if p==2:primes[key]={"prime":"2"};return
        factors=list((p-1).factor(proof=True))
        for z,e in factors:prime_certificate(z)
        g=next(a for a in range(2,10000)
               if pow(a,int(p-1),int(p))==1 and
               all(gcd(pow(a,int((p-1)//z),int(p))-1,p)==1 for z,e in factors))
        primes[key]={"prime":key,"p_minus_one":[[str(z),int(e)] for z,e in factors],
                     "full_order_witness":g}
    bad=ZZ(2*dq*dr*resultant);factorization=list(abs(bad).factor(proof=True))
    local=[]
    def ev(f,t):return sum(a*t**i for i,a in enumerate(f))
    def square_at(v,p):
        if not v:return True
        valuation=v.valuation(p);unit=v//p**valuation
        return valuation%2==0 and (unit%8==1 if p==2 else pow(int(unit%p),int((p-1)//2),int(p))==1)
    for p,e in factorization:
        prime_certificate(p)
        t=next(t for t in range(1000) if all(square_at(ev(f,t),p) for f in primitive))
        local.append({"prime":str(p),"exponent_in_bad_product":int(e),"integer_parameter":t})
    assert all(ev(f,0)>0 for f in primitive)
    rank=E.pari_curve().ellrank(0)
    rank_header=[int(rank[i]) for i in range(3)]
    assert rank_header==[2,2,0]
    points=[[str(z) for z in P] for P in rank[3]]
    assert points and all(E(list(map(QQ,P))) for P in points)
    covers=E.pari_curve().ell2cover()
    assert len(covers)==3
    R=PolynomialRing(QQ,"x");x=R.gen();K=R.fraction_field()
    S=PolynomialRing(K,"Y");Y=S.gen();coverrows=[]
    for f,P in covers:
        f=R(sage_eval(str(f),locals={"x":x}))
        L=S.quotient(Y**2-K(f),"z");z=L.gen()
        images=[sage_eval(str(v),locals={"x":L(x),"y":z}) for v in P]
        assert len(images)==2
        xx,yy=images
        assert yy**2==xx**3+E.a2()*xx**2+E.a4()*xx
        exported=[]
        for v in images:
            cs=list(v.list())+[K(0),K(0)]
            exported.append([{"num":rows(h.numerator()),"den":rows(h.denominator())}
                             for h in cs[:2]])
        coverrows.append({"quartic":rows(f),"map_x_y_quadratic_algebra":exported})
    reductions=[]
    for p in [53,59,61]:
        assert E.discriminant().valuation(p)==0
        prime_certificate(p)
        reductions.append({"prime":p,"cardinality":int(E.change_ring(GF(p)).cardinality())})
    assert gcd([z["cardinality"] for z in reductions])==2
    assert all(points[0][1]!="0" for _ in [0])
    return {"schema":"hill1152-first-pair-complete-descent-v1","status":"PASS",
        "selection":"lexicographically first of the frozen15 pairs; no rank-based selection",
        "pair":pair,"pilot_sha256":sha256(PILOT.read_bytes()).hexdigest(),
        "jacobians_sha256":sha256(JACS.read_bytes()).hexdigest(),
        "raw_radicands":[list(map(str,q)) for q in raw],
        "primitive_radicands":[list(map(str,q)) for q in primitive],
        "positive_square_contents":list(map(str,contents)),
        "source_jacobian_to_integral_scale":str(k),
        "integral_jacobian_a_invariants":list(map(str,E.a_invariants())),
        "bad_product":str(bad),"local_witnesses":local,"real_parameter":0,
        "prime_certificates":primes,"torsion_reductions":reductions,
        "pari_ellrank_header":rank_header,"pari_ellrank_points":points,
        "pari_ell2cover_dimension":len(covers),"selmer_basis_covers":coverrows,
        "rational_2torsion_dimension":1,"rank_exact":2,"Sha_2_dimension":0,
        "paired_base_has_rational_point":True,"paired_base_rational_points_infinite":True,
        "explicit_paired_base_point":"NOT_EXTRACTED",
        "software":{"sage":version(),"pari":str(pari("version()"))},
        "limits":{"cpu_seconds":30,"address_space_bytes":4*1024**3,
                  "local_parameter_candidates_per_bad_prime":1000,"ellrank_effort":0,
                  "random_seed":11520001},
        "cpu_seconds":time.process_time()-started,
        "boundary":"PARI-complete2-descent trusted; exact arithmetic replay does not independently prove Selmer completeness or the rank lower bound2. Genus-one paired base existence follows from local solubility, period dividing2 and Sha[2]=0. No explicit paired-base point, new specialized rank or single-quadratic two-gain result."}


if __name__=="__main__":
    from sage.all import GF
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--output",required=True,type=Path)
    a=p.parse_args();result=run()
    with a.output.open("x") as f:json.dump(result,f,indent=2,sort_keys=True);f.write("\n")
    print(json.dumps({k:result[k] for k in ["status","pair","rank_exact","pari_ell2cover_dimension",
                                           "paired_base_has_rational_point","cpu_seconds"]}))
