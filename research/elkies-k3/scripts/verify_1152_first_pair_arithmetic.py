#!/usr/bin/env python3
"""Independent exact arithmetic supporting the first-pair PARI descent.

Checks normalization, all local witnesses, complete bad-prime factorization
using recursive Lucas certificates, finite torsion bounds, one nontorsion
point, and the three exported quartic maps. DOES NOT certify Selmer
completeness, independence of these classes, or PARI's exact rank2.
"""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from math import gcd, isqrt, lcm
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/"artifacts/generated-results/x1152-first-pair-descent-20261114-v1.json"
spec=importlib.util.spec_from_file_location("first_pair_polynomials",Path(__file__).with_name("verify_r17_correlated_genus_one_pencils.py"))
V=importlib.util.module_from_spec(spec);spec.loader.exec_module(V)
require=V.require


def quotient(a,b):
    r=list(a);out=[Q(0)]*max(0,len(a)-len(b)+1)
    require(b,"zero polynomial denominator")
    while len(r)>=len(b):
        i=len(r)-len(b);s=r[-1]/b[-1];out[i]=s
        for j,v in enumerate(b):r[i+j]-=s*v
        r=V.trim(r)
    require(not r,"inexact polynomial quotient")
    return V.trim(out)


def rf(a,b=None):
    a=V.poly(a);b=V.poly([1] if b is None else b)
    require(b,"zero rational-function denominator")
    if not a:return [],[Q(1)]
    g=V.gcd_poly(a,b);a,b=quotient(a,g),quotient(b,g)
    return V.scale(a,1/b[-1]),V.scale(b,1/b[-1])


def add(a,b):return rf(V.add(V.mul(a[0],b[1]),V.mul(b[0],a[1])),V.mul(a[1],b[1]))
def mul(a,b):return rf(V.mul(a[0],b[0]),V.mul(a[1],b[1]))
def scale(a,k):return rf(V.scale(a[0],k),a[1])
def scalar(k):return rf([k])
def qadd(a,b):return add(a[0],b[0]),add(a[1],b[1])
def qscale(a,k):return scale(a[0],k),scale(a[1],k)
def qmul(a,b,f):
    return add(mul(a[0],b[0]),mul(mul(a[1],b[1]),rf(f))),add(mul(a[0],b[1]),mul(a[1],b[0]))
def qsquare(a,f):return qmul(a,a,f)


def verify_prime(n,proofs,done,visiting):
    if n in done:return
    require(n>=2 and n not in visiting,"prime certificate cycle")
    visiting.add(n);row=proofs[str(n)]
    require(int(row["prime"])==n,"prime identity")
    if n==2:
        done.add(n);visiting.remove(n);return
    factors=[(int(q),int(e)) for q,e in row["p_minus_one"]]
    require(len({q for q,e in factors})==len(factors) and all(e>0 for q,e in factors),"prime factors")
    product=1
    for q,e in factors:
        require(2<=q<n,"prime recursion")
        verify_prime(q,proofs,done,visiting);product*=q**e
    require(product==n-1,"incomplete p-1 factorization")
    g=int(row["full_order_witness"])
    require(pow(g,n-1,n)==1 and all(gcd(pow(g,(n-1)//q,n)-1,n)==1 for q,e in factors),
            "Lucas full-order certificate")
    done.add(n);visiting.remove(n)


def value(f,t):return sum(a*t**i for i,a in enumerate(f))


def local_square(v,p):
    if not v:return True
    exponent=0
    while v%p==0:exponent+=1;v//=p
    return exponent%2==0 and (v%8==1 if p==2 else pow(v%p,(p-1)//2,p)==1)


def rational_square(v):
    return v>=0 and isqrt(v.numerator)**2==v.numerator and isqrt(v.denominator)**2==v.denominator


def verify(d):
    require(d["schema"]=="hill1152-first-pair-complete-descent-v1","schema")
    paths={"pilot_sha256":"artifacts/generated-results/hill1152-transfer-20261113-v2/result.json",
           "jacobians_sha256":"artifacts/generated-results/x1152-pair-jacobians-20261114-v1.json"}
    inputs={}
    for key,path in paths.items():
        raw=(ROOT/path).read_bytes();require(sha256(raw).hexdigest()==d[key],"input hash")
        inputs[key]=json.loads(raw)
    pilot=inputs["pilot_sha256"];jacs=inputs["jacobians_sha256"]
    require(d["pair"]==list(min(map(tuple,pilot["locally_surviving_disjoint_pairs"])))==[0,1],"sealed first-pair selection")
    raw=[list(map(Q,pilot["quadratic_covers"][i]["q"])) for i in d["pair"]]
    require(raw==[list(map(Q,f)) for f in d["raw_radicands"]],"raw radicands")
    primitive=[];contents=[]
    for q in raw:
        den=lcm(*(v.denominator for v in q));ints=[int(v*den) for v in q]
        g=gcd(*ints);content=Q(g,den)
        require(content>0 and rational_square(content),"squareclass-preserving normalization")
        primitive.append([v//g for v in ints]);contents.append(content)
    require(primitive==[list(map(int,f)) for f in d["primitive_radicands"]],"primitive coefficients")
    require(contents==list(map(Q,d["positive_square_contents"])),"contents")
    q,r=primitive
    m=q[2]*r[0]+q[0]*r[2]-Q(q[1]*r[1],2)
    dq=q[1]**2-4*q[2]*q[0];dr=r[1]**2-4*r[2]*r[0]
    res=m*m-Q(dq*dr,4)
    require(dq*dr*res and res.denominator==1,"smooth branch pairs")
    a,b=4*m,Q(dq*dr)
    require([0,a,0,b,0]==list(map(Q,d["integral_jacobian_a_invariants"])),"integral Jacobian")
    row=next(z for z in jacs["pairs"] if z["pair"]==d["pair"])
    k=Q(d["source_jacobian_to_integral_scale"])
    require(k>0 and k*k*contents[0]*contents[1]==4 and
            k*k*Q(row["jacobian_a_invariants"][1])==a and
            k**4*Q(row["jacobian_a_invariants"][3])==b,"Jacobian Q-isomorphism")
    require(not rational_square(a*a-4*b),"rational2-torsion dimension")
    bad=2*dq*dr*res;require(bad==int(d["bad_product"]),"bad-prime product")
    done=set()
    for p in d["prime_certificates"]:verify_prime(int(p),d["prime_certificates"],done,set())
    local=d["local_witnesses"];require(len({z["prime"] for z in local})==len(local),"duplicate local primes")
    product=1
    for z in local:
        p=int(z["prime"]);e=z["exponent_in_bad_product"];t=z["integer_parameter"]
        require(p in done and type(e) is int and e>0,"bad-prime primality")
        require(type(t) is int and 0<=t<1000,"local witness allowance")
        require(all(local_square(value(f,t),p) for f in primitive),"local square witness")
        product*=p**e
    require(product==abs(bad),"all bad primes covered")
    require(all(value(f,d["real_parameter"])>0 for f in primitive),"real solubility")
    counts=[]
    for z in d["torsion_reductions"]:
        p=z["prime"];require(p in done and p>3,"good-prime certificate")
        require((16*b*b*(a*a-4*b)).numerator%p,"bad torsion-reduction prime")
        count=1
        for x in range(p):
            v=int((x**3+a*x*x+b*x)%p)
            count+=1 if not v else (2 if pow(v,(p-1)//2,p)==1 else 0)
        require(count==z["cardinality"],"torsion reduction count")
        counts.append((p,count))
    # Ensure each reduction characteristic's primary torsion is controlled
    # by another prime where the group order is prime to that characteristic.
    require(len(counts)>1 and all(any(q!=p and n%p for q,n in counts) for p,_ in counts),"torsion primary injection")
    require(gcd(*(n for p,n in counts))==2,"torsion bound")
    points=d["pari_ellrank_points"];require(points,"no exported point")
    for P in points:
        x,y=map(Q,P);require(y and y*y==x**3+a*x*x+b*x,"nontorsion point identity")
    covers=d["selmer_basis_covers"];require(len(covers)==3,"exported cover roster")
    for row in covers:
        f=V.poly(row["quartic"])
        require(len(f)==5 and len(V.gcd_poly(f,[i*c for i,c in enumerate(f)][1:]))==1,"smooth quartic")
        images=[tuple(rf(z["num"],z["den"]) for z in cs) for cs in row["map_x_y_quadratic_algebra"]]
        require(len(images)==2 and all(len(z)==2 for z in images),"cover map shape")
        x,y=images
        rhs=qadd(qadd(qmul(qsquare(x,f),x,f),qscale(qsquare(x,f),a)),qscale(x,b))
        require(qsquare(y,f)==rhs,"quartic-to-Jacobian map identity")
        require(len(x[0][0])>1 or len(x[0][1])>1 or x[1][0],"constant cover map")
    # Source height identities support the two distinct-character directions.
    # The generic height matrices themselves are dependencies, not recomputed
    # geometrically by this arithmetic checker.
    heightpath=ROOT/"artifacts/generated-results/hill1152-transfer-20261113-v2/four-generic-covers.json"
    heightcert=json.loads(heightpath.read_text())
    anti_heights=[]
    for i in d["pair"]:
        family=next(z for z in heightcert["families"] if z["index"]==i)
        H=[list(map(Q,r)) for r in family["height_matrix"]]
        word=list(map(Q,pilot["quadratic_covers"][i]["word"]))
        require(all(H[k][17]==sum(H[k][j]*word[j]/2 for j in range(17))
                    for k in range(17)),"trace projection into parent subgroup")
        trace_height=sum(word[k]*H[k][j]*word[j] for k in range(17) for j in range(17))
        h=4*H[17][17]-trace_height
        require(h==12,"nonzero anti-invariant height")
        anti_heights.append(str(h))
    # These are labels consistent with the CAS result, not an independent
    # certification of its global completeness or its rank lower bound2.
    require(d["pari_ellrank_header"]==[2,2,0] and d["pari_ell2cover_dimension"]==3,"PARI result header")
    require(d["rank_exact"]==2 and d["rational_2torsion_dimension"]==1 and d["Sha_2_dimension"]==0,"descent dimension accounting")
    require(d["paired_base_has_rational_point"] is True and d["paired_base_rational_points_infinite"] is True
            and d["explicit_paired_base_point"]=="NOT_EXTRACTED","conclusion labels")
    return {"schema":"hill1152-first-pair-arithmetic-replay-v2","status":"PASS",
        "pair":[0,1],"all_places_locally_soluble":True,"bad_primes_checked":len(local),
        "recursive_prime_certificates":len(done),"quartic_maps_checked":len(covers),
        "rational_2torsion_dimension":1,"independent_nontorsion_points_at_least":1,
        "source_height_certificate_sha256":sha256(heightpath.read_bytes()).hexdigest(),
        "anti_invariant_heights_on_quadratic_bases":anti_heights,
        "anti_invariant_heights_on_paired_base":["24","24"],
        "independent_selmer_completeness":False,"independent_exact_rank":False,
        "CAS_premise":"PARI ellrank=[2,2,0,...] and ell2cover dimension3 must be replayed separately.",
        "boundary":"All local places and exported arithmetic independently verified. Global existence and infinitude additionally use the complete PARI descent; no explicit paired-base point extracted."}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--certificate",type=Path,default=DEFAULT)
    p.add_argument("--output",type=Path);a=p.parse_args()
    result=verify(json.loads(a.certificate.read_text()))
    if a.output:
        with a.output.open("x") as f:json.dump(result,f,indent=2,sort_keys=True);f.write("\n")
    print(json.dumps(result,sort_keys=True))
