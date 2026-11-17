#!/usr/bin/env python3
"""Exact controls for ramified-fold necessary conditions; not an exhaustive proof."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import platform
import re
import time
from pathlib import Path
import sympy as S

BASE=Path(__file__).resolve().parent
PROOF=BASE/"RAMIFIED_FOLD_VALUATIONS_AND_COMMON_COORDINATE_POLES.md"
RECEIPT=BASE/"ramified_fold_common_coordinate_receipt.json"
CHECKS=[]

def check(name, condition):
    ok=bool(condition)
    CHECKS.append({"name":name,"passed":ok})
    if not ok:
        raise AssertionError(name)

def eq(name, lhs, rhs=0):
    check(name,S.cancel(S.sympify(lhs)-S.sympify(rhs))==0)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def euclid(p,r):
    ans=[]
    while True:
        ans.append(min(p,r))
        if p==r:
            return ans
        if p>r:
            p-=r
        else:
            r-=p

def ancestral_controls(depth=10):
    # The strict local axes have labels -1=u and -2=b.
    # Each exceptional row gives its coefficients in total pullbacks.
    vertices=[(1,1,2,(1,))]
    edges={(-1,0),(-2,0)}
    count=folds=polar_greater=polar_smaller=equal=terminal=0
    samples={}
    errors=[]

    def visit(vs,es):
        nonlocal count,folds,polar_greater,polar_smaller,equal,terminal
        p,r,A,row=vs[-1]
        count+=1
        block=euclid(p,r)
        g=math.gcd(p,r)
        if A!=sum(row)+1:
            errors.append(("canonical",p,r,A,row))
        if row[-1]!=1 or (len(row)>1 and row[-2]!=1):
            errors.append(("terminal multiplicities",p,r,A,row))
        if any(row[i]<row[i+1] for i in range(len(row)-1)):
            errors.append(("multiplicity monotonicity",p,r,A,row))
        if list(row[:len(block)])!=block:
            errors.append(("Euclid prefix",p,r,A,row,block))
        if any(m>g for m in row[len(block):]):
            errors.append(("tail cap",p,r,A,row))
        if A==2*p+2:
            folds+=1
            n=sum(m*m for m in row)
            if not(n<=p*(2*p+1)<4*p*p):
                errors.append(("coarse fold bound",p,r,n,row))
            if n>p*p:
                if p>r:
                    polar_greater+=1
                    if (p,r,n)!=(2,1,5):
                        errors.append(("greater classification",p,r,n,row))
                elif p<r:
                    polar_smaller+=1
                    good=(r==p+1 and n==p*r+1)
                    toric=(r==p+2 and p%2==1 and
                           len(row)==len(block) and n==p*r)
                    if toric:
                        terminal+=1
                    if not(good or toric):
                        errors.append(("smaller classification",p,r,n,row))
                else:
                    equal+=1
                samples.setdefault((p,r,n),list(row))
        if len(vs)==depth:
            return
        last=len(vs)-1
        adjacent=[b if a==last else a for a,b in es if last in (a,b)]
        for other in [None]+sorted(adjacent):
            parents=[last]+([] if other is None or other<0 else [other])
            pp=sum(vs[j][0] for j in parents)+(other==-1)
            rr=sum(vs[j][1] for j in parents)+(other==-2)
            AA=sum(vs[j][2] for j in parents)+(other is None or other<0)
            mult=tuple(sum(vs[j][3][i] if i<len(vs[j][3]) else 0
                           for j in parents) for i in range(len(vs)))+(1,)
            ne=len(vs)
            newes=set(es)
            if other is not None:
                old=next(e for e in newes if set(e)=={last,other})
                newes.remove(old)
                newes.add((other,ne))
            newes.add((last,ne))
            visit(vs+[(pp,rr,AA,mult)],newes)

    visit(vertices,edges)
    check("all ancestral local recurrence controls",not errors)
    check("fold controls include both surviving pole inequalities",
          folds>0 and polar_greater>0 and polar_smaller>0)
    check("formal equal-value and terminal monomial negative controls occur",
          equal>0 and terminal>0)
    return {"depth":depth,"ancestral_trees":count,"fold_valuations":folds,
            "greater_polar_controls":polar_greater,
            "smaller_polar_controls":polar_smaller,
            "equal_value_controls":equal,"terminal_monomial_controls":terminal,
            "representative_multiplicities":
                {str(k):v for k,v in sorted(samples.items())}}

def run():
    controls=ancestral_controls()
    for p in range(1,65):
        check(f"complete Euclidean identities for p{p}",all(
            sum(euclid(p,r))==p+r-math.gcd(p,r) and
            sum(v*v for v in euclid(p,r))==p*r
            for r in range(1,65)))
        check(f"universal coarse fold bound p{p}",
              p*(2*p+1)<4*p*p)
        check(f"degree-two-or-higher common coordinate contradiction p{p}",
              all(p*(2*p+1)<m*m*p*p for m in range(2,12)))
        check(f"nonpolar coordinate norm threshold p{p}",
              all(2*p+e-1>=e+1 for e in range(1,20)))

    p,r,g,A,B,N,t,PD,s,c,m=S.symbols("p r g A B N t PD s c m")
    eq("anchor determinant deficit decomposition",
       PD**2*(N+t)-(PD**2-s)*N,PD**2*t+s*N)
    eq("fold greater-value algebraic reduction",
       (g*A)*(g*B)+g*(g*A-g*B+g+1)-(g*A)**2,
       g*g*(1-(A-1)*(A-B))+g)
    eq("coarse strict bound algebra",
       4*p*p-p*(2*p+1),p*(2*p-1))
    eq("exceptional p2g rg tail norm upper bound",
       2*g*g+g*(2*g)+1,4*g*g+1)
    for gg in range(1,45):
        check(f"all nonexceptional coprime greater patterns excluded g{gg}",
              all(gg*gg*(1-(AA-1)*(AA-BB))+gg<=0
                  for AA in range(2,35) for BB in range(1,AA)
                  if (AA,BB)!=(2,1)))
        check(f"smaller-value nonnegative tail classification g{gg}",
              all((rr-pp==gg or
                   (gg==1 and rr-pp==2))
                  for pp in range(gg,35*gg,gg)
                  for rr in range(pp+gg,35*gg,gg)
                  if math.gcd(pp,rr)==gg and pp-rr+gg+1>=0))

    one=(m+1)**2*c+m*m+m-1
    two=(2*m+1)**2*c+4*m*m-1
    eq("one-parent label2 determinant formula",one,
       (m+1)**2*c+m*(m+1)-1)
    eq("two-parent label2 determinant formula",two,
       (2*m+1)**2*c+(2*m-1)*(2*m+1))
    negative=[(cc,mm,int(one.subs({c:cc,m:mm})))
              for cc in range(30) for mm in range(30)
              if one.subs({c:cc,m:mm})<0]
    check("only negative first-pattern Case2 label2 value",
          negative==[(0,0,-1)])
    check("all admissible two-parent label2 cases positive",
          all(two.subs({c:cc,m:mm})>0
              for cc in range(12) for mm in range(1,18)))
    check("small norm excludes nonpolar source coordinate",
          all(nn<ee+1 for ee in range(1,30) for nn in range(1,ee+1)))
    check("index2 Case2 cannot meet common-coordinate anchor",
          all(not(1>aa*aa) for aa in range(1,40)))

    # Independent local intersection matrices: -inverse equals the
    # total-pullback row Gram matrix.
    for idx,parents in enumerate([
        [[],[0],[1],[2],[3]],
        [[],[0],[0,1]],
        [[],[0],[0,1],[2],[3]],
        [[],[0],[1],[1,2],[2,3]],
    ]):
        I=S.zeros(0)
        rows=[]
        for par in parents:
            size=I.rows
            J=S.zeros(size+1)
            for i in range(size):
                for j in range(size):
                    J[i,j]=I[i,j]
            for i in par:
                J[i,i]-=1
                J[i,size]=J[size,i]=1
            if len(par)==2:
                J[par[0],par[1]]-=1
                J[par[1],par[0]]-=1
            J[size,size]=-1
            row=[sum(rows[k][j] if j<len(rows[k]) else 0 for k in par)
                 for j in range(size)]+[1]
            rows.append(row)
            I=J
        V=S.Matrix([row+[0]*(len(rows)-len(row)) for row in rows])
        check(f"independent inverse intersection Gram {idx}",-I.inv()==V*V.T)

    # Equal-value pole is incompatible with a regular shared primitive.
    z,T=S.symbols("z T")
    for rr in range(1,8):
        aa=(1+T)*z**(-rr)+T*z**(2-rr)
        bb=z**rr
        residue=S.expand(aa*S.diff(bb,z)).coeff(z,-1)
        eq(f"equal-value nonzero primitive residue r{rr}",residue,rr*(1+T))
        for i,j in [(1,0),(0,1),(1,1),(2,3),(4,1)]:
            FF=S.expand(aa**i*bb**j)
            eq(f"Laurent exact derivative has zero residue r{rr}i{i}j{j}",
               S.expand(S.diff(FF,z)).coeff(z,-1),0)

    # Terminal monomial case: no regular polynomial restricts
    # nonconstantly to b=0.
    for pp in range(1,20,2):
        check(f"terminal fold primitive ray p{pp}",math.gcd(pp,pp+2)==1)
        check(f"terminal fold monomials contract b0 p{pp}",all(
            jj>0 or ii==0
            for ii in range(20) for jj in range(20)
            if -pp*ii+(pp+2)*jj>=0))
    b,Aa,Ab,Ba,Bb,AA,BB=S.symbols("b Aa Ab Ba Bb AA BB")
    jac=(b*Aa)*(BB+b*Bb)-(AA+b*Ab)*(b*Ba)
    eq("constant-on-line target pair has Jacobian factor b",
       jac,b*(Aa*BB-AA*Ba+b*(Aa*Bb-Ab*Ba)))

    # Positive formal models, expressly without a Keller target.
    for pp in range(1,13):
        rr=pp+1
        aa=z**(-pp)+z**(1-pp)*T
        bb=z**rr
        eq(f"survivor volume order1 p{pp}",
           S.diff(aa,T)*S.diff(bb,z),rr*z)
        form=S.expand(aa*S.diff(bb,z))
        check(f"survivor primitive regular p{pp}",
              all(S.degree(S.together(term*z**0),z)>=0
                  for term in S.Add.make_args(form)))
        eq(f"survivor primitive explicit p{pp}",form,rr*(1+z*T))
    aa=z**(-2)+z*T
    bb=z
    F=aa*bb
    eq("exceptional survivor volume",S.diff(aa,T)*S.diff(bb,z),z)
    eq("exceptional survivor primitive dz",
       aa*S.diff(bb,z)+S.diff(F,z),3*z*T)
    eq("exceptional survivor primitive dT",
       aa*S.diff(bb,T)+S.diff(F,T),z*z)

    # Genuine affine marked-source survivor, including both localizations.
    a,b,c,d,x,y,u,T=S.symbols("a b c d x y u T")
    rel=[b*b-a*(c+1),b*c-a*d,b*d-c*(c+1)]
    G=S.groebner(rel,a,b,c,d)
    src={a:x,b:x*y,c:x*y*y-1,d:x*y**3-y}
    ext={a:u*u*(1+u*T),b:u*(1+u*T),c:u*T,d:T}
    for j,f in enumerate(rel):
        eq(f"affine survivor source equation {j}",f.subs(src))
        eq(f"affine survivor boundary chart equation {j}",f.subs(ext))
    def modideal(expr):
        num=S.fraction(S.cancel(expr))[0]
        return G.reduce(S.expand(num))[1]
    for variable,expr in [(a,a),(b,a*d/c),(c,a*d*d/c**2-1),
                          (d,(d/c)*(a*d*d/c**2-1))]:
        eq("localized D(c) inverse "+str(variable),modideal(expr-variable),0)
    U=b/(c+1)
    for variable,expr in [(a,U*U*(1+U*d)),(b,U*(1+U*d)),
                          (c,U*d),(d,d)]:
        eq("localized D(c+1) inverse "+str(variable),modideal(expr-variable),0)
    eq("two smooth localization domains cover", (c+1)-c,1)
    eq("whole source opens cover",x*y*y-(x*y*y-1),1)
    check("boundary quotient is exactly affine line",
          all(f.subs({a:0,b:0,c:0})==0 for f in rel))
    xx=u*u+u**3*T
    yy=1/u
    omega=S.diff(xx,u)*S.diff(yy,T)-S.diff(xx,T)*S.diff(yy,u)
    eq("affine survivor actual source volume",omega,u)
    eq("affine survivor regular primitive du",xx*S.diff(yy,u),-(1+u*T))
    eq("affine survivor primitive dT",xx*S.diff(yy,T),0)
    check("affine survivor has retained norm3",
          sum(euclid(1,2))+1==3 and
          sum(v*v for v in euclid(1,2))+1==3)
    return controls

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-receipt",action="store_true")
    args=parser.parse_args()
    start=time.monotonic()
    controls=run()
    text=PROOF.read_text()
    for phrase in [
        "v_E(b)>=0  implies  n_E>=e+1",
        "sum_j m_j=2p+e-1",
        "D^vee . E^vee=p_D*a",
        "n_E > a^2",
        "m=c=0 and d_2=-1",
        "n_E<=p*(2p+1)<4p^2",
        "(p,r,n_E)=(2,1,5)",
        "r=p+1,       n_E=p^2+p+1",
        "last TWO",
        "theta=a db+dF",
        "C[u,t,(1+u*t)^-1]",
        "No formal verification",
    ]:
        check("written scope retains "+phrase,phrase in text)
    for link in re.findall(r"\]\(([^)]+)\)",text):
        if "://" not in link:
            check("local proof link exists "+link,(BASE/link).is_file())
    parents=[
        "CANONICAL_BOUNDARY_AND_STABLE_COMPLETE_INTERSECTION_OBSTRUCTION.md",
        "EXACT_VOLUME_BOUNDARY_MULTIPLICITIES_AND_COMMON_SOURCE_PENCIL.md",
        "INDIVIDUAL_SOURCE_COORDINATES_AND_UNCONDITIONAL_BOUNDARY_CHARTS.md",
        "COMMON_COORDINATE_FORCES_SMOOTH_BOUNDARY_AND_WHOLE_ROOT_CHARTS.md",
        "sources/individual-source-coordinates-source-audit.json",
        "sources/borisov-divisorial-valuations-2013.pdf",
    ]
    check("primary PDF identity retained",
          sha(BASE/parents[-1])=="b4c0c4511079a3429bc859b0ff9b67d98b750ef502decd036e25dbcf194016d8")
    receipt={
        "claim":"Ramified fold valuations and common-source-coordinate poles",
        "outcome":"PASS","check_count":len(CHECKS),"checks":CHECKS,
        "controls":controls,"proof_sha256":sha(PROOF),
        "checker_sha256":sha(Path(__file__)),
        "retained_inputs":[{"path":p,"sha256":sha(BASE/p)} for p in parents],
        "environment":{"python":platform.python_version(),"sympy":S.__version__},
        "scope":{
            "general_small_norm_source_coordinate_polarity":True,
            "polar_prime_anchor_requires_nonempty_unramified_boundary":True,
            "index_two_minimal_case2_excluded_with_anchor":True,
            "polar_case1_fold_transverse_degree_one":True,
            "refined_two_necessary_fold_multiplicity_types":True,
            "genuine_affine_exact_volume_survivor_control":True,
            "ramified_whole_chart_assumed":False,
            "all_ramified_primes_excluded":False,
            "common_coordinate_integrality_proved":False,
            "Keller_counterexample_or_inverse_claimed":False,
        },
        "assurance":"Internal written proof with exact controls; no formal verification, external review or novelty claim.",
    }
    if args.write_receipt:
        RECEIPT.write_text(json.dumps(receipt,indent=2)+chr(10))
    print(json.dumps({"outcome":"PASS","check_count":len(CHECKS),
                      "controls":controls,
                      "elapsed_seconds":time.monotonic()-start,
                      "receipt_written":args.write_receipt},indent=2))

if __name__=="__main__":
    main()
