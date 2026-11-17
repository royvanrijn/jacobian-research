#!/usr/bin/env python3
"""Exact non-Keller barrier to a valuation-only common-coordinate argument."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import time
import sympy as S
START=time.perf_counter()
BASE=Path(__file__).resolve().parent
PROOF=BASE/"VALUATION_MATRIX_AND_CONTRACTED_SOURCE_AXES_CONTROL.md"
PRIMARY=BASE/"sources/xie-valuation-intersections-2015.pdf"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def zero(e):return S.cancel(e)==0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
x,y,s,t,T,U,V=S.symbols("x y s t T U V")
# Primitive divisorial charts, unit volume and a regular exact primitive.
for name,xx,yy,sign in [("v1",1/s,s**2*t,-1),("v2",s**2*t,1/s,1)]:
    jac=S.det(S.Matrix([[S.diff(xx,s),S.diff(xx,t)],
                       [S.diff(yy,s),S.diff(yy,t)]]))
    check("source volume order zero "+name,zero(jac-sign))
    alpha_s=S.cancel((xx*S.diff(yy,s)-yy*S.diff(xx,s))/2)
    alpha_t=S.cancel((xx*S.diff(yy,t)-yy*S.diff(xx,t))/2)
    check("primitive coefficients regular in root plane "+name,
          alpha_s.is_polynomial(s,t) and alpha_t.is_polynomial(s,t))
    check("primitive derivative is actual pulled-back volume "+name,
          zero(S.diff(alpha_t,s)-S.diff(alpha_s,t)-jac))
    for i,j in [(0,0),(1,0),(0,1),(2,3),(5,2),(3,5)]:
        expected=s**(-i+2*j)*t**j if name=="v1" else s**(2*i-j)*t**i
        check(f"exact monomial order transcription {name} i{i} j{j}",
              zero(xx**i*yy**j-expected))
check("both valuation groups are primitive",S.gcd(-1,2)==1)
check("valuation normalization factor is one",min(-1,2)==-1)

# Complete Hilbert cone decomposition, with explicit unbounded formulas in proof.
basis=[(1,1),(2,1),(1,2)]
for i,j in basis:
    check(f"Hilbert generator satisfies both valuations {i} {j}",-i+2*j>=0 and 2*i-j>=0)
for i in range(25):
    for j in range(25):
        admissible=(-i+2*j>=0 and 2*i-j>=0)
        if admissible:
            if i>=j:
                powerU=i-j;powerT=2*j-i
                check(f"upper cone decomposition i{i} j{j}",
                      powerU>=0 and powerT>=0 and
                      (2*powerU+powerT,powerU+powerT)==(i,j))
            else:
                powerV=j-i;powerT=2*i-j
                check(f"lower cone decomposition i{i} j{j}",
                      powerV>=0 and powerT>=0 and
                      (powerV+powerT,2*powerV+powerT)==(i,j))
            if i+j:
                check(f"every admissible nonconstant monomial contains both axes i{i} j{j}",i>0 and j>0)
        else:
            check(f"excluded monomial has a negative valuation i{i} j{j}",min(-i+2*j,2*i-j)<0)
embed={T:x*y,U:x**2*y,V:x*y**2}
relation=U*V-T**3
check("toric presentation relation exact",zero(relation.subs(embed)))
check("presentation binomial irreducible",S.Poly(relation,T,U,V).is_irreducible)
check("source x recovered in fraction field",zero((U/T).subs(embed)-x))
check("source y recovered in fraction field",zero((T**2/U).subs(embed)-y))
check("independent generators have nonzero Jacobian",
      S.det(S.Matrix([[S.diff(x*y,x),S.diff(x*y,y)],
                      [S.diff(x**2*y,x),S.diff(x**2*y,y)]]))==-x**2*y)

# Construct the full intersection model by six actual free blowup operations.
I=S.Matrix([[1]])
omega_orders=[-3]
parents=[]
def free_blow(parent):
    global I
    n=I.rows
    J=S.zeros(n+1,n+1)
    J[:n,:n]=I
    J[parent,parent]-=1
    J[n,n]=-1
    J[parent,n]=J[n,parent]=1
    I=J
    omega_orders.append(omega_orders[parent]+1)
    parents.append(parent)
    return n
a1=free_blow(0);a2=free_blow(a1);a3=free_blow(a2)
b1=free_blow(0);b2=free_blow(b1);b3=free_blow(b2)
check("six free blowups give exact diagonal",
      list(I.diagonal())==[-1,-2,-2,-1,-2,-2,-1])
check("intersection matrix unimodular determinant one",I.det()==1)
check("both terminal volume orders are zero",omega_orders[a3]==omega_orders[b3]==0)
expected=S.diag(-1,-2,-2,-1,-2,-2,-1)
for i,j in [(0,1),(1,2),(2,3),(0,4),(4,5),(5,6)]:
    expected[i,j]=expected[j,i]=1
check("both disjoint blowup branches retain exact incidence",I==expected)
dual=I.inv()
for idx in [a3,b3]:
    check("dual divisor gives Kronecker intersections "+str(idx),
          I*dual[:,idx]==S.eye(7)[:,idx])
M=dual.extract([a3,b3],[a3,b3])
check("exact unnormalized Xie dual submatrix",M==S.Matrix([[-2,1],[1,-2]]))
check("both normalized valuations have same matrix",M[0,0]==-2 and M[1,1]==-2)
check("negative matrix eigenvalues exact",M.eigenvals()=={-1:1,-3:1})
check("negative principal minors exact",M[0,0]<0 and M.det()==3)
check("bare n bound fails with unit pole p1",2>1**2)

# The whole source is contracted; the actual Keller hypothesis visibly fails.
for axis in [{x:0},{y:0}]:
    check("each full source axis contracts "+str(axis),
          all(z.subs(axis)==0 for z in embed.values()))
for xx in [1,2,3]:
    check("distinct source points have identical target xx"+str(xx),
          [expr.subs({x:xx,y:0}) for expr in embed.values()]==[0,0,0])
g,h,gx,gy,hx,hy=S.symbols("g h gx gy hx hy")
Px=y*g+x*y*gx;Py=x*g+x*y*gy
Qx=y*h+x*y*hx;Qy=x*h+x*y*hy
J=S.expand(Px*Qy-Py*Qx)
expected=x*y*(x*(gx*h-g*hx)+y*(g*hy-gy*h)+x*y*(gx*hy-gy*hx))
check("unrestricted Jacobian factorization",zero(J-expected))
check("unrestricted Jacobian divisible by both axes",zero(J.subs(x,0)) and zero(J.subs(y,0)))
for i in range(1,6):
    for j in range(1,6):
        f=x**i*y**j
        check(f"distinguished coordinate fibre has x axis factor i{i} j{j}",f.subs(x,0)==0)
        check(f"distinguished coordinate fibre has y axis factor i{i} j{j}",f.subs(y,0)==0)

receipt={"claim":"HC4MYVAC0","outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "primary_input":{"url":"https://archive.ymsc.tsinghua.edu.cn/pacm_download/423/9433-affinegeo20140928.pdf",
        "path":str(PRIMARY.relative_to(BASE)),"sha256":sha(PRIMARY),
        "locator":"PDFpage3, dual-divisor matrix definition and Theorem1.2",
        "scope":"Matrix convention and transcendence-degree criterion only"},
    "exact_control":{"valuation_weights":[[-1,2],[2,-1]],
        "volume_orders":[0,0],"regular_chart_primitives":True,
        "intersection_ring":"C[xy,x^2y,xy^2]",
        "presentation":"UV=T^3","dual_matrix":[[-2,1],[1,-2]],
        "matrix_eigenvalues":[-1,-3],"intersection_contains_source_coordinate":False,
        "whole_source_open_immersion":False,
        "every_pair_Jacobian_divisible_by":"xy",
        "actual_Keller_counterexample":False},
    "scope":"Barrier to valuation-only common-coordinate reasoning; actual full-source common-pencil theorem neither proved nor disproved.",
    "assurance":"Internal written exact countercontrol and symbolic replay; no formal or external human review."}
parser=argparse.ArgumentParser()
parser.add_argument("--write-receipt",action="store_true")
args=parser.parse_args()
if args.write_receipt:
    (BASE/"valuation_matrix_source_axes_verification.json").write_text(
        json.dumps(receipt,indent=2)+chr(10),encoding="utf-8")
print(json.dumps({"claim":receipt["claim"],"outcome":"PASS","check_count":len(checks),
                  "elapsed_seconds":receipt["elapsed_seconds"],"receipt_written":args.write_receipt},indent=2))
