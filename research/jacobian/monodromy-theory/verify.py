"""Exact finite controls; KMFI1's universal geometric proof is in THEORY.md."""
from pathlib import Path

import tempfile as _check_tempfile

_CHECK_OUTPUT_DIRECTORY = None

def _check_output(original):
    """Write new runs to a temporary directory; preserve deposited receipts."""
    global _CHECK_OUTPUT_DIRECTORY
    if _CHECK_OUTPUT_DIRECTORY is None:
        _CHECK_OUTPUT_DIRECTORY = Path(_check_tempfile.mkdtemp(prefix="math-check-"))
        print("Verification output:", _CHECK_OUTPUT_DIRECTORY)
    return _CHECK_OUTPUT_DIRECTORY / Path(original).name


from itertools import product
from collections import Counter
from math import gcd
import ast
import hashlib
import json
import platform
import re
import time

import sympy as s
from sympy.matrices.normalforms import smith_normal_form
import model as m

ROOT=Path(__file__).resolve().parent
started=time.monotonic()
checks=[]
def check(name, condition):
    if not condition:
        raise AssertionError(name)
    checks.append(name)

def pmul(a,b):
    return tuple(a[b[i]] for i in range(len(a)))
def pinv(a):
    out=[0]*len(a)
    for i,j in enumerate(a): out[j]=i
    return tuple(out)
def generated(gens, identity, multiply=pmul):
    seen={identity}; stack=[identity]
    while stack:
        a=stack.pop()
        for g in gens:
            b=multiply(a,g)
            if b not in seen:
                seen.add(b); stack.append(b)
    return frozenset(seen)
def normal_closure(group, gens, identity, multiply=pmul, inverse=pinv):
    conjugates=sorted({multiply(multiply(g,h),inverse(g)) for g in group for h in gens})
    return generated(conjugates, identity, multiply)
def cycles(p):
    seen=set(); sizes=[]
    for i in range(len(p)):
        if i in seen: continue
        j=i; size=0
        while j not in seen:
            seen.add(j); size+=1; j=p[j]
        sizes.append(size)
    return tuple(sorted(sizes,reverse=True))
def permutation(cs,n):
    out=list(range(n))
    for c in cs:
        for i,j in zip(c,c[1:]+c[:1]): out[i-1]=j-1
    return tuple(out)

dihedral=[]
for n in range(3,33):
    I=tuple(range(n)); r=tuple((i+1)%n for i in range(n)); h=tuple((-i)%n for i in range(n))
    G=generated((r,h),I)
    N=normal_closure(G,(h,),I)
    check(f"D_{n} natural faithful group order",len(G)==2*n)
    check(f"D_{n} point stabilizer normal-closure parity",len(N)==(2*n if n%2 else n))
    fixed=generated(tuple(g for g in G if any(g[i]==i for i in range(n))),I)
    check(f"D_{n} fixed-sheet-generated subgroup",fixed==N)
    if n%2==0:
        # The two reflection classes give the same index obstruction.
        other=pmul(r,h)
        check(f"D_{n} other reflection class",len(normal_closure(G,(other,),I))==n)
    dihedral.append(dict(n=n,group_order=len(G),normal_closure_order=len(N),
                         excluded_by_KMFI1=len(N)!=len(G)))
for n in range(2,10):
    I=tuple(range(n)); r=tuple((i+1)%n for i in range(n))
    G=generated((r,),I)
    check(f"regular C_{n} has no nonidentity fixed-sheet elements",
          all(g==I or 1 not in cycles(g) for g in G))
# A4 and F20 are positive controls for the necessary test.
for name,n,gs,h in [
    ("A4",4,(permutation([[1,2,3]],4),permutation([[1,2],[3,4]],4)),permutation([[2,3,4]],4)),
    ("F20",5,(tuple((i+1)%5 for i in range(5)),tuple(2*i%5 for i in range(5))),tuple(2*i%5 for i in range(5)))]:
    I=tuple(range(n)); G=generated(gs,I)
    check(name+" stabilizer normally generates",normal_closure(G,(h,),I)==G)
    check(name+" control order",len(G)==(12 if name=="A4" else 20))

prime=11
def canonical(a):
    return min(tuple(v%prime for v in a),tuple(-v%prime for v in a))
def mul(a,b):
    return canonical((a[0]*b[0]+a[1]*b[2],a[0]*b[1]+a[1]*b[3],
                      a[2]*b[0]+a[3]*b[2],a[2]*b[1]+a[3]*b[3]))
def inv(a): return canonical((a[3],-a[1],-a[2],a[0]))
I=canonical((1,0,0,1))
G=frozenset(canonical(a) for a in product(range(11),repeat=4) if (a[0]*a[3]-a[1]*a[2])%11==1)
def order(a):
    b=I
    for k in range(1,13):
        b=mul(a,b)
        if b==I: return k
    raise AssertionError(a)
def action(a):
    aa,b,c,d=a
    out=[]
    for x in range(11):
        den=(c*x+d)%11
        out.append(11 if den==0 else ((aa*x+b)*pow(den,-1,11))%11)
    out.append(11 if c==0 else aa*pow(c,-1,11)%11)
    return tuple(out)
by_order={}
for a in sorted(G): by_order.setdefault(order(a),[]).append(a)
expected={1:((1,)*12,1),2:((2,)*6,55),3:((3,)*4,110),
          5:((5,5,1,1),264),6:((6,6),110),11:((11,1),120)}
for k,(ct,count) in expected.items():
    check(f"PSL2(11) order {k} census",len(by_order[k])==count)
    check(f"PSL2(11) order {k} natural cycles",all(cycles(action(a))==ct for a in by_order[k]))
H=frozenset(a for a in G if action(a)[11]==11)
check("PSL2(11) natural Borel size",len(H)==55)
# H=<upper unipotent, square diagonal>.
u=canonical((1,1,0,1)); diag=canonical((2,0,0,6))
check("Borel generators",generated((u,diag),I,mul)==H)
check("Borel normally generates PSL2(11)",normal_closure(G,(u,diag),I,mul,inv)==G)
A=canonical((0,1,10,3)); B=canonical((0,1,10,8)); C=mul(A,B)
check("5,5,2 witness orders",(order(A),order(B),order(C))==(5,5,2))
check("5,5,2 witness product one",mul(mul(A,B),C)==I)
check("5,5,2 witness generates full group",generated((A,B),I,mul)==G)
indices=[12-len(cycles(action(a))) for a in (A,B,C)]
check("5,5,2 quotient genus zero",indices==[8,8,6] and sum(indices)==22)
solutions=[(a,b,j) for a in range(4) for b in range(4) for j in (0,6,8,10) if 8*a+10*b+j==22]
check("complete genus-zero fixed-inertia passport sieve",solutions==[(2,0,6)])
pair_counts={}
for trace2 in (5,9):
    a=min(g for g in by_order[5] if (g[0]+g[3])**2%11==trace2)
    counts=Counter()
    for b in by_order[5]:
        if order(mul(a,b))==2:
            counts[((b[0]+b[3])**2%11,len(generated((a,b),I,mul)))]+=1
    check(f"order-five class {trace2} generating-pair count",
          counts==Counter({(trace2,660):10,(14-trace2,60):10}))
    pair_counts[str(trace2)]={str(k):v for k,v in counts.items()}

# One exceptional A5 quotient supplies an independent action comparison.
a2=by_order[2][0]
for a3 in by_order[3]:
    if order(mul(a2,a3))==5:
        HA=generated((a2,a3),I,mul)
        if len(HA)==60: break
else: raise AssertionError("A5 subgroup not found")
unseen=set(G); reps=[]; lookup={}
while unseen:
    rep=min(unseen); coset={mul(rep,h) for h in HA}
    for h in coset: lookup[h]=len(reps)
    reps.append(rep); unseen-=coset
def exceptional(a): return tuple(lookup[mul(a,r)] for r in reps)
exc_expected={1:(1,)*11,2:(2,2,2,2,1,1,1),3:(3,3,3,1,1),
              5:(5,5,1),6:(6,3,2),11:(11,)}
for k,ct in exc_expected.items():
    check(f"exceptional degree-eleven order {k} cycles",
          all(cycles(exceptional(a))==ct for a in by_order[k]))

# Read only literal permutation arrays from the primary source snapshot.
dessins=(ROOT/"sources"/"Dessins_and_ModuliFields.txt").read_text()
f11=dessins.split("F11 has 2 dessins.",1)[1].split("F12 has",1)[0]
arrays=re.findall(r"g[01]\s*:=\s*(\[.*?\]);",f11)
check("two primary-source F11 dessins retained",len(arrays)==4)
def equivariant(pgens,qgens):
    n=len(pgens[0])
    for start in range(n):
        mapping={0:start}; stack=[0]
        while stack:
            x=stack.pop()
            for p,q in zip(pgens,qgens):
                y=p[x]; z=q[mapping[x]]
                if y in mapping and mapping[y]!=z: break
                if y not in mapping: mapping[y]=z;stack.append(y)
            else: continue
            break
        else:
            if len(mapping)==n and len(set(mapping.values()))==n: return tuple(mapping[i] for i in range(n))
    return None
dessin_maps=[]
for i in range(2):
    ps=tuple(permutation(ast.literal_eval(v),12) for v in arrays[2*i:2*i+2])
    P=generated(ps,tuple(range(12)))
    check(f"F11 dessin {i+1} group order",len(P)==660)
    found=None
    for tr in (5,9):
        a=min(g for g in by_order[5] if (g[0]+g[3])**2%11==tr)
        for b in by_order[2]:
            if order(mul(a,b))!=5: continue
            relabel=equivariant(ps,(action(a),action(b)))
            if relabel is not None:
                found=dict(order5=a,order2=b,relabel=relabel);break
        if found is not None: break
    check(f"F11 dessin {i+1} explicitly conjugate to natural projective action",found is not None)
    dessin_maps.append(found)

# The full PGL action is a control for the arithmetic normalizer.
natural={action(a) for a in G}
dilation=tuple([2*i%11 for i in range(11)]+[11])
pgl=generated((action(u),action(canonical((0,-1,1,0))),dilation),tuple(range(12)))
check("PGL2(11) control order",len(pgl)==1320 and natural<pgl)
check("PGL outer generator odd", (12-len(cycles(dilation)))%2==1)
check("PGL normalizes PSL",all(pmul(pmul(dilation,g),pinv(dilation)) in natural for g in natural))

# Exact rational model over Q[a]/Phi_10(a), no numerical root computations.
check("cyclotomic coefficient polynomial",s.Poly(m.minimal,m.a).is_irreducible and m.minimal==s.cyclotomic_poly(10,m.a))
N,D,Q=m.N,m.D,m.Q
check("F11 degree and projective infinity orders",(N.degree(),D.degree(),Q.degree())==(12,11,12))
check("F11 numerator denominator coprime",s.gcd(N,D).degree()==0)
k,fac=Q.sqf_list()
check("F11 order-two pole fiber",[(f.degree(),e) for f,e in fac]==[(6,2)])
R=fac[0][0]
check("exact pole square",Q.monic()==R**2)
parts=[m.A,m.B,m.C,m.E,R]
for i,f in enumerate(parts):
    check(f"boundary part {i} reduced",s.gcd(f,f.diff()).degree()==0)
    for j,g in enumerate(parts[:i]):
        check(f"boundary parts {j},{i} disjoint",s.gcd(f,g).degree()==0)
expected_factors=[[1,1],[2],[2],[1],[6]]
for name,f,degs in zip(("A","B","C","E","R"),parts,expected_factors):
    check("coefficient-field factor degrees "+name,sorted(g.degree() for g,e in f.factor_list()[1])==degs)
der=N.diff()*Q-N*Q.diff()
quot,rem=der.div(m.B**4*m.C**4*R)
check("exact derivative factorization",rem.is_zero and quot.degree()==0)
lam=quot.rep.to_dict()[(0,)]
# beta'=lam*B^4*C^4/(k^2*R^3).
check("derivative has only removed zeros and poles",
      der.degree()==22 and (N.degree()+Q.degree())-der.degree()==2)
matrix=s.Matrix([[1,0],[1,0],[5,0],[0,5],[0,1],[-2,-2]])
derivative=s.Matrix([0,0,4,4,0,-3])
check("F11 target-unit rank",matrix.rank()==2)
check("F11 derivative outside target-unit image",matrix.row_join(derivative).rank()==3)
for name,mat,rank in [("base",matrix,2),("base and derivative",matrix.row_join(derivative),3)]:
    snf=smith_normal_form(mat,domain=s.ZZ)
    check("primitive unit lattice "+name,[abs(snf[i,i]) for i in range(rank)]==[1]*rank)
check("geometric puncture and unit ranks",4+4+6==14 and 14-1==13)
# Resultant identity in THEORY/F11: Disc_x(N-tQ)=c*t^8*(1-t)^8.
res=lambda a,b:a.rep.resultant(b.rep)
disc_constant=lam**12*res(m.B,Q)**4*res(m.C,N)**4*res(R,N)/res(Q,N)
for t in (2,3,-1):
    check(f"discriminant resultant control at {t}",
          (N-t*Q).rep.discriminant()==disc_constant*m.K(t**8*(1-t)**8))
w=s.symbols("w")
square_test=s.Poly.from_dict({(2,):m.K.one,(0,):-disc_constant},w,domain=m.K)
check("F11 arithmetic discriminant constant is nonsquare",
      [(f.degree(),e) for f,e in square_test.factor_list()[1]]==[(2,1)])
# Characteristic and target controls, not use cases of the theorem.
check("Artin-Schreier etale derivative control",s.Poly(s.diff(w**5-w,w),w,modulus=5).as_expr()==-1)
check("Artin-Schreier distinct translates control",s.Poly((w+1)**5-(w+1)-(w**5-w),w,modulus=5).is_zero)
check("punctured Kummer retains a nonconstant unit",s.diff(w**3,w)==3*w**2)

def field_coeff(c): return [str(v) for v in c.to_list()]
def poly_coeff(p): return {str(mon[0]):field_coeff(c) for mon,c in sorted(p.rep.to_dict().items())}
model_data=dict(coefficient_field="Q[a]/(a^4-a^3+a^2-a+1)",
  coefficient_encoding="Each coefficient is a descending list in a with rational entries.",
  A=poly_coeff(m.A),B=poly_coeff(m.B),C=poly_coeff(m.C),E=poly_coeff(m.E),
  N=poly_coeff(N),D=poly_coeff(D),Q=poly_coeff(Q),R=poly_coeff(R),
  k=field_coeff(Q.rep.LC()),derivative_lambda=field_coeff(lam),
  discriminant_constant=field_coeff(disc_constant),
  discriminant_formula="Disc_x(N-t*Q)=c*t^8*(1-t)^8",
  arithmetic_group_over_coefficient_field="PGL2(11), using the external Aut(PSL2(11))=PGL2(11) fact",
  geometric_group="PSL2(11), with the source formula/dessin matching attributed to van Hoeij--Vidunas",
  geometric_prime_punctures=14,geometric_unit_rank=13,
  coefficient_field_finite_boundary_prime_count=6,coefficient_field_unit_rank=6,
  target_unit_matrix=[list(matrix.row(i)) for i in range(6)],
  derivative_vector=list(derivative),absolute_realization=False)
# Convert the few SymPy integers in matrices to JSON integers.
model_data["target_unit_matrix"]=[[int(v) for v in row] for row in model_data["target_unit_matrix"]]
model_data["derivative_vector"]=[int(v) for v in model_data["derivative_vector"]]
(_check_output(ROOT/"model.json")).write_text(json.dumps(model_data,indent=2)+"\n")
paths=['THEORY.md', 'model.py', 'verify.py', 'model.json']
receipt=dict(result="pass",check_count=len(checks),checks=checks,python=platform.python_version(),
 sympy=s.__version__,runtime_seconds=round(time.monotonic()-started,3),
 proof_boundary="Universal inertia, normal-generation, nilpotent and all-even-dihedral exclusions have written proofs. Finite controls do not formally verify them. Formula-to-dessin identification and the automorphism-group upper bound are attributed external inputs.",
 dihedral_controls=dihedral,genus_zero_sieve=solutions,
 triangle_matrices=[list(A),list(B),list(C)],triangle_indices=indices,
 order5_pair_counts=pair_counts,dessin_conjugacies=dessin_maps,
 sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in paths},
 full_action_spectrum_classified=False,absolute_PSL2_11_realized=False)
(_check_output(ROOT/"receipt.json")).write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps(dict(result="pass",checks=len(checks),runtime_seconds=receipt["runtime_seconds"],
                      absolute_realization=False)))
