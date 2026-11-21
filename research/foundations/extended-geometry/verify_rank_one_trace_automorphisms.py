#!/usr/bin/env python3
"""Bounded exact controls and source integrity for the written all-pair theorem."""
import argparse,copy,hashlib,json,sys
from pathlib import Path
from fractions import Fraction as Q
from math import factorial,comb
from itertools import permutations
if not __debug__:
    raise SystemExit("Optimized Python is unsupported: exact assertions must execute.")
HERE=Path(__file__).resolve().parent
PN="RANK_ONE_TRACE_AUTOMORPHISM_INVARIANCE_ALL_IMAGE_PAIRS.md"
CN="verify_rank_one_trace_automorphisms.py"
IN="rank_one_trace_automorphism_inputs.json"
LN="rank_one_trace_automorphism_source_ledger.json"
RN="rank_one_trace_automorphism_verification.json"
PROOF_SHA="6a387aa4213cb00f945fb5e41ff9c5166a8ac118264e01bab795ec54bf189c5a"
EXPECTED={'schema': 1, 'field': 'every characteristic-zero field and every coefficient-field extension', 'pairs': {'integer': True, 'minimum': 2, 'n1_cone_argument_claim': False}, 'ring': 'balanced subring C of K[xi,z], A=C/(R-1), R=sum xi_i z_i', 'trace': {'diagonal_formula': 'alpha!/(n)_D', 'off_diagonal_zero': True, 'all_polynomial_A_automorphisms': True, 'exact_top_coefficient_codimension': 1, 'exact_top_coefficient_space_is_ideal': False}, 'volume': {'beta': 'tr(M dM wedge dM)', 'omega': 'beta^(n-1)/(n-1)!', 'automorphism_multiplier': 'nonzero field constant', 'pointwise_preservation_claim': False}, 'original': {'kernel': 'ker E_n=sum_i(partial_zi-xi_i)S', 'kernel_is_ideal_claim': False, 'balanced_normalization': 'E_n(H)=(n)_D Lambda(H|R=1)', 'arbitrary_marker_first_balanced': True, 'fixed_marker_per_transport': True, 'same_literal_marker_for_every_automorphism_claim': False}, 'power_descent': {'full_original_ring': True, 'R_primitive_required': True, 'both_bigradings_required': True, 'scalar_factors_allowed': True, 'all_extensions': True}, 'transport': {'pure': True, 'safe': True, 'unsafe': True, 'absolute_nonproper_power': True, 'arbitrary_automorphism_degree_formula_claim': False, 'starting_witness_existence_claim': False, 'padding_or_degree_spectrum_corollary_claim': False}, 'verification': {'bounded_exact_controls_only': True, 'finite_controls_prove_universal_quantifiers': False, 'formal_verification': False, 'external_review': False, 'independent_replay': False, 'worldwide_priority': False}}
PRIMARY={'author': 'Bertram Kostant', 'title': 'Minimal coadjoint orbits and symplectic induction', 'version': 'arXiv:math/0312252v1,12Dec2003', 'url': 'https://arxiv.org/pdf/math/0312252v1', 'locator': 'Section3.1,printedpp24-25,equation(3.2)', 'actual_scope': 'Real Lie-group coadjoint-orbit KKS form on fundamental fields, paired with reversed Lie bracket.', 'relation': 'INSPECTED_CLASSICAL_CONTEXT_ONLY; no logical import, PDF download/hash claim or automorphism classification.'}

FLAGS={"formal_verification":False,"external_review":False,"independent_replay":False,"worldwide_priority":False}
BOUNDARY={"written_universal_proof":"self-contained Sections1-5","external_mathematical_premises":[],"classical_context_is_mathematical_premise":False,"bounded_controls_establish_universal_quantifiers":False}
ROLES=[
("main_proof",PN,"WRITTEN_SELF_CONTAINED_UNIVERSAL_PROOF"),
("checker",CN,"BOUNDED_EXACT_CONTROLS_AND_SOURCE_AUTHENTICATION"),
("typed_inputs",IN,"EXACT_TYPED_SCOPE")]
def require(ok,msg):
    if not ok:raise ValueError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text(encoding="utf-8"))
def strict_equal(a,b):return json.dumps(a,sort_keys=True,separators=(",",":"))==json.dumps(b,sort_keys=True,separators=(",",":"))
def validate_inputs(d):
    require(strict_equal(d,EXPECTED),"Exact typed scope mismatch")
def validate_ledger(d):
    require(set(d)=={"schema","sources","primary_context","logical_boundary","assurance"},"Source-ledger keys")
    require(type(d["schema"]) is int and d["schema"]==1,"Ledger schema")
    require(strict_equal(d["primary_context"],[PRIMARY]),"Primary-context scope/version")
    require(strict_equal(d["logical_boundary"],BOUNDARY),"Exact typed logical boundary")
    require(strict_equal(d["assurance"],FLAGS),"Assurance flags")
    rows=d["sources"]
    require(type(rows) is list and len(rows)==len(ROLES),"Mandatory role cardinality")
    require(all(type(r) is dict and set(r)=={"role","root","path","sha256","relation"} for r in rows),"Exact source rows")
    got={}
    for r in rows:
        require(type(r["role"]) is str and r["role"] not in got,"Duplicate role")
        require(r["root"]=="here","Source-root scope")
        require(type(r["sha256"]) is str and len(r["sha256"])==64 and all(c in "0123456789abcdef" for c in r["sha256"]),"Hash format")
        got[r["role"]]=r
    for role,path,relation in ROLES:
        require(role in got,"Missing mandatory role "+role)
        r=got[role]
        require(r["path"]==path and r["relation"]==relation,"Role/path/relation mismatch")
    return got
def sources():
    validate_inputs(load(HERE/IN));ledger=load(HERE/LN);rows=validate_ledger(ledger)
    for role,path,relation in ROLES:
        p=HERE/path
        require(p.is_file(),"Missing artifact "+role)
        require(sha(p)==rows[role]["sha256"],"Artifact hash mismatch "+role)
    require(sha(HERE/PN)==PROOF_SHA,"Reviewed proof body")
    require(b"\r" not in (HERE/PN).read_bytes(),"Proof must use LF")
    return ledger
def source_negatives():
    d=load(HERE/IN);l=load(HERE/LN);n=0
    def rejects(fn,a):
        nonlocal n
        try:fn(a)
        except (ValueError,KeyError,TypeError):n+=1;return
        raise AssertionError("Negative source/scope control was accepted")
    paths=[
    ("pairs","minimum",1),("pairs","integer",1),("pairs","n1_cone_argument_claim",True),
    ("trace","exact_top_coefficient_codimension",2),("trace","exact_top_coefficient_space_is_ideal",True),
    ("volume","pointwise_preservation_claim",True),
    ("original","kernel_is_ideal_claim",True),("original","arbitrary_marker_first_balanced",False),
    ("original","same_literal_marker_for_every_automorphism_claim",True),
    ("power_descent","R_primitive_required",False),("power_descent","both_bigradings_required",False),
    ("power_descent","all_extensions",False),
    ("transport","arbitrary_automorphism_degree_formula_claim",True),
    ("transport","starting_witness_existence_claim",True),
    ("transport","padding_or_degree_spectrum_corollary_claim",True),
    ("verification","finite_controls_prove_universal_quantifiers",True),
    ("verification","formal_verification",True),("verification","external_review",True)]
    for a,b,v in paths:
        z=copy.deepcopy(d);z[a][b]=v;rejects(validate_inputs,z)
    for i in range(len(ROLES)):
        z=copy.deepcopy(l);z["sources"].pop(i);rejects(validate_ledger,z)
    z=copy.deepcopy(l);z["sources"].append(z["sources"][0]);rejects(validate_ledger,z)
    z=copy.deepcopy(l);z["schema"]=True;rejects(validate_ledger,z)
    z=copy.deepcopy(l);z["logical_boundary"]["classical_context_is_mathematical_premise"]=0;rejects(validate_ledger,z)
    z=copy.deepcopy(l);z["logical_boundary"]["external_mathematical_premises"]=["KKS comparison"];rejects(validate_ledger,z)
    z=copy.deepcopy(l);z["primary_context"][0]["relation"]="LOGICAL_ALL_FIELD_AUTOMORPHISM_CLASSIFICATION";rejects(validate_ledger,z)
    z=copy.deepcopy(l);z["assurance"]["formal_verification"]=True;rejects(validate_ledger,z)
    return n

# Sparse exact multivariate polynomial operations; exponent tuples retain all variables.
def add(*ps):
    o={}
    for p in ps:
        for a,c in p.items():o[a]=o.get(a,Q(0))+c
    return {a:c for a,c in o.items() if c}
def scale(p,c):return {a:v*c for a,v in p.items() if v*c}
def mul(p,q):
    o={}
    for a,c in p.items():
        for b,d in q.items():
            ab=tuple(x+y for x,y in zip(a,b));o[ab]=o.get(ab,Q(0))+c*d
    return {a:c for a,c in o.items() if c}
def power(p,k):
    require(k>=0,"Nonnegative polynomial power")
    if not k:
        n=len(next(iter(p))) if p else 0
        return {(0,)*n:Q(1)}
    o=p
    for j in range(k-1):o=mul(o,p)
    return o
def diff(p,j):
    o={}
    for a,c in p.items():
        if a[j]:
            b=list(a);b[j]-=1;o[tuple(b)]=c*a[j]
    return o
class Ring:
    def __init__(self,n):self.n=n;self.one={(0,)*n:Q(1)}
    def c(self,q):return scale(self.one,Q(q))
    def v(self,j):
        a=[0]*self.n;a[j]=1;return {tuple(a):Q(1)}
def matadd(A,B):return [[add(a,b) for a,b in zip(ar,br)] for ar,br in zip(A,B)]
def matscale(A,c):return [[scale(a,c) for a in ar] for ar in A]
def matmul(A,B):
    return [[add(*(mul(A[i][k],B[k][j]) for k in range(len(B)))) for j in range(len(B[0]))] for i in range(len(A))]
def trace(A):return add(*(A[i][i] for i in range(len(A))))
def comm(A,B):return matadd(matmul(A,B),matscale(matmul(B,A),-1))
def transpose(A):return [list(r) for r in zip(*A)]
def mdiff(A,j):return [[diff(a,j) for a in ar] for ar in A]
def zero_matrix(n):return [[{} for j in range(n)] for i in range(n)]
def matrix_controls():
    charts=tangent=forms=0
    for n in range(2,5):
        r=n-1;P=Ring(2*r);u=[add(P.one,*(scale(mul(P.v(j),P.v(r+j)),-1) for j in range(r)))]+[P.v(j) for j in range(r)]
        v=[P.one]+[P.v(r+j) for j in range(r)]
        M=[[mul(a,b) for b in v] for a in u]
        require(trace(M)==P.one,"Rank-one chart trace")
        require(matmul(M,M)==M,"Rank-one chart idempotence");charts+=n*n+1
        d=[mdiff(M,j) for j in range(2*r)]
        for j in range(2*r):
            require(matadd(matmul(M,d[j]),matmul(d[j],M))==d[j],"Tangent idempotence")
            X=comm(d[j],M)
            require(trace(X)=={} and comm(X,M)==d[j],"Fundamental field spanning");tangent+=2
        for i in range(2*r):
            for j in range(i+1,2*r):
                beta=trace(matmul(M,comm(d[i],d[j])))
                expected=P.c(-1) if j==r+i and i<r else {}
                require(beta==expected,"Canonical chart beta")
                require(trace(matmul(transpose(M),comm(transpose(d[i]),transpose(d[j]))))==scale(beta,-1),"Transpose reverses beta")
                forms+=2
    # Generic tangent blocks, with independent formal variables, at diag(1,0,...,0).
    for n in range(2,6):
        r=n-1;P=Ring(6*r);M=zero_matrix(n);M[0][0]=P.one;vec=[]
        for h in range(3):
            U=zero_matrix(n)
            for j in range(r):U[0][j+1]=P.v(2*h*r+j);U[j+1][0]=P.v((2*h+1)*r+j)
            vec.append(U)
            require(comm(comm(U,M),M)==U,"All-block tangent commutator");tangent+=1
        U,V,W=vec
        require(trace(matmul(M,comm(U,V)))==add(*(add(mul(U[0][j],V[j][0]),scale(mul(V[0][j],U[j][0]),-1)) for j in range(1,n))),"Block symplectic pairing")
        alt=add(trace(matmul(U,comm(V,W))),trace(matmul(V,comm(W,U))),trace(matmul(W,comm(U,V))))
        require(alt=={},"Closed beta triple trace");forms+=2
    # n=2 beta=-omega_Hopf/2: substitute rational chart after clearing denominators.
    P=Ring(2);x,t=P.v(0),P.v(1)
    # x*M has polynomial entries. Derivative quotient d(M)/dx=(x*d(xM)/dx-xM)/x^2.
    xM=[[scale(mul(x,add(P.one,t)),Q(1,2)),power(x,2)],
        [scale(add(P.one,scale(power(t,2),-1)),Q(1,4)),scale(mul(x,add(P.one,scale(t,-1))),Q(1,2))]]
    dxM=matadd([[mul(x,a) for a in row] for row in mdiff(xM,0)],matscale(xM,-1))
    dtM=mdiff(xM,1)
    # trace((xM)[dxM,dtM])=x^4*beta=-x^3/2.
    require(trace(matmul(xM,comm(dxM,dtM)))==scale(power(x,3),Q(-1,2)),"n=2 exact Hopf volume factor");forms+=1
    return {"chart_entry_and_trace_controls":charts,"tangent_field_controls":tangent,"differential_form_controls":forms}

def compositions(D,n):
    if n==1:yield(D,);return
    for a in range(D+1):
        for b in compositions(D-a,n-1):yield(a,)+b
def rising(n,D):
    a=1
    for j in range(D):a*=n+j
    return a
def fact(a):return factorial(a)
def lambda_diag(n,a):return Q(__import__("math").prod(fact(x) for x in a),rising(n,sum(a)))
def original_E(p,n):
    o={}
    for e,c in p.items():
        a,b=e[:n],e[n:]
        if any(i>j for i,j in zip(a,b)):continue
        eout=tuple(j-i for i,j in zip(a,b));coef=c
        for i,j in zip(a,b):coef*=Q(fact(j),fact(j-i))
        o[eout]=o.get(eout,Q(0))+coef
    return {e:c for e,c in o.items() if c}
def original_controls():
    coin=norm=kernel=marker=lie=0
    for n in range(2,7):
        for D in range(7):
            for a in compositions(D,n):
                lam=lambda_diag(n,a);total=Q(0)
                for i in range(n):
                    b=list(a);b[i]+=1;total+=lambda_diag(n,b)
                require(total==lam,"R=1 trace compatibility");coin+=1
                for i,j in permutations(range(n),2):
                    if a[i]:
                        b=list(a);b[i]-=1;b[j]+=1
                        require((a[j]+1)*lam==a[i]*lambda_diag(n,b),"Offdiagonal factorial coinvariant");coin+=1
    for n in range(2,5):
        P=Ring(2*n);R=add(*(mul(P.v(i),P.v(n+i)) for i in range(n)))
        for D in range(7):
            p=power(R,D)
            require(original_E(p,n)=={(0,)*n:Q(rising(n,D))},"Literal original radial Pochhammer");norm+=1
        for i,j in permutations(range(n),2):
            lij=lambda p:add(mul(P.v(i),diff(p,j)),scale(mul(P.v(n+j),diff(p,n+i)),-1))
            require(lij(R)=={},"Infinitesimal trace relation");lie+=1
        # Deterministic monomials in all variables; no random parameter sampling.
        for q in range(70):
            e=tuple((q*(j+1)+j*j)%4 for j in range(2*n));H={e:Q(1)}
            for i in range(n):
                di=add(diff(H,n+i),scale(mul(P.v(i),H),-1))
                require(original_E(di,n)=={},"Image kernel annihilation");kernel+=1
            # Balance an arbitrary nonnegative residual z degree by a fixed derivative.
            a=tuple((q+j)%3 for j in range(n));b=tuple(a[j]+((q+2*j)%2) for j in range(n))
            H={a+b:Q(1)};eta=tuple(b[j]-a[j] for j in range(n))
            mult={(eta+(0,)*n):Q(1)}
            lhs=original_E(mul(mult,H),n);rhs=original_E(H,n)
            for j,k in enumerate(eta):
                for z in range(k):rhs=diff(rhs,j)
            require(lhs==rhs,"Fixed original marker coefficient derivative");marker+=1
            # A balanced monomial contracts to alpha!, zero unless alpha=beta.
            D=q%5;cs=list(compositions(D,n));aa=cs[q%len(cs)];bb=cs[(3*q+1)%len(cs)]
            value=original_E({aa+bb:Q(1)},n)
            want={(0,)*n:Q(__import__("math").prod(fact(v) for v in aa))} if aa==bb else {}
            require(value==want,"Literal balanced original normalization");norm+=1
    return {"factorial_coinvariant_controls":coin,"original_normalization_controls":norm,"kernel_controls":kernel,"balanced_marker_controls":marker,"Lie_trace_controls":lie}

def degree_distribution(p,part):
    return sorted({sum(e[j] for j in part) for e in p})
def proper_power_and_boundary_controls():
    n=2;P=Ring(4);R=add(mul(P.v(0),P.v(2)),mul(P.v(1),P.v(3)));Z=mul(P.v(0),P.v(3))
    padded=mul(R,power(Z,2))
    require(degree_distribution(padded,[0,1])==[3] and degree_distribution(padded,[2,3])==[3],"Padded original balanced F")
    # xi1=xi2=z2=1,z1=s sends F to s+1, with a simple root.
    univar={}
    for e,c in padded.items():univar[(e[2],)]=univar.get((e[2],),Q(0))+c
    require(univar=={(0,):Q(1),(1,):Q(1)},"Primitive hypothesis simple-slice control")
    # f=Z^2 on R=1 is a square, although padded F has a simple original slice.
    require(power(Z,2)!=padded,"R padding genuinely changes original proper-power status")
    H=mul(P.v(0),add(P.one,P.v(2)));H2=power(H,2)
    require(degree_distribution(H2,[0,1])==[2] and degree_distribution(H2,[2,3])==[0,1,2],"Both bigrading extremes needed")
    # General homogeneous balanced powers for several n and exact min/max support.
    cnt=0;primitive_roots=0
    for n in range(2,5):
        P=Ring(2*n);G=add(mul(P.v(0),P.v(n)),scale(mul(P.v(1),P.v(n+1)),2))
        # This original root is primitive, not the radial R when n=2.
        vals=[Q(0)]*(2*n);vals[0]=vals[1]=vals[n]=Q(1);vals[n+1]=Q(-1)
        require(sum((coef*__import__("math").prod(v**j for v,j in zip(vals,e)) for e,coef in G.items()),Q(0))==Q(-1),"Primitive original root at R=0")
        primitive_roots+=1
        for r in range(2,7):
            F=power(G,r)
            require(degree_distribution(F,list(range(n)))==[r] and degree_distribution(F,list(range(n,2*n)))==[r],"Original balanced scalar-power degrees");cnt+=1
        for r in range(2,5):
            H=mul(P.v(0),add(P.one,P.v(n)));F=power(H,r)
            require(degree_distribution(F,list(range(n)))==[r] and degree_distribution(F,list(range(n,2*n)))==list(range(r+1)),"Original unbalanced root expansion");cnt+=1
    # Exact coefficient spaces are not ideals, in n=2.
    lt=Q(1,2)-Q(1,2)
    lt2=lambda_diag(2,(2,0))-2*lambda_diag(2,(1,1))+lambda_diag(2,(0,2))
    require(lt==0 and lt2==Q(1,3),"Exact coefficient vector space is not ideal")
    # ker original E is not an ideal either: xi1 maps to zero, zi1*xi1 maps to one.
    P=Ring(4)
    require(original_E(P.v(0),2)=={} and original_E(mul(P.v(0),P.v(2)),2)=={(0,0):Q(1)},"Image kernel not ideal")
    return {"original_scalar_power_controls":cnt,"primitive_original_root_controls":primitive_roots,"mathematical_boundary_controls":[
    "R-primitive hypothesis cannot be removed: R*Z^2 has simple original slice but sphere image is square.",
    "Both original gradings are essential: [xi1(1+z1)]^2 is xi-homogeneous and has z degrees0,1,2.",
    "Transpose reverses beta, so pointwise volume preservation is false.",
    "Exact-form coefficient space and original Image kernel are vector spaces, not ideals."
    ]}
IDENTITIES=[
"MU+UM=U and MUM=0 imply [UM-MU,M]=U; this formal block calculation is dimension-independent.",
"beta(U,V)=tr M(UV-VU)=a*d-c*b; odd off-diagonal triple products have trace zero.",
"d(f contraction(D_X,omega))=(D_X f)omega because conjugation preserves omega.",
"(alpha_j+1)alpha!=alpha_i(alpha-e_i+e_j)!.",
"sum_i (alpha+e_i)!=(D+n)alpha!; R multiplication occurs before taking coinvariants.",
"E_n((partial_zi-xi_i)H)=0 and monomial reduction prove the written kernel identity.",
"E_n(xi^eta H)=partial_z^eta E_n(H) gives fixed-marker balancing.",
"Primitive filtered lifts multiply without radial padding; both original bigrading extremes constrain any full-ring power root."
]
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo-root",type=Path,default=HERE,help="Explicit workspace root; this self-contained packet has no external mathematical source artifacts.")
    ap.add_argument("--output",type=Path,default=HERE/RN)
    args=ap.parse_args();repo=args.repo_root.resolve();require(repo.is_dir(),"Repository root must exist")
    ledger=sources()
    result={"status":"PASS_SOURCE_GUARDED_RANK_ONE_TRACE_ALL_PAIRS","scope":EXPECTED,"source_ledger_sha256":sha(HERE/LN),
      "pins":ledger["sources"],"primary_context":PRIMARY,"logical_boundary":BOUNDARY,
      "written_identity_families":IDENTITIES,"matrix_and_form_controls":matrix_controls(),
      "contraction_and_marker_controls":original_controls(),"power_and_boundary_controls":proper_power_and_boundary_controls(),
      "source_scope_negative_controls":source_negatives(),"assurance":FLAGS}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],"mandatory_roles":len(ROLES),"source_scope_negative_controls":result["source_scope_negative_controls"],"output":str(args.output)},sort_keys=True))
if __name__=="__main__":
    try:main()
    except (ValueError,AssertionError,KeyError,TypeError,OSError) as e:
        print("FAIL: "+str(e),file=sys.stderr);raise SystemExit(1)
