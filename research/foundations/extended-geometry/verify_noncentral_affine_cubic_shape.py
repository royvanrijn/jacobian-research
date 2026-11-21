#!/usr/bin/env python3
"""Bounded exact controls and source integrity for constant affine cubic classification."""
import argparse,copy,hashlib,json,sys
from pathlib import Path
from fractions import Fraction as Q
from math import comb,factorial
from itertools import permutations
if not __debug__:
    raise SystemExit("Optimized Python is unsupported: exact assertions must execute.")
HERE=Path(__file__).resolve().parent
PN="NONCENTRAL_CONSTANT_AFFINE_CUBIC_SHAPE_CLASSIFICATION.md"
CN="verify_noncentral_affine_cubic_shape.py"
IN="noncentral_affine_cubic_shape_inputs.json"
LN="noncentral_affine_cubic_shape_source_ledger.json"
RN="noncentral_affine_cubic_shape_verification.json"
PROOF_SHA="20ff19e08cc7decb8564346c363e2a3657d846241267a462899934d481458e07"
EXPECTED={'schema': 1, 'field': 'every characteristic-zero field and every coefficient-field extension', 'shape': {'lambda_constant': True, 'excluded_lambda': '1/2', 'E_polynomial_unit_at_zero': True, 's0': -1, 'h0': 2, 'formal_s_h_allowed': True, 'requires_all_kappa_m_nonzero': False, 'purity_forces_constant_shape_claim': False}, 'purity': {'criterion': "pure iff h/2=1-x(E*h^3)'/(E*h^3)", 'universal_theorem_on_fixed_valid_H2_data': True, 'unconditional_Pure_f_seed': False}, 'scalar': {'recurrence_minimum_m': 1, 'm0_recurrence_claim': False, 'local_continuant_first_squared_index': 'm+2', 'actual_kappa_has_Jacobi_initial_data_claim': False, 'fixed_gap_finite_for_r_Delta_nonzero': True, 'all_exceptions_handled': ['0', '1', '-1', '2'], 'only_midpoint_excluded': True}, 'support': {'algebraic_required': True, 'arbitrarily_late_zero_blocks_required': True, 'infinitely_many_zero_coefficients_alone_suffices': False, 'second_marked_support_application_required': True}, 'rationality': {'h': '2/Q', 'Q0': 1, 'Q_squarefree': True, 'nonconstant_Q_degree_E': '3*degQ+1', 'u_prime_and_J_may_vanish_at_finite_fiber': True}, 'safety': {'safe_iff_Q_constant': True, 'safe_degrees': '2,5,8,...', 'unsafe_degrees': '6,9,12,...', 'every_nonconstant_Kx_marker_infinitely_bad': True, 'every_order_marked_nonzero_claim': False, 'all_arbitrary_A_markers_bad_claim': False}, 'nonclaims': {'generic_mixed_cubic_classification': False, 'global_degree5_resolved': False, 'global_degree7_resolved': False, 'general_original_absolute_nonpower': False, 'new_original_Image_witness_theorem': False}, 'source_boundary': {'self_contained_math': True, 'Stanley': 'inspected context only', 'Wimp': 'inspected local-associated-Jacobi context only', 'FPYZ_logical_premise': False, 'all_Aut_logical_premise': False, 'old429_logical_premise': False}, 'verification': {'bounded_exact_controls_only': True, 'finite_controls_prove_universal_quantifiers': False, 'formal_verification': False, 'external_review': False, 'independent_replay': False, 'worldwide_priority': False}}
PRIMARY=[{'author': 'Richard P. Stanley', 'title': 'Differentiably finite power series', 'year': 1980, 'journal': 'European Journal of Combinatorics1,175-188', 'url': 'https://math.mit.edu/~rstan/pubs/pubfiles/45.pdf', 'doi': '10.1016/S0195-6698(80)80051-5', 'locators': ['Theorem2.1 printedp178', 'Definition1.3 and Theorem1.5 printedp176'], 'relation': 'INSPECTED_CLASSICAL_CONTEXT_ONLY; printed statements over complex coefficients; direct all-field proof supplied.'}, {'author': 'Jet Wimp', 'title': 'Explicit formulas for the associated Jacobi polynomials and some applications', 'year': 1987, 'journal': 'Canadian Journal of Mathematics39(4),983-1000', 'url': 'https://doi.org/10.4153/CJM-1987-050-4', 'doi': '10.4153/CJM-1987-050-4', 'locator': 'equation12 printedp986; definitions printedpp983-985', 'parameters': {'alpha': '1/3', 'beta': '-1/3', 'gamma': '1', 'association': 'c=m+1', 'first_squared_index': 'm+2'}, 'identity': 'D_(m,h)(z)=27^(-h/2)*P_hat_h^(1/3,-1/3)(sqrt(27)*z;c=m+1)', 'relation': 'INSPECTED_CLASSICAL_CONTEXT_ONLY; local monic continuants, not scalar initial data or an imported association-shift zero theorem.'}]

FLAGS={"formal_verification":False,"external_review":False,"independent_replay":False,"worldwide_priority":False}
BOUNDARY={"external_mathematical_premises":[],"universal_mathematics":"self-contained written H1-H19 and all-field proofs","classical_sources_are_logical_premises":False,"actual_phase_properties_unconditionally_seeded":False,"bounded_controls_prove_universal_quantifiers":False}
ROLES=[("main_proof",PN,"WRITTEN_SELF_CONTAINED_UNIVERSAL_PROOF"),("checker",CN,"BOUNDED_EXACT_IDENTITIES_AND_SOURCE_AUTHENTICATION"),("typed_inputs",IN,"EXACT_TYPED_SCOPE")]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text(encoding="utf-8"))
def strict_equal(a,b):return json.dumps(a,sort_keys=True,separators=(",",":"))==json.dumps(b,sort_keys=True,separators=(",",":"))
def require(ok,msg):
    if not ok:raise ValueError(msg)
def validate_inputs(d):require(strict_equal(d,EXPECTED),"Exact typed theorem/scope mismatch")
def validate_ledger(d):
    require(set(d)=={"schema","sources","classical_context","logical_boundary","assurance"},"Ledger keys")
    require(type(d["schema"]) is int and d["schema"]==1,"Ledger schema")
    require(strict_equal(d["classical_context"],PRIMARY),"Exact classical context/parameter boundary")
    require(strict_equal(d["logical_boundary"],BOUNDARY),"Exact typed logical boundary")
    require(strict_equal(d["assurance"],FLAGS),"Assurance flags")
    rows=d["sources"];require(type(rows) is list and len(rows)==len(ROLES),"Mandatory role cardinality")
    got={}
    for r in rows:
        require(type(r) is dict and set(r)=={"role","root","path","sha256","relation"},"Exact source-row schema")
        require(type(r["role"]) is str and r["role"] not in got,"Duplicate source role")
        require(type(r["sha256"]) is str and len(r["sha256"])==64 and all(c in "0123456789abcdef" for c in r["sha256"]),"Hash format")
        got[r["role"]]=r
    for role,path,relation in ROLES:
        require(role in got,"Missing role "+role);r=got[role]
        require((r["root"],r["path"],r["relation"])==("here",path,relation),"Role/path/root/relation mismatch")
    return got
def check_sources():
    validate_inputs(load(HERE/IN));d=load(HERE/LN);rows=validate_ledger(d)
    for role,path,relation in ROLES:
        p=HERE/path;require(p.is_file(),"Missing artifact "+role);require(sha(p)==rows[role]["sha256"],"Artifact hash mismatch "+role)
    require(sha(HERE/PN)==PROOF_SHA,"Reviewed proof body")
    require(b"\r" not in (HERE/PN).read_bytes(),"Proof requires LF")
    return d
def source_negatives():
    n=0;di=load(HERE/IN);le=load(HERE/LN)
    def reject(fn,x):
        nonlocal n
        try:fn(x)
        except (ValueError,KeyError,TypeError):n+=1;return
        raise AssertionError("Negative source/scope control accepted")
    cases=[("shape","lambda_constant",False),("shape","excluded_lambda","none"),("shape","E_polynomial_unit_at_zero",False),("shape","s0",1),("shape","h0",1),("shape","requires_all_kappa_m_nonzero",True),("shape","purity_forces_constant_shape_claim",True),
    ("purity","unconditional_Pure_f_seed",True),("scalar","recurrence_minimum_m",0),("scalar","m0_recurrence_claim",True),("scalar","local_continuant_first_squared_index","m+1"),("scalar","actual_kappa_has_Jacobi_initial_data_claim",True),("scalar","only_midpoint_excluded",False),
    ("support","algebraic_required",False),("support","arbitrarily_late_zero_blocks_required",False),("support","infinitely_many_zero_coefficients_alone_suffices",True),("support","second_marked_support_application_required",False),
    ("rationality","Q0",True),("rationality","u_prime_and_J_may_vanish_at_finite_fiber",False),
    ("safety","every_order_marked_nonzero_claim",True),("safety","all_arbitrary_A_markers_bad_claim",True),("safety","safe_degrees","all"),("safety","unsafe_degrees","5,7"),
    ("nonclaims","generic_mixed_cubic_classification",True),("nonclaims","global_degree5_resolved",True),("nonclaims","global_degree7_resolved",True),("nonclaims","general_original_absolute_nonpower",True),("nonclaims","new_original_Image_witness_theorem",True),
    ("source_boundary","all_Aut_logical_premise",True),("verification","finite_controls_prove_universal_quantifiers",True),("verification","external_review",True)]
    for a,b,v in cases:
        x=copy.deepcopy(di);x[a][b]=v;reject(validate_inputs,x)
    for i in range(len(ROLES)):
        x=copy.deepcopy(le);x["sources"].pop(i);reject(validate_ledger,x)
    x=copy.deepcopy(le);x["sources"].append(x["sources"][0]);reject(validate_ledger,x)
    x=copy.deepcopy(le);x["schema"]=True;reject(validate_ledger,x)
    x=copy.deepcopy(le);x["logical_boundary"]["classical_sources_are_logical_premises"]=0;reject(validate_ledger,x)
    x=copy.deepcopy(le);x["logical_boundary"]["external_mathematical_premises"]=["spectral theorem"];reject(validate_ledger,x)
    x=copy.deepcopy(le);x["classical_context"][1]["parameters"]["association"]="c=m";reject(validate_ledger,x)
    x=copy.deepcopy(le);x["classical_context"][1]["relation"]="LOGICAL_FIXED_SHIFT_ZERO_THEOREM";reject(validate_ledger,x)
    x=copy.deepcopy(le);x["classical_context"][0]["relation"]="LOGICAL_ALL_FIELD_IMPORT";reject(validate_ledger,x)
    x=copy.deepcopy(le);x["assurance"]["formal_verification"]=True;reject(validate_ledger,x)
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


def coeff(p,j,k):
    o={}
    for e,c in p.items():
        if e[j]==k:
            z=list(e);z[j]=0;o[tuple(z)]=o.get(tuple(z),Q(0))+c
    return {e:c for e,c in o.items() if c}
def at(p,j,a):
    o={}
    for e,c in p.items():
        z=list(e);z[j]=0;o[tuple(z)]=o.get(tuple(z),Q(0))+c*Q(a)**e[j]
    return {e:c for e,c in o.items() if c}
def degree(p,j):return max((e[j] for e in p),default=-10**6)
def integrate_interval(p,j,lo,hi):
    o={}
    for e,c in p.items():
        z=list(e);z[j]=0;v=c*(Q(hi)**(e[j]+1)-Q(lo)**(e[j]+1))/Q(e[j]+1)
        o[tuple(z)]=o.get(tuple(z),Q(0))+v
    return {e:c for e,c in o.items() if c}
def falling(p,j,R):
    o=R.one
    for k in range(j):o=mul(o,add(p,R.c(-k)))
    return o
def reduce_univar(p,g):
    p=dict(p);dg=degree(g,0);cg=g[(dg,)]
    while p and degree(p,0)>=dg:
        d=degree(p,0);c=p[(d,)]/cg;shift={(d-dg,):c};p=add(p,scale(mul(shift,g),-1))
    return p

def universal_algebra():
    # Universal IBP identity H6 in independent variables t,p,r,m.
    R=Ring(4);t,p,r,m=[R.v(j) for j in range(4)]
    phi=add(power(t,3),mul(p,t),r);Delta=add(scale(power(p,3),-4),scale(power(r,2),-27))
    sm=add(mul(add(scale(m,9),R.c(6)),power(t,4)),
           mul(mul(add(scale(m,15),R.c(12)),p),power(t,2)),scale(mul(r,t),-3),
           scale(mul(add(m,R.one),power(p,2)),4))
    left=add(mul(diff(sm,0),phi),mul(mul(m,sm),diff(phi,0)))
    right=add(scale(mul(mul(add(scale(m,3),R.c(2)),add(scale(m,3),R.c(4))),power(phi,2)),3),
              scale(mul(mul(mul(add(m,R.one),add(scale(m,2),R.one)),r),phi),-27),
              scale(mul(mul(m,add(m,R.one)),Delta),-1))
    require(left==right,"H6 universal IBP certificate")
    # Depressed coefficients H14, in t,s,h,lambda.
    R=Ring(4);t,s,h,l=[R.v(j) for j in range(4)];a=add(t,scale(s,-1))
    phi=add(power(a,3),scale(mul(mul(add(R.one,l),h),power(a,2)),-1),mul(mul(l,power(h,2)),a))
    C,B,D=[coeff(phi,0,j) for j in range(3)]
    pl=scale(add(power(l,2),scale(l,-1),R.one),Q(-1,3))
    rl=scale(mul(mul(add(l,R.one),add(scale(l,2),R.c(-1))),add(l,R.c(-2))),Q(-1,27))
    require(add(scale(B,3),scale(power(D,2),-1))==scale(mul(pl,power(h,2)),3),"Depressed h^2 identity")
    require(add(scale(C,27),scale(mul(B,D),-9),scale(power(D,3),2))==scale(mul(rl,power(h,3)),27),"Depressed h^3 identity")
    require(add(scale(power(pl,3),-4),scale(power(rl,2),-27))==mul(power(l,2),power(add(R.one,scale(l,-1)),2)),"Cubic discriminant identity")
    # H8 rational normalization, after clearing all denominators in m.
    R=Ring(1);m=R.v(0)
    na=scale(mul(add(m,R.one),add(scale(m,2),R.one)),9);da=mul(add(scale(m,3),R.c(2)),add(scale(m,3),R.c(4)))
    nap=scale(mul(m,add(scale(m,2),R.c(-1))),9);dap=mul(add(scale(m,3),R.c(-1)),add(scale(m,3),R.one))
    nb=mul(m,add(m,R.one));db=scale(da,3)
    nt=add(power(m,2),R.c(Q(-1,9)));dt=scale(add(scale(power(m,2),4),R.c(-1)),27)
    require(mul(mul(mul(nb,da),dap),dt)==mul(mul(mul(db,na),nap),nt),"Jacobi normalization t_m")
    # Wimp Eq12 specialization/monic coefficient; context only.
    An=add(scale(m,2),R.one);Ad=add(m,R.one);Apn=add(scale(m,2),R.c(-1));Apd=m
    Cn=add(power(m,2),R.c(Q(-1,9)));Cd=mul(m,add(m,R.one))
    require(mul(mul(mul(Cn,Ad),Apd),dt)==scale(mul(mul(mul(Cd,An),Apn),nt),27),"Wimp monic local coefficient identification")
    return ["H6 universal polynomial IBP certificate","H14 depressed h^2 coefficient","H14 depressed h^3 coefficient","Discriminant Delta=lambda^2(1-lambda)^2","H8 exact rational normalization","Wimp Eq12 monic local coefficient identification (context)"]

def scalar_moment_controls():
    R=Ring(2);v,l=R.v(0),R.v(1);phi=mul(mul(v,add(v,R.c(-1))),add(v,scale(l,-1)))
    K=[]
    for m in range(17):
        k=integrate_interval(power(phi,m),0,0,1)
        K.append({(e[1],):c for e,c in k.items()})
    U=Ring(1);l=U.v(0);pl=scale(add(power(l,2),scale(l,-1),U.one),Q(-1,3))
    rl=scale(mul(mul(add(l,U.one),add(scale(l,2),U.c(-1))),add(l,U.c(-2))),Q(-1,27))
    Delta=mul(power(l,2),power(add(U.one,scale(l,-1)),2))
    require(K[1]==scale(add(scale(l,2),U.c(-1)),Q(1,12)),"First moment")
    n=0
    for m in range(1,16):
        lhs=scale(K[m+1],3*(3*m+2)*(3*m+4))
        rhs=add(scale(mul(rl,K[m]),27*(m+1)*(2*m+1)),scale(mul(Delta,K[m-1]),m*(m+1)))
        require(lhs==rhs,"Actual scalar recurrence");n+=1
    for m in range(1,17):
        require(degree(K[m],0)==m and K[m][(m,)]==Q(factorial(m)**2,factorial(2*m+1)),"Nonzero kappa leading beta coefficient")
        for val in [-1,0,1,2]:require(at(K[m],0,val)!={},"Allowed exceptional shape moment")
        mid=at(K[m],0,Q(1,2));require(mid=={} if m%2 else list(mid.values())[0]>0,"Midpoint parity exception")
    # The m=0 recurrence fails by the exact nonzero polynomial below.
    fail0=add(scale(K[1],24),scale(rl,-27))
    require(fail0==mul(mul(l,add(l,U.c(-1))),add(scale(l,2),U.c(-1))),"m0 failure polynomial")
    require(at(fail0,0,2)==U.c(6),"m0 recurrence genuinely false")
    # A noncentral complex algebraic shape has kappa_2=0.
    g=add(scale(power(l,2),7),scale(l,-7),U.c(2))
    require(K[2]==scale(g,Q(1,210)),"Exceptional kappa2 polynomial")
    require(reduce_univar(K[2],g)=={} and reduce_univar(K[1],g)!={} and reduce_univar(K[3],g)!={},"Allowed kappa2-zero shape")
    return K,{"scalar_recurrence_controls":n,"kappa_degree_and_exception_orders":16,"m0_failure": "lambda(lambda-1)(2lambda-1), value6 atlambda2","allowed_kappa2_shape":"7lambda^2-7lambda+2=0, kappa1/kappa3 nonzero"}

def local_continuants():
    det=idx=0
    for h in range(1,7):
        R=Ring(h+1);z=R.v(0);D=[R.one,z]
        for j in range(2,h+1):D.append(add(mul(z,D[-1]),scale(mul(R.v(j),D[-2]),-1)))
        vals=[{},R.one]
        for j in range(1,h+1):vals.append(add(mul(z,vals[-1]),scale(mul(R.v(j),vals[-2]),-1)))
        require(vals[h+1]==D[h],"Local zero-start indexing")
        require(all(e[1]==0 for e in D[h]),"First used coefficient is m+2, not m+1");idx+=2
        M=[[{} for j in range(h)] for i in range(h)]
        for i in range(h):
            M[i][i]=z
            if i<h-1:M[i][i+1]=R.c(-1)
            if i>0:M[i][i-1]=scale(R.v(i+1),-1)
        out={}
        for perm in permutations(range(h)):
            term=R.one
            for i,j in enumerate(perm):term=mul(term,M[i][j])
            sign=(-1)**sum(perm[i]>perm[j] for i in range(h) for j in range(i+1,h))
            out=add(out,scale(term,sign))
        require(out==D[h],"Literal path determinant");det+=1
    # H10 constant continuant squarefree certificate.
    R=Ring(2);z,T=R.v(0),R.v(1);D=[R.one,z]
    for h in range(2,19):D.append(add(mul(z,D[-1]),scale(mul(T,D[-2]),-1)))
    for h in range(1,19):
        lhs=mul(add(power(z,2),scale(T,-4)),diff(D[h],0))
        rhs=add(scale(mul(z,D[h]),h),scale(mul(T,D[h-1]),-2*(h+1)))
        require(lhs==rhs,"H10 constant continuant derivative certificate")
    # H11 formal common m^-2 perturbation, with x=1/m.
    R=Ring(2);x,z=R.v(0),R.v(1);T=Q(1,108)
    def trunc(p):return {e:c for e,c in p.items() if e[0]<=3}
    C=[R.one,z];L=[R.one,z]
    for h in range(2,17):
        C.append(add(mul(z,C[-1]),scale(C[-2],-T)))
        entry=scale(add(R.one,scale(power(x,2),Q(5,36)),scale(power(x,3),Q(-5*h,18))),T)
        L.append(trunc(add(mul(z,L[-1]),scale(mul(entry,L[-2]),-1))))
    for h in range(1,17):
        require(coeff(L[h],0,0)==C[h] and coeff(L[h],0,1)=={},"H11 constant/linear terms")
        target=scale(add(scale(C[h],h),scale(mul(z,diff(C[h],1)),-1)),Q(5,72))
        require(coeff(L[h],0,2)==target,"H11 universal first nontrivial perturbation")
    return {"local_path_determinants":det,"local_index_controls":idx,"constant_continuant_H10_controls":18,"formal_H11_Laurent_controls":16}

def extract_t(p,j):return coeff(p,1,j)
def build_phase(P,branch,j):
    x,t,l=P.v(0),P.v(1),P.v(2);Qp=add(P.one,x)
    s=P.c(-1) if j==0 else add(P.c(-1),power(x,j))
    u=add(t,scale(s,-1))
    if branch=="safe":
        return add(power(u,3),scale(mul(add(P.one,l),power(u,2)),-2),scale(mul(l,u),4))
    return add(mul(power(Qp,4),power(u,3)),scale(mul(mul(add(P.one,l),power(Qp,3)),power(u,2)),-2),scale(mul(mul(l,power(Qp,2)),u),4))
def realization_controls(K):
    R=Ring(3);x,t,l=R.v(0),R.v(1),R.v(2);Qp=add(R.one,x);endpoint=marked=degrees=0
    for branch in ["safe","unsafe"]:
        for j in range(5):
            P=build_phase(R,branch,j);C,B,D,E=[extract_t(P,k) for k in range(4)]
            require(add(at(C,0,0),at(D,0,0))=={} and add(at(B,0,0),at(E,0,0))=={},"Literal cubic regularity")
            d=max(degree(C,0)-1,degree(B,0),degree(D,0)+1,degree(E,0)+2)
            wanted=(2 if j==0 else 3*j-1) if branch=="safe" else (6 if j==0 else 3*j+3)
            require(d==wanted,"Exact realization degree");degrees+=1
            for m in range(1,5):
                avg=scale(integrate_interval(power(P,m),1,-1,1),Q(1,2))
                kap={(0,0,e[0]):c for e,c in K[m].items()}
                ref=scale(kap,8**m)
                if branch=="unsafe":ref=mul(ref,power(Qp,m-1))
                for k in range(m+1):
                    require(coeff(avg,0,k)==coeff(ref,0,k),"Actual generic-lambda endpoint coefficients");endpoint+=1
                for ell in [1,2,3]:
                    require(coeff(avg,0,m-ell)==coeff(ref,0,m-ell),"Actual polynomial-x marked coefficient");marked+=1
    # Integer degree comparisons under the proved constraints, not assumed phase existence.
    comp=0
    for q in range(1,11):
        e=3*q+1
        for j in range(11):
            a=e+2 if j==0 else max(e+3*j-1,e+2*j,e+j+1,e+2)
            require(a==3*(q+max(1,j)),"All unsafe degree-comparison branches");comp+=1
    for j in range(11):
        require((2 if j==0 else max(3*j-1,2*j,j+1,2))==(2 if j==0 else 3*j-1),"Safe degree comparison");comp+=1
    return {"literal_realization_degree_controls":degrees,"generic_lambda_endpoint_coefficients":endpoint,"generic_lambda_marked_coefficients":marked,"integer_degree_comparison_controls":comp}

def support_and_omission_controls(K):
    # Fixed-shift falling factorials have distinct degrees and leading coefficient1.
    R=Ring(1);n=R.v(0);tri=0
    for shift in range(-4,5):
        for j in range(7):
            p=falling(add(n,R.c(shift)),j,R)
            require(degree(p,0)==j and p[(j,)]==1,"Recurrence shift leading coefficient");tri+=1
    # Actual noncentral kappa2-zero diagonal: coefficient support alone need not force B=1.
    U=Ring(1);l=U.v(0);g=add(scale(power(l,2),7),scale(l,-7),U.c(2))
    require(reduce_univar(K[2],g)=={},"Nonconstant polynomial B=1+u^2 supported at allowed scalar zero")
    # Ramification is allowed in the omitted-value argument: E has a pole nonetheless.
    # u=x+x^2, b=-1/4, J=(1+2x)/(1+x), B(u)=1+4u=(1+2x)^2.
    x=U.v(0);q=add(U.one,x);r=add(U.one,scale(x,2))
    require(add(U.one,scale(add(x,power(x,2)),4))==power(r,2),"Ramified finite fiber B composition")
    require(diff(add(x,power(x,2)),0)==r,"Ramified rational-map derivative")
    require(at(power(q,2),0,Q(-1,2))==U.c(Q(1,4)) and at(power(r,9),0,Q(-1,2))=={},"E=(1+x)^2/[8(1+2x)^9] has a genuine finite pole")
    # Essential midpoint exclusion: central shape gives Bdiag=1-u and pure but ell!=J.
    # a=1+x,h=2/(1+x)^2,E=(1+x)^7/8,s=-h/2.
    require(add(q,scale(x,-1))==U.one,"Central Bdiag=1-u identity")
    require(add(U.one,scale(q,-1))==scale(x,-1),"Central ell-J=-x/(1+x)^2 nonzero")
    require(K[1] and at(K[1],0,Q(1,2))=={},"Central B linear coefficient can survive purity")
    V=Ring(2);xx,tt=V.v(0),V.v(1);qq=add(V.one,xx)
    ee=scale(power(qq,7),Q(1,8))
    central=scale(add(mul(power(qq,7),power(tt,3)),scale(mul(power(qq,3),tt),-1)),Q(1,8))
    require(mul(power(qq,4),central)==mul(ee,add(mul(power(qq,4),power(tt,3)),scale(tt,-1))),"Literal central polynomial-compatible cubic shape")
    central_count=0
    for m in range(1,7):
        avg=scale(integrate_interval(power(central,m),1,-1,1),Q(1,2))
        require(coeff(avg,0,m)=={} and coeff(avg,0,m-1)=={},"Literal midpoint pure and x-marked controls")
        central_count+=2
    # Algebraic series alone + infinitely many zeros is insufficient:
    # F=1/(1-u^2) and F=u/(1-u^2) are rational nonpolynomials.
    even=[Q(1 if j%2==0 else 0) for j in range(24)]
    odd=[Q(1 if j%2 else 0) for j in range(24)]
    require(all(even[j+2]==even[j] and odd[j+2]==odd[j] for j in range(22)),"Exact rational parity recurrences")
    require(all(at(K[m],0,Q(1,2))=={} or odd[m]==0 for m in range(1,17)),"Nonpolynomial marked series can be killed by midpoint kappas")
    return {"shift_triangular_controls":tri,"central_literal_boundary_controls":central_count,"mathematical_negative_controls":[
        "Allowed noncentral shape7lambda^2-7lambda+2=0 has kappa2=0: all-kappa nonvanishing and EVERY-order marked survival are false.",
        "The scalar recurrence cannot be extended to m0; its defect is lambda(lambda-1)(2lambda-1).",
        "Algebraic rational1/(1-u^2) has infinitely many zero coefficients but is nonpolynomial: late arbitrary-length zero blocks are essential.",
        "Nonpolynomial algebraic u/(1-u^2) is killed coefficientwise by midpoint kappas: marked-tail support must be reused.",
        "Polynomial Bdiag=1+u^2 meets scalar diagonal constraints at the kappa2 shape: the omitted-value geometry is indispensable.",
        "Central lambda1/2 admits polynomial-compatible a=1+x,h=2/(1+x)^2,E=(1+x)^7/8, Bdiag=1-u and ell!=J; it must remain excluded.",
        "At a ramified finite fiber J may vanish; E=a/h^3 still has a pole, so nonzero u' is not assumed."
    ]}
IDENTITIES=[
"H6 universal IBP polynomial certificate and exact scalar recurrence for m>=1 only.",
"H8 local Jacobi normalization; local zero-start determinant first uses t_(m+2).",
"H10 derivative certificate proves squarefree constant continuants algebraically.",
"H11 formal common m^-2 correction distinguishes every fixed nonzero spectral value.",
"Algebraic derivative span gives a finite polynomial coefficient recurrence with a nonzero highest shift.",
"H13 endpoint congruence and residue formula give Bdiag; both purity directions are scoped theorem propositions.",
"Two omitted finite values contradict a nonconstant rational map; polynomial E excludes every Bdiag root fiber.",
"H15 rational squarefreeQ degree data plus the SECOND support argument prove infinite polynomial-x marked survival.",
"H18-H19 degree comparisons and literal general-lambda safe/unsafe realizations."
]
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--repo-root",type=Path,default=HERE,help="Explicit workspace root; this self-contained packet imports no external theorem artifacts.");ap.add_argument("--output",type=Path,default=HERE/RN)
    a=ap.parse_args();require(a.repo_root.resolve().is_dir(),"Repository root must exist");ledger=check_sources()
    ids=universal_algebra();K,scalar=scalar_moment_controls()
    result={"status":"PASS_SOURCE_GUARDED_NONCENTRAL_AFFINE_CUBIC","scope":EXPECTED,"source_ledger_sha256":sha(HERE/LN),"pins":ledger["sources"],"classical_context":PRIMARY,"logical_boundary":BOUNDARY,"written_identity_families":IDENTITIES,"universal_polynomial_identity_controls":ids,"scalar_controls":scalar,"continuant_controls":local_continuants(),"realization_controls":realization_controls(K),"support_and_boundary_controls":support_and_omission_controls(K),"source_scope_negative_controls":source_negatives(),"assurance":FLAGS}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],"mandatory_roles":len(ROLES),"source_scope_negative_controls":result["source_scope_negative_controls"],"output":str(a.output)},sort_keys=True))
if __name__=="__main__":
    try:main()
    except (ValueError,AssertionError,KeyError,TypeError,OSError) as e:print("FAIL: "+str(e),file=sys.stderr);raise SystemExit(1)
