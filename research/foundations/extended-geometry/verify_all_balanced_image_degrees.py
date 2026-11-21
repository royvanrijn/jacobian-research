#!/usr/bin/env python3
"""Source-guarded bounded exact controls for the written all-degree Image witnesses."""
import argparse,copy,hashlib,json,sys
from pathlib import Path
from fractions import Fraction as Q
from math import factorial,comb
if not __debug__:
    raise SystemExit("Optimized Python is unsupported: exact assertions must execute.")
HERE=Path(__file__).resolve().parent
PN="ALL_BALANCED_DEGREES_IMAGE_WITNESSES_AT_LEAST_THREE_PAIRS.md"
CN="verify_all_balanced_image_degrees.py"
IN="all_balanced_image_degree_inputs.json"
LN="all_balanced_image_degree_source_ledger.json"
RN="all_balanced_image_degree_verification.json"
PROOF_SHA="cd00c7e229e0bb738726608e45525e8694c93131ce38021861b553886c92b918"
SEED_PATH="results/three-pair-image-counterexample/main.tex"
SEED_SHA="0c49d6bb76957f2bc8a2159fbd70e191b526b4552945cfcb9dab163a88a18dcb"
EXPECTED={'schema': 1, 'field': 'every characteristic-zero field and every coefficient-field extension', 'pairs': {'integer': True, 'minimum': 3, 'two_pair_coverage_claim': False}, 'degree': {'integer': True, 'balanced_minimum': 2, 'ordinary_total': '2D', 'least_matrix_degree': 'D', 'degree_one_claim': False}, 'seed': {'retained_preprint_credit': True, 'journal_publication_claim': False, 'all_order_identities_reproved': True}, 'family': {'D2': 'tau(t-y)(wz+vt)', 'D_at_least3': 'tau(t-y)[R^(D-2)(wz+vt)+(wy)^(D-2)w(t-y)]', 'same_marker': 'tau*y', 'marker_balanced_degree': 1, 'marker_same_for_all_n_D_m': True}, 'moments': {'all_positive_integer_orders': True, 'pure': 0, 'mixed': '(-1)^(m-1)(m+1)!m! (n+2m+1)_(m(D-2))', 'Image_kernel_is_vector_space_not_ideal': True}, 'power': {'full_original_ring': True, 'R_primitive': True, 'tau_valuation': 1, 'absolute_scalar_nonproper_power': True, 'all_extensions': True}, 'source_boundary': {'all_math_reproved': True, 'all_Aut_theorem_logical_premise': False, 'automorphism_group_classification_logical_premise': False, 'finite_prefix_logical_premise': False, 'seed_file_relation': 'AUTHENTICATED_RETAINED_PREPRINT_CREDIT_ONLY'}, 'nonclaims': {'all_polynomials_of_these_degrees_fail': False, 'two_pair_degree5_resolved': False, 'two_pair_degree7_resolved': False, 'Jacobian_consequence': False, 'historical_priority': False}, 'verification': {'bounded_exact_controls_only': True, 'finite_controls_prove_universal_quantifiers': False, 'formal_verification': False, 'external_review': False, 'independent_replay': False}}

FLAGS={"formal_verification":False,"external_review":False,"independent_replay":False,"worldwide_priority":False}
BOUNDARY={"external_mathematical_premises":[],"all_order_math":"self-contained N1-N12, normalization, kernel and tau-valuation proof","credit_file_is_logical_premise":False,"all_Aut_theorem_is_logical_premise":False,"finite_checks_prove_all_orders":False}
ROLES=[("main_proof","here",PN,"WRITTEN_SELF_CONTAINED_UNIVERSAL_PROOF"),
("checker","here",CN,"BOUNDED_EXACT_CONTROLS_AND_SOURCE_AUTHENTICATION"),
("typed_inputs","here",IN,"EXACT_TYPED_SCOPE"),
("retained_seed_credit","repo",SEED_PATH,"AUTHENTICATED_RETAINED_PREPRINT_CREDIT_ONLY")]
CREDIT={"file":SEED_PATH,"status":"deposited preprint; not a journal-publication claim","locators":["Theorem thm:main","Lemma lem:kernel","Proof binomial identities eq:pure and eq:mixed"],"relation":"Seed and moment cancellation credited; all identities needed by this theorem are reproved in the new note."}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text(encoding="utf-8"))
def strict_equal(a,b):return json.dumps(a,sort_keys=True,separators=(",",":"))==json.dumps(b,sort_keys=True,separators=(",",":"))
def require(ok,msg):
    if not ok:raise ValueError(msg)
def validate_inputs(d):require(strict_equal(d,EXPECTED),"Exact typed scope mismatch")
def validate_ledger(d):
    require(set(d)=={"schema","sources","credit","logical_boundary","assurance"},"Ledger keys")
    require(type(d["schema"]) is int and d["schema"]==1,"Ledger schema")
    require(strict_equal(d["credit"],CREDIT),"Exact retained credit boundary")
    require(strict_equal(d["logical_boundary"],BOUNDARY),"Exact typed logical boundary")
    require(strict_equal(d["assurance"],FLAGS),"Assurance flags")
    rows=d["sources"];require(type(rows) is list and len(rows)==len(ROLES),"Mandatory source cardinality")
    got={}
    for r in rows:
        require(type(r) is dict and set(r)=={"role","root","path","sha256","relation"},"Exact source row schema")
        require(type(r["role"]) is str and r["role"] not in got,"Duplicate source role")
        require(type(r["sha256"]) is str and len(r["sha256"])==64 and all(c in "0123456789abcdef" for c in r["sha256"]),"Hash format")
        got[r["role"]]=r
    for role,root,path,relation in ROLES:
        require(role in got,"Missing role "+role);r=got[role]
        require((r["root"],r["path"],r["relation"])==(root,path,relation),"Role/path/root/relation mismatch")
    return got
def check_sources(repo):
    validate_inputs(load(HERE/IN));d=load(HERE/LN);rows=validate_ledger(d)
    for role,root,path,relation in ROLES:
        p=(HERE if root=="here" else repo)/path
        require(p.is_file(),"Missing artifact "+role)
        require(sha(p)==rows[role]["sha256"],"Artifact hash mismatch "+role)
    require(sha(HERE/PN)==PROOF_SHA,"Reviewed proof body")
    require(sha(repo/SEED_PATH)==SEED_SHA,"Retained actual seed credit")
    require(b"\r" not in (HERE/PN).read_bytes(),"Proof uses LF")
    return d
def source_negatives():
    d=load(HERE/IN);l=load(HERE/LN);n=0
    def reject(fn,x):
        nonlocal n
        try:fn(x)
        except (ValueError,KeyError,TypeError):n+=1;return
        raise AssertionError("Negative source/scope control accepted")
    cases=[("pairs","minimum",2),("pairs","integer",1),("pairs","two_pair_coverage_claim",True),
    ("degree","balanced_minimum",1),("degree","ordinary_total","D"),("degree","degree_one_claim",True),
    ("seed","journal_publication_claim",True),("seed","all_order_identities_reproved",False),
    ("family","marker_same_for_all_n_D_m",False),("family","marker_balanced_degree",True),
    ("moments","all_positive_integer_orders",False),("power","full_original_ring",False),
    ("power","R_primitive",False),("power","tau_valuation",True),("power","all_extensions",False),
    ("source_boundary","all_Aut_theorem_logical_premise",True),
    ("nonclaims","two_pair_degree5_resolved",True),("nonclaims","two_pair_degree7_resolved",True),
    ("nonclaims","all_polynomials_of_these_degrees_fail",True),
    ("verification","finite_controls_prove_universal_quantifiers",True),
    ("verification","external_review",True)]
    for a,b,v in cases:
        z=copy.deepcopy(d);z[a][b]=v;reject(validate_inputs,z)
    for i in range(len(ROLES)):
        z=copy.deepcopy(l);z["sources"].pop(i);reject(validate_ledger,z)
    z=copy.deepcopy(l);z["sources"].append(z["sources"][0]);reject(validate_ledger,z)
    z=copy.deepcopy(l);z["schema"]=True;reject(validate_ledger,z)
    z=copy.deepcopy(l);z["logical_boundary"]["credit_file_is_logical_premise"]=0;reject(validate_ledger,z)
    z=copy.deepcopy(l);z["logical_boundary"]["external_mathematical_premises"]=["all-Aut theorem"];reject(validate_ledger,z)
    z=copy.deepcopy(l);z["credit"]["status"]="journal publication";reject(validate_ledger,z)
    z=copy.deepcopy(l);z["assurance"]["formal_verification"]=True;reject(validate_ledger,z)
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


def degree_distribution(p,inds):return sorted({sum(e[j] for j in inds) for e in p})
def evaluate(p,vals):
    return sum((c*__import__("math").prod(Q(v)**j for v,j in zip(vals,e)) for e,c in p.items()),Q(0))
def literal_seed_binomial_controls():
    counts=0
    for m in range(1,81):
        pure=0;mixed=0
        for s in range(m+1):
            a=(-1)**s*comb(m,s)**2*factorial(m)*factorial(s)*factorial(m-s)
            b=(-1)**s*comb(m,s)*factorial(m)**2
            require(a==b,"Literal pure seed factorial cancellation");pure+=a;counts+=1
        for s in range(m):
            a=(-1)**s*comb(m,s)*comb(m,s+1)*factorial(m+1)*factorial(s+1)*factorial(m-s-1)
            b=(-1)**s*comb(m,s)*factorial(m+1)*factorial(m)
            require(a==b,"Literal marked seed factorial cancellation");mixed+=a;counts+=1
        require(pure==0 and mixed==(-1)**(m-1)*factorial(m+1)*factorial(m),"Full/truncated seed binomial sums");counts+=2
    return counts

def original_family_controls():
    family=flow=primitive=nonpower=0;records=[]
    for n in range(3,6):
        P=Ring(2*n);tau,w,v=[P.v(i) for i in range(3)];t,z,y=[P.v(n+i) for i in range(3)]
        R=add(*(mul(P.v(i),P.v(n+i)) for i in range(n)))
        seed=mul(mul(tau,add(t,scale(y,-1))),add(mul(w,z),mul(v,t)))
        h=mul(w,y);marker=mul(tau,y)
        def L(a):return add(mul(w,diff(a,2)),scale(mul(y,diff(a,n+1)),-1))
        lf=mul(mul(tau,w),power(add(t,scale(y,-1)),2))
        require(L(seed)==lf and L(lf)=={} and L(R)=={} and L(h)=={} and L(marker)=={},"Exact replica identities");flow+=5
        for D in range(2,7):
            k=D-2
            F=seed if D==2 else add(mul(power(R,k),seed),mul(power(h,k),lf))
            require(degree_distribution(F,range(n))==[D] and degree_distribution(F,range(n,2*n))==[D],"Literal original balanced degree")
            require(max(sum(e) for e in F)==2*D,"Ordinary total degree");family+=2
            # Literal R=0 evaluations supply exact nondivisibility controls.
            vals=[Q(0)]*(2*n)
            if D==2:
                vals[0]=1;vals[2]=1;vals[n]=1;vals[n+2]=-1
                require(evaluate(R,vals)==0 and evaluate(F,vals)==2,"Seed primitive point")
            else:
                vals[0]=1;vals[1]=1;vals[n+2]=1
                require(evaluate(R,vals)==0 and evaluate(F,vals)==1,"Higher-degree primitive point")
            primitive+=1
            # F has tau factor and nonzero first tau coefficient, over every extension.
            require(all(e[0]>=1 for e in F),"tau divides F")
            vals=[Q(0)]*(2*n);vals[2]=1;vals[n]=2;vals[n+2]=1
            require(evaluate(diff(F,0),vals)==2,"Simple tau valuation");nonpower+=1
            for m in range(1,4):
                fpower=power(F,m);actual=original_E(mul(marker,fpower),n)
                target=(-1)**(m-1)*factorial(m+1)*factorial(m)*rising(n+2*m+1,m*(D-2))
                require(original_E(fpower,n)=={} and actual=={(0,)*n:Q(target)},"Actual original pure/mixed family contraction")
                records.append({"n":n,"D":D,"m":m,"mixed_integer":str(target),"terms":len(fpower)})
    return {"original_family_contractions":len(records),"records":records,"homogeneous_degree_controls":family,"exact_flow_controls":flow,"R_primitive_controls":primitive,"simple_tau_controls":nonpower}

def normalization_replica_kernel_controls():
    norm=replica=kernel=marker=0
    for n in range(3,6):
        P=Ring(2*n);w,v=P.v(1),P.v(2);z,y=P.v(n+1),P.v(n+2)
        R=add(*(mul(P.v(i),P.v(n+i)) for i in range(n)));h=mul(w,y)
        def L(a):return add(mul(w,diff(a,2)),scale(mul(y,diff(a,n+1)),-1))
        for D in range(5):
            cs=list(compositions(D,n))
            for ia,a in enumerate(cs):
                for ib,b in enumerate(cs):
                    H={a+b:Q(1)}
                    require(original_E(L(H),n)=={},"N7 literal infinitesimal factorial cancellation");replica+=1
                    if (ia+ib)%max(1,len(cs)//5)==0:
                        for k in [1,2,3]:
                            require(original_E(mul(power(h,k),L(H)),n)=={},"N8 replica factorial cancellation");replica+=1
                # Independent N4 normalization, including several radial powers.
                b=cs[(ia*3+1)%len(cs)];H={a+b:Q(1)};EH=original_E(H,n)
                for j in range(4):
                    require(original_E(mul(power(R,j),H),n)==scale(EH,rising(n+D,j)),"N4 radial Pochhammer normalization");norm+=1
        for q in range(60):
            e=tuple((q*(j+1)+j*j)%4 for j in range(2*n));H={e:Q(1)}
            for i in range(n):
                d=add(diff(H,n+i),scale(mul(P.v(i),H),-1))
                require(original_E(d,n)=={},"Literal original Image kernel");kernel+=1
            a=tuple((q+j)%3 for j in range(n));b=tuple(a[j]+((q+j*j)%2) for j in range(n))
            H={a+b:Q(1)};eta=tuple(b[j]-a[j] for j in range(n));lhs=original_E(mul({eta+(0,)*n:Q(1)},H),n);rhs=original_E(H,n)
            for j,order in enumerate(eta):
                for r in range(order):rhs=diff(rhs,j)
            require(lhs==rhs,"Balanced original marker derivative");marker+=1
    return {"radial_normalization_controls":norm,"literal_replica_controls":replica,"original_kernel_controls":kernel,"marker_balancing_controls":marker}

def meaningful_negatives():
    P=Ring(6);tau,w,v,t,z,y=[P.v(i) for i in range(6)]
    R=add(mul(tau,t),mul(w,z),mul(v,y))
    seed=mul(mul(tau,add(t,scale(y,-1))),add(mul(w,z),mul(v,t)))
    marker=mul(tau,y);h=mul(w,y);lf=mul(mul(tau,w),power(add(t,scale(y,-1)),2))
    F3=add(mul(R,seed),mul(h,lf))
    require(original_E(mul(marker,F3),3)=={(0,0,0):Q(12)} and original_E(mul(marker,seed),3)=={(0,0,0):Q(2)},"Radial Pochhammer omission is false")
    # Radial padding alone retains these moments but is not primitive.
    padded=mul(R,seed)
    require(original_E(mul(marker,padded),3)==original_E(mul(marker,F3),3),"Padding has same first marked moment")
    vals=[1,1,0,0,0,1]
    require(evaluate(R,vals)==0 and evaluate(padded,vals)==0 and evaluate(F3,vals)==1,"Moments alone do not establish primitivity")
    require(original_E(mul(y,seed),3)=={(1,0,0):Q(2),(0,0,1):Q(-1)},"Omitting tau leaves an unbalanced nonconstant contraction")
    require(original_E(tau,3)=={} and original_E(mul(tau,t),3)=={(0,0,0):Q(1)},"Original Image kernel is not an ideal")
    return ["Pochhammer normalization cannot be omitted.","Radial padding alone can share moments while failing primitivity.","The literal balanced marker is tau*y; y alone has nonconstant contraction.","The original Image subspace is not an ideal."]

IDENTITIES=[
"N2-N3: matching pure indices j=s and marked indices j=s+1 reduce to exact full and truncated binomial sums.",
"N4: E_n(R H)=(n+D) E_n(H) on balanced degree D; repeated radial powers give a rising factorial.",
"N7: alpha_v*beta!=beta_z*alpha! under beta=alpha+e_w-e_v.",
"N8: L(h)=0 implies h^k L(a)=L(h^k a); each derivative term has zero trace.",
"N10-N11: L(F2)=tau*w*(t-y)^2, L^2(F2)=0, and the balanced lift restricts to exp(h^k L)f2.",
"N1: (n)_(mD+1)/(n)_(2m+1)=(n+2m+1)_(m(D-2)).",
"N12: reduction modulo R leaves a nonzero product of displayed linear factors; literal R=0 points are additional controls.",
"Simple tau valuation one excludes every scalar proper power over all coefficient extensions.",
"Monomial reduction by partial_z-xi and E annihilation prove the original Image kernel identity."
]
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--repo-root",type=Path,required=True);ap.add_argument("--output",type=Path,default=HERE/RN)
    args=ap.parse_args();repo=args.repo_root.resolve();require(repo.is_dir(),"Repository root must exist")
    ledger=check_sources(repo)
    result={"status":"PASS_SOURCE_GUARDED_ALL_BALANCED_IMAGE_DEGREES","scope":EXPECTED,"source_ledger_sha256":sha(HERE/LN),"pins":ledger["sources"],"logical_boundary":BOUNDARY,"credit":CREDIT,
    "written_identity_families":IDENTITIES,"seed_binomial_controls":literal_seed_binomial_controls(),"family_controls":original_family_controls(),
    "normalization_and_replica_controls":normalization_replica_kernel_controls(),"mathematical_negative_controls":meaningful_negatives(),
    "source_scope_negative_controls":source_negatives(),"assurance":FLAGS}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],"mandatory_roles":len(ROLES),"output":str(args.output)},sort_keys=True))
if __name__=="__main__":
    try:main()
    except (ValueError,AssertionError,TypeError,KeyError,OSError) as e:print("FAIL: "+str(e),file=sys.stderr);raise SystemExit(1)
