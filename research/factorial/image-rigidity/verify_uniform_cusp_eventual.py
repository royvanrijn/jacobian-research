"""Exact/source support for uniform integer-cusp eventual Strong Factorial.

Finite algebra controls authenticate identities and inputs; analytic complete
contours and uniformity are proved in the accompanying written sources.
"""
if not __debug__:
    raise SystemExit("Optimized Python execution is refused.")
import argparse
import copy
from fractions import Fraction as F
import hashlib
import json
from math import comb, factorial
from pathlib import Path

HERE=Path(__file__).resolve().parent
EXPECTED_INPUT={
    "schema":1,
    "scope":{
        "one_cutoff_all_integer_cusps":True,
        "coefficient_and_characteristic_zero_field_uniform":True,
        "global_bad_pair_set_finite":True,
        "cutoff_evaluated":False,
        "finite_box_evaluated":False,
        "all_order_sfc_claimed":False,
        "exact_finite_gamma_circle_claimed":False,
        "small_m_over_s_assumed":False,
        "uses_407_main":False
    },
    "finite_controls":{
        "moment_s_max":6,"moment_order_max":6,
        "coefficient_pairs":[[-2,1],[1,1],[3,2],[0,1],[1,0]],
        "jet_s_max":20,"annulus_s_max":80,
        "box_N0_max":24,"box_s_max":64,
        "local_rho":["1/10","1/20","1/50"],
        "local_k":[20,100,1000]
    }
}
ASSURANCE={"independent_replay":False,"external_review":False,"formal_verification":False}
# Each expected role has one immutable kind/path/relation. The ledger pins the
# actual current bytes. All proof-level imports exclude the 407 main theorem.
SPECS={
    'main_proof': ('local', 'UNIFORM_INTEGER_CUSP_EVENTUAL_STRONG_FACTORIAL.md', 'written_argument'),
    'uniform_local_proof': ('local', 'UNIFORM_GROWING_QUARTIC_CUSP_NEIGHBORHOOD.md', 'mathematical_input_uniform_local'),
    'uniform_global_proof': ('local', 'UNIFORM_CUSP_CURIE_WEISS_REGULAR_COMPLEMENT.md', 'mathematical_input_uniform_regular'),
    'exact_checker': ('local', 'verify_uniform_cusp_eventual.py', 'exact_and_source_support'),
    'exact_inputs': ('local', 'uniform_cusp_eventual_inputs.json', 'exact_fixture_input'),
    'context_audit': ('local', 'UNIFORM_CUSP_PRIMARY_CONTEXT.md', 'classical_context_only'),
    'primary_kabluchko': ('local', 'uniform_cusp_sources/kabluchko_2203_05533v3.pdf', 'classical_context_only'),
    'compact_402_proof': ('repo', 'research/factorial/image-rigidity/UNIFORM_BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md', 'mathematical_input_compact_family_and_bounded_H'),
    'fixed_cusp_405_proof': ('repo', 'research/factorial/image-rigidity/CUSP_EVENTUAL_STRONG_FACTORIAL_THEOREM.md', 'mathematical_input_narrow_bounded_s_route'),
    'fixed_local_construction': ('repo', 'research/factorial/image-rigidity/CUSP_SHRINKING_WINDOW_LOG_RATIO_NONCOLLISION.md', 'mathematical_input_C1_to_C15_revalidated_uniformly'),
    'fixed_global_complement': ('repo', 'research/factorial/image-rigidity/CUSP_GLOBAL_REGULAR_COMPLEMENT_AND_STOKES_DETOUR.md', 'mathematical_input_narrow_405_closure'),
    'fixed_morse_proof': ('repo', 'research/factorial/image-rigidity/CUSP_TWO_THIMBLE_POSITIVE_CYCLE.md', 'mathematical_input_narrow_405_closure'),
    'cusp_stationary_geometry': ('repo', 'research/factorial/image-rigidity/CUSP_STATIONARY_BELYI_MONODROMY_AND_RADIAL_BARRIER.md', 'mathematical_input_narrow_405_closure'),
    'actual_coalescing_bridge': ('repo', 'research/factorial/image-rigidity/COALESCING_BINOMIAL_SADDLE_BRIDGE.md', 'mathematical_input_actual_phase_and_jets'),
    'fixed_bounded_proof': ('repo', 'research/factorial/image-rigidity/BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md', 'mathematical_input_402_model_mechanism'),
    'quartic_unfolding_proof': ('repo', 'research/factorial/image-rigidity/QUARTIC_CUSP_UNFOLDING_AND_ACTUAL_ORDER_DRIFT.md', 'mathematical_input_402_actual_jet_mechanism'),
    'quartic_model_proof': ('repo', 'research/factorial/image-rigidity/QUARTIC_LEE_YANG_ZERO_AND_ROBIN_BRIDGE.md', 'mathematical_input_model_location_and_simplicity_only'),
    'primary_newman': ('repo', 'research/factorial/image-rigidity/quartic_lee_yang_sources/newman_fourier_transforms_1976.pdf', 'mathematical_input_model_location_only'),
    'old_gamma_upper_soft': ('repo', 'research/factorial/image-rigidity/GAMMA_PRODUCT_ROOT_MODULUS.tex', 'mathematical_input_old_upper_and_local_catch_soft'),
    'sharp_annulus_soft_404': ('repo', 'research/factorial/image-rigidity/GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md', 'mathematical_input_compact_annulus_and_finite_bad_box'),
    'closed_disc_grace': ('repo', 'research/factorial/image-rigidity/GAP_TWO_GRACE_GLOBAL_LOCALIZATION.md', 'mathematical_input_404_closed_disc_mechanism_only'),
}
LOGICAL_ROLES=[
    'uniform_local_proof',
    'uniform_global_proof',
    'compact_402_proof',
    'fixed_cusp_405_proof',
    'fixed_local_construction',
    'fixed_global_complement',
    'fixed_morse_proof',
    'cusp_stationary_geometry',
    'actual_coalescing_bridge',
    'fixed_bounded_proof',
    'quartic_unfolding_proof',
    'quartic_model_proof',
    'primary_newman',
    'old_gamma_upper_soft',
    'sharp_annulus_soft_404',
    'closed_disc_grace',
]

def need(condition,message):
    if not condition:
        raise ValueError(message)

def exact_equal(actual,expected,label="input"):
    need(type(actual) is type(expected),label+" has the wrong exact type")
    if isinstance(expected,dict):
        need(set(actual)==set(expected),label+" has wrong keys")
        for k in expected: exact_equal(actual[k],expected[k],label+"."+k)
    elif isinstance(expected,list):
        need(len(actual)==len(expected),label+" has wrong length")
        for i,(a,b) in enumerate(zip(actual,expected)):
            exact_equal(a,b,label+"["+str(i)+"]")
    else:
        need(actual==expected,label+" has the wrong value")

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def check_header(ledger):
    need(type(ledger) is dict,"ledger must be an object")
    need(set(ledger)=={"schema","sources","logical_inputs","assurance_flags"},"ledger has wrong keys")
    need(type(ledger.get("schema")) is int and ledger["schema"]==1,"ledger schema must be integer1")
    exact_equal(ledger.get("assurance_flags"),ASSURANCE,"assurance")
    exact_equal(ledger.get("logical_inputs"),LOGICAL_ROLES,"logical_inputs")
    need(type(ledger.get("sources")) is list,"source records must be a list")
    roles=[]
    for rec in ledger["sources"]:
        need(type(rec) is dict and set(rec)=={"role","kind","path","relation","sha256"},"malformed source record")
        role=rec["role"]
        need(type(role) is str and role in SPECS,"unknown source role")
        need(role not in roles,"duplicate source role")
        roles.append(role)
        kind,path,relation=SPECS[role]
        exact_equal(rec["kind"],kind,"source kind")
        exact_equal(rec["path"],path,"source path")
        exact_equal(rec["relation"],relation,"source relation")
        digest=rec["sha256"]
        need(type(digest) is str and len(digest)==64 and set(digest)<=set("0123456789abcdef"),"malformed hash")
    need(set(roles)==set(SPECS),"every source role is mandatory")
    return ledger["sources"]

def validate_sources(ledger,repo):
    records=check_header(ledger)
    out={}
    for rec in records:
        base=HERE if rec["kind"]=="local" else repo
        path=(base/rec["path"]).resolve()
        need(path.is_relative_to(base.resolve()),"source path escapes its root")
        need(path.is_file(),"missing source: "+rec["role"])
        digest=sha(path)
        need(digest==rec["sha256"],"source hash mismatch: "+rec["role"])
        out[rec["role"]]=digest
    return out

# Small exact sparse polynomial/rational-function kernel. No CAS is imported.
VARS=("t","q","z","h","r","s","l","p","H","b","lam","y","x","d","m","rho","k","delta")
DIM=len(VARS); ZERO=(0,)*DIM
def pc(a):
    a=F(a)
    return {ZERO:a} if a else {}
def padd(a,b):
    c=dict(a)
    for e,v in b.items():
        c[e]=c.get(e,F(0))+v
        if not c[e]: del c[e]
    return c
def pneg(a): return {e:-v for e,v in a.items()}
def pmul(a,b):
    c={}
    for e,v in a.items():
        for f,w in b.items():
            g=tuple(x+y for x,y in zip(e,f))
            c[g]=c.get(g,F(0))+v*w
    return {e:v for e,v in c.items() if v}
def ppow(a,n):
    c=pc(1)
    while n:
        if n&1:c=pmul(c,a)
        a=pmul(a,a);n//=2
    return c
def pdiff(a,name):
    i=VARS.index(name);c={}
    for e,v in a.items():
        if e[i]:
            f=list(e);f[i]-=1;c[tuple(f)]=v*e[i]
    return c

class RF:
    def __init__(self,n=0,d=None):
        if isinstance(n,RF) and d is None:
            self.n,self.d=n.n,n.d
        else:
            self.n=n if isinstance(n,dict) else pc(n)
            self.d=pc(1) if d is None else d
            need(bool(self.d),"zero rational-function denominator")
    def __add__(self,o):
        o=RF(o);return RF(padd(pmul(self.n,o.d),pmul(o.n,self.d)),pmul(self.d,o.d))
    __radd__=__add__
    def __neg__(self):return RF(pneg(self.n),self.d)
    def __sub__(self,o):return self+-RF(o)
    def __rsub__(self,o):return RF(o)+-self
    def __mul__(self,o):
        o=RF(o);return RF(pmul(self.n,o.n),pmul(self.d,o.d))
    __rmul__=__mul__
    def __truediv__(self,o):
        o=RF(o);need(bool(o.n),"division by zero");return RF(pmul(self.n,o.d),pmul(self.d,o.n))
    def __rtruediv__(self,o):return RF(o)/self
    def __pow__(self,n):
        if n<0:return RF(ppow(self.d,-n),ppow(self.n,-n))
        return RF(ppow(self.n,n),ppow(self.d,n))
    def diff(self,name):
        return RF(padd(pmul(pdiff(self.n,name),self.d),pneg(pmul(self.n,pdiff(self.d,name)))),ppow(self.d,2))
    def iszero(self):return not self.n
    def subst(self,**values):
        def evalp(poly):
            ans=RF()
            for e,v in poly.items():
                term=RF(v)
                for name,power in zip(VARS,e):
                    if power:term*=RF(values.get(name,var(name)))**power
                ans+=term
            return ans
        return evalp(self.n)/evalp(self.d)
def var(name):
    e=list(ZERO);e[VARS.index(name)]=1
    return RF({tuple(e):F(1)})
def identity(label,expression,out):
    need(expression.iszero(),"universal identity failed: "+label)
    out.append(label)

# Truncated ordinary power series, exact over Q.
def sadd(a,b,n):
    return [(a[i] if i<len(a) else F(0))+(b[i] if i<len(b) else F(0)) for i in range(n+1)]
def smul(a,b,n):
    out=[F(0)]*(n+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<=n:out[i+j]+=x*y
    return out
def slog_unit(a,n):
    need(a[0]==1,"series log requires constant1")
    u=list(a);u[0]=0
    out=[F(0)]*(n+1);power=[F(1)]+[F(0)]*n
    for j in range(1,n+1):
        power=smul(power,u,n)
        out=sadd(out,[F((-1)**(j+1),j)*x for x in power],n)
    return out
def sinverse(a,n):
    need(a[0]!=0,"series inverse requires nonzero constant")
    out=[1/a[0]]+[F(0)]*n
    for j in range(1,n+1):
        out[j]=-sum(a[k]*out[j-k] for k in range(1,min(j,len(a)-1)+1))/a[0]
    return out


def universal_checks():
    t,q,z,h,r,s,l,p,H,b,lam,y,x,d,m,rho,k,delta=(var(n) for n in VARS)
    out=[]
    def dt(f):return f.diff("t")+r*f.diff("r")
    def dz(f):return f.diff("z")+2*r*f.diff("r")
    Q=r*r*(1-t)/(1+t)
    identity("limiting stationary-map derivative",dt(Q)+2*t*t*r*r/(1+t)**2,out)
    identity("critical factor at a stationary point",Q/r+r-2*r/(1+t),out)
    identity("reciprocal stationary map",Q.subst(t=-t,r=1/r)-1/Q,out)
    identity("reciprocal entire factor",(q/r+r).subst(t=-t,q=1/q,r=1/r)*q-(q/r+r),out)
    sh=(r-1/r)/2;ch=(r+1/r)/2
    E0=(1+2*z*z)*sh-2*z*ch
    identity("positive horizontal primitive derivative",dz(E0)-4*z*z*ch,out)
    Th=sh/(z-h)-ch
    identity("tilted horizontal monotonicity numerator",
             -dz(Th)*(z-h)**2-((1+2*(z-h)**2)*sh-2*(z-h)*ch),out)
    # V''=-1+2(1+C cosh)/(cosh+C)^2 at x=-z.
    V2=-1+2*(1+Th*ch)/(ch+Th)**2
    identity("horizontal Gaussian curvature at the favored maximum",
             V2-(z-h)*dz(Th)/(ch+Th),out)
    tail0=s*s-l*l*p+(2*s+1)*H
    tail1=-s-l*t+(2*s+1)*H
    identity("complete-tail numerator modulo L squared equals s(s+1)",
             (tail0-tail1).subst(t=l*(p-1))-(s*s+s-l*l),out)
    identity("exact normalized eta path",
             b*(lam*q-1)-(lam*b*(q-1)+b*(lam-1)),out)
    ratio=(lam*x-b*(lam-1))/(lam*y)
    identity("exact endpoint cone ratio",ratio-(x/y-b*(1-1/lam)/y),out)
    identity("endpoint ratio derivative sign formula",ratio.diff("lam")+b/(lam*lam*y),out)
    identity("actual adjacent bare-action envelope",
             -delta+m*(-l*x/(d*q))*(d*q/m)-(-delta-l*x),out)
    identity("first local transported error after exact m=k/rho^4",
             ((d/m)*(rho**-2)/(l*rho)).subst(m=k/rho**4)-(d/l)*rho/k,out)
    identity("second local transported error after exact m=k/rho^4",
             ((d/m)/(k*rho**3)/(l*rho)).subst(m=k/rho**4)-(d/l)/(k*k),out)
    identity("critical CW and unitary-Hermite logarithmic weights",
             (2*k-m)**2/(2*m)-(m/2-2*k*(m-k)/m),out)
    # Nonzero remainders detect real corruptions, not only mutated flags.
    need(not (dt(Q)+2*t*r*r/(1+t)**2).iszero(),"wrong derivative was accepted")
    need(not (V2+(z-h)*dz(Th)/(ch+Th)).iszero(),"wrong Gaussian curvature sign was accepted")
    need(not ((2*k-m)**2/(2*m)+2*k*(m-k)/m).iszero(),"missing CW scalar was accepted")
    need(not (q/r+r).subst(q=-1).iszero(),"global q+1 normalization corruption was accepted")
    positive=[]
    for j in range(16):
        coefficient=F(4**(j+1),factorial(2*j)*(2*j+3))
        need(coefficient>0,"horizontal positive coefficient is not positive")
        positive.append(str(coefficient))
    return out,positive,4

def local_jet_checks(max_s):
    records=[]
    n=6
    for s in range(1,max_s+1):
        a=s*s;d=2*s+1;C=s*(s+1);kap=F(s+1,s)
        power=[F(comb(d,j)) if j<=d else F(0) for j in range(n+1)]
        A=[F(0)]+[F(a*(-1)**(j+1),j) for j in range(1,n+1)]
        A[1]-=C
        denom=power[:];denom[0]+=kap
        normalized=[v/denom[0] for v in denom]
        phase=sadd(A,slog_unit(normalized,n),n)
        need(all(phase[j]==0 for j in [0,1,2,3]),"actual cusp lower jet failed")
        need(phase[4]==-F(C*C,12),"actual normalized quartic jet failed")
        inv=sinverse(denom,n)
        Fq=[kap*x for x in inv];Fq[0]-=kap/denom[0]
        need(Fq[1]==-F(C,d),"actual mixed field jet failed")
        need(Fq[2]==0,"raw mixed quadratic jet failed")
        records.append({"s":s,"fourth_z_coefficient":str(phase[4]),
                        "first_z_q_coefficient":str(Fq[1])})
    return records

def moment_checks(inp):
    count=0;coefficient_count=0;monomial_count=0
    for s in range(1,inp["moment_s_max"]+1):
        a=s*s;d=2*s+1
        for m in range(1,inp["moment_order_max"]+1):
            am=a*m
            coefficients=[]
            previous=1
            for j in range(m+1):
                direct=comb(m,j)*factorial(am+d*j)//factorial(am)
                if j:
                    for z in range(am+d*(j-1)+1,am+d*j+1):previous*=z
                sliding=comb(m,j)*previous
                need(direct==sliding,"actual factorial/rising coefficient mismatch")
                coefficients.append(direct);coefficient_count+=1
            for A,B in inp["coefficient_pairs"]:
                actual=sum(comb(m,j)*A**(m-j)*B**j*factorial(am+d*j) for j in range(m+1))
                if B:
                    value=sum(F(c)*F(A,B)**(m-j) for j,c in enumerate(coefficients))
                    normalized=F(B**m*factorial(am))*value
                else:
                    normalized=A**m*factorial(am);monomial_count+=1
                need(actual==normalized,"coefficient/field moment normalization failed")
                count+=1
    return coefficient_count,count,monomial_count

def finite_controls(inp):
    annulus=[]
    for s in range(1,inp["annulus_s_max"]+1):
        lo=F(s,s+1)**(2*s+2);hi=F(s+1,s)**(2*s)
        need(lo>=F(1,16) and hi<8,"cusp annulus control failed")
        need((2*s+1)**2<=F(9,2)*s*(s+1),"d/L control failed")
        annulus.append(s)
    box=[]
    for N0 in range(4,inp["box_N0_max"]+1):
        for N in range(2,N0):
            for s in range(1,inp["box_s_max"]+1):
                hard=s*s<(2*s+1)*(N-2)
                if hard:
                    need(3<=N<N0 and 1<=s<=2*N0-6,"unevaluated-box control failed")
                    box.append([s,N,N0])
    # A wrong bound s<N can fail even within the exact hard region.
    need(4*4<(2*4+1)*(4-2) and not 4<4,"wrong s<N box corruption did not fail")
    # An exact finite-Gamma circle law would require D2^3=D1^3 D3.
    circle=[]
    def rise(z,n):
        ans=1
        for j in range(n):ans*=z+j
        return ans
    for s in range(1,13):
        a=s*s;d=2*s+1
        D1,D2,D3=[rise(3*a+1,d*j) for j in [1,2,3]]
        need(D2**3>D1**3*D3,"finite-Gamma circle obstruction was lost")
        circle.append(s)
    # Local scale error budgets, with harmless dimensionless constants1,
    # support exact exponent bookkeeping; they do not certify analytic C.
    budgets=[]
    for rho in inp["local_rho"]:
        rho=F(rho)
        for k in inp["local_k"]:
            err=rho/F(k)+F(1,k*k)
            need(0<err<1,"local exact exponent budget control failed")
            budgets.append([str(rho),k,str(err)])
    return annulus,box,circle,budgets,2

def source_negative_controls(ledger,repo,inputs):
    count=0
    def reject(fn):
        nonlocal count
        try:fn()
        except (ValueError,KeyError,TypeError):count+=1
        else:raise ValueError("a deliberate source/scope corruption was accepted")
    for role in SPECS:
        mutant=copy.deepcopy(ledger)
        mutant["sources"]=[r for r in mutant["sources"] if r["role"]!=role]
        reject(lambda m=mutant:check_header(m))
        mutant=copy.deepcopy(ledger)
        next(r for r in mutant["sources"] if r["role"]==role)["sha256"]="0"*64
        reject(lambda m=mutant:validate_sources(m,repo))
    mutant=copy.deepcopy(ledger);mutant["schema"]=True
    reject(lambda:check_header(mutant))
    mutant=copy.deepcopy(ledger);mutant["sources"].append(copy.deepcopy(mutant["sources"][0]))
    reject(lambda:check_header(mutant))
    mutant=copy.deepcopy(ledger);mutant["sources"]=[]
    reject(lambda:check_header(mutant))
    for key in EXPECTED_INPUT["scope"]:
        mutant=copy.deepcopy(inputs);mutant["scope"][key]=not mutant["scope"][key]
        reject(lambda m=mutant:exact_equal(m,EXPECTED_INPUT))
    mutant=copy.deepcopy(inputs);mutant["schema"]=True
    reject(lambda:exact_equal(mutant,EXPECTED_INPUT))
    mutant=copy.deepcopy(inputs);mutant["scope"]["uses_407_main"]=0
    reject(lambda:exact_equal(mutant,EXPECTED_INPUT))
    mutant=copy.deepcopy(inputs);mutant["evaluated_N0"]=4
    reject(lambda:exact_equal(mutant,EXPECTED_INPUT))
    return count

def find_repo(given):
    if given:return Path(given).resolve()
    for p in HERE.parents:
        if (p/"research/factorial/image-rigidity/GAMMA_PRODUCT_ROOT_MODULUS.tex").is_file():
            return p
    raise ValueError("Use --repo-root to locate repository proof inputs.")

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root")
    parser.add_argument("--source-ledger",type=Path,default=HERE/"uniform_cusp_eventual_source_ledger.json")
    parser.add_argument("--output",type=Path,default=HERE/"uniform_cusp_eventual_verification.json")
    args=parser.parse_args()
    repo=find_repo(args.repo_root)
    ledger=json.loads(args.source_ledger.read_text())
    inputs=json.loads((HERE/"uniform_cusp_eventual_inputs.json").read_text())
    exact_equal(inputs,EXPECTED_INPUT)
    pins=validate_sources(ledger,repo)
    universal,positive,neg1=universal_checks()
    jets=local_jet_checks(inputs["finite_controls"]["jet_s_max"])
    coefficients,moments,monomials=moment_checks(inputs["finite_controls"])
    annulus,box,circle,budgets,neg2=finite_controls(inputs["finite_controls"])
    negatives=source_negative_controls(ledger,repo,inputs)
    result={
        "schema":1,"status":"PASS",
        "scope":"One unevaluated eventual cutoff for every integer cusp, all coefficient pairs and characteristic-zero fields; equivalent global bad-pair finiteness. No all-order theorem or evaluated finite box.",
        "analytic_proof_notice":"Exact/source checks support the written proofs; finite samples do not certify complete contours, analytic estimates or uniformity.",
        "source_ledger_sha256":sha(args.source_ledger),
        "source_pins":pins,
        "logical_inputs":LOGICAL_ROLES,
        "excluded_routes":["407_main","broad_405_square_root_gap","finite_terminal_prefix","exact_finite_gamma_circle_law"],
        "counts":{"source_pins":len(pins),"universal_identity_families":len(universal),
                  "positive_horizontal_coefficients":len(positive),"actual_cusp_jet_records":len(jets),
                  "two_source_actual_coefficients":coefficients,"actual_moment_normalizations":moments,
                  "monomial_axis_controls":monomials,"annulus_and_d_over_L_controls":len(annulus),
                  "finite_soft_box_controls":len(box),"finite_gamma_circle_obstructions":len(circle),
                  "local_error_exponent_controls":len(budgets),
                  "genuine_math_negative_controls":neg1+neg2,
                  "source_scope_negative_controls":negatives},
        "universal_identities":universal,"jet_records":jets,
        "positive_horizontal_coefficients":positive,
        "assurance_flags":ASSURANCE
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print(json.dumps({"status":"PASS","counts":result["counts"]},sort_keys=True))
if __name__=="__main__":
    main()

