#!/usr/bin/env python3
"""Exact support for actual endpoint/bulk simplicity and the high-ratio corollary.

The two retained written proofs supply the holomorphic contour estimates.
Finite controls prove only their stated algebra and rational inequalities.
The fresh parent replay authenticates the retained numerical input closure;
its coprimality conclusion is not a simplicity premise.
"""
from __future__ import annotations
if not __debug__:
    raise SystemExit("Refusing -O/-OO: exact checks must execute.")
from pathlib import Path
from fractions import Fraction as Q
from math import comb, factorial
from hashlib import sha256
import argparse, copy, json, subprocess, sys, tempfile
HERE = Path(__file__).resolve().parent
PREFIX = "research/factorial/image-rigidity/"
CHECKER = "verify_uniform_actual_high_ratio_simplicity.py"
INPUT = "uniform_actual_high_ratio_simplicity_inputs.json"
LEDGER = "uniform_actual_high_ratio_simplicity_source_ledger.json"
RECEIPT = "uniform_actual_high_ratio_simplicity_verification.json"
# Frozen reviewed proof hashes; a re-signed changed proof is not sufficient.
LOCAL_PROOFS = {
 "endpoint_proof": ("UNIFORM_ACTUAL_ENDPOINT_ROOT_SIMPLICITY.md",
  "799db033a60dfdd7bc6383252b894f8a4b0c001ce240dbc83858eb93956e97aa"),
 "bulk_proof": ("UNIFORM_ACTUAL_BULK_ROOT_SIMPLICITY.md",
  "4a2f09edc09fcfab3dfdcec9beab321c68e9d656ee67172249afc904f033633a")}
EXPECTED = {'schema': 1, 'scope': {'integer_a_min': 1, 'integer_d_min': 2, 'integer_m_min': 10000000000000000, 'derivative_delta_max': '1/10000000000', 'all_root_location_delta_max': '1/1000000000000', 'all_root_corollary_a_min': '10^24*d^2', 'endpoint_abs_real_max': '1/10000000', 'endpoint_abs_imag_min': '19/10', 'endpoint_abs_imag_max': '21/10', 'bulk_abs_real_max': '1/10000000', 'bulk_abs_imag_max': '19/10', 'actual_moment': 'sum binom(m,j)*(a*m+1)_(d*j)*x^(m-j)', 'all_characteristic_zero_fields': True, 'individual_root_simplicity': True, 'logical_phase_inputs': ['U9', 'U12'], 'logical_location_input': 'R6', 'parent_replay_role': 'authenticate complete numerical source boundary only', 'uses_parent_coprimality': False, 'uses_soft_order_as_simplicity': False, 'uses_HJN': False, 'small_order_simplicity': False, 'all_support_uniform_squarefreeness': False, 'all_order_SFC': False}, 'analytic_boundary': {'physical_real_open': '1/100000', 'physical_imag_open': ['18999/10000', '21001/10000'], 'common_zeta_collar': '3/1000000', 'bulk_holomorphic_real_open': '1/5000000', 'bulk_holomorphic_imag_open': '191/100', 'compact_R_constant': 100000000000000, 'large_R_constant': 1000000000}, 'assurance': {'formal_verification': False, 'external_review': False, 'independent_replay': False}}
REPO_PINS = {'actual_uniform_phase': {'root': 'repo', 'path': 'research/factorial/image-rigidity/UNIFORM_ACTUAL_PHASE_AND_CONTOUR_ESTIMATES.md', 'sha256': 'f2c08b978834b3ce9b3511b930814851cfed3719d7f232972d41febcd0e4774a'}, 'all_root_location': {'root': 'repo', 'path': 'research/factorial/image-rigidity/UNIFORM_QUADRATIC_STRONG_FACTORIAL_TAIL.md', 'sha256': '7e443fc6cd6743fc6885d4de047de0a43f5547526d40b872fa39fc2471d72d5a'}, 'numerical_boundary_checker': {'root': 'repo', 'path': 'research/factorial/image-rigidity/verify_uniform_quadratic_strong_factorial_completion.py', 'sha256': '68ddabfd88f71c9bdf1b1f91fa556d171b09d319dc5637169c5f943000884215'}, 'numerical_boundary_inputs': {'root': 'repo', 'path': 'research/factorial/image-rigidity/uniform_quadratic_strong_factorial_completion_inputs.json', 'sha256': 'c7daea01d23171451c898bb4e86f1e64170ddb0744f18db846d12f05166da2fe'}, 'numerical_boundary_ledger': {'root': 'repo', 'path': 'research/factorial/image-rigidity/uniform_quadratic_strong_factorial_completion_source_ledger.json', 'sha256': '27019c6e286efb8d79ce56068c31c73c75407caadc333e1ec377a2ad52f713a9'}, 'numerical_boundary_receipt': {'root': 'repo', 'path': 'research/factorial/image-rigidity/uniform_quadratic_strong_factorial_completion_verification.json', 'sha256': '3bae7d6483722dbe3246d6254aded998fe139fa7e5150abe71c6cea65cc4ee22'}, 'airy_ode': {'root': 'repo', 'path': 'research/factorial/image-rigidity/sources/dlmf-9-2-endpoint.html', 'sha256': 'a237add8efb2c719dc7cc74baf78b435b50beec368dbc1adc672645f792b7ac7'}}
def need(condition, message):
    if not condition:
        raise ValueError(message)
def digest(p):
    return sha256(p.read_bytes()).hexdigest()
def strict_equal(a, b):
    if type(a) is not type(b):
        return False
    if isinstance(a, dict):
        return set(a) == set(b) and all(strict_equal(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return len(a) == len(b) and all(strict_equal(x, y) for x, y in zip(a,b))
    return a == b
def read_json(p):
    return json.loads(p.read_text(encoding="utf-8"))
def source_specs():
    specs = {k:dict(v) for k,v in REPO_PINS.items()}
    for role,(name,sha) in LOCAL_PROOFS.items():
        specs[role] = dict(root="local",path=name,sha256=sha)
    specs["combined_checker"] = dict(root="local",path=CHECKER)
    specs["typed_inputs"] = dict(root="local",path=INPUT)
    return specs
def resolve_row(row, repo, local):
    need(row["root"] in ("repo","local"), "Invalid source root")
    rel = Path(row["path"])
    need(not rel.is_absolute() and ".." not in rel.parts, "Unsafe source path")
    return (repo if row["root"] == "repo" else local)/rel
def validate_sources(ledger, repo):
    need(type(ledger) is dict and type(ledger.get("schema")) is int
         and ledger["schema"] == 1, "Exact schema=1 required")
    need(strict_equal(ledger.get("assurance"),EXPECTED["assurance"]), "No assurance upgrade in source ledger")
    specs = source_specs()
    rows = ledger.get("files")
    need(type(rows) is list and len(rows) == len(specs), "All eleven roles required")
    need(all(type(v) is dict and set(v) == {"role","root","path","sha256"} for v in rows),
         "Exact source row fields")
    need({v["role"] for v in rows} == set(specs), "Exact source roles")
    for row in rows:
        spec = specs[row["role"]]
        need(row["root"] == spec["root"] and row["path"] == spec["path"],
             "Role-specific source path")
        h = row["sha256"]
        need(type(h) is str and len(h) == 64 and all(c in "0123456789abcdef" for c in h),
             "Exact SHA256 syntax")
        if "sha256" in spec:
            need(h == spec["sha256"], "Reviewed source digest: "+row["role"])
        p = resolve_row(row,repo,HERE)
        need(p.is_file() and digest(p) == h, "Missing/stale source: "+row["role"])
    need(strict_equal(read_json(HERE/INPUT),EXPECTED), "Typed analytical scope changed")
def authenticate_closure(repo):
    """Hash every descendant of the exact current parent ledger before replay."""
    start = repo/REPO_PINS["numerical_boundary_ledger"]["path"]
    stack, seen, artifacts = [start], set(), {}
    while stack:
        p = stack.pop().resolve()
        if p in seen:
            continue
        seen.add(p)
        data = read_json(p)
        regular = type(data.get("schema")) is int and data["schema"] == 1
        legacy = (p.name == "explicit_airy_contour_envelope_source_ledger.json"
                  and "schema" not in data and type(data.get("schema_version")) is int
                  and data["schema_version"] == 1)
        need(regular or legacy, "Exact descendant ledger schema")
        rows = data.get("files")
        need(type(rows) is list and len(rows)>0, "Nonempty descendant source ledger")
        keys = set()
        for row in rows:
            need(type(row) is dict and set(row) == {"root","path","sha256"},
                 "Exact descendant source row")
            key = (row["root"],row["path"])
            need(key not in keys, "Duplicate descendant source")
            keys.add(key)
            q = resolve_row(row,repo,p.parent)
            need(q.is_file() and digest(q) == row["sha256"],
                 "Stale descendant: "+row["path"])
            rel = q.resolve().relative_to(repo).as_posix()
            artifacts[rel] = digest(q)
            if q.name.endswith(("_source_ledger.json","_source_manifest.json")):
                stack.append(q)
    need(len(seen) == 10, "Expected current 401 plus nine retained source ledgers")
    return dict(ledger_count=len(seen),artifact_count=len(artifacts),
                files=[dict(path=k,sha256=artifacts[k]) for k in sorted(artifacts)])
def polynomial_engine(n):
    def var(i):
        e = [0]*n; e[i]=1
        return {tuple(e):Q(1)}
    def const(c):
        return {tuple([0]*n):Q(c)} if c else {}
    def add(*terms):
        out = {}
        for term in terms:
            for e,c in term.items():
                out[e] = out.get(e,Q(0))+c
        return {e:c for e,c in out.items() if c}
    def neg(p):
        return {e:-c for e,c in p.items()}
    def sub(a,b):
        return add(a,neg(b))
    def mul(*terms):
        out=const(1)
        for term in terms:
            nxt={}
            for e,c in out.items():
                for f,d in term.items():
                    k=tuple(x+y for x,y in zip(e,f))
                    nxt[k]=nxt.get(k,Q(0))+c*d
            out={e:c for e,c in nxt.items() if c}
        return out
    return var,const,add,neg,sub,mul
def universal_identities():
    X,C,add,neg,sub,mul=polynomial_engine(11)
    A,B,xi,eps,mu,Y,Yp,R,At,Bt,Rt=[X(i) for i in range(11)]
    derivative=add(mul(At,Y),mul(A,mu,Yp),mul(eps,Bt,Yp),
                   mul(eps,B,mu,xi,Y),Rt)
    # A*(P_tau - the stated expression) equals mu*eps*B*xi*P.
    lhs=add(mul(A,derivative),
            neg(mul(mu,sub(mul(A,A),mul(B,B,eps,eps,xi)),Yp)),
            mul(mu,eps,B,xi,R),neg(mul(A,At,Y)),
            neg(mul(A,eps,Bt,Yp)),neg(mul(A,Rt)))
    rhs=mul(mu,eps,B,xi,add(mul(A,Y),mul(eps,B,Yp),R))
    need(not sub(lhs,rhs),"Cleared exact endpoint derivative identity")
    need(bool(sub(sub(lhs,mul(mu,eps,B,xi,R)),rhs)),
         "Omitting the mixed remainder must give a nonzero defect")
    X,C,add,neg,sub,mul=polynomial_engine(8)
    z,delta,vp,vm,hp,hm,jp,jm=[X(i) for i in range(8)]
    pp=sub(add(mul(z,vp),mul(sub(vp,delta),hp)),C(1))
    pm=sub(add(mul(z,vm),mul(sub(vm,delta),hm)),C(1))
    rhs=sub(sub(mul(delta,z),C(1)),mul(sub(delta,vp),add(z,hp)))
    need(not sub(pp,rhs),"Cleared one-saddle stationarity")
    numerator=sub(mul(sub(C(1),mul(delta,z)),sub(hm,hp)),
                  mul(sub(vp,vm),add(z,hp),add(z,hm)))
    need(not add(numerator,mul(pp,add(z,hm)),neg(mul(pm,add(z,hp)))),
         "Cleared exact two-saddle action derivative")
    # At Jp+Jm=0, derivative of the sum is Jp*(log Jp-log Jm)'.
    need(not sub(sub(add(mul(jp,hp),mul(jm,hm)),mul(jp,sub(hp,hm))),
                 mul(add(jp,jm),hm)),"Exact cancellation derivative")
    return ["endpoint_differentiated_period_cleared",
            "single_saddle_stationarity_cleared",
            "two_saddle_action_derivative_cleared",
            "exact_period_cancellation_derivative"]
def rational_margins():
    rows=[]
    def check(name,left,right,relation="<"):
        need(left < right if relation=="<" else left == right,name)
        rows.append(dict(name=name,left=str(left),relation=relation,right=str(right)))
    rz=Q(1,125000); mu=4*10**10; eps=Q(1,200000); m=10**16
    check("endpoint real disk inclusion",Q(1,10**7)+rz,Q(1,10**5))
    check("endpoint imaginary disk inclusion",rz,Q(1,10000))
    check("co-Lipschitz image radius",Q(3,4)*rz,Q(6,10**6),"=")
    check("reserved open zeta collar",Q(3,8)*rz,Q(3,10**6),"=")
    check("A_zeta Cauchy bound",4/rz*Q(4,3),10**6)
    check("B_zeta Cauchy bound",Q(26,25)/rz*Q(4,3),10**6)
    check("epsilon power bound",200000**3,m)
    check("mu power bound",mu**3,m*m)
    check("D epsilon squared bound",Q(1,mu)+Q(1,8),Q(3,8)**2)
    check("D/m bound",Q(1,m)+Q(3,8*mu),Q(1,10**11))
    t0=(Q(13,25)*eps+Q(1,100))/(Q(49,50)-Q(1,100))
    check("zero Airy ratio",t0,Q(1,90))
    check("compact rescaled initial Y",Q(5,2)*t0,Q(1,4))
    check("large rescaled initial Y",
          (Q(13,25)*Q(3,8)+Q(1,100))/(Q(49,50)-Q(1,10**7)),Q(1,4))
    check("Gronwall rational majorant",Q(5,4)/(1-Q(5,16)),Q(20,11),"=")
    check("Gronwall norm below two",Q(20,11),2)
    check("Airy spacing disk inside common collar",Q(1,4*mu),Q(3,10**6))
    check("large xi spacing disk remains above four",4,5-Q(1,4))
    # Equality here is an upper constant; D<5/2 makes the analytic bound strict.
    check("compact Cauchy upper constant",Q(8*10**14,m)*Q(5,2),Q(1,5),"=")
    check("large Cauchy upper constant",Q(8*10**9,10**11),Q(1,5))
    lead=Q(49,50)-Q(13,25)**2/(8*Q(49,50))
    check("endpoint exact leading coefficient",lead,Q(4633,4900),"=")
    check("endpoint leading lower bound",Q(47,50),lead)
    check("endpoint amplitude derivative errors",
          Q(10**6,90*mu)+Q(10**6,m),Q(1,10**6))
    check("compact mixed remainder",eps*Q(26,49)*5*Q(1,100)*Q(91,90),Q(1,400))
    check("large mixed remainder",Q(26,49)*Q(10**9,8*mu)*Q(91,90),Q(1,400))
    check("endpoint positive derivative reserve",Q(7,10),
          Q(47,50)-Q(1,5)-Q(1,400)-Q(1,10**6))
    r=Q(1,10**6);delta=Q(1,10**10);xr=Q(2,10**7);b=Q(1,400000);rv=1+r
    check("Hermite central height gap",Q(1,600),Q(9,200)**2)
    check("midpoint denominator perturbation",xr+2*r+delta/(2*(1-delta)),b)
    check("midpoint U perturbation",2*r+delta,Q(3,10**6))
    check("midpoint log perturbation",Q(3,10**6)+b/(Q(1,10)-b),Q(1,25000))
    peakden=xr+r+delta*rv*rv/(2*(1-delta*rv))
    check("peak denominator perturbation",peakden,Q(13,10**7))
    check("peak logarithm perturbation",Q(13,10**7)/(1-Q(13,10**7)),Q(2,10**6))
    peakU=r*(rv+1)/2+delta*rv**3/(6*(1-delta*rv/2))
    check("peak U perturbation",peakU,Q(2,10**6))
    check("peak phase reserve",Q(4,10**6),Q(1,50000))
    check("both midpoint connector losses",Q(1,1000),Q(1,600)-Q(1,25000)-Q(1,50000))
    smallv=Q(1,10)+2*r
    check("small denominator U",smallv**2/2+delta*smallv**3/(6*(1-delta*smallv/2)),Q(1,100))
    check("small denominator modulus",Q(1,10)+b,Q(101,1000))
    e2upper=sum(Q(2**k,factorial(k)) for k in range(7))+Q(2**7,factorial(7))/(1-Q(1,4))
    check("elementary e squared majorant",e2upper,8)
    check("small denominator phase gap",Q(101,1000)*8,1)
    check("connector length budget",2*r,2)
    rho=Q(1,10**8)
    check("bulk Cauchy real reserve",Q(1,10**7)+rho,Q(2,10**7))
    check("bulk Cauchy imaginary reserve",Q(19,10)+rho,Q(191,100))
    check("curvature derivative bound",3/rho,3*10**8,"=")
    check("error derivative bound",Q(1,100)/rho,10**6,"=")
    check("action derivative correction",6*delta/(1-2*delta),Q(1,10**8))
    check("action positive real part",Q(1,2),Q(59,100)-Q(1,10**8))
    check("curvature logarithmic derivative error",Q(3*10**8)/Q(17,100),2*10**9)
    check("two errors logarithmic derivative",Q(2*10**6)/Q(99,100),3*10**6)
    check("bulk positive derivative",Q(m,3),Q(m,2)-2*10**9-3*10**6)
    # Complete errors: e^-x <= 2/x^2 and sqrt(m/2)<=sqrt(m).
    # After substitution, every displayed upper bound decreases for m>=m0.
    rootm=10**8
    budgets=[Q(360,rootm), Q(6*10000**2,m*rootm),
             Q(8,m*rootm),Q(4*1000**2,m*rootm)]
    for i,budget in enumerate(budgets):
        check("complete bulk packet error "+str(i+1),budget,Q(1,1000))
    check("sum complete relative errors",sum(budgets),Q(1,100))
    check("high ratio location delta squared",Q(1,10**24),Q(1,10**12)**2,"=")
    check("bulk endpoint root-window overlap",Q(19,10),Q(21,10))
    return rows
def exact_moment_controls():
    controls=0
    wrong_shape=0
    for a in range(1,5):
        for d in range(2,7):
            for m in range(1,7):
                for j in range(m+1):
                    # Expansion of the actual positive-cycle integral coefficient.
                    integral=Q(comb(m,j)*factorial(a*m+d*j),m**(a*m+d*j+1))
                    rising=1
                    for h in range(d*j):
                        rising*=a*m+1+h
                    # Gamma(am+1)/m^(Dm+1) times A_m(m^d w).
                    polynomial=Q(factorial(a*m)*comb(m,j)*rising*m**(d*(m-j)),
                                 m**((a+d)*m+1))
                    need(integral==polynomial,"Literal Gamma/moment normalization")
                    controls+=1
                    if j:
                        wrong=1
                        for h in range(d*j):wrong*=a*m+h
                        wrong_shape+=int(wrong!=rising)
    need(wrong_shape>0,"Wrong Gamma shape control")
    return controls,wrong_shape
def mathematical_negatives():
    # Honest abstract controls, not alleged counterexamples in the analytic domain.
    need((-1)**2-2*(-1)+1==4 and 1-2+1==0 and 2*1-2==0,
         "A coprime pair can contain a repeated component")
    # p_t=x^2-t preserves the degree/root count through t=0, yet gains a double root.
    need(0*0==0 and 2*0==0 and Q(1,2)**2-Q(1,4)==0,
         "Root-count homotopy is not a simplicity theorem")
    return ["omitted_endpoint_mixed_remainder_nonzero_polynomial_defect",
            "wrong_Gamma_shape_fails_literal_coefficients",
            "coprime_pair_does_not_force_component_simplicity",
            "constant_root_count_does_not_force_simplicity"]
def fresh_boundary(repo,python):
    script=repo/REPO_PINS["numerical_boundary_checker"]["path"]
    receipt=repo/REPO_PINS["numerical_boundary_receipt"]["path"]
    with tempfile.TemporaryDirectory(prefix="high_ratio_boundary_") as td:
        out=Path(td)/"parent.json"
        run=subprocess.run([str(python),"-B",str(script),"--repo-root",str(repo),
                            "--output",str(out)],capture_output=True,text=True)
        need(run.returncode==0,"Fresh numerical-boundary replay failed: "+run.stderr[-2500:])
        need(out.read_bytes()==receipt.read_bytes(),"Fresh numerical boundary must be byte-identical")
        data=read_json(out)
        need(data["status"]=="PASS_SOURCE_GUARDED_UNIFORM_QUADRATIC_SUPPORT"
             and data["fresh_replays"]==3 and data["mandatory_local_pins"]==8
             and data["mandatory_repo_pins"]==5 and data["checked_source_ledger_closure_count"]==9,
             "Exact current parent receipt structure")
        need(all(data.get(k) is False for k in
                 ("formal_verification","external_review","independent_replay")),
             "No assurance upgrade")
        return dict(receipt_sha256=digest(out),fresh_parent_replays=1,
                    fresh_parent_children=3,role="numerical source authentication only")
def rejection_controls(ledger,repo):
    count=0
    def reject(v):
        nonlocal count
        try:validate_sources(v,repo)
        except (ValueError,KeyError,TypeError,FileNotFoundError):
            count+=1
            return
        raise ValueError("Malformed ledger accepted")
    for i in range(len(ledger["files"])):
        bad=copy.deepcopy(ledger);del bad["files"][i];reject(bad)
        bad=copy.deepcopy(ledger);bad["files"][i]["sha256"]="0"*64;reject(bad)
    for v in [True,0,2,"1"]:
        bad=copy.deepcopy(ledger);bad["schema"]=v;reject(bad)
    bad=copy.deepcopy(ledger);bad["files"].append(copy.deepcopy(bad["files"][0]));reject(bad)
    for k,v in [("uses_parent_coprimality",True),("uses_soft_order_as_simplicity",True),
                ("uses_HJN",True),("small_order_simplicity",True),
                ("all_support_uniform_squarefreeness",True),
                ("integer_m_min",True),("all_root_location_delta_max","1/10000000000")]:
        bad=copy.deepcopy(EXPECTED);bad["scope"][k]=v
        need(not strict_equal(bad,EXPECTED),"False scope rejected");count+=1
    return count
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--repo-root",required=True,type=Path)
    ap.add_argument("--parent-python",type=Path)
    ap.add_argument("--output",type=Path,default=HERE/RECEIPT)
    args=ap.parse_args()
    repo=args.repo_root.resolve()
    candidates=[repo/".venv/bin/python",repo/".venv/Scripts/python.exe",Path(sys.executable)]
    python=args.parent_python.absolute() if args.parent_python else next(p for p in candidates if p.is_file())
    need(python.is_file(),"Parent interpreter is a file")
    ledger=read_json(HERE/LEDGER)
    validate_sources(ledger,repo)
    closure=authenticate_closure(repo)
    identities=universal_identities()
    margins=rational_margins()
    moments,wrongshape=exact_moment_controls()
    negatives=mathematical_negatives()
    controls=rejection_controls(ledger,repo)
    parent=fresh_boundary(repo,python)
    result=dict(schema=1,status="PASS_SOURCE_GUARDED_ACTUAL_HIGH_RATIO_SIMPLICITY",
      scope=EXPECTED["scope"],analytical_boundary=EXPECTED["analytic_boundary"],
      mathematical_premises={"actual_phase":["U9","U12"],"all_root_location":"R6",
       "numerical_replay_is_not_a_coprimality_premise":True,
       "simple_roots_from_derivatives_not_homotopy":True},
      finite_check_scope="Exact identities, rational budgets, literal normalizations and source freshness only; the retained written proofs supply analytic universal estimates.",
      universal_identities=identities,rational_guards=margins,
      exact_moment_coefficient_controls=moments,wrong_Gamma_shape_controls=wrongshape,
      mathematical_negative_controls=negatives,source_scope_controls=controls,
      mandatory_roles=len(ledger["files"]),pins=ledger["files"],
      authenticated_numerical_closure=closure,parent_replay=parent,
      source_ledger_sha256=digest(HERE/LEDGER),
      formal_verification=False,external_review=False,independent_replay=False)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(dict(status="PASS",roles=len(ledger["files"]),closure_ledgers=closure["ledger_count"],
      closure_artifacts=closure["artifact_count"],identities=len(identities),
      rational_guards=len(margins),moment_coefficients=moments,negatives=len(negatives),
      source_controls=controls,fresh_boundary_replays=1)))
if __name__=="__main__":
    main()
