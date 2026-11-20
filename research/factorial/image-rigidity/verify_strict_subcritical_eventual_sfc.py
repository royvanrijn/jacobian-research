#!/usr/bin/env python3
"""Exact support controls for the strict-subcritical eventual proof.

The analytic finite-level contour and Gaussian arguments are written proofs,
not consequences of the finite examples checked here. No effective cutoff,
all-order theorem, or assurance beyond this scoped internal replay is claimed.
"""
import argparse
from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import comb, factorial
from pathlib import Path
import sys

if not __debug__:
    raise SystemExit("Refusing optimized Python: proof-support assertions are required.")

HERE = Path(__file__).resolve().parent
PROOF = "SUBCRITICAL_ORDERED_PAIR_EVENTUAL_SFC.md"
CHECKER = "verify_strict_subcritical_eventual_sfc.py"
INPUT = "strict_subcritical_inputs.json"
LEDGER = "strict_subcritical_source_ledger.json"
EXPECTED = {
    "schema": 1,
    "fixed_ordinary_integer_supports": True,
    "strict_subcritical_only": True,
    "globally_ordered_physical_pair": True,
    "stokes_argument_revalidated": True,
    "spinodal_inflections_strictly_inactive": True,
    "coexistence_two_complete_packets": True,
    "eventual_coprimality": True,
    "effective_cutoff": False,
    "uniform_terminal_bound": False,
    "all_order_coprimality": False,
    "critical_equality_proved_by_this_note": False,
}
SOURCES = {
    "proof": ("local", PROOF),
    "checker": ("local", CHECKER),
    "inputs": ("local", INPUT),
    "two_thimble_proof": ("support", "CUSP_TWO_THIMBLE_POSITIVE_CYCLE.md"),
    "global_complement_proof": ("support", "CUSP_GLOBAL_REGULAR_COMPLEMENT_AND_STOKES_DETOUR.md"),
    "retained_gamma_upper_bound": (
        "repo", "research/factorial/image-rigidity/GAMMA_PRODUCT_ROOT_MODULUS.tex"),
}

def need(cond, msg):
    if not cond:
        raise ValueError(msg)

def digest(p):
    return sha256(p.read_bytes()).hexdigest()

def validate_inputs(doc):
    need(type(doc) is dict and set(doc) == set(EXPECTED), "Exact input schema/keys required")
    for k, v in EXPECTED.items():
        need(type(doc[k]) is type(v) and doc[k] == v, "Input type or scope mismatch: "+k)

def validate_ledger(doc, repo, support):
    need(type(doc) is dict and set(doc) == {"schema", "sources"}, "Exact ledger schema required")
    need(type(doc["schema"]) is int and doc["schema"] == 1, "Ledger schema must be integer1")
    need(type(doc["sources"]) is list and len(doc["sources"]) == len(SOURCES), "All six sources required")
    seen = set()
    for rec in doc["sources"]:
        need(type(rec) is dict and set(rec) == {"role", "root", "path", "sha256"}, "Malformed source")
        need(all(type(rec[k]) is str for k in rec), "Source fields must be strings")
        role = rec["role"]
        need(role in SOURCES and role not in seen, "Missing, duplicate, or unknown source role")
        seen.add(role)
        need((rec["root"], rec["path"]) == SOURCES[role], "Exact source selection required")
        h = rec["sha256"]
        need(len(h) == 64 and all(c in "0123456789abcdef" for c in h), "Malformed SHA256")
        root = {"local": HERE, "support": support, "repo": repo}[rec["root"]]
        p = root / rec["path"]
        need(p.is_file() and digest(p) == h, "Source pin mismatch: "+role)
    need(seen == set(SOURCES), "Incomplete source roles")

# Tiny sparse polynomial ring in six indeterminates over Q.
ZERO = (0,)*6
def C(c):
    return {} if not c else {ZERO: Q(c)}
def X(i):
    e = [0]*6
    e[i] = 1
    return {tuple(e): Q(1)}
def add(*ps):
    out = {}
    for p in ps:
        for e,c in p.items():
            out[e] = out.get(e,Q(0))+c
            if not out[e]:
                del out[e]
    return out
def scale(p,c):
    return {e:v*c for e,v in p.items() if v*c}
def mul(p,q):
    out = {}
    for e,c in p.items():
        for f,d in q.items():
            g = tuple(a+b for a,b in zip(e,f))
            out[g] = out.get(g,Q(0))+c*d
            if not out[g]:
                del out[g]
    return out
def power(p,n):
    out = C(1)
    for _ in range(n):
        out = mul(out,p)
    return out
def sub(p,q):
    return add(p,scale(q,-1))

def rising(x,n):
    out = 1
    for j in range(n):
        out *= x+j
    return out

def coeffs(s,m):
    a,d = s*s,2*s+1
    return [comb(m,k)*rising(a*m+1,d*k) for k in range(m+1)]

def core():
    # Denominator-cleared derivative of the physical stationary map, valid
    # symbolically without fixing the integer support.
    S,a,d = X(0),X(1),X(2)
    D = add(a,d)
    b = scale(add(d,C(-1)),Q(1,2))
    A = add(a,b)
    q2 = sub(power(b,2),a)
    left = sub(mul(sub(mul(d,D),mul(add(d,C(1)),S)),sub(S,a)),mul(S,sub(D,S)))
    expected = scale(mul(d,sub(power(sub(S,A),2),q2)),-1)
    need(not sub(left,expected), "General Omega derivative")
    identities=1
    # Exact stationary roots and value-factor identities in b,q.
    b,q = X(0),X(1)
    a = sub(power(b,2),power(q,2)); d=add(scale(b,2),C(1))
    D=add(a,d); A=add(a,b)
    for sign in (-1,1):
        P=add(A,scale(q,sign))
        x=add(b,scale(q,sign)); y=add(b,C(1),scale(q,-sign))
        need(not sub(P,mul(x,y)), "Spinodal product factorization"); identities+=1
        need(not sub(sub(D,P),y), "Spinodal numerator factor"); identities+=1
        need(not sub(sub(P,a),x), "Spinodal denominator factor"); identities+=1
    # Log(w_plus/w_minus)'=4*d*q^2/[a*(a+d)].
    der_num = sub(scale(mul(power(b,2),D),4),scale(mul(power(add(b,C(1)),2),a),4))
    need(not sub(der_num,scale(mul(d,power(q,2)),4)), "Critical-value ordering derivative")
    identities+=1
    # Full rational derivative identity for the vertical-divider angle.
    b,Y = X(0),X(1)
    d=add(scale(b,2),C(1)); sc=mul(b,add(b,C(1)))
    E0=add(C(1),power(Y,2))
    E1=add(C(1),mul(power(b,2),power(Y,2)))
    E2=add(C(1),mul(power(add(b,C(1)),2),power(Y,2)))
    numerator=add(mul(d,mul(E1,E2)),scale(mul(b,mul(E0,E2)),-1),
                  scale(mul(add(b,C(1)),mul(E0,E1)),-1))
    rhs=mul(mul(d,sc),mul(power(Y,2),add(C(1),mul(sub(sc,C(1)),power(Y,2)))))
    need(not sub(numerator,rhs), "Divider angle derivative"); identities+=1
    # The stationary envelope, retaining S^d as an independent formal T.
    S,a,d,w,T = X(0),X(1),X(2),X(3),X(4)
    D=add(a,d)
    need(not sub(sub(mul(d,w),mul(sub(D,S),add(w,T))),
                 sub(mul(sub(S,a),add(w,T)),mul(d,T))), "Envelope"); identities+=1
    # Exact square-root-gap algebra, with x=sqrt(a), h=the positive gap.
    x,h=X(0),X(1)
    d=sub(power(add(x,h),2),power(x,2))
    need(not sub(add(d,C(-1),scale(x,-2)),mul(sub(h,C(1)),add(scale(x,2),h,C(1)))),
         "Reversible square-root-gap factor"); identities+=1

    pairs=[(1,4),(1,5),(2,4),(2,5),(3,5),(4,6),(6,6),(9,8)]
    coefficient_controls=0
    factor_controls=0
    for av,dv in pairs:
        need(4*av<(dv-1)**2 and av>=1 and dv>=3, "Actual strict-support fixture")
        for mv in range(1,6):
            for k in range(mv+1):
                ck=comb(mv,k)*rising(av*mv+1,dv*k)
                literal=Q(comb(mv,k)*factorial(av*mv+dv*k),mv**(av*mv+dv*k+1))
                via=Q(factorial(av*mv)*ck*mv**(dv*(mv-k)),mv**((av+dv)*mv+1))
                need(literal==via,"Actual own-order normalization")
                coefficient_controls+=1
            for j in range(dv):
                need(av*mv+1+dv*(mv-1)+j<=(av+dv)*mv,"Old upper-root factor")
                factor_controls+=1

    # Wrongly replacing the strict derivative by the cusp square loses q^2.
    sv,av,dv=Q(2),Q(1),4
    Nv=sv**dv*(av+dv-sv)
    Np=dv*sv**(dv-1)*(av+dv-sv)-sv**dv
    actual=(Np*(sv-av)-Nv)/(sv-av)**2
    wrong=-dv*sv**(dv-1)*(sv-(av+Q(dv-1,2)))**2/(sv-av)**2
    need(actual==32 and wrong==-8 and actual!=wrong, "Corrupt cusp-square derivative")
    # A fixed nonclosed field cannot be used for the root-existence converse:
    # p=x^2+1 has no rational root although gcd_Q(p,p)=p.
    # A primitive integral equation n^2+q^2=0 has no solution:
    # modulo3 both residues must vanish, contradicting primitivity.
    need(all((n*n+q*q)%3 != 0 for n in range(3) for q in range(3)
             if (n,q)!=(0,0)), "Exact primitive mod3 nonclosed-field control")
    # Wrong lower own-order scaling changes a genuine ordinary polynomial.
    av,dv,mv,N=1,4,2,3
    cc=[comb(mv,k)*rising(av*mv+1,dv*k) for k in range(mv+1)]
    good=[cc[k]*mv**(dv*(mv-k)) for k in range(mv+1)]
    bad=[cc[k]*N**(dv*(mv-k)) for k in range(mv+1)]
    need(good!=bad and good[-1]==bad[-1],"Wrong lower own-order scale")
    # Explicit equality instances belong to the independent cusp theorem.
    equality_controls=0
    for sv in range(1,9):
        av,dv=sv*sv,2*sv+1
        need(4*av==(dv-1)**2 and (av+dv)==(sv+1)**2,"Critical equality")
        equality_controls+=1
    return {"universal_identity_families":identities,
            "strict_integer_support_fixtures":len(pairs),
            "actual_moment_coefficient_controls":coefficient_controls,
            "old_upper_factor_controls":factor_controls,
            "independent_critical_equality_controls":equality_controls,
            "genuine_false_statement_controls":3}

def rejected(fn,doc):
    try:
        fn(doc)
    except (ValueError,KeyError,TypeError,OSError):
        return True
    return False

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo-root",type=Path,required=True)
    ap.add_argument("--support-root",type=Path,default=HERE,
                    help="Directory holding both pinned Morse/complement proofs; defaults to this package directory")
    ap.add_argument("--input",type=Path,default=HERE/INPUT)
    ap.add_argument("--ledger",type=Path,default=HERE/LEDGER)
    ap.add_argument("--output",type=Path)
    args = ap.parse_args()
    repo,support = args.repo_root.resolve(),args.support_root.resolve()
    data = json.loads(args.input.read_text())
    ledger = json.loads(args.ledger.read_text())
    validate_inputs(data)
    validate_ledger(ledger,repo,support)
    tests = core()
    checks = 0
    for i in range(6):
        q = deepcopy(ledger); del q["sources"][i]
        need(rejected(lambda j:validate_ledger(j,repo,support),q), "Missing-source negative failed")
        checks += 1
        q = deepcopy(ledger); q["sources"][i]["sha256"] = "0"*64
        need(rejected(lambda j:validate_ledger(j,repo,support),q), "Bad-source negative failed")
        checks += 1
    for replacement in (True,0):
        q=deepcopy(ledger); q["schema"]=replacement
        need(rejected(lambda j:validate_ledger(j,repo,support),q), "Ledger schema type failed")
        checks += 1
    q=deepcopy(ledger); q["sources"][1]=deepcopy(q["sources"][0])
    need(rejected(lambda j:validate_ledger(j,repo,support),q), "Duplicate-source failed")
    checks += 1
    for key in ("all_order_coprimality","uniform_terminal_bound","critical_equality_proved_by_this_note","effective_cutoff"):
        q=deepcopy(data); q[key]=True
        need(rejected(validate_inputs,q), "False stronger scope accepted")
        checks += 1
    q=deepcopy(data); q["fixed_ordinary_integer_supports"]=1
    need(rejected(validate_inputs,q), "Boolean/integer scope confusion")
    checks += 1
    tests["source_and_scope_negative_controls"]=checks
    out = {
        "schema":1,"status":"PASS",
        "scope":"Written strict-subcritical eventual proof support for fixed ordinary integer supports; no effective cutoff, uniform bound or all-order conclusion",
        "proof_sha256":digest(HERE/PROOF),
        "checker_sha256":digest(Path(__file__).resolve()),
        "inputs_sha256":digest(args.input),
        "source_ledger_sha256":digest(args.ledger),
        "mandatory_source_records":len(SOURCES),
        "checks":tests,
        "assurance":{"formal_verification":False,"independent_replay":False,"external_review":False},
    }
    text=json.dumps(out,indent=2,sort_keys=True)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(text)
    else:
        print(text,end="")

if __name__=="__main__":
    try:
        main()
    except Exception as e:
        print("FAIL: "+str(e),file=sys.stderr)
        raise SystemExit(1)

