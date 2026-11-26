#!/usr/bin/env python3
"""Exact support controls for the fixed-integer-cusp compact-complement proof.

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
PROOF = "CUSP_GLOBAL_REGULAR_COMPLEMENT_AND_STOKES_DETOUR.md"
CHECKER = "verify_cusp_global_complement.py"
INPUT = "cusp_global_complement_inputs.json"
LEDGER = "cusp_global_complement_source_ledger.json"
EXPECTED = {
    "schema": 1,
    "fixed_integer_cusps": True,
    "compact_excludes_positive_cusp": True,
    "complete_stokes_detour": True,
    "nonnegative_boundary_includes_zero": True,
    "negative_boundary_pair_exclusion_only": True,
    "retained_upper_bound_only": True,
    "effective_cutoff": False,
    "uniform_in_s": False,
    "all_order_coprimality": False,
    "all_root_localization": False,
}
SOURCES = {
    "proof": ("local", PROOF),
    "checker": ("local", CHECKER),
    "inputs": ("local", INPUT),
    "two_thimble_proof": ("support", "CUSP_TWO_THIMBLE_POSITIVE_CYCLE.md"),
    "distinguished_sector_proof": (
        "repo", "research/factorial/image-rigidity/CUSP_STATIONARY_BELYI_MONODROMY_AND_RADIAL_BARRIER.md"),
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
    # Universal cusp identities after a=s^2, d=2s+1.
    s,S = X(0),X(1)
    a = power(s,2)
    d = add(scale(s,2),C(1))
    D = add(a,d)
    sc = mul(s,add(s,C(1)))
    need(not sub(add(scale(a,2),d,C(-1)),scale(sc,2)), "Cusp double-root center")
    need(not sub(mul(a,D),power(sc,2)), "Cusp double-root discriminant")
    # Denominator-cleared Omega' identity, no fixed value of s is used.
    left = sub(mul(sub(mul(d,D),mul(add(d,C(1)),S)),sub(S,a)),mul(S,sub(D,S)))
    need(not add(left,mul(d,power(sub(S,sc),2))), "Exact Omega derivative factor")
    # Universal envelope identity with T denoting S^d independently.
    ss,aa,dd,w,T = X(0),X(1),X(2),X(3),X(4)
    DD = add(aa,dd)
    env = sub(mul(dd,w),mul(sub(DD,ss),add(w,T)))
    stationarity = sub(mul(sub(ss,aa),add(w,T)),mul(dd,T))
    need(not sub(env,stationarity), "Envelope modulo exact stationarity")
    # Radial chain rule after multiplying by m; Dw*d*w=-Delta.
    mm,H,Delta = X(0),X(1),X(2)
    chain = add(scale(mul(mm,H),-1),scale(mul(mm,Delta),-1))
    need(not sub(chain,scale(mul(mm,add(H,Delta)),-1)), "Exact radial action drop")
    identities = 5

    normalizations = 0
    upper_factor_controls = 0
    for sv in range(1,6):
        av,dv = sv*sv,2*sv+1
        for mv in range(1,7):
            cc = coeffs(sv,mv)
            for k,ck in enumerate(cc):
                literal = Q(comb(mv,k)*factorial(av*mv+dv*k),mv**(av*mv+dv*k+1))
                via_A = Q(factorial(av*mv)*ck*mv**(dv*(mv-k)),mv**((av+dv)*mv+1))
                need(literal == via_A, "Actual integral/moment normalization")
                normalizations += 1
            gamma = av*mv+1
            for j in range(dv):
                need(gamma+dv*(mv-1)+j <= (av+dv)*mv, "Retained upper-factor bound")
                upper_factor_controls += 1

    # Genuine algebra controls against tempting false strengthened statements.
    wrong_sc = add(sc,C(1))
    need(add(left,mul(d,power(sub(S,wrong_sc),2))), "Corrupt cusp center was not detected")
    # Using the upper N scale for the lower order changes the actual polynomial,
    # and the vectors are not proportional (constant terms already agree).
    cv = coeffs(1,2)
    correct = [cv[k]*2**(3*(2-k)) for k in range(3)]
    wrong = [cv[k]*3**(3*(2-k)) for k in range(3)]
    need(correct[-1] == wrong[-1] and correct != wrong, "False common own-order scale not detected")
    # Every odd-degree actual A_m has a negative real root; do not claim global
    # eventual zero-freeness on all negative-real neighborhoods.
    odd_controls = 0
    for sv in range(1,4):
        for mv in (1,3,5):
            cc = coeffs(sv,mv)
            R = 1+sum(cc)
            at_zero = cc[-1]
            at_minus_R = sum(ck*(-R)**(mv-k) for k,ck in enumerate(cc))
            need(at_zero > 0 and at_minus_R < 0, "Actual odd-m negative-axis sign control")
            odd_controls += 1
    return {"universal_identities":identities,
            "actual_moment_coefficient_controls":normalizations,
            "retained_upper_factor_controls":upper_factor_controls,
            "negative_axis_sign_controls":odd_controls,
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
                    help="Directory holding the pinned two-thimble proof; defaults to this package directory")
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
    for key in ("all_order_coprimality","uniform_in_s","all_root_localization","effective_cutoff"):
        q=deepcopy(data); q[key]=True
        need(rejected(validate_inputs,q), "False stronger scope accepted")
        checks += 1
    q=deepcopy(data); q["fixed_integer_cusps"]=1
    need(rejected(validate_inputs,q), "Boolean/integer scope confusion")
    checks += 1
    tests["source_and_scope_negative_controls"]=checks
    out = {
        "schema":1,"status":"PASS",
        "scope":"Written analytic compact-complement proof support for each fixed integer cusp; no effective cutoff or all-order conclusion",
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

