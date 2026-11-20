#!/usr/bin/env python3
"""Exact support controls for the above-critical regular compact complement.

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
PROOF = "ABOVE_CRITICAL_REGULAR_COMPACT_COMPLEMENT.md"
CHECKER = "verify_above_critical_regular_complement.py"
INPUT = "above_critical_regular_inputs.json"
LEDGER = "above_critical_regular_source_ledger.json"
EXPECTED = {
    "schema":1,
    "integer_a_min":1,
    "integer_d_min":2,
    "strict_abovecritical_support":True,
    "fixed_support":True,
    "compact_avoids_both_turning_values":True,
    "a_one_d_two_included":True,
    "shrunk_wedge_for_d_two":True,
    "equal_height_stokes_complete_packets":True,
    "novertical_input_only_at_equal_moduli":True,
    "global_fixed_branch_divider":False,
    "turning_neighborhood_excluded_by_this_module":False,
    "global_eventual_sfc_claimed_by_this_module":False,
    "all_order_coprimality":False,
    "effective_cutoff":False,
    "uniform_support_cutoff":False,
}
SOURCES = {
    "proof":("local",PROOF),
    "checker":("local",CHECKER),
    "inputs":("local",INPUT),
    "novertical_proof":("support","NO_VERTICAL_PHYSICAL_EQUIMODULAR_PAIR.md"),
    "coalescing_actual_bridge":("repo","research/factorial/image-rigidity/COALESCING_BINOMIAL_SADDLE_BRIDGE.md"),
    "morse_mechanisms":("repo","research/factorial/image-rigidity/CUSP_TWO_THIMBLE_POSITIVE_CYCLE.md"),
    "finite_level_stokes_detour":("repo","research/factorial/image-rigidity/CUSP_GLOBAL_REGULAR_COMPLEMENT_AND_STOKES_DETOUR.md"),
    "retained_gamma_upper_bound":("repo","research/factorial/image-rigidity/GAMMA_PRODUCT_ROOT_MODULUS.tex"),
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
    need(type(doc["sources"]) is list and len(doc["sources"]) == len(SOURCES), "All eight sources required")
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


def core():
    # Full formal parameter identities; the six-variable ring is over Q.
    S,b,tau = X(0),X(1),X(2)
    a=add(power(b,2),power(tau,2));d=add(scale(b,2),C(1))
    D=add(a,d);A=add(a,b);c=mul(b,add(b,C(1)))
    Np=sub(mul(d,D),mul(add(d,C(1)),S))
    left=sub(mul(Np,sub(S,a)),mul(S,sub(D,S)))
    need(not add(left,mul(d,add(power(sub(S,A),2),power(tau,2)))),"Abovecritical Omega derivative")
    need(not sub(add(power(A,2),power(tau,2)),mul(a,D)),"Turning |S| square")
    need(not sub(add(power(b,2),power(tau,2)),a),"Turning |S-a| square")
    need(not sub(add(power(add(b,C(1)),2),power(tau,2)),D),"Turning |D-S| square")
    rhs=power(sub(a,scale(mul(d,sub(d,C(1))),Q(1,2))),2)
    need(not sub(sub(power(A,2),mul(sub(power(d,2),C(1)),power(tau,2))),rhs),"Physical-sector angle bound")
    # Exact upper critical-value argument derivative, after multiplying aD.
    phase_num=add(scale(mul(d,sub(c,power(tau,2))),-1),mul(b,D),mul(add(b,C(1)),a))
    need(not sub(phase_num,scale(mul(d,power(tau,2)),2)),"Critical argument derivative")
    # Imag[z/(w+z)] times |w+z|^2 = Im[z conjugate(w)].
    ux,uy,wx,wy=X(0),X(1),X(2),X(3)
    imag_product=sub(mul(uy,add(ux,wx)),mul(ux,add(uy,wy)))
    need(not sub(imag_product,sub(mul(uy,wx),mul(ux,wy))),"Exact angular boundary sign")
    # Stationary envelope, with T denoting S^d.
    S,a,d,w,T=X(0),X(1),X(2),X(3),X(4)
    D=add(a,d)
    need(not sub(sub(mul(d,w),mul(sub(D,S),add(w,T))),
                 sub(mul(sub(S,a),add(w,T)),mul(d,T))),"Exact stationary envelope")
    identities=8
    fixtures=[(1,2),(2,2),(2,3),(3,2),(4,2),(5,3),(7,4),(10,5),(20,7)]
    norms=0;upper=0;slopes=0
    for av,dv in fixtures:
        need(av>=1 and dv>=2 and 4*av>(dv-1)**2,"Actual abovecritical fixture")
        A=Q(2*av+dv-1,2);t2=Q(4*av-(dv-1)**2,4);D=av+dv
        for sv in (Q(1,3),Q(av,2),Q(2*av+dv,2),Q(D+1)):
            if sv==av:continue
            num=-dv*sv**(dv-1)*((sv-A)**2+t2)
            need(num<0 and (sv-av)**2>0,"Strict positive-real stationary monotonicity")
            slopes+=1
        for mv in range(1,6):
            for k in range(mv+1):
                ck=comb(mv,k)*rising(av*mv+1,dv*k)
                literal=Q(comb(mv,k)*factorial(av*mv+dv*k),mv**(av*mv+dv*k+1))
                via=Q(factorial(av*mv)*ck*mv**(dv*(mv-k)),mv**((av+dv)*mv+1))
                need(literal==via,"Actual own-order factorial identity")
                norms+=1
            for j in range(dv):
                need(av*mv+1+dv*(mv-1)+j <= (av+dv)*mv,"Old global upper factor")
                upper+=1
    # Exact cap and equal-height exit budgets; these support the explicit
    # general triangle inequalities, not a numerical global cutoff.
    budgets=0
    for eta in (Q(1,1000),Q(1,20),Q(1),Q(7,3)):
        for diff in (-eta/4,Q(0),eta/4):
            heights=(Q(2),Q(2)+diff)
            need(max(h-eta for h in heights)<min(heights),"Both local exits lie below both heights")
            budgets+=1
    for cosa in (Q(1,10000),Q(1,7),Q(1,2)):
        need(cosa-cosa/2>0,"Shrunk infinity radial stability")
        budgets+=1
    # Three actual countercontrols to false strengthened statements.
    # A full d2 lower boundary S=-ir has no exponential loss:
    # at a1,w=i, |psi|^2=r^2*(r^4+1), increasing to infinity.
    need(all((r+1)**2*((r+1)**4+1)>r*r*(r**4+1) for r in range(1,9)),"False unshrunk d2 decay not detected")
    def cmul(z,v):return(z[0]*v[0]-z[1]*v[1],z[0]*v[1]+z[1]*v[0])
    def cdiv(z,v):
        den=v[0]*v[0]+v[1]*v[1]
        return((z[0]*v[0]+z[1]*v[1])/den,(z[1]*v[0]-z[0]*v[1])/den)
    # The proposed global vertical divider ReS=a+(d-1)/2 FAILS abovecritical:
    # a2,d2,S=5/2-i is physical and maps to upper w=139/20+47i/5.
    z=(Q(5,2),Q(-1));value=cmul(cmul(z,z),cdiv((Q(3,2),Q(1)),(Q(1,2),Q(-1))))
    need(value==(Q(139,20),Q(47,5)) and value[1]>0,"False fixed global divider not rejected")
    av,dv,mv,N=2,2,2,3
    cc=[comb(mv,k)*rising(av*mv+1,dv*k) for k in range(mv+1)]
    good=[cc[k]*mv**(dv*(mv-k)) for k in range(mv+1)]
    bad=[cc[k]*N**(dv*(mv-k)) for k in range(mv+1)]
    need(good!=bad and good[-1]==bad[-1],"False shared own-order scale")
    return {"universal_identity_families":identities,
      "actual_abovecritical_support_fixtures":len(fixtures),
      "actual_own_order_coefficient_controls":norms,
      "old_global_upper_factor_controls":upper,
      "strict_real_stationary_slope_controls":slopes,
      "exit_and_infinity_budget_controls":budgets,
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
                    help="Directory holding the pinned no-vertical proof; defaults to this package directory")
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
    for i in range(len(SOURCES)):
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
    for key in ("all_order_coprimality","uniform_support_cutoff","global_eventual_sfc_claimed_by_this_module","turning_neighborhood_excluded_by_this_module","global_fixed_branch_divider","effective_cutoff"):
        q=deepcopy(data); q[key]=True
        need(rejected(validate_inputs,q), "False stronger scope accepted")
        checks += 1
    q=deepcopy(data); q["integer_a_min"]=True
    need(rejected(validate_inputs,q), "Boolean/integer scope confusion")
    checks += 1
    tests["source_and_scope_negative_controls"]=checks
    out = {
        "schema":1,"status":"PASS",
        "scope":"Written actual abovecritical regular compact-complement support; no turning-neighborhood, global SFC, effective or uniform-cutoff conclusion",
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

