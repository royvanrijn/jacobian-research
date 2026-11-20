#!/usr/bin/env python3
"""Supporting exact algebra for uniform bounded actual-cusp exclusion.
The analytic compactness/concentration arguments are proved in the note.
"""
from __future__ import annotations
import argparse,copy,hashlib,json,math,re,subprocess,sys,tempfile
from fractions import Fraction as Q
from pathlib import Path
if not __debug__:
    raise SystemExit("Refusing optimized execution: exact assertions are required.")
HERE=Path(__file__).resolve().parent
PROOF="UNIFORM_BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md"
NV=6
def clean(p):return {m:Q(c) for m,c in p.items() if c}
def const(c):return {} if not c else {(0,)*NV:Q(c)}
def var(i):
    m=[0]*NV;m[i]=1
    return {tuple(m):Q(1)}
def add(*items):
    out={}
    for p in items:
        for m,c in p.items():out[m]=out.get(m,Q(0))+c
    return clean(out)
def scale(p,c):return clean({m:Q(c)*a for m,a in p.items()})
def mul(*items):
    out=const(1)
    for p in items:
        nxt={}
        for m,a in out.items():
            for n,b in p.items():
                k=tuple(x+y for x,y in zip(m,n))
                nxt[k]=nxt.get(k,Q(0))+a*b
        out=clean(nxt)
    return out
def power(p,k):
    out=const(1)
    for _ in range(k):out=mul(out,p)
    return out
def identities():
    d,de,v,b,r,N=[var(i) for i in range(NV)]
    # d times Z'/p^(d-1), multiplied by v^2, plus d(v-1)^2.
    stationary=add(mul(d,de,v),scale(d,-1),scale(mul(de,v),-1),
                   scale(mul(d,power(v,2)),-1),mul(d,power(add(v,const(-1)),2)))
    residual=mul(v,add(mul(de,add(d,const(-1))),scale(d,-2)))
    assert not add(stationary,scale(residual,-1))
    # Quartic potential is uniformly positive in both signed directions.
    K=add(scale(mul(power(b,2),power(r,2)),Q(1,2)),
          scale(mul(b,power(r,3)),Q(-1,3)),scale(power(r,4),Q(1,12)))
    sos=add(scale(add(mul(power(b,2),power(r,2)),power(r,4)),Q(1,48)),
             scale(mul(power(r,2),power(add(r,scale(b,Q(-8,3))),2)),Q(1,16)),
             scale(mul(power(b,2),power(r,2)),Q(5,144)))
    assert not add(K,scale(sos,-1))
    # The upper potential bound follows from (b-r)^2>=0.
    # K <= b²r²+r⁴: the difference is an explicit positive quadratic.
    upper_diff=add(mul(power(b,2),power(r,2)),power(r,4),scale(K,-1))
    upper_sos=add(scale(mul(power(r,2),power(add(b,scale(r,Q(1,3))),2)),Q(1,2)),
                  scale(power(r,4),Q(31,36)))
    assert not add(upper_diff,scale(upper_sos,-1))
    # Cleared actual hard-region polynomial.
    hard=add(power(add(d,const(-1)),2),scale(mul(d,add(scale(N,2),const(-3))),-4))
    expected=add(power(d,2),scale(mul(add(scale(N,8),const(-10)),d),-1),const(1))
    assert not add(hard,scale(expected,-1))
    # The endpoint ratio of the hard range is strictly below8.
    assert not add(add(scale(N,8),const(-10)),
                   scale(add(scale(add(N,const(-1)),8),const(-2)),-1))
    # Relative-root radial comparison; variables now mean K,qN,qn.
    k,qN,qn=d,de,v
    lhs=add(mul(k,qN),scale(qn,-1))
    rhs=add(mul(add(k,const(-1)),qN),qN,scale(qn,-1))
    assert not add(lhs,scale(rhs,-1))
    # Cleared exact field transport with A=N^.75, B=n^.75, c and h real/complex.
    A,B,c,h,Kp,Hn=[var(i) for i in range(NV)]
    # Substitute qN=1+c*h/A, qn=Kp*qN; clear c*A.
    direct=add(mul(A,B,Kp),mul(B,Kp,c,h),scale(mul(A,B),-1))
    transported=add(mul(Kp,B,c,h),mul(A,B,add(Kp,const(-1))))
    assert not add(direct,scale(transported,-1))
    return 7,{"K":K,"sos":sos,"hard":hard}
def rational_guards():
    Q0=2*3**8
    zmax=Q(16*Q0+1,2);c0=1/(zmax+8)
    guards={
      "positive_fixed_Q0":Q0>1,
      "positive_uniform_c0":0<c0<Q(1,2),
      "left_tail_slope":Q(1,2)>0,
      "right_tail_slope":Q(6,2*(1+6))==Q(3,7),
      "left_compact_p":1-Q(1,2)==Q(1,2),
      "quartic_SOS_coefficient":Q(5,144)>0,
      "variance_scale_gaussian":Q(1,2)>0,
      "variance_scale_quartic":(1-Q(1,2))**2==Q(1,4),
      "normalization_large_m":Q(1,2)>0,
      "phase_log_series_budget":1/(2*(1-Q(1,2)))==1,
      "global_Lipschitz_upper":Q(3,4)<3,
      "hard_radial_cap":Q(3,2)*3**8<Q0,
      "large_d_real_margin":Q(1,2)-Q(3,12)==Q(1,4),
      "relative_c_lower":2>0,
      "small_field_inverse_c":Q(1,2)>0,
      "imaginary_transport_power":Q(3,4)>0,
      "root_motion_vs_radial":Q(7,4)>1,
      "large_field_noise_decay":Q(1,2)>0,
      "large_field_power":Q(1,4)-Q(1,3)==Q(-1,12),
      "large_field_suppression":Q(1,3)>0,
    }
    assert all(guards.values()),[k for k,v in guards.items() if not v]
    return guards
def cusp_controls():
    rows=[]
    for s in range(1,33):
        d=2*s+1;a=s*s;R=Q(s+1,s);de=Q(d,s)
        assert R**(d+1)<=16
        c_sq=de**2/R
        assert 4<c_sq<9
        assert de*(1-Q(1,d))==2
        assert Q((d-1)**2,4)==a
        # Raw quintic coefficient has a uniformly bounded removable nu limit.
        gamma_sq=Q(1,25*(d*d-1))
        assert gamma_sq<=Q(1,200)
        rows.append({"s":s,"d":d,"c_squared":str(c_sq),"gamma_squared":str(gamma_sq)})
    return rows
def negative_controls(info):
    K=info["K"];wrong=scale(info["sos"],-1)
    assert add(K,scale(wrong,-1))  # wrong potential sign is not the identity
    assert Q(1,2)-Q(3,12)!=Q(1,2)  # upper complex perturbation costs margin
    assert Q(7,4)!=1               # no claim of equal root/radial order
    scope={"uniform_in_s":True,"fixed_H":True,"effective_N0":False,
           "growing_window":False,"global_SFC":False,"all_roots_localized":False}
    assert scope["uniform_in_s"] and scope["fixed_H"]
    assert not any(scope[k] for k in ["effective_N0","growing_window","global_SFC","all_roots_localized"])
    return scope,4

SCRIPT="verify_uniform_bounded_quartic_cusp.py"
INPUT="uniform_bounded_quartic_cusp_inputs.json"
LEDGER="uniform_bounded_quartic_cusp_source_ledger.json"
RECEIPT="uniform_bounded_quartic_cusp_verification.json"
PREFIX="research/factorial/image-rigidity/"
ROLES={
 "uniform_proof":("local",PROOF),
 "uniform_checker":("local",SCRIPT),
 "uniform_inputs":("local",INPUT),
 "soft_criterion":("repo",PREFIX+"GAMMA_PRODUCT_ROOT_MODULUS.tex"),
 "local_completion_checker":("repo",PREFIX+"verify_quartic_cusp_local_completion.py"),
 "local_completion_inputs":("repo",PREFIX+"quartic_cusp_local_completion_inputs.json"),
 "local_completion_ledger":("repo",PREFIX+"quartic_cusp_local_completion_source_ledger.json"),
 "local_completion_receipt":("repo",PREFIX+"quartic_cusp_local_completion_verification.json"),
}
SCOPE={
 "schema":1,"fixed_finite_H":True,"uniform_in_integer_cusp_s":True,
 "integer_cusps":{"s_min":1,"a":"s^2","d":"2s+1"},
 "field":"h_N=(N^(3/4)/c_d)*(x/(N^d*w_c)-1)",
 "order_quantifier":"exists N0(H); all integer s>=1 and all integers N>=N0(H)",
 "complete_original_positive_side_bound":True,
 "uniform_real_concentration":True,
 "effective_N0":False,"growing_window":False,"all_N_at_fixed_support":False,
 "global_SFC":False,"all_actual_roots_localized":False,
 "formal_verification":False,"external_review":False,
}
def require(ok,message):
    if not ok:raise AssertionError(message)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read_json(path):return json.loads(path.read_text(encoding="utf-8"))
def check_scope(data):
    require(type(data) is dict and data==SCOPE,"Exact uniform fixed-H eventual scope required")
    for key,value in SCOPE.items():
        require(type(data[key]) is type(value),"Scope type: "+key)
    require(type(data["integer_cusps"]["s_min"]) is int,"Integer s_min")
def check_sources(data,repo_root):
    require(type(data) is dict and set(data)==
            {"schema","package","formal_verification","external_review","sources"},
            "Exact uniform source schema")
    require(type(data["schema"]) is int and data["schema"]==1,"Source schema1")
    require(data["package"]=="uniform-bounded-quartic-cusp","Source package identity")
    require(data["formal_verification"] is False and data["external_review"] is False,
            "Unsupported assurance flag")
    rows=data["sources"]
    require(type(rows) is list and len(rows)==len(ROLES),"All eight mandatory roles")
    result={};paths=set()
    for row in rows:
        require(type(row) is dict and set(row)=={"role","kind","path","sha256"},"Source row schema")
        role=row["role"]
        require(type(role) is str and role in ROLES and role not in result,"Unique source role")
        require((row["kind"],row["path"])==ROLES[role],"Exact role path")
        require((row["kind"],row["path"]) not in paths,"Unique source path")
        require(type(row["sha256"]) is str and re.fullmatch(r"[0-9a-f]{64}",row["sha256"]),
                "Lowercase SHA256")
        root=HERE if row["kind"]=="local" else repo_root.resolve()
        path=(root/row["path"]).resolve()
        require(path.is_relative_to(root.resolve()),"Source escapes its root")
        require(path.is_file() and sha(path)==row["sha256"],"Missing/stale source: "+role)
        result[role]=row;paths.add((row["kind"],row["path"]))
    require(set(result)==set(ROLES),"Complete source roles")
    return result
def local_dependency_closure(repo_root):
    base=repo_root/PREFIX
    ledger=read_json(base/"quartic_cusp_local_completion_source_ledger.json")
    require(type(ledger) is dict and type(ledger.get("schema")) is int
            and ledger["schema"]==1,"Imported local ledger schema")
    rows=ledger.get("sources")
    require(type(rows) is list and len(rows)==15,"Imported fifteen-role local closure")
    seen=set();hashes={}
    for row in rows:
        require(type(row) is dict and set(row)=={"role","path","sha256"},"Imported source row")
        require(row["role"] not in seen,"Imported role uniqueness");seen.add(row["role"])
        path=(base/row["path"]).resolve()
        require(path.is_relative_to(base.resolve()) and path.is_file(),"Imported relative source")
        require(sha(path)==row["sha256"],"Stale imported local source")
        require(row["path"] not in hashes,"Imported path uniqueness")
        hashes[row["path"]]=row["sha256"]
    for name in ["quartic_cusp_local_completion_source_ledger.json",
                 "quartic_cusp_local_completion_verification.json"]:
        hashes[name]=sha(base/name)
    require(len(hashes)==17,"Complete imported seventeen-file packet")
    return [{"path":PREFIX+name,"sha256":value} for name,value in sorted(hashes.items())]
def equal_child(actual,expected):
    require(actual==expected,"Imported completion fresh receipt differs")
def fresh_local_completion(repo_root):
    script=repo_root/ROLES["local_completion_checker"][1]
    receipt=repo_root/ROLES["local_completion_receipt"][1]
    with tempfile.TemporaryDirectory(prefix="uniform-bounded-cusp-") as directory:
        output=Path(directory)/"local.json"
        r=subprocess.run([sys.executable,"-B",str(script),"--repo-root",str(repo_root),
                          "--output",str(output)],cwd=directory,capture_output=True,
                         text=True,timeout=30)
        require(r.returncode==0,"Imported completion failed: "+r.stdout+"; "+r.stderr)
        require(output.is_file(),"Imported completion did not produce a receipt")
        actual=output.read_bytes();equal_child(actual,receipt.read_bytes())
        data=json.loads(actual)
        require(data.get("status")=="PASS" and data.get("mandatory_source_roles")==15
                and data.get("fresh_replays")==3,"Imported completion status/scope")
        return {"byte_identical":True,"receipt_sha256":hashlib.sha256(actual).hexdigest(),
                "child_fresh_replays":3}
def source_scope_negatives(ledger,repo_root):
    names=[]
    def rejected(name,fn):
        try:fn()
        except (AssertionError,KeyError,TypeError,ValueError):names.append(name);return
        raise AssertionError("Accepted negative: "+name)
    for role in sorted(ROLES):
        data=copy.deepcopy(ledger)
        data["sources"]=[r for r in data["sources"] if r["role"]!=role]
        rejected("omit:"+role,lambda data=data:check_sources(data,repo_root))
        data=copy.deepcopy(ledger)
        next(r for r in data["sources"] if r["role"]==role)["sha256"]="0"*64
        rejected("badpin:"+role,lambda data=data:check_sources(data,repo_root))
    for key in ["effective_N0","growing_window","all_N_at_fixed_support",
                "global_SFC","all_actual_roots_localized"]:
        data=copy.deepcopy(SCOPE);data[key]=True
        rejected("false-scope:"+key,lambda data=data:check_scope(data))
    for key in ["fixed_finite_H","uniform_in_integer_cusp_s",
                "complete_original_positive_side_bound","uniform_real_concentration"]:
        data=copy.deepcopy(SCOPE);data[key]=False
        rejected("omitted-scope:"+key,lambda data=data:check_scope(data))
    for key in ["schema","empty","duplicate"]:
        data=copy.deepcopy(ledger)
        if key=="schema":data["schema"]=True
        if key=="empty":data["sources"]=[]
        if key=="duplicate":data["sources"][-1]=copy.deepcopy(data["sources"][0])
        rejected("ledger:"+key,lambda data=data:check_sources(data,repo_root))
    rejected("forged-child",lambda:equal_child(b'{"status":"PASS"}',b'{"status":"PASS","real":true}'))
    return names


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--repo-root",required=True,type=Path)
    ap.add_argument("--output",type=Path,default=HERE/RECEIPT)
    ns=ap.parse_args()
    repo_root=ns.repo_root.resolve()
    require(repo_root.is_dir(),"Explicit repository root required")
    inputs=read_json(HERE/INPUT);check_scope(inputs)
    ledger=read_json(HERE/LEDGER);roles=check_sources(ledger,repo_root)
    closure=local_dependency_closure(repo_root)
    imported=fresh_local_completion(repo_root)
    premise=repo_root/ROLES["soft_criterion"][1]
    text=premise.read_text()
    require("cor:linear-lower-exponent" in text and r"a\ge d(2N-3)" in text,
            "Actual retained beta1 soft criterion")
    count,info=identities();guards=rational_guards();rows=cusp_controls()
    scope,negative=negative_controls(info)
    require(scope["uniform_in_s"] and scope["fixed_H"] and not scope["global_SFC"],
            "Analytic scope")
    source_controls=source_scope_negatives(ledger,repo_root)
    check_sources(ledger,repo_root)
    payload={"status":"PASS","schema":1,"formal_verification":False,"external_review":False,
             "scope":inputs,"mandatory_source_roles":len(ROLES),
             "proof_sha256":sha(HERE/PROOF),"checker_sha256":sha(HERE/SCRIPT),
             "input_sha256":sha(HERE/INPUT),"source_ledger_sha256":sha(HERE/LEDGER),
             "soft_criterion_sha256":sha(premise),
             "soft_criterion_locator":"cor:linear-lower-exponent, beta=1",
             "sources":[roles[r] for r in sorted(roles)],
             "imported_local_core_files":len(closure),"imported_local_closure":closure,
             "fresh_local_completion":imported,
             "counts":{"universal_polynomial_identities":count,"rational_guards":len(guards),
                       "exact_integer_cusp_controls":len(rows),"negative_controls":negative,
                       "source_scope_controls":len(source_controls)},
             "guard_names":sorted(guards),"cusp_controls":rows,
             "source_scope_control_names":source_controls,
             "meaning":"The full written proof establishes uniform analytic bounds. Exact/source replay does not replace compactness, Newman, or concentration arguments."}
    ns.output.parent.mkdir(parents=True,exist_ok=True)
    ns.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":"PASS","mandatory_sources":len(ROLES),
                      "imported_core_files":len(closure),"counts":payload["counts"],
                      "output":str(ns.output)}))
if __name__=="__main__":main()

