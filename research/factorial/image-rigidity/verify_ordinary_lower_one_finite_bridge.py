#!/usr/bin/env python3
"""Complete full-degree modular verification of the lower-one finite bridge."""
if not __debug__:
    raise SystemExit("Refusing -O/-OO: exact verification checks are required.")
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor, as_completed
from hashlib import sha256
from math import comb, factorial, isqrt
from pathlib import Path
import copy, json, struct, sys
import flint
from flint import nmod_poly
HERE=Path(__file__).resolve().parent
PROOF="ORDINARY_LOWER_ONE_STRONG_FACTORIAL.md"
INPUTS="ordinary_lower_one_finite_bridge_inputs.json"
LEDGER="ordinary_lower_one_finite_bridge_source_ledger.json"
RECEIPT="ordinary_lower_one_finite_bridge_verification.json"
def require(test,message):
    if not test:raise ValueError(message)
def digest(path):return sha256(path.read_bytes()).hexdigest()
def domain(inputs):
    required={"schema_version","prime","minimum_gap","maximum_gap","minimum_terminal",
              "tail_minimum","tail_linear_factor","expected_cases","gaps"}
    require(set(inputs)==required,"invalid finite-domain schema")
    require([inputs[k]for k in ["schema_version","prime","minimum_gap","maximum_gap",
             "minimum_terminal","tail_minimum","tail_linear_factor","expected_cases"]]
            ==[1,1000003,10,43,2,5000,128,171812],"wrong mathematical domain")
    p=inputs["prime"]
    require(p>1 and all(p%j for j in range(2,isqrt(p)+1)),"coefficient modulus is not prime")
    require(len(inputs["gaps"])==34,"missing or duplicate gap family")
    cases=[]
    for d,row in zip(range(10,44),inputs["gaps"]):
        last=max(5000,128*(d+1))-1
        require(set(row)=={"gap","first_terminal","last_terminal","tested","case_stream_sha256"},
                "invalid expected gap record")
        require([row[k]for k in ["gap","first_terminal","last_terminal","tested"]]
                ==[d,2,last,last-1],"incomplete or reordered gap domain")
        require(isinstance(row["case_stream_sha256"],str) and len(row["case_stream_sha256"])==64,
                "invalid expected case digest")
        require((d+1)*last-1<p,"full-degree modular specialization can lose its endpoint")
        cases.append((d,last,p,row["case_stream_sha256"]))
    require(sum(last-1 for d,last,p,h in cases)==171812,"wrong complete case count")
    require(max((d+1)*last-1 for d,last,p,h in cases)==247763,"wrong factorial endpoint")
    return cases
def tables(p,limit):
    require(0<=limit<p,"factorial denominator is not a modular unit")
    fact=[1]*(limit+1)
    for j in range(1,limit+1):fact[j]=fact[j-1]*j%p
    inv=[1]*(limit+1);inv[limit]=pow(fact[limit],-1,p)
    for j in range(limit,0,-1):inv[j-1]=inv[j]*j%p
    return fact,inv
def coefficients(N,d,p,t):
    # Ascending reciprocal T_N coefficients, before reversing to monic A_N.
    fact,inv=t
    return [N*fact[N+d*k-1]*inv[k]*inv[N-k]%p for k in range(N+1)]
def squarefree(coefficients_ascending,p,degree):
    f=nmod_poly(coefficients_ascending,p)
    require(f.degree()==degree,"original polynomial degree was lost")
    df=f.derivative()
    require(df.degree()==degree-1,"derivative degree was lost")
    g=f.gcd(df)
    require(g.degree()==0 and int(g[0])==1,"finite case is not certified squarefree")
    return f
def gap_worker(args):
    d,last,p,expected=args
    tab=tables(p,(d+1)*last-1);stream=sha256()
    for N in range(2,last+1):
        c=coefficients(N,d,p,tab)
        require(c[0]==1 and c[-1]!=0,"reciprocal endpoint is zero")
        monic=list(reversed(c))
        require(monic[-1]==1,"A_N is not monic")
        squarefree(monic,p,N)
        stream.update(struct.pack("<IIII",d,N,p,0))
    actual=stream.hexdigest()
    require(actual==expected,"case stream differs from the complete retained bridge")
    return {"gap":d,"first_terminal":2,"last_terminal":last,"tested":last-1,
            "prime":p,"all_full_degree_gcds_equal_one":True,"case_stream_sha256":actual}
def controls(inputs):
    p=inputs["prime"];identities=0
    for d in [1,2,6,10,43]:
        tab=tables(p,(d+1)*9)
        for N in range(2,9):
            c=coefficients(N,d,p,tab)
            literal=[comb(N,k)*factorial(N+d*k-1)//factorial(N-1)%p for k in range(N+1)]
            require(c==literal,"small exact factorial normalization failed")
            f=nmod_poly(list(reversed(c)),p)
            manual=list(reversed([(N-k)*c[k]%p for k in range(N)]))
            require(list(map(int,f.derivative()))==manual,"literal derivative control failed")
            identities+=1
    require((factorial(2+3-1)//factorial(2-1))!=(factorial(2+3)//factorial(2)),
            "Gamma-shift negative control is invalid")
    repeated=nmod_poly([1,2,1],p)
    require(repeated.gcd(repeated.derivative()).degree()==1,"repeated-root negative control failed")
    # A full-degree polynomial corrupted to have the repeated root z=1.
    c=list(reversed(coefficients(4,10,p,tables(p,44))))
    c[1]=-sum(j*c[j]for j in range(2,len(c)))%p
    c[0]=-(c[1]+sum(c[2:]))%p
    try:squarefree(c,p,4)
    except ValueError:pass
    else:raise ValueError("actual repeated-root corruption was accepted")
    rejected=0
    variants=[]
    bad=copy.deepcopy(inputs);bad["gaps"].pop();variants.append(bad)
    bad=copy.deepcopy(inputs);bad["gaps"][-1]["last_terminal"]-=1;variants.append(bad)
    bad=copy.deepcopy(inputs);bad["prime"]=1000005;variants.append(bad)
    for bad in variants:
        try:domain(bad)
        except ValueError:rejected+=1
        else:raise ValueError("malformed finite domain was accepted")
    return {"exact_factorial_and_derivative_controls":identities,
            "malformed_domains_rejected":rejected,
            "negative_controls":["wrong_Gamma_shift","literal_repeated_root",
                                 "full_degree_repeated_root_corruption"]}
def root_path(explicit):
    if explicit:return explicit.resolve()
    for path in [HERE,*HERE.parents]:
        if(path/"research/MATH_STATUS.json").is_file():return path
    raise SystemExit("Use --repo-root when replaying a relocated packet.")
def main():
    ap=ArgumentParser(description=__doc__)
    ap.add_argument("--repo-root",type=Path)
    ap.add_argument("--jobs",type=int,default=4)
    ap.add_argument("--output","--write-receipt",dest="output",type=Path)
    args=ap.parse_args();require(1<=args.jobs<=32,"jobs must be between one and32")
    root=root_path(args.repo_root)
    source=json.loads((HERE/LEDGER).read_text())
    mandatory={("packet",PROOF),("packet",Path(__file__).name),("packet",INPUTS),
               ("repository","research/factorial/image-rigidity/GAP_SIX_RESULTANT_STRICT_POSITIVITY.md")}
    keys={(v["location"],v["path"])for v in source["files"]}
    require(mandatory<=keys and len(keys)==len(source["files"]),"missing essential source pins")
    pins=[]
    for item in source["files"]:
        path=(root if item["location"]=="repository"else HERE)/item["path"]
        require(path.is_file() and digest(path)==item["sha256"],"changed source:"+str(path))
        pins.append({"path":item["path"],"sha256":item["sha256"]})
    inputs=json.loads((HERE/INPUTS).read_text())
    cases=domain(inputs);checks=controls(inputs);rows=[]
    with ProcessPoolExecutor(max_workers=args.jobs)as pool:
        futures={pool.submit(gap_worker,c):c[0]for c in cases}
        for future in as_completed(futures):
            row=future.result();rows.append(row)
            print("PASS gap",row["gap"],"all",row["tested"],"full-degree cases",flush=True)
    rows.sort(key=lambda r:r["gap"])
    require(len(rows)==34 and sum(r["tested"]for r in rows)==171812,"finite bridge did not complete")
    result={"status":"PASS","scope":"Completed exact full-degree modular squarefreeness for EVERY required pair in10<=d<=43,2<=N<max(5000,128(d+1)); unrestricted theorem also requires the written analytic tail and small-gap inputs.",
            "arithmetic":"F_1000003, with preserved full polynomial and derivative degrees",
            "completed_cases":171812,"complete_gap_families":34,
            "maximum_terminal":5631,"maximum_factorial_argument":247763,
            "gaps":rows,"controls":checks,
            "versions":{"python":sys.version.split()[0],"python_flint":flint.__version__,
                        "FLINT":flint.__FLINT_VERSION__},
            "source_sha256":{"proof":digest(HERE/PROOF),"checker":digest(Path(__file__)),
                             "inputs":digest(HERE/INPUTS),"ledger":digest(HERE/LEDGER)},
            "verified_source_pins":pins,
            "external_review_or_formal_verification":False}
    output=args.output or HERE/RECEIPT
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("PASS complete finite bridge:171812 cases,34 gap families",flush=True)
if __name__=="__main__":main()
