#!/usr/bin/env sage -python
"""
Recover the exact H92 q24 section by reconstructing WHOLE modular sections
computed at many independent primes.

Strategy
========
1. Reuse the proven fresh-prime worker:
     recover_h92_q24_degree46_direct_global_modp_reuse_sign.sage
   It computes the normalized q24 section
       Z (monic deg 24), X (deg 52), Y (deg 78)
   over GF(p) from the degree-46 direct trace.

2. Run that worker for several primes in parallel.

3. CRT the 156 unknown normalized coefficients:
       Z[0..23], X[0..52], Y[0..78].

4. Rational-reconstruct coefficient-wise.  Continue adding primes until every
   coefficient reconstructs and the FULL exact identity
       Y^2 = X^3 + A X Z^4 + B Z^6
   holds over QQ[U].

5. Require literal reduction to the independently certified GF(100003)
   degree-46 section.

Optional acceleration
=====================
If artifacts/local/elkies-k3/q24-direct-hensel-p512.json exists, its p-adic
checkpoint can seed the CRT with thousands of bits.  It is used ONLY after its
mod-p reduction is independently reproduced by the fresh-prime global worker,
so a stale/wrong branch cannot silently contaminate the lift.

This intentionally differs from the slow exact-specialization path: all heavy
degree-46 traces happen over machine finite fields; characteristic zero appears
only in CRT/rational reconstruction and the final exact identity.
"""

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, gcd, matrix


def locate_repo(explicit=None):
    candidates=[]
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd=Path.cwd().resolve()
    candidates += [cwd,*cwd.parents]
    h=Path.home()
    candidates += [
        h/"Documents"/"jacobian-research",
        h/"jacobian-research",
        h/"src"/"jacobian-research",
        h/"git"/"jacobian-research",
        h/"projects"/"jacobian-research",
    ]
    seen=set()
    for c in candidates:
        try:
            c=c.resolve()
        except Exception:
            continue
        if c in seen:
            continue
        seen.add(c)
        if (c/"elkies-k3/scripts").is_dir() and (c/"artifacts/generated-results").is_dir():
            return c
    raise SystemExit("Could not locate jacobian-research")


parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo",type=Path)
parser.add_argument("--workers",type=int,default=min(8,max(1,os.cpu_count() or 1)))
parser.add_argument("--batch-size",type=int,default=None)
parser.add_argument("--max-fresh-primes",type=int,default=192)
parser.add_argument("--prime-start",type=int,default=1000000007)
parser.add_argument("--samples-per-prime",type=int,default=105)
parser.add_argument("--scan-limit",type=int,default=500)
parser.add_argument("--no-hensel-seed",action="store_true")
parser.add_argument("--work-dir",type=Path)
parser.add_argument("--output",type=Path)
args=parser.parse_args()

if args.workers<1:
    raise ValueError("workers must be positive")
if args.samples_per_prime<105:
    raise ValueError("the current modular worker requires at least 105 samples")

ROOT=locate_repo(args.repo)
LOCAL=ROOT/"artifacts/local/elkies-k3"
GEN=ROOT/"artifacts/generated-results"
SCRIPTS=ROOT/"elkies-k3/scripts"

WORKER=SCRIPTS/"recover_h92_q24_degree46_direct_global_modp_reuse_sign.sage"
SEED=LOCAL/"q24-degree46-direct-global-mod-100003.json"
HENSEL=LOCAL/"q24-direct-hensel-p512.json"

q8_candidates=[
    LOCAL/"q8-corrected2cover-qq-child.json",
    GEN/"elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
]
CHILD=next((
    p for p in q8_candidates
    if p.exists()
    and "child" in json.loads(p.read_text())
    and "minimal_A_coefficients_low_to_high"
        in json.loads(p.read_text()).get("child",{})
),None)

for p in (WORKER,SEED):
    if not p.exists():
        raise SystemExit(f"Missing prerequisite: {p}")
if CHILD is None:
    raise SystemExit("No complete corrected q8 child artifact")

WORKDIR=(
    args.work_dir.resolve()
    if args.work_dir
    else LOCAL/"q24-multiprime-sections"
)
WORKDIR.mkdir(parents=True,exist_ok=True)

OUTPUT=(
    args.output.resolve()
    if args.output
    else LOCAL/"q8-q24-horizontal-section-qq-multiprime.json"
)

sage_exe=shutil.which("sage")
if not sage_exe:
    raise SystemExit("Could not find `sage` executable in PATH")

seed=json.loads(SEED.read_text())
child=json.loads(CHILD.read_text())
assert seed["status"]=="PASS_MODULAR_Q24_FROM_DIRECT_DEGREE46_BRIDGE"
assert int(seed["prime"])==100003
assert child["status"]=="PASS_EXACT_CORRECTED_Q8_D13_CHILD"

global_sign=int(seed["direct_global_map"]["global_scale_root_sign"])
assert global_sign in (-1,1)

A_QQ=[QQ(v) for v in child["child"]["minimal_A_coefficients_low_to_high"]]
B_QQ=[QQ(v) for v in child["child"]["minimal_B_coefficients_low_to_high"]]

print(
    "Q24MULTI_SETUP|"
    f"workers={args.workers}|samples_per_prime={args.samples_per_prime}|"
    f"global_sign={global_sign:+d}|seed_prime=100003|status=PASS",
    flush=True,
)


def section_vector(artifact):
    sec=artifact["section_mod_p"]
    z=[ZZ(v) for v in sec["Z_coefficients_low_to_high"]]
    x=[ZZ(v) for v in sec["X_coefficients_low_to_high"]]
    y=[ZZ(v) for v in sec["Y_coefficients_low_to_high"]]
    if not (len(z)==25 and len(x)==53 and len(y)==79):
        raise ArithmeticError(
            f"bad coefficient lengths Z/X/Y={len(z)}/{len(x)}/{len(y)}"
        )
    p=ZZ(artifact["prime"])
    if z[-1]%p != 1:
        raise ArithmeticError("Z is not monic")
    return [v%p for v in z[:-1]] + [v%p for v in x] + [v%p for v in y]


def load_mod_artifact(path, expected_prime=None):
    try:
        a=json.loads(path.read_text())
        if a.get("status")!="PASS_MODULAR_Q24_FROM_DIRECT_DEGREE46_BRIDGE":
            return None
        p=ZZ(a["prime"])
        if expected_prime is not None and p!=ZZ(expected_prime):
            return None
        if int(a["direct_global_map"]["global_scale_root_sign"])!=global_sign:
            return None
        v=section_vector(a)
        if len(v)!=156:
            return None
        return a
    except Exception:
        return None


def worker_path(p):
    return WORKDIR/f"q24-degree46-direct-global-mod-{int(p)}.json"


def worker_log_path(p):
    return WORKDIR/f"q24-degree46-direct-global-mod-{int(p)}.log"


def run_prime(p):
    p=ZZ(p)
    out=worker_path(p)
    cached=load_mod_artifact(out,p)
    if cached is not None:
        return {
            "prime":p,
            "artifact":cached,
            "seconds":0.0,
            "cached":True,
            "ok":True,
        }

    cmd=[
        sage_exe,"-python",str(WORKER),
        "--repo",str(ROOT),
        "--prime",str(p),
        "--global-sign",str(global_sign),
        "--samples",str(args.samples_per_prime),
        "--scan-limit",str(args.scan_limit),
        "--output",str(out),
    ]
    started=time.monotonic()
    proc=subprocess.run(
        cmd,
        cwd=str(ROOT),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    elapsed=time.monotonic()-started
    worker_log_path(p).write_text(proc.stdout)

    if proc.returncode:
        return {
            "prime":p,
            "artifact":None,
            "seconds":elapsed,
            "cached":False,
            "ok":False,
            "reason":f"exit_{proc.returncode}",
        }

    artifact=load_mod_artifact(out,p)
    if artifact is None:
        return {
            "prime":p,
            "artifact":None,
            "seconds":elapsed,
            "cached":False,
            "ok":False,
            "reason":"invalid_artifact",
        }

    return {
        "prime":p,
        "artifact":artifact,
        "seconds":elapsed,
        "cached":False,
        "ok":True,
    }


def run_batch(primes,label="fresh"):
    results=[]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures={pool.submit(run_prime,p):ZZ(p) for p in primes}
        for fut in as_completed(futures):
            p=futures[fut]
            try:
                r=fut.result()
            except Exception as exc:
                r={
                    "prime":p,
                    "artifact":None,
                    "seconds":0.0,
                    "cached":False,
                    "ok":False,
                    "reason":type(exc).__name__+":"+str(exc),
                }
            results.append(r)
            if r["ok"]:
                print(
                    "Q24MULTI_WORKER|"
                    f"kind={label}|prime={p}|seconds={r['seconds']:.3f}|"
                    f"cached={int(r['cached'])}|status=PASS",
                    flush=True,
                )
            else:
                print(
                    "Q24MULTI_WORKER|"
                    f"kind={label}|prime={p}|seconds={r['seconds']:.3f}|"
                    f"reason={r.get('reason','unknown')}|status=SKIP",
                    flush=True,
                )
    return results


def crt_merge_vec(residues,M,vec,p):
    p=ZZ(p)
    if M%p==0:
        raise ArithmeticError(f"duplicate CRT prime {p}")
    inv=ZZ(pow(int(M%p),-1,int(p)))
    out=[]
    for r,a in zip(residues,vec):
        t=((ZZ(a)-ZZ(r)%p)*inv)%p
        out.append(ZZ(r)+M*t)
    return out,M*p


def rational_reconstruct(a,m):
    """Standard symmetric rational reconstruction, uniqueness bound sqrt(m/2)."""
    a=ZZ(a)%ZZ(m)
    m=ZZ(m)
    if not a:
        return QQ.zero()

    B=ZZ(math.isqrt(int(m//2)))
    r0,r1=m,a
    t0,t1=ZZ.zero(),ZZ.one()

    while r1 and abs(r1)>B:
        q=r0//r1
        r0,r1=r1,r0-q*r1
        t0,t1=t1,t0-q*t1

    if not r1 or not t1 or abs(r1)>B or abs(t1)>B:
        return None

    n,d=ZZ(r1),ZZ(t1)
    if d<0:
        n,d=-n,-d
    g=n.gcd(d)
    if g!=1:
        n//=g
        d//=g
    if d<=0 or (a*d-n)%m:
        return None
    return QQ(n)/QQ(d)


def _centered_mod(r,M):
    r=ZZ(r)%M
    return r-M if r>M//2 else r


def _srr_candidate_from_row(row, original, M):
    d=ZZ(row[0])
    if not d or gcd(d,M)!=1:
        return None
    nums=[ZZ(x) for x in row[1:]]
    if d<0:
        d=-d
        nums=[-x for x in nums]

    g=abs(gcd([d]+nums))
    variants=[]
    if g>1:
        variants.append((d//g,[a//g for a in nums]))
    variants.append((d,nums))

    for dd,aa in variants:
        if not dd or gcd(dd,M)!=1:
            continue
        if all((aa[i]-dd*ZZ(original[i]))%M==0 for i in range(len(aa))):
            return [QQ(aa[i])/QQ(dd) for i in range(len(aa))]
    return None


def _srr_chunk(values,M,label):
    n=len(values)
    sym=[_centered_mod(v,M) for v in values]

    B=matrix(ZZ,n+1,n+1)
    B[0,0]=1
    for j,r in enumerate(sym):
        B[0,j+1]=r
    for i in range(n):
        B[i+1,i+1]=M

    started_srr=time.monotonic()
    L=B.LLL(delta=0.99)
    elapsed=time.monotonic()-started_srr

    candidates=[]
    seen=set()
    for row in L.rows():
        q=_srr_candidate_from_row(row,values,M)
        if q is None:
            continue
        key=tuple(q)
        if key in seen:
            continue
        seen.add(key)
        candidates.append(q)

    if not candidates:
        print(
            "Q24MULTI_SRR|"
            f"label={label}|n={n}|seconds={elapsed:.3f}|"
            "status=NO_VECTOR",
            flush=True,
        )
        return None

    def height(qs):
        return max(
            max(abs(ZZ(q.numerator())).nbits(),
                abs(ZZ(q.denominator())).nbits())
            for q in qs
        )

    candidates.sort(key=height)
    q=candidates[0]
    print(
        "Q24MULTI_SRR|"
        f"label={label}|n={n}|seconds={elapsed:.3f}|"
        f"candidates={len(candidates)}|max_bits={height(q)}|"
        "status=CANDIDATE",
        flush=True,
    )
    return q


_SRR_BLOCKS={
    "Z":(0,24),
    "X":(24,77),
    "Y":(77,156),
}


def _srr_reconstruct_chunks(residues,M,chunk_size):
    out=[None]*156
    for block,(lo,hi) in _SRR_BLOCKS.items():
        pos=lo
        local=0
        while pos<hi:
            end=min(pos+chunk_size,hi)
            q=_srr_chunk(
                residues[pos:end],
                M,
                f"{block}[{local}:{local+(end-pos)}]/c{chunk_size}"
            )
            if q is None:
                return None
            out[pos:end]=q
            local+=end-pos
            pos=end
    return out


def _best_product_reconstruct(a,m):
    # Asymmetric rational reconstruction via Euclidean convergents.
    a=ZZ(a)%ZZ(m)
    m=ZZ(m)
    if not a:
        return QQ.zero(),0,None

    r0,r1=m,a
    t0,t1=ZZ.zero(),ZZ.one()
    cand=[]

    while r1:
        n=ZZ(r1)
        d=ZZ(t1)
        if d:
            if d<0:
                n,d=-n,-d
            if gcd(d,m)==1 and gcd(abs(n),d)==1 and (a*d-n)%m==0:
                prod=abs(n)*d
                if prod and 2*prod < m:
                    cand.append((prod,QQ(n)/QQ(d)))

        q=r0//r1
        r0,r1=r1,r0-q*r1
        t0,t1=t1,t0-q*t1

    if not cand:
        return None,None,None

    uniq={}
    for prod,q in cand:
        if q not in uniq or prod<uniq[q]:
            uniq[q]=prod
    ranked=sorted((prod,q) for q,prod in uniq.items())
    best_prod,best=ranked[0]

    gap=None
    if len(ranked)>1:
        gap=max(0,ranked[1][0].nbits()-best_prod.nbits())
    return best,best_prod.nbits(),gap


def candidate_from_crt(residues,M):
    # First strict symmetric reconstruction.
    coeffs=[]
    unresolved=[]
    for i,r in enumerate(residues):
        q=rational_reconstruct(r,M)
        if q is None:
            coeffs.append(None)
            unresolved.append(i)
        else:
            coeffs.append(q)

    if not unresolved:
        return coeffs,[]

    # Then allow highly asymmetric numerator/denominator sizes.
    recovered=0
    products=[]
    gaps=[]
    for i in list(unresolved):
        q,pbits,gap=_best_product_reconstruct(residues[i],M)
        if q is None:
            continue
        coeffs[i]=q
        recovered+=1
        products.append(pbits)
        if gap is not None:
            gaps.append(gap)

    unresolved=[i for i,q in enumerate(coeffs) if q is None]
    print(
        "Q24MULTI_ASYM|"
        f"crt_bits={M.nbits()}|recovered={recovered}|"
        f"remaining={len(unresolved)}|"
        f"max_product_bits={max(products) if products else -1}|"
        f"min_gap_bits={min(gaps) if gaps else -1}|status=TRY",
        flush=True,
    )

    # Exact verification of the full section remains the acceptance test.
    return coeffs,unresolved

def split_coeffs(coeffs):
    if len(coeffs)!=156:
        raise ArithmeticError("wrong coefficient count")
    z=list(coeffs[:24])+[QQ.one()]
    x=list(coeffs[24:77])
    y=list(coeffs[77:156])
    return z,x,y


def reduce_q(q,F,p):
    q=QQ(q)
    d=ZZ(q.denominator())
    if d%p==0:
        raise ZeroDivisionError
    return F(ZZ(q.numerator()))/F(d)


def exact_verify(coeffs,artifacts):
    if any(q is None for q in coeffs):
        return None

    RQ=PolynomialRing(QQ,"U")
    Zq,Xq,Yq=split_coeffs(coeffs)
    Z=RQ(Zq)
    X=RQ(Xq)
    Y=RQ(Yq)
    A=RQ(A_QQ)
    B=RQ(B_QQ)

    if not (
        Z.degree()==24 and Z.leading_coefficient()==1
        and X.degree()==52 and Y.degree()==78
    ):
        return None

    identity=Y**2-X**3-A*X*Z**4-B*Z**6
    if identity:
        return None

    # Mandatory independent seed-prime replay.
    for art in [seed]+list(artifacts):
        p=ZZ(art["prime"])
        F=GF(p)
        RF=PolynomialRing(F,"U")
        try:
            Zr=RF([reduce_q(v,F,p) for v in Z.list()])
            Xr=RF([reduce_q(v,F,p) for v in X.list()])
            Yr=RF([reduce_q(v,F,p) for v in Y.list()])
        except ZeroDivisionError:
            # A prime dividing an exact coefficient denominator is an unlucky
            # normalization prime.  Such a prime must never be part of CRT.
            return None

        sec=art["section_mod_p"]
        Za=RF([F(int(v)) for v in sec["Z_coefficients_low_to_high"]])
        Xa=RF([F(int(v)) for v in sec["X_coefficients_low_to_high"]])
        Ya=RF([F(int(v)) for v in sec["Y_coefficients_low_to_high"]])
        if Zr!=Za or Xr!=Xa or Yr!=Ya:
            return None

    max_num_bits=max(
        abs(ZZ(q.numerator())).nbits()
        for q in list(Z)+list(X)+list(Y)
    )
    max_den_bits=max(
        abs(ZZ(q.denominator())).nbits()
        for q in list(Z)+list(X)+list(Y)
    )
    return Z,X,Y,max_num_bits,max_den_bits


def write_exact(result,artifacts,M,hensel_used):
    Z,X,Y,max_num_bits,max_den_bits=result
    payload={
        "schema":"elkies-k3.h92-q8-q24-horizontal-section-qq.multiprime.v1",
        "status":"PASS_EXACT_Q24_HORIZONTAL_SECTION",
        "zero":"II*_E8_1_branch_anchor",
        "formula":"Q24 = AJ(Qmap-S3) + 2*G1",
        "method":"parallel whole-section multimodular CRT + rational reconstruction",
        "bridge":{
            "formula":"Qmap-S3",
            "q8_degree":46,
        },
        "profile":{
            "P_dot_O":24,
            "height":"52",
            "D13_local_correction":"0",
            "Z_degree":24,
            "X_degree":52,
            "Y_degree":78,
            "x_degrees":[52,48],
            "y_degrees":[78,72],
        },
        "section":{
            "Z_coefficients_low_to_high":[str(v) for v in Z.list()],
            "X_coefficients_low_to_high":[str(v) for v in X.list()],
            "Y_coefficients_low_to_high":[str(v) for v in Y.list()],
            "x_numerator_coefficients_low_to_high":[str(v) for v in X.list()],
            "x_denominator_coefficients_low_to_high":[str(v) for v in (Z**2).list()],
            "y_numerator_coefficients_low_to_high":[str(v) for v in Y.list()],
            "y_denominator_coefficients_low_to_high":[str(v) for v in (Z**3).list()],
            "x":"X/Z^2",
            "y":"Y/Z^3",
        },
        "verification":{
            "exact_weierstrass_identity":True,
            "seed_mod_100003_identical":True,
            "all_fresh_prime_sections_replayed":True,
        },
        "reconstruction":{
            "crt_modulus_bits":int(M.nbits()),
            "fresh_primes":[int(a["prime"]) for a in artifacts],
            "fresh_prime_count":len(artifacts),
            "hensel_seed_used":bool(hensel_used),
            "max_numerator_bits":int(max_num_bits),
            "max_denominator_bits":int(max_den_bits),
        },
        "inputs":{
            "seed_modular_section":str(SEED.relative_to(ROOT)),
            "modular_worker":str(WORKER.relative_to(ROOT)),
            "child":str(CHILD.relative_to(ROOT)),
        },
        "next":(
            "Use this exact q24 section to construct the exact q24 isotropic "
            "divisor / resolved Riemann-Roch pencil and continue to D12/MW5."
        ),
    }
    OUTPUT.parent.mkdir(parents=True,exist_ok=True)
    OUTPUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(f"OUTPUT|{OUTPUT}",flush=True)
    print(
        "Q24MULTI_RESULT|"
        f"fresh_primes={len(artifacts)}|crt_bits={M.nbits()}|"
        f"max_num_bits={max_num_bits}|max_den_bits={max_den_bits}|"
        "identity=PASS|mod100003=IDENTICAL|status=PASS_EXACT_Q24",
        flush=True,
    )


# ---------------------------------------------------------------------------
# Start CRT with independently certified GF(100003) whole section.
# ---------------------------------------------------------------------------
seed_vec=section_vector(seed)
crt_res=[ZZ(v) for v in seed_vec]
crt_M=ZZ(seed["prime"])
fresh_artifacts=[]
used_primes={ZZ(seed["prime"])}
hensel_used=False

# ---------------------------------------------------------------------------
# Optional p-adic seed: validate its mod-p section independently first.
# ---------------------------------------------------------------------------
if not args.no_hensel_seed and HENSEL.exists():
    try:
        hen=json.loads(HENSEL.read_text())
        valid=(
            hen.get("schema")=="elkies-k3.h92-q24-direct-hensel-lift.v1"
            and hen.get("status")=="PASS_Q24_DIRECT_HENSEL"
            and int(hen.get("jacobian_rank",-1))==156
            and int(hen.get("final_residual_valuation",0))>=int(hen.get("precision",1))
            and len(hen.get("residues",[]))==156
        )
        if valid:
            hp=ZZ(hen["prime"])
            hN=int(hen["precision"])
            hM=hp**hN

            print(
                "Q24MULTI_HENSEL|"
                f"prime={hp}|precision={hN}|modulus_bits={hM.nbits()}|"
                "stage=VALIDATE_MOD_P|status=BEGIN",
                flush=True,
            )
            hv=run_prime(hp)
            if hv["ok"]:
                modvec=section_vector(hv["artifact"])
                hres=[ZZ(v)%hM for v in hen["residues"]]
                if all(hres[i]%hp==modvec[i]%hp for i in range(156)):
                    # The initial CRT seed is already modulo hp=100003.
                    # Upgrade that residue to the compatible hp^N checkpoint
                    # instead of trying to CRT two non-coprime moduli.
                    if crt_M == hp:
                        if any((hres[i]-crt_res[i])%hp for i in range(156)):
                            raise ArithmeticError("Hensel seed disagrees with base-prime CRT seed")
                        crt_res=list(hres)
                        crt_M=hM
                        hensel_used=True
                        used_primes.add(hp)
                        print(
                            "Q24MULTI_HENSEL|"
                            f"prime={hp}|precision={hN}|"
                            f"combined_modulus_bits={crt_M.nbits()}|"
                            "modp_match=1|status=ACCEPTED_SEED",
                            flush=True,
                        )
                    elif crt_M%hp:
                        # Generic fallback if a future invocation reaches this
                        # point with a coprime accumulated modulus.
                        inv=ZZ(pow(int(crt_M%hM),-1,int(hM)))
                        new=[]
                        for r,a in zip(crt_res,hres):
                            t=((a-r)%hM*inv)%hM
                            new.append(r+crt_M*t)
                        crt_res=new
                        crt_M*=hM
                        hensel_used=True
                        used_primes.add(hp)
                        print(
                            "Q24MULTI_HENSEL|"
                            f"prime={hp}|precision={hN}|"
                            f"combined_modulus_bits={crt_M.nbits()}|"
                            "modp_match=1|status=ACCEPTED_SEED",
                            flush=True,
                        )
                else:
                    print(
                        "Q24MULTI_HENSEL|modp_match=0|status=IGNORED",
                        flush=True,
                    )
            else:
                print(
                    "Q24MULTI_HENSEL|"
                    f"reason={hv.get('reason','worker_failed')}|status=IGNORED",
                    flush=True,
                )
    except Exception as exc:
        print(
            "Q24MULTI_HENSEL|"
            f"reason={type(exc).__name__}:{exc}|status=IGNORED",
            flush=True,
        )

# See whether the seed(s) alone already suffice.
coeffs,unresolved=candidate_from_crt(crt_res,crt_M)
print(
    "Q24MULTI_RECON|"
    f"fresh_primes=0|crt_bits={crt_M.nbits()}|"
    f"resolved={156-len(unresolved)}|unresolved={len(unresolved)}|status=TRY",
    flush=True,
)
verified=exact_verify(coeffs,[])
if verified:
    write_exact(verified,[],crt_M,hensel_used)
    raise SystemExit(0)


# ---------------------------------------------------------------------------
# Include any previously completed fresh-prime artifacts from WORKDIR.
# ---------------------------------------------------------------------------
cached=[]
for path in sorted(WORKDIR.glob("q24-degree46-direct-global-mod-*.json")):
    art=load_mod_artifact(path)
    if art is None:
        continue
    p=ZZ(art["prime"])
    if p in used_primes:
        continue
    cached.append(art)

for art in cached:
    p=ZZ(art["prime"])
    crt_res,crt_M=crt_merge_vec(crt_res,crt_M,section_vector(art),p)
    fresh_artifacts.append(art)
    used_primes.add(p)

if cached:
    coeffs,unresolved=candidate_from_crt(crt_res,crt_M)
    print(
        "Q24MULTI_RECON|"
        f"fresh_primes={len(fresh_artifacts)}|crt_bits={crt_M.nbits()}|"
        f"resolved={156-len(unresolved)}|unresolved={len(unresolved)}|"
        "source=CACHE|status=TRY",
        flush=True,
    )
    verified=exact_verify(coeffs,fresh_artifacts)
    if verified:
        write_exact(verified,fresh_artifacts,crt_M,hensel_used)
        raise SystemExit(0)


# ---------------------------------------------------------------------------
# Generate fresh prime sections in parallel, reconstruct after every batch.
# ---------------------------------------------------------------------------
batch_size=args.batch_size or args.workers
candidate=ZZ(args.prime_start)
launched=0

while launched<args.max_fresh_primes:
    batch=[]
    while len(batch)<batch_size and launched+len(batch)<args.max_fresh_primes:
        candidate=candidate.next_prime()
        if candidate in used_primes:
            continue
        batch.append(candidate)

    if not batch:
        break

    print(
        "Q24MULTI_BATCH|"
        f"count={len(batch)}|first={batch[0]}|last={batch[-1]}|"
        f"fresh_before={len(fresh_artifacts)}|status=BEGIN",
        flush=True,
    )
    results=run_batch(batch)
    launched+=len(batch)

    good=0
    for r in sorted(results,key=lambda x:int(x["prime"])):
        if not r["ok"]:
            continue
        art=r["artifact"]
        p=ZZ(art["prime"])
        if p in used_primes:
            continue
        crt_res,crt_M=crt_merge_vec(crt_res,crt_M,section_vector(art),p)
        fresh_artifacts.append(art)
        used_primes.add(p)
        good+=1

    coeffs,unresolved=candidate_from_crt(crt_res,crt_M)
    print(
        "Q24MULTI_RECON|"
        f"fresh_primes={len(fresh_artifacts)}|batch_good={good}|"
        f"crt_bits={crt_M.nbits()}|"
        f"resolved={156-len(unresolved)}|unresolved={len(unresolved)}|status=TRY",
        flush=True,
    )

    if not unresolved:
        verified=exact_verify(coeffs,fresh_artifacts)
        if verified:
            write_exact(verified,fresh_artifacts,crt_M,hensel_used)
            raise SystemExit(0)
        print(
            "Q24MULTI_RECON|all_coefficients_reconstructed=1|"
            "exact_identity=FAIL|action=ADD_MORE_PRIMES",
            flush=True,
        )

raise RuntimeError(
    f"No exact q24 after {len(fresh_artifacts)} fresh prime sections; "
    f"CRT modulus has {crt_M.nbits()} bits. "
    f"Increase --max-fresh-primes and rerun; completed prime artifacts are cached."
)
