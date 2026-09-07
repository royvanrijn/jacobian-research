#!/usr/bin/env sage -python
import argparse, json
from pathlib import Path
from sage.all import QQ, ZZ, Zmod

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3"

parser=argparse.ArgumentParser()
parser.add_argument("--primes",required=True)
parser.add_argument("--output",type=Path)
args=parser.parse_args()

primes=[ZZ(x) for x in args.primes.split(",") if x.strip()]
if len(primes)<2:
    raise ValueError("need at least two signatures")

sigs=[]
for p in primes:
    path=LOCAL/f"q32-signature-mod-{p}.json"
    if not path.exists():
        raise SystemExit(f"missing {path}")
    d=json.loads(path.read_text())
    assert d["status"]=="PASS_Q32_MODP_SIGNATURE"
    assert ZZ(d["prime"])==p
    sigs.append(d)

ref=sigs[0]
for d in sigs[1:]:
    assert d["plane_pivots"]==ref["plane_pivots"]
    assert d["quartic_degree"]==ref["quartic_degree"]==4
    assert d["jacobian_A"]["num_degree"]==ref["jacobian_A"]["num_degree"]
    assert d["jacobian_A"]["den_degree"]==ref["jacobian_A"]["den_degree"]
    assert d["jacobian_B"]["num_degree"]==ref["jacobian_B"]["num_degree"]
    assert d["jacobian_B"]["den_degree"]==ref["jacobian_B"]["den_degree"]
    assert (d["child_root_rank"],d["child_root_det"])==(12,4)

def crt_scalar(residues,mods):
    x=ZZ(0)
    M=ZZ(1)
    for rr,pp in zip(residues,mods):
        pp=ZZ(pp)
        rr=ZZ(rr)%pp
        t=((rr-x)%pp)*((M%pp).inverse_mod(pp)) % pp
        x += M*t
        M *= pp
        x %= M
    return x,M

def rr_scalar(residues,mods):
    x,M=crt_scalar(residues,mods)
    try:
        q=QQ(Zmod(M)(x).rational_reconstruction())
    except (ValueError, ArithmeticError, ZeroDivisionError):
        return None,M
    return q,M

def reduce_q(q,p):
    q=QQ(q); p=ZZ(p)
    den=ZZ(q.denominator())%p
    if den==0:
        return None
    return int((ZZ(q.numerator())%p)*den.inverse_mod(p)%p)

# We deliberately reconstruct from all but the final prime and reserve the
# final prime as an independent held-out test.
train_primes=primes[:-1]
hold_p=primes[-1]
train=sigs[:-1]
hold=sigs[-1]

def reconstruct_array(getter):
    shape0=getter(train[0])
    vals=[]
    total=0
    recovered=0
    valid=0
    max_nb=0
    max_db=0

    # getter returns a flat list
    for j in range(len(shape0)):
        residues=[getter(d)[j] for d in train]
        q,M=rr_scalar(residues,train_primes)
        total+=1
        if q is None:
            vals.append(None)
            continue
        recovered+=1
        max_nb=max(max_nb,abs(ZZ(q.numerator())).nbits())
        max_db=max(max_db,abs(ZZ(q.denominator())).nbits())
        hv=reduce_q(q,hold_p)
        if hv is not None and hv==int(getter(hold)[j])%int(hold_p):
            valid+=1
        vals.append(q)
    return {
        "values":vals,
        "total":total,
        "recovered":recovered,
        "heldout_valid":valid,
        "max_num_bits":max_nb,
        "max_den_bits":max_db,
    }

def flat_plane(d):
    return [int(v) for row in d["plane_rref_2x56"] for v in row]

def flat_A(d):
    return [int(v) for v in d["jacobian_A"]["num"]] + [int(v) for v in d["jacobian_A"]["den"]]

def flat_B(d):
    return [int(v) for v in d["jacobian_B"]["num"]] + [int(v) for v in d["jacobian_B"]["den"]]

rp=reconstruct_array(flat_plane)
ra=reconstruct_array(flat_A)
rb=reconstruct_array(flat_B)

def stat(name,r):
    print(
        "Q32CRT_PROGRESS|"
        f"object={name}|train={len(train_primes)}|holdout={hold_p}|"
        f"recovered={r['recovered']}/{r['total']}|"
        f"heldout={r['heldout_valid']}/{r['total']}|"
        f"max_num_bits={r['max_num_bits']}|max_den_bits={r['max_den_bits']}|"
        "status="+(
            "PASS_HELDOUT" if r["heldout_valid"]==r["total"]
            else "PARTIAL"
        ),
        flush=True,
    )

stat("plane",rp)
stat("jacA",ra)
stat("jacB",rb)

complete=all(
    r["heldout_valid"]==r["total"]
    for r in (rp,ra,rb)
)

payload={
    "schema":"elkies-k3.h3-q32-crt-qq-candidate.v1",
    "status":(
        "PASS_Q32_QQ_CANDIDATE_HELDOUT_VALIDATED"
        if complete else
        "PARTIAL_Q32_CRT_RECONSTRUCTION"
    ),
    "training_primes":[int(p) for p in train_primes],
    "heldout_prime":int(hold_p),
    "plane_pivots":ref["plane_pivots"],
    "progress":{
        "plane":{k:v for k,v in rp.items() if k!="values"},
        "jacobian_A":{k:v for k,v in ra.items() if k!="values"},
        "jacobian_B":{k:v for k,v in rb.items() if k!="values"},
    },
}

if complete:
    pv=rp["values"]
    plane=[
        [str(pv[56*i+j]) for j in range(56)]
        for i in range(2)
    ]
    Avals=ra["values"]
    An=len(ref["jacobian_A"]["num"])
    Bvals=rb["values"]
    Bn=len(ref["jacobian_B"]["num"])
    payload["qq_candidate"]={
        "plane_rref_2x56":plane,
        "jacobian_A":{
            "num":[str(q) for q in Avals[:An]],
            "den":[str(q) for q in Avals[An:]],
            "num_degree":ref["jacobian_A"]["num_degree"],
            "den_degree":ref["jacobian_A"]["den_degree"],
        },
        "jacobian_B":{
            "num":[str(q) for q in Bvals[:Bn]],
            "den":[str(q) for q in Bvals[Bn:]],
            "num_degree":ref["jacobian_B"]["num_degree"],
            "den_degree":ref["jacobian_B"]["den_degree"],
        },
    }

OUT=args.output.resolve() if args.output else LOCAL/"q32-crt-qq-candidate.json"
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q32CRT_RESULT|"
    f"primes={','.join(map(str,primes))}|"
    f"train={len(train_primes)}|holdout={hold_p}|"
    f"status={payload['status']}",
    flush=True,
)
