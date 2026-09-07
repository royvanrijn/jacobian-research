#!/usr/bin/env python3
"""Bounded CT construction from the already certified projection quartics."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import projection_fibres as pf

PROTOCOL = pf.HERE/"PROJECTION_PAIR_PROTOCOL.json"
OUTPUT = r.OUT/"rank_jump_projection_pairings_v1.json"
WORK = r.ROOT/"artifacts/local/rank-jump-projection-pairings-v1"


def bindings():
    paths = [PROTOCOL,Path(__file__),pf.OUTPUT,pf.HERE/"verify_projection_fibres.py",
             pf.old.CAS/"research_runtime/fisher.py",pf.old.CAS/"research_runtime/sage_subspace.py"]
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}


def worker(index,path):
    import sage.version
    from sage.all import QQ,ZZ,matrix,vector,prod,pari
    from verify_projection_fibres import RetainedBackend,ArithmeticContext
    from research_runtime.fisher import fisher_gamma,hilbert_symbol
    from research_runtime.sage_subspace import quartic_local_witness
    data = r.read(pf.OUTPUT)
    assert data["bindings"] == pf.bindings()
    masks = r.read(pf.PROTOCOL)["pair_class_masks"][index]
    rows = [next(c for c in data["covers"] if c["mask"] == m) for m in masks]
    backend = RetainedBackend(ArithmeticContext.from_record(rows[0]["context"]))
    quartics = [backend.R(list(map(QQ,c["quartic"]))) for c in rows]
    reps = [backend.K(list(map(QQ,b))) for b in rows[0]["classes"]]
    root = prod(backend.K(list(map(QQ,c["cubic_invariant_over_beta_square_root"]))) for c in rows)
    exponents = [sum(m>>i&1 for m in masks) for i in range(len(reps))]
    assert all(e%2 == 0 for e in exponents)
    root *= prod(b**(e//2) for b,e in zip(reps,exponents))
    phi = -3*backend.alpha-backend.a2
    cubic_product = prod((4*q[4]*phi+3*q[3]**2-8*q[4]*q[2])/3 for q in quartics)
    assert root*root == cubic_product
    coordinates = matrix(QQ,3,3,lambda i,j:(phi**j).lift()[i])
    sqrt_phi = list(map(str,coordinates.inverse()*vector(QQ,[root.lift()[i] for i in range(3)])))
    print("ROOT verified by multiplication",flush=True)
    gamma = backend.R(fisher_gamma(quartics,sqrt_phi,backend.I,backend.J))
    known = sorted(set([2,3]+[p for p,e in rows[0]["context"]["discriminant_factorization"]]))
    def factor(value):
        n = ZZ(value)
        facts = []
        for p in known:
            e = n.valuation(p)
            if e:
                facts.append([p,int(e)])
                n //= ZZ(p)**e
        if n > 1:
            print("FACTOR remaining bits",n.nbits(),flush=True)
            facts.extend([int(p),int(e)] for p,e in n.factor(proof=True))
        assert prod(ZZ(p)**e for p,e in facts) == value
        return sorted(facts)
    q,other,_ = quartics
    values = [abs(ZZ(q.discriminant().numerator())),ZZ(q.discriminant().denominator()),
              ZZ(q.denominator()),abs(ZZ(other[4].numerator())),ZZ(other[4].denominator())]
    factors = [factor(v) for v in values]
    places,_ = backend._pair_support(quartics,retained=factors)
    print("SUPPORT verified",len(places),flush=True)
    terms = []
    answer = 1
    for place in places:
        print("LOCAL",place,flush=True)
        xx = quartic_local_witness(q,gamma,place,node_cap=r.read(PROTOCOL)["limits"]["local_witness_node_cap"])
        qv,gv = q(xx),gamma(xx)
        symbol = int(hilbert_symbol(other[4],gv,place))
        answer *= symbol
        terms.append({"place":place,"x":str(xx),"q_value":str(qv),"gamma_value":str(gv),"hilbert_symbol":symbol})
    result = {"masks":masks,"square_root_phi_coefficients":sqrt_phi,
              "gamma":[str(gamma[i]) for i in range(3)],"support_factors":factors,
              "local_terms":terms,"value":int(answer == -1)}
    assert backend._pair(None,masks,rows,retained=result) == result
    old = next(c["matrix"] for c in r.read(pf.lc.INPUT)["ct"] if c["u"] == -1)
    expected = old[{1:0,2:1,4:2}[masks[0]]][{1:0,2:1,4:2}[masks[1]]]
    assert result["value"] == expected
    r.write_new(path,{"bindings":bindings(),"index":index,"pair":result,
                     "expected_retained_CT_value":expected,
                     "software":{"sage":sage.version.version,"pari":str(pari.version())}})
    print("DONE",index,result["value"],flush=True)


def capture():
    WORK.mkdir(parents=True,exist_ok=True)
    records = []
    for i in range(3):
        path = WORK/f"pair-{i}.json"
        if not path.exists():
            with (WORK/f"pair-{i}.log").open("x") as log:
                try:
                    proc = subprocess.run(["sage","-python",str(Path(__file__).resolve()),"worker",
                                           "--index",str(i),"--destination",str(path)],
                                          cwd=r.ROOT,stdout=log,stderr=log,timeout=60)
                    failure = None if proc.returncode == 0 else "worker failure"
                except subprocess.TimeoutExpired:
                    failure = "60-second timeout"
                if failure:
                    log.write("\nUNKNOWN: "+failure+"\n")
                    r.write_new(path,{"bindings":bindings(),"index":i,"status":"UNKNOWN","reason":failure})
        row = r.read(path)
        assert row["bindings"] == bindings()
        records.append(row)
        print("checkpoint",i,row.get("status",row.get("pair",{}).get("value")),flush=True)
    r.write_new(OUTPUT,{"schema":"rank-jump.projection-pairings.v1","bindings":bindings(),"pairs":records})


def verify():
    from sage.all import QQ,pari
    from verify_projection_fibres import RetainedBackend,ArithmeticContext,verify as verify_covers
    verify_covers()
    data = r.read(OUTPUT)
    assert data["bindings"] == bindings()
    covers = r.read(pf.OUTPUT)["covers"]
    backend = RetainedBackend(ArithmeticContext.from_record(covers[0]["context"]))
    assert len(data["pairs"]) == 3
    for i,row in enumerate(data["pairs"]):
        assert row["bindings"] == bindings() and row["index"] == i
        if row.get("status") == "UNKNOWN":
            print("UNKNOWN retry",i,row["reason"])
            continue
        masks = r.read(pf.PROTOCOL)["pair_class_masks"][i]
        selected = [next(c for c in covers if c["mask"] == m) for m in masks]
        pair = backend._pair(None,masks,selected,retained=row["pair"])
        for term in pair["local_terms"]:
            a,b = QQ(selected[1]["quartic"][4]),QQ(term["gamma_value"])
            place = term["place"]
            symbol = (-1 if a < 0 and b < 0 else 1) if place == "infinity" else int(pari.hilbert(a,b,place))
            assert symbol == term["hilbert_symbol"]
        old = next(c["matrix"] for c in r.read(pf.lc.INPUT)["ct"] if c["u"] == -1)
        assert old[{1:0,2:1,4:2}[masks[0]]][{1:0,2:1,4:2}[masks[1]]] == row["expected_retained_CT_value"] == pair["value"]
        print("PASS projection pairing",masks,pair["value"],len(pair["local_terms"]),"Hilbert witnesses")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode",choices=("capture","worker","verify"))
    parser.add_argument("--index",type=int)
    parser.add_argument("--destination",type=Path)
    args = parser.parse_args()
    if args.mode == "capture":
        capture()
    elif args.mode == "verify":
        verify()
    else:
        worker(args.index,args.destination)

