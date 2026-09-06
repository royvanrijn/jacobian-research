#!/usr/bin/env python3
"""Fresh bounded replay without reusing or overwriting local checkpoints."""
import argparse
from pathlib import Path
import subprocess
import tempfile
import retrospective as r
import verify_equal_class_picard as verifier


def replay():
    retained=r.read(verifier.OUTPUT)
    assert retained["status"]=="PASS"
    assert retained["verifier_sha256"]==r.digest(Path(verifier.__file__).read_bytes())
    assert retained["analysis_sha256"]==r.digest(verifier.ex.OUTPUT.read_bytes())
    assert retained["counts_sha256"]==r.digest(verifier.ex.RAW.read_bytes())
    with tempfile.TemporaryDirectory(prefix="rank-jump-equal-class-replay-") as name:
        work=Path(name)
        for record in retained["records"]:
            case,p=record["case"],record["p"]
            subprocess.run(["python3",str(Path(__file__).resolve()),"--case",str(case),"--prime",str(p),
                            "--work",str(work)],cwd=r.ROOT,check=True,timeout=60)
            assert r.read(work/f"prime-{case}-{p}.json")==record
    assert verifier.geometry_check()==retained["frobenius_checks"]
    print("PASS fresh independent traces, input provenance, and Frobenius geometry")


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--case",type=int)
    parser.add_argument("--prime",type=int);parser.add_argument("--work",type=Path)
    args=parser.parse_args()
    if args.work:
        verifier.ex.WORK=args.work
        r.write_new(args.work/f"prime-{args.case}-{args.prime}.json",verifier.prime_check(args.case,args.prime))
    else:replay()
