#!/usr/bin/env python3
"""One frozen high-control prime, reusing the sealed construction and verifier."""
import argparse
from pathlib import Path
import subprocess
import tempfile
import retrospective as r
import equal_class_picard as ex
import verify_equal_class_picard as verifier

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/"EQUAL_CLASS_HIGH_FOLLOWUP_PROTOCOL.json"
OLD_RAW=ex.RAW
OLD_OUTPUT=ex.OUTPUT
OLD_VERIFICATION=verifier.OUTPUT
RAW=r.OUT/"rank_jump_equal_class_high_followup_counts_v1.json"
OUTPUT=r.OUT/"rank_jump_equal_class_high_followup_v1.json"
VERIFICATION=r.OUT/"rank_jump_equal_class_high_followup_verification_v1.json"
WORK=r.ROOT/"artifacts/local/rank-jump-equal-class-high-followup-v1"


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),PROTOCOL,ex.INPUT,OLD_RAW,OLD_OUTPUT,OLD_VERIFICATION,
             Path(ex.__file__),Path(ex.previous.__file__),Path(verifier.__file__))}


def setup():
    # Runtime adapters only: the sealed source files and old outputs stay intact.
    ex.RAW=RAW;ex.OUTPUT=OUTPUT;ex.bindings=bindings;ex.WORK=WORK


def worker():
    protocol=r.read(PROTOCOL);p=protocol["prime"];row=r.read(ex.INPUT)["cases"][0]
    tested=[{"p":q,"good":ex.previous.good(row,q,7)} for q in r.primes(p) if q>17]
    assert tested==protocol["eligibility"] and tested[-1]["good"]
    assert not any(v["good"] for v in tested[:-1])
    assert row["model"]==r.read(ex.previous.INPUT)["cases"][0]["model"]
    ex.previous.WORK=WORK
    record=ex.previous.new_count(0,p)
    record["producer_bindings"]=record.pop("bindings");record["bindings"]=bindings()
    r.write_new(WORK/"new-prime.json",record)


def capture():
    WORK.mkdir(parents=True,exist_ok=True);path=WORK/"new-prime.json"
    if not path.exists():
        with (WORK/"new-prime.log").open("x") as log:
            try:
                proc=subprocess.run(["sage","-python",str(Path(__file__).resolve()),"worker"],cwd=r.ROOT,
                                    stdout=log,stderr=log,timeout=r.read(PROTOCOL)["limits"]["seconds_count"])
                error=None if proc.returncode==0 else "worker failure"
            except subprocess.TimeoutExpired:error="40-second cap"
        if error:r.write_new(path,{"case":0,"p":29,"status":"UNKNOWN","reason":error,"bindings":bindings(),
                                  "transcript":(WORK/"new-prime.log").read_text()})
    record=r.read(path);assert record["bindings"]==bindings()
    r.write_new(RAW,{"schema":"rank-jump.equal-class-high-followup-counts.v1","bindings":bindings(),
                     "records":r.read(OLD_RAW)["records"]+[record]})
    print("checkpoint",record["status"],flush=True)


def independent(check=False):
    setup();counts=r.read(RAW);report=r.read(OUTPUT);old=r.read(OLD_VERIFICATION)
    assert counts["bindings"]==report["bindings"]==bindings()
    assert old["status"]=="PASS" and old["counts_sha256"]==r.digest(OLD_RAW.read_bytes())
    assert counts["records"][:-1]==r.read(OLD_RAW)["records"]
    assert report["reductions"][:-1]==r.read(OLD_OUTPUT)["reductions"]
    with tempfile.TemporaryDirectory(prefix="rank-jump-high-followup-") as tmp:
        subprocess.run(["python3",str(Path(__file__).resolve()),"verify-prime","--work",tmp],cwd=r.ROOT,
                       check=True,timeout=r.read(PROTOCOL)["limits"]["seconds_independent_replay"])
        new=r.read(Path(tmp)/"new-replay.json")
    geometry=verifier.geometry_check()
    result={"schema":"rank-jump.equal-class-high-followup-verification.v1","status":"PASS",
            "bindings":bindings(),"counts_sha256":r.digest(RAW.read_bytes()),"analysis_sha256":r.digest(OUTPUT.read_bytes()),
            "new_record":new,"reused_verification_sha256":r.digest(OLD_VERIFICATION.read_bytes()),
            "frobenius_checks":geometry}
    if check:assert result==r.read(VERIFICATION)
    else:r.write_new(VERIFICATION,result)
    print("PASS independent new-prime traces and all five Frobenius reductions")


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("mode",choices=["worker","capture","build","check","verify","replay","verify-prime"])
    parser.add_argument("--work",type=Path);args=parser.parse_args()
    if args.mode=="worker":worker()
    elif args.mode=="capture":capture()
    elif args.mode=="verify-prime":
        setup();ex.WORK=args.work
        r.write_new(args.work/"new-replay.json",verifier.prime_check(0,r.read(PROTOCOL)["prime"]))
    elif args.mode in ("verify","replay"):independent(args.mode=="replay")
    else:setup();ex.build(args.mode=="check")
