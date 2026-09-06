#!/usr/bin/env python3
"""Bounded descents on the two fixed small elliptic quotients."""
import argparse
from pathlib import Path
import subprocess
import retrospective as r
import nonscalar_cup_control as control

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/"SMALL_QUOTIENT_SOLUBILITY_PROTOCOL.json"
OUTPUT=r.OUT/"rank_jump_small_quotient_descents_v1.json"
WORK=r.ROOT/"artifacts/local/rank-jump-small-quotient-solubility-v1"


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),PROTOCOL,control.OUTPUT)}


def worker(index):
    from sage.all import pari
    from sage.version import version
    spec=r.read(PROTOCOL);model=spec["models"][index]
    E=pari.ellinit(model)
    result=pari.ellrank(E,spec["limits"]["PARI_ellrank_effort"],spec["supplied_points"][index])
    out={"bindings":bindings(),"index":index,"model":model,"status":"COMPUTED",
         "rank_interval":[int(result[0]),int(result[1])],"Sha2_mod_twice_Sha4_dimension":int(result[2]),
         "points":[list(map(str,P)) for P in result[3]],"raw_result":str(result),
         "software":{"sage":version,"pari":str(pari.version())}}
    r.write_new(WORK/f"curve-{index}.json",out)


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for index in range(2):
        path=WORK/f"curve-{index}.json"
        if not path.exists():
            with (WORK/f"curve-{index}.log").open("x") as log:
                try:
                    proc=subprocess.run(["sage","-python",str(Path(__file__).resolve()),"worker","--index",str(index)],
                                        cwd=r.ROOT,stdout=log,stderr=log,timeout=30)
                    error=None if proc.returncode==0 else "worker failed"
                except subprocess.TimeoutExpired:error="30-second timeout"
            if error:r.write_new(path,{"bindings":bindings(),"index":index,"status":"UNKNOWN","reason":error,
                                      "transcript":(WORK/f"curve-{index}.log").read_text()})
        row=r.read(path);assert row["bindings"]==bindings();rows.append(row)
        print("checkpoint",index,row.get("rank_interval"),row.get("Sha2_mod_twice_Sha4_dimension"),row.get("points"),flush=True)
    r.write_new(OUTPUT,{"schema":"rank-jump.small-quotient-descents.v1","bindings":bindings(),"records":rows})


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("mode",choices=["worker","capture"])
    parser.add_argument("--index",type=int);args=parser.parse_args()
    if args.mode=="worker":worker(args.index)
    else:capture()
