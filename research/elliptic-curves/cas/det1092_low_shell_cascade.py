#!/usr/bin/env python3
"""A target-blind low-shell adaptive cascade on a frozen determinant-1092 cohort."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import det1092_record_scale_selection as selection
import low_shell_adaptive_engine as engine
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint, digest
from research_runtime.supervisor import Limits, run


ROOT, CAS, ART, SAGE = selection.ROOT, selection.CAS, selection.ART, selection.SAGE
cert = selection.cert
BATCH = ROOT / "artifacts/local/elliptic-curves/det1092-low-shell-cascade-v3"
READS: set[str] = set()


def install_guard() -> None:
    allowed = (BATCH.resolve(), selection.D.resolve(), selection.SOURCE.resolve())
    def guard(event, args):
        if event != "open" or not args or not isinstance(args[0], (str, bytes)): return
        path = Path(args[0]).resolve()
        if not path.is_relative_to(ROOT / "artifacts"): return
        if path.is_relative_to(allowed[0]) or path.is_relative_to(allowed[1]) or path == allowed[2]:
            READS.add(str(path.relative_to(ROOT))); return
        raise PermissionError("low-shell execution input rejected: " + str(path.relative_to(ROOT)))
    sys.addaudithook(guard)


def sources() -> dict[str,str]:
    names=["det1092_low_shell_cascade.py","low_shell_adaptive_engine.py","prepare_det1092_low_shell_seed.sage","prepare_det1092_low_shell_maps.sage","certify_low_shell_factor_free_exposure.py","verify_low_shell_factor_free_exposure.sage","verify_factor_free_rank.sage"]
    return {**engine.sources(), **{str((CAS/name).relative_to(ROOT)):cert.hashed(CAS/name) for name in names}}


def frozen_rows() -> list[dict]:
    result=cert.read(selection.D / "selection-result.json")
    if result["status"] != "PASS" or len(result["selected"]) != 48: raise ArithmeticError("the score-only determinant-1092 selection changed")
    rows=[]
    for stratum in ("strong","moderate","lower_fixed"):
        choices=[row for row in result["selected"] if row["stratum"]==stratum]
        if not choices: raise ArithmeticError("missing frozen score stratum")
        # Score and ID are pre-point-exposure data.  One row per stratum avoids
        # using the completed deep-centre exposure as an outcome filter.
        rows.append(min(choices,key=lambda row:(-row["score_units"],row["id"])))
    return rows


def freeze() -> None:
    if BATCH.exists(): raise FileExistsError("preserve frozen low-shell cascade")
    rows=frozen_rows(); BATCH.mkdir(parents=True)
    inputs=[selection.D / "protocol.json",selection.D / "selection-result.json",selection.SOURCE]
    checkpoint(BATCH / "protocol.json",{"schema":"elliptic-curves.det1092-low-shell-cascade.v1","sources":sources(),"inputs":{str(path.relative_to(ROOT)):cert.hashed(path) for path in inputs},"rows":rows,"maximum_curves":3,"maximum_waves_per_curve":5,"charts":49,"sample_size":2048,"sample_domain":"det1092-intrinsic-low-shell-v1","height":125000,"seconds_per_chart":10,"target_rank":32,"rank_stop":False,"rss_bytes":2147483648,"seed_seconds":180,"geometry_wall_seconds":180,"worker_wall_seconds":1200,"replay_wall_seconds":1200,"gp_sha256":cert.hashed(selection.scalar.GP),"selection":"Exactly one maximum frozen score row in each of the strong, moderate, and lower-fixed strata. The source contains no catalogue rank, exceptional-point coordinate, or prior point outcome. No outcome-dependent replacement.","adaptive_policy":"At every wave recompute 2048 SHA-derived nonzero parities in the current certified subgroup, retain the 49 smallest exact rounded canonical-height coset representatives (ties by parity), map all 49 before point execution, and continue only after an independently certified rank gain. Later waves require a parity nonzero in a newly added coordinate. Stop at no gain, rank at least32, or five waves.","basis_invariance":"The score is an intrinsic coset norm in the rounded canonical-height lattice. The geometry stage records a unimodular LLL transport and rechecks each original-coordinate parity/representative identity; a separate basis-randomization audit is required before any invariance claim.","execution_blindness":{"allowed_artifact_inputs":[str(path.relative_to(ROOT)) for path in inputs],"forbidden_inputs":["known exceptional points","catalogue ranks","prior determinant-1092 point-exposure outputs","curve-302 residual diagnostics"]},"boundary":"Prospective bounded point search. Complete misses are retained; no miss is an upper-bound statement and rank at least32 needs 32 independently certified points."})
    print("FROZEN DET1092 LOW-SHELL CASCADE|curves=3|waves<=5|target=32",flush=True)


def campaign() -> dict:
    p=cert.read(BATCH / "protocol.json")
    if p["sources"] != sources() or any(cert.hashed(ROOT/path)!=value for path,value in p["inputs"].items()): raise ArithmeticError("frozen low-shell protocol changed")
    return p


def configure(index:int,wave:int) -> None:
    global D, SEED, ROW
    candidate=campaign()["rows"][index]; D=BATCH/candidate["id"]/("wave-"+str(wave+1).zfill(2)); SEED=D/"seed.json"; ROW=cert.read(D/"protocol.json")["rows"][0]; engine.bind(sys.modules[__name__])


def protocol() -> dict:
    p=cert.read(D/"protocol.json"); top=campaign()
    if p["campaign_sha256"]!=cert.hashed(BATCH/"protocol.json") or p["seed_sha256"]!=cert.hashed(SEED) or p["sources"]!=top["sources"]: raise ArithmeticError("wave protocol changed")
    return p


def masks(p:dict) -> list[int]:
    values=[]; counter=0
    while len(values)<p["sample_size"]:
        mask=int(digest([p["sample_domain"],counter]),16)%(1<<ROW["initial_rank"]); counter+=1
        if mask>>ROW["mask_floor"] and mask not in values: values.append(mask)
    return values


def _stage(ledger,entry,name,command,seconds,directory):
    status=run(command,limits=Limits(seconds,ledger["rss_bytes"]),log_path=directory/(name+".log"),checkpoint_path=directory/(name+".supervisor.json"),cwd=ROOT)
    ok=status["outcome"]=="completed" and status["returncode"]==0
    entry["stages"].append({"name":name,"status":"PASS" if ok else "FAILED_OR_CENSORED","supervision":status}); checkpoint(BATCH/"ledger.json",ledger)
    if not ok: ledger["status"]=entry["status"]="FAILED_OR_CENSORED"; checkpoint(BATCH/"ledger.json",ledger); raise ArithmeticError("preserve failed low-shell stage")


def launch() -> None:
    p=campaign()
    if (BATCH/"ledger.json").exists(): raise FileExistsError("preserve low-shell cascade ledger")
    ledger={"status":"RUNNING","rss_bytes":p["rss_bytes"],"rows":[],"completed_boxes":0};checkpoint(BATCH/"ledger.json",ledger)
    for index,candidate in enumerate(p["rows"]):
        folder=BATCH/candidate["id"];folder.mkdir(); entry={"id":candidate["id"],"stratum":candidate["stratum"],"status":"RUNNING","waves":[],"stages":[]};ledger["rows"].append(entry);checkpoint(BATCH/"ledger.json",ledger)
        _stage(ledger,entry,"initial-seed",[SAGE,str(CAS/"prepare_det1092_low_shell_seed.sage"),"--index",str(index)],p["seed_seconds"],folder)
        seed=cert.read(folder/"initial-seed.json");original=seed;floor=0
        for wave in range(p["maximum_waves_per_curve"]):
            directory=folder/("wave-"+str(wave+1).zfill(2));directory.mkdir();rank=len(seed["points"]);proof=seed["rank_certificate"]
            actual=checked_rank(tuple(map(cert.F,seed["curve"])),[tuple(map(cert.F,point)) for point in seed["points"]],[row["prime"] for row in proof["signatures"]],proof["no_rational_2_torsion_prime"])
            if digest(actual)!=digest(proof) or seed["points"][:17]!=original["points"] or rank<=floor: raise ArithmeticError("certified cascade seed changed")
            checkpoint(directory/"seed.json",seed);checkpoint(directory/"seed-cloud.json",{"curve":seed["curve"],"points":seed["points"],"signatures":proof["signatures"],"rank_certificate":proof,"rank_lower_bound":rank,"independent_column_indices":list(range(rank))})
            checkpoint(directory/"protocol.json",dict(p,campaign_sha256=cert.hashed(BATCH/"protocol.json"),seed_sha256=cert.hashed(directory/"seed.json"),candidate_id=candidate["id"],wave=wave+1,rows=[{"id":directory.name,"initial_rank":rank,"generic_dimension":17,"mask_floor":floor}]))
            configure(index,wave);wave_entry={"id":directory.name,"initial_rank":rank,"mask_floor":floor,"status":"RUNNING","stages":[]};entry["waves"].append(wave_entry);checkpoint(BATCH/"ledger.json",ledger)
            _stage(ledger,wave_entry,"geometry",[SAGE,str(CAS/"prepare_det1092_low_shell_maps.sage"),"--index",str(index),"--wave",str(wave)],p["geometry_wall_seconds"],directory)
            for stage in ("worker","replay"):_stage(ledger,wave_entry,stage,[sys.executable,str(Path(__file__).resolve()),stage,"--index",str(index),"--wave",str(wave)],p[stage+"_wall_seconds"],directory)
            prefix="det1092_low_shell_v3_"+candidate["id"].replace("-","_")+"_"+directory.name.replace("-","_")
            _stage(ledger,wave_entry,"certificates",[sys.executable,str(CAS/"certify_low_shell_factor_free_exposure.py"),"--run",str(directory),"--protocol",str(directory/"protocol.json"),"--prefix",prefix],1800,directory)
            certified=cert.read(directory/"certification-ledger.json"); cloud=cert.read(ART/(prefix+"_mod2_v1.json"))
            if certified["status"]!="PASS" or any(value!=certified["rank_lower_bound"] for value in certified["odd_modulus_ranks"].values()) or cloud["independent_points"][:rank]!=seed["points"]: raise ArithmeticError("low-shell certificate differs")
            wave_entry.update(status="PASS",rank_lower_bound=certified["rank_lower_bound"],completed_boxes=49,cloud_path=str((ART/(prefix+"_mod2_v1.json")).relative_to(ROOT)));ledger["completed_boxes"]+=49;entry["rank_lower_bound"]=certified["rank_lower_bound"];checkpoint(BATCH/"ledger.json",ledger)
            if certified["rank_lower_bound"]==rank:entry["stop_reason"]="NO_CERTIFIED_GAIN";break
            if certified["rank_lower_bound"]>=32:entry["stop_reason"]="CERTIFIED_AT_LEAST32";break
            floor=rank;seed=dict(original,points=cloud["independent_points"],rank_certificate=cloud["rank_certificate"])
        else:entry["stop_reason"]="FIVE_WAVE_LIMIT"
        entry["status"]="PASS";checkpoint(BATCH/"ledger.json",ledger)
    ledger.pop("rss_bytes");ledger["status"]="PASS";checkpoint(BATCH/"ledger.json",ledger);print("LOW-SHELL DET1092 CASCADE COMPLETE",flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("stage",choices=("freeze","launch","worker","replay"));parser.add_argument("--index",type=int);parser.add_argument("--wave",type=int);args=parser.parse_args()
    if args.stage in ("worker","replay"):
        install_guard();configure(args.index,args.wave);getattr(engine,args.stage)();checkpoint(D/(args.stage+"-data-access.json"),sorted(READS))
    else:getattr(sys.modules[__name__],args.stage)()
