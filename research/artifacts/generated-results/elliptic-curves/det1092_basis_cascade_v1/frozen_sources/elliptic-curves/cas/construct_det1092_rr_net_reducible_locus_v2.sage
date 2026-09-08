#!/usr/bin/env sage-python
"""Correct the v1 claim boundary; retain its exact enumeration unchanged.

A reducible member C_min+F_tau can inherit useful points from C_min at
other fibres. Only the22 section pairs are everywhere old-subgroup points;
the already nonsplit C_min contributes no rational point at302.
"""
import json,runpy,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves'
if __name__=='__main__':
    original=runpy.run_path(str(CAS/'construct_det1092_rr_net_reducible_locus.sage'))
    d=original['build']()
    d['conclusion']='The22 section-pair members supply only known-subgroup specializations. The remaining reducible pencil C_min+F_tau repeats the prior bisection, nonsplit at302. The vertical component at tau=0 is the entire302 fibre and does not construct a point; no exclusion at other fibres is claimed.'
    d['v1_boundary_rejection']='The v1 conclusion incorrectly excluded useful specializations on every smooth fibre. Only this conclusion is rejected; the exact enumeration and equations are unchanged.'
    old=ART/'det1092_rr_net_reducible_locus_v1.json'
    d['inputs'][str(old.relative_to(ROOT))]=hashlib.sha256(old.read_bytes()).hexdigest()
    d['inputs'][str(Path(__file__).relative_to(ROOT))]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out=ART/'det1092_rr_net_reducible_locus_v2.json'
    if out.exists():raise FileExistsError(out)
    out.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(d['status'],'corrected302-specific scope',flush=True)
