#!/usr/bin/env python3
"""Independent exact no-intersection proof for the inherited node probes.

Recompute x(2Q)-x(T) directly from the generic coordinates. A retained finite
projective no-root witness for its numerator excludes every rational finite
node position for this word pair. This bypasses Sage's T-2Q factorizations.
"""
import argparse
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/"artifacts/generated-results/elkies-k3-r17-inherited-halving-nodes-v1"
SPEC=importlib.util.spec_from_file_location("node_fraction_replay",ROOT/"elkies-k3/scripts/verify_r17_one_node_correlated.py")
N=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(N)
V=N.V;require=V.require


def intersection_polynomial(trace,point,A):
    h,Nx=[V.poly(trace[k]) for k in ("h","Nx")]
    x,y=[V.poly(point[k]) for k in ("x","y")]
    tangent_numerator=V.add(V.scale(V.power(x,2),3),A)
    denominator=V.scale(y,2)
    x2_numerator=V.sub(V.power(tangent_numerator,2),V.scale(V.mul(x,V.power(denominator,2)),2))
    return V.sub(V.mul(x2_numerator,V.power(h,2)),V.mul(Nx,V.power(denominator,2)))


def verify(path,export=False):
    packet=json.loads((path/"input.json").read_text())
    result=json.loads((path/"result.json").read_text())
    require(packet["schema"]=="r17-inherited-halving-node-input-v1","input schema")
    require(result["input_sha256"]==sha256((path/"input.json").read_bytes()).hexdigest(),"input binding")
    prior_path=N.DEFAULT/"input.json"
    prior=json.loads(prior_path.read_text())
    require(packet["source_input_sha256"]==sha256(prior_path.read_bytes()).hexdigest(),"source input binding")
    require(packet["generic_source"]==prior["generic_source"] and packet["traces"]==prior["traces"],"frozen generic fields")
    source=packet["generic_source"];V.source_projection(source)
    require(packet["probe_words"]==[{"index":i,"sign":s} for i in range(17) for s in (1,-1)],"probe bank")
    A=V.poly(source["A"])
    witnesses={} if export else json.loads((path/"root-witnesses.json").read_text())["witnesses"]
    used=set()
    for index,trace in enumerate(packet["traces"]):
        for j,point in enumerate(source["basis"]):
            F=intersection_polynomial(trace,point,A)
            key=f"{index}:{j}"
            if export:witnesses[key]=N.root_witness(F)
            require(key in witnesses and not N.has_projective_root(N.primitive(F),witnesses[key]),"generic intersection root exclusion")
            used.add(key)
    require(set(witnesses)==used and len(used)==1139,"root witness coverage")
    records=N.read_traces(path)
    for k,row in enumerate(records):
        require(row["trace_index"]==k and [r["probe"] for r in row["records"]]==packet["probe_words"],"checkpoint coverage")
        require(row["carriers"]==[] and all(r["rational_roots"]==[] for r in row["records"]),"checkpoint conclusion")
    require(result["pairs"]==2278 and result["rational_root_incidence"]==result["distinct_carriers"]==0 and result["genus_one_carriers"]==[],"result reconciliation")
    if export:
        with (path/"root-witnesses.json").open("x") as f:
            json.dump({"schema":"r17-inherited-node-root-witnesses-v1","witnesses":witnesses},f,indent=2,sort_keys=True);f.write("\n")
    return {"status":"PASS_NO_FINITE_RATIONAL_INHERITED_NODE_IN_FIXED_WORD_BANK",
            "traces":67,"signed_generic_probes":34,"word_pairs":2278,
            "independent_no_root_polynomials":1139,"rational_fibre_bound":None,
            "scope":"No finite rational smooth fibre has T=2Q for these 67 traces and 34 signed published basis sections. This does not exclude halves outside this fixed generic word bank, or other covers."}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--input-dir",type=Path,default=DEFAULT);p.add_argument("--export",action="store_true");p.add_argument("--output",type=Path)
    a=p.parse_args();started=time.monotonic();result=verify(a.input_dir,a.export);result["elapsed_seconds"]=round(time.monotonic()-started,3)
    if a.output:
        with a.output.open("x") as f:json.dump(result,f,indent=2,sort_keys=True);f.write("\n")
    print(json.dumps(result,sort_keys=True))


if __name__=="__main__":main()
