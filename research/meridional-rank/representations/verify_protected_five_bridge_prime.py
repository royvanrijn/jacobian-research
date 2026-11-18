#!/usr/bin/env python3
"""Literal braid/PD certificate, with optional rigorous Sage hyperbolicity replay.

Default mode uses only stdlib. --sage additionally uses pinned Sage 10.9,
SnapPy 3.3.2 and Spherogram 2.4.1, and reruns interval hyperbolicity.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from copy import deepcopy
import hashlib
import importlib
import json
from pathlib import Path
import platform
import random
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
INPUT = HERE / "protected_five_bridge_prime_inputs.json"
from verify_new_seventeen_crossing_knot import diagram, group_order
from verify_unbounded_prime_augmentation_gap_audit import (
    compose, cycle, inverse, perm_hurwitz, free_hurwitz, word_mul, word_inv)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def exact_certificate(data):
    prefix = [1,2,-3,-3,2,-3,-3,2,1,1,2,2,2,-3]
    word = prefix + [2,3,4,-3]*4
    require(data["literal_braid_word"] == word and data["braid_strands"] == 5,
            "literal thirty-letter protected braid")
    pd = data["pd"]
    require(len(pd) == len(word) == 30, "thirty presented crossings")
    signed, components, arc = diagram(pd)
    require([list(r) for r in signed] == data["signed_wirtinger_relations"],
            "oriented Wirtinger relation binding")
    require(components == data["arc_edge_components"], "Wirtinger arc binding")
    seed = [cycle(10,(0,1),(2,3)), cycle(10,(1,4),(2,3)),
            cycle(10,(0,5),(1,6)), cycle(10,(1,2),(6,7)), cycle(10,(0,8),(7,9))]
    require([list(p) for p in seed] == data["seed_permutations"], "literal A10 seed")
    require(perm_hurwitz(seed,word) == seed, "literal Hurwitz closure")
    old = [p[:8] for p in seed[:4]]
    require(group_order(old) == 20160, "first four generators give A8")
    require(seed[4] == cycle(10,(0,8),(7,9)), "two-new-symbol alternating extension")
    images = [tuple(p) for p in data["arc_permutations"]]
    require(len(images) == 30, "one image per Wirtinger arc")
    identity = tuple(range(10))
    for p in images:
        require(sorted(p) == list(identity) and compose(p,p) == identity
                and sum(p[i] != i for i in range(10)) == 4, "double transposition image")
    for o,a,b,e in signed:
        by = images[o] if e == 1 else inverse(images[o])
        require(compose(compose(by, images[a]), inverse(by)) == images[b],
                "literal signed Wirtinger relation")
    current=seed[:];wire=[None]*5;top=[None]*5;tags=list(range(5))
    seededges=[None]*5;edgecolors=defaultdict(set);strand=list(range(5))
    for letter,row in zip(word,pd):
        j=abs(letter)-1;a,b=current[j:j+2]
        if letter>0:
            c=compose(compose(b,a),b);ports=[a,b,c,b]
            incoming=[0,3];outgoing=[1,2];new=[b,c]
        else:
            c=compose(compose(a,b),a);ports=[b,a,c,a]
            incoming=[1,0];outgoing=[2,3];new=[c,a]
        for col,port in [(j,incoming[0]),(j+1,incoming[1])]:
            if wire[col] is None: wire[col]=row[port];top[col]=row[port]
            require(wire[col] == row[port], "literal braid incoming wiring")
            if tags[col] is not None: seededges[tags[col]]=row[port]
        wire[j:j+2]=[row[k] for k in outgoing];current[j:j+2]=new
        tags[j]=tags[j+1]=None;strand[j],strand[j+1]=strand[j+1],strand[j]
        for edge,value in zip(row,ports): edgecolors[edge].add(value)
    require(wire == top and current == seed, "all five closure identifications")
    require([arc[e] for e in seededges] == data["seed_arcs"], "top arc labeling")
    require(all(len(values) == 1 for values in edgecolors.values()), "consistent edge coloring")
    for j,block in enumerate(components):
        require({next(iter(edgecolors[e])) for e in block} == {images[j]}, "edge-to-arc coloring")
    visited=set();at=0
    while at not in visited: visited.add(at);at=strand[at]
    require(at == 0 and len(visited) == 5, "one-component five-strand closure")
    # Universal protected-tail identity, not a sample matrix identity.
    free=[(j,) for j in range(1,6)];C=free[3]
    def frame(v): return word_mul(v[1],v[2],C,v[4],word_inv(C))
    for generator in ([2],[3,4,-3]):
        after=free_hurwitz(free,generator)
        require(after[0] == free[0] and after[3] == C and frame(after) == frame(free),
                "protected outer entries and framed product")
    # Non-example: the original W5^2 family has a full rank-five reflection quotient.
    return {"literal_braid_and_PD_wiring":True,"spherical_one_component_PD":True,
            "presented_crossings":30,"full_A8_order":20160,
            "A10_extension_lemma_inputs":True,"full_A10_order_by_extension_lemma":1814400,
            "double_transposition_lower_bound":5,"closed_braid_bridge_upper":5,
            "protected_tail_free_group_identity":True,
            "augmentation_upper_via_paired_span_theorem":4,
            "stored_hyperbolicity_flag_is_not_used_as_a_proof":True}


def sage_certificate(data):
    # Match the retained project-local runtime without modifying the default env.
    for p in (ROOT/".scratch/mrc-sage-snappy",
              Path("/tmp/mrc-knot-audit-venv/lib/python3.12/site-packages")):
        if p.exists(): sys.path.append(str(p))
    import snappy, spherogram
    from sage.version import version as sage_version
    from sage.all import ComplexIntervalField, QQ
    from snappy.verify.hyperbolicity import check_logarithmic_gluing_equations_and_positively_oriented_tets
    g=data["geometry"]
    versions={"sage":sage_version,"snappy":snappy.__version__,
              "spherogram":spherogram.__version__,"python":sys.version.split()[0]}
    require(versions == g["runtime_versions"], "pinned analytic runtime versions")
    source_hashes={name:sha(importlib.import_module(name).__file__)
                   for name in g["runtime_analytic_source_hashes"]}
    require(source_hashes == g["runtime_analytic_source_hashes"], "analytic source hashes")
    binary=Path(snappy.__file__).parent/g["runtime_binary"]["file"]
    require(sha(binary) == g["runtime_binary"]["sha256"], "analytic SnapPy kernel hash")
    link=spherogram.Link(braid_closure=data["literal_braid_word"])
    require([list(row) for row in link.PD_code()] == data["pd"], "fresh literal braid gives retained PD")
    require(json.loads(json.dumps(link.DT_code())) == data["DT"], "fresh DT binding")
    random.seed(g["random_seed"]);snappy.set_rand_seed(g["random_seed"])
    fresh=spherogram.Link(data["pd"]).exterior()
    for _ in range(g["randomizations"]): fresh.randomize()
    require(fresh.triangulation_isosig() == g["decorated_triangulation_isosig"],
            "PD exterior triangulation binding")
    M=snappy.Manifold(g["triangulation_snappy"])
    require(M.triangulation_isosig() == g["decorated_triangulation_isosig"],
            "retained complete triangulation binding")
    require(M.num_tetrahedra() == g["tetrahedra"] == 30 and M.num_cusps() == g["cusps"] == 1,
            "thirty-tetrahedron one-cusp triangulation")
    require(json.loads(json.dumps(M._get_tetrahedra_gluing_data())) == g["face_pairings"],
            "literal face pairing data")
    require(Counter(tuple(map(int,row)) for row in M.gluing_equations()) == Counter(map(tuple,g["log_gluing_equations"])),
            "retained logarithmic gluing equations")
    CIF=ComplexIntervalField(g["bits_prec"]);RIF=CIF.real_field()
    retained=[]
    for a,b,c,d in g["shape_intervals_rational"]:
        retained.append(CIF(RIF(QQ(a),QQ(b)),RIF(QQ(c),QQ(d))))
    check_logarithmic_gluing_equations_and_positively_oriented_tets(M,retained)
    okay,shapes=M.verify_hyperbolicity(bits_prec=g["bits_prec"])
    require(okay, "fresh Krawczyk interval existence verification")
    check_logarithmic_gluing_equations_and_positively_oriented_tets(M,shapes)
    # Replay on the exterior constructed from the PD as a second topological binding.
    fresh_okay,_=fresh.verify_hyperbolicity(bits_prec=g["bits_prec"])
    require(fresh_okay, "fresh PD-exterior interval verification")
    def endpoints(x): return [str(x.lower().exact_rational()),str(x.upper().exact_rational())]
    shape_intervals=[endpoints(z.real())+endpoints(z.imag()) for z in shapes]
    return {"runtime_versions":versions,"analytic_source_hashes":source_hashes,
            "kernel_sha256":sha(binary),"literal_braid_PD_exterior_binding":True,
            "retained_interval_logarithmic_and_positive_orientation_checks":True,
            "fresh_Krawczyk_hyperbolicity":True,"fresh_PD_exterior_hyperbolicity":True,
            "tetrahedra":30,"cusps":1,"bits_prec":g["bits_prec"],
            "decorated_triangulation_isosig":M.triangulation_isosig(),
            "fresh_shape_intervals_rational":shape_intervals,
            "prime_via_hyperbolic_knot_theorem":True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sage",action="store_true")
    parser.add_argument("--write-receipt",type=Path)
    args=parser.parse_args()
    data=json.loads(INPUT.read_text())
    report={"schema_version":1,"python":platform.python_version(),
            "checker_sha256":sha(__file__),"input_sha256":sha(INPUT),
            "dependency_sha256":{name:sha(HERE/name) for name in (
                "verify_new_seventeen_crossing_knot.py",
                "verify_unbounded_prime_augmentation_gap_audit.py")},
            "finite_certificate":exact_certificate(data)}
    rejected=0
    for field in ("literal_braid_word","seed_permutations"):
        bad=deepcopy(data)
        if field == "literal_braid_word": bad[field][0] *= -1
        else: bad[field][4] = list(range(10))
        try: exact_certificate(bad)
        except AssertionError: rejected += 1
        else: raise AssertionError("corrupted input accepted: "+field)
    report["negative_controls_rejected"] = rejected
    if args.sage: report["rigorous_geometry_certificate"]=sage_certificate(data)
    report["result"]="PASS"
    if args.write_receipt: args.write_receipt.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({k:v for k,v in report.items() if k != "rigorous_geometry_certificate"},indent=2))
    if args.sage: print("Rigorous geometry: PASS (fresh Krawczyk and fresh PD exterior replay).")

if __name__ == "__main__": main()