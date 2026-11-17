#!/usr/bin/env python3
"""Known P(2,3,3,5) control:23 A5 kernels improve the raw Hom bound."""

from collections import Counter
import hashlib
from itertools import product
import json
from math import gcd
from pathlib import Path
import time

from joint_a5 import group_model,core_extension,generated,certify_distinct_same_source
from parabolic_lifts import oriented_diagram
from search_lift_panel import seed_certificate
from verify_diagram_audit import Tangle,plane_map,walk

HERE=Path(__file__).resolve().parent


def main():
    integers=[2,3,3,5]
    tangle=Tangle.integer(0)
    for a in integers:
        tangle=tangle.add(Tangle.rational(1,a))
    raw_pd=tangle.closure_pd()
    visits=walk(raw_pd)
    edge_labels={raw_pd[d//4][d%4]:i+1 for i,d in enumerate(visits)}
    rotations={d//4:d%4 for d in visits if d%2==0}
    assert len(edge_labels)==26 and len(rotations)==13
    pd=[[edge_labels[row[(rotations[c]+h)%4]] for h in range(4)] for c,row in enumerate(raw_pd)]
    plane_map(pd)
    assert len(pd)==13 and len(walk(pd))==26
    crossings=[c[:3] for c in oriented_diagram(pd)]
    seeds,trace,trials=seed_certificate(13,crossings,4,715,time.monotonic()+10)
    assert seeds is not None
    catalogue=json.loads((HERE/"joint_a5_catalogue.json").read_text())
    elements,table,inverses,orders,automorphisms,core=group_model()
    source="P(2,3,3,5):double-cover:arc"+str(seeds[0])
    representatives=[]
    for index,values in enumerate(catalogue["generating_triple_orbits"]):
        labels=core_extension(13,crossings,seeds,trace,values,core)
        if labels is not None:
            representatives.append({"source_group":source,"orbit_index":index,"generator_images":values,"arc_images":labels})
    assert len(representatives)==23
    assert certify_distinct_same_source(representatives[:20],source,table,automorphisms)==20
    raw_count=sum(core_extension(13,crossings,seeds,trace,values,core) is not None for values in product(range(60),repeat=3))
    assert raw_count==2901
    # Independently specified Seifert fibre quotient, rather than a diagram.
    # Every nonidentity possible fibre image is ruled out by one x_i^a_i=h^-1.
    for h in range(1,60):
        exponent=orders[h]
        assert exponent in integers
        assert not any(orders[x]//gcd(orders[x],exponent)==orders[h] for x in range(60))
    orbifold_hom=0
    image_orders=Counter()
    allowed=[[i for i in range(60) if exponent%orders[i]==0] for exponent in integers]
    for a,b,c in product(*allowed[:3]):
        d=inverses[table[table[a][b]][c]]
        if d in allowed[3]:
            orbifold_hom+=1
            image_orders[len(generated((a,b,c),table))]+=1
    assert orbifold_hom==raw_count and image_orders=={1:1,3:20,12:120,60:2760}
    inputs=["build_joint_a5_pretzel_control.py","joint_a5.py","joint_a5_catalogue.json",
            "verify_diagram_audit.py","parabolic_lifts.py","search_lift_panel.py"]
    result={"status":"PASS","knot":"P(2,3,3,5)","definition":"Numerator closure of the sum1/2+1/3+1/3+1/5",
            "pd":pd,"raw_tangle_pd":raw_pd,"edge_labels":sorted(edge_labels.items()),"crossing_rotations":rotations,
            "seeds":seeds,"seed_trace":trace,"normalized_arc":seeds[0],"source_group":source,
            "A5_quotient_kernels":23,"epimorphism_representatives":representatives,
            "Hom_to_A5":raw_count,"Seifert_orbifold_Hom_count":orbifold_hom,"Hom_image_orders":dict(image_orders),
            "raw_Hom_meridional_lower_bound":3,"joint_kernel_meridional_lower_bound":4,"bridge_upper_bound":4,
            "known_family":"Four-noninteger-tangle Montesinos, already covered by Boileau-Zieschang1985; no novel knot case",
            "input_sha256":{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs}}
    (HERE/"joint_a5_pretzel_control.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":"PASS","knot":result["knot"],"kernels":23,"Hom":raw_count,
                      "raw_Hom_lower":3,"joint_kernel_lower":4,"bridge_upper":4,"image_orders":dict(image_orders)},indent=2))


if __name__=="__main__":
    main()
