#!/usr/bin/env python3
"""Exact actual-double-cover topology of24 A6 sharp-four knots.

Requires SnapPy3.3.2/Spherogram2.4.1/Regina7.4 and pinned named census.
No numerical volume or numerical isometry is an acceptance condition.
"""
from copy import deepcopy
from importlib.metadata import version
from pathlib import Path
import hashlib,json,random,warnings
warnings.filterwarnings("ignore",message="Plink failed to import tkinter.*")
import snappy,spherogram,regina,snappy_15_knots

HERE=Path(__file__).resolve().parent
INPUT=HERE/"a6_cover_topology_inputs.json"
EXPECTED=["13a_304","13a_362","13a_623","13a_629","13a_3232","13a_3259","13a_3374","13a_3375",
          "13n_834","13n_838","13n_839","13n_991","13n_992","13n_993","13n_1367","13n_1370",
          "13n_1371","13n_1374","13n_1375","13n_1376","13n_3814","13n_3840","13n_3911","13n_3914"]
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def edge_pairs(pd):
    ends={}
    for i,row in enumerate(pd):
        for h,e in enumerate(row):ends.setdefault(e,[]).append(4*i+h)
    assert all(len(v)==2 for v in ends.values())
    return {a:b for pair in ends.values() for a,b in [pair,pair[::-1]]}
def named_binding(row):
    link=spherogram.Link("K"+row["name"].replace("_",""))
    actual=[list(r) for r in link.PD_code()]
    assert actual==row["census_pd"] and json.loads(json.dumps(link.DT_code()))==row["census_DT"]
    a=edge_pairs(row["pd"]);b=edge_pairs(actual)
    witness=row["census_map_witness"];f=witness["dart_images"]
    epsilon,swap=witness["orientation"],witness["over_under_swap"]
    assert epsilon in [-1,1] and swap in [0,1]
    assert len(row["pd"])==len(actual)==13 and sorted(f)==list(range(52))
    for x,y in enumerate(f):
        assert f[a[x]]==b[y]
        assert (x%2)^(y%2)==swap
        assert f[4*(x//4)+(x+1)%4]==4*(y//4)+(y+epsilon)%4
def census_table():
    table={}
    for line in (HERE/"sources/burton-hypcensus-5-or.dat").read_text().splitlines():
        fields=line.split()
        if fields:
            assert int(fields[5])==len(fields[6:])
            for sig in fields[6:]:table[sig]=fields[0]
    return table
def close(row):
    seed=row["random_seed"];random.seed(seed);snappy.set_rand_seed(seed);regina.RandomEngine.reseedWithDefault()
    E=spherogram.Link(row["pd"]).exterior()
    assert E.num_cusps()==1 and E.triangulation_isosig()==row["exterior_isosig_with_peripheral"]
    E.dehn_fill((2,0));covers=E.covers(2);assert len(covers)==1
    C=covers[0];assert C.cover_info()["degree"]==2 and C.cover_info()["type"]=="cyclic"
    assert C.num_cusps()==1 and C.cusp_info("filling")==[(1.0,0.0)]
    T=regina.Triangulation3(C.filled_triangulation()._to_string());T.simplify()
    assert T.isoSig()==row["closed_cover_isosig"] and T.size()==row["closed_tetrahedra"]
    assert T.isValid() and T.isClosed() and T.isConnected() and T.isOrientable()
    return T
def check_coordinates(surface,row):
    vector=[int(surface.quads(j,k).longValue()) for j in range(row["closed_tetrahedra"]) for k in range(3)]
    assert vector==row["normal_quad_coordinates"]
def check_piece_name(row,ideal_sig,table):
    assert ideal_sig==row["ideal_piece_isosig"]
    assert table[ideal_sig]==row["hyperbolic_census_name"]
def check_case(row,table):
    named_binding(row)
    T=close(row);assert T.isIrreducible()
    surfaces=regina.NormalSurfaces(T,regina.NormalCoords.Quad,regina.NormalList.Vertex)
    surface=surfaces[row["normal_surface_index"]];check_coordinates(surface,row)
    assert surface.isConnected() and surface.isOrientable() and surface.isTwoSided() and surface.eulerChar()==0
    assert surface.isIncompressible()
    regina.RandomEngine.reseedWithDefault()
    cut=surface.cutAlong();cut.simplify();parts=cut.triangulateComponents()
    assert len(parts)==2
    matches=[]
    for U in parts:
        U.simplify()
        assert U.isValid() and U.isOrientable() and U.countBoundaryComponents()==1
        assert U.boundaryComponent(0).eulerChar()==0
        V=regina.Triangulation3(U);assert V.finiteToIdeal();V.simplify()
        if U.isoSig()==row["finite_piece_isosig"]:
            check_piece_name(row,V.isoSig(),table);matches.append(V.isoSig())
    assert matches==[row["ideal_piece_isosig"]]
    return {"name":row["name"],"closed_cover_isosig":T.isoSig(),"closed_tetrahedra":T.size(),
            "irreducible":True,"incompressible_separating_normal_torus":True,
            "hyperbolic_piece_census_name":row["hyperbolic_census_name"],"hyperbolic_piece_isosig":matches[0],
            "outside_BJW18_graph_cover_hypothesis":True}
def rejected(call):
    try:call()
    except (AssertionError,IndexError,KeyError,ValueError,RuntimeError):return True
    raise AssertionError("corrupted topology assertion accepted")
def main():
    data=json.loads(INPUT.read_text());assert [r["name"] for r in data["cases"]]==EXPECTED
    versions={p:version(p) for p in ["snappy","spherogram","snappy_15_knots","snappy_manifolds"]}
    versions["regina"]=regina.versionString()
    assert versions=={"snappy":"3.3.2","spherogram":"2.4.1","snappy_15_knots":"1.2.1","snappy_manifolds":"1.4","regina":"7.4"}
    for filename,pin in data["source_sha256"].items():assert digest(HERE/filename)==pin
    provenance=data["census_provenance"];package=Path(snappy_15_knots.__file__).parent
    assert provenance["version"]==versions["snappy_15_knots"] and provenance["spherogram_version"]==versions["spherogram"]
    for pin in provenance["source_files"]:
        raw=(package/pin["path"]).read_bytes()
        assert len(raw)==pin["bytes"] and hashlib.sha256(raw).hexdigest()==pin["sha256"]
    panel={r["name"]:r["pd"] for r in json.loads((HERE/"joint_a5_panel.json").read_text())["cases"]}
    assert all(r["pd"]==panel[r["name"]] for r in data["cases"])
    table=census_table();rows=[]
    for row in data["cases"]:
        rows.append(check_case(row,table))
        print("PASS",row["name"],row["hyperbolic_census_name"],flush=True)
    bad=deepcopy(data["cases"][0]);bad["name"]="13a_373"
    bad2=deepcopy(data["cases"][0]);bad2["hyperbolic_census_name"]="m004"
    bad3=deepcopy(data["cases"][0]);bad3["normal_quad_coordinates"][0]+=1
    T=close(data["cases"][0]);assert T.isIrreducible()
    surfaces=regina.NormalSurfaces(T,regina.NormalCoords.Quad,regina.NormalList.Vertex)
    controls={"wrong_named_knot":rejected(lambda:named_binding(bad)),
              "wrong_hyperbolic_census_name":rejected(lambda:check_piece_name(bad2,bad2["ideal_piece_isosig"],table)),
              "changed_normal_quad_vector":rejected(lambda:check_coordinates(surfaces[bad3["normal_surface_index"]],bad3))}
    files=[INPUT,Path(__file__)]+[HERE/p for p in data["source_sha256"]]
    receipt={"status":"PASS","versions":versions,"scope":data["scope"],"cases":rows,"negative_controls":controls,
             "input_sha256":{p.relative_to(HERE).as_posix():digest(p) for p in files},
             "assurance":"Exact internal actual-cover construction, current named-census binding, irreducibility and incompressible-normal-torus replay. Burton certified census hyperbolicity and trusted SnapPy/Regina algorithms are imported; no numerical volume/isometry, formal verification or external specialist review."}
    (HERE/"a6_cover_topology_verification.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"status":"PASS","outside_BJW18":len(rows),"m143":sum(r["hyperbolic_piece_census_name"]=="m143" for r in rows),
                      "m239":sum(r["hyperbolic_piece_census_name"]=="m239" for r in rows),"negative_controls":controls},indent=2))
if __name__=="__main__":main()
