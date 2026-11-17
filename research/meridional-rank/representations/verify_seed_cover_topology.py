#!/usr/bin/env python3
"""Exact finite topology checks; no floating-point geometric acceptance.

Replay: topology Python with SnapPy 3.3.2, Spherogram 2.4.1, Regina 7.4.
Hyperbolicity is imported from Burton's certified census, Corollary 1.2,
and its exact author data. This does not replay the census-wide proof.
"""
from copy import deepcopy
from importlib.metadata import version
from pathlib import Path
import hashlib
import json
import random
import warnings

warnings.filterwarnings("ignore", message="Plink failed to import tkinter.*")
import snappy
import spherogram
import regina
from audit_diagram_identity import normalize, pd_walk

HERE = Path(__file__).resolve().parent
INPUT = HERE / "seed_cover_topology_inputs.json"
CENSUS_FILES = ["burton-hypcensus-4-or.dat", "burton-hypcensus-5-or.dat"]
SOURCE_PDF = "burton-cusped-census-1405.2695v1.pdf"
EXPECTED = ['13a_3145', '13a_3147', '13a_3148', '13a_3149', '13a_373', '13a_375', '13a_391', '13a_392', '13a_3264', '13a_3266', '13a_3275', '13a_3276', '13n_1005', '13n_1006', '13n_1007', '13n_1011', '13n_1012', '13n_1013', '13n_1043', '13n_1044', '13n_1045', '13n_1046', '13n_1047', '13n_1048', '13n_3842', '13n_3844', '13n_3855', '13n_3856']


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def census_table():
    result = {}
    for filename in CENSUS_FILES:
        for line in (HERE / "sources" / filename).read_text().splitlines():
            fields = line.split()
            if not fields:
                continue
            assert int(fields[5]) == len(fields[6:])
            for sig in fields[6:]:
                assert sig not in result
                assert regina.Triangulation3.fromIsoSig(sig).isValid()
                result[sig] = fields[0]
    return result


def named_diagram_matches(row):
    lookup = "K" + row["name"].replace("_", "")
    census_pd = [list(c) for c in spherogram.Link(lookup).PD_code()]
    assert normalize(pd_walk(row["pd"]), True)[0] == normalize(pd_walk(census_pd), True)[0]
    return lookup


def claimed_piece_matches(row, ideal_sig, census):
    assert census.get(ideal_sig) == row["hyperbolic_census_name"]
    assert ideal_sig == row["ideal_piece_isosig"]


def check_case(row, census):
    lookup = named_diagram_matches(row)
    seed = row["random_seed"]
    random.seed(seed)
    snappy.set_rand_seed(seed)
    regina.RandomEngine.reseedWithDefault()
    exterior = spherogram.Link(row["pd"]).exterior()
    assert exterior.num_cusps() == 1
    assert exterior.triangulation_isosig() == row["exterior_isosig_with_peripheral"]
    # Meridians supplied by the diagram exterior; the unique twofold cover of
    # the order-two orbifold filling lifts the branching slope primitively.
    exterior.dehn_fill((2, 0))
    covers = exterior.covers(2)
    assert len(covers) == 1
    cover = covers[0]
    assert cover.cover_info()["degree"] == 2
    assert cover.cover_info()["type"] == "cyclic"
    assert cover.num_cusps() == 1
    assert cover.cusp_info("filling") == [(1.0, 0.0)]
    filled = cover.filled_triangulation()
    tri = regina.Triangulation3(filled._to_string())
    tri.simplify()
    assert tri.isoSig() == row["closed_cover_isosig"]
    assert tri.size() == row["closed_tetrahedra"]
    assert tri.isValid() and tri.isClosed() and tri.isOrientable() and tri.isConnected()
    assert tri.isIrreducible()
    surfaces = regina.NormalSurfaces(tri, regina.NormalCoords.Quad, regina.NormalList.Vertex)
    surface = surfaces[row["normal_surface_index"]]
    vector = [int(surface.quads(t, q).longValue()) for t in range(tri.size()) for q in range(3)]
    assert vector == row["normal_quad_coordinates"]
    assert surface.isConnected() and surface.isOrientable() and surface.isTwoSided()
    assert surface.eulerChar() == 0
    assert surface.isIncompressible()
    regina.RandomEngine.reseedWithDefault()
    cut = surface.cutAlong()
    cut.simplify()
    pieces = cut.triangulateComponents()
    assert len(pieces) == 2
    found = []
    for piece in pieces:
        piece.simplify()
        assert piece.isValid() and piece.isOrientable()
        assert piece.countBoundaryComponents() == 1
        assert piece.boundaryComponent(0).eulerChar() == 0
        ideal = regina.Triangulation3(piece)
        assert ideal.finiteToIdeal()
        ideal.simplify()
        if piece.isoSig() == row["finite_piece_isosig"]:
            claimed_piece_matches(row, ideal.isoSig(), census)
            found.append(ideal.isoSig())
    assert found == [row["ideal_piece_isosig"]]
    return {"name": row["name"], "independent_census_lookup": lookup,
            "closed_cover_isosig": tri.isoSig(), "irreducible": True,
            "incompressible_normal_torus": True,
            "hyperbolic_piece_census_name": row["hyperbolic_census_name"],
            "hyperbolic_piece_isosig": found[0],
            "outside_BJW18_graph_cover_class": True}


def rejection_test(call):
    try:
        call()
    except (AssertionError, ValueError, RuntimeError, IndexError):
        return True
    raise AssertionError("Deliberately invalid claim was accepted")


def main():
    versions = {p: version(p) for p in
                ["snappy", "spherogram", "snappy_15_knots", "snappy_manifolds"]}
    versions["regina"] = regina.versionString()
    assert versions == {"snappy": "3.3.2", "spherogram": "2.4.1",
                        "snappy_15_knots": "1.2.1", "snappy_manifolds": "1.4",
                        "regina": "7.4"}
    data = json.loads(INPUT.read_text())
    assert [r["name"] for r in data["cases"]] == EXPECTED
    assert set(data["primary_source_sha256"]) == {"sources/" + x for x in CENSUS_FILES + [SOURCE_PDF]}
    assert all(digest(HERE / filename) == expected for filename, expected in data["primary_source_sha256"].items())
    panel = {r["name"]: r["pd"] for r in
             json.loads((HERE / "joint_a5_panel.json").read_text())["cases"]}
    assert all(r["pd"] == panel[r["name"]] for r in data["cases"])
    census = census_table()
    rows = [check_case(row, census) for row in data["cases"]]
    bad = deepcopy(data["cases"][0])
    bad["name"] = "13a_373"
    changed_name_rejected = rejection_test(lambda: named_diagram_matches(bad))
    bad = deepcopy(data["cases"][0])
    bad["hyperbolic_census_name"] = "m004"
    wrong_census_rejected = rejection_test(
        lambda: claimed_piece_matches(bad, bad["ideal_piece_isosig"], census))
    pins = [INPUT, Path(__file__), HERE / "audit_diagram_identity.py",
            HERE / "joint_a5_panel.json", HERE / "sources" / SOURCE_PDF]
    pins.extend(HERE / "sources" / name for name in CENSUS_FILES)
    report = {"status": "PASS", "versions": versions,
              "command": "python verify_seed_cover_topology.py",
              "scope": "Twenty-eight exact diagrams; cover irreducibility and an incompressible torus cutting off a certified-census hyperbolic piece. No global novelty or full census replay.",
              "input_sha256": {str(p.relative_to(HERE)): digest(p) for p in pins},
              "cases": rows,
              "negative_controls": {"wrong_named_diagram": changed_name_rejected,
                                    "wrong_hyperbolic_census_name": wrong_census_rejected}}
    (HERE / "seed_cover_topology_verification.json").write_text(
        json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "outside_graph_cover_class": EXPECTED,
                      "negative_controls": report["negative_controls"]}, indent=2))


if __name__ == "__main__":
    main()
