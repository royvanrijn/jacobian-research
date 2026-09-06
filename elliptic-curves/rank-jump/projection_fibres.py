#!/usr/bin/env python3
"""Jacobian projection fibres as explicit two-quadrics and quartics."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import retrospective as r
import local_collision as lc
import affine_ct as old

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE/"PROJECTION_FIBRE_PROTOCOL.json"
OUTPUT = r.OUT/"rank_jump_projection_fibres_v1.json"
WORK = r.ROOT/"artifacts/local/rank-jump-projection-fibres-v1"


def bindings():
    paths = (PROTOCOL, Path(__file__), lc.INPUT, HERE/"affine_ct.py",
             old.CAS/"research_runtime/sage_subspace.py", old.CAS/"research_runtime/fisher.py")
    return {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes()) for p in paths}


def setup():
    from sage.all import QQ, EllipticCurve
    sys.path.insert(0, str(old.CAS))
    from research_runtime.subspace import GlobalSquareclasses
    old.WORK = WORK/"arithmetic"
    _, _, _, context, _, backend, _ = old.setup(0)
    backend.cap = r.read(PROTOCOL)["limits"]["local_witness_node_cap"]
    inp = r.read(lc.INPUT)
    A, B = map(QQ, inp["anchor"]["short_model_ainvariants"][3:])
    E = EllipticCurve([A, B])
    basis = [E(list(map(QQ, P))) for P in inp["anchor"]["known_points_on_short_model"]]
    points = {}
    for mask in r.read(PROTOCOL)["cover_anchor_masks"]:
        P = sum((P for i, P in enumerate(basis) if mask >> i & 1), E(0))
        assert not P.is_zero()
        points[mask] = P
    reps = [backend.coordinates(points[m][0]-backend.theta)
            for m in r.read(PROTOCOL)["anchor_basis_masks"]]
    classes = GlobalSquareclasses(context.two_torsion.key, reps,
                                  r.digest(lc.INPUT.read_bytes())+":three independent anchor masks")
    return A, B, context, backend, classes, points


def verify_geometry(A, B, p, q, u=-1):
    from sage.all import QQ, PolynomialRing
    assert q*q == p**3+A*p+B
    R = PolynomialRing(QQ, ["s", "t", "v"], order="lex")
    s, t, v = R.gens()
    kappa = 1+u*p+u*u*(A+p*p)
    Q1 = v*v+u*s*s-2*u*t-u*(u*p+2)
    Q2 = u*t*t+(1-u*p)*v*v+2*u*u*q*v-u*kappa
    ideal = R.ideal([Q1, Q2])
    Z = PolynomialRing(R, "z")
    z = Z.gen()
    D = 1+A*u*u+B*u**3
    g = u*D-(3*u+A*u**3)*z*z+3*u*z**4-u*z**6
    h = z*z-s*z+t
    n = v*(u*p-1-t)-u*u*q
    remainder = (g-(s*v*z+n)**2) % h
    assert all(ideal.reduce(c) == 0 for c in remainder)
    # The norm on the anchor has chord slope -v/u and ordinate q.
    anchor_x = v*v/(u*u)-(2-s*s+2*t)/u
    anchor_y = -n/(u*u)-v*(1+t)/(u*u)+v*anchor_x/u
    assert ideal.reduce(R(anchor_x-p)) == ideal.reduce(R(anchor_y-q)) == 0
    # The other norm is the usual degree-four map to Eu.
    F = R.fraction_field()
    e = v*(u*p-1)-u*u*q
    slope = F(v*t*t+e*(s*s-t))/(u*s*t)
    x = slope*slope+QQ(2)/u-D*(s*s-2*t)/(u*t*t)
    y = D*e*s/(u*u*t*t)-slope*(D/(u*t)+A*u+QQ(1)/u+x)
    equation = y*y-(x**3+2*A*u*x*x+(A+3*B*u+A*A*u*u)*x+B+A*B*u*u-B*B*u**3)
    assert ideal.reduce(equation.numerator()) == 0
    return {"Mumford_remainder_zero": True, "anchor_norm_verified": True,
            "specialized_norm_verified": True}


def cover_worker(index, destination):
    from sage.all import QQ, ZZ, matrix, vector, pari
    A, B, context, backend, classes, points = setup()
    from research_runtime.fisher import invariants
    protocol = r.read(PROTOCOL)
    anchor_mask = protocol["cover_anchor_masks"][index]
    mask = lc.coordinates(anchor_mask, protocol["anchor_basis_masks"])
    p, q = points[anchor_mask][0:2]
    geometry = verify_geometry(A, B, p, q)
    kappa = p*p-p+A+1
    # Coordinates (t,v,w) on the rank-three quadric.
    G = matrix(QQ, [[1,0,0], [0,-p-1,-q], [0,-q,-kappa]])
    assert G.det() == 1+A-B
    integer = matrix(ZZ, G/QQ(pari.content(G)))
    factors = pari.factor(abs(integer.det()))
    witness = pari.qfsolve([pari(integer), factors])
    assert witness.type() == "t_COL"
    witness_vector = vector(QQ, list(witness))
    assert witness_vector*G*witness_vector == 0
    parameter = matrix(QQ, pari.qfparam(pari(integer), witness, 1))
    X = backend.R.gen()
    tvw = parameter*vector([X*X, X, 1])
    t, v, w = tvw
    raw = backend.R(v*v+2*t*w+(2-p)*w*w)
    assert t*t-(p+1)*v*v-2*q*v*w-kappa*w*w == 0
    current = raw
    change = matrix.identity(QQ, 2)
    for name in ("hyperellred", "hyperellminimalmodel", "hyperellred"):
        current *= current.denominator()**2
        function = pari(f"(q)->{{my(m);my(z={name}(q,&m));[z,m]}}")
        z, m = function(current)
        current = backend.R(z[0])+backend.R(z[1])**2/4
        change *= matrix(QQ, m[1])
    if current.degree() < 4:
        k = next(k for k in range(5) if current(k))
        extra = matrix(QQ, [[k,1],[1,0]])
        a, b = extra*vector([X, 1])
        current = backend.R(sum(current[i]*a**i*b**(4-i) for i in range(5)))
        change *= extra
    I, J = invariants(current)
    scaling = QQ(backend.I/I).nth_root(4)
    current *= scaling**2
    assert invariants(current) == (backend.I, backend.J)
    a, b = change*vector([X,1])
    pull = backend.R(sum(raw[i]*a**i*b**(4-i) for i in range(5)))
    ratio = QQ(pull[4]/current[4])
    assert ratio.is_square() and pull == ratio*current
    beta = backend.beta(classes, mask)
    phi = -3*backend.alpha-backend.a2
    cubic = (4*current[4]*phi+3*current[3]**2-8*current[4]*current[2])/3
    root = backend.arithmetic.square_root(context.two_torsion, backend.coordinates(cubic/beta), discover=True)
    assert root is not None
    record = {"bindings": bindings(), "index": index, "anchor_mask": anchor_mask,
              "mask": mask, "anchor_point": [str(p),str(q)], "geometry": geometry,
              "conic_matrix": [[str(c) for c in row] for row in G.rows()],
              "conic_point": list(map(str, witness)),
              "conic_parameter_matrix": [[str(c) for c in row] for row in parameter.rows()],
              "raw_quartic": [str(raw[i]) for i in range(5)],
              "parameter_change": [[str(c) for c in row] for row in change.rows()],
              "raw_y_over_final_y": str(ratio.sqrt()),
              "quartic": [str(current[i]) for i in range(5)],
              "cubic_invariant_over_beta_square_root": list(root),
              "classes": list(classes.representatives), "context": context.record()}
    r.write_new(destination, record)
    print("DONE cover", anchor_mask, flush=True)


def pair_worker(index, destination):
    _, _, _, backend, classes, _ = setup()
    masks = r.read(PROTOCOL)["pair_class_masks"][index]
    records = [r.read(WORK/f"cover-{i}.json") for i in range(6)]
    selected = [next(c for c in records if c["mask"] == mask) for mask in masks]
    result = backend._pair(classes, masks, selected)
    old_matrix = next(x["matrix"] for x in r.read(lc.INPUT)["ct"] if x["u"] == -1)
    expected = old_matrix[{1:0,2:1,4:2}[masks[0]]][{1:0,2:1,4:2}[masks[1]]]
    assert result["value"] == expected
    r.write_new(destination, {"bindings": bindings(), "index": index, "pair": result,
                              "expected_retained_CT_value": expected})
    print("DONE pair", masks, result["value"], flush=True)


def capture():
    WORK.mkdir(parents=True, exist_ok=True)
    records = {}
    for kind, count in (("cover",6), ("pair",3)):
        rows = []
        for i in range(count):
            path = WORK/f"{kind}-{i}.json"
            if not path.exists():
                with (WORK/f"{kind}-{i}.log").open("x") as log:
                    try:
                        process = subprocess.run(["sage","-python",str(Path(__file__).resolve()),kind,
                                                  "--index",str(i),"--destination",str(path)],
                                                 cwd=r.ROOT,stdout=log,stderr=log,timeout=60)
                        failure = None if process.returncode == 0 else "worker failed"
                    except subprocess.TimeoutExpired:
                        failure = "60-second timeout"
                    if failure:
                        log.write("\nUNKNOWN: "+failure+"\n")
                        r.write_new(path, {"bindings":bindings(), "index":i,
                                           "status":"UNKNOWN", "reason":failure})
            row = r.read(path)
            assert row["bindings"] == bindings()
            rows.append(row)
            print("checkpoint",kind,i,flush=True)
        if kind == "cover" and any(c.get("status") == "UNKNOWN" for c in rows):
            raise RuntimeError("incomplete cover construction; pairing inputs unavailable")
        records[kind+"s"] = rows
    r.write_new(OUTPUT, {"schema":"rank-jump.projection-fibres.v1","bindings":bindings(),**records})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode",choices=("capture","cover","pair"))
    parser.add_argument("--index",type=int)
    parser.add_argument("--destination",type=Path)
    args = parser.parse_args()
    if args.mode == "capture":
        capture()
    elif args.mode == "cover":
        cover_worker(args.index,args.destination)
    else:
        pair_worker(args.index,args.destination)
