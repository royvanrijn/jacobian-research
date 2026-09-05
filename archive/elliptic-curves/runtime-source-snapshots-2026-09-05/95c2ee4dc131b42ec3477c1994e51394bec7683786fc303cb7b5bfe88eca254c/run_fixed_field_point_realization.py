#!/usr/bin/env sage
"""Checkpointed searches of fixed-field covers, with exact maps and replay.

Run with sage -python. A cover is converted to a ternary conic and then a
binary quartic. PARI minimizes the conic and quartic and reduces the quartic.
Every coordinate change and every returned point is checked over QQ. Search
misses and timeouts never decide Sha or give a Mordell--Weil upper bound.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import subprocess
import sys
import time

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, gcd, lcm, matrix, pari, vector
from sage.version import version as sage_version

from fixed_cubic_field_curve_family import field_multiply, field_product, f2_rank

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "artifacts/generated-results/elliptic-curves/fixed_cubic_field_fermigier_rank20_local_kummer_u2_v1.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/fixed_field_u_minus1_point_realization_v1.json"
SCHEMA = "elliptic-curves.fixed-field-point-realization.v1"
R = PolynomialRing(QQ, "x")
x = R.gen()


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def rows(M):
    return [[str(v) for v in row] for row in M.rows()]


def coeffs(f, length=5):
    f = R(f)
    return [str(f[i]) for i in range(length)]


def poly(values):
    return R([QQ(v) for v in values])


def primitive(M):
    M = matrix(QQ, M)
    scale = lcm(v.denominator() for v in M.list())
    scale /= gcd(ZZ(v * scale) for v in M.list())
    return scale * M, scale


def class_masks(basis, max_weight):
    answer = []
    for weight in range(1, max_weight + 1):
        for indices in combinations(range(len(basis)), weight):
            mask = 0
            for i in indices:
                mask ^= basis[i]["mask"]
            answer.append((mask, [i + 1 for i in indices]))
    if len({mask for mask, _ in answer}) != len(answer):
        raise ArithmeticError("dependent input basis")
    return answer


def context(source, u):
    data = json.loads(Path(source).read_text())
    if data["status"] != "PASS_EXACT_FULL_SPAN_LOCAL_INTERSECTIONS_NO_CLASS_GROUP":
        raise ArithmeticError("local certificate is incomplete")
    run = next(r for r in data["runs"] if QQ(r["parameter_u"]) == u)
    assert run["all_local_kummer_images_complete"]
    B, A = map(QQ, data["anchor"]["base_polynomial_ascending"][:2])
    E = EllipticCurve(QQ, list(map(QQ, run["raw_curve_ainvariants"])))
    return data, run, A, B, E


def class_input(data, mask, A, B):
    points = data["anchor"]["known_points_on_short_model"]
    indices = [i for i in range(len(points)) if mask >> i & 1]
    beta = field_product([(QQ(points[i][0]), QQ(-1), QQ(0)) for i in indices], A, B)
    norm_root = QQ(1)
    for i in indices:
        norm_root *= QQ(points[i][1])
    assert norm(beta, A, B) == norm_root**2
    return beta, norm_root


def norm(v, A, B):
    return matrix(QQ, [field_multiply(v, e, A, B) for e in [(1,0,0),(0,1,0),(0,0,1)]]).det()


def divide(left, right, A, B):
    M = matrix(QQ, [field_multiply(right, e, A, B) for e in [(1,0,0),(0,1,0),(0,0,1)]]).transpose()
    return tuple(M.solve_right(vector(QQ,left)))


def universal_point_certificate(data, A, B, E):
    """At u=-1, a rational point and its independence from all of W."""
    Q = E(A+1,A-B+1)
    assert 2*Q == E(QQ(1)/4, B+A/2+QQ(1)/8)
    eta = (A+1,QQ(-1),QQ(1))
    assert norm(eta,A,B) == (A-B+1)**2
    f = pari("y^3+("+str(A)+")*y+("+str(B)+")")
    nf = pari.nfinit(f)
    theta = pari.Mod(pari("y"),f)
    primes = pari.idealprimedec(nf,19)
    beta = theta*theta-theta+A+1
    valuations = [int(pari.nfeltval(nf,beta,P)) for P in primes]
    anchors = [[int(pari.nfeltval(nf,pari(QQ(pt[0]))-theta,P)) for P in primes]
               for pt in data["anchor"]["known_points_on_short_model"]]
    assert all(all(v%2==0 for v in row) for row in anchors)
    assert any(v%2 for v in valuations)
    assert (x**3+A*x+B).is_irreducible()
    return {"point":list(map(str,Q[:2])), "double":list(map(str,(2*Q)[:2])),
            "kummer_beta":list(map(str,eta)), "norm_square_root":str(A-B+1),
            "separating_rational_prime":19,
            "prime_ideal_residue_degrees":[int(P[3]) for P in primes],
            "point_kummer_valuations":valuations,"anchor_kummer_valuations":anchors,
            "independent_from_entire_anchor_Kummer_span":True,
            "certified_rank_lower_bound":1,
            "reason":"Nonzero valuation parity separates eta from W; irreducibility excludes rational 2-torsion."}


def translate_back(point, beta_original, norm_original, beta_total, A, B, E):
    """Realize the original mask by subtracting Q, with an explicit square root."""
    P = E(*map(QQ,point["raw_point"]))
    Q = E(A+1,A-B+1)
    T = P-Q
    assert P[0] != Q[0] and not T.is_zero()
    slope = (P[1]+Q[1])/(P[0]-Q[0])
    intercept = P[1]-slope*P[0]
    line_at_alpha = (intercept,slope,-slope)
    g = tuple(map(QQ,point["kummer_ratio_square_root"]))
    h = divide(line_at_alpha,field_multiply(beta_total,g,A,B),A,B)
    if norm_original*norm(h,A,B) != T[1]:
        h = tuple(-v for v in h)
    witness = verify_point(beta_original,norm_original,h,1,A,B,-1,E)
    assert witness["raw_point"] == list(map(str,T[:2]))
    return witness


def quadric_matrices(beta, A, B, u):
    S = PolynomialRing(QQ, "a,b,c")
    vals = field_multiply(beta, field_multiply(S.gens(), S.gens(), A, B), A, B)
    matrices = [matrix(QQ, 3, 3, lambda i,j: f.derivative(S.gen(i)).derivative(S.gen(j))/2) for f in vals]
    return matrices, primitive(matrices[2] - u * matrices[1])[0]


def model_record(C):
    return [coeffs(R(C[0])), coeffs(R(C[1]), 3)]


def change_record(m):
    return {"e": str(m[0]), "matrix": rows(matrix(QQ, m[1])), "H": coeffs(R(m[2]), 3)}


def homogeneous(f, s, t, degree):
    return sum(QQ(v) * s**i * t**(degree-i) for i, v in enumerate(f))


def verify_change(old, new, change):
    e = QQ(change["e"])
    a,b,c,d = matrix(QQ, change["matrix"]).list()
    H = poly(change["H"])
    assert e and a*d-b*c
    pp = homogeneous(old[0], a*x+b, c*x+d, 4)
    qp = homogeneous(old[1], a*x+b, c*x+d, 2)
    assert pp - H**2 - qp*H == e**2 * poly(new[0])
    assert 2*H + qp == e * poly(new[1])


def inverse_change(point, change):
    s,t,v = map(QQ, point)
    a,b,c,d = matrix(QQ, change["matrix"]).list()
    return (a*s+b*t, c*s+d*t, QQ(change["e"])*v + homogeneous(change["H"], s,t,2))


def verify_point(beta, norm_root, gamma, d, A, B, u, E):
    """Compute the actual class: (x-alpha)/beta=(gamma/d)^2 in K."""
    gamma = tuple(map(QQ, gamma))
    d = QQ(d)
    if not d:
        # A nonzero class cannot lie over O; never silently discard this chart.
        raise ArithmeticError("nontrivial cover point over infinity requires separate audit")
    values = field_multiply(beta, field_multiply(gamma, gamma, A, B), A, B)
    assert values[1] == -d**2 and values[2] == -u*d**2
    px = values[0] / d**2
    py = norm_root * norm(gamma, A, B) / d**3
    P = E(px, py)
    square_root = [v/d for v in gamma]
    actual = field_multiply(beta, field_multiply(square_root, square_root, A, B), A, B)
    assert actual == (px, -QQ(1), -QQ(u))
    return {"raw_point": [str(P[0]), str(P[1])],
            "gamma": [str(v) for v in gamma], "d": str(d),
            "actual_kummer_power_coefficients": [str(v) for v in actual],
            "kummer_ratio_square_root": [str(v) for v in square_root],
            "exact_curve_and_kummer_identity_verified": True}


def worker(source, u, mask, height, checkpoint, translate=False):
    started = time.monotonic()
    data, run, A, B, E = context(source, u)
    beta_original, norm_original = class_input(data, mask, A, B)
    beta, norm_root = beta_original, norm_original
    if translate:
        assert u == -1
        beta = field_multiply(beta,(A+1,-1,1),A,B)
        norm_root *= A-B+1
    record = {"mask": mask, "parameter_u": str(u), "height_bound": height,
              "beta": [str(v) for v in beta], "norm_square_root": str(norm_root),
              "source_sha256": digest(source), "translated_by_universal_point": translate,
              "points": [], "stage": "conic_minimization"}
    def checkpoint_now(stage):
        record["stage"] = stage
        record["elapsed_seconds"] = round(time.monotonic()-started, 6)
        save(checkpoint, record)
    checkpoint_now("conic_minimization")
    pari.allocatemem(128_000_000, silent=True)
    pari.set_real_precision(100)
    pari.setrand(1)
    Qs, G = quadric_matrices(beta, A, B, u)
    record["coefficient_quadrics"] = [rows(q) for q in Qs]
    record["conic_matrix"] = rows(G)
    H, U, scale = pari(G).qfminimize()
    H, U, scale = matrix(QQ,H), matrix(QQ,U), QQ(scale)
    assert H == scale * U.transpose()*G*U and U.det()
    sol = pari(H).qfsolve()
    if sol.type() != "t_COL":
        raise ArithmeticError("local-surviving class has no conic point")
    sol = vector(QQ,sol)
    assert sol*H*sol == 0
    M, mscale = primitive(U * matrix(QQ, pari(H).qfparam(pari(sol), 3)))
    vv = M * vector([x*x,x,1])
    assert vv*G*vv == 0 and M.det()
    quartic = R(-(vv*Qs[1]*vv))
    L = lcm(v.denominator() for v in quartic.list())
    C0 = [coeffs(L**2*quartic), ["0"]*3]
    record.update({"minimized_conic_matrix": rows(H), "conic_change": rows(U),
                   "conic_scale": str(scale), "conic_point": [str(v) for v in sol],
                   "conic_parameter_matrix": rows(M), "quartic_ordinate_scale": str(L),
                   "raw_quartic_model": C0})
    checkpoint_now("quartic_minimization")
    pari("C0=["+str(poly(C0[0])).replace("**","^")+",0]")
    out = pari("C1=hyperellminimalmodel(C0,&m1); [C1,m1]")
    C1, m1 = model_record(out[0]), change_record(out[1])
    verify_change(C0,C1,m1)
    record.update({"minimal_quartic_model": C1, "minimization_change": m1})
    checkpoint_now("quartic_reduction")
    out = pari("C2=hyperellred(C1,&m2); [C2,m2]")
    C2, m2 = model_record(out[0]), change_record(out[1])
    verify_change(C1,C2,m2)
    record.update({"reduced_quartic_model": C2, "reduction_change": m2,
                   "reduced_max_coefficient_bits": max(abs(QQ(v).numerator()).nbits() for f in C2 for v in f),
                   "model_changes_exactly_verified": True})
    checkpoint_now("point_search")
    found = [[QQ(p[0]), QQ(1), QQ(p[1])] for p in pari("hyperellratpoints(C2,"+str(height)+")")]
    # The search is affine. Also test every rational point over the parameter infinity.
    infinity_poly = x*x + QQ(C2[1][2])*x - QQ(C2[0][4])
    found.extend([[QQ(1),QQ(0),v] for v,multiplicity in infinity_poly.roots(QQ)])
    record["searched_affine_point_count"] = sum(p[1] != 0 for p in found)
    record["rational_parameter_infinity_point_count"] = sum(p[1] == 0 for p in found)
    for pt in found:
        s,t,v = pt
        assert v*v + homogeneous(C2[1],s,t,2)*v == homogeneous(C2[0],s,t,4)
        old = inverse_change(inverse_change(pt,m2),m1)
        s,t,v = old
        gamma = M*vector([s*s,s*t,t*t])
        point = verify_point(beta,norm_root,gamma,v/L,A,B,u,E)
        point.update({"reduced_quartic_point": list(map(str,pt)), "mask": mask})
        if translate:
            point["original_W_class_realization"] = translate_back(point,beta_original,norm_original,beta,A,B,E)
        record["points"].append(point)
        # A nonzero class is immediately a rank-one certificate (E[2](Q)=0).
        record["realized_class_rank"] = 1
        checkpoint_now("point_verification")
    record["status"] = "POINTS_CERTIFIED" if record["points"] else "NO_POINT_WITHIN_DECLARED_BOUND"
    checkpoint_now("complete")


def replay_row(row, data, A, B, u, E):
    beta, root = class_input(data,row["mask"],A,B)
    beta_original, root_original = beta, root
    if row.get("translated_by_universal_point",False):
        assert u == -1
        beta = field_multiply(beta,(A+1,-1,1),A,B)
        root *= A-B+1
    assert list(map(QQ,row["beta"])) == list(beta) and QQ(row["norm_square_root"]) == root
    Qs,G = quadric_matrices(beta,A,B,u)
    if "conic_parameter_matrix" not in row:
        return
    assert [rows(q) for q in Qs] == row["coefficient_quadrics"]
    assert G == matrix(QQ,row["conic_matrix"])
    H,U = matrix(QQ,row["minimized_conic_matrix"]),matrix(QQ,row["conic_change"])
    assert H == QQ(row["conic_scale"])*U.transpose()*G*U and U.det()
    sol = vector(QQ,row["conic_point"])
    assert sol*H*sol == 0 and any(sol)
    M = matrix(QQ,row["conic_parameter_matrix"])
    vv = M*vector([x*x,x,1])
    assert vv*G*vv == 0 and M.det()
    L = QQ(row["quartic_ordinate_scale"])
    assert -L**2*(vv*Qs[1]*vv) == poly(row["raw_quartic_model"][0])
    if "minimal_quartic_model" in row:
        verify_change(row["raw_quartic_model"],row["minimal_quartic_model"],row["minimization_change"])
    if "reduced_quartic_model" in row:
        verify_change(row["minimal_quartic_model"],row["reduced_quartic_model"],row["reduction_change"])
    for point in row["points"]:
        pt = tuple(map(QQ,point["reduced_quartic_point"]))
        s,t,v = pt
        C2 = row["reduced_quartic_model"]
        assert v*v + homogeneous(C2[1],s,t,2)*v == homogeneous(C2[0],s,t,4)
        s,t,v = inverse_change(inverse_change(pt,row["reduction_change"]),row["minimization_change"])
        gamma = M*vector([s*s,s*t,t*t])
        assert list(gamma) == list(map(QQ,point["gamma"])) and v/L == QQ(point["d"])
        checked = verify_point(beta,root,gamma,v/L,A,B,u,E)
        for key,value in checked.items():
            assert point[key] == value
        if row.get("translated_by_universal_point",False):
            assert point["original_W_class_realization"] == translate_back(point,beta_original,root_original,beta,A,B,E)


def summarize(document):
    records = document["covers"]
    realized = [r["mask"] for r in records if r.get("points")]
    rank = f2_rank([[mask >> j & 1 for j in range(20)] for mask in realized])
    document["summary"] = {
        "completed_searches": sum(r.get("stage")=="complete" for r in records),
        "reduced_covers": sum("reduced_quartic_model" in r for r in records),
        "timeouts": sum(r.get("status")=="TIMEOUT" for r in records),
        "errors": sum(r.get("status")=="ERROR" for r in records),
        "verified_points_including_signs": sum(len(r.get("points",[])) for r in records),
        "realized_class_masks": realized, "certified_realized_subspace_dimension": rank,
        "unresolved_dimension_in_W": document["W_dimension"]-rank,
        "new_rank_lower_bound": rank + (1 if document.get("universal_point") else 0) or None,
        "Sha_classification": "UNKNOWN", "exact_rank": "UNKNOWN"}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source",type=Path,default=SOURCE)
    ap.add_argument("--output",type=Path,default=DEFAULT_OUTPUT)
    ap.add_argument("--u",type=int,default=-1)
    ap.add_argument("--max-weight",type=int,default=2)
    ap.add_argument("--height",type=int,default=100000)
    ap.add_argument("--timeout",type=float,default=60)
    ap.add_argument("--jobs",type=int,default=2)
    ap.add_argument("--checkpoint-dir",type=Path,default=ROOT/"artifacts/local/fixed-field-u-minus1")
    ap.add_argument("--worker",type=int)
    ap.add_argument("--translate",action="store_true")
    ap.add_argument("--check",action="store_true")
    args = ap.parse_args()
    if args.worker is not None:
        worker(args.source,args.u,args.worker,args.height,args.output,args.translate)
        return
    data,run,A,B,E = context(args.source,args.u)
    if args.check:
        doc = json.loads(args.output.read_text())
        assert doc["schema"] == SCHEMA
        assert doc["source_sha256"] == digest(args.source)
        assert doc["parameter_u"] == args.u
        if doc.get("universal_point"):
            assert doc["universal_point"] == universal_point_certificate(data,A,B,E)
        expected = class_masks(run["W_u_basis"], doc["policy"]["max_basis_weight"])
        assert [(r["mask"],r["basis_indices"]) for r in doc["covers"]] == expected
        minimal_pari,_ = pari.ellinit(list(E.a_invariants())).ellminimalmodel()
        minimal = EllipticCurve(QQ,[QQ(minimal_pari[i]) for i in range(5)])
        iso = E.isomorphism_to(minimal)
        assert doc["global_minimal_curve_ainvariants"] == list(map(str,minimal.a_invariants()))
        assert doc["raw_to_minimal_urst"] == list(map(str,iso.tuple()))
        for row in doc["covers"]:
            assert row["height_bound"] == doc["policy"]["quartic_height_bound"]
            assert row.get("translated_by_universal_point",False) == doc["policy"].get("translate_by_universal_point",False)
            replay_row(row,data,A,B,args.u,E)
            if row["stage"] == "complete":
                assert row["status"] == ("POINTS_CERTIFIED" if row["points"] else "NO_POINT_WITHIN_DECLARED_BOUND")
            for point in row["points"]:
                assert point["minimal_point"] == list(map(str,iso(E(*map(QQ,point["raw_point"])))[:2]))
        previous = doc["summary"]
        summarize(doc)
        assert doc["summary"] == previous
        print("PASS_EXACT_COVER_MAP_AND_POINT_REPLAY",doc["summary"],flush=True)
        return
    assert 1 <= args.max_weight <= len(run["W_u_basis"]) and args.height > 0 and args.timeout > 0 and args.jobs > 0
    policy = class_masks(run["W_u_basis"],args.max_weight)
    # PARI's model minimization only factors the required invariant gcd;
    # Sage's default path factors the much larger full discriminant first.
    minimal_pari, _ = pari.ellinit(list(E.a_invariants())).ellminimalmodel()
    minimal = EllipticCurve(QQ, [QQ(minimal_pari[i]) for i in range(5)])
    iso = E.isomorphism_to(minimal)
    doc = {"schema": SCHEMA, "status": "RUNNING", "parameter_u": args.u,
           "source_sha256": digest(args.source), "runner_sha256": digest(__file__),
           "software": {"sage":sage_version,"pari":str(pari.version())},
           "W_dimension": run["W_u_dimension"], "raw_curve_ainvariants":list(map(str,E.a_invariants())),
           "global_minimal_curve_ainvariants":list(map(str,minimal.a_invariants())),
           "raw_to_minimal_urst":list(map(str,iso.tuple())),
           "policy": {"max_basis_weight":args.max_weight,"cover_count":len(policy),
                      "quartic_height_bound":args.height,"per_cover_timeout_seconds":args.timeout,
                      "workers":args.jobs,"pari_random_seed":1,"real_precision_decimal_digits":100,
                      "parameter_infinity_checked":True,"translate_by_universal_point":args.translate},
           "independence_method":"Exact (x-alpha)/beta square identity, then GF(2) rank in the certified independent anchor Kummer basis; irreducibility gives E(Q)[2]=0.",
           "claim_boundary":["Equivalent quartic models are minimized and reduced; this is not a degree-four-model minimality theorem.",
                             "A completed bounded search miss does not prove global insolubility or a nontrivial Sha class.",
                             "Only the selected masks and the declared reduced-coordinate height are searched.",
                             "Realized mask rank is a lower bound; dependent masks can conceal further independent points."],
           "covers":[]}
    if args.u == -1:
        doc["universal_point"] = universal_point_certificate(data,A,B,E)
    args.checkpoint_dir.mkdir(parents=True,exist_ok=True)
    def launch(entry):
        mask,indices = entry
        suffix = "_translated" if args.translate else ""
        target = args.checkpoint_dir/f"u{args.u}_mask{mask}_h{args.height}{suffix}.json"
        log = target.with_suffix(".log")
        if target.exists():
            old = json.loads(target.read_text())
            if (old.get("stage")=="complete" and old.get("source_sha256")==digest(args.source)
                and old.get("translated_by_universal_point",False)==args.translate
                and old.get("height_bound")==args.height):
                old["basis_indices"] = indices
                return old
        cmd = [sys.executable,str(Path(__file__).resolve()),"--source",str(args.source),"--u",str(args.u),
               "--worker",str(mask),"--height",str(args.height),"--output",str(target)]
        if args.translate:
            cmd.append("--translate")
        status = None
        with log.open("w") as stream:
            try:
                result = subprocess.run(cmd,stdout=stream,stderr=stream,timeout=args.timeout)
                if result.returncode:
                    status = "ERROR"
            except subprocess.TimeoutExpired:
                status = "TIMEOUT"
        row = json.loads(target.read_text()) if target.exists() else {"mask":mask,"stage":"startup","points":[]}
        if status:
            row["status"] = status
            row["diagnostic"] = log.read_text()[-2000:]
        row["basis_indices"] = indices
        return row
    # Submit the basis first and finish it before the pair/higher-weight wave.
    indexed = {}
    for weight in range(1,args.max_weight+1):
        wave = [entry for entry in policy if len(entry[1])==weight]
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            for future in as_completed([pool.submit(launch,entry) for entry in wave]):
                row = future.result()
                # Sage arithmetic stays on the main thread; workers here only
                # supervise separate CAS processes and read checkpoint files.
                if "beta" in row:
                    replay_row(row,data,A,B,args.u,E)
                for pt in row.get("points",[]):
                    P = E(*map(QQ,pt["raw_point"]))
                    pt["minimal_point"] = list(map(str,iso(P)[:2]))
                indexed[row["mask"]] = row
                doc["covers"] = [indexed[mask] for mask,_ in policy if mask in indexed]
                summarize(doc)
                save(args.output,doc)
                print(f"COVER {len(indexed)}/{len(policy)} mask={row['mask']} weight={weight} status={row.get('status')} points={len(row.get('points',[]))} rank={doc['summary']['certified_realized_subspace_dimension']}",flush=True)
    doc["status"] = "BOUNDED_CAMPAIGN_COMPLETE" if not (doc["summary"]["errors"] or doc["summary"]["timeouts"]) else "BOUNDED_CAMPAIGN_WITH_UNRESOLVED_JOBS"
    save(args.output,doc)
    print(json.dumps(doc["summary"],sort_keys=True),flush=True)


if __name__ == "__main__":
    main()
