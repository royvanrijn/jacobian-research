#!/usr/bin/env sage-python
"""Witness-hidden, post-search calibration of degree-four covers on E302.

The preparation phase may use the displayed rank-31 group only to make the
two test panels.  The solver phase runs in a fresh process and reads one
sealed pair of quadrics only.  In particular it never receives alpha, the
strict word, the public points, an expected point, or the control/missing
label.  Returned points are checked only after that process has exited.

This is a calibration gate, not a new point discovery or a Selmer/class-group
calculation.  A missing-panel solution is explicitly excluded from the next
principal-relation constructor's inputs.
"""

from __future__ import annotations

import argparse
import ast
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import subprocess
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

from sage.all import GF, QQ, ZZ, EllipticCurve, PolynomialRing, lcm, matrix, vector
from sage.version import version as SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
ART = ROOT / "artifacts" / "generated-results" / "elliptic-curves"
FILTER = ART / "curve302_recovered_quotient_local_filtration_v1.json"
SPAN = ART / "curve302_recovered_public_span_v1"
PARENT = ART / "curve302_recovered_mw17_parent_v1.json"
PUBLIC_SOURCE = CAS / "icarm_curve302.py"
BUILDER = CAS / "build_bnf_free_two_covers.py"
PROTOCOL = CAS / "CURVE302_STRICT_COVER_HIDDEN_BENCHMARK_PROTOCOL.json"
OUTPUT = ART / "curve302_strict_cover_hidden_benchmark_v1.json"
WORK = ROOT / "artifacts" / "local" / "elliptic-curves" / "curve302-strict-cover-hidden-benchmark-v1"

sys.path.insert(0, str(CAS))
from build_bnf_free_two_covers import multiply_mod_cubic  # noqa: E402
import icarm_curve302 as curve  # noqa: E402


def read(path: Path):
    return json.loads(path.read_text())


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def put_new(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def rational_vector(values):
    return [QQ(value) for value in values]


def xor(*rows):
    answer = [0] * 31
    for row in rows:
        assert len(row) == 31
        answer = [left ^ int(right) for left, right in zip(answer, row)]
    return answer


def binary_combinations(basis):
    assert len(basis) == 3
    answer = []
    for mask in range(8):
        answer.append(xor(*(basis[i] for i in range(3) if (mask >> i) & 1)) if mask else [0] * 31)
    assert len({tuple(row) for row in answer}) == 8
    return answer


def point_beta(point, f, ring):
    """Integral [4x-theta] with a declared square-root of its norm."""
    x, y = map(QQ, point.xy())
    X, Y = 4 * x, 8 * y + 4 * x + 4
    assert Y * Y == f(X)
    den = ZZ(X.denominator()).sqrt()
    assert den in ZZ and den * den == X.denominator()
    beta = ring([ZZ(X * den**2), -den**2])
    return [QQ(beta[i]) for i in range(3)], QQ(Y * den**3)


def primitive_integral(poly):
    values = poly.coefficients()
    scale = lcm([value.denominator() for value in values])
    integral = poly * scale
    content = 0
    for value in integral.coefficients():
        content = gcd(content, abs(int(value)))
    assert content
    return integral / content


def norm_form(beta, coefficients, ring):
    u, v, w, _ = ring.gens()
    basis = ([1, 0, 0], [0, 1, 0], [0, 0, 1])
    return matrix(ring, [multiply_mod_cubic(beta, row, coefficients) for row in basis]).det()


def build_panels():
    filtration = read(FILTER)
    protocol = read(PROTOCOL)
    assert filtration["status"] == "PASS"
    local = filtration["local_filtration_mod_2"]
    assert filtration["displayed_lattices"]["M24_rank"] == 24
    assert local["strict_kernel_dimensions"] == {"M17": 0, "M24": 3, "D": 10}
    h = [[int(value) for value in row] for row in local["M24_strict_kernel_public_words"]]
    assert len(h) == 3 and matrix(GF(2), h).rank() == 3
    missing = [int(value) for value in filtration["explicit_strict_witness"]["public_word_mod_2"]]
    m24_words = matrix(ZZ, [row["word"] for row in read(SPAN / "result.json")["relations"]]).transpose()
    m24 = matrix(GF(2), m24_words).column_space()
    assert all(vector(GF(2), row) in m24 for row in h)
    assert vector(GF(2), missing) not in m24

    E = EllipticCurve(QQ, list(map(QQ, curve.GENERAL_WEIERSTRASS_COEFFICIENTS)))
    public = [E(point) for point in curve.POINTS]
    ring = PolynomialRing(QQ, "theta")
    a1, a2, a3, a4, a6 = map(QQ, curve.GENERAL_WEIERSTRASS_COEFFICIENTS)
    assert (a1, a2, a3) == (1, 1, 1)
    f = ring([64 * a6 + 16, 16 * a4 + 8, 5, 1])
    coefficients = [QQ(value) for value in f.list()]
    betas, roots = zip(*(point_beta(point, f, ring) for point in public))
    variables = PolynomialRing(QQ, names=("u", "v", "w", "z"))

    def make_case(token, cohort, member, word):
        alpha = [QQ(1), QQ(0), QQ(0)]
        root = QQ(1)
        for index, bit in enumerate(word):
            if bit:
                alpha = multiply_mod_cubic(alpha, betas[index], coefficients)
                root *= roots[index]
        alpha = [QQ(value) for value in alpha]
        u, v, w, z = variables.gens()
        products = multiply_mod_cubic(alpha, multiply_mod_cubic([u, v, w], [u, v, w], coefficients), coefficients)
        first, second = primitive_integral(products[1] + z**2), primitive_integral(products[2])
        # determinant over the power basis agrees with the declared product norm.
        Kalpha = matrix(QQ, [multiply_mod_cubic(alpha, row, coefficients) for row in ([1, 0, 0], [0, 1, 0], [0, 0, 1])]).det()
        assert Kalpha == root * root
        return {
            "token": token,
            "cohort": cohort,
            "member": member,
            "public_word_mod_2": word,
            "alpha_coefficients": [str(value) for value in alpha],
            "norm_square_root": str(root),
            "integral_quadrics": [str(first), str(second)],
        }

    h_members = binary_combinations(h)
    cases = []
    for member, row in enumerate(h_members):
        cases.append(make_case(f"case-{member:02d}", "recovered_control", member, row))
    for member, row in enumerate(h_members):
        cases.append(make_case(f"case-{member + 8:02d}", "missing_strict_coset", member, xor(missing, row)))
    assert len(cases) == 16
    return {
        "protocol": protocol,
        "filtration": filtration,
        "m24_words": m24_words,
        "curve": E,
        "f": f,
        "coefficients": coefficients,
        "ring": variables,
        "cases": cases,
    }


def magma_program(quadrics, height, preparation="minimise_reduce"):
    first, second = quadrics
    if preparation == "reduce_only":
        return f'''SetColumns(0); SetSeed(1); Q:=Rationals();
S<u,v,w,z>:=PolynomialRing(Q,4);
g:=GenusOneModel([{first},{second}]);
gr,rr:=Reduce(g);
print "QI_MATRICES",[Eltseq(m):m in Matrices(gr)];
print "EQ_TRANS",Eltseq(Tuple(rr)[1]);
print "VAR_TRANS",Eltseq(Tuple(rr)[2]);
SetVerbose("QISearch",1); print "SEARCH_START";
pts:=PointsQI(Curve(gr),{height} : OnlyOne:=true);
print "POINTS",[Eltseq(pt):pt in pts]; print "DONE";
'''
    assert preparation == "minimise_reduce"
    return f'''SetColumns(0); SetSeed(1); Q:=Rationals();
S<u,v,w,z>:=PolynomialRing(Q,4);
g:=GenusOneModel([{first},{second}]);
gm,tr:=Minimise(g); gr,rr:=Reduce(gm);
print "QI_MATRICES",[Eltseq(m):m in Matrices(gr)];
print "EQ_TRANS",Eltseq(Tuple(rr*tr)[1]);
print "VAR_TRANS",Eltseq(Tuple(rr*tr)[2]);
SetVerbose("QISearch",1); print "SEARCH_START";
pts:=PointsQI(Curve(gr),{height} : OnlyOne:=true);
print "POINTS",[Eltseq(pt):pt in pts]; print "DONE";
'''


def parse_xml(raw):
    root = ET.fromstring(raw)
    text = "\n".join("".join(row.itertext()) for row in root.findall(".//results/line"))
    bad = ("Runtime error", "User error", "System Error", "Syntax error", "memory limit", "Internal error")
    complete = not root.findall(".//warning") and "DONE" in text and not any(item in text for item in bad)
    return complete, text


def literal(text, label):
    start = text.index(label) + len(label)
    start = text.index("[", start)
    depth = 0
    for end in range(start, len(text)):
        depth += (text[end] == "[") - (text[end] == "]")
        if depth == 0:
            source = text[start:end + 1]
            assert __import__("re").fullmatch(r"[\[\],\s0-9+\-/]+", source)
            tree = ast.parse(source, mode="eval").body
            def decode(node):
                if isinstance(node, ast.List):
                    return [decode(child) for child in node.elts]
                if isinstance(node, ast.Constant) and type(node.value) is int:
                    return QQ(node.value)
                if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
                    return -decode(node.operand)
                if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
                    return decode(node.left) / decode(node.right)
                raise ValueError("unexpected Magma numeric syntax")
            return decode(tree)
    raise ValueError("unterminated Magma list")


def solver(sealed: Path) -> None:
    """The isolated worker: it deliberately reads only sealed cover equations."""
    job = read(sealed)
    assert set(job) == {"schema", "token", "integral_quadrics", "solver"}
    assert job["schema"] == "elliptic-curves.sealed-two-quadric-job.v1"
    assert len(job["integral_quadrics"]) == 2
    request = urllib.request.Request(
        "https://magma.maths.usyd.edu.au/xml/calculator.xml",
        data=urllib.parse.urlencode({"input": magma_program(
            job["integral_quadrics"], job["solver"]["pointsqi_height"],
            job["solver"].get("model_preparation", "minimise_reduce"),
        )}).encode(),
        headers={"Referer": "https://magma.maths.usyd.edu.au/calc/", "Content-Type": "application/x-www-form-urlencoded"},
    )
    started = time.monotonic()
    result = {"schema": "elliptic-curves.sealed-two-quadric-solver-result.v1", "token": job["token"],
              "sealed_input_sha256": digest(sealed), "engine": job["solver"]["engine"],
              "status": "UNKNOWN", "wall_seconds": None}
    try:
        raw = urllib.request.urlopen(request, timeout=job["solver"]["remote_request_timeout_seconds"]).read().decode()
        complete, text = parse_xml(raw)
        result["magma_xml"] = raw
        result["magma_text_sha256"] = sha256(text.encode()).hexdigest()
        result["status"] = "COMPLETED" if complete else "SOLVER_ERROR_OR_INCOMPLETE"
        if complete:
            result["quadric_matrices"] = [[str(value) for value in row] for row in literal(text, "QI_MATRICES")]
            result["equation_transform"] = [str(value) for value in literal(text, "EQ_TRANS")]
            result["variable_transform"] = [str(value) for value in literal(text, "VAR_TRANS")]
            result["returned_points"] = [[str(value) for value in row] for row in literal(text, "POINTS")]
    except Exception as error:  # error is retained; it is never an arithmetic conclusion.
        result["status"] = "TRANSPORT_OR_PARSE_FAILURE"
        result["error"] = repr(error)
    result["wall_seconds"] = time.monotonic() - started
    put_new(sealed.with_name(sealed.stem + "-result.json"), result)
    print(job["token"], result["status"], flush=True)


def make_sealed(case, protocol):
    return {
        "schema": "elliptic-curves.sealed-two-quadric-job.v1",
        "token": case["token"],
        "integral_quadrics": case["integral_quadrics"],
        "solver": {
            "engine": protocol["solver"]["engine"],
            "pointsqi_height": protocol["solver"]["pointsqi_height"],
            "remote_request_timeout_seconds": protocol["solver"]["remote_request_timeout_seconds"],
        },
    }


def map_solution(case, solved, panel):
    """Post-worker exact verification and M24 escape test."""
    if solved["status"] != "COMPLETED":
        return {"solver_status": solved["status"], "exact_rational_cover_points": []}
    ring = panel["ring"]
    u, v, w, z = ring.gens()
    first, second = [ring(text) for text in case["integral_quadrics"]]
    transform = matrix(QQ, 4, [QQ(value) for value in solved["variable_transform"]])
    alpha = [QQ(value) for value in case["alpha_coefficients"]]
    root = QQ(case["norm_square_root"])
    # This is N(u+v*theta+w*theta^2), not N(alpha).  The latter is only
    # the fixed square that labels the cover.
    norm = norm_form([u, v, w], panel["coefficients"], ring)
    result = []
    for raw in solved["returned_points"]:
        point = vector(QQ, [QQ(value) for value in raw]) * transform
        assert any(point) and first(*point) == 0 and second(*point) == 0
        if point[3] == 0:
            result.append({"primitive_cover_point": [str(value) for value in point], "map": "AT_INFINITY"})
            continue
        X = multiply_mod_cubic(alpha, multiply_mod_cubic([u, v, w], [u, v, w], panel["coefficients"]), panel["coefficients"])[0](*point) / point[3]**2
        Y = root * norm(*point) / point[3]**3
        assert Y * Y == panel["f"](X)
        original = panel["curve"](X / 4, (Y - X - 4) / 8)
        result.append({"primitive_cover_point": [str(value) for value in point],
                       "cubic_point": [str(X), str(Y)],
                       "point_on_E302": [str(original[0]), str(original[1])]})
    return {"solver_status": solved["status"], "exact_rational_cover_points": result}


def capture():
    panel = build_panels()
    protocol = panel["protocol"]
    WORK.mkdir(parents=True, exist_ok=True)
    records = []
    for case in panel["cases"]:
        sealed = WORK / f"{case['token']}.json"
        result_path = sealed.with_name(sealed.stem + "-result.json")
        if not sealed.exists():
            put_new(sealed, make_sealed(case, protocol))
        else:
            assert read(sealed) == make_sealed(case, protocol)
        if not result_path.exists():
            try:
                subprocess.run(["sage", "-python", str(Path(__file__)), "--solver", "--sealed", str(sealed)],
                               check=False, timeout=protocol["solver"]["subprocess_timeout_seconds"])
            except subprocess.TimeoutExpired:
                put_new(result_path, {"schema": "elliptic-curves.sealed-two-quadric-solver-result.v1", "token": case["token"],
                                       "sealed_input_sha256": digest(sealed), "engine": protocol["solver"]["engine"],
                                       "status": "SUBPROCESS_TIMEOUT", "wall_seconds": protocol["solver"]["subprocess_timeout_seconds"]})
        solved = read(result_path)
        assert solved["token"] == case["token"] and solved["sealed_input_sha256"] == digest(sealed)
        verified = map_solution(case, solved, panel)
        records.append({key: case[key] for key in ("token", "cohort", "member", "public_word_mod_2", "alpha_coefficients", "norm_square_root", "integral_quadrics")} |
                       {"sealed_solver_input": str(sealed.relative_to(ROOT)), "sealed_solver_input_sha256": digest(sealed),
                        "solver_result": str(result_path.relative_to(ROOT)), "solver_result_sha256": digest(result_path), **verified})
        print(case["token"], solved["status"], len(verified["exact_rational_cover_points"]), flush=True)
    recovered_nontrivial = [row for row in records if row["cohort"] == "recovered_control" and row["member"] != 0 and row["exact_rational_cover_points"]]
    missing = [row for row in records if row["cohort"] == "missing_strict_coset" and row["exact_rational_cover_points"]]
    status = "PASS" if recovered_nontrivial and missing else "PARTIAL_NO_ADVANCE"
    output = {
        "schema": "elliptic-curves.curve302-strict-cover-hidden-benchmark.v1", "status": status,
        "advance_gate_passed": status == "PASS", "baseline": "M24", "software": {"sage": SAGE_VERSION},
        "protocol": str(PROTOCOL.relative_to(ROOT)), "protocol_sha256": digest(PROTOCOL),
        "bindings": {str(path.relative_to(ROOT)): digest(path) for path in (FILTER, SPAN / "result.json", PARENT, PUBLIC_SOURCE, BUILDER, Path(__file__))},
        "cases": records,
        "comparison": {"representatives_per_cohort": 8,
                       "recovered_control_nontrivial_solution_count": len(recovered_nontrivial),
                       "missing_strict_coset_solution_count": len(missing)},
        "boundary": protocol["boundaries"],
        "next_action": "RUN_GENERIC_ONLY_PRINCIPAL_RELATION_CONSTRUCTOR" if status == "PASS" else "STOP_NO_CONSTRUCTOR",
    }
    put_new(OUTPUT, output)
    print(status, output["comparison"], flush=True)


def check():
    panel = build_panels()
    stored = read(OUTPUT)
    assert stored["bindings"] == {str(path.relative_to(ROOT)): digest(path) for path in (FILTER, SPAN / "result.json", PARENT, PUBLIC_SOURCE, BUILDER, Path(__file__))}
    assert stored["baseline"] == "M24" and len(stored["cases"]) == 16
    for expected, got in zip(panel["cases"], stored["cases"]):
        assert all(got[key] == expected[key] for key in ("token", "cohort", "member", "public_word_mod_2", "alpha_coefficients", "norm_square_root", "integral_quadrics"))
        sealed = ROOT / got["sealed_solver_input"]
        result = ROOT / got["solver_result"]
        assert digest(sealed) == got["sealed_solver_input_sha256"] and digest(result) == got["solver_result_sha256"]
        assert map_solution(expected, read(result), panel) == {key: got[key] for key in ("solver_status", "exact_rational_cover_points")}
    print("PASS exact replay of hidden E302 strict-cover benchmark", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capture", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--solver", action="store_true")
    parser.add_argument("--sealed", type=Path)
    args = parser.parse_args()
    assert sum((args.capture, args.check, args.solver)) == 1
    if args.solver:
        assert args.sealed is not None
        solver(args.sealed)
    elif args.capture:
        capture()
    else:
        check()


if __name__ == "__main__":
    main()
